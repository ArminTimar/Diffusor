"""Forward models and time fits for diffusion-driven isotope fractionation.

Two models are offered, because the physics differs:

:class:`DiluteIsotopeModel`
    A trace element (Li) diffusing in a host of fixed major-element composition.
    Each isotope diffuses on its own with D_m = D_ref (m_ref/m)^beta, using any
    scalar law of the registry for the reference isotope. This is the
    single-species model; Richter et al. (2014, 2017) needed a two-site model to
    reproduce the shape of their Li isotope profiles, and the beta of a
    two-site fit is smaller than the beta a single-site model needs (2014,
    p. 363), so match the beta to the model.

:class:`CoupledFeMgIsotopeModel`
    Fe-Mg interdiffusion in olivine with the three Mg and four Fe isotopes as
    seven components of one ideal ionic solution (Oeser et al. 2026, eqs 3-4).
    The isotopes are coupled through the exchange: modelling 54Fe/56Fe and
    24Mg/26Mg as separate binaries gives different beta values (their Fig. 11).
    The element tracer coefficients come either from the Oeser et al. (2026)
    tracer laws (reference isotopes 57Fe and 25Mg, the doped isotopes of the
    experiments) or from any olivine Fe-Mg interdiffusion law together with a
    ratio D*Fe/D*Mg (reference isotopes 56Fe and 24Mg), using eq. 6 backwards.

Both fit one duration to the concentration and delta profiles together; every
profile needs its own uncertainties because the units differ.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import least_squares, minimize_scalar

from ..coefficients.base import Conditions, DiffusionCoefficient
from ..coefficients.isotopes import ABUNDANCES, DELTA_LABEL, delta_value, mass_factor, split_isotopes
from ..coefficients.transport import ideal_ionic_matrix
from ..solvers.boundary import BoundaryCondition
from ..solvers.convolution import gaussian_convolve
from ..solvers.geometry import Geometry, suggest_grid
from ..solvers.history import ThermalHistory
from ..solvers.initial import InitialCondition
from ..solvers.multicomponent import solve_multicomponent
from ..thermo.buffers import log_fo2_from_delta
from ..thermo.units import human_time
from .fit import T_MAX_DEFAULT, T_MIN_DEFAULT
from .multicomponent import resolution_warnings
from .model import DiffusionModel, M2_PER_S_TO_UM2_PER_S
from .objective import FitStatistics, statistics

OFFSET = "__log10_D_offset__"


def _with_offset(overrides, extra):
    ov = dict(overrides or {})
    ov[OFFSET] = ov.get(OFFSET, 0.0) + extra
    return ov


# ---------------------------------------------------------------------------
@dataclass
class DiluteIsotopeModel:
    """Isotopes of a dilute element, each diffusing with its own scalar D."""
    model: DiffusionModel
    element: str
    beta: float
    delta_left: float = 0.0          # left plateau (core for a plateau-rim initial state)
    delta_right: float = 0.0         # right plateau (rim)
    reference_mass: Optional[int] = None

    def __post_init__(self):
        if self.element not in ABUNDANCES:
            raise ValueError(f"no isotope data for {self.element}; available: {sorted(ABUNDANCES)}")
        if self.model.composition_dependent and self.model.comp_key:
            raise ValueError("splitting into isotopes needs a D that does not depend on the profile's own "
                             "concentration; this law makes D a function of the modelled composition")
        if not np.isfinite(self.beta) or self.beta < 0:
            raise ValueError("beta must be finite and nonnegative")
        ab = ABUNDANCES[self.element]
        self.reference_mass = self.reference_mass or max(ab, key=ab.get)
        if self.reference_mass not in ab:
            raise ValueError(f"{self.element} has no isotope of mass {self.reference_mass}")
        kind = self.model.initial.kind
        if kind not in ("step", "plateau_rim") and self.delta_left != self.delta_right:
            raise ValueError("different left and right delta values need a step or plateau-rim initial state")

    def _isotope_model(self, mass: int) -> DiffusionModel:
        fl = split_isotopes(self.element, 1.0, self.delta_left)[mass]
        fr = split_isotopes(self.element, 1.0, self.delta_right)[mass]
        m = copy.copy(self.model)
        m.initial = copy.deepcopy(self.model.initial)
        p = m.initial.params
        if m.initial.kind == "step":
            p["C_left"], p["C_right"] = p["C_left"] * fl, p["C_right"] * fr
        elif m.initial.kind == "plateau_rim":
            p["C_core"], p["C_rim"] = p["C_core"] * fl, p["C_rim"] * fr
        elif m.initial.kind == "multi_step":
            p["values"] = [v * fl for v in p["values"]]
        elif m.initial.kind == "table":
            p["C_table"] = np.asarray(p["C_table"], dtype=float) * fl
        else:
            raise ValueError(f"initial state {m.initial.kind!r} is not supported for isotopes")

        def scaled(bc: BoundaryCondition, f):
            if bc.kind != "dirichlet":
                return bc
            v = bc.value
            return BoundaryCondition("dirichlet", (lambda t, v=v: v(t) * f) if callable(v) else v * f)
        m.bc_left, m.bc_right = scaled(self.model.bc_left, fl), scaled(self.model.bc_right, fr)
        return m

    def isotope_profiles(self, t_seconds, x_out, overrides=None) -> Dict[int, np.ndarray]:
        out = {}
        for mass in ABUNDANCES[self.element]:
            off = np.log10(mass_factor(self.reference_mass, mass, self.beta))
            out[mass] = self._isotope_model(mass).profile(t_seconds, x_out, _with_offset(overrides, off))
        return out

    def profiles(self, t_seconds, x_out, overrides=None) -> Tuple[np.ndarray, np.ndarray]:
        """Total concentration and principal delta value (per mil) at ``x_out``."""
        iso = self.isotope_profiles(t_seconds, x_out, overrides)
        return sum(iso.values()), delta_value(self.element, iso)

    def warnings(self, t_seconds=None) -> List[str]:
        return self.model.warnings(t_seconds) + [
            f"Single-species isotope model with beta = {self.beta:g} for {self.element}; D of "
            f"{self.reference_mass}{self.element} from {self.model.coefficient.key}. beta depends on the "
            "diffusion model it was fitted with (Richter et al. 2014)."]


# ---------------------------------------------------------------------------
MG = (24, 25, 26)
FE = (54, 56, 57, 58)
COMPONENTS = tuple(f"{m}Mg" for m in MG) + tuple(f"{m}Fe" for m in FE)


@dataclass
class CoupledFeMgIsotopeModel:
    """Seven-isotope Fe-Mg interdiffusion in olivine (Oeser et al. 2026, eqs 3-4).

    Concentrations are site fractions on the octahedral (M) site, X_Fe = Fe/(Fe+Mg).
    The initial state is given for X_Fe; the delta values change at the same
    position as X_Fe for a step or plateau-rim initial state.
    """
    conditions: Conditions
    initial_XFe: InitialCondition
    beta_fe: float
    beta_mg: float
    tracer_fe: Optional[DiffusionCoefficient] = None
    tracer_mg: Optional[DiffusionCoefficient] = None
    interdiffusion: Optional[DiffusionCoefficient] = None
    ratio_fe_mg: Optional[float] = None
    delta_fe: Tuple[float, float] = (0.0, 0.0)     # delta56Fe left, right
    delta_mg: Tuple[float, float] = (0.0, 0.0)     # delta26Mg left, right
    geometry: Geometry = field(default_factory=Geometry)
    left_fixed: bool = False
    right_fixed: bool = False
    history: Optional[ThermalHistory] = None
    fo2_buffer: Optional[Tuple[str, float]] = None
    beam_sigma_um: float = 0.0
    n_nodes: int = 201
    min_steps: int = 150
    log10_offset: float = 0.0
    reference_masses: Tuple[int, int] = (0, 0)      # (Fe, Mg); set from the source if 0

    def __post_init__(self):
        if (self.tracer_fe is None) == (self.interdiffusion is None):
            raise ValueError("give either the Fe and Mg tracer laws or an interdiffusion law with a ratio")
        if self.tracer_fe is not None:
            if self.tracer_mg is None or self.tracer_fe.kind != "tracer" or self.tracer_mg.kind != "tracer":
                raise ValueError("tracer_fe and tracer_mg must both be tracer laws")
            default_ref = (57, 25)        # the doped isotopes the tracer laws were measured on
        else:
            if self.interdiffusion.kind != "interdiffusion" or self.interdiffusion.species != "Fe-Mg":
                raise ValueError("the interdiffusion law must be an Fe-Mg interdiffusion law")
            if self.ratio_fe_mg is None or not self.ratio_fe_mg > 0:
                raise ValueError("an interdiffusion law needs a positive ratio D*Fe/D*Mg")
            default_ref = (56, 24)
        if self.reference_masses == (0, 0):
            self.reference_masses = default_ref
        if self.initial_XFe.kind not in ("step", "plateau_rim"):
            raise ValueError("the coupled isotope model needs a step or plateau-rim initial state")

    # -- coefficients ---------------------------------------------------------
    def conditions_at(self, T_K):
        cond = self.conditions.replace(T_K=T_K)
        if self.fo2_buffer is not None and T_K != self.conditions.T_K:
            cond = cond.replace(log_fo2_bar=log_fo2_from_delta(self.fo2_buffer[0], self.fo2_buffer[1], T_K, cond.P_Pa))
        return cond

    def element_tracers(self, XFe, T_K) -> Tuple[np.ndarray, np.ndarray]:
        """D*Fe and D*Mg (m2/s) of the reference isotopes at X_Fe."""
        XFe = np.asarray(XFe, dtype=float)
        cond = self.conditions_at(T_K)
        cond = cond.replace(X=dict(cond.X, XFe=XFe))
        scale = 10.0 ** self.log10_offset
        if self.tracer_fe is not None:
            return (np.asarray(self.tracer_fe.D(cond), float) * scale,
                    np.asarray(self.tracer_mg.D(cond), float) * scale)
        D = np.asarray(self.interdiffusion.D(cond), float) * scale
        r = self.ratio_fe_mg
        D_mg = D * ((1.0 - XFe) + XFe * r) / r          # eq. 6 solved for D*Mg at D*Fe = r D*Mg
        return r * D_mg, D_mg

    def _beta_masses(self):
        ref_fe, ref_mg = self.reference_masses
        return ([mass_factor(ref_mg, m, self.beta_mg) for m in MG]
                + [mass_factor(ref_fe, m, self.beta_fe) for m in FE])

    def _D_func(self):
        factors = np.array(self._beta_masses())[:, None]

        def D_func(C, T):
            X = np.vstack([C, 1.0 - C.sum(axis=0)])
            X = np.clip(X, 1e-12, None)
            XFe = X[3:].sum(axis=0) / X.sum(axis=0)
            dfe, dmg = self.element_tracers(XFe, T)
            Dt = np.vstack([np.repeat(dmg[None], 3, 0), np.repeat(dfe[None], 4, 0)]) * factors
            return ideal_ionic_matrix(Dt, X, dependent=-1) * M2_PER_S_TO_UM2_PER_S
        return D_func

    # -- state ----------------------------------------------------------------
    def _side(self, x):
        p = self.initial_XFe.params
        x0 = p.get("x0", p.get("rim_start"))
        return np.asarray(x, dtype=float) > x0

    def _plateaus(self):
        p = self.initial_XFe.params
        if self.initial_XFe.kind == "step":
            return p["C_left"], p["C_right"]
        return p["C_core"], p["C_rim"]

    def _isotopes(self, XFe, d_fe, d_mg) -> np.ndarray:
        fe = split_isotopes("Fe", XFe, d_fe)
        mg = split_isotopes("Mg", 1.0 - np.asarray(XFe, dtype=float), d_mg)
        return np.array([mg[m] for m in MG] + [fe[m] for m in FE])

    def initial_state(self, x) -> np.ndarray:
        XFe = self.initial_XFe.evaluate(x)
        right = self._side(x)
        d_fe = np.where(right, self.delta_fe[1], self.delta_fe[0])
        d_mg = np.where(right, self.delta_mg[1], self.delta_mg[0])
        return self._isotopes(XFe, d_fe, d_mg)

    def profiles(self, t_seconds, x_out) -> Dict[str, np.ndarray]:
        """X_Fe, delta56Fe and delta26Mg (per mil) at ``x_out``."""
        x_out = np.asarray(x_out, dtype=float)
        x = suggest_grid(x_out, self.n_nodes)
        full = self.initial_state(x)
        lo, hi = self._plateaus()
        left = self._isotopes(lo, self.delta_fe[0], self.delta_mg[0])[:-1] if self.left_fixed else None
        right = self._isotopes(hi, self.delta_fe[1], self.delta_mg[1])[:-1] if self.right_fixed else None
        hist = (self.history.shifted_to_end(t_seconds) if self.history is not None
                else ThermalHistory.isothermal(self.conditions.T_K, t_seconds))
        res = solve_multicomponent(x, full[:-1], self._D_func(), t_seconds, m=self.geometry.m,
                                   left=left, right=right, history=hist, min_steps=self.min_steps)
        C = np.vstack([res.C_final, 1.0 - res.C_final.sum(axis=0)])
        if self.beam_sigma_um > 0:
            C = np.array([gaussian_convolve(x, c, self.beam_sigma_um) for c in C])
        C = np.array([np.interp(x_out, x, c) for c in C])
        mg = {m: C[k] for k, m in enumerate(MG)}
        fe = {m: C[3 + k] for k, m in enumerate(FE)}
        XFe = sum(fe.values()) / (sum(fe.values()) + sum(mg.values()))
        return {"XFe": XFe, "d56Fe": delta_value("Fe", fe), "d26Mg": delta_value("Mg", mg)}

    def warnings(self) -> List[str]:
        w = []
        cond = self.conditions.replace(X=dict(self.conditions.X, XFe=np.array(self._plateaus(), float)))
        if self.tracer_fe is not None:
            w += self.tracer_fe.check_conditions(cond) + self.tracer_mg.check_conditions(cond)
        else:
            w += self.interdiffusion.check_conditions(cond)
            w.append(f"D*Fe/D*Mg = {self.ratio_fe_mg:g} is an input; Oeser et al. (2026, Table 3) found 0.5-5 "
                     "at 1150-1300 C. The isotope profiles depend on it (their Fig. 11).")
        w.append(f"Coupled seven-isotope model (Oeser et al. 2026) with beta_Fe = {self.beta_fe:g} and beta_Mg = "
                 f"{self.beta_mg:g}; reference isotopes {self.reference_masses[0]}Fe and {self.reference_masses[1]}Mg.")
        return list(dict.fromkeys(w))


# ---------------------------------------------------------------------------
@dataclass
class IsotopeFitResult:
    t_seconds: float
    x: np.ndarray
    observed: Dict[str, np.ndarray]
    sigma: Dict[str, np.ndarray]
    predicted: Dict[str, np.ndarray]
    stats: FitStatistics
    per_profile: Dict[str, FitStatistics]
    beta: Optional[float]
    scan_times: np.ndarray
    scan_chi2: np.ndarray
    warnings: List[str] = field(default_factory=list)
    labels: Dict[str, str] = field(default_factory=dict)
    model: object = None

    def summary(self) -> str:
        lines = [f"t = {human_time(self.t_seconds)} ({self.t_seconds:.4g} s), "
                 f"fitted to {', '.join(self.labels.get(k, k) for k in self.observed)} jointly",
                 self.stats.describe()]
        if self.beta is not None:
            lines.append(f"  beta fitted = {self.beta:.3g}")
        for k, s in self.per_profile.items():
            lines.append(f"  {self.labels.get(k, k)}: chi2 = {s.chi2:.3g}, RMSE = {s.rmse:.3g}")
        lines.extend(f"  ! {w}" for w in self.warnings)
        return "\n".join(lines)


def _check(observed, sigma, x):
    obs, sig = {}, {}
    for k, v in observed.items():
        v = np.asarray(v, dtype=float)
        s = sigma.get(k) if sigma else None
        if s is None:
            raise ValueError(f"{k} needs uncertainties: concentration and delta values have different units")
        s = np.broadcast_to(np.asarray(s, dtype=float), v.shape).copy()
        if v.shape != x.shape or np.any(~np.isfinite(v)) or np.any(~np.isfinite(s)) or np.any(s <= 0):
            raise ValueError(f"{k} needs finite values and positive uncertainties on every point")
        obs[k], sig[k] = v, s
    if not obs:
        raise ValueError("nothing to fit")
    return obs, sig


def _fit(predict, x, obs, sig, *, t_min, t_max, scan_points, progress=None):
    """Log-time scan then bounded refinement of the joint chi2; returns the best log10 t."""
    def resid(t, beta=None):
        pred = predict(t, beta)
        return np.concatenate([(obs[k] - pred[k]) / sig[k] for k in obs]), pred

    logts = np.linspace(np.log10(t_min), np.log10(t_max), scan_points)
    chis = np.empty(scan_points)
    for i, lt in enumerate(logts):
        if progress is not None and progress(0.8 * i / scan_points):
            raise InterruptedError("fit cancelled")
        try:
            chis[i] = float(np.sum(resid(10.0 ** lt)[0] ** 2))
        except Exception:
            chis[i] = np.inf
    i0 = int(np.argmin(chis))
    a, b = logts[max(i0 - 1, 0)], logts[min(i0 + 1, scan_points - 1)]
    opt = minimize_scalar(lambda lt: float(np.sum(resid(10.0 ** lt)[0] ** 2)),
                          bounds=(a, b), method="bounded", options={"xatol": 1e-6})
    logt = float(opt.x) if opt.fun <= chis[i0] else float(logts[i0])
    return logt, logts, chis, resid


def fit_dilute_isotopes(model: DiluteIsotopeModel, x, concentration, sigma_concentration,
                        delta, sigma_delta, *, fit_beta: bool = False, t_min=T_MIN_DEFAULT,
                        t_max=T_MAX_DEFAULT, scan_points: int = 60, progress=None) -> IsotopeFitResult:
    """Fit one duration to a concentration profile and its delta profile (either may be None)."""
    x = np.asarray(x, dtype=float)
    observed = {k: v for k, v in (("C", concentration), ("delta", delta)) if v is not None}
    obs, sig = _check(observed, {"C": sigma_concentration, "delta": sigma_delta}, x)

    def predict(t, beta=None):
        m = model if beta is None else _with_beta(model, beta)
        C, d = m.profiles(t, x)
        return {"C": C, "delta": d}

    logt, logts, chis, resid = _fit(predict, x, obs, sig, t_min=t_min, t_max=t_max,
                                    scan_points=scan_points, progress=progress)
    beta = None
    if fit_beta:
        sol = least_squares(lambda v: resid(10.0 ** v[0], v[1])[0], [logt, model.beta],
                            bounds=([np.log10(t_min), 0.0], [np.log10(t_max), 1.0]), xtol=1e-8, max_nfev=200)
        logt, beta = float(sol.x[0]), float(sol.x[1])
    t = 10.0 ** logt
    r, pred = resid(t, beta)
    stats = statistics(np.concatenate([obs[k] for k in obs]), np.concatenate([pred[k] for k in obs]),
                       np.concatenate([sig[k] for k in obs]), n_params=1 + (beta is not None))
    per = {k: statistics(obs[k], pred[k], sig[k], n_params=0) for k in obs}
    used = model if beta is None else _with_beta(model, beta)
    labels = {"C": f"{model.element} concentration", "delta": DELTA_LABEL[model.element]}
    if progress is not None:
        progress(1.0)
    return IsotopeFitResult(t, x, obs, sig, pred, stats, per, beta, 10.0 ** logts, chis,
                            used.warnings(t) + resolution_warnings(logts, chis, logt, stats.chi2)
                            + ([stats.dof_warning()] if stats.dof_warning() else []), labels, used)


def _with_beta(model, beta):
    m = copy.copy(model)
    m.beta = float(beta)
    return m


def fit_coupled_isotopes(model: CoupledFeMgIsotopeModel, x, observed: Dict[str, Sequence[float]],
                         sigma: Dict[str, Sequence[float]], *, t_min=T_MIN_DEFAULT, t_max=T_MAX_DEFAULT,
                         scan_points: int = 40, progress=None) -> IsotopeFitResult:
    """Fit one duration to any of ``XFe``, ``d56Fe`` and ``d26Mg`` profiles.

    A profile measured as forsterite content is converted to X_Fe by the
    caller (X_Fe = 1 - Fo for a fraction).
    """
    x = np.asarray(x, dtype=float)
    unknown = set(observed) - {"XFe", "d56Fe", "d26Mg"}
    if unknown:
        raise ValueError(f"unknown profiles {sorted(unknown)}; use XFe, d56Fe, d26Mg")
    obs, sig = _check(observed, sigma, x)

    def predict(t, beta=None):
        return model.profiles(t, x)

    logt, logts, chis, resid = _fit(predict, x, obs, sig, t_min=t_min, t_max=t_max,
                                    scan_points=scan_points, progress=progress)
    t = 10.0 ** logt
    r, pred = resid(t)
    stats = statistics(np.concatenate([obs[k] for k in obs]), np.concatenate([pred[k] for k in obs]),
                       np.concatenate([sig[k] for k in obs]), n_params=1)
    per = {k: statistics(obs[k], pred[k], sig[k], n_params=0) for k in obs}
    if progress is not None:
        progress(1.0)
    return IsotopeFitResult(t, x, obs, sig, pred, stats, per, None, 10.0 ** logts, chis,
                            model.warnings() + resolution_warnings(logts, chis, logt, stats.chi2)
                            + ([stats.dof_warning()] if stats.dof_warning() else []),
                            {"XFe": "X_Fe", "d56Fe": "δ56Fe", "d26Mg": "δ26Mg"}, model)
