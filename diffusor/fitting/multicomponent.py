"""Forward model and time fit for multicomponent (garnet Fe-Mg-Mn-Ca) diffusion.

:class:`MulticomponentModel` turns a duration into the profiles of every
component, using a :class:`~diffusor.coefficients.families.TracerFamily` for
the tracer coefficients, the ideal ionic matrix for their coupling and
:func:`~diffusor.solvers.multicomponent.solve_multicomponent` for transport.
:func:`fit_multicomponent_time` fits one duration to all measured components
at once by weighted least squares, so a component whose profile the cross terms
make non-monotonic (uphill diffusion of Ca, Carlson 2006 Fig. 4) constrains
the time together with the others.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass, field
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import least_squares, minimize_scalar

from ..coefficients.families import TracerFamily
from ..solvers.convolution import gaussian_convolve
from ..solvers.geometry import Geometry, suggest_grid
from ..solvers.history import ThermalHistory
from ..solvers.initial import InitialCondition
from ..solvers.multicomponent import solve_multicomponent
from ..thermo.units import human_time
from .fit import T_MAX_DEFAULT, T_MIN_DEFAULT
from .objective import FitStatistics, statistics

UM2 = 1.0e12        # m2/s -> um2/s


@dataclass
class MulticomponentModel:
    family: TracerFamily
    T_K: float
    initial: Dict[str, InitialCondition]          # one per component of the family
    P_Pa: float = 1.0e5
    options: Dict[str, float] = field(default_factory=dict)
    geometry: Geometry = field(default_factory=Geometry)
    left_fixed: Optional[Dict[str, float]] = None  # fixed rim composition, else closed
    right_fixed: Optional[Dict[str, float]] = None
    history: Optional[ThermalHistory] = None
    beam_sigma_um: float = 0.0
    n_nodes: int = 201
    min_steps: int = 150
    x_grid: Optional[np.ndarray] = None
    dependent: str = ""
    log10_offsets: Dict[str, float] = field(default_factory=dict)   # Monte Carlo perturbations

    def __post_init__(self):
        comps = self.family.components
        if set(self.initial) != set(comps):
            raise ValueError(f"an initial condition is needed for each of {', '.join(comps)}")
        self.dependent = self.dependent or self.family.default_dependent or comps[-1]
        if self.dependent not in comps:
            raise ValueError(f"unknown dependent component {self.dependent!r}")
        for fixed in (self.left_fixed, self.right_fixed):
            if fixed is not None and set(fixed) != set(comps):
                raise ValueError("a fixed boundary needs a value for every component")
        self.family._options(self.options)       # validates option names

    @property
    def independent(self) -> List[str]:
        return [c for c in self.family.components if c != self.dependent]

    def grid(self, x_data) -> np.ndarray:
        return self.x_grid if self.x_grid is not None else suggest_grid(x_data, self.n_nodes)

    def initial_state(self, x) -> Dict[str, np.ndarray]:
        raw = {c: np.asarray(self.initial[c].evaluate(x), dtype=float) for c in self.family.components}
        total = sum(raw.values())
        if np.any(total <= 0):
            raise ValueError("initial mole fractions must have a positive sum everywhere")
        return {c: v / total for c, v in raw.items()}

    def _fixed(self, values):
        if values is None:
            return None
        total = sum(values.values())
        return [values[c] / total for c in self.independent]

    def _matrix(self, T_K):
        def D_func(C, T):
            X = {c: C[i] for i, c in enumerate(self.independent)}
            X[self.dependent] = np.clip(1.0 - C.sum(axis=0), 0.0, None)
            X = {c: np.clip(v, 1e-9, None) for c, v in X.items()}
            return self.family.matrix(T, self.P_Pa, X, self.dependent, self.options,
                                      self.log10_offsets) * UM2
        return D_func

    def profiles(self, t_seconds: float, x_out) -> Dict[str, np.ndarray]:
        """Mole fraction of every component at ``x_out`` after ``t_seconds``."""
        x_out = np.asarray(x_out, dtype=float)
        x = self.grid(x_out)
        X0 = self.initial_state(x)
        C0 = np.array([X0[c] for c in self.independent])
        hist = (self.history.shifted_to_end(t_seconds) if self.history is not None
                else ThermalHistory.isothermal(self.T_K, t_seconds))
        res = solve_multicomponent(x, C0, self._matrix(self.T_K), t_seconds, m=self.geometry.m,
                                   left=self._fixed(self.left_fixed), right=self._fixed(self.right_fixed),
                                   history=hist, min_steps=self.min_steps)
        out = {c: res.C_final[i] for i, c in enumerate(self.independent)}
        out[self.dependent] = 1.0 - res.C_final.sum(axis=0)
        if self.beam_sigma_um > 0:
            out = {c: gaussian_convolve(x, v, self.beam_sigma_um) for c, v in out.items()}
        return {c: np.interp(x_out, x, out[c]) for c in self.family.components}

    def representative_composition(self) -> Dict[str, float]:
        x = np.linspace(-1, 1, 3)
        p = self.initial_state(x)
        return {c: float(np.mean(v)) for c, v in p.items()}

    def warnings(self) -> List[str]:
        w = list(self.family.check(self.T_K, self.P_Pa, self.representative_composition()))
        if self.history is not None:
            for T in np.unique(self.history.temps_K):
                w.extend(self.family.check(float(T), self.P_Pa))
        if self.geometry.kind == "plane":
            w.append("1-D modelling can be biased by 3-D geometry and section orientation "
                     "(Shea et al. 2015, Krimer & Costa 2017).")
        return list(dict.fromkeys(w))

    def describe(self) -> str:
        opts = ", ".join(f"{k} = {v:g}" for k, v in self.family._options(self.options).items())
        return (f"{self.family.label}; T = {self.T_K - 273.15:.1f} C, P = {self.P_Pa / 1e9:.3g} GPa"
                + (f", {opts}" if opts else "") + f"; dependent component {self.dependent}; "
                f"{self.geometry.kind} geometry; left end "
                + ("fixed" if self.left_fixed else "closed") + ", right end "
                + ("fixed" if self.right_fixed else "closed"))


@dataclass
class MulticomponentFitResult:
    t_seconds: float
    model: MulticomponentModel
    x: np.ndarray
    data: Dict[str, np.ndarray]
    sigma: Dict[str, Optional[np.ndarray]]
    fitted: Tuple[str, ...]
    prediction: Dict[str, np.ndarray]
    stats: FitStatistics
    per_component: Dict[str, FitStatistics]
    scan_times: np.ndarray
    scan_chi2: np.ndarray
    free: Dict[str, float] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)
    mc_times: Optional[np.ndarray] = None

    def summary(self) -> str:
        lines = [f"t = {human_time(self.t_seconds)} ({self.t_seconds:.4g} s), fitted to "
                 f"{', '.join(self.fitted)} jointly", self.stats.describe()]
        for c, s in self.per_component.items():
            lines.append(f"  {c}: chi2 = {s.chi2:.3g}, RMSE = {s.rmse:.3g}")
        if self.mc_times is not None and self.mc_times.size:
            lo, med, hi = np.percentile(self.mc_times, [16, 50, 84])
            lines.append(f"  Monte Carlo ({self.mc_times.size} draws): median {human_time(med)}, "
                         f"68 % interval {human_time(lo)} to {human_time(hi)}")
        lines.extend(f"  ! {w}" for w in self.warnings)
        return "\n".join(lines)


def _prepare(data, sigma, fitted):
    data = {c: np.asarray(v, dtype=float) for c, v in data.items()}
    sig = {}
    for c in fitted:
        s = (sigma or {}).get(c)
        sig[c] = None if s is None else np.broadcast_to(np.asarray(s, dtype=float), data[c].shape).copy()
        if sig[c] is not None and (np.any(~np.isfinite(sig[c])) or np.any(sig[c] <= 0)):
            raise ValueError(f"uncertainties of {c} must be finite and positive")
    return data, sig


def _residuals(model, t, x, data, sig, fitted, free=None):
    m = model
    if free:
        m = copy.copy(model)
        m.initial = {c: InitialCondition(ic.kind, dict(ic.params, **({"x0": free["x0"]} if "x0" in ic.params else
                                                                     {"rim_start": free["x0"]} if "rim_start" in ic.params else {})),
                                         ic.description)
                     for c, ic in model.initial.items()}
    pred = m.profiles(t, x)
    r = [(data[c] - pred[c]) / (sig[c] if sig[c] is not None else 1.0) for c in fitted]
    return np.concatenate(r), pred


def fit_multicomponent_time(model: MulticomponentModel, x, data: Mapping[str, Sequence[float]],
                            sigma: Optional[Mapping[str, object]] = None,
                            fitted: Optional[Sequence[str]] = None, *, fit_x0: bool = False,
                            t_min: float = T_MIN_DEFAULT, t_max: float = T_MAX_DEFAULT,
                            scan_points: int = 50, t_guess: Optional[float] = None,
                            progress=None) -> MulticomponentFitResult:
    """Fit one duration (and optionally the interface position) to several components.

    ``data`` maps component names to measured mole fractions on ``x`` (um).
    ``sigma`` gives one-sigma uncertainties per component; without them all
    components are weighted equally in mole-fraction units and chi2 is only a
    relative measure. ``progress(fraction)`` may return True to cancel.
    """
    x = np.asarray(x, dtype=float)
    fitted = tuple(fitted or [c for c in model.family.components if c in data])
    if not fitted or any(c not in model.family.components for c in fitted):
        raise ValueError(f"choose components to fit from {model.family.components}")
    data, sig = _prepare(data, sigma, fitted)
    if any(data[c].shape != x.shape for c in fitted):
        raise ValueError("every fitted component needs one value per distance")

    def chi2(logt):
        r, _ = _residuals(model, 10.0 ** logt, x, data, sig, fitted)
        return float(np.sum(r ** 2))

    if t_guess is not None and t_guess > 0:
        lo, hi = max(np.log10(t_min), np.log10(t_guess) - 1.5), min(np.log10(t_max), np.log10(t_guess) + 1.5)
        scan_points = max(9, scan_points // 4)
    else:
        lo, hi = np.log10(t_min), np.log10(t_max)
    logts = np.linspace(lo, hi, scan_points)
    chis = np.empty(scan_points)
    for i, lt in enumerate(logts):
        if progress is not None and progress(0.8 * i / scan_points):
            raise InterruptedError("fit cancelled")
        try:
            chis[i] = chi2(lt)
        except Exception:
            chis[i] = np.inf
    i0 = int(np.argmin(chis))
    a, b = logts[max(i0 - 1, 0)], logts[min(i0 + 1, scan_points - 1)]
    opt = minimize_scalar(chi2, bounds=(a, b), method="bounded", options={"xatol": 1e-6})
    logt = float(opt.x) if opt.fun <= chis[i0] else float(logts[i0])
    free: Dict[str, float] = {}
    if fit_x0:
        ic = next(iter(model.initial.values()))
        x0 = float(ic.params.get("x0", ic.params.get("rim_start", np.mean(x))))

        def fun(v):
            r, _ = _residuals(model, 10.0 ** v[0], x, data, sig, fitted, {"x0": v[1]})
            return r
        sol = least_squares(fun, [logt, x0], bounds=([np.log10(t_min), x.min()], [np.log10(t_max), x.max()]),
                            xtol=1e-8, ftol=1e-8, max_nfev=200)
        logt, free = float(sol.x[0]), {"x0": float(sol.x[1])}
    t = 10.0 ** logt
    r, pred = _residuals(model, t, x, data, sig, fitted, free or None)
    if free:
        model = copy.copy(model)
        model.initial = {c: InitialCondition(ic.kind, dict(ic.params, **({"x0": free["x0"]} if "x0" in ic.params
                                                                          else {"rim_start": free["x0"]})), ic.description)
                         for c, ic in model.initial.items()}
    obs = np.concatenate([data[c] for c in fitted])
    mod = np.concatenate([pred[c] for c in fitted])
    s_all = None if any(sig[c] is None for c in fitted) else np.concatenate([sig[c] for c in fitted])
    stats = statistics(obs, mod, s_all, n_params=1 + len(free))
    per = {c: statistics(data[c], pred[c], sig[c], n_params=0) for c in fitted}
    w = model.warnings()
    if any(sig[c] is None for c in fitted):
        w.append("No uncertainties were given for every component, so all are weighted equally in mole "
                 "fraction and chi2 is only a relative measure.")
    if abs(logt - np.log10(t_min)) < 0.05 or abs(logt - np.log10(t_max)) < 0.05:
        w.append("The fitted time is at the edge of the search window.")
    if all(sig[c] is not None for c in fitted):      # the test is in units of sigma
        w.extend(resolution_warnings(logts, chis, logt, stats.chi2))
    if progress is not None:
        progress(1.0)
    return MulticomponentFitResult(t, model, x, data, sig, fitted, pred, stats, per,
                                   10.0 ** logts, chis, {"t": t, **free}, w)


def resolution_warnings(logts, chis, logt, chi2_best) -> List[str]:
    """Warn when the scan shows that the duration is not resolved.

    If chi2 stays within 1 of the minimum at every scanned time longer than
    the best fit, the profiles may have relaxed fully: any longer time fits as
    well, and the result is only a lower limit. The same test is applied to
    shorter times (an upper limit).
    """
    w = []
    finite = np.isfinite(chis)
    later = finite & (logts > logt + 0.5)
    earlier = finite & (logts < logt - 0.5)
    if later.any() and np.all(chis[later] <= chi2_best + 1.0):
        w.append("chi2 stays within 1 of its minimum at every longer time scanned: the profiles may have "
                 "homogenised, and the fitted duration is only a lower limit.")
    if earlier.any() and np.all(chis[earlier] <= chi2_best + 1.0):
        w.append("chi2 stays within 1 of its minimum at every shorter time scanned: the profiles show no "
                 "resolvable diffusion, and the fitted duration is only an upper limit.")
    return w


def monte_carlo_multicomponent(result: MulticomponentFitResult, n_draws: int = 50, *,
                               T_sigma_K: float = 0.0, seed: int = 0, progress=None) -> np.ndarray:
    """Refit the duration for draws of temperature and of the tracer coefficients.

    Each draw shifts the temperature (and a cooling path with it) by a normal
    deviate of ``T_sigma_K`` and every component's log10 D* by an independent
    normal deviate of the family's ``sigma_logD`` (no draw if the family has
    none). The measured profiles are refitted with the nominal initial and
    boundary conditions. Returns the fitted times in seconds.
    """
    rng = np.random.default_rng(seed)
    fam = result.model.family
    times = []
    for k in range(n_draws):
        if progress is not None and progress(k / n_draws):
            break
        m = copy.copy(result.model)
        dT = float(rng.normal(0.0, T_sigma_K)) if T_sigma_K > 0 else 0.0
        m.T_K = result.model.T_K + dT
        if result.model.history is not None:
            h = result.model.history
            m.history = ThermalHistory(h.times.copy(), h.temps_K + dT, h.label)
        m.log10_offsets = ({c: float(rng.normal(0.0, fam.sigma_logD)) for c in fam.components}
                           if fam.sigma_logD else {})
        try:
            r = fit_multicomponent_time(m, result.x, result.data, result.sigma, result.fitted,
                                        t_guess=result.t_seconds)
            times.append(r.t_seconds)
        except Exception:
            continue
    out = np.array(times)
    result.mc_times = out
    if progress is not None:
        progress(1.0)
    return out
