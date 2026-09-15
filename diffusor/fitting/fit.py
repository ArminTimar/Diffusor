"""Least-squares fitting of the diffusion time (and optional nuisance parameters).

Strategy: a coarse logarithmic scan in time to locate the global minimum,
then a bounded refinement with :func:`scipy.optimize.least_squares`.  The scan
matters because chi2(t) is very flat at long times and has a single broad
minimum -- a local optimiser started in the wrong decade converges to the edge.

Free parameters
---------------
``t``            always
``x0``           interface position (only if the initial condition is a step)
``C_left``, ``C_right``   the plateau compositions
``beam_sigma``   analytical resolution, if the user wants it fitted

Weighting is ``1/sigma`` per point (Bevington-style weighted least squares);
with no uncertainty column the weights are uniform and the reported chi2 is
then only a relative measure.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import least_squares

from ..thermo.units import human_time
from .model import DiffusionModel
from .objective import FitStatistics, residuals, statistics

FREE_PARAMETERS = ("t", "x0", "C_left", "C_right", "beam_sigma")


@dataclass
class FitResult:
    t_seconds: float
    model: DiffusionModel
    x_data: np.ndarray
    C_data: np.ndarray
    sigma: Optional[np.ndarray]
    C_model: np.ndarray
    stats: FitStatistics
    free: Dict[str, float] = field(default_factory=dict)
    t_bounds_scan: Tuple[float, float] = (0.0, 0.0)
    success: bool = True
    message: str = ""
    warnings: List[str] = field(default_factory=list)
    route: str = ""

    @property
    def t_years(self) -> float:
        from ..constants import SEC_PER_YEAR
        return self.t_seconds / SEC_PER_YEAR

    def summary(self) -> str:
        lines = [f"t = {human_time(self.t_seconds)} ({self.t_seconds:.4g} s)",
                 f"solver: {self.route}",
                 self.stats.describe()]
        for k, v in self.free.items():
            if k != "t":
                lines.append(f"  fitted {k} = {v:.4g}")
        for w in self.warnings:
            lines.append(f"  ! {w}")
        return "\n".join(lines)


def _apply_free(model: DiffusionModel, values: Dict[str, float]) -> DiffusionModel:
    """Return a shallow copy of the model with nuisance parameters applied."""
    import copy
    if not any(k in values for k in ("x0", "C_left", "C_right", "beam_sigma")):
        return model
    m = copy.copy(model)
    m.initial = copy.deepcopy(model.initial)
    for k in ("x0", "C_left", "C_right"):
        if k in values:
            m.initial.params[k] = float(values[k])
    if "beam_sigma" in values:
        m.beam_sigma_um = max(0.0, float(values["beam_sigma"]))
    return m


def predict(model: DiffusionModel, t_seconds: float, x_data, free=None, overrides=None):
    m = _apply_free(model, free or {})
    return m.profile(max(t_seconds, 0.0), x_data, overrides)


def scan_time(model: DiffusionModel, x_data, C_data, sigma=None,
              t_min: float = 1.0e2, t_max: float = 3.2e12, n: int = 60,
              free=None, overrides=None) -> Tuple[float, np.ndarray, np.ndarray]:
    """Coarse logarithmic scan of chi2 against time.

    Default range is 100 s to about 100 kyr, which brackets everything from
    syn-eruptive ascent to long crustal residence (Costa et al. 2020).
    """
    ts = np.logspace(np.log10(t_min), np.log10(t_max), n)
    chi = np.empty(n)
    for i, t in enumerate(ts):
        try:
            Cm = predict(model, float(t), x_data, free, overrides)
            chi[i] = float(np.sum(residuals(C_data, Cm, sigma) ** 2))
        except Exception:
            chi[i] = np.inf
    return float(ts[int(np.argmin(chi))]), ts, chi


def fit_time(model: DiffusionModel, x_data, C_data, sigma=None,
             free_parameters: Sequence[str] = ("t",),
             t_min: float = 1.0e2, t_max: float = 3.2e12,
             overrides=None, verbose: bool = False,
             t_guess: Optional[float] = None, scan_points: int = 60) -> FitResult:
    """Fit the diffusion time, optionally with nuisance parameters.

    ``t_guess`` skips the wide logarithmic scan and searches only two decades
    around the guess.  Monte Carlo draws use this (seeded from the best fit),
    which cuts the number of forward solves per draw by roughly five times
    without changing the answer, because each draw perturbs the conditions only
    slightly.
    """
    x_data = np.asarray(x_data, dtype=float)
    C_data = np.asarray(C_data, dtype=float)
    sigma = None if sigma is None else np.asarray(sigma, dtype=float)
    for p in free_parameters:
        if p not in FREE_PARAMETERS:
            raise ValueError(f"'{p}' is not a fittable parameter; choose from {FREE_PARAMETERS}")

    if t_guess is not None and np.isfinite(t_guess) and t_guess > 0:
        lo_s = max(t_min, t_guess / 30.0)
        hi_s = min(t_max, t_guess * 30.0)
        t0, ts, chis = scan_time(model, x_data, C_data, sigma, lo_s, hi_s,
                                 n=max(9, scan_points // 5), overrides=overrides)
    else:
        t0, ts, chis = scan_time(model, x_data, C_data, sigma, t_min, t_max,
                                 n=scan_points, overrides=overrides)

    # --- build the parameter vector (log time + nuisance parameters) ---------
    names = [p for p in free_parameters if p != "t"]
    p0 = [np.log10(t0)]
    lo = [np.log10(t_min)]
    hi = [np.log10(t_max)]
    span = float(x_data.max() - x_data.min())
    c_lo, c_hi = float(np.min(C_data)), float(np.max(C_data))
    c_pad = 0.5 * (c_hi - c_lo) + 1e-12
    for nme in names:
        if nme == "x0":
            p0.append(float(model.initial.params.get("x0", np.mean(x_data))))
            lo.append(float(x_data.min())); hi.append(float(x_data.max()))
        elif nme in ("C_left", "C_right"):
            p0.append(float(model.initial.params.get(nme, c_lo if nme == "C_left" else c_hi)))
            lo.append(c_lo - c_pad); hi.append(c_hi + c_pad)
        elif nme == "beam_sigma":
            p0.append(max(model.beam_sigma_um, 0.5))
            lo.append(0.0); hi.append(max(0.1 * span, 1.0))

    def unpack(v):
        t = 10.0 ** v[0]
        free = {n: float(val) for n, val in zip(names, v[1:])}
        return t, free

    def fun(v):
        t, free = unpack(v)
        try:
            Cm = predict(model, t, x_data, free, overrides)
        except Exception:
            return np.full(x_data.size, 1e6)
        return residuals(C_data, Cm, sigma)

    try:
        sol = least_squares(fun, np.array(p0), bounds=(np.array(lo), np.array(hi)),
                            xtol=1e-10, ftol=1e-10, max_nfev=400)
        t_best, free_best = unpack(sol.x)
        success, message = bool(sol.success), str(sol.message)
    except Exception as exc:                      # pragma: no cover - defensive
        t_best, free_best = t0, {}
        success, message = False, f"refinement failed ({exc}); using the scan minimum"

    m_final = _apply_free(model, free_best)
    C_model = m_final.profile(t_best, x_data, overrides)
    stats = statistics(C_data, C_model, sigma, n_params=1 + len(names))
    ok, why = model.can_use_analytical()
    free_best["t"] = t_best
    res = FitResult(
        t_seconds=float(t_best), model=m_final, x_data=x_data, C_data=C_data, sigma=sigma,
        C_model=C_model, stats=stats, free=free_best,
        t_bounds_scan=(float(ts[0]), float(ts[-1])), success=success, message=message,
        warnings=m_final.warnings(t_best),
        route=("analytical (Crank 1975) -- " + why) if ok else ("numerical Crank-Nicolson -- " + why))
    if t_best <= 1.05 * t_min or t_best >= 0.95 * t_max:
        res.warnings.append(
            f"the fitted time is at the edge of the search window "
            f"({human_time(t_min)} to {human_time(t_max)}); widen it or check the data")
    if verbose:
        print(res.summary())
    return res
