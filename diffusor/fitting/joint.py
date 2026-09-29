"""Shared-duration inference for independent scalar diffusion profiles.

The objective is sum(((C_observed-C_model)/sigma)**2) across profiles.
Each profile retains its own geometry, boundary conditions and coefficient.
This is a joint likelihood, not a multicomponent flux or reaction model.
"""
from dataclasses import dataclass
from typing import Sequence

import numpy as np
from scipy.optimize import minimize_scalar

from .fit import FitResult, T_MIN_DEFAULT, T_MAX_DEFAULT
from .model import DiffusionModel
from .objective import statistics


@dataclass
class ProfileConstraint:
    name: str
    model: DiffusionModel
    x: np.ndarray
    concentration: np.ndarray
    sigma: np.ndarray

    def __post_init__(self):
        self.x = np.asarray(self.x, dtype=float).copy()
        self.concentration = np.asarray(self.concentration, dtype=float).copy()
        if self.sigma is None:
            raise ValueError("Joint fitting requires measurement uncertainties for every profile; different concentration units cannot be weighted equally.")
        self.sigma = np.broadcast_to(np.asarray(self.sigma, dtype=float), self.x.shape).copy()
        if (self.x.ndim != 1 or self.x.size < 3 or self.concentration.shape != self.x.shape
                or not np.all(np.isfinite(self.x)) or not np.all(np.isfinite(self.concentration))
                or not np.all(np.isfinite(self.sigma)) or np.any(self.sigma <= 0)
                or np.ptp(self.x) == 0):
            raise ValueError("Each profile requires at least three finite paired measurements and strictly positive finite uncertainties.")
        if not self.name.strip():
            raise ValueError("Each profile needs a name.")


@dataclass
class JointFitResult:
    t_seconds: float
    names: list[str]
    profiles: list[FitResult]
    chi2: float
    dof: int
    scan_times: np.ndarray
    scan_chi2: np.ndarray
    warnings: list[str]

    @property
    def reduced_chi2(self):
        return self.chi2 / self.dof


def fit_joint_time(profiles: Sequence[ProfileConstraint], *,
                   t_min=T_MIN_DEFAULT, t_max=T_MAX_DEFAULT, scan_points=80,
                   progress=None) -> JointFitResult:
    """Fit one duration to two or more independent measured profiles.

    Known Gaussian errors are required in each profile's own concentration
    units. All parameters other than time remain fixed. Coarse scanning and
    refinement of every sampled local minimum reduce sensitivity to initial
    guesses. This does not claim a global optimum below scan resolution.
    ``progress(fraction)`` can return True to cancel.
    """
    profiles = list(profiles)
    if len(profiles) < 2:
        raise ValueError("A joint study requires at least two profiles.")
    if len({p.name for p in profiles}) != len(profiles):
        raise ValueError("Profile names must be unique.")
    if not (np.isfinite(t_min) and np.isfinite(t_max) and 0 < t_min < t_max):
        raise ValueError("Time bounds must be finite, positive and increasing.")
    if scan_points < 5:
        raise ValueError("At least five scan points are required.")

    def objective(logt):
        total = 0.0
        for p in profiles:
            pred = p.model.profile(10.0**logt, p.x)
            if not np.all(np.isfinite(pred)):
                raise ValueError(f"Non-finite prediction in profile {p.name!r}.")
            total += float(np.sum(((p.concentration - pred) / p.sigma)**2))
        return total

    logtimes = np.linspace(np.log10(t_min), np.log10(t_max), scan_points)
    chi = np.empty(scan_points)
    for i, logt in enumerate(logtimes):
        if progress is not None and progress(i / scan_points):
            raise InterruptedError("Joint fit cancelled.")
        chi[i] = objective(logt)
    candidates = [(chi[0], logtimes[0]), (chi[-1], logtimes[-1])]
    for i in range(1, scan_points-1):
        if chi[i] <= chi[i-1] and chi[i] <= chi[i+1]:
            if progress is not None and progress(.95):
                raise InterruptedError("Joint fit cancelled.")
            opt = minimize_scalar(objective, bounds=(logtimes[i-1], logtimes[i+1]),
                                  method="bounded", options={"xatol": 1e-8})
            if not opt.success:
                raise RuntimeError(f"Joint time refinement failed: {opt.message}")
            candidates.append((float(opt.fun), float(opt.x)))
    chi_best, log_best = min(candidates)
    t_best = 10.0**log_best
    warnings = ["Independent scalar profiles share one fitted duration. Cross-element fluxes, site exchange and reactions are not represented.",
                "The likelihood assumes independent Gaussian measurement errors. Temperature, diffusion coefficients, initial states and boundaries are held fixed; systematic uncertainty is not propagated.",
                "A shared duration is a scientific hypothesis; inspect every profile's residuals for incompatible histories."]
    if abs(log_best-logtimes[0]) < 1e-5 or abs(log_best-logtimes[-1]) < 1e-5:
        warnings.append("The optimum is at a time bound; the duration is not resolved within the selected range.")
    results = []
    for p in profiles:
        ok, why = p.model.can_use_analytical()
        prediction = p.model.profile(t_best, p.x)
        caveats = p.model.warnings(t_best) + ["Time estimated jointly with other profiles. Per-profile degrees of freedom are descriptive; use the joint N-1 degrees of freedom for inference."]
        results.append(FitResult(t_best, p.model, p.x, p.concentration, p.sigma,
                                 prediction, statistics(p.concentration, prediction, p.sigma, n_params=0),
                                 free={"t": t_best}, t_bounds_scan=(t_min, t_max), warnings=caveats,
                                 route=("analytical: " if ok else "Crank-Nicolson: ") + why))
    if progress is not None:
        progress(1.)
    return JointFitResult(t_best, [p.name for p in profiles], results, chi_best,
                          sum(p.x.size for p in profiles)-1, 10**logtimes, chi, warnings)
