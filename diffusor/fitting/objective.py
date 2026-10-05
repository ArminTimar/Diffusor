"""Goodness-of-fit statistics."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import warnings

import numpy as np


@dataclass
class FitStatistics:
    chi2: float
    reduced_chi2: float
    dof: int
    rmse: float
    r_squared: float
    n_points: int
    n_params: int

    def describe(self) -> str:
        return (f"chi2 = {self.chi2:.3g}, reduced chi2 = {self.reduced_chi2:.3g} "
                f"({self.dof} dof), RMSE = {self.rmse:.4g}, R2 = {self.r_squared:.4f}")

    def dof_warning(self) -> Optional[str]:
        """Text to show the user when the fit has no degrees of freedom, else None."""
        if self.dof > 0:
            return None
        return (f"The fit has {self.n_points} points for {self.n_params} fitted parameters "
                f"({self.dof} degrees of freedom), so the reduced chi-squared is undefined (NaN) "
                "and the fit cannot be tested against the uncertainties. Add data or fix a parameter.")


def weights_from_sigma(sigma: Optional[np.ndarray], n: int) -> np.ndarray:
    """1/sigma weights; uniform weights if no uncertainties were supplied."""
    if sigma is None:
        return np.ones(n)
    s = np.asarray(sigma, dtype=float)
    good = np.isfinite(s) & (s > 0)
    if not np.all(good):
        # a NaN, infinite, zero or negative uncertainty is not a usable weight: use the
        # median of the usable ones (or 1 if there are none) and say so
        fill = float(np.median(s[good])) if np.any(good) else 1.0
        warnings.warn(f"{int((~good).sum())} of {s.size} uncertainties are not finite and positive; "
                      f"they were replaced by {fill:.3g} (the median of the valid ones, or 1 if none) "
                      "when weighting the fit.", UserWarning, stacklevel=2)
        s = np.where(good, s, fill)
    return 1.0 / s


def residuals(observed, model, sigma=None) -> np.ndarray:
    """Weighted residuals ``(obs - model)/sigma`` as used by least_squares."""
    observed = np.asarray(observed, dtype=float)
    model = np.asarray(model, dtype=float)
    return (observed - model) * weights_from_sigma(sigma, observed.size)


def statistics(observed, model, sigma=None, n_params: int = 1) -> FitStatistics:
    observed = np.asarray(observed, dtype=float)
    model = np.asarray(model, dtype=float)
    r = residuals(observed, model, sigma)
    chi2 = float(np.sum(r ** 2))
    n = observed.size
    # true degrees of freedom: zero or negative when the fit is exactly or over-determined,
    # in which case the reduced chi-squared is undefined (NaN), not a finite number
    dof = n - n_params
    ss_res = float(np.sum((observed - model) ** 2))
    ss_tot = float(np.sum((observed - np.mean(observed)) ** 2))
    return FitStatistics(
        chi2=chi2, reduced_chi2=chi2 / dof if dof > 0 else float("nan"), dof=dof,
        rmse=float(np.sqrt(ss_res / n)),
        r_squared=1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan"),
        n_points=n, n_params=n_params)
