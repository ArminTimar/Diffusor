"""Goodness-of-fit statistics."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

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


def weights_from_sigma(sigma: Optional[np.ndarray], n: int) -> np.ndarray:
    """1/sigma weights; uniform weights if no uncertainties were supplied."""
    if sigma is None:
        return np.ones(n)
    s = np.asarray(sigma, dtype=float)
    if np.any(s <= 0):
        s = np.where(s > 0, s, np.nanmedian(s[s > 0]) if np.any(s > 0) else 1.0)
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
    dof = max(n - n_params, 1)
    ss_res = float(np.sum((observed - model) ** 2))
    ss_tot = float(np.sum((observed - np.mean(observed)) ** 2))
    return FitStatistics(
        chi2=chi2, reduced_chi2=chi2 / dof, dof=dof,
        rmse=float(np.sqrt(ss_res / n)),
        r_squared=1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan"),
        n_points=n, n_params=n_params)
