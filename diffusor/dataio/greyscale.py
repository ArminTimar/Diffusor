"""Greyscale (BSE) profiles calibrated to composition.

Back-scattered electron intensity varies with mean atomic number, so for a
binary substitution such as Fe-Mg it is a linear proxy for composition over a
narrow compositional range.  This is the basis of the greyscale approach used
by NIDIS (Petrone et al. 2016, Methods, "Rationale of working with greyscale
values of BSE images") and by Morgan et al. (2004); it buys a spatial
resolution far better than an electron microprobe traverse, at the cost of
needing a calibration.

Diffusor never fits raw grey values.  A greyscale profile must be calibrated
against at least two microprobe anchor points, and the calibration residuals
are reported so a non-linear response is visible.  If the calibration is poor,
the timescale that follows is not trustworthy.

The uncertainty on each calibrated point combines the scatter of the grey
values across the averaged raster lines with the standard error of the
calibration itself.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

import numpy as np


@dataclass
class GreyscaleCalibration:
    """Linear or quadratic map from grey value to composition."""
    coefficients: np.ndarray            # highest power first, as numpy.polyfit
    degree: int
    grey_anchors: np.ndarray
    comp_anchors: np.ndarray
    residuals: np.ndarray
    rmse: float
    r_squared: float
    covariance: Optional[np.ndarray] = None
    notes: List[str] = field(default_factory=list)

    def apply(self, grey) -> np.ndarray:
        return np.polyval(self.coefficients, np.asarray(grey, dtype=float))

    def sigma_from_calibration(self, grey) -> np.ndarray:
        """1-sigma prediction uncertainty from the calibration fit alone."""
        g = np.asarray(grey, dtype=float)
        if self.covariance is None:
            return np.full(g.shape, self.rmse)
        powers = np.vstack([g ** (self.degree - i) for i in range(self.degree + 1)])
        var = np.einsum("ij,ik,jk->k", self.covariance, powers, powers)
        return np.sqrt(np.maximum(var, 0.0) + self.rmse ** 2)

    def describe(self) -> str:
        terms = []
        for i, c in enumerate(self.coefficients):
            p = self.degree - i
            terms.append(f"{c:+.6g}" + ("" if p == 0 else f"*g" if p == 1 else f"*g^{p}"))
        return (f"C = {' '.join(terms)}  (degree {self.degree}, {self.grey_anchors.size} anchors, "
                f"RMSE = {self.rmse:.4g}, R2 = {self.r_squared:.4f})")


def calibrate(grey_anchors, comp_anchors, degree: int = 1) -> GreyscaleCalibration:
    """Fit grey value -> composition on microprobe anchor points.

    Source of the approach: Petrone et al. (2016) Nature Communications 7:12946,
    Methods; Morgan et al. (2004) EPSL 222:933-946.
    """
    g = np.asarray(grey_anchors, dtype=float)
    c = np.asarray(comp_anchors, dtype=float)
    if g.size != c.size:
        raise ValueError("grey and composition anchors must have the same length")
    if g.size < degree + 1:
        raise ValueError(f"a degree-{degree} calibration needs at least {degree+1} anchor points")
    notes: List[str] = []
    cov = None
    if g.size > degree + 1:
        coeffs, cov = np.polyfit(g, c, degree, cov=True)
    else:
        coeffs = np.polyfit(g, c, degree)
        notes.append("exactly determined fit: no degrees of freedom, so the calibration "
                     "uncertainty cannot be estimated from the anchors")
    pred = np.polyval(coeffs, g)
    resid = c - pred
    ss_res = float(np.sum(resid ** 2))
    ss_tot = float(np.sum((c - c.mean()) ** 2))
    rmse = float(np.sqrt(ss_res / g.size))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    if degree == 1 and g.size >= 4 and r2 < 0.95:
        notes.append(f"the linear calibration explains only R2 = {r2:.3f} of the anchor "
                     "variance. Check for a non-linear BSE response or for a second element "
                     "(e.g. Ca) varying along the traverse")
    return GreyscaleCalibration(coeffs, degree, g, c, resid, rmse, r2, cov, notes)


def apply_calibration(x_um, grey, calibration: GreyscaleCalibration,
                      grey_scatter=None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert a grey profile to composition with a propagated uncertainty.

    ``grey_scatter`` is the standard deviation of the grey values across the
    averaged raster lines at each point (NIDIS exports this as the standard
    error of the mean of 200-600 lines).
    """
    x = np.asarray(x_um, dtype=float)
    g = np.asarray(grey, dtype=float)
    C = calibration.apply(g)
    sigma = calibration.sigma_from_calibration(g)
    if grey_scatter is not None:
        slope = np.polyval(np.polyder(calibration.coefficients), g)
        sigma = np.sqrt(sigma ** 2 + (np.abs(slope) * np.asarray(grey_scatter, dtype=float)) ** 2)
    return x, C, sigma


def anchors_from_microprobe(x_grey, grey, x_probe, comp_probe,
                            window_um: float = 1.0) -> Tuple[np.ndarray, np.ndarray]:
    """Average the grey profile around each microprobe spot to build anchor pairs.

    ``window_um`` should be comparable to the microprobe interaction volume so
    the two measurements sample the same material.
    """
    x_grey = np.asarray(x_grey, dtype=float)
    grey = np.asarray(grey, dtype=float)
    out_g, out_c = [], []
    for xp, cp in zip(np.asarray(x_probe, dtype=float), np.asarray(comp_probe, dtype=float)):
        sel = np.abs(x_grey - xp) <= window_um
        if not np.any(sel):
            sel = np.array([int(np.argmin(np.abs(x_grey - xp)))])
        out_g.append(float(np.mean(grey[sel])))
        out_c.append(float(cp))
    return np.array(out_g), np.array(out_c)
