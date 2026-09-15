"""Analytical spatial-resolution (convolution) correction.

A microbeam measures a concentration averaged over its interaction volume.
For a 1-D traverse this is modelled as the convolution of the true profile
with a Gaussian of standard deviation ``sigma`` (Ganguly et al. 1988, Am
Mineral 73:901-909; Bradshaw & Kent 2017, Chem Geol 466:667-677; DMG Short
Course 2025 Practical 5).  The model profile is convolved *before* it is
compared with the data, so the fitted time is corrected for beam smearing.

Rule of thumb (Bradshaw & Kent 2017): if the apparent diffusion length
``2 sqrt(Dt)`` is smaller than roughly 3 sigma the retrieved time is dominated
by the beam and should be reported as an upper bound only.
"""
from __future__ import annotations

import numpy as np


def gaussian_convolve(x: np.ndarray, C: np.ndarray, sigma: float) -> np.ndarray:
    """Convolve profile C(x) (uniform grid) with a Gaussian of std ``sigma`` (same units as x).

    Edges are padded with the end values (no artificial gradient at the ends).
    """
    x = np.asarray(x, dtype=float)
    C = np.asarray(C, dtype=float)
    if sigma is None or sigma <= 0.0:
        return C.copy()
    dx = float(np.mean(np.diff(x)))
    half = int(np.ceil(4.0 * sigma / dx))
    if half < 1:
        return C.copy()
    k = np.arange(-half, half + 1) * dx
    kern = np.exp(-0.5 * (k / sigma) ** 2)
    kern /= kern.sum()
    padded = np.concatenate([np.full(half, C[0]), C, np.full(half, C[-1])])
    return np.convolve(padded, kern, mode="valid")


def resolution_warning(Dt: float, sigma: float, factor: float = 3.0):
    """Return a warning string if 2*sqrt(Dt) < factor*sigma, else None."""
    if sigma is None or sigma <= 0:
        return None
    L = 2.0 * np.sqrt(max(Dt, 0.0))
    if L < factor * sigma:
        return (f"Diffusion length 2*sqrt(Dt) = {L:.3g} is < {factor:g} x beam sigma ({sigma:.3g}); "
                "the timescale is at the analytical resolution limit and should be treated as an "
                "upper bound (Bradshaw & Kent 2017).")
    return None
