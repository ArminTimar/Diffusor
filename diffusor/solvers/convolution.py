"""Analytical spatial-resolution (convolution) correction.

A microbeam measures a concentration averaged over its interaction volume.
For a 1-D traverse this is modelled as the convolution of the true profile
with a Gaussian of standard deviation ``sigma``:

    f(x) = integral C(x + sigma z) phi(z) dz,    phi = standard normal density

which is Ganguly, Bhattacharya & Chakraborty (1988, Am Mineral 73:901-909)
eq. 9 (printed p. 903; their Gaussian standard deviation is called epsilon,
and their eq. 10 is the discrete sum over +-4 epsilon/Delta grid steps used
here).  Bradshaw & Kent (2017, Chem Geol 466:667-677) and the DMG Short Course
2025 Practical 5 treat the same effect.  The model profile is convolved
*before* it is compared with the data, so the fitted time is corrected for beam
smearing.

Resolution warning (Diffusor's own rule of thumb, not a published criterion):
if the apparent diffusion length ``2 sqrt(Dt)`` is smaller than ``factor`` x
sigma (default 3) the retrieved time is dominated by the beam and should be
reported as an upper bound only.  The published criterion of Bradshaw & Kent
(2017) is different in form: a modelled timescale is accurate to about 20 %
when the gradient width (length between the 5th and 95th percentile of the
step, which is 4*erfinv(0.9)*sqrt(Dt) = 2.33 x 2 sqrt(Dt) for an erf step)
exceeds twice the spot size, and at least three measurements fall inside the
gradient.  They resample by averaging over an interval equal to the spot size,
so a Gaussian sigma maps onto their spot size only through an assumption; for a
box of equal variance (width sqrt(12) sigma) their rule corresponds to
2 sqrt(Dt) > 2.98 sigma, close to the default factor of 3, whereas if the spot
size were the Gaussian FWHM (2.355 sigma) it would be 2 sqrt(Dt) > 2.02 sigma.
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
    """Return a warning string if 2*sqrt(Dt) < factor*sigma, else None.

    The default ``factor = 3`` is Diffusor's own choice; see the module
    docstring for how it compares with the Bradshaw & Kent (2017) criterion.
    """
    if sigma is None or sigma <= 0:
        return None
    L = 2.0 * np.sqrt(max(Dt, 0.0))
    if L < factor * sigma:
        return (f"Diffusion length 2*sqrt(Dt) = {L:.3g} is below {factor:g} x beam sigma "
                f"({sigma:.3g}). The timescale is at the analytical resolution limit. Treat it "
                "as an upper bound (the 3-sigma threshold is Diffusor's rule of thumb; compare the "
                "gradient width/spot size > 2 criterion of Bradshaw & Kent 2017).")
    return None
