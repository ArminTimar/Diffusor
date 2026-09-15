"""Closed-form solutions of the linear diffusion equation (constant D).

All functions take the product ``Dt`` (m^2 or um^2 -- any consistent unit with
``x``) so that non-isothermal problems can be handled by substituting the
effective integral from :mod:`diffusor.solvers.history`.

Equation numbers refer to Crank, J. (1975) *The Mathematics of Diffusion*,
2nd ed., Oxford University Press (citekey ``crank1975``).  Where a formula is
also given in Costa et al. (2008) RiMG 69 it is cited as a secondary source.

Solutions implemented
---------------------
step_infinite       infinite medium, initial step at x0        Crank section 2.3.1, eq. 2.14
                    (C = C_left for x < x0, C_right for x > x0)
semi_infinite_fixed_surface   x >= 0, C(0,t) = C_s, C(x,0) = C_0   Crank section 2.4, eq. 2.45
band_infinite       band of half-width h at C_in inside C_out   Crank section 2.3.2, eq. 2.15
plane_sheet         -l < x < l, C(x,0)=C_0, C(+-l,t)=C_1        Crank section 4.3.2, eq. 4.17
cylinder            0 < r < a, C(r,0)=C_0, C(a,t)=C_1          Crank section 5.3.1, eq. 5.22
sphere              0 < r < a, C(r,0)=C_0, C(a,t)=C_1          Crank section 6.3.1, eq. 6.18
fraction_*          fractional uptake M_t/M_inf                Crank eqs 4.18, 5.23, 6.20
"""
from __future__ import annotations

import numpy as np
from scipy.special import erf, erfc, j0, j1, jn_zeros

_SQRT_EPS = 1e-300


def _L(Dt):
    return 2.0 * np.sqrt(max(float(Dt), _SQRT_EPS))


# ---------------------------------------------------------------------------
def step_infinite(x, x0: float, C_left: float, C_right: float, Dt: float):
    """Infinite medium with an initial step at x0.

    C(x,t) = C_left + (C_right - C_left) * 1/2 * erfc((x0 - x)/(2 sqrt(Dt)))

    Source: Crank (1975) eq. 2.14 written for a step between two plateaus
    (the classic "diffusion couple"); Costa et al. (2008) eq. 9.  This is the
    form fitted by NIDIS (Petrone et al. 2016, Methods) with free x0 and 2sqrt(Dt).
    """
    x = np.asarray(x, dtype=float)
    return C_left + (C_right - C_left) * 0.5 * erfc((x0 - x) / _L(Dt))


def semi_infinite_fixed_surface(x, C_s: float, C_0: float, Dt: float):
    """Semi-infinite medium x >= 0, surface held at C_s, initial C_0.

    (C - C_0)/(C_s - C_0) = erfc(x/(2 sqrt(Dt)))      Crank (1975) eq. 2.45
    """
    x = np.asarray(x, dtype=float)
    return C_0 + (C_s - C_0) * erfc(x / _L(Dt))


def band_infinite(x, x_center: float, half_width: float, C_in: float, C_out: float, Dt: float):
    """Band |x - x_center| < h initially at C_in in an infinite medium at C_out.

    (C - C_out)/(C_in - C_out) = 1/2 [ erf((h - X)/(2 sqrt(Dt))) + erf((h + X)/(2 sqrt(Dt))) ],
    X = x - x_center.                                   Crank (1975) eq. 2.15
    """
    X = np.asarray(x, dtype=float) - x_center
    L = _L(Dt)
    f = 0.5 * (erf((half_width - X) / L) + erf((half_width + X) / L))
    return C_out + (C_in - C_out) * f


def plane_sheet(x, half_thickness: float, C_0: float, C_1: float, Dt: float, n_terms: int = 200):
    """Plane sheet -l < x < l, uniform initial C_0, surfaces held at C_1.

    (C - C_0)/(C_1 - C_0) = 1 - (4/pi) sum_{n=0}^inf (-1)^n/(2n+1)
                            exp(-D(2n+1)^2 pi^2 t / (4 l^2)) cos((2n+1) pi x / (2 l))
    Source: Crank (1975) eq. 4.17.
    """
    x = np.asarray(x, dtype=float)
    l = float(half_thickness)
    n = np.arange(n_terms)[:, None]
    k = (2 * n + 1)
    terms = ((-1.0) ** n / k) * np.exp(-Dt * k ** 2 * np.pi ** 2 / (4 * l ** 2)) * np.cos(k * np.pi * x[None, :] / (2 * l))
    f = 1.0 - (4.0 / np.pi) * terms.sum(axis=0)
    return C_0 + (C_1 - C_0) * f


def cylinder(r, radius: float, C_0: float, C_1: float, Dt: float, n_terms: int = 100):
    """Infinite cylinder 0 <= r <= a, uniform initial C_0, surface held at C_1.

    (C - C_0)/(C_1 - C_0) = 1 - (2/a) sum_n exp(-D alpha_n^2 t) J0(r alpha_n) / (alpha_n J1(a alpha_n)),
    alpha_n the positive roots of J0(a alpha) = 0.    Source: Crank (1975) eq. 5.22.
    """
    r = np.asarray(r, dtype=float)
    a = float(radius)
    alphas = jn_zeros(0, n_terms) / a
    terms = (np.exp(-Dt * alphas[:, None] ** 2) * j0(r[None, :] * alphas[:, None])
             / (alphas[:, None] * j1(a * alphas[:, None])))
    f = 1.0 - (2.0 / a) * terms.sum(axis=0)
    return C_0 + (C_1 - C_0) * f


def sphere(r, radius: float, C_0: float, C_1: float, Dt: float, n_terms: int = 200):
    """Sphere 0 <= r <= a, uniform initial C_0, surface held at C_1.

    (C - C_0)/(C_1 - C_0) = 1 + (2a/(pi r)) sum_{n=1}^inf ((-1)^n / n) sin(n pi r / a) exp(-D n^2 pi^2 t / a^2)
    Source: Crank (1975) eq. 6.18.  At r = 0 the limit sin(n pi r/a)/r -> n pi/a is used.
    """
    r = np.asarray(r, dtype=float)
    a = float(radius)
    n = np.arange(1, n_terms + 1)[:, None]
    expo = np.exp(-Dt * n ** 2 * np.pi ** 2 / a ** 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(r[None, :] > 0, np.sin(n * np.pi * r[None, :] / a) / np.where(r[None, :] > 0, r[None, :], 1.0),
                         n * np.pi / a)
    f = 1.0 + (2.0 * a / np.pi) * (((-1.0) ** n / n) * ratio * expo).sum(axis=0)
    return C_0 + (C_1 - C_0) * f


# ---- fractional uptake / loss ------------------------------------------------
def fraction_plane_sheet(half_thickness: float, Dt: float, n_terms: int = 200) -> float:
    """M_t/M_inf for a plane sheet, Crank (1975) eq. 4.18."""
    l = float(half_thickness)
    n = np.arange(n_terms)
    k = 2 * n + 1
    return float(1.0 - (8.0 / np.pi ** 2) * np.sum(np.exp(-Dt * k ** 2 * np.pi ** 2 / (4 * l ** 2)) / k ** 2))


def fraction_cylinder(radius: float, Dt: float, n_terms: int = 100) -> float:
    """M_t/M_inf for a cylinder, Crank (1975) eq. 5.23."""
    a = float(radius)
    alphas = jn_zeros(0, n_terms) / a
    return float(1.0 - (4.0 / a ** 2) * np.sum(np.exp(-Dt * alphas ** 2) / alphas ** 2))


def fraction_sphere(radius: float, Dt: float, n_terms: int = 200) -> float:
    """M_t/M_inf for a sphere, Crank (1975) eq. 6.20."""
    a = float(radius)
    n = np.arange(1, n_terms + 1)
    return float(1.0 - (6.0 / np.pi ** 2) * np.sum(np.exp(-Dt * n ** 2 * np.pi ** 2 / a ** 2) / n ** 2))


ANALYTICAL_MODELS = {
    "step_infinite": step_infinite,
    "semi_infinite_fixed_surface": semi_infinite_fixed_surface,
    "band_infinite": band_infinite,
    "plane_sheet": plane_sheet,
    "cylinder": cylinder,
    "sphere": sphere,
}
