"""Initial conditions.

The initial condition is the single largest source of systematic error in
diffusion chronometry: a profile produced partly by growth and partly by
diffusion will give a spuriously long time if the whole of it is attributed to
diffusion (Costa et al. 2008 section "Initial conditions"; Shea et al. 2015).
Diffusor therefore keeps the initial condition explicit and plots it alongside
the fit so the assumption is always visible.

Available forms
---------------
``step``            two plateaus meeting at x0 (the classic diffusion couple)
``multi_step``      several plateaus (growth zoning with sharp boundaries)
``plateau_rim``     interior plateau plus a rim of different composition
``table``           an arbitrary user-supplied profile
``equilibrium_plag``  the quasi-steady-state trace-element profile implied by a
                    frozen anorthite profile (Dohmen et al. 2017 App. eq. A13)
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, List, Optional, Sequence, Tuple

import numpy as np


def step(x, x0: float, C_left: float, C_right: float, smooth: float = 0.0):
    """Sharp step at x0, optionally smoothed over ``smooth`` (same units as x).

    A little smoothing is numerically helpful (the DMG Short Course 2025 script
    ``Diff_Model_Sr_in_Plag_implicit.m`` applies 10 explicit smoothing sweeps
    for exactly this reason), but it must be small compared with the diffusion
    length or it will bias the fitted time.
    """
    x = np.asarray(x, dtype=float)
    if smooth and smooth > 0:
        from scipy.special import erf
        return C_left + (C_right - C_left) * 0.5 * (1.0 + erf((x - x0) / (smooth * np.sqrt(2.0))))
    # A node sitting exactly on the interface gets the cell average of the two
    # plateaus. Without this the discretised step is offset by half a cell and
    # the numerical solution is shifted by dx/2 relative to the Crank solution.
    return np.where(x < x0, float(C_left),
                    np.where(x > x0, float(C_right),
                             0.5 * (float(C_left) + float(C_right)))).astype(float)


def multi_step(x, edges: Sequence[float], values: Sequence[float], smooth: float = 0.0):
    """Piecewise-constant profile: ``values`` has one more entry than ``edges``."""
    x = np.asarray(x, dtype=float)
    if len(values) != len(edges) + 1:
        raise ValueError("multi_step needs len(values) == len(edges) + 1")
    C = np.full_like(x, float(values[0]))
    for e, v_next in zip(edges, values[1:]):
        C = C + (float(v_next) - C) * (
            0.5 * (1.0 + np.tanh((x - e) / max(smooth, 1e-12))) if smooth > 0
            else (x >= e).astype(float))
    return C


def plateau_rim(x, rim_start: float, C_core: float, C_rim: float, smooth: float = 0.0):
    """Interior plateau with a rim beyond ``rim_start``."""
    return step(x, rim_start, C_core, C_rim, smooth)


def from_table(x, x_table, C_table):
    """Interpolate a user-supplied initial profile onto the model grid."""
    return np.interp(np.asarray(x, dtype=float), np.asarray(x_table, dtype=float),
                     np.asarray(C_table, dtype=float))


@dataclass
class InitialCondition:
    """A named, re-evaluable initial condition."""
    kind: str = "step"
    params: dict = field(default_factory=dict)
    description: str = ""

    def evaluate(self, x) -> np.ndarray:
        k = self.kind
        p = self.params
        if k == "step":
            return step(x, p["x0"], p["C_left"], p["C_right"], p.get("smooth", 0.0))
        if k == "multi_step":
            return multi_step(x, p["edges"], p["values"], p.get("smooth", 0.0))
        if k == "plateau_rim":
            return plateau_rim(x, p["rim_start"], p["C_core"], p["C_rim"], p.get("smooth", 0.0))
        if k == "table":
            return from_table(x, p["x_table"], p["C_table"])
        if k == "equilibrium_plag":
            from ..coefficients.plagioclase import equilibrium_profile
            X_An = np.interp(np.asarray(x, dtype=float), p["x_an_x"], p["x_an_values"])
            return equilibrium_profile(X_An, p["T_K"], p["species"], p["C_ref"], p.get("X_An_ref"),
                                       p.get("activity_set", "dohmen_blundy2014"))
        raise ValueError(f"unknown initial condition '{k}'")

    def describe(self) -> str:
        if self.description:
            return self.description
        if self.kind == "step":
            p = self.params
            return (f"step at x = {p['x0']:g}: {p['C_left']:g} -> {p['C_right']:g}"
                    + (f", smoothed over {p['smooth']:g}" if p.get("smooth") else ""))
        return self.kind


def guess_step_from_data(x, C, plateau_fraction: float = 0.15) -> InitialCondition:
    """First guess for a step initial condition from the measured profile.

    The plateaus are the medians of the outer ``plateau_fraction`` of the
    points at each end, and the interface is placed where the profile crosses
    the midpoint between them.  The midpoint crossing is used rather than the
    steepest gradient because on a broad, coarsely sampled profile the noisiest
    single interval can easily out-gradient the real step and drag the
    interface to the edge of the traverse.  The search is additionally confined
    to the interior of the profile, outside the plateau windows.

    This is only a starting point: inspect it, and let the fit refine x0 if the
    guess is not obviously right.
    """
    x = np.asarray(x, dtype=float)
    C = np.asarray(C, dtype=float)
    order = np.argsort(x)
    x, C = x[order], C[order]
    n = max(2, int(plateau_fraction * x.size))
    C_left = float(np.median(C[:n]))
    C_right = float(np.median(C[-n:]))

    lo, hi = n - 1, x.size - n
    if hi <= lo:
        lo, hi = 0, x.size - 1
    mid = 0.5 * (C_left + C_right)

    x0 = None
    if abs(C_right - C_left) > 1e-12:
        # first crossing of the midpoint inside the interior window
        seg = C[lo:hi + 1] - mid
        sign_change = np.where(np.sign(seg[:-1]) * np.sign(seg[1:]) <= 0)[0]
        if sign_change.size:
            i = lo + int(sign_change[sign_change.size // 2])
            c0, c1 = C[i], C[min(i + 1, x.size - 1)]
            if abs(c1 - c0) > 1e-12:
                frac = (mid - c0) / (c1 - c0)
                x0 = float(x[i] + frac * (x[min(i + 1, x.size - 1)] - x[i]))
            else:
                x0 = float(x[i])
    if x0 is None:
        # fall back to the steepest gradient of a lightly smoothed profile,
        # still restricted to the interior
        k = max(3, x.size // 20) | 1
        kern = np.ones(k) / k
        Cs = np.convolve(np.pad(C, k // 2, mode="edge"), kern, mode="valid")
        g = np.abs(np.gradient(Cs, x))
        g[:lo] = -np.inf
        g[hi + 1:] = -np.inf
        x0 = float(x[int(np.argmax(g))])

    return InitialCondition(
        "step", {"x0": x0, "C_left": C_left, "C_right": C_right},
        f"step at the midpoint crossing (x0 = {x0:.2f}), plateaus from the median of the "
        f"outer {n} points at each end ({C_left:.4g} and {C_right:.4g})")
