"""Temperature-time histories and the effective diffusion integral.

For a time-dependent diffusion coefficient D(t) (e.g. through T(t)), the
solution of the linear diffusion equation with constant-D form is recovered by
replacing ``D t`` with the *effective* integral

    Dt_eff = integral_0^t D(T(t')) dt'

(Crank 1975, section 7.1 "Time-dependent diffusion coefficients", eqs 7.2-7.3,
substitution ``dT = D(t) dt``, which reduces eq. 7.1 to dC/dT = d2C/dx2 (7.4) and
is exact only when D depends on time alone; used in geospeedometry by Lasaga 1983 and reviewed by Costa
et al. 2008.)  The numerical solver instead evaluates D at every time step, so
both routes give identical answers for constant initial/boundary geometry.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, List, Sequence, Tuple

import numpy as np


@dataclass
class ThermalHistory:
    """Piecewise-linear T(t) path.

    ``times`` in seconds (increasing, starting at 0) and ``temps_K`` in K.
    An isothermal history has a single segment.
    """
    times: np.ndarray
    temps_K: np.ndarray
    label: str = "isothermal"

    @classmethod
    def isothermal(cls, T_K: float, duration_s: float = np.inf) -> "ThermalHistory":
        end = duration_s if np.isfinite(duration_s) else 1.0
        return cls(np.array([0.0, end]), np.array([T_K, T_K]), "isothermal")

    @classmethod
    def linear(cls, T_start_K: float, T_end_K: float, duration_s: float) -> "ThermalHistory":
        return cls(np.array([0.0, duration_s]), np.array([T_start_K, T_end_K]),
                   f"linear {T_start_K:.0f}->{T_end_K:.0f} K")

    @classmethod
    def piecewise(cls, points: Sequence[Tuple[float, float]]) -> "ThermalHistory":
        pts = sorted(points)
        t = np.array([p[0] for p in pts], dtype=float)
        T = np.array([p[1] for p in pts], dtype=float)
        if t[0] != 0.0:
            raise ValueError("history must start at t = 0")
        return cls(t, T, "piecewise")

    @property
    def is_isothermal(self) -> bool:
        return bool(np.allclose(self.temps_K, self.temps_K[0]))

    @property
    def duration(self) -> float:
        return float(self.times[-1])

    def T(self, t):
        """Temperature at time t (held constant beyond the last node)."""
        return np.interp(np.asarray(t, dtype=float), self.times, self.temps_K,
                         left=self.temps_K[0], right=self.temps_K[-1])

    def shifted_to_end(self, t_total: float) -> "ThermalHistory":
        """Return a copy whose time axis is rescaled so the path ends at t_total.

        Used when the user fixes the *shape* of a cooling path but fits its
        duration.
        """
        if self.is_isothermal:
            return ThermalHistory(np.array([0.0, t_total]), self.temps_K[:2].copy(), self.label)
        scale = t_total / self.duration
        return ThermalHistory(self.times * scale, self.temps_K.copy(), self.label)


def effective_Dt(D_of_T: Callable[[np.ndarray], np.ndarray], history: ThermalHistory,
                 t_total: float, n: int = 2001) -> float:
    """integral_0^t D(T(t')) dt' by the trapezoidal rule on a fine time grid.

    ``D_of_T`` must accept an array of temperatures (K) and return D (m^2/s).
    Source: Crank (1975) eqs 7.2-7.3 (T = integral of D(t') dt'); Lasaga (1983).
    """
    if history.is_isothermal:
        return float(D_of_T(np.array([history.temps_K[0]]))[0]) * t_total
    tt = np.linspace(0.0, t_total, n)
    D = D_of_T(history.T(tt))
    _f = getattr(np, "trapezoid", None) or np.trapz
    return float(_f(D, tt))
