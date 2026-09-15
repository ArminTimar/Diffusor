"""Finite-difference solver for the 1-D diffusion equation.

Governing equation (Crank 1975 eq. 1.7 for planar; section 5.1 eq. 5.4 cylinder;
section 6.1 eq. 6.3 sphere; compact form with geometry index m):

    dC/dt = (1/x^m) d/dx [ x^m ( D(C,x,T) dC/dx  -  theta D C dXAn/dx ) ]

with m = 0 (plane), 1 (cylinder, x = r), 2 (sphere, x = r).  The second flux
term is the activity (non-ideality) term of Costa et al. (2003) eq. 7 with
theta = A_i/(R T) (Dohmen, Faak & Blundy 2017, RiMG 83, Appendix eqs A7-A8);
it is only active for plagioclase trace elements (``an_profile`` given).

Discretisation: conservative finite volumes on a uniform grid with D at
half-nodes, ``D_{i+1/2} = (D_i + D_{i+1})/2`` (Dohmen et al. 2017 App. eqs
A17-A19), integrated in time with the theta-scheme

    (I - theta_t dt L^{j}) C^{j+1} = (I + (1 - theta_t) dt L^{j}) C^{j}

theta_t = 1/2 gives Crank-Nicolson (Crank 1975 section 8.4, eq. 8.35; Dohmen et al.
2017 App. eq. A21), theta_t = 0 the explicit scheme (Crank eq. 8.31) and
theta_t = 1 fully implicit.  The operator L is lagged (evaluated at C^j); for
composition-dependent D the time step is kept small enough that this is
accurate (default Courant-like number dt D/dx^2 <= 0.5, the explicit stability
limit, Crank eq. 8.33).

Boundary rows: Dirichlet -> identity row with the boundary value at t^{j+1};
zero flux at a node with x > 0 -> mirror ghost node; symmetry at x = 0 for
m > 0 -> limit  dC/dt = (m+1) D d2C/dx2  (L'Hopital, Crank section 8.5 eq. 8.45).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence

import numpy as np
from scipy.linalg import solve_banded

from .boundary import BoundaryCondition, zero_flux
from .history import ThermalHistory


@dataclass
class NumericalResult:
    x: np.ndarray
    C_final: np.ndarray
    t_final: float
    snapshots: Dict[float, np.ndarray] = field(default_factory=dict)
    n_steps: int = 0
    dt_used: float = 0.0
    mass_initial: float = 0.0
    mass_final: float = 0.0
    warnings: List[str] = field(default_factory=list)


def _operator(x, D_nodes, m, an_profile=None, theta_act=0.0):
    """Tridiagonal operator L (lower, diag, upper) for interior nodes."""
    n = x.size
    dx = x[1] - x[0]
    Dh = 0.5 * (D_nodes[:-1] + D_nodes[1:])            # D at i+1/2, length n-1
    xh = 0.5 * (x[:-1] + x[1:])                         # x at i+1/2
    lower = np.zeros(n)
    diag = np.zeros(n)
    upper = np.zeros(n)
    i = np.arange(1, n - 1)
    if m == 0:
        gp = np.ones(n - 1)
        gm = np.ones(n - 1)
    else:
        gp = np.empty(n - 1)
        gm = np.empty(n - 1)
        # factors (x_{i+1/2}/x_i)^m and (x_{i-1/2}/x_i)^m applied per interior node
        gp[:] = (xh / np.where(x[:-1] > 0, x[:-1], 1.0)) ** m       # for node i using xh[i]
        gm[:] = (xh / np.where(x[1:] > 0, x[1:], 1.0)) ** m         # for node i using xh[i-1]
    # coefficients for interior nodes
    cp = Dh[i] * (gp[i] if m else 1.0) / dx ** 2          # multiplies C_{i+1}
    cm = Dh[i - 1] * (gm[i - 1] if m else 1.0) / dx ** 2  # multiplies C_{i-1}
    upper[i] = cp
    lower[i] = cm
    diag[i] = -(cp + cm)
    if an_profile is not None and theta_act != 0.0:
        dAn = np.diff(an_profile)                        # An_{i+1} - An_i, length n-1
        # Dohmen et al. 2017 App. eq. A20 (explicit form) coefficients of the activity term
        ap = Dh[i] * (gp[i] if m else 1.0) * theta_act * 0.5 * dAn[i] / dx ** 2
        am = Dh[i - 1] * (gm[i - 1] if m else 1.0) * theta_act * 0.5 * dAn[i - 1] / dx ** 2
        upper[i] -= ap
        diag[i] -= (ap - am)
        lower[i] += am
    return lower, diag, upper


def _apply_bc_rows(lower, diag, upper, x, D_nodes, m, bc_left, bc_right):
    n = x.size
    dx = x[1] - x[0]
    # left node 0
    if bc_left.kind == "dirichlet":
        lower[0] = diag[0] = upper[0] = 0.0
    else:
        Dh = 0.5 * (D_nodes[0] + D_nodes[1])
        if m > 0 and x[0] <= 0.0:
            c = 2.0 * (m + 1) * Dh / dx ** 2          # symmetry at r = 0
        else:
            c = 2.0 * Dh / dx ** 2                     # mirror ghost node
        upper[0] = c
        diag[0] = -c
        lower[0] = 0.0
    # right node n-1
    if bc_right.kind == "dirichlet":
        lower[-1] = diag[-1] = upper[-1] = 0.0
    else:
        Dh = 0.5 * (D_nodes[-2] + D_nodes[-1])
        c = 2.0 * Dh / dx ** 2
        if m > 0:
            xh = 0.5 * (x[-2] + x[-1])
            c *= (xh / x[-1]) ** m
        lower[-1] = c
        diag[-1] = -c
        upper[-1] = 0.0
    return lower, diag, upper


def _trapz(y, x):
    """numpy >= 2.0 renamed trapz to trapezoid; support both."""
    f = getattr(np, "trapezoid", None) or np.trapz
    return f(y, x)


def _mass(x, C, m):
    w = x ** m if m else np.ones_like(x)
    return float(_trapz(C * w, x))


def solve_1d(x: np.ndarray, C0: np.ndarray, D_func: Callable[[np.ndarray, float], np.ndarray],
             t_total: float, *, m: int = 0,
             bc_left: BoundaryCondition = None, bc_right: BoundaryCondition = None,
             history: Optional[ThermalHistory] = None, T_K: Optional[float] = None,
             theta_time: float = 0.5, courant: float = 0.5, max_steps: int = 200_000,
             min_steps: int = 400, dt_growth: float = 1.07,
             an_profile: Optional[np.ndarray] = None, theta_activity: float = 0.0,
             snapshot_times: Optional[Sequence[float]] = None,
             progress: Optional[Callable[[float], bool]] = None) -> NumericalResult:
    """Integrate the diffusion equation from C0 to time t_total.

    Parameters
    ----------
    x : uniform grid (any length unit; D must use the same unit squared)
    C0 : initial profile on x
    D_func : callable (C_array, T_K) -> D array on nodes (length x)
    t_total : seconds
    m : 0 plane, 1 cylinder, 2 sphere (x is then the radius, x[0] may be 0)
    bc_left, bc_right : BoundaryCondition (default zero flux both sides)
    history : ThermalHistory giving T(t); if None, ``T_K`` is used (isothermal)
    theta_time : 0.5 Crank-Nicolson, 1 implicit, 0 explicit
    courant : first time step as a fraction of the explicit stability limit
    min_steps : lower bound on the number of steps (sets the dt cap for implicit schemes)
    dt_growth : geometric growth factor applied to dt each step (implicit schemes only)
    an_profile, theta_activity : plagioclase activity term (theta = A/RT)
    snapshot_times : store the profile at these times (seconds) exactly
    progress : optional callback receiving fraction done; return True to abort
    """
    x = np.asarray(x, dtype=float)
    C = np.asarray(C0, dtype=float).copy()
    n = x.size
    dx = x[1] - x[0]
    if not np.allclose(np.diff(x), dx, rtol=1e-6):
        raise ValueError("numerical solver requires a uniform grid")
    bc_left = bc_left or zero_flux()
    bc_right = bc_right or zero_flux()
    if history is None:
        if T_K is None:
            raise ValueError("either history or T_K must be given")
        history = ThermalHistory.isothermal(T_K, t_total)
    snaps = sorted(set(float(s) for s in (snapshot_times or []) if 0.0 < s <= t_total))
    result = NumericalResult(x=x, C_final=C, t_final=0.0, mass_initial=_mass(x, C, m))
    if t_total <= 0:
        result.C_final = C
        return result

    # --- time-step plan ----------------------------------------------------
    # The explicit scheme is stable only for dt <= dx^2/(2 D) (Crank 1975 eq.
    # 8.33), but Crank-Nicolson (theta >= 1/2) is unconditionally stable, so
    # for the implicit schemes dt is limited by *accuracy*, not stability.
    # A sharp initial step needs small steps while the front is unresolved,
    # after which the profile is smooth and large steps are accurate. Diffusor
    # therefore starts at the explicit limit and grows dt geometrically up to
    # a cap set by ``min_steps``, which keeps a 100 kyr run to a few hundred
    # steps instead of millions.
    Dmax0 = float(np.max(D_func(C, float(history.T(0.0)))))
    if Dmax0 <= 0:
        raise ValueError("D must be positive")
    dt_explicit = courant * dx ** 2 / Dmax0
    if theta_time < 0.5:
        dt_start = dt_cap = dt_explicit
        growth = 1.0
    else:
        dt_start = min(dt_explicit, t_total)
        dt_cap = max(dt_start, t_total / max(min_steps, 2))
        growth = float(dt_growth)
    if t_total / dt_cap > max_steps:
        result.warnings.append(
            f"{t_total/dt_cap:.3g} steps would be needed at the stability limit "
            f"(max_steps = {max_steps}); the time step was enlarged, so check the "
            "result against a finer grid.")
        dt_cap = t_total / max_steps
    dt_nominal = dt_cap

    t = 0.0
    step = 0
    dt_step = dt_start
    ab = np.zeros((3, n))
    next_snap_idx = 0
    while t < t_total * (1 - 1e-12):
        T_now = float(history.T(t))
        D_nodes = np.asarray(D_func(C, T_now), dtype=float)
        if D_nodes.shape != (n,):
            D_nodes = np.full(n, float(np.ravel(D_nodes)[0]))
        Dmax = float(D_nodes.max())
        if theta_time < 0.5:
            dt = courant * dx ** 2 / Dmax if Dmax > 0 else dt_nominal
        else:
            dt = min(dt_step * growth, dt_cap)
        dt_step = dt
        dt = min(dt, t_total - t)
        if next_snap_idx < len(snaps):
            dt = min(dt, snaps[next_snap_idx] - t)
        lower, diag, upper = _operator(x, D_nodes, m, an_profile, theta_activity)
        lower, diag, upper = _apply_bc_rows(lower, diag, upper, x, D_nodes, m, bc_left, bc_right)
        # explicit part: rhs = C + (1-theta) dt L C
        LC = diag * C
        LC[:-1] += upper[:-1] * C[1:]
        LC[1:] += lower[1:] * C[:-1]
        rhs = C + (1.0 - theta_time) * dt * LC
        # implicit part matrix (I - theta dt L) in banded form
        ab[0, 1:] = -theta_time * dt * upper[:-1]
        ab[1, :] = 1.0 - theta_time * dt * diag
        ab[2, :-1] = -theta_time * dt * lower[1:]
        t_new = t + dt
        if bc_left.kind == "dirichlet":
            rhs[0] = bc_left.value_at(t_new)
            ab[1, 0] = 1.0
            ab[0, 1] = 0.0
        if bc_right.kind == "dirichlet":
            rhs[-1] = bc_right.value_at(t_new)
            ab[1, -1] = 1.0
            ab[2, -2] = 0.0
        if theta_time > 0:
            C = solve_banded((1, 1), ab, rhs)
        else:
            C = rhs
        t = t_new
        step += 1
        if next_snap_idx < len(snaps) and abs(t - snaps[next_snap_idx]) <= 1e-9 * max(t, 1.0):
            result.snapshots[snaps[next_snap_idx]] = C.copy()
            next_snap_idx += 1
        if progress is not None and step % 200 == 0:
            if progress(t / t_total):
                result.warnings.append("aborted by user")
                break
    result.C_final = C
    result.t_final = t
    result.n_steps = step
    result.dt_used = dt_nominal
    result.mass_final = _mass(x, C, m)
    return result


def make_grid(x_min: float, x_max: float, n: int) -> np.ndarray:
    return np.linspace(x_min, x_max, int(n))
