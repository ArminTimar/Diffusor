"""Finite-difference solver for the 1-D diffusion equation.

Governing equation (Crank 1975: eq. 1.5 for plane flow with variable D, eq. 5.1
for a cylinder, and the radial reduction of eq. 1.8 for a sphere -- eq. 6.1 is
the constant-D sphere form -- written in a compact form with geometry index m):

    dC/dt = (1/x^m) d/dx [ x^m ( D(C,x,T) dC/dx  -  theta D C dXAn/dx ) ]

with m = 0 (plane), 1 (cylinder, x = r), 2 (sphere, x = r).  The second flux
term is the activity (non-ideality) term of Costa et al. (2003) eq. 7 with
theta = A_i/(R T) (general form: Dohmen, Faak & Blundy 2017, RiMG 83, main-text
eqs 6-7; their Electronic Appendix eqs A7-A8 were not available to check this
implementation against);
it is only active for plagioclase trace elements (``an_profile`` given).

Discretisation: conservative finite volumes on a uniform grid with D at
half-nodes, ``D_{i+1/2} = (D_i + D_{i+1})/2`` (an arithmetic face mean, the usual
conservative finite-volume choice; Dohmen et al. 2017 describe their scheme only
in the Electronic Appendix, which was not available to check, so this is
Diffusor's own implementation), integrated in time with the theta-scheme

    (I - theta_t dt L^{j}) C^{j+1} = (I + (1 - theta_t) dt L^{j}) C^{j}

theta_t = 1/2 gives Crank-Nicolson (Crank 1975 section 8.5, eqs 8.16-8.17, stable
for all r; Dohmen et al. 2017 also use it, but their appendix equations could not
be checked), theta_t = 0 the explicit scheme (Crank eq. 8.12) and
theta_t = 1 fully implicit.  The operator L is lagged (evaluated at C^j); for
composition-dependent D the time step is kept small enough that this is
accurate (default Courant-like number dt D/dx^2 <= 0.5, the explicit stability
limit r = dT/dX^2 <= 1/2, which Crank states in words in sections 8.4.1 and 8.11
without an equation number).

Boundary rows use exact half-cell volumes and zero total external flux for
closed boundaries. Dirichlet rows are replaced by the prescribed boundary value.
Radial cell volumes are integrals of r**m; mass diagnostics use the same weights.
This retains the 2(m+1)D/dx**2 centre coefficient and conserves mass including
activity-driven fluxes at boundary-adjacent faces.
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


def _cell_volumes(x, m):
    """Exact control volumes with node-centred half cells at both ends.

    Common geometric factors (2*pi for cylinders, 4*pi for spheres) cancel.
    The same weights are used by the flux operator and mass diagnostics.
    """
    edges = np.r_[x[0], .5*(x[:-1]+x[1:]), x[-1]]
    return np.diff(edges**(m+1))/(m+1)


def _operator(x, D_nodes, m, an_profile=None, theta_act=0.0):
    """Conservative face-flux divergence, including closed boundary cells."""
    n = x.size
    dx = x[1]-x[0]
    volumes = _cell_volumes(x, m)
    face_area = (.5*(x[:-1]+x[1:]))**m
    conductance = .5*(D_nodes[:-1]+D_nodes[1:])*face_area/dx
    # F = D*dC/dx - theta*D*C*dAn/dx, with centred face concentration.
    drift = np.zeros(n-1) if an_profile is None else .5*theta_act*np.diff(an_profile)
    left = conductance*(1+drift)
    right = conductance*(1-drift)
    lower, diag, upper = np.zeros(n), np.zeros(n), np.zeros(n)
    upper[:-1] = right/volumes[:-1]
    diag[:-1] -= left/volumes[:-1]
    lower[1:] = left/volumes[1:]
    diag[1:] -= right/volumes[1:]
    return lower, diag, upper


def _apply_bc_rows(lower, diag, upper, x, D_nodes, m, bc_left, bc_right):
    # A closed boundary has zero *total* external flux, including activity.
    # Interior-face contributions already appear in the half-cell rows.
    for i, bc in ((0, bc_left), (-1, bc_right)):
        if bc.kind == "dirichlet":
            lower[i] = diag[i] = upper[i] = 0.0
    return lower, diag, upper


def _trapz(y, x):
    """numpy >= 2.0 renamed trapz to trapezoid; support both."""
    f = getattr(np, "trapezoid", None) or np.trapz
    return f(y, x)


def _mass(x, C, m):
    return float(np.dot(C, _cell_volumes(x, m)))


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
    an_profile, theta_activity : plagioclase activity term, theta = A/RT; a callable
        theta_activity(T_K) is evaluated at the temperature of every step
    snapshot_times : store the profile at these times (seconds) exactly
    progress : optional callback receiving fraction done; return True to abort
    """
    x = np.asarray(x, dtype=float)
    C = np.asarray(C0, dtype=float).copy()
    if (x.ndim != 1 or x.size < 3 or C.shape != x.shape
            or not np.all(np.isfinite(x)) or not np.all(np.isfinite(C))
            or np.any(np.diff(x) <= 0)):
        raise ValueError("x and C0 must be finite paired arrays on an increasing grid with at least three nodes")
    if m not in (0, 1, 2) or (m and x[0] < 0):
        raise ValueError("geometry index must be 0, 1 or 2; radial coordinates must be nonnegative")
    if not np.isfinite(t_total) or t_total < 0 or not 0 <= theta_time <= 1 or courant <= 0:
        raise ValueError("invalid duration, theta_time or courant")
    theta_of_T = theta_activity if callable(theta_activity) else (lambda T, th=theta_activity: th)
    T_first = float(history.T(0.0)) if history is not None else float(T_K or 1.0)
    theta0 = float(theta_of_T(T_first))
    if theta0 and an_profile is None:
        raise ValueError("activity-driven transport requires an anorthite profile")
    if an_profile is not None:
        an_profile = np.asarray(an_profile, dtype=float)
        if an_profile.shape != x.shape or not np.all(np.isfinite(an_profile)):
            raise ValueError("anorthite profile must be finite and match the grid")
        if np.max(np.abs(theta0*np.diff(an_profile))) >= 2:
            raise ValueError("activity gradient is unresolved; refine the spatial grid")
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
        result.mass_final = result.mass_initial
        result.C_final = C
        return result

    # --- time-step plan ----------------------------------------------------
    # The explicit scheme is stable only for dt <= dx^2/(2 D) (Crank 1975,
    # r <= 1/2, stated in words in section 8.4.1), but Crank-Nicolson (theta >= 1/2) is unconditionally stable, so
    # for the implicit schemes dt is limited by *accuracy*, not stability.
    # A sharp initial step needs small steps while the front is unresolved,
    # after which the profile is smooth and large steps are accurate. Diffusor
    # therefore starts at the explicit limit and grows dt geometrically up to
    # a cap set by ``min_steps``, which keeps a 100 kyr run to a few hundred
    # steps instead of millions.
    Dmax0 = float(np.max(D_func(C, float(history.T(0.0)))))
    if not np.isfinite(Dmax0) or Dmax0 <= 0:
        raise ValueError("D must be positive")
    dt_explicit = courant * dx ** 2 / Dmax0
    if theta_time < 0.5:
        dt_start = dt_cap = dt_explicit
        growth = 1.0
    else:
        dt_start = min(dt_explicit, t_total)
        dt_cap = max(dt_start, t_total / max(min_steps, 2))
        growth = float(dt_growth)
    if t_total / dt_cap > max_steps and theta_time < 0.5:
        raise ValueError("Explicit stability requires more than max_steps; use an implicit method or increase max_steps")
    if t_total / dt_cap > max_steps:
        result.warnings.append(
            f"{t_total/dt_cap:.3g} steps would be needed at the stability limit "
            f"(max_steps = {max_steps}). The time step was enlarged. Check the "
            "result against a finer grid.")
        dt_cap = t_total / max_steps
    dt_nominal = dt_cap

    t = 0.0
    step = 0
    dt_step = dt_start
    ab = np.zeros((3, n))
    next_snap_idx = 0
    while t < t_total * (1 - 1e-12):
        if step >= max_steps:
            raise RuntimeError("max_steps reached before the requested duration; no partial solution returned")
        T_now = float(history.T(t))
        D_nodes = np.asarray(D_func(C, T_now), dtype=float)
        if D_nodes.size == 1:
            D_nodes = np.full(n, float(D_nodes.item()))
        if D_nodes.shape != (n,) or not np.all(np.isfinite(D_nodes)) or np.any(D_nodes < 0):
            raise ValueError("D must be finite, nonnegative and scalar or match the grid")
        Dmax = float(D_nodes.max())
        if theta_time < 0.5:
            dt = courant * dx ** 2 / Dmax if Dmax > 0 else dt_nominal
        else:
            dt = min(dt_step * growth, dt_cap)
        dt_step = dt
        dt = min(dt, t_total - t)
        if next_snap_idx < len(snaps):
            dt = min(dt, snaps[next_snap_idx] - t)
        lower, diag, upper = _operator(x, D_nodes, m, an_profile, float(theta_of_T(T_now)))
        lower, diag, upper = _apply_bc_rows(lower, diag, upper, x, D_nodes, m, bc_left, bc_right)
        if theta_time < 0.5:
            # The radial centre and drift can have a larger exit rate than
            # the planar dx²/(2D) estimate. Keep the explicit RHS nonnegative.
            exit_rate = float(np.max(-diag))
            if exit_rate > 0:
                dt = min(dt, .95/((1-theta_time)*exit_rate))
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
