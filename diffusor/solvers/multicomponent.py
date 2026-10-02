"""Finite-volume solver for 1-D multicomponent diffusion with a composition-dependent D matrix.

Governing equation (Chakraborty & Ganguly 1992 eq. 1; Oeser et al. 2026 eq. 3),
in plane (m = 0), cylinder (m = 1) or sphere (m = 2) geometry:

    dC/dt = (1/x^m) d/dx [ x^m D(C) dC/dx ]

C is the vector of the n-1 independent mole fractions and D the (n-1)x(n-1)
Fick matrix, which need not be symmetric and can have negative off-diagonal
terms (uphill diffusion). The discretisation is the one of
:mod:`diffusor.solvers.numerical` applied to a vector: one flux per face with
the face matrix ``(D_i + D_{i+1})/2``, exact radial control volumes, closed
boundaries with zero total flux, and the theta scheme with the matrix lagged
at the start of each step. Each component is conserved to round-off with closed
boundaries, and the dependent component follows from the sum.

A closed-form check for a constant matrix (:func:`step_constant_matrix`)
diagonalises D and lets each eigen-component diffuse on its own (Toor 1964,
as used by Chakraborty & Ganguly 1992, p. 78).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from typing import Callable, List, Optional, Sequence, Union

import numpy as np
from scipy.linalg import solve_banded

from .analytical import step_infinite
from .history import ThermalHistory
from .numerical import _cell_volumes

Fixed = Union[None, Sequence[float], Callable[[float], Sequence[float]]]


@dataclass
class MulticomponentResult:
    x: np.ndarray
    C_final: np.ndarray                      # (n-1, N)
    t_final: float
    n_steps: int = 0
    mass_initial: np.ndarray = field(default_factory=lambda: np.zeros(0))
    mass_final: np.ndarray = field(default_factory=lambda: np.zeros(0))
    min_value: float = 0.0
    warnings: List[str] = field(default_factory=list)


def _fixed(value: Fixed, t: float, n: int) -> Optional[np.ndarray]:
    if value is None:
        return None
    v = np.asarray(value(t) if callable(value) else value, dtype=float)
    if v.shape != (n,) or not np.all(np.isfinite(v)):
        raise ValueError(f"a fixed boundary needs {n} finite values")
    return v


@lru_cache(maxsize=32)
def _pattern(n: int, N: int, fixed: tuple):
    """Index pattern of the operator: for each of four contribution groups the
    (row, column, face, i, j, cell) arrays, with rows of fixed nodes removed."""
    i, j, f = (a.ravel() for a in np.meshgrid(np.arange(n), np.arange(n), np.arange(N - 1), indexing="ij"))
    out = []
    for rows, cols, cell in (((f * n + i), (f + 1) * n + j, f), (f * n + i, f * n + j, f),
                             ((f + 1) * n + i, (f + 1) * n + j, f + 1), ((f + 1) * n + i, f * n + j, f + 1)):
        keep = ~np.isin(rows, np.array(fixed, dtype=int))
        out.append((rows[keep], cols[keep], f[keep], i[keep], j[keep], cell[keep]))
    return out


_SIGNS = (1.0, -1.0, -1.0, 1.0)


def _operator(x, D_nodes, m, fixed_rows):
    """Divergence operator L as groups of COO triplets over the interleaved unknowns
    k*n + i (node k, component i). Within a group no (row, column) pair repeats.
    Rows of fixed nodes are dropped."""
    n, N = D_nodes.shape[0], x.size
    dx = x[1] - x[0]
    vol = _cell_volumes(x, m)
    area = (0.5 * (x[:-1] + x[1:])) ** m
    G = 0.5 * (D_nodes[:, :, :-1] + D_nodes[:, :, 1:]) * (area / dx)       # (n, n, N-1)
    fixed = tuple(int(r) for r in (fixed_rows if fixed_rows is not None else ()))
    return [(rows, cols, sign * G[i, j, f] / vol[cell])
            for sign, (rows, cols, f, i, j, cell) in zip(_SIGNS, _pattern(n, N, fixed))]


def solve_multicomponent(x: np.ndarray, C0: np.ndarray,
                         D_func: Callable[[np.ndarray, float], np.ndarray], t_total: float, *,
                         m: int = 0, left: Fixed = None, right: Fixed = None,
                         history: Optional[ThermalHistory] = None, T_K: Optional[float] = None,
                         theta_time: float = 0.5, courant: float = 0.5, min_steps: int = 400,
                         dt_growth: float = 1.07, max_steps: int = 200_000,
                         progress: Optional[Callable[[float], bool]] = None) -> MulticomponentResult:
    """Integrate the vector diffusion equation from ``C0`` (shape (n-1, N)) to ``t_total``.

    ``D_func(C, T_K)`` returns the matrix at every node, shape (n-1, n-1, N), in
    the square of the length unit of ``x`` per second. ``left`` and ``right``
    are ``None`` (closed, zero flux) or the fixed independent compositions at
    that end (a sequence, or a callable of time).
    """
    x = np.asarray(x, dtype=float)
    C = np.array(C0, dtype=float, copy=True)
    if C.ndim != 2 or x.ndim != 1 or C.shape[1] != x.size or x.size < 3:
        raise ValueError("C0 must have shape (components, nodes) on a grid of at least three nodes")
    if not (np.all(np.isfinite(x)) and np.all(np.isfinite(C))) or np.any(np.diff(x) <= 0):
        raise ValueError("x and C0 must be finite and x increasing")
    dx = x[1] - x[0]
    if not np.allclose(np.diff(x), dx, rtol=1e-6):
        raise ValueError("the solver requires a uniform grid")
    if m not in (0, 1, 2) or (m and x[0] < 0):
        raise ValueError("geometry index must be 0, 1 or 2 with nonnegative radii")
    if not np.isfinite(t_total) or t_total < 0 or not 0.5 <= theta_time <= 1:
        raise ValueError("invalid duration or theta_time (0.5 to 1)")
    n, N = C.shape
    if history is None:
        if T_K is None:
            raise ValueError("either history or T_K must be given")
        history = ThermalHistory.isothermal(T_K, t_total)
    vol = _cell_volumes(x, m)
    res = MulticomponentResult(x=x, C_final=C, t_final=0.0, mass_initial=C @ vol)
    if t_total == 0:
        res.mass_final = res.mass_initial
        return res

    def matrices(C, T):
        D = np.asarray(D_func(C, T), dtype=float)
        if D.shape != (n, n, N) or not np.all(np.isfinite(D)):
            raise ValueError(f"D_func must return finite matrices of shape {(n, n, N)}")
        return D

    D0 = matrices(C, float(history.T(0.0)))
    rate0 = float(np.max(np.sum(np.abs(D0), axis=1)))           # Gershgorin bound
    if rate0 <= 0:
        raise ValueError("the diffusion matrix must not vanish")
    dt_step = min(courant * dx ** 2 / rate0, t_total)
    dt_cap = max(dt_step, t_total / max(min_steps, 2))
    if t_total / dt_cap > max_steps:
        dt_cap = t_total / max_steps
        res.warnings.append("the step count limit enlarged the time step; check against a finer grid")
    fixed_nodes = [k for k, b in ((0, left), (N - 1, right)) if b is not None]
    fixed_rows = np.array([k * n + i for k in fixed_nodes for i in range(n)], dtype=int)
    bw = 2 * n - 1                                              # half bandwidth, interleaved order
    size = n * N
    t, step = 0.0, 0
    flat = C.T.ravel()                                          # node-major: k*n + i
    while t < t_total * (1 - 1e-12):
        if step >= max_steps:
            raise RuntimeError("max_steps reached before the requested duration")
        D = matrices(flat.reshape(N, n).T, float(history.T(t)))
        dt = min(dt_step * dt_growth, dt_cap, t_total - t)
        dt_step = dt
        # A = I - theta dt L in banded storage; rows of fixed nodes stay identity rows
        Lflat = np.zeros(size)
        ab = np.zeros((2 * bw + 1, size))
        for rows, cols, vals in _operator(x, D, m, fixed_rows):
            Lflat += np.bincount(rows, vals * flat[cols], minlength=size)
            ab[bw + rows - cols, cols] -= theta_time * dt * vals
        ab[bw] += 1.0
        rhs = flat + (1.0 - theta_time) * dt * Lflat
        t_new = t + dt
        for k, b in ((0, left), (N - 1, right)):
            v = _fixed(b, t_new, n)
            if v is not None:
                rhs[k * n:(k + 1) * n] = v
        flat = solve_banded((bw, bw), ab, rhs)
        t = t_new
        step += 1
        if progress is not None and step % 50 == 0 and progress(t / t_total):
            res.warnings.append("aborted by user")
            break
    C = flat.reshape(N, n).T.copy()
    res.C_final, res.t_final, res.n_steps = C, t, step
    res.mass_final = C @ vol
    res.min_value = float(C.min())
    return res


def step_constant_matrix(x, x0: float, C_left, C_right, D: np.ndarray, t: float) -> np.ndarray:
    """Infinite-medium step solution for a constant matrix, via its eigen-components.

    Returns shape (n-1, len(x)). The eigenvalues must be real and positive,
    which holds for the matrices of an ideal ionic solution in the tests.
    """
    D = np.asarray(D, dtype=float)
    lam, P = np.linalg.eig(D)
    if np.any(np.abs(lam.imag) > 1e-12 * np.max(np.abs(lam))) or np.any(lam.real <= 0):
        raise ValueError("the matrix has complex or nonpositive eigenvalues")
    lam, P = lam.real, P.real
    Pinv = np.linalg.inv(P)
    uL, uR = Pinv @ np.asarray(C_left, float), Pinv @ np.asarray(C_right, float)
    u = np.array([step_infinite(np.asarray(x, float), x0, uL[k], uR[k], lam[k] * t) for k in range(len(lam))])
    return P @ u
