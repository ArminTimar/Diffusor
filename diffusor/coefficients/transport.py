"""Relations between tracer coefficients and exchange (interdiffusion) coefficients.

A tracer coefficient D*_i describes how fast species i moves in a host whose
chemical composition does not change. A zoning profile of major components
relaxes by exchange: on a site occupied by n equally charged ions only n-1
compositions are independent, and the flux of each depends on the gradients
of all of them. For an ideal ionic solution the Fick matrix follows from the
tracer coefficients (Lasaga 1979, as printed by Chakraborty & Ganguly 1992,
eq. 2, and by Oeser et al. 2026, eq. 4):

    D_ij = D*_i delta_ij - D*_i X_i (D*_j - D*_n) / sum_k(X_k D*_k)

with component n dependent. For two components this reduces to the
binary interdiffusion coefficient (Oeser et al. 2026, eq. 6; Borinski et al.
2012, eq. 8, for the ideal binary, with their eq. 1 as the ideal multicomponent
matrix; their eq. 2 is the non-ideal matrix)

    D_AB = D*_A D*_B / (X_A D*_A + X_B D*_B)

which lies between the two tracer coefficients and approaches the tracer
coefficient of the minor species as that species becomes dilute. A non-ideal
solution multiplies the binary form by the thermodynamic factor
1 + dln(gamma)/dln(X) (Schaffer et al. 2014, eq. 3).

Only equal valences are handled: the printed forms of the general equation
differ in where the valences enter, and every application here (Fe, Mg, Mn
and Ca in garnet; Fe and Mg isotopes in olivine) exchanges divalent ions.

A law built here is a calculation from two fitted tracer laws, not a fit of its
own. Its composition range (``X_range``) is the mathematical domain of the
formula unless the caller passes the range the tracer laws were fitted over; the
tracer coefficients are assumed independent of composition.

The tracer coefficients of a family must have been extracted with the same
formalism they are used in. Chakraborty & Ganguly (1992, p. 79) make this
point for garnet: their D* were fitted with the ideal model, so inverting them
with the ideal model returns a consistent D matrix.
"""
from __future__ import annotations

from typing import Dict, Optional, Sequence

import numpy as np

from .base import DiffusionCoefficient, Parameter, Range


def binary_interdiffusion(D_A, D_B, X_A, thermodynamic_factor=1.0):
    """Interdiffusion coefficient of an equal-valence binary exchange.

    ``D_A`` and ``D_B`` are the tracer coefficients, ``X_A`` the mole fraction
    of A on the exchange site. Arrays broadcast. The thermodynamic factor is 1
    for an ideal solution.
    """
    D_A = np.asarray(D_A, dtype=float)
    D_B = np.asarray(D_B, dtype=float)
    X_A = np.asarray(X_A, dtype=float)
    if np.any(D_A <= 0) or np.any(D_B <= 0):
        raise ValueError("tracer coefficients must be positive")
    if np.any((X_A < 0) | (X_A > 1)):
        raise ValueError("X_A must lie between 0 and 1")
    return D_A * D_B / (X_A * D_A + (1.0 - X_A) * D_B) * thermodynamic_factor


def ideal_ionic_matrix(D_tracer, X, dependent: int = -1) -> np.ndarray:
    """Fick diffusion matrix of an ideal ionic solution of equally charged ions.

    ``D_tracer`` and ``X`` have shape ``(n, ...)``: one row per component,
    trailing dimensions (for example the nodes of a grid) broadcast. The
    result has shape ``(n-1, n-1, ...)`` over the independent components in
    their original order with ``dependent`` removed. Mole fractions are
    normalised to sum to one.
    """
    D = np.asarray(D_tracer, dtype=float)
    X = np.asarray(X, dtype=float)
    if D.shape != X.shape or D.ndim < 1 or D.shape[0] < 2:
        raise ValueError("D_tracer and X need the same shape (n, ...) with n >= 2")
    if np.any(D <= 0) or not np.all(np.isfinite(D)):
        raise ValueError("tracer coefficients must be finite and positive")
    n = D.shape[0]
    dep = dependent % n
    total = X.sum(axis=0)
    if np.any(total <= 0):
        raise ValueError("mole fractions must have a positive sum")
    X = X / total
    idx = [k for k in range(n) if k != dep]
    weight = np.sum(D * X, axis=0)
    out = np.empty((n - 1, n - 1) + D.shape[1:])
    for a, i in enumerate(idx):
        for b, j in enumerate(idx):
            out[a, b] = (D[i] if i == j else 0.0) - D[i] * X[i] * (D[j] - D[dep]) / weight
    return out


def effective_binary(D_matrix: np.ndarray, gradients: np.ndarray, i: int) -> np.ndarray:
    """Effective binary diffusion coefficient of component ``i`` (Chakraborty &
    Ganguly 1992, eq. 4.3): ``D_i(EB) = sum_j D_ij (dC_j/dx) / (dC_i/dx)``.

    It depends on the local gradients, so it is not transferable between
    crystals (their p. 81); Diffusor uses it only for reporting.
    """
    D_matrix = np.asarray(D_matrix, dtype=float)
    g = np.asarray(gradients, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.einsum("j...,j...->...", D_matrix[i], g) / g[i]


# ---------------------------------------------------------------------------
def interdiffusion_from_tracers(key: str, *, tracer_a: DiffusionCoefficient,
                                tracer_b: DiffusionCoefficient, species: str, comp_key: str,
                                label: str, citation: str, equation_number: str,
                                X_range: Optional[Range] = None, notes: str = "",
                                calibration_notes: Sequence[str] = (),
                                verified_from: str = "", recommended: bool = False) -> DiffusionCoefficient:
    """A scalar interdiffusion law computed from two tracer laws.

    ``cond.X[comp_key]`` is the mole fraction of the species of ``tracer_a``
    on the exchange site. Both tracer laws are evaluated at the same
    conditions (temperature, pressure, fO2, composition and direction), so a
    composition term inside a tracer law sees the same ``comp_key``. Where both
    tracer laws have independent principal-axis laws, the derived law has them
    too, and a direction between the axes is projected from the derived
    principal values.
    """
    if tracer_a.kind != "tracer" or tracer_b.kind != "tracer":
        raise ValueError("both source laws must be tracer coefficients")
    if tracer_a.mineral != tracer_b.mineral:
        raise ValueError("the two tracer laws must be for the same mineral")

    def split(p: Dict[str, float]):
        a = {k[2:]: v for k, v in p.items() if k.startswith("A:")}
        b = {k[2:]: v for k, v in p.items() if k.startswith("B:")}
        return a, b

    def combine(fa, fb):
        def evaluate(dc, cond, p):
            pa, pb = split(p)
            return binary_interdiffusion(fa(tracer_a, cond, pa), fb(tracer_b, cond, pb),
                                         np.asarray(cond.x(comp_key), dtype=float))
        return evaluate

    params = {}
    for prefix, src in (("A:", tracer_a), ("B:", tracer_b)):
        for name, p in src.params.items():
            params[prefix + name] = Parameter(prefix + name, p.value, p.sigma, p.unit,
                                              p.sigma_level, f"{src.species} tracer: {p.description}")
    principal = {}
    if tracer_a.principal_funcs and tracer_b.principal_funcs:
        for axis in sorted(set(tracer_a.principal_funcs) & set(tracer_b.principal_funcs)):
            principal[axis] = combine(tracer_a.principal_funcs[axis], tracer_b.principal_funcs[axis])
    requires = tuple(dict.fromkeys((comp_key, *tracer_a.requires, *tracer_b.requires)))
    T_lo = max(v for v in (tracer_a.T_range.lo, tracer_b.T_range.lo) if v is not None)
    T_hi = min(v for v in (tracer_a.T_range.hi, tracer_b.T_range.hi) if v is not None)
    eq = (f"D_{species} = D*_A D*_B / (X_A D*_A + X_B D*_B), A = {tracer_a.species}, "
          f"B = {tracer_b.species}, X_A = {comp_key}; D*_A: {tracer_a.equation_text}; "
          f"D*_B: {tracer_b.equation_text}")
    return DiffusionCoefficient(
        key=key, mineral=tracer_a.mineral, species=species, label=label, citation=citation,
        kind="interdiffusion", transported_variable=f"{comp_key}, mole fraction of {tracer_a.species}",
        equation_text=eq, equation_number=equation_number,
        func=combine(tracer_a.func, tracer_b.func), principal_funcs=principal,
        params=params, requires=requires,
        needs_fo2=tracer_a.needs_fo2 or tracer_b.needs_fo2, fo2_unit=tracer_a.fo2_unit,
        T_range=Range(T_lo, T_hi, tracer_a.T_range.unit), P_range=tracer_a.P_range,
        fo2_range=tracer_a.fo2_range, X_range=X_range or tracer_a.X_range,
        reference_axis=tracer_a.reference_axis, allowed_axes=tracer_a.allowed_axes,
        orientation_required=tracer_a.orientation_required or tracer_b.orientation_required,
        verified=tracer_a.verified and tracer_b.verified,
        verified_from=verified_from or f"computed from {tracer_a.key} and {tracer_b.key}",
        reference_state=tracer_a.reference_state, calibration_notes=tuple(calibration_notes),
        uncertainty_note=" ".join(dict.fromkeys(n for n in (tracer_a.uncertainty_note,
                                                            tracer_b.uncertainty_note) if n)),
        notes=notes, recommended=recommended, derived_from=(tracer_a.key, tracer_b.key),
    )
