"""Families of tracer coefficients that together define a multicomponent D matrix.

A :class:`TracerFamily` holds one published set of tracer laws, one per
exchanging component, that were fitted together with one diffusion model. The
multicomponent solver (:mod:`diffusor.solvers.multicomponent`) turns them into
the Fick matrix with :func:`diffusor.coefficients.transport.ideal_ionic_matrix`
at every node and time step, because the matrix changes with composition.

Families are kept apart from the scalar registry: a single tracer coefficient
from a family is not a law for a single zoning profile, and mixing the laws of
two families breaks the consistency that the families were fitted under
(Chakraborty & Ganguly 1992, p. 79).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from ..references import cite, get as get_reference
from .base import Range
from .transport import ideal_ionic_matrix

# law(T_K, P_Pa, X, options) -> D* in m2/s; X maps component -> mole fraction array
TracerLaw = Callable[[float, float, Mapping[str, np.ndarray], Mapping[str, float]], np.ndarray]


@dataclass
class TracerFamily:
    key: str
    label: str
    mineral: str
    components: Tuple[str, ...]
    citation: str
    equation_text: str
    equation_number: str
    laws: Dict[str, TracerLaw]
    T_range: Range = field(default_factory=Range)
    P_range: Range = field(default_factory=Range)
    verified: bool = False
    verified_from: str = ""
    sigma_logD: Optional[float] = None          # 1 sigma of log10 D* of every component
    uncertainty_note: str = ""
    calibration_notes: Sequence[str] = ()
    notes: str = ""
    secondary_citations: Sequence[str] = ()
    options: Dict[str, float] = field(default_factory=dict)        # name -> default
    option_help: Dict[str, str] = field(default_factory=dict)
    composition_labels: Dict[str, str] = field(default_factory=dict)
    default_dependent: str = ""
    recommended: bool = False
    # extra composition checks: (description, test(X) -> bool for "inside")
    composition_checks: Sequence[Tuple[str, Callable[[Mapping[str, np.ndarray]], bool]]] = ()

    def __post_init__(self):
        if set(self.laws) != set(self.components):
            raise ValueError(f"{self.key}: one tracer law per component is required")
        if self.default_dependent and self.default_dependent not in self.components:
            raise ValueError(f"{self.key}: unknown dependent component {self.default_dependent!r}")

    # -- evaluation ---------------------------------------------------------
    def _options(self, options: Optional[Mapping[str, float]]) -> Dict[str, float]:
        out = dict(self.options)
        for k, v in (options or {}).items():
            if k not in out:
                raise KeyError(f"{self.key} has no option {k!r} (has {sorted(out)})")
            out[k] = float(v)
        return out

    def _composition(self, X: Mapping[str, object]) -> Dict[str, np.ndarray]:
        missing = [c for c in self.components if c not in X]
        if missing:
            raise KeyError(f"{self.key} needs the mole fractions of {', '.join(missing)}")
        arr = {c: np.asarray(X[c], dtype=float) for c in self.components}
        total = sum(arr.values())
        if np.any(total <= 0):
            raise ValueError("mole fractions must have a positive sum")
        return {c: v / total for c, v in arr.items()}

    def tracer(self, T_K: float, P_Pa: float, X: Mapping[str, object],
               options: Optional[Mapping[str, float]] = None,
               log10_offsets: Optional[Mapping[str, float]] = None) -> Dict[str, np.ndarray]:
        """Tracer coefficient of every component, m2/s."""
        if not np.isfinite(T_K) or T_K <= 0 or not np.isfinite(P_Pa) or P_Pa < 0:
            raise ValueError("temperature must be positive and pressure nonnegative")
        Xn = self._composition(X)
        opts = self._options(options)
        out = {}
        for c in self.components:
            D = np.asarray(self.laws[c](T_K, P_Pa, Xn, opts), dtype=float)
            off = (log10_offsets or {}).get(c, 0.0)
            out[c] = D * 10.0 ** off if off else D
        return out

    def matrix(self, T_K: float, P_Pa: float, X: Mapping[str, object], dependent: Optional[str] = None,
               options=None, log10_offsets=None) -> np.ndarray:
        """Ideal ionic Fick matrix over the independent components (dependent removed)."""
        dep = dependent or self.default_dependent or self.components[-1]
        D = self.tracer(T_K, P_Pa, X, options, log10_offsets)
        Xn = self._composition(X)
        order = list(self.components)
        shape = np.broadcast(*[np.asarray(v) for v in D.values()], *[np.asarray(v) for v in Xn.values()]).shape
        Dm = np.stack([np.broadcast_to(D[c], shape) for c in order])
        Xm = np.stack([np.broadcast_to(Xn[c], shape) for c in order])
        return ideal_ionic_matrix(Dm, Xm, dependent=order.index(dep))

    def latex(self) -> str:
        from .latex import family_latex
        return family_latex(self)

    # -- reporting ----------------------------------------------------------
    def check(self, T_K: float, P_Pa: float, X: Optional[Mapping[str, object]] = None) -> List[str]:
        w: List[str] = []
        if not self.T_range.contains(T_K):
            w.append(f"T = {T_K - 273.15:.0f} C is outside the calibration range ({self.T_range}) "
                     f"of {cite(self.citation)}")
        if not self.P_range.contains(P_Pa):
            w.append(f"P = {P_Pa / 1e9:.2f} GPa is outside the calibration range ({self.P_range}) "
                     f"of {cite(self.citation)}")
        if X is not None:
            Xn = self._composition(X)
            for text, inside in self.composition_checks:
                if not inside(Xn):
                    w.append(text)
        if not self.verified:
            w.append(f"{self.key}: tracer laws NOT verified against the primary publication.")
        w.extend(self.calibration_notes)
        if self.uncertainty_note:
            w.append(self.uncertainty_note)
        return w

    def describe(self) -> str:
        lines = [self.label, f"  mineral         : {self.mineral}",
                 f"  components      : {', '.join(self.components)}",
                 f"  source          : {get_reference(self.citation).full()}",
                 f"  equation        : {self.equation_number}", f"     {self.equation_text}",
                 f"  LaTeX           : ${self.latex()}$",
                 "  D matrix        : ideal ionic solution of equally charged ions "
                 f"({cite('lasaga1979')}; {cite('chakraborty_ganguly1992', 'eq. 2')})",
                 f"  calibration     : T {self.T_range}, P {self.P_range}",
                 f"  verified        : {'YES -- ' + self.verified_from if self.verified else 'NO'}"]
        if self.sigma_logD is not None:
            lines.append(f"  scatter         : {self.sigma_logD:g} log10 units (1 sigma) for each D*")
        for k, v in self.options.items():
            lines.append(f"  option {k} = {v:g}: {self.option_help.get(k, '')}")
        lines.extend(f"  limitation      : {s}" for s in self.calibration_notes)
        if self.uncertainty_note:
            lines.append(f"  uncertainty     : {self.uncertainty_note}")
        if self.secondary_citations:
            lines.append("  see also        : " + ", ".join(cite(c) for c in self.secondary_citations))
        if self.notes:
            lines.append(f"  notes           : {self.notes}")
        return "\n".join(lines)

    def methods_sentence(self) -> str:
        return (f"Multicomponent diffusion of {', '.join(self.components)} in {self.mineral} was modelled "
                f"with the tracer coefficients of {get_reference(self.citation).short()} "
                f"({self.equation_number}), combined into a composition-dependent diffusion matrix for an "
                f"ideal ionic solution ({get_reference('lasaga1979').short()}; "
                f"{get_reference('chakraborty_ganguly1992').short()}, eq. 2).")


FAMILIES: Dict[str, TracerFamily] = {}


def register(family: TracerFamily) -> TracerFamily:
    if family.key in FAMILIES:
        raise RuntimeError(f"duplicate tracer family {family.key}")
    FAMILIES[family.key] = family
    return family


def get_family(key: str) -> TracerFamily:
    try:
        return FAMILIES[key]
    except KeyError as exc:
        raise KeyError(f"Unknown tracer family '{key}'. Available: {sorted(FAMILIES)}") from exc


def families_for(mineral: str) -> List[TracerFamily]:
    return sorted((f for f in FAMILIES.values() if f.mineral == mineral),
                  key=lambda f: (not f.recommended, f.key))
