"""Boundary conditions for the 1-D solvers.

Types (Crank 1975, section 8.6 "Other boundary conditions"; Costa et al. 2008 section "Boundary conditions"):

* ``dirichlet``  -- fixed concentration at the boundary (crystal rim in
  contact with an infinite melt reservoir; "open system").  The value may be
  a constant or a function of time ``f(t_seconds) -> C``.
* ``neumann``    -- zero flux (closed system rim, or the symmetry plane /
  centre of a crystal).  Equivalent to a mirror ghost node.
* ``symmetry``   -- alias of zero flux used at r = 0 for cylinder / sphere.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional, Union

Number = float
ValueOrFunc = Union[Number, Callable[[float], Number]]


@dataclass
class BoundaryCondition:
    kind: str = "neumann"            # 'dirichlet' | 'neumann' | 'symmetry'
    value: Optional[ValueOrFunc] = None

    def __post_init__(self):
        if self.kind not in ("dirichlet", "neumann", "symmetry"):
            raise ValueError(f"Unknown BC kind {self.kind}")
        if self.kind == "dirichlet" and self.value is None:
            raise ValueError("Dirichlet BC needs a value or a function of time")

    def value_at(self, t: float) -> float:
        if callable(self.value):
            return float(self.value(t))
        return float(self.value)

    @property
    def is_flux_free(self) -> bool:
        return self.kind in ("neumann", "symmetry")

    def describe(self) -> str:
        if self.kind == "dirichlet":
            v = "f(t)" if callable(self.value) else f"{self.value:g}"
            return f"fixed concentration = {v} (open system)"
        return "zero flux (closed system / symmetry)"


def dirichlet(value: ValueOrFunc) -> BoundaryCondition:
    return BoundaryCondition("dirichlet", value)


def zero_flux() -> BoundaryCondition:
    return BoundaryCondition("neumann")


def symmetry() -> BoundaryCondition:
    return BoundaryCondition("symmetry")
