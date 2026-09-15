"""Geometry of the modelled domain.

The 1-D diffusion equation is solved in the general form

    dC/dt = (1/x^m) d/dx [ x^m D dC/dx ]

with the geometry index ``m``:

==========  ===  ===========================================================
plane       0    a traverse across a slab / diffusion couple (Crank section 2, section 4)
cylinder    1    radius of an infinite cylinder (Crank section 5)
sphere      2    radius of a sphere (Crank section 6)
==========  ===  ===========================================================

Choosing the geometry is a petrological decision, not a numerical one.  A
traverse from rim to core of a prismatic crystal that is long compared with
the diffusion length is a plane-sheet problem; a small equant crystal that has
equilibrated from all sides is better described as a sphere.  Modelling a 3-D
crystal as 1-D always *over*-estimates the time, and sectioning effects add
further bias (Costa et al. 2008; Shea et al. 2015; Krimer & Costa 2017).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict

import numpy as np

GEOMETRIES: Dict[str, int] = {"plane": 0, "cylinder": 1, "sphere": 2}

GEOMETRY_NOTES = {
    "plane": ("Plane sheet / semi-infinite medium. Use for a traverse across a zone boundary "
              "in a crystal that is large compared with the diffusion length "
              "(Crank 1975 section 2 and section 4)."),
    "cylinder": ("Infinite cylinder, x is the radius. Use for a prismatic crystal equilibrating "
                 "radially with negligible transport along the prism axis (Crank 1975 section 5)."),
    "sphere": ("Sphere, x is the radius. Use for a small equant crystal equilibrating from all "
               "sides (Crank 1975 section 6)."),
}


@dataclass
class Geometry:
    kind: str = "plane"

    def __post_init__(self):
        if self.kind not in GEOMETRIES:
            raise ValueError(f"Unknown geometry '{self.kind}'. Choose from {list(GEOMETRIES)}.")

    @property
    def m(self) -> int:
        return GEOMETRIES[self.kind]

    @property
    def is_radial(self) -> bool:
        return self.m > 0

    def describe(self) -> str:
        return GEOMETRY_NOTES[self.kind]


def make_grid(x_min: float, x_max: float, n_nodes: int) -> np.ndarray:
    """Uniform grid; the numerical solver requires uniform spacing."""
    if n_nodes < 5:
        raise ValueError("use at least 5 nodes")
    return np.linspace(float(x_min), float(x_max), int(n_nodes))


def suggest_grid(x_data: np.ndarray, n_nodes: int = 401, pad: float = 0.0) -> np.ndarray:
    """Grid spanning the measured traverse, optionally padded at both ends."""
    x_data = np.asarray(x_data, dtype=float)
    lo, hi = float(x_data.min()), float(x_data.max())
    span = hi - lo
    return make_grid(lo - pad * span, hi + pad * span, n_nodes)


def stability_dt(dx: float, D_max: float, courant: float = 0.5) -> float:
    """Explicit-scheme stability limit ``dt <= dx^2 / (2 D)`` (Crank 1975 eq. 8.33).

    ``courant`` scales it; the DMG Short Course 2025 scripts use 0.2-0.4.
    """
    return courant * dx ** 2 / D_max
