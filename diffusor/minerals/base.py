"""Mineral descriptions: crystallographic axes, composition variables, species.

A :class:`Mineral` tells the rest of Diffusor

* which crystallographic axes exist and how a traverse is described
  (direction cosines to a, b, c -- Costa & Chakraborty 2004; DMG Short Course
  2025 Lecture 6, ``D_V = D_a cos^2(alpha) + D_b cos^2(beta) + D_c cos^2(gamma)``),
* which composition variable the diffusion coefficients depend on
  (X_Fe, X_An, X_Fo, x_Ti ...), and how to build it from two measured columns,
* which species (diffusing components) can be modelled.

Composition conventions follow Deer, Howie & Zussman (1992) and the usage of
the diffusion papers themselves (see :mod:`diffusor.coefficients`).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np


@dataclass(frozen=True)
class CompositionVariable:
    """How to turn the measured columns into the variable D depends on."""
    key: str                      # 'XFe', 'XAn', 'XFo', 'xTi', 'none'
    label: str                    # human readable
    definition: str               # e.g. 'Fe/(Fe+Mg) molar'
    default_mode: str = "A/(A+B)"  # composition_variable() mode
    citation: str = "deer1992"

    def __str__(self) -> str:
        return f"{self.label} = {self.definition}"


@dataclass(frozen=True)
class Species:
    """A diffusing component of a mineral."""
    key: str                      # 'Fe-Mg', 'Ni', 'Sr', 'Mg', 'Fe-Ti'
    label: str
    kind: str = "interdiffusion"  # 'interdiffusion' | 'tracer' | 'trace'
    note: str = ""


@dataclass(frozen=True)
class Mineral:
    key: str
    name: str
    formula: str
    system: str                                   # crystal system
    axes: Sequence[str]                           # e.g. ('a', 'b', 'c')
    species: Dict[str, Species]
    composition_variable: CompositionVariable
    typical_habit: str = ""
    notes: str = ""
    isotropic: bool = False

    def species_keys(self) -> List[str]:
        return list(self.species)

    def get_species(self, key: str) -> Species:
        if key not in self.species:
            raise KeyError(f"{self.name} has no species '{key}'. Available: {self.species_keys()}")
        return self.species[key]


def direction_factor(D_a, D_b, D_c, alpha_deg: float, beta_deg: float, gamma_deg: float):
    """Diffusivity along a traverse making angles alpha/beta/gamma with a, b, c.

    ``D_V = D_a cos^2(alpha) + D_b cos^2(beta) + D_c cos^2(gamma)``

    Source: Costa & Chakraborty (2004) EPSL 227, eq. for diffusion along an
    arbitrary direction in an orthorhombic crystal; DMG Short Course 2025
    Lecture 6 "Crystal orientation" (also used in the course script
    ``MCdiff_OlFo_MO.m``).  The angles must satisfy
    ``cos^2 a + cos^2 b + cos^2 g = 1`` for orthogonal axes.
    """
    ca = np.cos(np.radians(alpha_deg)) ** 2
    cb = np.cos(np.radians(beta_deg)) ** 2
    cg = np.cos(np.radians(gamma_deg)) ** 2
    s = ca + cb + cg
    if not np.isclose(s, 1.0, atol=0.05):
        raise ValueError(
            f"direction cosines do not close: cos^2 sum = {s:.3f} (should be 1.0). "
            "Check the measured angles between the traverse and the a, b, c axes.")
    return D_a * ca + D_b * cb + D_c * cg
