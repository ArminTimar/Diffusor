"""Lookup of all diffusion coefficients available in Diffusor."""
from __future__ import annotations

from typing import Dict, List, Optional

from .base import Conditions, DiffusionCoefficient, Parameter, Range
from . import cpx, magnetite, olivine, opx, plagioclase

_MODULES = (olivine, opx, cpx, plagioclase, magnetite)

REGISTRY: Dict[str, DiffusionCoefficient] = {}
for _m in _MODULES:
    for _c in _m.COEFFICIENTS:
        if _c.key in REGISTRY:
            raise RuntimeError(f"duplicate coefficient key {_c.key}")
        REGISTRY[_c.key] = _c


def get(key: str) -> DiffusionCoefficient:
    try:
        return REGISTRY[key]
    except KeyError as exc:
        raise KeyError(f"Unknown diffusion coefficient '{key}'. "
                       f"Use list_coefficients() to see the options.") from exc


def list_coefficients(mineral: Optional[str] = None, species: Optional[str] = None,
                      verified_only: bool = False) -> List[DiffusionCoefficient]:
    out = [c for c in REGISTRY.values()
           if (mineral is None or c.mineral == mineral)
           and (species is None or c.species == species)
           and (not verified_only or c.verified)]
    # recommended first, then verified, then alphabetical
    return sorted(out, key=lambda c: (not c.recommended, not c.verified, c.key))


def species_for(mineral: str) -> List[str]:
    return sorted({c.species for c in REGISTRY.values() if c.mineral == mineral})


def minerals() -> List[str]:
    return sorted({c.mineral for c in REGISTRY.values()})


def summary_table() -> str:
    rows = [("key", "mineral", "species", "source", "verified")]
    for c in list_coefficients():
        from ..references import get as get_ref
        rows.append((c.key, c.mineral, c.species, get_ref(c.citation).short(),
                     "yes" if c.verified else "NO"))
    w = [max(len(r[i]) for r in rows) for i in range(5)]
    lines = ["  ".join(r[i].ljust(w[i]) for i in range(5)) for r in rows]
    lines.insert(1, "  ".join("-" * w[i] for i in range(5)))
    return "\n".join(lines)
