"""Unit and composition conversions.

Temperature / pressure / time helpers plus the mineral composition variables
used by the diffusion coefficients:

* ``X_Fe = Fe / (Fe + Mg)``  (molar; Fe = total Fe as Fe2+ unless stated)
  -- convention of Dohmen & Chakraborty (2007), Dohmen et al. (2016),
  Mueller et al. (2013).
* ``Fo = 100 * Mg / (Mg + Fe)`` mol%; ``Mg# = 100 * Mg / (Mg + Fe)``.
* ``X_An = Ca / (Ca + Na + K)`` (molar), Deer, Howie & Zussman (1992).
* ``X_Usp`` (ulvoespinel fraction) for titanomagnetite from Ti apfu on a
  3-cation basis: ``X_Usp = 3 Ti / (Ti + Fe)`` i.e. Ti per formula unit
  (Fe3-xTixO4 => X_Usp = x); Stormer (1983)-style recalculation is *not*
  applied here -- users should supply Ti apfu or X_Usp directly for Fe-Ti work.

wt% oxide -> cation mole conversions use the molar masses in
:mod:`diffusor.constants` (IUPAC 2021).
"""
from __future__ import annotations

import numpy as np

from ..constants import (T_KELVIN_OFFSET, BAR_TO_PA, GPA_TO_PA, KBAR_TO_PA,
                         SEC_PER_MIN, SEC_PER_HOUR, SEC_PER_DAY, SEC_PER_YEAR,
                         OXIDE_MOLAR_MASS, OXIDE_CATIONS)


# --- temperature ---------------------------------------------------------------
def c_to_k(T_C):
    return np.asarray(T_C, dtype=float) + T_KELVIN_OFFSET


def k_to_c(T_K):
    return np.asarray(T_K, dtype=float) - T_KELVIN_OFFSET


# --- pressure ------------------------------------------------------------------
_P_FACTORS = {"Pa": 1.0, "bar": BAR_TO_PA, "kbar": KBAR_TO_PA, "MPa": 1.0e6, "GPa": GPA_TO_PA}


def pressure_to_pa(value, unit: str) -> float:
    try:
        return float(value) * _P_FACTORS[unit]
    except KeyError as exc:
        raise ValueError(f"Unknown pressure unit {unit}; use one of {list(_P_FACTORS)}") from exc


# --- time ----------------------------------------------------------------------
_T_FACTORS = {"s": 1.0, "min": SEC_PER_MIN, "h": SEC_PER_HOUR, "d": SEC_PER_DAY,
              "yr": SEC_PER_YEAR, "kyr": 1.0e3 * SEC_PER_YEAR, "Myr": 1.0e6 * SEC_PER_YEAR}


def time_to_seconds(value, unit: str):
    return np.asarray(value, dtype=float) * _T_FACTORS[unit]


def seconds_to(value, unit: str):
    return np.asarray(value, dtype=float) / _T_FACTORS[unit]


def human_time(seconds: float) -> str:
    """Pick a readable unit for a duration in seconds."""
    s = float(seconds)
    for unit in ("Myr", "kyr", "yr", "d", "h", "min"):
        v = s / _T_FACTORS[unit]
        if v >= 1.0:
            return f"{v:.3g} {unit}"
    return f"{s:.3g} s"


# --- length ----------------------------------------------------------------------
_L_FACTORS = {"m": 1.0, "mm": 1e-3, "um": 1e-6, "µm": 1e-6, "nm": 1e-9, "cm": 1e-2}


def length_to_m(value, unit: str):
    return np.asarray(value, dtype=float) * _L_FACTORS[unit]


# --- compositions --------------------------------------------------------------
def oxide_wt_to_cation_moles(wt: np.ndarray, oxide: str) -> np.ndarray:
    """Cation moles per 100 g from oxide wt% (IUPAC 2021 molar masses)."""
    return np.asarray(wt, dtype=float) / OXIDE_MOLAR_MASS[oxide] * OXIDE_CATIONS[oxide]


def ratio_A_over_AplusB(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return a / (a + b)


def x_fe_from_oxides(feo_wt, mgo_wt):
    """X_Fe = Fe/(Fe+Mg) molar from FeO(total) and MgO wt%."""
    fe = oxide_wt_to_cation_moles(feo_wt, "FeO")
    mg = oxide_wt_to_cation_moles(mgo_wt, "MgO")
    return ratio_A_over_AplusB(fe, mg)


def x_an_from_oxides(cao_wt, na2o_wt, k2o_wt=0.0):
    """X_An = Ca/(Ca+Na+K) molar from oxide wt%."""
    ca = oxide_wt_to_cation_moles(cao_wt, "CaO")
    na = oxide_wt_to_cation_moles(na2o_wt, "Na2O")
    k = oxide_wt_to_cation_moles(k2o_wt, "K2O")
    return ca / (ca + na + k)


def composition_variable(a, b, mode: str, a_oxide: str | None = None, b_oxide: str | None = None):
    """Build the modelled variable from columns A and B.

    mode:
      'A'            -> use column A as is (single-element profile)
      'A/(A+B)'      -> molar ratio if oxides are named, else raw ratio of the columns
      'B/(A+B)'      -> complement
      'A-B'          -> difference (rarely used)
    """
    a = np.asarray(a, dtype=float)
    if mode == "A":
        return a
    b = np.asarray(b, dtype=float)
    if a_oxide and b_oxide:
        a = oxide_wt_to_cation_moles(a, a_oxide)
        b = oxide_wt_to_cation_moles(b, b_oxide)
    if mode == "A/(A+B)":
        return a / (a + b)
    if mode == "B/(A+B)":
        return b / (a + b)
    if mode == "A-B":
        return a - b
    raise ValueError(f"Unknown composition mode {mode}")
