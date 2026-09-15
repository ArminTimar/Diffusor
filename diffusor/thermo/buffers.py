"""Oxygen fugacity buffers.

All buffers are returned as ``log10 fO2`` (fO2 in bar) for a temperature in K
and pressure in Pa.

Two parameterisations are offered:

``frost1991`` (default)
    Compilation form ``log fO2 = A/T + B + C (P - 1)/T`` with *P* in bar,
    Frost (1991), Reviews in Mineralogy 25, Table 1 (citekey ``frost1991``).
    Coefficients (A, B, C):

    ==== ========== ======= =====
    IW   -27489     6.702   0.055
    WM   -32807     13.012  0.083
    FMQ  -25096.3   8.735   0.110
    NNO  -24930     9.36    0.046
    HM   -25700.6   14.558  0.019
    ==== ========== ======= =====

``oneill``
    Free-energy expressions of the primary calibrations, 1 bar:

    * NNO: O'Neill & Pownceby (1993), CMP 114, eq. for 2 Ni + O2 = 2 NiO,
      ``dG(J/mol) = -478967 + 248.514 T - 9.7961 T ln T`` (citekey
      ``oneill_pownceby1993``), so ``log fO2 = dG / (R T ln 10)``.
    * FMQ: O'Neill (1987), Am Mineral 72, eq. for 3 Fe2SiO4 + O2 = 2 Fe3O4 + 3 SiO2,
      ``dG(J/mol) = -587474 + 1584.427 T - 203.3164 T ln T + 0.092710 T^2``
      (citekey ``oneill1987``).
    The Frost (1991) pressure term ``C (P-1)/T`` is added to these so that both
    parameterisations agree on the pressure dependence.

Note on correlation: in Monte Carlo runs the user supplies fO2 as an offset
from a buffer (e.g. NNO+1).  The buffer is re-evaluated at every sampled
temperature, so fO2 and T are correlated by construction; this is the
dependence that NIDIS (Petrone et al. 2016) ignores when it propagates errors
as if the variables were independent.
"""
from __future__ import annotations

import math
from typing import Dict, Tuple

from ..constants import R_GAS, BAR_TO_PA

LN10 = math.log(10.0)

# Frost (1991) Table 1 coefficients: A, B, C
_FROST: Dict[str, Tuple[float, float, float]] = {
    "IW":  (-27489.0, 6.702, 0.055),
    "WM":  (-32807.0, 13.012, 0.083),
    "FMQ": (-25096.3, 8.735, 0.110),
    "NNO": (-24930.0, 9.36, 0.046),
    "HM":  (-25700.6, 14.558, 0.019),
}
_ALIASES = {"QFM": "FMQ", "IQF": "IW"}

BUFFER_CITATIONS = {
    "frost1991": "frost1991",
    "NNO_oneill": "oneill_pownceby1993",
    "FMQ_oneill": "oneill1987",
}


def available_buffers():
    return list(_FROST) + ["QFM"]


def _canon(name: str) -> str:
    n = name.upper()
    n = _ALIASES.get(n, n)
    if n not in _FROST:
        raise ValueError(f"Unknown buffer '{name}'. Choose from {available_buffers()}.")
    return n


def _pressure_term(name: str, T_K: float, P_Pa: float) -> float:
    _, _, C = _FROST[name]
    P_bar = P_Pa / BAR_TO_PA
    return C * (P_bar - 1.0) / T_K


def log_fo2_buffer(name: str, T_K: float, P_Pa: float = 1.0e5, parameterisation: str = "frost1991") -> float:
    """log10 fO2 (bar) of a solid buffer at T and P.

    Source: Frost (1991) Table 1 [frost1991]; or O'Neill & Pownceby (1993)
    [oneill_pownceby1993] for NNO and O'Neill (1987) [oneill1987] for FMQ
    when ``parameterisation='oneill'``.
    """
    n = _canon(name)
    if parameterisation == "frost1991":
        A, B, _ = _FROST[n]
        return A / T_K + B + _pressure_term(n, T_K, P_Pa)
    if parameterisation == "oneill":
        if n == "NNO":
            dG = -478967.0 + 248.514 * T_K - 9.7961 * T_K * math.log(T_K)
        elif n == "FMQ":
            dG = -587474.0 + 1584.427 * T_K - 203.3164 * T_K * math.log(T_K) + 0.092710 * T_K ** 2
        else:
            raise ValueError("O'Neill parameterisation only available for NNO and FMQ")
        return dG / (R_GAS * T_K * LN10) + _pressure_term(n, T_K, P_Pa)
    raise ValueError(f"Unknown parameterisation '{parameterisation}'")


def log_fo2_from_delta(buffer: str, delta: float, T_K: float, P_Pa: float = 1.0e5,
                       parameterisation: str = "frost1991") -> float:
    """Absolute log10 fO2 from an offset ``delta`` (log units) relative to a buffer."""
    return log_fo2_buffer(buffer, T_K, P_Pa, parameterisation) + delta


def delta_from_log_fo2(buffer: str, log_fo2: float, T_K: float, P_Pa: float = 1.0e5,
                       parameterisation: str = "frost1991") -> float:
    return log_fo2 - log_fo2_buffer(buffer, T_K, P_Pa, parameterisation)
