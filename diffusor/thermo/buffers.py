"""Oxygen fugacity buffers.

All buffers are returned as ``log10 fO2`` (fO2 in bar) for a temperature in K
and pressure in Pa.

Two parameterisations are offered:

``frost1991`` (default)
    Compilation form ``log fO2 = A/T + B + C (P - 1)/T`` with *P* in bar,
    Frost (1991), Reviews in Mineralogy 25, Table 1 (citekey ``frost1991``).
    The table prints one expression per temperature range; every row used
    here is listed with its printed range in degrees C:

    ===== ========== ======= ===== ========= ============================
    IW    -27489     6.702   0.055 565-1200
    WM    -32807     13.012  0.083 565-1200
    FMQ   -26455.3   10.344  0.092 400-573   alpha-quartz row (FMaQ)
    FMQ   -25096.3   8.735   0.110 573-1200  beta-quartz row (FMbQ)
    NNO   -24930     9.36    0.046 600-1200  (NiNiO)
    HM    -25497.5   14.330  0.019 300-573   (MH)
    HM    -26452.6   15.455  0.019 573-682
    HM    -25700.6   14.558  0.019 682-1100
    QIF   -29435.7   7.391   0.044 150-573   alpha-quartz row (aQIF)
    QIF   -29520.8   7.492   0.050 573-1200  beta-quartz row (bQIF)
    ===== ========== ======= ===== ========= ============================

    Quartz-bearing buffers (FMQ, QIF) switch between the alpha- and beta-quartz
    row at the alpha-beta transition, T(C) = 573 + 0.025 P(bar), the
    approximation printed under the table.  The HM rows change at 573 and
    682 C.  Table 1 states that the expressions "may be strongly curved at low
    pressures" and advises against extrapolating below the limits indicated.
    When the requested temperature is outside the printed range of the buffer
    (IW and WM below 565 C, where the table stops them and lists the
    iron-magnetite buffer IM for 300-565 C instead, which Diffusor does not
    implement; above 1200 C, or 1100 C for HM; below 400 C for FMQ,
    300 C for HM, 150 C for QIF, 600 C for NNO),
    :class:`BufferRangeWarning` is issued and the nearest printed row is
    extrapolated.  :func:`buffer_range_status` tests a temperature without
    computing anything.  The table prints no pressure limit.

    Two features of the printed table are kept as printed: the 573-682 C and
    682-1100 C HM expressions differ by 0.11 log units at 682 C, and the alpha-
    and beta-quartz rows meet within 0.003 log units at 573 C and 1 bar.

``oneill``
    Free-energy expressions of the primary calibrations, 1 bar:

    * NNO: O'Neill & Pownceby (1993), CMP 114, eq. for 2 Ni + O2 = 2 NiO,
      ``dG(J/mol) = -478967 + 248.514 T - 9.7961 T ln T``, printed for
      700 < T < 1700 K (citekey ``oneill_pownceby1993``), so
      ``log fO2 = dG / (R T ln 10)``.
    * FMQ: O'Neill (1987), Am Mineral 72, eq. for 3 Fe2SiO4 + O2 = 2 Fe3O4 + 3 SiO2,
      ``dG(J/mol) = -587474 + 1584.427 T - 203.3164 T ln T + 0.092710 T^2``,
      printed for 900 < T < 1420 K (citekey ``oneill1987``).
    The Frost (1991) pressure term ``C (P-1)/T`` is added to these so that both
    parameterisations agree on the pressure dependence.  A temperature outside
    the printed range of the O'Neill expression issues
    :class:`BufferRangeWarning` and the expression is extrapolated.

Note on correlation: in Monte Carlo runs the user supplies fO2 as an offset
from a buffer (e.g. NNO+1).  The buffer is re-evaluated at every sampled
temperature, so fO2 and T are correlated by construction; this is the
dependence that NIDIS (Petrone et al. 2016) ignores when it propagates errors
as if the variables were independent.
"""
from __future__ import annotations

import math
import warnings
from typing import Dict, List, Tuple

import numpy as np

from ..constants import R_GAS, BAR_TO_PA, T_KELVIN_OFFSET

LN10 = math.log(10.0)


class BufferRangeWarning(UserWarning):
    """A buffer was evaluated outside the temperature range its source prints.

    The value is still returned (the nearest printed expression is
    extrapolated), so existing callers keep working.  Callers that must not
    extrapolate can turn the warning into an error with
    ``warnings.simplefilter("error", BufferRangeWarning)`` or test the
    temperature first with :func:`buffer_range_status`.
    """


# Frost (1991) Table 1 rows in order of temperature: (A, B, C, T_low, T_high), T in degrees C.
_ROWS: Dict[str, List[Tuple[float, float, float, float, float]]] = {
    "IW":  [(-27489.0, 6.702, 0.055, 565.0, 1200.0)],
    "WM":  [(-32807.0, 13.012, 0.083, 565.0, 1200.0)],
    "FMQ": [(-26455.3, 10.344, 0.092, 400.0, 573.0),      # FMaQ, alpha-quartz
            (-25096.3, 8.735, 0.110, 573.0, 1200.0)],     # FMbQ, beta-quartz
    "NNO": [(-24930.0, 9.36, 0.046, 600.0, 1200.0)],      # NiNiO
    "HM":  [(-25497.5, 14.330, 0.019, 300.0, 573.0),      # MH
            (-26452.6, 15.455, 0.019, 573.0, 682.0),
            (-25700.6, 14.558, 0.019, 682.0, 1100.0)],
    "QIF": [(-29435.7, 7.391, 0.044, 150.0, 573.0),       # aQIF, alpha-quartz
            (-29520.8, 7.492, 0.050, 573.0, 1200.0)],     # bQIF, beta-quartz
}
# Buffers whose two rows are the alpha- and beta-quartz branches; they switch at the
# alpha-beta transition, T(C) = 573 + 0.025 P(bar) (Frost 1991, Table 1 note).
_QUARTZ_BUFFERS = ("FMQ", "QIF")
# "IQF" (iron-quartz-fayalite) names the equilibrium Fe2SiO4 = 2 Fe + SiO2 + O2, which
# Table 1 lists as QIF.  It is not the iron-wuestite buffer (IW), which lies 0.7-1.1 log
# units higher between 800 and 1100 C.
_ALIASES = {"QFM": "FMQ", "IQF": "QIF"}

# Printed validity ranges of the O'Neill expressions, K.
_ONEILL_RANGE_K = {"NNO": (700.0, 1700.0), "FMQ": (900.0, 1420.0)}
_ONEILL_SOURCE = {"NNO": "O'Neill & Pownceby (1993)", "FMQ": "O'Neill (1987)"}

BUFFER_CITATIONS = {
    "frost1991": "frost1991",
    "NNO_oneill": "oneill_pownceby1993",
    "FMQ_oneill": "oneill1987",
}


def available_buffers():
    return list(_ROWS) + ["QFM"]


def _canon(name: str) -> str:
    n = name.upper()
    n = _ALIASES.get(n, n)
    if n not in _ROWS:
        raise ValueError(f"Unknown buffer '{name}'. Choose from {available_buffers()}.")
    return n


def _row_index(name: str, T_C, P_bar):
    """Index of the printed row used at T (degrees C) and P (bar); arrays allowed."""
    rows = _ROWS[name]
    T_C = np.asarray(T_C, dtype=float)
    if len(rows) == 1:
        return np.zeros(T_C.shape, dtype=int)
    if name in _QUARTZ_BUFFERS:
        return (T_C >= 573.0 + 0.025 * np.asarray(P_bar, dtype=float)).astype(int)
    idx = np.zeros(T_C.shape, dtype=int)
    for k in range(len(rows) - 1):          # a row ends where the next one starts
        idx = idx + (T_C >= rows[k][4]).astype(int)
    return idx


def _frost_coefficients(name: str, T_K, P_Pa):
    T_K = np.asarray(T_K, dtype=float)
    idx = _row_index(name, T_K - T_KELVIN_OFFSET, np.asarray(P_Pa, dtype=float) / BAR_TO_PA)
    table = np.array([r[:3] for r in _ROWS[name]])
    return table[idx, 0], table[idx, 1], table[idx, 2]


def buffer_range_status(name: str, T_K, parameterisation: str = "frost1991") -> Tuple[bool, str]:
    """``(inside, message)`` for a buffer at a temperature in K (scalar or array).

    ``inside`` is True when every temperature lies within the range the source
    prints for the buffer; ``message`` names that range and is the text of the
    :class:`BufferRangeWarning` issued by :func:`log_fo2_buffer` when it is False.
    """
    n = _canon(name)
    T_K = np.asarray(T_K, dtype=float)
    if parameterisation == "oneill":
        if n not in _ONEILL_RANGE_K:
            raise ValueError("O'Neill parameterisation only available for NNO and FMQ")
        lo, hi = _ONEILL_RANGE_K[n]
        inside = bool(np.all((T_K >= lo) & (T_K <= hi)))
        return inside, (f"{n} from {_ONEILL_SOURCE[n]} is printed for {lo:.0f}-{hi:.0f} K. A requested "
                        f"temperature lies outside, so the expression is extrapolated")
    if parameterisation != "frost1991":
        raise ValueError(f"Unknown parameterisation '{parameterisation}'")
    lo, hi = _ROWS[n][0][3], _ROWS[n][-1][4]
    T_C = T_K - T_KELVIN_OFFSET
    inside = bool(np.all((T_C >= lo) & (T_C <= hi)))
    return inside, (f"{n} from Frost (1991) Table 1 is printed for {lo:.0f}-{hi:.0f} C. A requested "
                    f"temperature lies outside, so the nearest printed expression is extrapolated "
                    f"(the table advises against extrapolating below its limits)")


def _pressure_term(C, T_K, P_Pa):
    P_bar = np.asarray(P_Pa, dtype=float) / BAR_TO_PA
    return C * (P_bar - 1.0) / T_K


def _scalar_if_scalar(x):
    x = np.asarray(x, dtype=float)
    return float(x) if x.ndim == 0 else x


def log_fo2_buffer(name: str, T_K: float, P_Pa: float = 1.0e5, parameterisation: str = "frost1991") -> float:
    """log10 fO2 (bar) of a solid buffer at T and P.

    Source: Frost (1991) Table 1 [frost1991]; or O'Neill & Pownceby (1993)
    [oneill_pownceby1993] for NNO and O'Neill (1987) [oneill1987] for FMQ
    when ``parameterisation='oneill'``.  A temperature outside the range the
    source prints issues :class:`BufferRangeWarning`; the value is still
    returned by extrapolation.  ``T_K`` may be an array.
    """
    n = _canon(name)
    if parameterisation not in ("frost1991", "oneill"):
        raise ValueError(f"Unknown parameterisation '{parameterisation}'")
    if parameterisation == "oneill" and n not in _ONEILL_RANGE_K:
        raise ValueError("O'Neill parameterisation only available for NNO and FMQ")
    inside, message = buffer_range_status(n, T_K, parameterisation)
    if not inside:
        warnings.warn(message, BufferRangeWarning, stacklevel=2)
    T_K = np.asarray(T_K, dtype=float)
    A, B, C = _frost_coefficients(n, T_K, P_Pa)
    if parameterisation == "frost1991":
        return _scalar_if_scalar(A / T_K + B + _pressure_term(C, T_K, P_Pa))
    if n == "NNO":
        dG = -478967.0 + 248.514 * T_K - 9.7961 * T_K * np.log(T_K)
    else:
        dG = -587474.0 + 1584.427 * T_K - 203.3164 * T_K * np.log(T_K) + 0.092710 * T_K ** 2
    return _scalar_if_scalar(dG / (R_GAS * T_K * LN10) + _pressure_term(C, T_K, P_Pa))


def log_fo2_from_delta(buffer: str, delta: float, T_K: float, P_Pa: float = 1.0e5,
                       parameterisation: str = "frost1991") -> float:
    """Absolute log10 fO2 from an offset ``delta`` (log units) relative to a buffer."""
    return log_fo2_buffer(buffer, T_K, P_Pa, parameterisation) + delta


def delta_from_log_fo2(buffer: str, log_fo2: float, T_K: float, P_Pa: float = 1.0e5,
                       parameterisation: str = "frost1991") -> float:
    return log_fo2 - log_fo2_buffer(buffer, T_K, P_Pa, parameterisation)
