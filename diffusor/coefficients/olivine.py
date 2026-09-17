"""Diffusion coefficients for olivine.

IMPORTANT -- verification status
--------------------------------
The Dohmen & Chakraborty (2007) law below is the standard formulation used by
essentially every olivine diffusion study, but the primary PDF (and its
erratum) was not available offline when this module was written, so the
numerical coefficients were transcribed from secondary sources and are marked
``verified=False``.  The parts that *are* independently confirmed from sources
in the local library are noted per entry.  Check the transcription against
Dohmen & Chakraborty (2007) Phys Chem Minerals 34:409-430 together with the
erratum (34:597-598) before using these numbers in a publication.
"""
from __future__ import annotations

import numpy as np

from ..constants import R_GAS
from .base import Conditions, DiffusionCoefficient, Parameter, Range, LN10, arrhenius

COEFFICIENTS = []


def _add(c):
    COEFFICIENTS.append(c)
    return c


# ---------------------------------------------------------------------------
# Dohmen & Chakraborty (2007), corrected by the 2007 erratum
# ---------------------------------------------------------------------------
def _dohmen_chakraborty_tamed(dc, cond: Conditions, p):
    """TaMED regime (transition-metal extrinsic diffusion), D//[001]:

    log D = -9.21 - [201000 + (P - 1e5) * 7e-6] / (2.303 R T)
            + 1/6 log(fO2 / 1e-7) + 3 (X_Fe - 0.1)

    fO2 in Pa, P in Pa, D in m^2/s.
    """
    XFe = np.asarray(cond.x("XFe"), dtype=float)
    fo2_Pa_log = cond.log_fo2_Pa
    logD = (p["c0"]
            - (p["Q"] * 1.0e3 + (cond.P_Pa - 1.0e5) * p["dV"]) / (LN10 * R_GAS * cond.T_K)
            + p["n_fo2"] * (fo2_Pa_log - p["log_fo2_ref"])
            + p["m"] * (XFe - p["XFe_ref"]))
    return 10.0 ** logD


def _dohmen_chakraborty_ped(dc, cond: Conditions, p):
    """PED regime (pure extrinsic diffusion, fO2-independent), D//[001]:

    log D = -8.91 - [220000 + (P - 1e5) * 7e-6] / (2.303 R T) + 3 (X_Fe - 0.1)
    """
    XFe = np.asarray(cond.x("XFe"), dtype=float)
    logD = (p["c0"]
            - (p["Q"] * 1.0e3 + (cond.P_Pa - 1.0e5) * p["dV"]) / (LN10 * R_GAS * cond.T_K)
            + p["m"] * (XFe - p["XFe_ref"]))
    return 10.0 ** logD


_OLIVINE_ANISO = {"c": 1.0, "b": 1.0 / 6.0, "a": 1.0 / 6.0}

_add(DiffusionCoefficient(
    key="ol_FeMg_dohmen_chakraborty2007_tamed",
    mineral="olivine", species="Fe-Mg",
    label="Olivine Fe-Mg // [001], Dohmen & Chakraborty (2007), TaMED regime",
    citation="dohmen_chakraborty2007",
    equation_number="model equation of Dohmen & Chakraborty (2007) as corrected by the erratum",
    equation_text=("log D_Fe-Mg [m2/s] = -9.21 - (201000 + (P - 1e5) * 7e-6) / (2.303 R T) "
                   "+ 1/6 log(fO2 / 1e-7) + 3 (X_Fe - 0.1),  fO2 and P in Pa"),
    func=_dohmen_chakraborty_tamed,
    params={
        "c0": Parameter("c0", -9.21, 0.0, "log10(m2/s)", "1s", "intercept at the reference state"),
        "Q": Parameter("Q", 201.0, 0.0, "kJ/mol", "1s", "activation energy (TaMED)"),
        "dV": Parameter("dV", 7.0e-6, 0.0, "m3/mol", "1s", "activation volume (Holzapfel et al. 2007)"),
        "n_fo2": Parameter("n_fo2", 1.0 / 6.0, 0.0, "-", "1s", "fO2 exponent"),
        "log_fo2_ref": Parameter("log_fo2_ref", -7.0, 0.0, "log10 Pa", "1s", "reference fO2"),
        "m": Parameter("m", 3.0, 0.0, "-", "1s", "X_Fe coefficient"),
        "XFe_ref": Parameter("XFe_ref", 0.1, 0.0, "-", "1s", "reference X_Fe of the formulation"),
    },
    sigma_logD=0.21,
    axis_factors=_OLIVINE_ANISO,
    reference_axis="c",
    requires=("XFe",),
    needs_fo2=True, fo2_unit="Pa",
    T_range=Range(973.15, 1473.15, "K (700-1200 C)"),
    P_range=Range(1.0e5, 1.2e10, "Pa"),
    fo2_range=Range(-12.0, -5.0, "log10 Pa"),
    X_range=Range(0.0, 0.5, "XFe"),
    verified=False,
    verified_from=("coefficients transcribed from secondary sources, NOT from the primary PDF. "
                   "Independently confirmed from the local library: the fO2 exponent of 1/6 and "
                   "the TaMED/PED mechanism change (Dohmen et al. 2016, p. 2216 and p. 2219). "
                   "The ~6x anisotropy of D//[001] over D//[100] and D//[010] (Hartley et al. "
                   "2016, p. 60. Mutch et al. 2021 DFENS source code uses aniso = 6.0). The "
                   "activation volume of 7e-6 m3/mol (Costa et al. 2008, p. 571, citing "
                   "Holzapfel et al. 2007)."),
    secondary_citations=("dohmen2007", "holzapfel2007", "chakraborty2010", "hartley2016"),
    recommended=True,
    notes=("Valid where transition-metal extrinsic diffusion dominates, i.e. at relatively "
           "oxidising conditions. Below about the IW buffer the mechanism changes to PED and "
           "D stops depending on fO2 -- use the PED entry there. D//[001] is about 6 times "
           "D//[100] and D//[010]. The sigma_logD of 0.21 reflects the roughly 0.2-log-unit "
           "scatter of the experimental database and should be replaced with the published "
           "covariance once the primary paper is checked."),
))

_add(DiffusionCoefficient(
    key="ol_FeMg_dohmen_chakraborty2007_ped",
    mineral="olivine", species="Fe-Mg",
    label="Olivine Fe-Mg // [001], Dohmen & Chakraborty (2007), PED regime",
    citation="dohmen_chakraborty2007",
    equation_number="model equation of Dohmen & Chakraborty (2007) as corrected by the erratum",
    equation_text=("log D_Fe-Mg [m2/s] = -8.91 - (220000 + (P - 1e5) * 7e-6) / (2.303 R T) "
                   "+ 3 (X_Fe - 0.1),  P in Pa. Independent of fO2"),
    func=_dohmen_chakraborty_ped,
    params={
        "c0": Parameter("c0", -8.91, 0.0, "log10(m2/s)", "1s", "intercept"),
        "Q": Parameter("Q", 220.0, 0.0, "kJ/mol", "1s", "activation energy (PED)"),
        "dV": Parameter("dV", 7.0e-6, 0.0, "m3/mol", "1s", "activation volume"),
        "m": Parameter("m", 3.0, 0.0, "-", "1s", "X_Fe coefficient"),
        "XFe_ref": Parameter("XFe_ref", 0.1, 0.0, "-", "1s", "reference X_Fe"),
    },
    sigma_logD=0.21,
    axis_factors=_OLIVINE_ANISO,
    reference_axis="c",
    requires=("XFe",),
    needs_fo2=False,
    T_range=Range(973.15, 1473.15, "K (700-1200 C)"),
    X_range=Range(0.0, 0.5, "XFe"),
    verified=False,
    verified_from="as for the TaMED entry -- transcribed from secondary sources",
    notes=("Pure extrinsic diffusion: point defect concentrations are fixed by aliovalent "
           "impurities, so D does not depend on fO2. Applies at reducing "
           "conditions (around and below the IW buffer). Diffusor does not switch between "
           "the two regimes automatically -- choose the one appropriate to your fO2, or run "
           "both and compare."),
))


# ---------------------------------------------------------------------------
# Older comparison calibration
# ---------------------------------------------------------------------------
def _chakraborty1997(dc, cond: Conditions, p):
    return arrhenius(p["D0"], p["Q"] * 1.0e3, cond.T_K)


_add(DiffusionCoefficient(
    key="ol_FeMg_chakraborty1997",
    mineral="olivine", species="Fe-Mg",
    label="Olivine Fe-Mg // [001], Chakraborty (1997) -- older calibration",
    citation="chakraborty1997",
    equation_number="(Arrhenius fit at fO2 = 1e-7 Pa, Fo86)",
    equation_text="D = D0 exp(-Q / R T) m2/s with D0 = 1.0e-9 m2/s and Q = 226 kJ/mol",
    func=_chakraborty1997,
    params={
        "D0": Parameter("D0", 1.0e-9, 0.0, "m2/s", "1s", "pre-exponential factor"),
        "Q": Parameter("Q", 226.0, 0.0, "kJ/mol", "1s", "activation energy"),
    },
    sigma_logD=0.3,
    axis_factors=_OLIVINE_ANISO,
    reference_axis="c",
    needs_fo2=False,
    T_range=Range(1253.15, 1573.15, "K (980-1300 C)"),
    X_range=Range(0.12, 0.16, "XFe (around Fo86)"),
    verified=False,
    verified_from="approximate transcription from secondary summaries. NOT checked against the primary paper",
    notes=("Superseded by Dohmen & Chakraborty (2007) but retained so that older published "
           "timescales can be reproduced and compared. Treat the numbers as provisional."),
))
