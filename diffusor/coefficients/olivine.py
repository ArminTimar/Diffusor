"""Diffusion coefficients for olivine.

Verification status
--------------------------------
Dohmen & Chakraborty (2007) equations 27 (TaMED) and 28 (PED) were checked on
29 September 2026 against the printed paper (Phys. Chem. Minerals 34:409-430,
pp. 424-425) and its erratum (34:597-598, doi:10.1007/s00269-007-0185-3).
The printed equations end in ``+ 3 X_Fe``; the erratum corrects that to
``+ 3 (X_Fe - 0.1)``, consistent with the data being normalised to Fo90. Every
other constant (log D0, Q, dV, the 1/6 fO2 exponent, fO2 reference 1e-7 Pa,
the 1e-10 Pa regime boundary, log 6 for [100] and [010]) is unchanged by the
erratum and matches the code.

Between 20 and 29 September 2026 Diffusor used the printed ``+ 3 X_Fe``, which
makes D 10**0.3 (about 2x) too large and olivine Fe-Mg times about half too
short. Refit olivine results made in that window.
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
# Dohmen & Chakraborty (2007), eqs 27-28 as corrected by the erratum
# ---------------------------------------------------------------------------
def _dohmen_chakraborty_tamed(dc, cond: Conditions, p):
    """TaMED regime (transition-metal extrinsic diffusion), D//[001]:

    log D = -9.21 - [201000 + (P - 1e5) * 7e-6] / (2.303 R T)
            + 1/6 log(fO2 / 1e-7) + 3 (X_Fe - 0.1)

    fO2 in Pa, P in Pa, D in m^2/s. Eq. 27 with the erratum's composition term.
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

    Eq. 28 with the erratum's composition term.
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
    kind="interdiffusion", transported_variable="XFe = Fe/(Fe+Mg), mole fraction",
    calibration_notes=("Eq. 27 applies at fO2 above 1e-10 Pa (TaMED); below it use the PED "
                       "entry (eq. 28).",),
    label="Olivine Fe-Mg // [001], Dohmen & Chakraborty (2007), TaMED regime",
    citation="dohmen_chakraborty2007",
    equation_number="27 (TaMED), p. 424, as corrected by the erratum (PCM 34:597-598)",
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
        "XFe_ref": Parameter("XFe_ref", 0.1, 0.0, "-", "1s",
                             "reference X_Fe (Fo90), from the erratum"),
    },
    sigma_logD=0.21,
    axis_factors=_OLIVINE_ANISO,
    reference_axis="c",
    requires=("XFe",),
    needs_fo2=True, fo2_unit="Pa",
    T_range=Range(973.15, 1473.15, "K (700-1200 C)"),
    P_range=Range(1.0e5, 1.2e10, "Pa"),
    fo2_range=Range(-10.0, None, "log10 Pa: TaMED regime above 1e-10 Pa (p. 424), no upper limit given"),
    X_range=Range(0.0, 0.5, "XFe"),
    verified=True,
    verified_from=("Checked 29 September 2026 against the primary paper, Phys. Chem. Minerals "
                   "34:409-430, pp. 424-425 (eqs 27, 28), and its erratum, 34:597-598 "
                   "(doi:10.1007/s00269-007-0185-3), which corrects the composition term of "
                   "eqs 27 and 28 from 3 X_Fe to 3 (X_Fe - 0.1). "
                   "Also confirmed from the local library: the fO2 exponent of 1/6 and "
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
           "scatter of the experimental database; the paper gives no parameter covariance, "
           "only that eqs 27-28 reproduce all 113 data points within 0.5 log units."),
))

_add(DiffusionCoefficient(
    key="ol_FeMg_dohmen_chakraborty2007_ped",
    mineral="olivine", species="Fe-Mg",
    kind="interdiffusion", transported_variable="XFe = Fe/(Fe+Mg), mole fraction",
    label="Olivine Fe-Mg // [001], Dohmen & Chakraborty (2007), PED regime",
    citation="dohmen_chakraborty2007",
    equation_number="28 (PED), p. 425, as corrected by the erratum (PCM 34:597-598)",
    equation_text=("log D_Fe-Mg [m2/s] = -8.91 - (220000 + (P - 1e5) * 7e-6) / (2.303 R T) "
                   "+ 3 (X_Fe - 0.1),  P in Pa. Independent of fO2"),
    func=_dohmen_chakraborty_ped,
    params={
        "c0": Parameter("c0", -8.91, 0.0, "log10(m2/s)", "1s", "intercept"),
        "Q": Parameter("Q", 220.0, 0.0, "kJ/mol", "1s", "activation energy (PED)"),
        "dV": Parameter("dV", 7.0e-6, 0.0, "m3/mol", "1s", "activation volume"),
        "m": Parameter("m", 3.0, 0.0, "-", "1s", "X_Fe coefficient"),
        "XFe_ref": Parameter("XFe_ref", 0.1, 0.0, "-", "1s", "reference X_Fe (Fo90), from the erratum"),
    },
    sigma_logD=0.21,
    axis_factors=_OLIVINE_ANISO,
    reference_axis="c",
    requires=("XFe",),
    needs_fo2=False,
    T_range=Range(973.15, 1473.15, "K (700-1200 C)"),
    X_range=Range(0.0, 0.5, "XFe"),
    verified=True,
    verified_from=("Checked 29 September 2026 against the primary paper, Phys. Chem. Minerals "
                   "34:409-430, pp. 424-425 (eqs 27, 28), and its erratum, 34:597-598 "
                   "(doi:10.1007/s00269-007-0185-3), which corrects the composition term of "
                   "eqs 27 and 28 from 3 X_Fe to 3 (X_Fe - 0.1). Eq. 28 applies at fO2 below 1e-10 Pa."),
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
    equation_number="Arrhenius fit at Fo86, // [001], fO2 = 1e-12 bar (abstract p. 12,317 and p. 12,325)",
    equation_text=("D = D0 exp(-Q / R T) with D0 = (5.38 +/- 0.89) x 10^-9 m2/s and "
                   "Q = 226 +/- 18.5 kJ/mol (54 +/- 4.4 kcal/mol), for Fo86 // [001] at "
                   "fO2 = 1e-12 bar (1e-7 Pa)"),
    func=_chakraborty1997,
    params={
        "D0": Parameter("D0", 5.38e-9, 0.89e-9, "m2/s", "1s", "pre-exponential factor"),
        "Q": Parameter("Q", 226.0, 18.5, "kJ/mol", "1s", "activation energy"),
    },
    sigma_logD=0.3,
    reference_axis="c",
    allowed_axes=("c",),
    reference_state="Fo86, // [001], fO2 = 1e-12 bar, 1 atm",
    needs_fo2=False,
    fo2_unit="Pa",
    fo2_range=Range(-7.0, -7.0, "log10 Pa (fO2 = 1e-12 bar)"),
    T_range=Range(1253.15, 1573.15, "K (980-1300 C)"),
    X_range=Range(0.12, 0.16, "XFe (around Fo86)"),
    verified=True,
    verified_from=("read from the paper PDF: abstract (p. 12,317) and p. 12,325. The abstract "
                   "rounds the Q uncertainty to 18 kJ/mol. Before 1 October 2026 this entry used "
                   "D0 = 1.0e-9 m2/s from a secondary summary, which made D 5.4 times too low"),
    notes=("Superseded by Dohmen & Chakraborty (2007) but retained so that older published "
           "timescales can be reproduced and compared. The fit holds for one composition "
           "(Fo86), one direction ([001]) and one oxygen fugacity (1e-12 bar). The paper "
           "measured no other direction, so no anisotropy is applied and other axes are "
           "refused. It gives no fO2 term either, so other fO2 values raise a range warning. "
           "For Fo92 the paper gives D0 = 6.59e-9 m2/s and Q = 229 +/- 18 kJ/mol. sigma_logD "
           "is Diffusor's judgement, not a published scatter."),
))
