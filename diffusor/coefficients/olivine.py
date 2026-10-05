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

What the paper does and does not state for eqs 27-28 (so that the ranges below
are not mistaken for printed validity limits):

* The paper prints no temperature, pressure or composition limits for the
  equations. They are fitted to 113 data points (largest deviation 0.5 log
  units, normalised to 1e5 Pa, fO2 = 1e-7 Pa, Fo90; p. 425). Diffusor's T range
  (700-1200 C) is the range of the Part I experiments (Dohmen et al. 2007,
  abstract); the point-defect model is stated for 600-1300 C (p. 3). The P range
  up to 12 GPa is the range of the high-pressure data of Holzapfel et al. (2007;
  6-12 GPa) that is included in the comparison of Fig. 7b. The X_Fe range 0-0.5
  is Diffusor's choice.
* The regime boundary is fO2 = 1e-10 Pa (and "consequently temperatures below
  900 C" for buffered conditions); the paper does not tie it to the IW buffer.
* The paper's single compromise equation (eq. 29; abstract eq. 3), which in the
  erratum reads log D = -8.27 - [226000 + (P - 1e5) 7e-6]/(2.303 R T)
  + 3 (X_Fe - 0.14), is not implemented. Use the TaMED or PED law.
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
    calibration_notes=("Eq. 27 applies at fO2 above 1e-10 Pa (TaMED). At or below it use the PED "
                       "entry (eq. 28).",
                       "The paper prints no temperature, pressure or composition limits for eq. 27. The "
                       "700-1200 C range is that of the Part I experiments (Dohmen et al. 2007), the "
                       "pressure range reaches the 12 GPa data of Holzapfel et al. (2007) used in Fig. 7b, "
                       "and the X_Fe range of 0 to 0.5 is Diffusor's choice (the m = 3 term derives from "
                       "data in Chakraborty 1997, p. 424).",
                       "The paper's single compromise equation (eq. 29, erratum: -8.27, 226 kJ/mol, "
                       "3 (X_Fe - 0.14)) is not implemented."),
    label="Olivine Fe-Mg // [001], Dohmen & Chakraborty (2007), TaMED regime",
    citation="dohmen_chakraborty2007",
    equation_number="27 (TaMED), p. 424, as corrected by the erratum (PCM 34:597-598)",
    equation_text=("log D_Fe-Mg [m2/s] = -9.21 - (201000 + (P - 1e5) * 7e-6) / (2.303 R T) "
                   "+ 1/6 log(fO2 / 1e-7) + 3 (X_Fe - 0.1),  fO2 and P in Pa"),
    func=_dohmen_chakraborty_tamed,
    params={
        "c0": Parameter("c0", -9.21, 0.0, "log10(m2/s)", "1s", "intercept at the reference state"),
        "Q": Parameter("Q", 201.0, 0.0, "kJ/mol", "1s", "activation energy (TaMED)"),
        "dV": Parameter("dV", 7.0e-6, 0.0, "m3/mol", "1s",
                        "activation volume at constant fO2 (Holzapfel et al. 2007: ~7 cm3/mol)"),
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
    T_range=Range(973.15, 1473.15, "K (700-1200 C, range of the Part I experiments)"),
    P_range=Range(1.0e5, 1.2e10, "Pa (1 atm to 12 GPa, Holzapfel et al. 2007 data)"),
    fo2_range=Range(-10.0, None, "log10 Pa: TaMED regime above 1e-10 Pa (p. 424), no upper limit given"),
    X_range=Range(0.0, 0.5, "XFe (Diffusor's choice, not stated in the paper)"),
    verified=True,
    verified_from=("Checked 29 September 2026 against the primary paper, Phys. Chem. Minerals "
                   "34:409-430, pp. 424-425 (eqs 27, 28), and its erratum, 34:597-598 "
                   "(doi:10.1007/s00269-007-0185-3), which corrects the composition term of "
                   "eqs 27 and 28 from 3 X_Fe to 3 (X_Fe - 0.1). "
                   "The fO2 exponent of 1/6 (n ~ 6), the 1e-10 Pa boundary between the TaMED and "
                   "PED regimes and the subtraction of log 6 for [100] and [010] are printed in the "
                   "paper itself (abstract, pp. 424-425). "
                   "The ~6x anisotropy of D//[001] over D//[100] and D//[010] (Hartley et al. "
                   "2016, p. 61. Mutch et al. 2021 DFENS source code uses aniso = 6.0). The "
                   "activation volume of 7e-6 m3/mol (Costa et al. 2008, p. 571, citing "
                   "Holzapfel et al. 2007)."),
    secondary_citations=("dohmen2007", "holzapfel2007", "chakraborty2010", "hartley2016"),
    recommended=True,
    notes=("Valid where transition-metal extrinsic diffusion dominates, i.e. at fO2 above 1e-10 Pa. "
           "At or below 1e-10 Pa (for buffered conditions the paper says 'consequently, temperatures "
           "below 900 C') the mechanism changes to PED and D stops depending on fO2 -- use the PED "
           "entry there. For orientation, the Frost (1991) IW buffer is at log fO2 = -11.7 Pa at 900 C, "
           "so 1e-10 Pa lies about 1.7 log units above IW. The paper does not name IW as the boundary. "
           "D//[001] is about 6 times D//[100] and D//[010]. The sigma_logD of 0.21 is Diffusor's "
           "assumption, not a published scatter: the paper gives no parameter covariance and states "
           "only that eqs 27-28 reproduce all 113 data points within a maximum deviation of 0.5 log units."),
))

_add(DiffusionCoefficient(
    key="ol_FeMg_dohmen_chakraborty2007_ped",
    mineral="olivine", species="Fe-Mg",
    kind="interdiffusion", transported_variable="XFe = Fe/(Fe+Mg), mole fraction",
    label="Olivine Fe-Mg // [001], Dohmen & Chakraborty (2007), PED regime",
    citation="dohmen_chakraborty2007",
    equation_number="28 (PED), p. 425, as corrected by the erratum (PCM 34:597-598)",
    calibration_notes=("Eq. 28 was fitted to the data at fO2 equal to or below 1e-10 Pa (p. 425). The paper "
                       "ties the PED regime to fO2 < 1e-10 Pa and, for buffered conditions, to temperatures "
                       "below about 900 C, and it prints no temperature, pressure or composition limits for the "
                       "equation. The 700-1200 C range is that of the Part I experiments and 0 to 0.5 XFe is "
                       "Diffusor's choice. At 900-1200 C the PED law applies only if fO2 is at or below 1e-10 Pa. "
                       "The data behind eq. 28 are 1 atm experiments, so the pressure term (dV) is carried over "
                       "from eq. 27 and no pressure range is enforced.",
                       "The paper's single compromise equation (eq. 29) is not implemented."),
    equation_text=("log D_Fe-Mg [m2/s] = -8.91 - (220000 + (P - 1e5) * 7e-6) / (2.303 R T) "
                   "+ 3 (X_Fe - 0.1),  P in Pa. Independent of fO2"),
    func=_dohmen_chakraborty_ped,
    params={
        "c0": Parameter("c0", -8.91, 0.0, "log10(m2/s)", "1s", "intercept"),
        "Q": Parameter("Q", 220.0, 0.0, "kJ/mol", "1s", "activation energy (PED)"),
        "dV": Parameter("dV", 7.0e-6, 0.0, "m3/mol", "1s",
                        "activation volume at constant fO2 (Holzapfel et al. 2007: ~7 cm3/mol)"),
        "m": Parameter("m", 3.0, 0.0, "-", "1s", "X_Fe coefficient"),
        "XFe_ref": Parameter("XFe_ref", 0.1, 0.0, "-", "1s", "reference X_Fe (Fo90), from the erratum"),
    },
    sigma_logD=0.21,
    axis_factors=_OLIVINE_ANISO,
    reference_axis="c",
    requires=("XFe",),
    needs_fo2=False,
    T_range=Range(973.15, 1473.15, "K (700-1200 C, range of the Part I experiments)"),
    fo2_range=Range(None, -10.0, "log10 Pa: PED regime at or below 1e-10 Pa (p. 425)"),
    X_range=Range(0.0, 0.5, "XFe (Diffusor's choice, not stated in the paper)"),
    verified=True,
    verified_from=("Checked 29 September 2026 against the primary paper, Phys. Chem. Minerals "
                   "34:409-430, pp. 424-425 (eqs 27, 28), and its erratum, 34:597-598 "
                   "(doi:10.1007/s00269-007-0185-3), which corrects the composition term of "
                   "eqs 27 and 28 from 3 X_Fe to 3 (X_Fe - 0.1). Eq. 28 is fitted to data at fO2 equal to "
                   "or below 1e-10 Pa (p. 425)."),
    notes=("Pure extrinsic diffusion: point defect concentrations are fixed by aliovalent "
           "impurities, so D does not depend on fO2. Applies at fO2 at or below 1e-10 Pa "
           "('consequently, temperatures below 900 C' for buffered conditions, and the paper does not "
           "name the IW buffer). The sigma_logD of 0.21 is Diffusor's assumption, as for eq. 27. "
           "Diffusor does not switch between the two regimes automatically -- choose the one "
           "appropriate to your fO2, or run both and compare."),
))


# ---------------------------------------------------------------------------
# Older comparison calibration
# ---------------------------------------------------------------------------
def _chakraborty1997(dc, cond: Conditions, p):
    return arrhenius(p["D0"], p["Q"] * 1.0e3, cond.T_K)


_add(DiffusionCoefficient(
    key="ol_FeMg_chakraborty1997",
    kind="interdiffusion",
    mineral="olivine", species="Fe-Mg",
    label="Olivine Fe-Mg // [001], Chakraborty (1997) -- older calibration",
    citation="chakraborty1997",
    equation_number="Arrhenius fit at Fo86, // [001], fO2 = 1e-12 bar (abstract p. 12,317 and p. 12,325)",
    equation_text=("D = D0 exp(-Q / R T) with D0 = (5.38 +/- 0.89) x 10^-9 m2/s and "
                   "Q = 226 +/- 18.5 kJ/mol (54 +/- 4.4 kcal/mol), for Fo86 // [001] at "
                   "fO2 = 1e-12 bar (1e-7 Pa)"),
    func=_chakraborty1997,
    params={
        "D0": Parameter("D0", 5.38e-9, 0.89e-9, "m2/s", "unstated", "pre-exponential factor"),
        "Q": Parameter("Q", 226.0, 18.5, "kJ/mol", "unstated", "activation energy"),
    },
    sigma_logD=0.3,
    reference_axis="c",
    allowed_axes=("c",),
    reference_state=("Fo86, // [001], fO2 = 1e-12 bar, gas-mixing furnace (the pressure is not stated, "
                     "ambient pressure of about 1 atm is assumed)"),
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
           "reports no values for other directions (its few exploratory tests of orientation and fO2 "
           "were only said to be consistent with earlier studies, p. 12,325), so no anisotropy is "
           "applied and other axes are refused. It gives no fO2 term either, so other fO2 values raise "
           "a range warning. [001] is taken as the crystallographic c axis, which assumes the Pbnm "
           "setting that the paper does not state. "
           "For Fo92 the paper gives D0 = 6.59e-9 m2/s and Q = 229 +/- 18 kJ/mol. sigma_logD "
           "is Diffusor's judgement, not a published scatter."),
))


# ---------------------------------------------------------------------------
# Oeser, Dohmen & Weyer (2026): Fe and Mg tracer coefficients and the Fe-Mg
# interdiffusion coefficient computed from them
# ---------------------------------------------------------------------------
# Table 4 (p. 57): Arrhenius fits for 1100-1250 C at fO2 ~ 1e-5 Pa, uncertainties
# are 95 % confidence bounds. The fits are for X_Fe = 0.085 (Fig. 6 caption). The
# composition dependence D*(X) = D*(0.1) 10^(n (X - 0.1)), n = 3, is the one used in their
# multicomponent model (p. 51; given there for D*25 as an example, with reference
# composition X_Fe = 0.1). Diffusor re-references the same exponential law to
# X_Fe = 0.085, the composition of the Table 4 values: D(X) = D(0.085) 10^(3 (X - 0.085)),
# which equals D(0.1) 10^(3 (X - 0.1)) with D(0.1) = 1.109 D(0.085).
_OESER_TRACERS = {
    "Fe": {"a": (287.0, 52.0, -6.29, 1.88), "b": (307.0, 49.0, -5.22, 1.74), "c": (328.0, 36.0, -4.18, 1.32)},
    "Mg": {"a": (231.0, 34.0, -8.69, 1.24), "b": (292.0, 45.0, -6.11, 1.60), "c": (329.0, 75.0, -4.39, 2.74)},
}
_OESER_XFE_REF = 0.085


def _oeser_axis(axis):
    def evaluate(dc, cond: Conditions, p):
        XFe = np.asarray(cond.x("XFe"), dtype=float)
        logD = (p[f"logD0_{axis}"] - p[f"Ea_{axis}"] * 1.0e3 / (LN10 * R_GAS * cond.T_K)
                + p["n"] * (XFe - p["XFe_ref"]))
        return 10.0 ** logD
    return evaluate


_OESER_NOTE = ("Table 4 gives 95 % confidence bounds on Ea and log D0 without their covariance. "
               "Sampling the two independently would greatly overstate the spread of D, so the "
               "coefficient is held fixed in the Monte Carlo.")

for _el, _fits in _OESER_TRACERS.items():
    _funcs = {ax: _oeser_axis(ax) for ax in _fits}
    _params = {"n": Parameter("n", 3.0, 0.0, "", "1s", "composition exponent"),
               "XFe_ref": Parameter("XFe_ref", _OESER_XFE_REF, 0.0, "", "1s",
                                    "X_Fe of the Table 4 fits (Fig. 6 caption), the paper's "
                                    "model reference is X_Fe = 0.1")}
    for _ax, (_Ea, _sEa, _lD0, _slD0) in _fits.items():
        # bounds are kept for reporting with sigma = 0 so they are not sampled
        _params[f"Ea_{_ax}"] = Parameter(f"Ea_{_ax}", _Ea, 0.0, "kJ/mol", "2s",
                                         f"activation energy // {_ax} (95 % bound +/- {_sEa:g})")
        _params[f"logD0_{_ax}"] = Parameter(f"logD0_{_ax}", _lD0, 0.0, "log10(m2/s)", "2s",
                                            f"log10 D0 // {_ax} (95 % bound +/- {_slD0:g})")
    _add(DiffusionCoefficient(
        key=f"ol_{_el}_oeser2026_tracer",
        kind="tracer",
        mineral="olivine", species=_el,
        label=f"Olivine {_el} tracer, Oeser, Dohmen & Weyer (2026), independent a/b/c laws",
        citation="oeser2026",
        equation_number=("Table 4 (p. 57) with eq. 5 (p. 53), composition term p. 51, "
                         "re-referenced from X_Fe = 0.1 to 0.085"),
        equation_text=("log D*_" + _el + "//axis = log D0_axis - Ea_axis / (2.303 R T) + 3 (X_Fe - 0.085), "
                       "Ea in J/mol here (eq. 5 is printed with Ea in kJ/mol and 10000/T); "
                       + "; ".join(f"{ax}: Ea = {v[0]:g} +/- {v[1]:g} kJ/mol, log D0 = {v[2]:g} +/- {v[3]:g}"
                                   for ax, v in _fits.items())
                       + " (95 % bounds), fO2 ~ 1e-5 Pa"),
        func=_funcs["c"], principal_funcs=_funcs, reference_axis="c", orientation_required=True,
        params=_params, requires=("XFe",), needs_fo2=False, fo2_unit="Pa",
        fo2_range=Range(-5.0, -5.0, "log10 Pa (Table 4 fits at fO2 ~ 1e-5 Pa)"),
        T_range=Range(1373.15, 1523.15, "K (1100-1250 C)"),
        P_range=Range(1.0e5, 101325.0, "Pa (atmospheric; 1 bar to 1 atm)"),
        X_range=Range(0.085, 0.085, "XFe (San Carlos olivine, Fo91.5)"),
        transported_variable=f"{_el} isotope tracer (25Mg and 57Fe doped powder source)",
        reference_state="San Carlos olivine Fo91.5, fO2 ~ 1e-5 Pa, high silica activity (SiO2 + opx in the source)",
        verified=True,
        verified_from="read from the rendered Table 4 (p. 57) and the model description on p. 51 of the open-access paper",
        uncertainty_note=_OESER_NOTE,
        calibration_notes=(
            "Fitted for 1100-1250 C. The high activation energies (230-370 kJ/mol) and the change of "
            "anisotropy at about 1150 C (D//b = D//a below, D//b > D//a from 1150 C) suggest a change of "
            "diffusion mechanism around 1150 C (abstract). Do not extrapolate to magmatic temperatures "
            "below 1100 C.",
            "Measured at high silica activity (SiO2 and orthopyroxene in the powder source). Jollands et "
            "al. (2020) showed a strong silica-activity effect on Mg diffusion in pure forsterite. "
            "Earlier work suggests the same for Ni, Co, Cr and Zr in San Carlos olivine (p. 57).",
            "The printed eq. 5 writes the Arrhenius term as Ea/(ln10 R) x 10000/T with Ea in kJ/mol. "
            "Diffusor uses Ea x 1000/(ln10 R T), which reproduces the log D values of the paper's Table 2.",
            "The Table 4 laws are stated for X_Fe = 0.085, the initial crystal composition. The rim "
            "composition changes during each run (core minus rim up to about 0.28 in Fo, Table 2), so the "
            "profiles sample X_Fe above 0.085 and the composition term of the multicomponent model applies.",
        ),
        notes=("The isotopic profiles (delta 25Mg, 57Fe, 26Mg, 56Fe) were fitted together with the isotope "
               "fractionation factors beta_Fe and beta_Mg using a seven-isotope multicomponent model. The "
               "chemical (Fo) profile is largely a prediction of that model and the fits are visual "
               "(p. 51). In the Table 4 regressions Fe tracer diffusion is faster than Mg along all three "
               "axes between 1100 and 1250 C, but individual 1250 C runs along c have D*Fe/D*Mg below 1 "
               "(Table 2)."),
    ))


from .transport import interdiffusion_from_tracers  # noqa: E402  (needs the tracer laws above)

_add(interdiffusion_from_tracers(
    "ol_FeMg_oeser2026",
    tracer_a=COEFFICIENTS[-2], tracer_b=COEFFICIENTS[-1], species="Fe-Mg", comp_key="XFe",
    label="Olivine Fe-Mg from the Fe and Mg tracer laws, Oeser, Dohmen & Weyer (2026), a/b/c",
    citation="oeser2026", equation_number="eq. 6 with the Table 4 tracer laws",
    X_range=Range(0.085, 0.085, "XFe (San Carlos olivine, Fo91.5)"),
    verified_from=("eq. 6 (p. 54) applied to the Table 4 tracer laws. At X_Mg = 0.915 it reproduces "
                   "the paper's own D_Fe-Mg regression (Table 4) within 0.2 log units from 1100 to 1250 C"),
    calibration_notes=(
        "Fitted for 1100-1250 C at fO2 ~ 1e-5 Pa. The authors read the high activation energies "
        "(230-370 kJ/mol) and the change of anisotropy at 1150 C as a possible change of mechanism "
        "around 1150 C (abstract). For magmatic temperatures below 1100 C Diffusor recommends "
        "Dohmen & Chakraborty (2007).",),
    notes=("The interdiffusion coefficient that the tracer laws imply for an ideal Fe-Mg binary. "
           "The paper computes its D_Fe-Mg values the same way (p. 54). Use it to compare with the "
           "Dohmen & Chakraborty (2007) laws, or as the chemical part of an isotope model."),
))
