"""Fe-Mg and Mg diffusion coefficients for orthopyroxene."""
from __future__ import annotations

import numpy as np

from ..constants import R_GAS
from .base import Conditions, DiffusionCoefficient, Parameter, Range, LN10, arrhenius

COEFFICIENTS = []


def _add(c: DiffusionCoefficient):
    COEFFICIENTS.append(c)
    return c


# ---------------------------------------------------------------------------
# Dohmen, ter Heege, Becker & Chakraborty (2016) Am Mineral 101, 2210-2221
# ---------------------------------------------------------------------------
def _dohmen2016_fs9(dc, cond: Conditions, p):
    """log D_FeMg//[001] = log D0 - Q/(ln10 R T) + n log fO2(Pa)   [eq. 1]

    with the compositional correction recommended on p. 2218:
        D(XFe, T, fO2) = D(T, fO2) * 10^(m (XFe - 0.09)),  m = 1
    """
    logD = (p["logD0"]
            - p["Q"] * 1.0e3 / (LN10 * R_GAS * cond.T_K)
            + p["n"] * cond.log_fo2_Pa)
    XFe = cond.X.get("XFe", None)
    if XFe is not None:
        logD = logD + p["m"] * (np.asarray(XFe, dtype=float) - 0.09)
    return 10.0 ** logD


_add(DiffusionCoefficient(
    key="opx_FeMg_dohmen2016",
    kind="interdiffusion",
    mineral="opx", species="Fe-Mg",
    label="Opx Fe-Mg, Dohmen et al. (2016) -- Fs9, fO2-dependent",
    citation="dohmen2016",
    equation_number="1 (+ compositional correction, p. 2218)",
    equation_text=("log D_Fe-Mg^c = log D0 - {Q / ln(10) R T} + n log fO2   with "
                   "Q = 308 +/- 23 kJ/mol, log(D0 [m2/s]) = -5.95 +/- 0.83, n = 0.053 +/- 0.027. "
                   "Abstract form D_Fe-Mg [m2/s] = 1.12e-6 (fO2 [Pa])^0.053 exp[-308 kJ/mol /(R T)]. "
                   "D(XFe,T,fO2) = D(T,fO2) * 10^(m (XFe - 0.09)) with m = 1"),
    func=_dohmen2016_fs9,
    params={
        "logD0": Parameter("logD0", -5.95, 0.83, "log10(m2/s)", "unstated", "pre-exponential factor"),
        "Q": Parameter("Q", 308.0, 23.0, "kJ/mol", "unstated", "activation energy"),
        "n": Parameter("n", 0.053, 0.027, "-", "unstated", "fO2 exponent (fO2 in Pa)"),
        "m": Parameter("m", 1.0, 0.0, "-", "1s", "composition exponent, 10^(m(XFe-0.09))"),
    },
    sigma_logD=0.20,
    axis_factors={"c": 1.0, "b": 1.0, "a": 1.0 / 3.5},
    reference_axis="c",
    requires=("XFe",),
    needs_fo2=True, fo2_unit="Pa",
    T_range=Range(1143.15, 1373.15, "K (870-1100 C, as in the abstract; the run table lists 872-1104 C)"),
    P_range=Range(1.0e5, 101325.0, "Pa (atmospheric; 1 bar to 1 atm)"),
    fo2_range=Range(-11.0, -7.0, "log10 Pa"),
    X_range=Range(0.09, 0.5, "XFe (0.09-0.5 is the range for which the authors recommend the composition term; the fit itself is Fs9 only)"),
    verified=True,
    verified_from="read from the paper PDF: abstract, eq. 1 (p. 2215), anisotropy p. 2215, composition p. 2218, run table (Table 2) p. 2212",
    secondary_citations=("sato2022", "polo_sanchez2023", "dias_dohmen2024"),
    recommended=False,
    superseded_by="dias2025",
    superseded_note=("Dias & Dohmen (2024) and Dias, Dohmen & Behrens (2025) refitted these "
                     "experiments with a temperature-dependent composition exponent m instead "
                     "of the fixed m = 1 used here, added new data for XFe up to 0.4, and give "
                     "separate laws above and below log fO2 = -10 Pa. Between 900 and 1100 C "
                     "the two agree to within about 0.3 log units for Fs10 at log fO2 = -7 Pa, "
                     "but differ by up to 0.9 log units for Fe-rich opx at 1100 C and by up "
                     "to 1 log unit under reducing conditions. Prefer opx_FeMg_dias2025."),
    notes=("Reference axis is [001]. D//[100] = D//[001]/3.5 (abstract and p. 2215: Da is slower "
           "than Dc by a factor of 3.5). D//[010] is 'similar' to D//[001] in the abstract, and "
           "p. 2215 gives Dc >= Db > Da from sparse b-axis data (two runs, Table 2). Diffusor sets "
           "D//[010] = D//[001], which is an upper bound for b, not a measured ratio. The "
           "compositional term is only recommended for XFe = 0.09-0.5 (p. 2218): the Eq. 1 fit "
           "itself comes from one crystal, Opx7 (XFe = 0.09), so every other XFe is "
           "extrapolation by that recommendation. "
           "Test case: at 950 C and log fO2 = -7 Pa the equation gives log D = -19.47 "
           "(-19.474), against the measured -19.49 +/- 0.07 (run OPXD_14, sample 7D_20). "
           "Uncertainties: the printed +/- on log D0, Q and n carry no stated confidence level "
           "(neither 1 sigma nor 2 sigma nor standard error is given), so they are labelled "
           "'level not stated'. The paper states overall experimental uncertainties of "
           "+/-0.2 log units on individual D values from reproducibility and time-series "
           "experiments (p. 2214, level also not stated). sigma_logD = 0.2 is that value, "
           "used as the scatter about the line; no scatter of the fit itself is given."),
))


def _dohmen2016_fs1(dc, cond: Conditions, p):
    """Opx8 (XFe = 0.01): single Arrhenius, no resolvable fO2 dependence."""
    return arrhenius(10.0 ** p["logD0"], p["Q"] * 1.0e3, cond.T_K)


_add(DiffusionCoefficient(
    key="opx_FeMg_dohmen2016_fs1",
    kind="interdiffusion",
    mineral="opx", species="Fe-Mg",
    label="Opx Fe-Mg, Dohmen et al. (2016) -- Fs1, fO2-independent",
    citation="dohmen2016",
    equation_number="1 with n = 0 (Opx8 fit, p. 2215)",
    equation_text=("D_Fe-Mg [m2/s] = 1.66e-4 exp[-377 +/- 30 kJ/mol /(R T)]. "
                   "log(D0 [m2/s]) = -3.78 +/- 1.26, Q = 377 +/- 30 kJ/mol"),
    func=_dohmen2016_fs1,
    params={
        "logD0": Parameter("logD0", -3.78, 1.26, "log10(m2/s)", "unstated", "pre-exponential factor"),
        "Q": Parameter("Q", 377.0, 30.0, "kJ/mol", "unstated", "activation energy"),
    },
    sigma_logD=0.20,
    axis_factors={"c": 1.0, "b": 1.0, "a": 1.0 / 3.5},
    reference_axis="c",
    needs_fo2=False,
    T_range=Range(1143.15, 1373.15, "K (870-1100 C, as in the abstract; the run table lists 872-1104 C)"),
    P_range=Range(1.0e5, 101325.0, "Pa (atmospheric; 1 bar to 1 atm)"),
    X_range=Range(0.0, 0.05, "XFe (Diffusor assumption around the single calibrated composition XFe = 0.01; not stated by the authors)"),
    verified=True,
    verified_from="read from the paper PDF: abstract and p. 2215",
    calibration_notes=("Anisotropy of this entry is Diffusor's assumption: the factors D//[100] = D//[001]/3.5 "
                       "and D//[010] = D//[001] were measured only on the Fs9 crystal (Opx7); no [100] or "
                       "[010] experiment exists for the Fs1 crystal (Opx8). Only the [001] direction is "
                       "calibrated.",),
    notes=("For near-end-member enstatite only (XFe ~ 0.01, crystal Opx8). No fO2 dependence resolved. "
           "The XFe range shown (0 to 0.05) is not from the paper, which calibrates the single "
           "composition XFe = 0.01; because this entry has no composition term it is not "
           "checked against the host composition. The a-axis factor 1/3.5 and the b-axis factor 1 are "
           "taken from the Fs9 crystal (Opx7) of the same study and applied here without data for Opx8. "
           "Uncertainties on log D0 and Q have no stated confidence level. sigma_logD = 0.2 is the "
           "overall experimental uncertainty the paper states for individual D values (p. 2214)."),
))


# ---------------------------------------------------------------------------
# Dias, Dohmen & Behrens (2025) GCA 395, 195-211 -- two fO2 regimes
# ---------------------------------------------------------------------------
LOG_FO2_PA_SWITCH = -10.0   # Dias et al. (2025) section 4.3: eq. 22 above, eq. 23 at or below


def _dias2025_m(p, T_K, regime: int):
    """Composition exponent m(T), eqs 24 and 25 (note the opposite signs)."""
    if regime == 1:
        return p["m1_slope"] * 1.0e4 / T_K + p["m1_int"]
    return p["m2_slope"] * 1.0e4 / T_K + p["m2_int"]


def _dias2025(dc, cond: Conditions, p):
    """D = D0 (fO2/1e-7 Pa)^n exp(-Q/RT) 10^(m (XFe - 0.1)), regime by log fO2 [Pa]."""
    lf = cond.log_fo2_Pa
    XFe = np.asarray(cond.x("XFe"), dtype=float)
    if lf > LOG_FO2_PA_SWITCH:
        logD = (np.log10(p["D0_1"]) + p["n_1"] * (lf + 7.0)
                - p["Q_1"] * 1.0e3 / (LN10 * R_GAS * cond.T_K)
                + _dias2025_m(p, cond.T_K, 1) * (XFe - 0.1))
    else:
        logD = (np.log10(p["D0_2"])
                - p["Q_2"] * 1.0e3 / (LN10 * R_GAS * cond.T_K)
                + _dias2025_m(p, cond.T_K, 2) * (XFe - 0.1))
    return 10.0 ** logD


def dias2025_regime(log_fo2_Pa: float) -> str:
    """Which of the two Dias et al. (2025) equations applies at this fO2."""
    return ("eq. 22 (fO2-dependent, log fO2 > -10 Pa)" if log_fo2_Pa > LOG_FO2_PA_SWITCH
            else "eq. 23 (fO2-independent, log fO2 <= -10 Pa)")


_add(DiffusionCoefficient(
    key="opx_FeMg_dias2025",
    kind="interdiffusion",
    mineral="opx", species="Fe-Mg",
    label="Opx Fe-Mg, Dias, Dohmen & Behrens (2025) -- T, XFe and fO2",
    citation="dias2025",
    equation_number="22-25",
    equation_text=(
        "log fO2 [Pa] > -10:  D [m2/s] = 3.085e-8 (fO2[Pa]/1e-7)^0.25 exp[-(284 +/- 19 kJ/mol)/(R T)] "
        "10^(m1 (XFe - 0.1)),  m1 = -2.37(+/-0.45) (1e4/T) + 21.09   (eqs 22, 24). "
        "log fO2 [Pa] <= -10:  D [m2/s] = 1.93e-10 exp[-(246 +/- 78 kJ/mol)/(R T)] "
        "10^(m2 (XFe - 0.1)),  m2 = 2.96(+/-0.49) (1e4/T) - 21.08   (eqs 23, 25). T in K"),
    func=_dias2025,
    params={
        "D0_1": Parameter("D0_1", 3.085e-8, 0.0, "m2/s", "1s", "pre-exponential factor, eq. 22"),
        "n_1": Parameter("n_1", 0.25, 0.0, "-", "1s", "fO2 exponent, eq. 22 (fO2 in Pa, ref 1e-7 Pa)"),
        "Q_1": Parameter("Q_1", 284.0, 19.0, "kJ/mol", "unstated", "activation energy, eq. 22"),
        "m1_slope": Parameter("m1_slope", -2.37, 0.45, "-", "unstated", "slope of m1 on 1e4/T, eq. 24"),
        "m1_int": Parameter("m1_int", 21.09, 0.0, "-", "1s", "intercept of m1, eq. 24"),
        "D0_2": Parameter("D0_2", 1.93e-10, 0.0, "m2/s", "1s", "pre-exponential factor, eq. 23"),
        "Q_2": Parameter("Q_2", 246.0, 78.0, "kJ/mol", "unstated", "activation energy, eq. 23"),
        "m2_slope": Parameter("m2_slope", 2.96, 0.49, "-", "unstated", "slope of m2 on 1e4/T, eq. 25"),
        "m2_int": Parameter("m2_int", -21.08, 0.0, "-", "1s", "intercept of m2, eq. 25"),
    },
    sigma_logD=0.34,     # p. 208: average misfit of the regressions
    axis_factors={"c": 1.0, "b": 1.0, "a": 1.0 / 3.5},
    reference_axis="c",
    requires=("XFe",),
    needs_fo2=True, fo2_unit="Pa",
    T_range=Range(1173.15, 1373.15, "K (900-1100 C, the interval Section 4.3 gives for both equations. The new experiments ran at 950-1100 C, and Section 4.3 fits them together with the refitted Dohmen et al. (2016) data, which reach lower temperatures (the sentence introducing eq. 23 names 870 C))"),
    P_range=Range(1.0e5, 101325.0, "Pa (atmospheric; 1 bar to 1 atm)"),
    fo2_range=Range(-11.0, -7.0, "log10 Pa"),
    X_range=Range(0.1, 0.4, "XFe"),
    verified=True,
    verified_from=("read from the paper PDF: abstract, eqs 22-25 (p. 208) and Table 3. "
                   "Cross-checked against Table 1 of Dias & Dohmen (2024), which it reproduces "
                   "to within 0.1 log units at 950, 1050 and 1100 C and 0.3 at 1000 C"),
    secondary_citations=("dias_dohmen2024", "dohmen2016"),
    recommended=True,
    notes=(
        "Measured along [001] only. Both anisotropy factors are Diffusor's carry-over from "
        "Dohmen et al. (2016), not from this paper, which did not measure anisotropy and only "
        "says (introduction) that diffusion along c is faster than along a or b: the "
        "a-axis factor of 1/3.5 is from Dohmen et al. (abstract and p. 2215), and the b-axis "
        "factor of 1 follows their statement that [010] is 'similar' to [001] (Dc >= Db > Da, "
        "sparse b data), so it is an upper bound. The two equations meet at log fO2 = -10 "
        "Pa but do not join smoothly there: Dias et al. found a change of diffusion mechanism, "
        "with D independent of fO2 under more reducing conditions and the composition effect "
        "m decreasing with temperature. The composition exponent m "
        "falls from about 3.8 at 1100 C to 0.9 at 900 C in the oxidised regime (eq. 24). "
        "That linear m1(T) reaches zero at 850.6 C and is negative below, so it must not be "
        "extrapolated downward: below the experimental range the authors recommend the corrected "
        "calibration of Kroll et al. (1997) (their eq. 7) instead, and Diffusor does not implement "
        "that and does not clip m. The printed +/- on Q and on the slopes of m (eqs 22-25) have no "
        "stated confidence level (only the Table 3 values of D are stated as 1 sigma), so they are "
        "labelled 'level not stated'. The pre-exponential factors carry no published uncertainty, "
        "so the Monte Carlo samples log D at T with a 0.34 log-unit scatter, the average misfit "
        "with which the regressions reproduce the experiments (p. 208). Natural arc magmas usually "
        "sit at log fO2 above -7 Pa, above the calibrated range (this statement is Diffusor's, not "
        "the paper's), so the fO2^(1/4) term is extrapolated there. Diffusor's comparison with "
        "Dohmen et al. (2016) at 950 C and NNO+1 (log fO2 about -5.0 Pa, buffer from Diffusor's "
        "own NNO law): this law is faster by 0.21 log units at XFe = 0.1, 0.36 at 0.3 and 0.43 "
        "at 0.4, so the difference is about 0.2-0.4 log units depending on XFe. The authors "
        "advise against extrapolating to conditions (fO2 and/or T) outside the range explored."),
))


# ---------------------------------------------------------------------------
# Dias & Dohmen (2024) CMP 179:36 -- first calibration of the Fe effect
# ---------------------------------------------------------------------------
def _dias2024(dc, cond: Conditions, p):
    """D = 3.8e-9 exp(-261 kJ/mol / RT) 10^(m (XFe - 0.1)), m = -2.711e4/T + 23.5408 (eqs 7, 12, 13)."""
    XFe = np.asarray(cond.x("XFe"), dtype=float)
    m = p["m_slope"] / cond.T_K + p["m_int"]
    logD = (np.log10(p["D0"]) - p["Q"] * 1.0e3 / (LN10 * R_GAS * cond.T_K)
            + m * (XFe - 0.1))
    return 10.0 ** logD


_add(DiffusionCoefficient(
    key="opx_FeMg_dias_dohmen2024",
    kind="interdiffusion",
    mineral="opx", species="Fe-Mg",
    label="Opx Fe-Mg, Dias & Dohmen (2024) -- log fO2 = -7 Pa only",
    citation="dias_dohmen2024",
    equation_number="12-13",
    equation_text=("D_Fe-Mg [m2/s] = 3.8e-9 exp[-(261.07 +/- 24 kJ/mol)/(R T)] at XFe = 0.1 and "
                   "log fO2 = -7 Pa (eq. 12). D(XFe) = D(XFe = 0.1) 10^(m (XFe - 0.1)) (eq. 7) with "
                   "m = -2.711e4 / T[K] + 23.5408 (eq. 13)"),
    func=_dias2024,
    params={
        "D0": Parameter("D0", 3.8e-9, 0.0, "m2/s", "1s", "pre-exponential factor"),
        "Q": Parameter("Q", 261.07, 24.0, "kJ/mol", "unstated", "activation energy"),
        "m_slope": Parameter("m_slope", -2.711e4, 0.0, "K", "1s", "slope of m on 1/T"),
        "m_int": Parameter("m_int", 23.5408, 0.0, "-", "1s", "intercept of m"),
    },
    sigma_logD=0.2,
    axis_factors={"c": 1.0, "b": 1.0, "a": 1.0 / 3.5},
    reference_axis="c",
    requires=("XFe",),
    needs_fo2=False,
    T_range=Range(1223.15, 1373.15, "K (950-1100 C)"),
    P_range=Range(1.0e5, 101325.0, "Pa (atmospheric; 1 bar to 1 atm)"),
    X_range=Range(0.09, 0.5, "XFe (Diffusor's inference, not stated by the authors: Opx7 crystal 0.09 to the Fs50 target film)"),
    verified=True,
    verified_from=("read from the paper PDF: abstract, eqs 12-13 (pp. 16-17) and Table 1. The "
                   "typeset m(T) expression is ambiguous in the text layer. The form -2.711e4/T "
                   "+ 23.5408 is the one that reproduces the fitted m values of Table 1 (3.7, "
                   "3.0, 2.4 and 1.1 at 1102, 1050, 1000 and 950 C) and matches the 1/T form of "
                   "Dias et al. (2025) eq. 2"),
    secondary_citations=("dohmen2016",),
    recommended=False,
    superseded_by="dias2025",
    superseded_note=("Dias, Dohmen & Behrens (2025) extended these experiments to XFe up to 0.4 "
                     "and log fO2 down to -11 Pa and give the general equations that replace "
                     "this one. At log fO2 = -7 Pa the two agree to within about 0.1 log units "
                     "for Fs9."),
    notes=("Combined regression of the new data with the Dohmen et al. (2016) Fs9 data refitted "
           "with the temperature-dependent m. Only valid at log fO2 = -7 Pa (all experiments were "
           "run at that nominal fO2, which is about FMQ-2.5 at 1102 C and about FMQ at 950 C, "
           "Table 1); the restriction to -7 Pa is Diffusor's wording, the paper only states that "
           "an fO2 dependence of m cannot be excluded. It has no fO2 term. Measured along [001] "
           "only. The anisotropy factors are not from this paper, which says only that diffusion "
           "along c is faster than along other directions: D//[100] = D//[001]/3.5 and D//[010] = "
           "D//[001] are carried over from Dohmen et al. (2016) (a: abstract and p. 2215; b: "
           "'similar', Dc >= Db > Da from sparse data, so the b-axis factor of 1 is an upper "
           "bound). The XFe range shown is not stated by the authors; it is inferred from the "
           "crystal (XFe = 0.09), the Fs50 target film and the compositions plotted in Figs 9b "
           "and 10b (XFe = 0.1, 0.3, 0.5). The quoted +/-24 kJ/mol on Q has no stated confidence "
           "level."),
))


# ---------------------------------------------------------------------------
# REE: Dias, Dohmen & Hartmann (2025) GCA 410, 85-100
# ---------------------------------------------------------------------------
def _dias2025ree(dc, cond: Conditions, p):
    logD = np.log10(p["D0"]) - p["Q"] * 1.0e3 / (LN10 * R_GAS * cond.T_K)
    if p.get("inv_n", 0.0):
        from ..thermo.buffers import log_fo2_buffer
        d_iw = cond.log_fo2_bar - log_fo2_buffer("IW", cond.T_K, cond.P_Pa)
        logD = logD + p["inv_n"] * d_iw
    return 10.0 ** logD


for _el, _D0, _Q, _sQ, _eq, _inv_n in [
        ("Lu", 1.51e-9, 263.0, 52.0, "6", 1.0 / 7.0),
        ("Ce", 5.75e-14, 166.0, 40.0, "7", 0.0),
        ("Eu", 1.90e-14, 147.0, 22.0, "8", 0.0)]:
    _is_lu = _inv_n > 0
    _add(DiffusionCoefficient(
        key=f"opx_{_el}_dias2025",
        kind="chemical",
        mineral="opx", species=_el,
        label=f"Opx {_el}, Dias, Dohmen & Hartmann (2025)",
        citation="dias2025ree",
        equation_number=_eq,
        equation_text=(f"D_{_el} [m2/s] = {_D0:g}"
                       + (" (fO2[Pa]/fO2[Pa]_0)^(1/7)" if _is_lu else "")
                       + f" exp[-({_Q:g} +/- {_sQ:g} kJ/mol)/(R T)]"),
        func=_dias2025ree,
        params={
            "D0": Parameter("D0", _D0, 0.0, "m2/s", "1s", "pre-exponential factor"),
            "Q": Parameter("Q", _Q, _sQ, "kJ/mol", "unstated", "activation energy"),
            "inv_n": Parameter("inv_n", _inv_n, 0.0, "-", "1s",
                               "fO2 exponent 1/n relative to the IW buffer (Lu only)"),
        },
        sigma_logD=0.3,
        reference_axis="c",
        needs_fo2=_is_lu, fo2_unit="bar",
        T_range=Range(1223.15, 1373.15, "K (950-1100 C as printed; the IW-buffered experiments are 1000-1100 C)"),
        P_range=Range(1.0e5, 101325.0, "Pa (atmospheric; not stated by the authors, inferred from the CO-CO2 gas-mixing furnaces; 1 bar to 1 atm)"),
        fo2_range=Range(None, None, "near the IW buffer"),
        verified=True,
        verified_from=(
            f"read from the paper PDF: abstract and eq. {_eq} (p. 93)"
            + (". The equation is transcribed as printed. The reference fugacity "
               "fO2[Pa]_0 is not defined in the text; the paper states the laws for 'fO2 close to the "
               "iron-wustite (IW) buffer' and the abstract calls them laws 'along the IW "
               "buffer', so Diffusor takes fO2_0 as the IW buffer at T. That choice is an "
               "interpretation, not part of the source. Check against the 12 Lu values of "
               "Table 3 (fO2 from IW-2 to IW+4): with fO2_0 = IW(T) eq. 6 reproduces them "
               "with an rms misfit of 0.10 log units, whereas the best fixed reference "
               "fugacity (log fO2_0 = -9 or -8 Pa) gives 0.16 and no fO2 term gives 0.23" if _is_lu else "")),
        notes=(("Along [001] in natural opx (crystals Opx7 and Opx31, same provenance and "
                "composition), measured by TOF-SIMS on thin-film diffusion couples. The laws are "
                "stated for 950-1100 C and fO2 close to the IW buffer, but only the 1000-1100 C "
                "runs are near IW (the 950 C run is at about IW+4), so the IW-level Arrhenius fit "
                "rests on 1000-1100 C data. The printed +/- on Q has no stated confidence level. "
                "The transport kind 'chemical' (dilute trace element in a host of fixed major-element "
                "composition) is Diffusor's classification; the authors only speak of diffusion "
                "coefficients. " )
               + ("Lu diffusion appears to depend on fO2 (exponent about 1/7 at 1080-1100 C) and "
                  "possibly on Lu concentration, suggesting a different mechanism from Ce and "
                  "Eu. " if _is_lu else
                  f"No resolvable dependence of {_el} diffusion on fO2 between 1e-7 and 1e-11 "
                  "Pa. " + ("Eu is assumed to be trivalent. " if _el == "Eu" else ""))
               + "Faster than earlier opx and diopside REE data (Cherniak & Liang 2007, Van "
               "Orman et al. 2001), and with lower activation energies. Relevant to Lu-Hf "
               "systematics and to diffusive fractionation of REE between opx and melt."),
    ))


# ---------------------------------------------------------------------------
# Ganguly & Tazzoli (1994) without an fO2 term, as used by Ostorero et al. (2022)
# ---------------------------------------------------------------------------
def _ganguly_tazzoli_no_fo2(dc, cond: Conditions, p):
    XFe = np.asarray(cond.x("XFe"), dtype=float)
    return 10.0 ** (p["c0"] + p["a"] * XFe - p["b"] / cond.T_K)


_add(DiffusionCoefficient(
    key="opx_FeMg_ganguly_tazzoli1994_nofo2",
    kind="interdiffusion",
    mineral="opx", species="Fe-Mg",
    label="Opx Fe-Mg, Ganguly & Tazzoli (1994) without fO2 term (Ostorero et al. 2022)",
    citation="ganguly_tazzoli1994",
    equation_number="Ostorero et al. (2022) eq. 1",
    equation_text="log D [cm2/s] = -5.54 + 2.6 XFe - 12530/T[K]  (= -9.54 + ... in m2/s)",
    func=_ganguly_tazzoli_no_fo2,
    params={
        "c0": Parameter("c0", -9.54, 0.0, "log10(m2/s)", "1s", "intercept (-5.54 in cm2/s)"),
        "a": Parameter("a", 2.6, 0.0, "-", "1s", "XFe coefficient"),
        "b": Parameter("b", 12530.0, 0.0, "K", "1s", "Q/(ln10 R) = 239.9 kJ/mol"),
    },
    sigma_logD=1.0,       # Ganguly & Tazzoli (1994) p. 934: "approximately one order of magnitude"
    axis_factors={"c": 1.0, "b": 1.0, "a": 1.0 / 3.5},
    reference_axis="c",
    requires=("XFe",),
    needs_fo2=False,
    T_range=Range(773.15, 1073.15, "K (500-800 C, extrapolated above)"),
    X_range=Range(0.10, 0.50, "XFe"),
    verified=True,
    verified_from=("read from the PDF of Ostorero et al. (2022), Methods, eq. 1, which states "
                   "the fO2 correction was deliberately omitted"),
    secondary_citations=("ostorero2022", "dohmen2016", "dias_dohmen2024"),
    calibration_notes=("Ganguly & Tazzoli (1994) state the law only for fO2 from the iron-wustite (IW) buffer "
                       "to 0.8 log units above IW; the fO2 window is not checked by Diffusor for this "
                       "entry (the law has no fO2 term), so make sure the crystal's fO2 lies in it or "
                       "accept the extrapolation.",),
    notes=("The same law as opx_FeMg_ganguly_tazzoli1994 with the hypothesised fO2^(1/6) term "
           "removed. Ostorero et al. (2022) used it this way for the Kizimen 2010 "
           "orthopyroxenes, modelled across the c-axis parallel to b, so select this entry to "
           "reproduce their timescales. Ganguly & Tazzoli (1994) formulate the law for 500-800 C, "
           "fO2 from IW to 0.8 log units above IW, and XFe 0.1-0.5. The law is the average "
           "D(Fe-Mg) along c and b, hence reference axis c with b equal to c. The a-axis factor "
           "1/3.5 is NOT from Ganguly & Tazzoli, who only say qualitatively that diffusion is "
           "slowest along a; the value 3.5 is the measured ratio of Dohmen et al. (2016), "
           "applied here by Diffusor. Dias "
           "& Dohmen (2024) conclude that it overestimates the composition dependence below "
           "about 1000 C, and that correcting for this lengthened the Mount St Helens "
           "timescales of Saunders et al. (2012) from about 5 to 11 weeks."),
))


# ---------------------------------------------------------------------------
# Ganguly & Tazzoli (1994), in the form used by Allan et al. (2013)
# ---------------------------------------------------------------------------
def _ganguly_tazzoli(dc, cond: Conditions, p):
    """log D [m2/s] = -9.54 + 2.6 XFe - 12530/T + (1/6) log(fO2_sample/fO2_IW)

    Quoted as Equation 3 of Dohmen et al. (2016), who attribute the general
    expression to Allan et al. (2013) derived from Ganguly & Tazzoli (1994).
    The fO2 term is optional: Ostorero et al. (2022) eq. 1 use the same
    expression without it (they write log D [cm2/s] = -5.54 + 2.6 XFe - 12530/T,
    which is the same law, 1 cm2/s = 1e-4 m2/s).
    """
    XFe = np.asarray(cond.x("XFe"), dtype=float)
    logD = p["c0"] + p["a"] * XFe - p["b"] / cond.T_K
    if p.get("use_fo2", 1.0) >= 0.5:
        from ..thermo.buffers import log_fo2_buffer
        d_iw = cond.log_fo2_bar - log_fo2_buffer("IW", cond.T_K, cond.P_Pa)
        logD = logD + p["n_fo2"] * d_iw
    return 10.0 ** logD


_add(DiffusionCoefficient(
    key="opx_FeMg_ganguly_tazzoli1994",
    kind="interdiffusion",
    mineral="opx", species="Fe-Mg",
    label="Opx Fe-Mg, Ganguly & Tazzoli (1994) / Allan et al. (2013) form",
    citation="ganguly_tazzoli1994",
    equation_number="Dohmen et al. (2016) eq. 3. Ostorero et al. (2022) eq. 1",
    equation_text=("log D_Fe-Mg [m2/s] = -9.54 + 2.6 XFe - 12530/T[K] "
                   "+ (1/6) log( fO2(sample,T) / fO2(IW,T) ). The first three terms are "
                   "Ganguly & Tazzoli (1994) eq. 12 (-5.54 in cm2/s); the fO2 term is the "
                   "addition of Allan et al. (2013) quoted by Dohmen et al. (2016) eq. 3"),
    func=_ganguly_tazzoli,
    params={
        "c0": Parameter("c0", -9.54, 0.0, "log10(m2/s)", "1s", "intercept"),
        "a": Parameter("a", 2.6, 0.0, "-", "1s", "XFe coefficient"),
        "b": Parameter("b", 12530.0, 0.0, "K", "1s", "Q/(ln10 R) = 239.9 kJ/mol"),
        "n_fo2": Parameter("n_fo2", 1.0 / 6.0, 0.0, "-", "1s", "fO2 exponent relative to IW"),
        "use_fo2": Parameter("use_fo2", 1.0, 0.0, "-", "1s", "1 = apply fO2 term, 0 = omit it"),
    },
    sigma_logD=1.0,       # Ganguly & Tazzoli (1994) p. 934: "approximately one order of magnitude"
    axis_factors={"c": 1.0, "b": 1.0, "a": 1.0 / 3.5},
    reference_axis="c",
    requires=("XFe",),
    needs_fo2=True, fo2_unit="bar",
    T_range=Range(773.15, 1073.15, "K (500-800 C, extrapolated above)"),
    fo2_range=Range(None, None, "relative to IW"),
    X_range=Range(0.10, 0.50, "XFe"),
    verified=True,
    verified_from=("Ganguly & Tazzoli (1994) eq. 12 (p. 934), read from the rendered page on 2 October "
                   "2026, and Ostorero et al. (2022) eq. 1. The fO2 term is the form of Allan et al. "
                   "(2013) as quoted by Dohmen et al. (2016) eq. 3; Ganguly & Tazzoli only speculate "
                   "an fO2^(1/6) dependence (p. 935), and Diffusor did not consult Allan et al. (2013) itself"),
    secondary_citations=("dohmen2016", "ostorero2022"),
    calibration_notes=("The first three terms are Ganguly & Tazzoli (1994) eq. 12, valid for 500-800 C, "
                       "fO2 from the iron-wustite (IW) buffer to 0.8 log units above IW, and XFe 0.10-0.50. "
                       "The (1/6) log(fO2/fO2_IW) term is NOT in that paper: Ganguly & Tazzoli only "
                       "speculate (p. 935, by analogy with olivine and garnet) that D varies roughly as "
                       "fO2^(1/6), and give no normalisation. The form with fO2 relative to IW(T) is "
                       "the one of Allan et al. (2013) as quoted in eq. 3 of Dohmen et al. (2016). "
                       "Diffusor does not check the fO2 window of this entry (a buffer-relative window "
                       "cannot be expressed as an absolute range): at fO2 above IW + 0.8 the fO2 term "
                       "is an extrapolation of a speculated exponent (about +0.8 log units of D at "
                       "NNO, 800-1000 C).",),
    notes=("Retrieved indirectly from Fe-Mg order-disorder kinetics. It comes from no "
           "direct diffusion measurement. The law is the average D(Fe-Mg) along c and b "
           "(reference axis c, b equal to c). The a-axis factor 1/3.5 is not from Ganguly & Tazzoli, "
           "who only say qualitatively that diffusion is slowest along a; the value 3.5 is the "
           "measured ratio of Dohmen et al. (2016), applied by Diffusor. "
           "The fO2 term is optional and on by default because it is the form used in the "
           "literature (Allan et al. 2013) that this entry is named after; it is zero at IW and "
           "at most +0.13 in log D within the IW to IW+0.8 window of the paper, well below the "
           "stated uncertainty of one order of magnitude, but it grows outside it. Dohmen et al. "
           "(2016, p. 2216) show that this strong fO2 dependence (exponent 1/6, against their "
           "measured n = 0.053) gives an artificially low activation energy (150 kJ/mol) at "
           "constant fO2, against their measured 308 kJ/mol, and deviates from their data at many "
           "temperatures. Set use_fo2 = 0 to drop the fO2 term, as Ostorero et al. (2022) did "
           "(or select opx_FeMg_ganguly_tazzoli1994_nofo2). Kept for comparison with older studies."),
))


# ---------------------------------------------------------------------------
# Schwandt, Cygan & Westrich (1998) -- Mg self-diffusion in orthoenstatite
# ---------------------------------------------------------------------------
def _schwandt(dc, cond: Conditions, p):
    return arrhenius(10.0 ** p["logD0"], p["Q"] * 1.0e3, cond.T_K)


for _axis, _label, _logD0, _slog, _Q, _sQ in [
        ("a", "(100)", -3.96, 2.48, 360.0, 52.0),
        ("b", "(010)", -5.16, 3.69, 339.0, 77.0),
        ("c", "(001)", -8.36, 3.27, 265.0, 66.0)]:
    _add(DiffusionCoefficient(
        key=f"opx_Mg_schwandt1998_{_axis}",
        kind="tracer",
        mineral="opx", species="Mg",
        label=f"Opx Mg self-diffusion // {_axis}-axis {_label}, Schwandt et al. (1998)",
        citation="schwandt1998",
        equation_number="Table 3",
        equation_text=(f"D = D0 exp(-Ea/RT). {_label}: Ea = {_Q:g} +/- {_sQ:g} kJ/mol, "
                       f"log D0 [m2/s] = {_logD0:g} +/- {_slog:g}"),
        func=_schwandt,
        params={
            "logD0": Parameter("logD0", _logD0, _slog, "log10(m2/s)", "1s", "pre-exponential factor"),
            "Q": Parameter("Q", _Q, _sQ, "kJ/mol", "1s", "activation energy"),
        },
        reference_axis=_axis,
        needs_fo2=False,
        T_range=Range(1023.15, 1173.15, "K (750-900 C)"),
        P_range=Range(1.0e5, 101325.0, "Pa (atmospheric; 1 bar to 1 atm)"),
        X_range=Range(0.10, 0.10, "En90Fs10 (Table 2 and 3 captions; abstract and Table 1 give En88Fs12)"),
        verified=True,
        verified_from="read from the paper PDF: abstract and Table 3",
        notes=("25Mg self-diffusion by SIMS depth profiling at the IW buffer on orthoenstatite. "
               "The abstract and the "
               "microprobe Table 1 give the host as En88Fs12 (Fs 0.11-0.13 for the three surfaces), "
               "while the captions of Tables 2 and 3 give En90Fs10. Diffusor shows the Table 3 caption value "
               "(XFs = 0.10) because the parameters come from that table; the difference does not "
               "enter the calculation. The anisotropy in D is very small (the three axes lie "
               "within 0.35 log units of each other at 900 C, p. 394), but the c-axis activation "
               "energy is reported as statistically smaller than the a-axis value, with "
               "D(001) > D(010) > D(100). Dohmen et al. (2016, p. 2216) report that these tracer "
               "data are very similar to their Eq. 1 interdiffusion coefficients at fO2 = IW "
               "(Dohmen et al. is a later paper; this sentence is not in Schwandt et al.). "
               "Uncertainties are the maximum 1 sigma values of Table 3."),
    ))
