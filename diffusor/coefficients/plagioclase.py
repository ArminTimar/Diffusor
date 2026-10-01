"""Diffusion coefficients for plagioclase.

Trace-element diffusion in plagioclase is coupled to the anorthite gradient:
the driving force is the chemical-potential gradient, not the concentration
gradient, so the flux carries an extra term (Costa et al. 2003 eq. 7; Dohmen,
Faak & Blundy 2017 Appendix eq. A7-A8):

    J_i = -D_i dC_i/dx + (D_i C_i / R T) A_i dX_An/dx

with ``A_i`` the slope of ``-R T ln(gamma_i)`` against X_An.  The A values are
in :data:`ACTIVITY_A` and are applied by the solver as
``theta = A_i / (R T)``.  X_An is treated as frozen because NaSi-CaAl
interdiffusion is orders of magnitude slower (Grove et al. 1984).
"""
from __future__ import annotations

import numpy as np

from ..constants import R_GAS
from .base import LN10, Conditions, DiffusionCoefficient, Parameter, Range, arrhenius

COEFFICIENTS = []


def _add(c):
    COEFFICIENTS.append(c)
    return c


# --- activity (non-ideality) slopes A_i, kJ/mol -------------------------------
# Source: Dohmen, Faak & Blundy (2017) RiMG 83, Appendix Fig. A1 caption and
# Table 1, calculated at 1200 C from the lattice-strain model of Dohmen &
# Blundy (2014).  -R T ln(gamma_i) = A_i X_An + B_i  (Appendix eq. A6).
ACTIVITY_A = {
    "Mg": 15.8,
    "Sr": -17.4,
    "Ba": -35.1,
    "Li": -1.7,
    "K": -8.0,
    "Rb": -15.9,
}
ACTIVITY_A_CITATION = "dohmen2017"
ACTIVITY_A_NOTE = (
    "A_i values (kJ/mol) at 1200 C from Dohmen, Faak & Blundy (2017) RiMG 83, Appendix "
    "Fig. A1 caption / Table 1, derived from the lattice-strain model of Dohmen & Blundy "
    "(2014). Table 1 of that appendix also lists values for 900 C. Costa et al. (2003) "
    "used A_Mg = -RT ln(gamma) slope of the same form. The DMG Short Course 2025 script "
    "Diff_Model_Sr_in_Plag_implicit.m uses A_Sr = -15.1 kJ/mol."
)


def activity_theta(species: str, T_K: float) -> float:
    """theta = A_i / (R T), the coefficient of the dX_An/dx flux term."""
    if species not in ACTIVITY_A:
        return 0.0
    return ACTIVITY_A[species] * 1.0e3 / (R_GAS * T_K)


def equilibrium_profile(X_An, T_K: float, species: str, C_ref: float, X_An_ref: float = None):
    """Quasi-steady-state ("equilibrated") profile for a frozen X_An(x).

    ``C_eq(x) = C0 exp(A_i X_An(x) / (R T))``  -- Dohmen et al. (2017) Appendix
    eq. A13, with C0 fixed by the rim condition eq. A14 when ``X_An_ref`` is
    given (C_ref is then the rim concentration).

    This is the state a fast trace element relaxes to, and is the physically
    correct initial condition when the crystal grew in equilibrium with one melt.
    """
    X_An = np.asarray(X_An, dtype=float)
    th = ACTIVITY_A.get(species, 0.0) * 1.0e3 / (R_GAS * T_K)
    if X_An_ref is None:
        X_An_ref = float(X_An[0])
    return C_ref * np.exp(th * (X_An - X_An_ref))


# ---------------------------------------------------------------------------
# Mg
# ---------------------------------------------------------------------------
def _vanorman2014(dc, cond: Conditions, p):
    """ln D = lnD0 + a X_An - Q/(R T)  (Van Orman et al. 2014 eq. 4)."""
    XAn = np.asarray(cond.x("XAn"), dtype=float)
    lnD = p["lnD0"] + p["a"] * XAn - p["Q"] * 1.0e3 / (R_GAS * cond.T_K)
    return np.exp(lnD)


_add(DiffusionCoefficient(
    key="plag_Mg_vanorman2014",
    mineral="plagioclase", species="Mg",
    label="Plagioclase Mg, Van Orman et al. (2014)",
    citation="vanorman2014",
    equation_number="4",
    equation_text=("ln D = (-6.06 +/- 1.10) - (7.96 +/- 0.42) x_An - (287 +/- 10 kJ/mol)/(R T), "
                   "uncertainties 2 sigma, D in m2/s"),
    func=_vanorman2014,
    params={
        "lnD0": Parameter("lnD0", -6.06, 1.10, "ln(m2/s)", "2s", "intercept"),
        "a": Parameter("a", -7.96, 0.42, "-", "2s", "x_An coefficient"),
        "Q": Parameter("Q", 287.0, 10.0, "kJ/mol", "2s", "activation energy"),
    },
    sigma_logD=0.25,
    requires=("XAn",),
    needs_fo2=False,
    T_range=Range(1073.15, 1423.15, "K (800-1150 C)"),
    P_range=Range(1.0e5, 1.0e5, "Pa (1 atm)"),
    X_range=Range(0.23, 0.95, "x_An"),
    verified=True,
    verified_from="read from the paper PDF: abstract and eq. 4 on p. 84",
    recommended=True,
    notes=("Fitted to anorthite (An93), labradorite (An67), andesine (An43) and oligoclase "
           "(An23) plus the An95 data of LaTourrette & Wasserburg (1998). Little anisotropy "
           "was found between the b and c directions, and the authors state that treating Mg "
           "diffusion in plagioclase as isotropic is adequate for most applications. Reduced "
           "chi-squared of the fit was 2.88."),
))


def _costa2003_mg(dc, cond: Conditions, p):
    XAn = np.asarray(cond.x("XAn"), dtype=float)
    lnD = p["lnD0"] + p["a"] * XAn - p["Q"] * 1.0e3 / (R_GAS * cond.T_K)
    return np.exp(lnD)


_add(DiffusionCoefficient(
    key="plag_Mg_costa2003",
    mineral="plagioclase", species="Mg",
    label="Plagioclase Mg, Costa et al. (2003) -- older estimate",
    citation="costa2003",
    equation_number="as re-written by Van Orman et al. (2014), p. 84",
    equation_text="ln D = -6.07 - 9.44 x_An - (266 kJ/mol)/(R T), D in m2/s",
    func=_costa2003_mg,
    params={
        "lnD0": Parameter("lnD0", -6.07, 0.0, "ln(m2/s)", "1s", "intercept"),
        "a": Parameter("a", -9.44, 0.0, "-", "1s", "x_An coefficient"),
        "Q": Parameter("Q", 266.0, 0.0, "kJ/mol", "1s", "activation energy"),
    },
    sigma_logD=0.3,
    requires=("XAn",),
    needs_fo2=False,
    T_range=Range(1073.15, 1473.15, "K"),
    X_range=Range(0.0, 1.0, "x_An"),
    verified=True,
    verified_from=("transcribed verbatim from Van Orman et al. (2014) p. 84, who re-write the "
                   "Costa et al. (2003) expression in this form. Not read from Costa et al. (2003)"),
    secondary_citations=("vanorman2014", "latourrette_wasserburg1998"),
    notes=("A preliminary estimate built on the compositional dependence of Sr diffusion. "
           "Van Orman et al. (2014) show it over-predicts D at low temperature and low An "
           "(a factor of 10 for albite at 850 C) and note that timescales based on it, "
           "including those of Druitt et al. (2012) for Santorini, are too short by factors "
           "of 2.25-5.5. Kept for reproducing published results. Avoid it for new work."),
))


# ---------------------------------------------------------------------------
# Mg -- Audetat, Grocolas & Mutch (2026), silica-activity dependent
# ---------------------------------------------------------------------------
def _audetat2026_mg(dc, cond: Conditions, p):
    """log10 D = a X_An + b - Q/(2.303 R T) + c (1 - aSiO2)   (Audetat et al. 2026 eq. 1)."""
    XAn = np.asarray(cond.x("XAn"), dtype=float)
    logD = (p["a"] * XAn + p["b"] - p["Q"] * 1.0e3 / (LN10 * R_GAS * cond.T_K)
            + p["c"] * (1.0 - p["aSiO2"]))
    return 10.0 ** logD


_add(DiffusionCoefficient(
    key="plag_Mg_audetat2026",
    mineral="plagioclase", species="Mg",
    label="Plagioclase Mg, Audetat et al. (2026) -- with silica activity",
    citation="audetat2026",
    equation_number="1",
    equation_text=("log10 D_Mg [m2/s] = -2.99(+/-0.35) X_An - 4.03(+/-0.63) "
                   "- [262,914(+/-15,529) / (2.303 R T)] - 1.87(+/-0.45)(1 - aSiO2)"),
    func=_audetat2026_mg,
    params={
        "a": Parameter("a", -2.99, 0.35, "-", "1s", "X_An coefficient"),
        "b": Parameter("b", -4.03, 0.63, "log10(m2/s)", "1s", "intercept"),
        "Q": Parameter("Q", 262.914, 15.529, "kJ/mol", "1s", "activation energy"),
        "c": Parameter("c", -1.87, 0.45, "-", "1s", "coefficient of (1 - aSiO2)"),
        "aSiO2": Parameter("aSiO2", 1.0, 0.0, "-", "1s",
                           "silica activity relative to quartz. 1 = quartz saturated"),
    },
    sigma_logD=0.25,
    requires=("XAn",),
    needs_fo2=False,
    T_range=Range(1073.15, 1473.15, "K (800-1200 C, the combined range of the two data sets)"),
    P_range=Range(1.0e5, 1.0e5, "Pa (1 atm)"),
    X_range=Range(0.23, 0.93, "X_An"),
    verified=True,
    verified_from=("read from the accepted-manuscript PDF of Audetat et al. (2026), eq. 1 "
                   "(manuscript lines 271-278). The confidence level of the +/- values is not "
                   "stated there. They are treated as 1 sigma"),
    secondary_citations=("faak2013", "vanorman2014"),
    recommended=True,
    notes=("A joint fit of the Faak et al. (2013) and Van Orman et al. (2014) experiments. "
           "Mg diffusion is faster at higher silica activity (Faak et al. 2013). Diffusor "
           "evaluates it at aSiO2 = 1 by default, which corresponds to quartz saturation. For "
           "a silica-undersaturated melt D is lower by up to 1.87 log units at aSiO2 = 0. The "
           "review recommends this law, together with the Mg partitioning model of Mutch et al. "
           "(2022), as the Mg calibration that gives timescales most consistent with other "
           "chronometers. The numbers come from the accepted manuscript. "
           "The Monte Carlo samples log D at T with a 0.25 log-unit scatter, because no "
           "covariance matrix is published."),
))


# ---------------------------------------------------------------------------
# Sr
# ---------------------------------------------------------------------------
def _grocolas_form(dc, cond: Conditions, p):
    """log10 D = a X_An + b - Q/(2.303 R T)   (Grocolas et al. 2025 eqs 7, 8, 12-14)."""
    XAn = np.asarray(cond.x("XAn"), dtype=float)
    return 10.0 ** (p["a"] * XAn + p["b"] - p["Q"] * 1.0e3 / (LN10 * R_GAS * cond.T_K))


def _compensated_covariance(s_a, s_b, s_Q):
    """Covariance over (a, b, Q) with log D0 and Q perfectly correlated.

    Grocolas et al. (2025, section 4.3) state that the uncertainties of log10 D0
    and Ea are strongly covariant (citing Mutch et al. 2021) and that sampling
    them independently overestimates the total uncertainty. In their own Monte
    Carlo they "assume that log10 D0 and Ea follow a linear trend without any
    uncertainty envelope", which they note slightly underestimates it. This
    matrix encodes exactly that assumption; the X_An slope a is independent.
    """
    return np.array([[s_a ** 2, 0.0, 0.0],
                     [0.0, s_b ** 2, s_b * s_Q],
                     [0.0, s_b * s_Q, s_Q ** 2]])


_add(DiffusionCoefficient(
    key="plag_Sr_grocolas2025",
    mineral="plagioclase", species="Sr",
    label="Plagioclase Sr, Grocolas et al. (2025)",
    citation="grocolas2025",
    equation_number="7",
    equation_text=("log10 D_Sr [m2/s] = -1.65(+/-0.24) X_An - 3.03(+/-1.16) "
                   "- [368,142(+/-27,141) / (2.303 R T)]"),
    func=_grocolas_form,
    params={
        "a": Parameter("a", -1.65, 0.24, "-", "1s", "X_An coefficient"),
        "b": Parameter("b", -3.03, 1.16, "log10(m2/s)", "1s", "intercept"),
        "Q": Parameter("Q", 368.142, 27.141, "kJ/mol", "1s", "activation energy"),
    },
    covariance=_compensated_covariance(0.24, 1.16, 27.141),
    cov_order=("a", "b", "Q"),
    sigma_logD=0.25,
    requires=("XAn",),
    needs_fo2=False,
    T_range=Range(1173.15, 1473.15, "K (900-1200 C)"),
    P_range=Range(1.0e5, 1.0e5, "Pa (1 atm, anhydrous)"),
    X_range=Range(0.28, 0.67, "X_An (oligoclase An28 and labradorite An67)"),
    verified=True,
    verified_from=("read from the paper PDF: abstract and eq. 7 (p. 9). Uncertainties are one "
                   "standard deviation from 2000 Monte Carlo resamplings (Fig. 7 caption)"),
    secondary_citations=("audetat2026",),
    recommended=True,
    notes=("Sr in-diffusion into oriented oligoclase and labradorite with the silica activity "
           "buffered by sol-gel mullite-cristobalite or mullite-corundum powders. No resolvable "
           "dependence on aSiO2 or crystal orientation, so the law is isotropic. The values are "
           "1-2 orders of magnitude slower than Giletti & Casserly (1994) and Cherniak & "
           "Watson (1994), and because the activation energy is higher the gap widens as "
           "temperature falls: about 1.6 log units at 1100 C, 2.2 at 900 C and 2.8 at 750 C for "
           "An36 against Giletti & Casserly. Grocolas et al. attribute the older, faster values to SrO-plagioclase "
           "reaction fronts because Sr-feldspar was not stable in those source materials. "
           "Applied to a Santorini plagioclase from Druitt et al. (2012) this gives about 99 kyr "
           "(+62/-42 kyr), and it is the Sr calibration Audetat et al. (2026) recommend. "
           "The Monte Carlo uses the correlated (a, b, Q) covariance that reproduces the "
           "authors' own sampling assumption. See _compensated_covariance. Diffusion runs were "
           "anhydrous, and the effect of water is untested."),
))


def _giletti_casserly1994(dc, cond: Conditions, p):
    """D = D0' * exp(-Q/RT) * 10^(-4.1 X_An), D0' = 8.3176e-5 m2/s."""
    XAn = np.asarray(cond.x("XAn"), dtype=float)
    return (p["D0"] * np.exp(-p["Q"] * 1.0e3 / (R_GAS * cond.T_K))
            * 10.0 ** (p["a"] * XAn))


_add(DiffusionCoefficient(
    key="plag_Sr_giletti_casserly1994",
    mineral="plagioclase", species="Sr",
    label="Plagioclase Sr, Giletti & Casserly (1994)",
    citation="giletti_casserly1994",
    equation_number="(An-dependent form as implemented in the DMG Short Course 2025 script)",
    equation_text=("D_Sr = 8.3176e-5 exp(-276000 / (R T)) x 10^(-4.1 X_An) m2/s. "
                   "Equivalently log D0 = -(4.1 X_An + 4.08)"),
    func=_giletti_casserly1994,
    params={
        "D0": Parameter("D0", 8.3176e-5, 0.0, "m2/s", "1s", "pre-exponential factor at X_An = 0"),
        "Q": Parameter("Q", 276.0, 0.0, "kJ/mol", "1s", "activation energy"),
        "a": Parameter("a", -4.1, 0.0, "-", "1s", "X_An coefficient of log10 D0"),
    },
    sigma_logD=0.3,
    requires=("XAn",),
    needs_fo2=False,
    T_range=Range(973.15, 1373.15, "K (about 700-1100 C)"),
    P_range=Range(1.0e5, 1.0e5, "Pa (1 atm)"),
    X_range=Range(0.2, 0.9, "X_An"),
    verified=True,
    verified_from=("read from the paper PDF: eq. 1 on p. 3791, log10 D0 [m2/s] = -8.18 + "
                   "0.041 [Ab] with Q = 276 kJ/g-atom, and Table 4 (p. 3790). With [Ab] = "
                   "100 (1 - X_An), which neglects Or, this is log10 D0 = -4.08 - 4.1 X_An, "
                   "the form used here. The DMG Short Course 2025 script "
                   "Diff_Model_Sr_in_Plag_implicit.m implements the same expression "
                   "(8.3176e-5*exp(-276000/(temp*8.314))*10^(-4.1*an_conc)). Grocolas et al. "
                   "(2025) eq. 12 re-fit the same data independently as log10 D = -3.76 X_An "
                   "- 4.52 - 270,607/(2.303 R T). The test suite checks that the two agree"),
    secondary_citations=("dmg2025", "cherniak2010", "grocolas2025"),
    recommended=False,
    superseded_by="grocolas2025",
    superseded_note=(
        "Grocolas et al. (2025) measured Sr diffusion in oligoclase and labradorite at "
        "900-1200 C with the silica activity buffered, and found it to be 1.5-2 ORDERS OF "
        "MAGNITUDE SLOWER than this calibration, which they attribute to a SrCl2-plagioclase "
        "reaction in the earlier experiments. Timescales from "
        "this entry are therefore likely to be too short by a factor of roughly 40 at 1100 C "
        "and 100-600 at 750-900 C. At 750 C and An36 the difference reaches 2.8 log units "
        "(Audetat et al. 2026). Use the "
        "plag_Sr_grocolas2025 entry for new work."),
    notes=("Use together with the activity term (A_Sr = -17.4 kJ/mol): because A_Sr is "
           "negative, the equilibrium Sr distribution is inversely correlated with An "
           "content (Dohmen et al. 2017 Appendix, Zellmer et al. 1999, Costa et al. 2003). "
           "Retained so that published timescales built on it can be reproduced."),
))


_add(DiffusionCoefficient(
    key="plag_Sr_cherniak_watson1994",
    mineral="plagioclase", species="Sr",
    label="Plagioclase Sr, Cherniak & Watson (1994), as re-fitted by Grocolas et al. (2025)",
    citation="cherniak_watson1994",
    equation_number="Grocolas et al. (2025) eq. 13",
    equation_text=("log10 D_Sr [m2/s] = -2.30(+/-0.02) X_An - 4.12(+/-0.08) "
                   "- [304,382(+/-1,773) / (2.303 R T)]"),
    func=_grocolas_form,
    params={
        "a": Parameter("a", -2.30, 0.02, "-", "1s", "X_An coefficient"),
        "b": Parameter("b", -4.12, 0.08, "log10(m2/s)", "1s", "intercept"),
        "Q": Parameter("Q", 304.382, 1.773, "kJ/mol", "1s", "activation energy"),
    },
    sigma_logD=0.3,
    requires=("XAn",),
    needs_fo2=False,
    T_range=Range(998.15, 1348.15, "K (725-1075 C)"),
    P_range=Range(1.0e5, 1.0e5, "Pa (1 atm)"),
    X_range=Range(0.23, 0.93, "X_An"),
    verified=True,
    verified_from=("read from Grocolas et al. (2025) eq. 13 (p. 13), their Monte Carlo "
                   "parameterisation of the Cherniak & Watson (1992, 1994) data. Not read from "
                   "the 1994 primary paper"),
    secondary_citations=("grocolas2025",),
    recommended=False,
    superseded_by="grocolas2025",
    superseded_note=("Grocolas et al. (2025) find Sr diffusion 1.5-2 orders of magnitude slower "
                     "than this calibration and attribute the difference to SrO-plagioclase "
                     "reaction fronts, because Sr-feldspar was not stable in the SrO-Al2O3-SiO2 "
                     "source powders. Timescales from it come out about 1.5-2 orders of "
                     "magnitude too short (Grocolas et al. 2025, section 4.3)."),
    notes=("Measured by Rutherford backscattering at 1 atm. Cherniak & Watson (1994) reported "
           "diffusion parallel to b about 0.7 log units slower than parallel to c in some "
           "compositions. The re-fit treats the data as isotropic. This entry replaces an "
           "earlier, unverified interpolation between per-composition fits."),
))


# ---------------------------------------------------------------------------
# Ba
# ---------------------------------------------------------------------------
_add(DiffusionCoefficient(
    key="plag_Ba_grocolas2025",
    mineral="plagioclase", species="Ba",
    label="Plagioclase Ba, Grocolas et al. (2025)",
    citation="grocolas2025",
    equation_number="8",
    equation_text=("log10 D_Ba [m2/s] = -1.43(+/-0.20) X_An - 4.65(+/-0.96) "
                   "- [337,037(+/-22,969) / (2.303 R T)]"),
    func=_grocolas_form,
    params={
        "a": Parameter("a", -1.43, 0.20, "-", "1s", "X_An coefficient"),
        "b": Parameter("b", -4.65, 0.96, "log10(m2/s)", "1s", "intercept"),
        "Q": Parameter("Q", 337.037, 22.969, "kJ/mol", "1s", "activation energy"),
    },
    covariance=_compensated_covariance(0.20, 0.96, 22.969),
    cov_order=("a", "b", "Q"),
    sigma_logD=0.25,
    requires=("XAn",),
    needs_fo2=False,
    T_range=Range(1173.15, 1473.15, "K (900-1200 C)"),
    P_range=Range(1.0e5, 1.0e5, "Pa (1 atm, anhydrous)"),
    X_range=Range(0.28, 0.67, "X_An (oligoclase An28 and labradorite An67)"),
    verified=True,
    verified_from=("read from the paper PDF: abstract and eq. 8 (p. 9). Uncertainties are one "
                   "standard deviation from 2000 Monte Carlo resamplings"),
    secondary_citations=("audetat2026", "cherniak2002"),
    recommended=True,
    notes=("Measured in the same experiments as the Sr law. Unlike Sr, Ba agrees with the "
           "earlier data of Cherniak (2002) to within about 0.5 log units, which Grocolas et al. "
           "explain by Ba-feldspar having been stable in that source powder. Ba diffuses "
           "slightly more slowly than Sr, consistent with its larger ionic radius. Isotropic "
           "within uncertainty. The Monte Carlo covariance follows the authors' own assumption "
           "of perfectly correlated log D0 and Q (see the Sr entry)."),
))

_add(DiffusionCoefficient(
    key="plag_Ba_cherniak2002",
    mineral="plagioclase", species="Ba",
    label="Plagioclase Ba, Cherniak (2002), as re-fitted by Grocolas et al. (2025)",
    citation="cherniak2002",
    equation_number="Grocolas et al. (2025) eq. 14",
    equation_text=("log10 D_Ba [m2/s] = -1.27(+/-0.14) X_An - 5.47(+/-0.44) "
                   "- [331,858(+/-10,077) / (2.303 R T)]"),
    func=_grocolas_form,
    params={
        "a": Parameter("a", -1.27, 0.14, "-", "1s", "X_An coefficient"),
        "b": Parameter("b", -5.47, 0.44, "log10(m2/s)", "1s", "intercept"),
        "Q": Parameter("Q", 331.858, 10.077, "kJ/mol", "1s", "activation energy"),
    },
    sigma_logD=0.3,
    requires=("XAn",),
    needs_fo2=False,
    T_range=Range(1048.15, 1397.15, "K (775-1124 C)"),
    P_range=Range(1.0e5, 1.0e5, "Pa (1 atm)"),
    X_range=Range(0.23, 0.67, "X_An"),
    verified=True,
    verified_from=("read from Grocolas et al. (2025) eq. 14 (p. 13), their Monte Carlo "
                   "parameterisation of the plagioclase data of Cherniak (2002). Not read from "
                   "the 2002 primary paper"),
    secondary_citations=("grocolas2025", "audetat2026"),
    recommended=False,
    notes=("Rutherford backscattering on An23 and An67 at 1 atm. Grocolas et al. (2025) find "
           "their own Ba diffusivities within about 0.5 log units of these, and Audetat et al. "
           "(2026) conclude Ba diffusion in plagioclase is relatively well established, so "
           "this remains a valid alternative to the 2025 law. "
           "Grocolas et al. note it gives timescales about 3 times longer for the Cerro Galan "
           "crystals. Ba is much slower than Mg, which is why Ba zoning survives where Mg has "
           "relaxed (Chamberlain et al. 2014). This entry replaces an earlier, unverified "
           "sanidine relation that Diffusor had listed under plagioclase."),
))


# ---------------------------------------------------------------------------
# Li
# ---------------------------------------------------------------------------
def _pohl2024(dc, cond: Conditions, p):
    return 10.0 ** (p["logD0"] - p["Q"] * 1.0e3 / (LN10 * R_GAS * cond.T_K))


_POHL_NOTE = (
    "Pohl et al. fitted a multispecies model (interstitial and A1-site Li in exchange). "
    "Diffusor applies each mechanism as a single effective diffusion coefficient, which is "
    "an approximation. Chemical diffusion of Li is charge balanced by Na and is 1.5-2 orders "
    "of magnitude slower than the tracer diffusion of Giletti & Shanahan (1997), so earlier "
    "Li timescales may be too short by up to a factor of about 20. Calibrated on An61 only: "
    "a dependence on An content is not yet known. No resolvable pressure effect between 0.1 "
    "and 400 MPa.")

for _mech, _sym, _logD0, _slog, _Q, _sQ, _eq, _rec, _lead in [
        ("interstitial", "i", -3.76, 0.58, 180.0, 12.0, "21", True,
         "The faster of the two mechanisms, and the one that sets the length of a Li profile "
         "(Pohl et al. 2024, section 4.6), so it is the entry to use for timescales from Li "
         "zoning. "),
        ("vacancy", "A", -5.53, 0.16, 151.7, 3.2, "22", False,
         "The slower, vacancy-mediated mechanism, within 0.2-1 orders of magnitude of the "
         "interstitial one over the experimental range. ")]:
    _add(DiffusionCoefficient(
        key=f"plag_Li_pohl2024_{_mech}",
        mineral="plagioclase", species="Li",
        label=f"Plagioclase Li, Pohl et al. (2024) -- {_mech} mechanism",
        citation="pohl2024",
        equation_number=_eq,
        equation_text=(f"D_Li^{_sym} = 10^({_logD0:g} +/- {_slog:g}) "
                       f"exp[-({_Q:g} +/- {_sQ:g} kJ/mol)/(R T)] m2/s"),
        func=_pohl2024,
        params={
            "logD0": Parameter("logD0", _logD0, _slog, "log10(m2/s)", "1s", "pre-exponential factor"),
            "Q": Parameter("Q", _Q, _sQ, "kJ/mol", "1s", "activation energy"),
        },
        sigma_logD=0.3,
        needs_fo2=False,
        T_range=Range(879.15, 1387.15, "K (606-1114 C)"),
        P_range=Range(1.0e5, 4.0e8, "Pa (mostly 200 MPa)"),
        X_range=Range(0.61, 0.61, "X_An (labradorite An61 only)"),
        verified=True,
        verified_from=(f"read from the paper PDF: abstract and eq. {_eq} (p. 997). The "
                       "confidence level of the +/- values is not stated. They are treated as "
                       "1 sigma"),
        secondary_citations=("giletti_shanahan1997",),
        recommended=_rec,
        notes=_lead + _POHL_NOTE,
    ))


# ---------------------------------------------------------------------------
# NaSi-CaAl (the anorthite profile itself)
# ---------------------------------------------------------------------------
def _grove1984(dc, cond: Conditions, p):
    return arrhenius(p["D0"], p["Q"] * 1.0e3, cond.T_K)


_add(DiffusionCoefficient(
    key="plag_NaSiCaAl_grove1984",
    mineral="plagioclase", species="NaSi-CaAl",
    label="Plagioclase coupled NaSi-CaAl interdiffusion, Grove et al. (1984)",
    citation="grove1984",
    equation_number="abstract (p. 2113) and regression on p. 2116",
    equation_text=("D = 10.99 cm2/s exp(-123.4 kcal/mol / R T) = 1.099 x 10^-3 m2/s "
                   "exp(-516.3 kJ/mol / R T), An80-81, 1 atm anhydrous, normal to the (03-1) "
                   "lamellae. Regression error on Q: +/- 5 kcal/mol (20.9 kJ/mol)"),
    func=_grove1984,
    params={
        "D0": Parameter("D0", 1.099e-3, 0.0, "m2/s", "1s", "pre-exponential factor"),
        "Q": Parameter("Q", 123.4 * 4.184, 5.0 * 4.184, "kJ/mol", "1s", "activation energy"),
    },
    sigma_logD=0.7,
    needs_fo2=False,
    T_range=Range(1373.15, 1673.15, "K (1100-1400 C)"),
    P_range=Range(1.0e5, 1.0e5, "Pa (1 atm)"),
    verified=True,
    verified_from=("read from the paper PDF: abstract (p. 2113), p. 2116 and Fig. 3, which the "
                   "line reproduces (ln D = -34.7 at 1400 C, cm2/s). Before 1 October 2026 this "
                   "entry was a placeholder (1.1e-4 m2/s, 520 kJ/mol) about 15 times too slow"),
    notes=("Average CaAl-NaSi interdiffusion coefficient from homogenising exsolution "
           "lamellae in An80-81 Stillwater bytownite, dry, at 1 atm. The paper assumes a "
           "composition-independent binary coefficient and says each D should be viewed as "
           "correct to within a factor of 2 to 5, hence sigma_logD = 0.7. Useful to check "
           "the 'is X_An frozen?' assumption: compare the NaSi-CaAl diffusion length with the "
           "trace-element diffusion length over the fitted time. Hydrous and more sodic "
           "plagioclase differ (compare Liu & Yund 1992)."),
))
