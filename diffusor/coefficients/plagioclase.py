"""Diffusion coefficients for plagioclase.

Trace-element diffusion in plagioclase is coupled to the anorthite gradient:
the driving force is the chemical-potential gradient, not the concentration
gradient, so the flux carries an extra term (Costa et al. 2003 eq. 7; Dohmen,
Faak & Blundy 2017 eq. 6, their Electronic Appendix has the numerical scheme):

    J_i = -D_i dC_i/dx + (D_i C_i / R T) A_i dX_An/dx

with ``A_i`` the slope of ``-R T ln(gamma_i)`` against X_An.  The A values are
in :data:`ACTIVITY_SETS` (two published sets, see below) and are applied by the
solver as ``theta = A_i / (R T)`` at the temperature of each time step.

X_An is treated as frozen.  Dohmen et al. (2017, p. 555) describe the host
anorthite zoning as effectively immobile compared with the trace elements, and
Costa et al. (2003, p. 2193) justify it by noting that Mg diffuses approximately
seven orders of magnitude faster than NaSi-CaAl interdiffusion under dry
conditions (Grove et al. 1984).  The Grove law was measured at 1100-1400 C and
is extrapolated to magmatic temperatures here, and it is a dry law, so use the
``plag_NaSiCaAl_grove1984`` entry to check the assumption for each application
rather than taking it as established.
"""
from __future__ import annotations

import warnings

import numpy as np

from ..constants import R_GAS
from .base import LN10, Conditions, DiffusionCoefficient, Parameter, Range, arrhenius

COEFFICIENTS = []


def _add(c):
    COEFFICIENTS.append(c)
    return c


# --- activity (non-ideality) slopes A_i, kJ/mol -------------------------------
# A_i is defined through the plagioclase-melt partition coefficient,
# R T ln K = A_i X_An + B_i (Dohmen, Faak & Blundy 2017 eq. 7). With a melt activity
# coefficient that does not depend on X_An this is the slope of -R T ln(gamma_i) against
# X_An, the form written as eq. A6 in their Electronic Appendix (the appendix itself was
# not consulted; the sign convention above is the one eq. 6-7 of the main text imply).
# Dohmen, Faak & Blundy (2017) Table 1 (p. 556) lists three sets, read from the
# rendered table on 2 October 2026. They disagree in sign for Mg: Dohmen & Blundy
# (2014) give a positive A_Mg from the lattice-strain model, Bindeman et al.
# (1998) a negative one, and the second was used by Costa et al. (2003) and
# Druitt et al. (2012). Dohmen et al. (2017, p. 556) call the Bindeman values
# potentially incorrect. The set is therefore a choice that every run records.
ACTIVITY_SPECIES = ("Mg", "Sr", "Ba", "Li", "K", "Rb")
ACTIVITY_SETS = {
    "dohmen_blundy2014": {
        "label": "Dohmen & Blundy (2014), lattice-strain model",
        "citation": "dohmen_blundy2014",
        # columns at 900 and 1200 C; the column nearer the run temperature is used
        "columns": {
            1173.15: {"Mg": 13.7, "Sr": -15.1, "Ba": -30.5, "Li": -2.5, "K": -8.8, "Rb": -16.7},
            1473.15: {"Mg": 15.8, "Sr": -17.4, "Ba": -35.1, "Li": -1.7, "K": -8.0, "Rb": -15.9},
        },
    },
    "bindeman1998": {
        "label": "Bindeman et al. (1998), as used by Costa et al. (2003) and Druitt et al. (2012)",
        "citation": "bindeman1998",
        "columns": {
            None: {"Mg": -26.1, "Sr": -30.4, "Ba": -55.0, "Li": -6.9, "K": -25.5, "Rb": -40.0},
        },
    },
}
DEFAULT_ACTIVITY_SET = "dohmen_blundy2014"
ACTIVITY_A_CITATION = "dohmen2017"
ACTIVITY_A_NOTE = (
    "A_i (kJ/mol) from Dohmen, Faak & Blundy (2017) Table 1. Dohmen & Blundy (2014) give values at "
    "900 and 1200 C; Diffusor uses the column nearer the run temperature and does not interpolate. "
    "Bindeman et al. (1998) give one value per element. The two sets differ in sign for Mg.")


def _column(name: str, T_K: float) -> dict:
    try:
        cols = ACTIVITY_SETS[name]["columns"]
    except KeyError as exc:
        raise KeyError(f"unknown activity set {name!r}; choose from {sorted(ACTIVITY_SETS)}") from exc
    if None in cols:
        return cols[None]
    return cols[min(cols, key=lambda T: abs(T - T_K))]


def activity_A(species: str, T_K: float, activity_set: str = DEFAULT_ACTIVITY_SET) -> float:
    """A_i in kJ/mol for the chosen set at temperature T_K.

    An element without a tabulated value gets 0, which means ideal mixing (no
    activity term). That is Diffusor's default, not a published statement that
    the element mixes ideally, so a warning is issued.
    """
    col = _column(activity_set, T_K)
    if species not in col:
        warnings.warn(f"no activity slope A_i is tabulated for {species!r} in the activity set "
                      f"{activity_set!r}; using 0 (ideal mixing), which is an assumption, not a "
                      "published value", UserWarning, stacklevel=2)
    return float(col.get(species, 0.0))


def activity_note(activity_set: str, T_K: float) -> str:
    s = ACTIVITY_SETS[activity_set]
    cols = s["columns"]
    if None in cols:
        return f"A_i from {s['label']}."
    T_col = min(cols, key=lambda T: abs(T - T_K))
    return f"A_i from {s['label']}, the {T_col - 273.15:.0f} C column (nearest the run temperature)."


def activity_theta(species: str, T_K: float, activity_set: str = DEFAULT_ACTIVITY_SET) -> float:
    """theta = A_i / (R T), the coefficient of the dX_An/dx flux term."""
    return activity_A(species, T_K, activity_set) * 1.0e3 / (R_GAS * T_K)


def equilibrium_profile(X_An, T_K: float, species: str, C_ref: float, X_An_ref: float = None,
                        activity_set: str = DEFAULT_ACTIVITY_SET):
    """Quasi-steady-state ("equilibrated") profile for a frozen X_An(x).

    ``C_eq(x) = C0 exp(A_i X_An(x) / (R T))``  -- Dohmen et al. (2017) Appendix
    eq. A13, with C0 fixed by the rim condition eq. A14 when ``X_An_ref`` is
    given (C_ref is then the rim concentration).

    This is the state a fast trace element relaxes to. It is the final state of
    diffusion, so it is an initial condition only for a crystal that had already
    equilibrated with one melt before the event being timed.
    """
    X_An = np.asarray(X_An, dtype=float)
    th = activity_theta(species, T_K, activity_set)
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
    kind="chemical",
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
    T_range=Range(1073.15, 1423.15, "K (800-1150 C, the range of the new experiments)"),
    P_range=Range(1.0e5, 101325.0, "Pa (atmospheric; 1 bar to 1 atm)"),
    X_range=Range(0.23, 0.95, "x_An"),
    verified=True,
    verified_from="read from the paper PDF: abstract and eq. 4 on p. 84",
    recommended=True,
    notes=("Fitted to anorthite (An93), labradorite (An67), andesine (An43) and oligoclase "
           "(An23) plus the An95 data of LaTourrette & Wasserburg (1998), which were measured "
           "at 1200-1400 C, so the fit itself extends to 1400 C for the most calcic end. Little "
           "anisotropy was found between the b and c directions (the anisotropy experiments are "
           "called preliminary and cover labradorite and, from the earlier data, anorthite), and "
           "the authors write that it is probably adequate for most applications to treat Mg "
           "diffusion in plagioclase as isotropic (section 3.3). Reduced chi-squared of the fit "
           "was 2.88. The paper labels the labradorite An66 in the experimental section and An67 "
           "elsewhere, and its Table 2 labels the An43 and An23 rows 'Olig' and 'And', the "
           "reverse of the text; Diffusor follows the abstract and text (An23 oligoclase, An43 "
           "andesine). The paper states that the uncertainties of eq. 4 are 2 sigma. The +/-10 kJ/mol "
           "on Q is printed in the abstract; the typeset eq. 4 on p. 84 leaves the number out."),
))


def _costa2003_mg(dc, cond: Conditions, p):
    XAn = np.asarray(cond.x("XAn"), dtype=float)
    lnD = p["lnD0"] + p["a"] * XAn - p["Q"] * 1.0e3 / (R_GAS * cond.T_K)
    return np.exp(lnD)


_add(DiffusionCoefficient(
    key="plag_Mg_costa2003",
    kind="effective",
    mineral="plagioclase", species="Mg",
    label="Plagioclase Mg, Costa et al. (2003) -- older estimate",
    citation="costa2003",
    equation_number="eq. 8, p. 2193 (as re-written by Van Orman et al. 2014, p. 84)",
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
    T_range=Range(1073.15, 1473.15, "K (800-1200 C; assumed by Diffusor, Costa et al. state no calibration range)"),
    X_range=Range(0.0, 1.0, "x_An (full range of the formula; assumed by Diffusor, see notes)"),
    verified=True,
    verified_from=("Costa et al. (2003) eq. 8 (p. 2193), D = 2.92 x 10^(-4.1 X_An - 3.1) "
                   "exp(-266000/RT) m2/s, read from the rendered page on 2 October 2026. It equals the "
                   "ln D = -6.07 - 9.44 X_An - 266/RT form of Van Orman et al. (2014, p. 84) and eq. 1 "
                   "of Druitt et al. (2012)"),
    secondary_citations=("vanorman2014", "latourrette_wasserburg1998"),
    notes=("A preliminary estimate: Costa et al. (2003, p. 2193) cite the Mg data for "
           "anorthite of LaTourrette & Wasserburg (1998) and estimate the albite end by assuming a "
           "compositional dependence similar to that of Sr; they say the estimate does not remove "
           "the need for direct measurements. It is not a fit to Mg experiments, so the paper gives no "
           "temperature range, no composition range and no uncertainty; the 800-1200 C and "
           "x_An 0-1 ranges and the 'effective' transport kind are Diffusor's own assumptions "
           "(the law is an estimate applied outside measured conditions, and the An-gradient "
           "term of Costa et al. eq. 7 is separate from D). "
           "Van Orman et al. (2014, p. 84, not Costa et al.) report that it over-predicts D at "
           "low temperature and low An (a factor of 10 for albite at 850 C) and that the "
           "expression as used by Druitt et al. (2012) over-estimates Mg diffusivities by a "
           "factor of 2.25-5.5 relative to their results, so Santorini timescales may need to be "
           "revised upward by a similar factor. Kept for reproducing published results. Avoid it "
           "for new work."),
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
    kind="chemical",
    mineral="plagioclase", species="Mg",
    label="Plagioclase Mg, Audetat et al. (2026) -- with silica activity",
    citation="audetat2026",
    equation_number="1",
    equation_text=("log10 D_Mg [m2/s] = -2.99(+/-0.35) X_An - 4.03(+/-0.63) "
                   "- [262,914(+/-15,529) / (2.303 R T)] - 1.87(+/-0.45)(1 - aSiO2)"),
    func=_audetat2026_mg,
    params={
        "a": Parameter("a", -2.99, 0.35, "-", "unstated", "X_An coefficient"),
        "b": Parameter("b", -4.03, 0.63, "log10(m2/s)", "unstated", "intercept"),
        "Q": Parameter("Q", 262.914, 15.529, "kJ/mol", "unstated", "activation energy"),
        "c": Parameter("c", -1.87, 0.45, "-", "unstated", "coefficient of (1 - aSiO2)"),
        "aSiO2": Parameter("aSiO2", 1.0, 0.0, "-", "1s",
                           "silica activity relative to quartz. 1 = quartz saturated"),
    },
    sigma_logD=0.25,
    requires=("XAn",),
    needs_fo2=False,
    T_range=Range(1073.15, 1473.15, "K (800-1200 C, the combined range of the two data sets)"),
    P_range=Range(1.0e5, 101325.0, "Pa (atmospheric; 1 bar to 1 atm)"),
    X_range=Range(0.23, 0.93, "X_An"),
    verified=True,
    verified_from=("read from the accepted-manuscript PDF of Audetat et al. (2026), eq. 1 "
                   "(manuscript lines 271-278). The paper does not say what the +/- values are "
                   "(standard deviation, standard error or a confidence bound), so the level is "
                   "recorded as not stated and the values are sampled as if they were 1 sigma. "
                   "The paper prints 2.303 for ln 10; Diffusor uses the exact ln 10 = 2.302585, "
                   "which lowers log D by 0.002 (1473 K) to 0.0023 (1073 K) log units"),
    secondary_citations=("faak2013", "vanorman2014"),
    recommended=True,
    notes=("Parameterised from the Faak et al. (2013) and Van Orman et al. (2014) experiments "
           "(LaTourrette & Wasserburg 1998 is excluded because aSiO2 of their glass is hard to "
           "estimate). Mg diffusion is faster at higher silica activity (Faak et al. 2013). "
           "Diffusor evaluates it at aSiO2 = 1 by default, which corresponds to quartz "
           "saturation (aSiO2 is defined relative to quartz). For a silica-undersaturated melt D "
           "is lower by up to 1.87 log units at aSiO2 = 0. The review recommends this law, "
           "together with the Mg partitioning model of Mutch et al. (2022) (p. 47), and notes "
           "that Mg-in-plagioclase equilibration times from three studies fit within a factor "
           "of two with Fe-Mg interdiffusion times in olivine from the same samples (p. 11). "
           "The numbers come from the accepted manuscript. The paper does not classify the "
           "coefficient as tracer or chemical; the 'chemical' kind is Diffusor's choice. Van Orman "
           "et al. (2014, p. 2) note that Mg chemical diffusion rates are the ones relevant to "
           "the relaxation of Mg zoning, and find that D does not depend significantly on the "
           "Mg concentration gradient (abstract). "
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
    matrix is Diffusor's encoding of that sentence: the paper prints no
    covariance matrix, no slope of the log D0-Ea trend and no sign. The positive
    sign (a larger Ea goes with a larger log D0, compensation) is Diffusor's
    choice. With perfect correlation the spread of log D from b and Q vanishes at
    T = s_Q / (ln10 R s_b), 949 C for Sr and 977 C for Ba, and the X_An slope a is
    independent.
    """
    return np.array([[s_a ** 2, 0.0, 0.0],
                     [0.0, s_b ** 2, s_b * s_Q],
                     [0.0, s_b * s_Q, s_Q ** 2]])


_add(DiffusionCoefficient(
    key="plag_Sr_grocolas2025",
    kind="chemical",
    mineral="plagioclase", species="Sr",
    label="Plagioclase Sr, Grocolas et al. (2025)",
    citation="grocolas2025",
    equation_number="7",
    equation_text=("log10 D_Sr [m2/s] = -1.65(+/-0.24) X_An - 3.03(+/-1.16) "
                   "- [368,142(+/-27,141) / (2.303 R T)]"),
    func=_grocolas_form,
    params={
        "a": Parameter("a", -1.65, 0.24, "-", "unstated", "X_An coefficient"),
        "b": Parameter("b", -3.03, 1.16, "log10(m2/s)", "unstated", "intercept"),
        "Q": Parameter("Q", 368.142, 27.141, "kJ/mol", "unstated", "activation energy"),
    },
    covariance=_compensated_covariance(0.24, 1.16, 27.141),
    cov_order=("a", "b", "Q"),
    sigma_logD=0.25,
    requires=("XAn",),
    needs_fo2=False,
    T_range=Range(1173.15, 1423.15, "K (900-1150 C, the runs used in the fit; the experiments "
                  "span 900-1200 C but the degraded 1200 C run and the 1100 C run OHSC_8 were "
                  "excluded)"),
    P_range=Range(1.0e5, 101325.0, "Pa (atmospheric, anhydrous; 1 bar to 1 atm)"),
    X_range=Range(0.28, 0.67, "X_An (oligoclase An28 and labradorite An67)"),
    verified=True,
    verified_from=("read from the paper PDF: abstract and eq. 7 (p. 9). The paper does not say "
                   "what the +/- printed in eq. 7 are; the Fig. 7 caption calls the fitted "
                   "envelopes one standard deviation from 2000 Monte Carlo resamplings. The level "
                   "is therefore recorded as not stated and the values are sampled as 1 sigma. "
                   "The paper prints 2.303 for ln 10; Diffusor uses the exact ln 10 = 2.302585, "
                   "which lowers log D by about 0.003 log units over 900-1200 C"),
    secondary_citations=("audetat2026",),
    recommended=True,
    notes=("Sr in-diffusion into oriented oligoclase and labradorite with the silica activity "
           "buffered by sol-gel mullite-cristobalite or mullite-corundum powders. No resolvable "
           "dependence on aSiO2 or crystal orientation, so the law is isotropic. The values are "
           "1-2 orders of magnitude slower than Giletti & Casserly (1994) and Cherniak & "
           "Watson (1994) (the paper says ~1-2 log units on p. 10 and 1.5-2 orders in the "
           "abstract), and because the activation energy is higher the gap widens as "
           "temperature falls. Audetat et al. (2026, p. 9) give 2.8 log units at 750 C and An36 "
           "against Giletti & Casserly; Diffusor's own calculation from eq. 7 and the Giletti & "
           "Casserly law as coded here gives about 1.6 log units at 1100 C, 2.2 at 900 C and 2.8 "
           "at 750 C for An36 (750 C is below the 900 C lower limit of eq. 7, so that value is an "
           "extrapolation). Grocolas et al. (pp. 10-11) attribute the older, faster values to "
           "SrO-plagioclase reaction fronts because Sr-feldspar was not stable in those source "
           "materials, and to a reaction with the evaporated Sr chloride in the Giletti & Casserly "
           "experiments. Applied to a Santorini plagioclase from Druitt et al. (2012) this gives "
           "99.4 kyr (+61.8/-41.5 kyr) (Grocolas et al., p. 13). Audetat et al. (2026, p. 47) "
           "recommend this Sr calibration and method for plagioclase, while also stating that a "
           "consensus on Sr in plagioclase has not yet been reached. "
           "The paper does not classify the coefficient as chemical or tracer; 'chemical' is "
           "Diffusor's classification of in-diffusion from a source powder. "
           "The Monte Carlo uses a correlated (a, b, Q) covariance that is Diffusor's encoding "
           "of the authors' own sampling assumption (log D0 and Ea on a linear trend without an "
           "envelope, section 4.3, p. 13); the paper prints no covariance matrix, and the authors "
           "note that the assumption slightly underestimates the total uncertainty. See "
           "_compensated_covariance. Diffusion runs were anhydrous, and the effect of water is "
           "untested."),
))


def _giletti_casserly1994(dc, cond: Conditions, p):
    """D = D0' * exp(-Q/RT) * 10^(-4.1 X_An), D0' = 8.3176e-5 m2/s."""
    XAn = np.asarray(cond.x("XAn"), dtype=float)
    return (p["D0"] * np.exp(-p["Q"] * 1.0e3 / (R_GAS * cond.T_K))
            * 10.0 ** (p["a"] * XAn))


_add(DiffusionCoefficient(
    key="plag_Sr_giletti_casserly1994",
    kind="tracer",
    mineral="plagioclase", species="Sr",
    label="Plagioclase Sr, Giletti & Casserly (1994)",
    citation="giletti_casserly1994",
    equation_number="eq. 1 (p. 3791) with Q = 276 kJ/g-atom; [Ab] = 100 (1 - X_An)",
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
    T_range=Range(823.15, 1573.15, "K (550-1300 C, the span of the experiments in the four "
                  "specimens; albite 550-1080, oligoclase 750-1100, labradorite 800-1300, "
                  "anorthite 900-1300 C, so only 900-1080 C is covered by all four)"),
    P_range=Range(1.0e5, 101325.0, "Pa (atmospheric, in air; 1 bar to 1 atm)"),
    X_range=Range(0.006, 0.956, "X_An (the four specimens, An0.6 to An95.6; the paper says "
                  "the relation can be used for any normal plagioclase composition)"),
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
        "reaction in the earlier experiments (Grocolas et al. 2025, p. 11). Timescales from "
        "this entry are therefore likely to be too short; Diffusor's own calculation for An36 "
        "from the two laws gives a factor of roughly 40 at 1100 C and 150-600 at 750-900 C. "
        "At 750 C and An36 the difference reaches 2.8 log units (Audetat et al. 2026, p. 9). "
        "Use the plag_Sr_grocolas2025 entry for new work."),
    notes=("Giletti & Casserly (1994) say nothing about an activity term. The Sr "
           "chemical-potential term of Dohmen et al. (2017, eq. 6; Table 1 gives A_Sr = -15.1 "
           "kJ/mol at 900 C and -17.4 kJ/mol at 1200 C in the default set) is a later "
           "construction: with a negative A_Sr, eq. 6 gives an equilibrium Sr distribution "
           "that falls with increasing An content. The temperature and composition ranges "
           "are the experimental spans stated in the paper; the paper warns that using the "
           "relation beyond the measured temperatures is done at the user's risk and notes "
           "that the four specimens avoided the peristerite, Boggild and Huttenlocher "
           "intergrowth ranges, so interpolation between them is not tested (p. 3791). "
           "Retained so that published timescales built on it can be reproduced."),
))


_add(DiffusionCoefficient(
    key="plag_Sr_cherniak_watson1994",
    kind="chemical",
    mineral="plagioclase", species="Sr",
    label="Plagioclase Sr, Cherniak & Watson (1994), as re-fitted by Grocolas et al. (2025)",
    citation="cherniak_watson1994",
    equation_number="Grocolas et al. (2025) eq. 13",
    equation_text=("log10 D_Sr [m2/s] = -2.30(+/-0.02) X_An - 4.12(+/-0.08) "
                   "- [304,382(+/-1,773) / (2.303 R T)]"),
    func=_grocolas_form,
    params={
        "a": Parameter("a", -2.30, 0.02, "-", "unstated", "X_An coefficient"),
        "b": Parameter("b", -4.12, 0.08, "log10(m2/s)", "unstated", "intercept"),
        "Q": Parameter("Q", 304.382, 1.773, "kJ/mol", "unstated", "activation energy"),
    },
    sigma_logD=0.3,
    requires=("XAn",),
    needs_fo2=False,
    T_range=Range(998.15, 1348.15, "K (725-1075 C, the range of Cherniak & Watson 1994)"),
    P_range=Range(1.0e5, 101325.0, "Pa (atmospheric, anhydrous; 1 bar to 1 atm)"),
    X_range=Range(0.23, 0.67, "X_An (An23 to An67, the compositions of Cherniak & Watson 1994)"),
    verified=True,
    verified_from=("read from Grocolas et al. (2025) eq. 13 (p. 13), their Monte Carlo "
                   "parameterisation of the Cherniak & Watson (1992, 1994) data. Grocolas et al. "
                   "do not state the temperature, pressure or composition range of that "
                   "data set; the ranges here are those of the abstract of Cherniak & Watson "
                   "(1994): 725-1075 C, anhydrous 1 atm, three compositions between An23 and "
                   "An67. Grocolas et al. print the +/- of eq. 13 without a confidence level. "
                   "The 1992 paper, which Grocolas et al. also cite for this fit, was not "
                   "available; an upper X_An of 0.93 used earlier could not be confirmed and "
                   "was removed"),
    secondary_citations=("grocolas2025",),
    recommended=False,
    superseded_by="grocolas2025",
    superseded_note=("Grocolas et al. (2025) find Sr diffusion 1.5-2 orders of magnitude slower "
                     "than this calibration and attribute the difference to SrO-plagioclase "
                     "reaction fronts, because Sr-feldspar was not stable in the SrO-Al2O3-SiO2 "
                     "source powders. Timescales from it come out about 1.5-2 orders of "
                     "magnitude too short (Grocolas et al. 2025, section 4.3)."),
    notes=("Measured by Rutherford backscattering at 1 atm, anhydrous (Cherniak & Watson 1994 "
           "abstract, which also calls the quantity Sr chemical diffusion). That paper "
           "reported diffusion normal to (010) in labradorite about 0.7 log units slower than "
           "normal to (001), and none in andesine; Grocolas et al. (p. 8) repeat the 0.7 "
           "for labradorite and no anisotropy for oligoclase, while p. 10 of the same paper says "
           "about 0.5 log units. Eq. 13 has a single set of coefficients, so the re-fit is "
           "isotropic by construction (Grocolas et al. do not state which orientations went into "
           "it). This entry replaces an earlier, unverified interpolation between "
           "per-composition fits."),
))


# ---------------------------------------------------------------------------
# Ba
# ---------------------------------------------------------------------------
_add(DiffusionCoefficient(
    key="plag_Ba_grocolas2025",
    kind="chemical",
    mineral="plagioclase", species="Ba",
    label="Plagioclase Ba, Grocolas et al. (2025)",
    citation="grocolas2025",
    equation_number="8",
    equation_text=("log10 D_Ba [m2/s] = -1.43(+/-0.20) X_An - 4.65(+/-0.96) "
                   "- [337,037(+/-22,969) / (2.303 R T)]"),
    func=_grocolas_form,
    params={
        "a": Parameter("a", -1.43, 0.20, "-", "unstated", "X_An coefficient"),
        "b": Parameter("b", -4.65, 0.96, "log10(m2/s)", "unstated", "intercept"),
        "Q": Parameter("Q", 337.037, 22.969, "kJ/mol", "unstated", "activation energy"),
    },
    covariance=_compensated_covariance(0.20, 0.96, 22.969),
    cov_order=("a", "b", "Q"),
    sigma_logD=0.25,
    requires=("XAn",),
    needs_fo2=False,
    T_range=Range(1173.15, 1423.15, "K (900-1150 C, the runs used in the fit; the experiments "
                  "span 900-1200 C but the degraded 1200 C run and the 1100 C run OHSC_8 were "
                  "excluded)"),
    P_range=Range(1.0e5, 101325.0, "Pa (atmospheric, anhydrous; 1 bar to 1 atm)"),
    X_range=Range(0.28, 0.67, "X_An (oligoclase An28 and labradorite An67)"),
    verified=True,
    verified_from=("read from the paper PDF: abstract and eq. 8 (p. 9). The paper does not say "
                   "what the +/- printed in eq. 8 are; the Fig. 7 caption calls the fitted "
                   "envelopes one standard deviation from 2000 Monte Carlo resamplings. The level "
                   "is therefore recorded as not stated and the values are sampled as 1 sigma. "
                   "The paper prints 2.303 for ln 10; Diffusor uses the exact ln 10 = 2.302585, "
                   "which lowers log D by about 0.003 log units over 900-1200 C"),
    secondary_citations=("audetat2026", "cherniak2002"),
    recommended=True,
    notes=("Measured in the same experiments as the Sr law. Unlike Sr, Ba agrees with the "
           "earlier data of Cherniak (2002) to within about 0.5 log units, which Grocolas et al. "
           "explain by Ba-feldspar having been stable in that source powder. Ba diffuses "
           "slightly more slowly than Sr, consistent with its larger ionic radius (Audetat et "
           "al. 2026, p. 11). Isotropic within uncertainty. The paper does not classify the "
           "coefficient as chemical or tracer; 'chemical' is Diffusor's classification. The "
           "Monte Carlo covariance is Diffusor's encoding of the authors' own assumption of "
           "perfectly correlated log D0 and Q (see the Sr entry); the paper prints no "
           "covariance matrix."),
))

_add(DiffusionCoefficient(
    key="plag_Ba_cherniak2002",
    kind="chemical",
    mineral="plagioclase", species="Ba",
    label="Plagioclase Ba, Cherniak (2002), as re-fitted by Grocolas et al. (2025)",
    citation="cherniak2002",
    equation_number="Grocolas et al. (2025) eq. 14",
    equation_text=("log10 D_Ba [m2/s] = -1.27(+/-0.14) X_An - 5.47(+/-0.44) "
                   "- [331,858(+/-10,077) / (2.303 R T)]"),
    func=_grocolas_form,
    params={
        "a": Parameter("a", -1.27, 0.14, "-", "unstated", "X_An coefficient"),
        "b": Parameter("b", -5.47, 0.44, "log10(m2/s)", "unstated", "intercept"),
        "Q": Parameter("Q", 331.858, 10.077, "kJ/mol", "unstated", "activation energy"),
    },
    sigma_logD=0.3,
    requires=("XAn",),
    needs_fo2=False,
    T_range=Range(1048.15, 1397.15, "K (775-1124 C, the plagioclase runs of Cherniak 2002)"),
    P_range=Range(1.0e5, 101325.0, "Pa (atmospheric, dry; 1 bar to 1 atm)"),
    X_range=Range(0.23, 0.67, "X_An (An23 and An67, the compositions of Cherniak 2002)"),
    verified=True,
    verified_from=("read from Grocolas et al. (2025) eq. 14 (p. 13), their Monte Carlo "
                   "parameterisation of the plagioclase data of Cherniak (2002). Grocolas et al. "
                   "do not state the temperature, pressure or composition range of that data "
                   "set. The ranges here come from Cherniak (2002): Table 1 plagioclase runs "
                   "775 C (oligoclase) to 1124 C (labradorite), An23 and An67, dry 1 atm "
                   "(abstract). Which of the runs went into eq. 14 is not stated. Grocolas et "
                   "al. print the +/- of eq. 14 without a confidence level"),
    secondary_citations=("grocolas2025", "audetat2026"),
    recommended=False,
    notes=("Rutherford backscattering on An23 and An67 at 1 atm, in air (Cherniak 2002 "
           "abstract and Grocolas et al. 2025, p. 10). Grocolas et al. (2025, p. 10) find "
           "their own Ba diffusivities similar to slightly faster (about 0.5 log units) than "
           "these, and Audetat et al. (2026, pp. 10 and 47) conclude from the good match that "
           "Ba diffusion in plagioclase is relatively well established and that results similar "
           "to the Grocolas et al. calibration are obtained with this one; 'a valid "
           "alternative' is Diffusor's wording. "
           "Grocolas et al. (p. 13) note it gives timescales about 3 times longer for the Cerro "
           "Galan crystals. Cherniak (2002, p. 1644) reports anisotropy in oligoclase, with "
           "diffusion normal to (010) slower and a larger activation energy, and Audetat et al. "
           "(2026, p. 9) quote up to about 0.6 orders of magnitude; none in labradorite. The "
           "re-fit of eq. 14 has one set of coefficients, so Diffusor treats it as isotropic. "
           "This entry replaces an earlier, unverified sanidine relation that Diffusor had "
           "listed under plagioclase."),
))


# ---------------------------------------------------------------------------
# Li
# ---------------------------------------------------------------------------
def _pohl2024(dc, cond: Conditions, p):
    return 10.0 ** (p["logD0"] - p["Q"] * 1.0e3 / (LN10 * R_GAS * cond.T_K))


_POHL_NOTE = (
    "Pohl et al. fitted a multispecies model (interstitial and A1-site Li in exchange). "
    "Diffusor applies each mechanism as a single effective diffusion coefficient, which is "
    "an approximation; the paper calls the Li coefficients chemical diffusion (interdiffusion "
    "with Na), and the 'effective' kind is Diffusor's label for a mechanism-specific "
    "coefficient taken out of a multispecies model. Chemical diffusion of Li is charge "
    "balanced by Na and is 1.5-2 orders of magnitude slower than the tracer diffusion of "
    "Giletti & Shanahan (1997); the authors conclude that using that tracer coefficient "
    "underestimates the time of their in-diffusion experiments by a factor of 20 (p. 1001; "
    "1.5-2 orders would be a factor of 30-100, so the paper is not consistent with itself), "
    "and note that Li typically diffuses out of crystals during magma ascent, a case not "
    "tested. Calibrated on An61 only, on crystals cut with the polished surface perpendicular "
    "to [001] (section 2.2), so a dependence on An content or on direction is not known. "
    "At 750 C the authors do not consider the slight positive trend of D with pressure "
    "(1 atm, 50, 200 and 400 MPa) a resolvable pressure effect (section 4.10).")

for _mech, _sym, _logD0, _slog, _Q, _sQ, _eq, _rec, _lead in [
        ("interstitial", "i", -3.76, 0.58, 180.0, 12.0, "21", True,
         "The faster of the two mechanisms, and the one that sets the length of a Li profile "
         "(Pohl et al. 2024, section 4.5, p. 997; the authors suggest using it with the "
         "profile length for natural samples in section 4.11, p. 1001), so it is the entry "
         "to use for timescales from Li zoning. "),
        ("vacancy", "A", -5.53, 0.16, 151.7, 3.2, "22", False,
         "The slower, vacancy-mediated mechanism, within 0.2-1 orders of magnitude of the "
         "interstitial one over the experimental range. ")]:
    _add(DiffusionCoefficient(
        key=f"plag_Li_pohl2024_{_mech}",
        kind="effective",
        mineral="plagioclase", species="Li",
        label=f"Plagioclase Li, Pohl et al. (2024) -- {_mech} mechanism",
        citation="pohl2024",
        equation_number=_eq,
        equation_text=(f"D_Li^{_sym} = 10^({_logD0:g} +/- {_slog:g}) "
                       f"exp[-({_Q:g} +/- {_sQ:g} kJ/mol)/(R T)] m2/s"),
        func=_pohl2024,
        params={
            "logD0": Parameter("logD0", _logD0, _slog, "log10(m2/s)", "unstated", "pre-exponential factor"),
            "Q": Parameter("Q", _Q, _sQ, "kJ/mol", "unstated", "activation energy"),
        },
        sigma_logD=0.3,
        needs_fo2=False,
        T_range=Range(879.15, 1387.15, "K (606-1114 C, the experimental span; the model "
                      "misfits the steep crystal gradient at 1000-1114 C and the 700-750 C "
                      "runs were fitted to the isotope data only)"),
        P_range=Range(1.0e5, 4.0e8, "Pa (mostly 200 MPa)"),
        X_range=Range(0.61, 0.61, "X_An (labradorite An61 only)"),
        verified=True,
        verified_from=(f"read from the paper PDF: abstract and eq. {_eq} (p. 997). The "
                       "paper says only that the equations were determined by linear "
                       "regression and does not state what the +/- values are, so the level is "
                       "recorded as not stated and the values are sampled as 1 sigma"),
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
    kind="interdiffusion",
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
        "Q": Parameter("Q", 123.4 * 4.184, 5.0 * 4.184, "kJ/mol", "unstated",
                       "activation energy (the +/- 5 kcal/mol is the regression error; level not stated)"),
    },
    sigma_logD=0.7,
    needs_fo2=False,
    T_range=Range(1373.15, 1673.15, "K (1100-1400 C)"),
    P_range=Range(1.0e5, 101325.0, "Pa (atmospheric; 1 bar to 1 atm)"),
    verified=True,
    verified_from=("read from the paper PDF: abstract (p. 2113), p. 2116 and Fig. 3, which the "
                   "line reproduces (ln D = -34.7 at 1400 C, cm2/s). Before 1 October 2026 this "
                   "entry was a placeholder (1.1e-4 m2/s, 520 kJ/mol) about 15 times too slow"),
    notes=("Average CaAl-NaSi interdiffusion coefficient from homogenising exsolution "
           "lamellae in An80-81 Stillwater bytownite, dry, at 1 atm. The paper assumes a "
           "composition-independent binary coefficient and says each D should be viewed as "
           "correct to within a factor of 2 to 5; sigma_logD = 0.7 is log10(5), the upper end "
           "of that range, applied by Diffusor as a conservative 1 sigma (the paper's own "
           "estimate for each D is about +/-0.5 in ln D, 0.22 log10 units). Useful to check "
           "the 'is X_An frozen?' assumption: compare the NaSi-CaAl diffusion length with the "
           "trace-element diffusion length over the fitted time. Grove et al. (p. 2119) say "
           "that the effect of water fugacity on the rate is unknown but could be important, "
           "and that they have no knowledge of a concentration dependence of D; "
           "Cherniak & Watson (2020, p. 1044) add that CaAl-NaSi interdiffusion may be "
           "enhanced by hydrous species (Yund 1986; Yund & Snow 1989; Liu & Yund 1992). This "
           "coefficient is for dry, calcic plagioclase only."),
))
