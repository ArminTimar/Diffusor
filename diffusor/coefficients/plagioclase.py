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
from .base import Conditions, DiffusionCoefficient, Parameter, Range, arrhenius

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
    "used A_Mg = -RT ln(gamma) slope of the same form; the DMG Short Course 2025 script "
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
                   "Costa et al. (2003) expression in this form; not read from Costa et al. (2003)"),
    secondary_citations=("vanorman2014", "latourrette_wasserburg1998"),
    notes=("A preliminary estimate built on the compositional dependence of Sr diffusion. "
           "Van Orman et al. (2014) show it over-predicts D at low temperature and low An "
           "(a factor of 10 for albite at 850 C) and note that timescales based on it, "
           "including those of Druitt et al. (2012) for Santorini, are too short by factors "
           "of 2.25-5.5. Kept for reproducing published results, not recommended for new work."),
))


# ---------------------------------------------------------------------------
# Sr
# ---------------------------------------------------------------------------
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
    equation_text=("D_Sr = 8.3176e-5 exp(-276000 / (R T)) x 10^(-4.1 X_An) m2/s; "
                   "equivalently log D0 = -(4.1 X_An + 4.08)"),
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
    verified_from=("transcribed from the DMG Short Course 2025 supplementary script "
                   "Diff_Model_Sr_in_Plag_implicit.m, which implements it as "
                   "8.3176e-5*exp(-276000/(temp*8.314))*10^(-4.1*an_conc); the -4.1 X_An form "
                   "is confirmed by secondary summaries of the paper"),
    secondary_citations=("dmg2025", "cherniak2010"),
    recommended=False,
    superseded_by="grocolas2025",
    superseded_note=(
        "Grocolas et al. (2025) measured Sr diffusion in oligoclase and labradorite at "
        "900-1200 C with the silica activity buffered, and found it to be 1.5-2 ORDERS OF "
        "MAGNITUDE SLOWER than this calibration, which they attribute to feldspar stability "
        "not having been controlled in the earlier experiments. Timescales from this entry "
        "are therefore likely to be too short by a factor of roughly 30-100. Their Arrhenius "
        "parameters are not yet implemented in Diffusor because the paper could not be "
        "retrieved offline; add them to diffusor/coefficients/plagioclase.py from the PDF "
        "(EPSL 651:119141, open access) before using Sr for new work."),
    notes=("Use together with the activity term (A_Sr = -17.4 kJ/mol): because A_Sr is "
           "negative, the equilibrium Sr distribution is inversely correlated with An "
           "content (Dohmen et al. 2017 Appendix; Zellmer et al. 1999; Costa et al. 2003). "
           "Retained so that published timescales built on it can be reproduced."),
))


def _cherniak_watson1994(dc, cond: Conditions, p):
    """log D = logD0 - Q/(ln10 R T) with composition-dependent parameters.

    Cherniak & Watson (1994) fitted each composition separately; Diffusor
    interpolates log D0 and Q linearly in X_An between their An23, An43 and
    An67 fits.
    """
    XAn = np.asarray(cond.x("XAn"), dtype=float)
    xs = np.array([0.23, 0.43, 0.67])
    logD0 = np.array([-6.07, -6.75, -7.03])
    Q = np.array([273.0, 265.0, 268.0])
    l0 = np.interp(XAn, xs, logD0)
    q = np.interp(XAn, xs, Q) * 1.0e3
    from .base import LN10
    return 10.0 ** (l0 - q / (LN10 * R_GAS * cond.T_K))


_add(DiffusionCoefficient(
    key="plag_Sr_cherniak_watson1994",
    mineral="plagioclase", species="Sr",
    label="Plagioclase Sr, Cherniak & Watson (1994) -- per-composition fits",
    citation="cherniak_watson1994",
    equation_number="(Arrhenius relations normal to (001) for An23, An43, An67)",
    equation_text=("oligoclase An23: log D = (-6.07 +/- 0.58) + (-273 +/- 13 kJ/mol)/RT; "
                   "andesine An43: log D = (-6.75 +/- 0.33) + (-265 +/- 8 kJ/mol)/RT; "
                   "labradorite An67: log D = (-7.03 +/- 0.37) + (-268 +/- 8 kJ/mol)/RT "
                   "(D in m2/s, 725-1075 C)"),
    func=_cherniak_watson1994,
    params={},
    sigma_logD=0.3,
    requires=("XAn",),
    needs_fo2=False,
    T_range=Range(998.15, 1348.15, "K (725-1075 C)"),
    X_range=Range(0.23, 0.67, "X_An"),
    verified=False,
    verified_from="values from secondary summaries of the abstract; primary PDF not available offline",
    superseded_by="grocolas2025",
    superseded_note=("Grocolas et al. (2025) find Sr diffusion 1.5-2 orders of magnitude slower "
                     "than the 1990s calibrations; see the note on the Giletti & Casserly entry."),
    notes=("Measured by Rutherford backscattering under anhydrous 1 atm conditions. "
           "Diffusor interpolates log D0 and Q linearly in X_An between the three published "
           "fits and holds them constant outside 0.23-0.67, so the composition dependence is "
           "weaker than in the Giletti & Casserly (1994) parameterisation."),
))


# ---------------------------------------------------------------------------
# Ba
# ---------------------------------------------------------------------------
def _cherniak2002_ba(dc, cond: Conditions, p):
    return arrhenius(p["D0"], p["Q"] * 1.0e3, cond.T_K)


_add(DiffusionCoefficient(
    key="plag_Ba_cherniak2002",
    mineral="plagioclase", species="Ba",
    label="Feldspar Ba, Cherniak (2002) -- alkali feldspar relation",
    citation="cherniak2002",
    equation_number="(abstract)",
    equation_text="D = 2.9 x 10^-1 exp(-455 +/- 20 kJ/mol / R T) m2/s",
    func=_cherniak2002_ba,
    params={
        "D0": Parameter("D0", 2.9e-1, 0.0, "m2/s", "1s", "pre-exponential factor"),
        "Q": Parameter("Q", 455.0, 20.0, "kJ/mol", "1s", "activation energy"),
    },
    sigma_logD=0.3,
    needs_fo2=False,
    T_range=Range(1073.15, 1423.15, "K"),
    verified=False,
    verified_from="values from secondary summaries of the abstract; primary PDF not available offline",
    secondary_citations=("grocolas2025",),
    notes=("Grocolas et al. (2025) report that Ba diffusion in plagioclase is SIMILAR to the "
           "earlier determinations, unlike Sr, so this entry is not thought to be badly wrong. "
           "This is the relation quoted for alkali feldspar (sanidine). Cherniak (2002) also "
           "reports plagioclase compositions; those were not recoverable offline, so use this "
           "entry with care for plagioclase and check the primary paper. Ba is much slower "
           "than Sr, which is why Ba zoning survives where Sr has relaxed "
           "(Chamberlain et al. 2014)."),
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
    equation_number="(cooling-rate speedometry calibration)",
    equation_text="D = 1.1 x 10^-4 exp(-520 kJ/mol / R T) m2/s (order-of-magnitude entry)",
    func=_grove1984,
    params={
        "D0": Parameter("D0", 1.1e-4, 0.0, "m2/s", "1s", "pre-exponential factor"),
        "Q": Parameter("Q", 520.0, 0.0, "kJ/mol", "1s", "activation energy"),
    },
    sigma_logD=0.7,
    needs_fo2=False,
    T_range=Range(1373.15, 1673.15, "K"),
    verified=False,
    verified_from="NOT verified: order-of-magnitude placeholder from secondary summaries",
    notes=("Included only so that the 'is X_An frozen?' assumption can be checked "
           "quantitatively: compare the NaSi-CaAl diffusion length with the trace-element "
           "diffusion length over the fitted time. Do not report timescales from this entry "
           "without first checking the numbers against Grove et al. (1984) and "
           "Liu & Yund (1992)."),
))
