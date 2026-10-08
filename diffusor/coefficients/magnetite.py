"""Diffusion coefficients for magnetite and titanomagnetite.

Cation diffusion in magnetite is unusual: a vacancy mechanism (D proportional
to fO2^(2/3)) competes with an interstitial(cy) mechanism (D proportional to
fO2^(-2/3)), so D passes through a *minimum* with respect to both fO2 and
temperature.  Van Orman & Crispin (2010) RiMG 72, Table 12 tabulate the two
branches:

    D* = D_V,0 exp(-Q_V / R T) a_O2^(2/3) + D_I,0 exp(-Q_I / R T) a_O2^(-2/3)

with ``a_O2 = fO2 / fO2_0`` and ``fO2_0 = 1 atm``.  Note that the vacancy
activation energies Q_V are *negative*, which is what produces the minimum.

This is the formulation Tomiya et al. (2013) used for the 2011 Shinmoedake
eruption.  Diffusor reproduces their published Ti and Fe values to within about
5 % when the two branches are interpolated separately in x_Ti (see the notes of
the Table 12 entries and the test suite).
"""
from __future__ import annotations

import numpy as np

from ..constants import R_GAS, BAR_TO_PA
from .base import Conditions, DiffusionCoefficient, Parameter, Range, arrhenius

COEFFICIENTS = []
ATM_IN_BAR = 1.01325     # 1 atm = 1.01325 bar (SI definition)


def _add(c):
    COEFFICIENTS.append(c)
    return c


# --- Van Orman & Crispin (2010) Table 12 --------------------------------------
# species -> (D_V0 [m2/s], Q_V [kJ/mol], D_I0 [m2/s], Q_I [kJ/mol]) for
# pure magnetite (x_Ti = 0) and for x_Ti = 0.2.
TABLE12_PURE = {
    "Cr": (5.12e-13, -7.3, 3.84e5, 752.5),
    "Al": (3.24e-13, -64.9, 6.92e3, 681.0),   # D_I0 printed "6.92x10" in Table 12; see below
    "Fe": (1.68e-14, -123.1, 9.79e3, 618.2),
    "Co": (2.20e-15, -144.8, 8.22e3, 612.0),
    "Mn": (2.02e-16, -176.7, 2.48e3, 604.6),
    "Ti": (1.37e-14, -107.6, 2.96e4, 688.7),
}
TABLE12_XTI02 = {
    "Fe": (1.36e-13, -167.5, 2.33e1, 562.3),
    "Co": (5.55e-13, -148.0, 1.47e2, 584.8),
    "Mn": (5.64e-13, -149.9, 2.47e1, 568.3),
    "Ti": (4.07e-13, -130.4, 3.20e9, 840.4),
}


# Temperature windows (K) of the data behind each law.  Table 12 itself prints no
# temperature range.  For Fe, Co, Mn and Ti the laws are fits "from all data
# reported in Table 5 of Aggarwal & Dieckmann (2002b)" (Table 12 footnote), and
# that table is now in hand: it lists the normalised partial coefficients
# D_Me[V] and D_Me[I] at 1100, 1200, 1300 and 1400 C.  Fitting ln D = ln D0 - Q/RT
# to those values reproduces the Table 12 rows (D0 within 1 %, Q within 0.2 kJ/mol,
# i.e. within the three printed digits; for x_Ti = 0 the literature entries of the
# same table, marked a, b and c, are part of the fit).
# The temperatures of the values used are therefore the real basis of the laws:
#   x_Ti = 0   Fe: 1100-1400 C (1373-1673 K)
#              Co, Mn: 1100-1300 C (1373-1573 K); these are literature values that
#                      Aggarwal & Dieckmann tabulate (Co partly flagged interpolated
#                      or extrapolated by them; Mn from Franke & Dieckmann, pers.
#                      comm., and Lu et al. 1993), not their own measurements
#              Ti: 1200-1400 C (1473-1673 K); no 1100 C value exists for pure
#                      magnetite (Table 10 prints 1373-1573 K for Ti, which the
#                      data of Table 5 do not support)
#   x_Ti = 0.2 Fe, Co, Mn, Ti: 1100-1300 C (1373-1573 K), measured in that paper
# A law that interpolates between the two compositions is limited to the overlap
# of the two windows, because the interpolation needs both end members.  Cr and
# Al are not Aggarwal & Dieckmann results: Table 12 cites Dieckmann, Hilton & Mason
# (1987), whose Table 1 (p. 61) gives the Cr data at 1210-1410 C and the Al values
# at 1280-1500 C, the windows printed in Table 10.
#
# Al interstitial D0.  Table 12 prints D_I,0 = "6.92x10" m2/s, the only entry without
# a visible exponent.  Dieckmann et al. (1987, p. 61, summarizing expressions, read
# from the rendered page on 8 October 2026) print D0_Al[I] = 6.92e7 exp(-81900 K/T)
# cm2/s, i.e. 6.92e3 m2/s and 681.0 kJ/mol; their other three Cr and Al expressions
# convert exactly to the Table 12 values.  Diffusor uses 6.92e3 (6.92e1 until
# 8 October 2026, which put the interstitial branch 2 log units low).
TABLE12_T_WINDOWS = {
    #        pure magnetite              x_Ti = 0.2
    "Cr": ((1483.0, 1683.0), None),
    "Al": ((1553.0, 1773.0), None),
    "Fe": ((1373.0, 1673.0), (1373.0, 1573.0)),
    "Co": ((1373.0, 1573.0), (1373.0, 1573.0)),
    "Mn": ((1373.0, 1573.0), (1373.0, 1573.0)),
    "Ti": ((1473.0, 1673.0), (1373.0, 1573.0)),
}


def _table12_T_range(species: str) -> Range:
    pure, ti02 = TABLE12_T_WINDOWS[species]
    lo, hi = pure
    text = f"pure magnetite {pure[0]:.0f}-{pure[1]:.0f} K"
    if ti02 is not None:
        lo, hi = max(lo, ti02[0]), min(hi, ti02[1])
        text += f", x_Ti = 0.2 {ti02[0]:.0f}-{ti02[1]:.0f} K; overlap used"
    basis = ("Table 5 of Aggarwal & Dieckmann 2002b, from which the Table 12 rows can be refitted"
             if species in ("Fe", "Co", "Mn", "Ti")
             else "Table 1 of Dieckmann, Hilton & Mason 1987, the data this row is fitted to")
    return Range(lo, hi, f"K (temperatures of the data behind the law: {text}; {basis}. "
                         "Table 12 prints no range, so the law is an extrapolation outside it)")


def _table12_branches(species: str, T_K: float, a_O2, table):
    """Vacancy and interstitial terms (m2/s) of one Table 12 row."""
    DV0, QV, DI0, QI = table[species]
    vac = DV0 * np.exp(-QV * 1.0e3 / (R_GAS * T_K)) * a_O2 ** (2.0 / 3.0)
    inter = DI0 * np.exp(-QI * 1.0e3 / (R_GAS * T_K)) * a_O2 ** (-2.0 / 3.0)
    return vac, inter


def _table12_D(species: str, T_K: float, a_O2, table) -> float:
    vac, inter = _table12_branches(species, T_K, a_O2, table)
    return vac + inter


# --- Composition dependence of the Table 12 laws (x_Ti between 0 and 0.2) -----------
# What the sources say.  Table 12 gives parameters at x_Ti = 0 and 0.2 only, and no
# source gives D for the intermediate compositions as a function of x_Ti:
#  * Aggarwal & Dieckmann (2002b, Phys Chem Minerals 29, 707-718) measured the
#    tracer coefficients at x = 0, 0.1, 0.2 and 0.3 (their x is Ti per cation site, the
#    x of (Ti_x Fe_1-x)_3 O_4, i.e. x_Ti here) and tabulate the partial coefficients
#    D*_Me[V] and D*_Me[I] separately for each x in Table 5 (p. 713): 1200 and
#    1300 C at x = 0, 0.1, 0.2 and 0.3, 1100 C at x = 0.2 and 1400 C at x = 0.  They
#    show them against x in Figs. 6-7 and fit no function of x.  Their Part I
#    (2002a, Phys Chem Minerals 29, 695-706, eqs. 17-19, p. 703) fits the normalised
#    defect concentrations as log10 [def] = A1 + A2 x - (C1 + C2 x)/T, i.e.
#    log-linear in x with an enthalpy linear in x; Van Orman & Crispin (2010,
#    p. 794) conclude that Ti enhances vacancy diffusion almost entirely through the
#    larger vacancy concentration.  Aragon et al. (1984, p. 177) find log D of the
#    Fe-Ti interdiffusion linear in x and relate it to the exponential increase of
#    the vacancy concentration with x.
#  * Van Orman & Crispin give no rule for intermediate x_Ti.  Tomiya et al. (2013,
#    p. 14) say only that they used the Van Orman & Crispin equation at X_Usp = 0.3
#    and print the four values reproduced in the Ti and Fe notes.
# Diffusor's choice, therefore, is its own: the vacancy and the interstitial branch,
# D = D_V,0 exp(-Q_V/RT) a_O2^(2/3) + D_I,0 exp(-Q_I/RT) a_O2^(-2/3), are each
# interpolated log-linearly between the two printed rows (ln D0 and Q of each
# branch linear in x_Ti, so ln D_branch is linear in x_Ti at every T and fO2).
# Why branch-wise: the two terms are separate mechanisms with opposite fO2
# exponents that respond oppositely to Ti, each is an Arrhenius law of its own, and
# the rows are two points on a log-linear curve of the defect-concentration type
# above.  The assumption made is that for each mechanism ln(D) is linear in x_Ti at
# fixed T and fO2 (the mobility of the defect being independent of x_Ti, or
# varying log-linearly).  Evidence and limits, all against Aggarwal & Dieckmann
# (2002b) Tables 2, 3 and 5 (checked in tests/test_magnetite_sources.py):
#  * At x_Ti = 0 and 0.2 the Table 12 rows reproduce their 65 raw tracer values at
#    1200 and 1300 C to an rms of 0.13 log units (0.11 with the Part I correction
#    C_V, which Table 12 omits and which is at most 0.1 log units here).
#  * At x_Ti = 0.1 (1200, 1300 C) the branch-wise interpolation is within 0.3 log
#    units of their vacancy coefficients and 0.26-0.94 too low for the interstitial
#    ones (their D_I is flat in x up to 0.1 and falls afterwards; Ti is not even
#    monotonic).  On their 39 raw values at x_Ti = 0.1 it has mean -0.09, rms
#    0.30 and maximum 0.76 log units; interpolating the total log D instead gives
#    +0.12, 0.29 and 0.52.  The data do not decide between the two.
#  * At x_Ti = 0.3 (extrapolation, fraction 1.5; 35 raw values) rms 0.34 with mean
#    +0.31 for branch-wise and 0.33 / +0.22 for total-D; the measured minimum of D
#    in fO2 is reproduced to 0.0-0.9 log units by the branch-wise form and 0.9-1.7
#    by the total-D form.  At x_Ti = 0.1 both put the minimum at the same fO2, 0.3-0.8
#    log units below the measured one; the depth of the minimum is wrong by 0.04-0.15
#    log units (Fe) for the branch-wise form but 0.35-0.45 (too fast) for the total-D
#    form, whereas for Ti the total-D form is closer (0.17-0.20 against 0.20-0.49).
#  * The two forms differ by up to 0.5 log units for Fe, Co and Mn (x_Ti = 0.1,
#    1300 C, log fO2 near -6) and 0.8 for Ti (850 C, log fO2 -16), so the choice
#    matters little compared with the 0.3 log unit scatter above.
#  * Tomiya et al.'s four values (x_Ti = 0.1) are reproduced within about 5 % by the
#    branch-wise form; total-D reproduces Ti within 1 %, Fe at 900 C within 10 %
#    and Fe at 950 C only within a factor 2.1.  Linear interpolation of the printed
#    parameters (D0 and Q themselves) is 0.2-0.5 log units off and is not used.
#  * Limits: no measurement at intermediate x_Ti exists below 1200 C, where natural
#    magnetite lives (800-1000 C), and Aggarwal & Dieckmann's activation energies
#    at x = 0.1 and 0.3 rest on two temperatures only (error unknown, some
#    nonphysical, e.g. Q_V(Ti, x = 0.1) = +25 kJ/mol in their Table 6).  Extrapolating
#    their x = 0.1 values with those energies to 950 C gives D 0.35 (Fe) to 1.35 (Ti)
#    log units below this interpolation, an indication of the true uncertainty of
#    the intermediate compositions at magmatic temperature.
# Node scheme, evaluated and NOT adopted.  The measured coefficients can be turned into
# nodes at x_Ti = 0.1 and 0.3 by the very fit that reproduces the Table 12 rows,
# ln D = ln D0 - Q/RT through the Table 5 values, and each branch interpolated
# log-linearly between adjacent nodes (0, 0.1, 0.2, 0.3), the form that Part I
# eqs. 17-19 give for the defect concentrations.  Table 5 has x = 0.1 and 0.3 at 1200
# and 1300 C only (not at 1100 or 1400 C, and not for Cr or Al), so every node is a
# two-temperature fit, (D_V,0 m2/s, Q_V kJ/mol, D_I,0 m2/s, Q_I kJ/mol):
#   x = 0.1  Fe 1.68e-13 -124.7 4.81e5 667.7   Co 6.87e-13 -105.6 3.52e4 631.8
#            Mn 2.28e-13 -120.7 1.80e5 659.3   Ti 4.89e-9 +25.3 2.25e11 867.8
#   x = 0.3  Fe 2.26e-23 -493.8 2.90e4 670.4   Co 1.98e-24 -523.5 1.40e4 661.5
#            Mn 7.59e-23 -477.8 2.55e5 705.0   Ti 1.47e-23 -478.3 9.31e7 813.3
# (the energies agree with A&D's Table 6 to 0.3 kJ/mol except Co at x = 0.3, where they
# print -409.7; Table 6 marks them "errors unknown").
#  * In sample the nodes reproduce the 39 (x = 0.1) and 35 (x = 0.3) raw values with an
#    rms of 0.06 and 0.07 log units, against 0.30 and 0.34 for the 0/0.2 interpolation.
#    That is circular: the nodes are those data, and it says nothing about other T.
#  * Away from 1200-1300 C the nodes rest on two temperatures and on energies that A&D
#    call unreliable: Q_V(Ti, x = 0.1) is +25 kJ/mol between -108 (x = 0) and -130 (x = 0.2)
#    and Q_V at x = 0.3 is -410 to -520 kJ/mol, "suspicious" in their words (pp. 712-713).
#    At 900-1000 C and log fO2 -11 to -13 the nodes lie 0.57 below to 0.14 above the 0/0.2
#    laws for Fe, Co, Mn and 1.6 below to 0.1 above for Ti at x = 0.1, and 1.9-4.0 log
#    units above them at x = 0.3.  They also miss Tomiya et al.'s values (Ti by a factor
#    22-41, Fe by 0.35-0.38 log units), which the 0/0.2 scheme matches.
#  * Out of sample in temperature (rescale the 0/0.2 branches to the measured
#    coefficient at 1200 C and predict the 1300 C values, and the reverse) the measured
#    amplitude helps at x = 0.1 (rms 0.17-0.18 against 0.26-0.34 for both unanchored
#    forms) and does not at x = 0.3 (0.27-0.55 against 0.32-0.34).
#  Decision: the 0/0.2 branch-wise scheme stays the default, because the node scheme can
#  only be checked on the data it was made from, its temperature dependence at
#  x = 0.1 and 0.3 is unconstrained and in part unphysical, and it changes D at magmatic
#  temperature by up to a factor of 40 relative to the laws that Tomiya et al. used.  The
#  measured coefficients at 1200 and 1300 C are the better description of x_Ti = 0.1 and
#  0.3 in that window; a user working there should expect the 0/0.2 laws to be off by
#  the 0.3 log units given above.
def _make_table12_func(species: str):
    def _f(dc, cond: Conditions, p):
        """Vacancy + interstitial sum, each branch log-interpolated in x_Ti."""
        a_O2 = 10.0 ** (cond.log_fo2_bar - np.log10(ATM_IN_BAR))   # fO2 in atm
        if species not in TABLE12_XTI02:
            return _table12_D(species, cond.T_K, a_O2, TABLE12_PURE)
        xTi = np.asarray(cond.X.get("xTi", 0.0), dtype=float)
        v0, i0 = _table12_branches(species, cond.T_K, a_O2, TABLE12_PURE)
        v1, i1 = _table12_branches(species, cond.T_K, a_O2, TABLE12_XTI02)
        # Composition dependence: Diffusor's own construction (see the block above
        # _make_table12_func for the evidence and its limits).  Each branch is
        # interpolated log-linearly in x_Ti (ln D0 and Q of each branch linear in
        # x_Ti), then the branches are summed; the fraction is capped at 1.5
        # (x_Ti = 0.3).
        f = np.clip(xTi / 0.2, 0.0, 1.5)
        vac = np.exp((1.0 - f) * np.log(v0) + f * np.log(v1))
        inter = np.exp((1.0 - f) * np.log(i0) + f * np.log(i1))
        return vac + inter
    return _f


_COMPOSITION_NOTE = (
    "Table 12 has rows for x_Ti = 0 and x_Ti = 0.2 only, and no source gives D between them as a "
    "function of x_Ti, so the composition dependence is Diffusor's own construction: the vacancy "
    "branch and the interstitial branch are each interpolated log-linearly in x_Ti (ln D0 and Q of "
    "each branch linear in x_Ti), the branches are summed afterwards, and the interpolation "
    "fraction is capped at 1.5 (x_Ti = 0.3); the code warns outside x_Ti = 0-0.2. Aggarwal & "
    "Dieckmann (2002b) also measured both branches at x_Ti = 0.1 and 0.3 (their Table 5, 1200 and "
    "1300 C) but fit no function of x_Ti; for the defect concentrations their Part I (2002a, "
    "eqs. 17-19) uses the same log-linear form, and Van Orman & Crispin (2010, p. 794) attribute "
    "the faster vacancy diffusion with Ti almost entirely to the larger vacancy concentration. "
    "Against the Aggarwal & Dieckmann values at x_Ti = 0.1 the interpolation is within 0.3 log "
    "units for the vacancy branch but 0.26-0.94 log units too low for the interstitial branch, "
    "which is flat in x_Ti up to 0.1 and falls afterwards (Ti is not monotonic). On their 39 "
    "tracer values at x_Ti = 0.1 the rms error is 0.30 log units (largest 0.76), the same as "
    "interpolating the total D log-linearly (0.29, largest 0.52), so those data do not decide "
    "between the two forms; the branch-wise form is kept because each branch is a separate "
    "mechanism with its own fO2 exponent and Arrhenius law, because it places the minimum of D "
    "nearer the measured one at x_Ti = 0.3 and because it reproduces Tomiya et al.'s values. The "
    "two forms differ by up to 0.5 log units (Fe, Co, Mn) and 0.8 (Ti, reducing, 850 C). At "
    "x_Ti = 0.3 (extrapolation) the rms error is 0.34 log units and D is over-predicted by 0.3 on "
    "average. No measurement at an intermediate x_Ti exists below 1200 C, so at magmatic "
    "temperatures the intermediate compositions are an extrapolation that can be wrong by more "
    "than the 0.3 log units found at 1200-1300 C. Using their Table 5 values at x_Ti = 0.1 and 0.3 "
    "as extra interpolation nodes was tested and not adopted: each node has only two "
    "temperatures (1200 and 1300 C), A&D mark the resulting activation energies as of unknown "
    "error (Q_V of Ti at x_Ti = 0.1 is positive, Q_V at x_Ti = 0.3 is -410 to -520 kJ/mol), the "
    "nodes fit their own data with an rms of 0.06-0.07 only because they are those data, and at "
    "900-1000 C they differ from the laws used here by 0.6 log units (Fe, Co, Mn) to 1.6 (Ti) at "
    "x_Ti = 0.1 and by 2-4 at x_Ti = 0.3, and miss Tomiya et al.'s Ti values by a factor 22-41.")


def _table12_notes(sp: str) -> str:
    """Species-specific notes of the Table 12 entries (numbers from this implementation)."""
    cites = {
        "Ti": "Aggarwal & Dieckmann (2002b)", "Fe": "Aggarwal & Dieckmann (2002b)",
        "Mn": "Aggarwal & Dieckmann (2002b)", "Co": "Aggarwal & Dieckmann (2002b)",
        "Cr": "Dieckmann et al. (1987)", "Al": "Dieckmann et al. (1987)"}
    out = [f"Table 12 cites {cites[sp]} for this species."]
    basis = {
        "Fe": "the pure-magnetite row on 1100-1400 C values (Aggarwal & Dieckmann's own at 1200-1400 C "
              "and literature values they tabulate)",
        "Co": "the pure-magnetite row on 1100-1300 C literature values they tabulate (partly flagged "
              "by them as interpolated or extrapolated)",
        "Mn": "the pure-magnetite row on 1100-1300 C values they take from Franke & Dieckmann "
              "(personal communication) and Lu et al. (1993)",
        "Ti": "the pure-magnetite row on their own 1200-1400 C values only (no 1100 C value for pure "
              "magnetite exists in that table, although Table 10 prints 1373-1573 K for Ti)"}
    if sp in basis:
        out.append(
            "Refitting ln D = ln D0 - Q/RT to the partial coefficients D_V and D_I of Aggarwal & "
            "Dieckmann (2002b, Table 5) reproduces the Table 12 rows of this species (D0 within 1 %, "
            f"Q within 0.2 kJ/mol), so the temperature window used here is that of those values: {basis[sp]}; "
            "the x_Ti = 0.2 row on their own 1100-1300 C measurements.")
    if sp == "Co":
        out.append("The pure-magnetite Co tracer data of Table 10 (1179-1483 K) are credited there "
                   "to Dieckmann et al. (1978).")
    if sp == "Cr":
        out.append("The text (p. 794) credits the Cr tracer data to Dieckmann et al. (1978) and "
                   "Hodge (1978), and Table 10 to Dieckmann et al. (1978). The row is the "
                   "summarizing expression of Dieckmann, Hilton & Mason (1987, p. 61) converted "
                   "from cm2/s and K to m2/s and kJ/mol, fitted to their Table 1 values at "
                   "1210-1410 C.")
    if sp == "Al":
        out.append("Unlike the other rows, the Al law is not a radiotracer result: the text "
                   "(p. 794) describes it as a re-analysis by Dieckmann et al. (1987) of "
                   "interdiffusion data of Petuskey (1977), extrapolated to pure magnetite "
                   "(Table 10 gives the Al window 1553-1773 K). Table 12 writes it as D* like the others. "
                   "Table 12 prints the interstitial D0 as '6.92x10' m2/s with no visible exponent; "
                   "Dieckmann et al. (1987, p. 61) give 6.92e7 cm2/s = 6.92e3 m2/s, which Diffusor "
                   "uses. With it the row reproduces their Table 1 values at 1280-1500 C within "
                   "0.06 log units.")
    if sp in TABLE12_XTI02:
        out.append(_COMPOSITION_NOTE)
    else:
        out.append("Table 12 has no x_Ti = 0.2 row for this species, so the law is for pure "
                   "magnetite and x_Ti is ignored.")
    if sp == "Ti":
        out.append(
            "Tomiya et al. (2013) evaluated this law at 950 C, log fO2 = -11 and X_Usp = 0.3 "
            "(x_Ti = X_Usp/3 = 0.1) and print 4.3e-16 m2/s for Ti (6.9e-16 at 900 C). Diffusor gives "
            "4.15e-16 and 6.83e-16 (3 % and 1 % lower). Tomiya et al. do not state the fO2 unit or "
            "how they obtained the x_Ti dependence. Solving for the x_Ti that reproduces each of "
            "their four published Ti and Fe values gives 0.100-0.102 with the branch-wise "
            "interpolation used here; if the total D is interpolated instead it gives 0.100 for Ti "
            "but 0.066 (950 C) and 0.097 (900 C) for Fe. That favours the branch-wise form only "
            "through one of the four numbers, which is weak evidence. "
            "Tomiya et al. took D_Al = D_Ti, so use this entry for Al to follow them; "
            "the mt_Al entry is a separate pure-magnetite law about 50 times slower there. "
            "They report a minimum of D(T) at about 980 C for Ti and conclude that the error from "
            "the temperature is below a factor of 2. This implementation, at log fO2 = -11 and "
            "x_Ti = 0.1, puts the Ti minimum at 1023 C (a Diffusor calculation; the minimum is "
            "flat, D is within a factor of 2 of it from about 930 to 1070 C).")
    elif sp == "Fe":
        out.append(
            "Tomiya et al. (2013) took D_Mg = D_Fe and give 4.4e-15 m2/s at 950 C and 6.6e-15 m2/s "
            "at 900 C (log fO2 = -11, X_Usp = 0.3, i.e. x_Ti = 0.1, with the Van Orman & Crispin "
            "equation); Diffusor gives 4.17e-15 and 6.54e-15 (5 % and 1 % lower). They place the "
            "Fe minimum of D(T) at about 950 C and conclude that the error from the temperature "
            "is below a factor of 2. This implementation, at log fO2 = -11 and x_Ti = 0.1, puts the "
            "Fe minimum at 955 C (a Diffusor calculation), with D within a factor of 2 of the "
            "minimum from about 880 to 1000 C.")
    elif sp in ("Co", "Mn"):
        out.append("Tomiya et al. (2013) did not use this species.")
    out.append("The minimum of D with temperature (negative Q_V) is the mechanism described by Van "
               "Orman & Crispin (2010, p. 794); its position depends on fO2 and composition.")
    return " ".join(out)


for _sp in ("Ti", "Fe", "Mn", "Co", "Cr", "Al"):
    _has02 = _sp in TABLE12_XTI02
    _add(DiffusionCoefficient(
        key=f"mt_{_sp}_vanorman_crispin2010",
        mineral="magnetite", species=_sp,
        kind="chemical" if _sp == "Al" else "tracer",
        transported_variable=("Al concentration (interdiffusion with Fe, extrapolated to pure magnetite)"
                              if _sp == "Al" else f"{_sp} tracer concentration"),
        label=(f"Magnetite Al diffusion (interdiffusion extrapolated to pure magnetite), "
               "Van Orman & Crispin (2010) Table 12" if _sp == "Al" else
               f"Magnetite {_sp} tracer diffusion, Van Orman & Crispin (2010) Table 12"),
        citation="vanorman_crispin2010",
        equation_number="Table 12",
        equation_text=("D* = D_V,0 exp(-Q_V/RT) a_O2^(2/3) + D_I,0 exp(-Q_I/RT) a_O2^(-2/3), "
                       "a_O2 = fO2 / fO2_0 with fO2_0 = 1 atm. "
                       f"{_sp} (x_Ti = 0): D_V,0 = {TABLE12_PURE[_sp][0]:g} m2/s, "
                       f"Q_V = {TABLE12_PURE[_sp][1]:g} kJ/mol, "
                       f"D_I,0 = {TABLE12_PURE[_sp][2]:g} m2/s, Q_I = {TABLE12_PURE[_sp][3]:g} kJ/mol"
                       + (f". {_sp} (x_Ti = 0.2): D_V,0 = {TABLE12_XTI02[_sp][0]:g} m2/s, "
                          f"Q_V = {TABLE12_XTI02[_sp][1]:g} kJ/mol, "
                          f"D_I,0 = {TABLE12_XTI02[_sp][2]:g} m2/s, "
                          f"Q_I = {TABLE12_XTI02[_sp][3]:g} kJ/mol" if _has02 else "")),
        func=_make_table12_func(_sp),
        params={},
        sigma_logD=0.3,
        requires=("xTi",) if _has02 else (),
        needs_fo2=True, fo2_unit="atm",
        T_range=_table12_T_range(_sp),
        X_range=(Range(0.0, 0.2, "x_Ti (Ti per cation site, the x of (TixFe1-x)3O4; the Table 12 rows are "
                                 "x_Ti = 0 and 0.2, Aggarwal & Dieckmann also measured 0.1 and 0.3 at "
                                 "1200 and 1300 C. X_Ti = X_Usp / 3 is Diffusor's conversion)")
                 if _has02 else
                 Range(0.0, 0.0, "x_Ti (pure magnetite only: Table 12 has no x_Ti = 0.2 row for "
                                 f"{_sp}, and the law ignores x_Ti)")),
        verified=True,
        verified_from=("read from the PDF of Van Orman & Crispin (2010) RiMG 72, Table 12 "
                       "(p. 821) and its footnote giving the vacancy/interstitial sum"
                       + ("; the rows were also reproduced by refitting the partial coefficients of "
                          "Aggarwal & Dieckmann (2002b) Table 5 (p. 713) read from the rendered page"
                          if _sp in ("Fe", "Co", "Mn", "Ti") else
                          "; the row was checked against the summarizing expressions and Table 1 of "
                          "Dieckmann, Hilton & Mason (1987, p. 61) and Table 1 of Dieckmann et al. "
                          "(1978, p. 779), read from the rendered pages on 8 October 2026")),
        secondary_citations=(("tomiya2013", "aggarwal_dieckmann2002", "aggarwal_dieckmann2002a")
                             if _sp in ("Fe", "Co", "Mn", "Ti")
                             else ("dieckmann1987", "dieckmann1978", "tomiya2013")),
        recommended=(_sp in ("Ti", "Fe")),
        notes=_table12_notes(_sp),
    ))


# --- Sievwright et al. (2020): modern data at 1150 C ----------------------------
# Table 5: log D_V1 and log D_I1 (m2/s) of the fit to eq. 5,
#     D = D_V1 a_O2^(2/3) + D_I1 a_O2^(-2/3),   T = 1150 C.
# The last two numbers are the published location of the minimum (log fO2, log D),
# kept so the test suite can check the transcription.
#
# fO2 unit.  Eq. 5 is written with the oxygen activity a_O2 and the paper does not
# state its standard state.  The Fig. 6 axes are labelled "log fO2 (atm)", the
# experiments were run at 1 bar total pressure, and the paper's own dFMQ values were
# calculated with O'Neill (1987), whose reference pressure is 1 bar.  Diffusor takes
# a_O2 = fO2 / (1 bar), i.e. the log10 fO2 in bar.  Reading the axis as atm instead
# (1 atm = 1.01325 bar) lowers log fO2 by 0.0057 and log D by at most 0.004.  The
# Van Orman & Crispin (2010) Table 12 entries state fO2_0 = 1 atm and use atm.
SIEVWRIGHT_T_K = 1423.15
# log10 fO2 (bar) of the experiments: FMQ-1 to FMQ+4.89 (Table 2), dFMQ from O'Neill
# (1987) as in the paper, at 1150 C and 1 bar: -10.017 to -4.127 (computed with
# diffusor.thermo.buffers.log_fo2_from_delta(..., parameterisation="oneill")).  The
# Frost (1991) FMQ used elsewhere in Diffusor would give -9.90 to -4.01.  Rounded
# outward to two decimals.  1423.15 K is 3 K above the 1420 K upper limit printed for
# the O'Neill (1987) expression, so that evaluation issues a BufferRangeWarning.
SIEVWRIGHT_FO2_RANGE = (-10.02, -4.12)
SIEVWRIGHT_TABLE5 = {
    "Mn": (-9.07, -18.8, -7.22, -13.61),
    "Co": (-9.13, -18.4, -6.94, -13.47),
    "Ni": (-9.68, -18.5, -6.58, -13.79),
    "Mg": (-9.27, -18.5, -6.87, -13.58),
    "Zn": (-9.18, -18.8, -7.20, -13.71),
    "Sc": (-9.60, -20.2, -7.91, -14.60),
    "Al": (-10.3, -21.0, -8.01, -15.33),
    "Ga": (-9.48, -20.3, -8.07, -14.58),
    "In": (-8.97, -20.1, -8.28, -14.22),
    "Y": (-9.04, -19.6, -7.89, -14.03),
    "Cr": (-12.0, -21.6, -7.18, -16.47),
    "Lu": (-9.15, -19.6, -7.81, -14.07),
    "V3+": (-10.5, -21.2, -7.94, -15.56),
    "Ti": (-10.1, -21.3, -8.40, -15.41),
    "V4+": (-10.5, -21.2, -7.94, -15.56),
    "Zr": (-9.17, -20.6, -8.56, -14.60),
    "Hf": (-9.86, -21.1, -8.40, -15.18),
    "U": (-8.23, -19.8, -8.63, -13.71),
    "Nb": (-9.50, -21.4, -8.90, -15.17),
    "Ta": (-10.3, -22.0, -8.77, -15.87),
    "Mo": (-10.3, -22.2, -8.89, -15.91),
}


def sievwright_D_1150(species: str, log_fo2_bar):
    """D (m2/s) at 1150 C from Sievwright et al. (2020) Table 5 and eq. 5."""
    lv, li = SIEVWRIGHT_TABLE5[species][:2]
    lf = np.asarray(log_fo2_bar, dtype=float)
    return 10.0 ** (lv + lf * 2.0 / 3.0) + 10.0 ** (li - lf * 2.0 / 3.0)


def _make_sievwright_func(species: str, scale_with_table12: bool):
    def _f(dc, cond: Conditions, p):
        lf = cond.log_fo2_bar
        vac = 10.0 ** (SIEVWRIGHT_TABLE5[species][0] + lf * 2.0 / 3.0)
        inter = 10.0 ** (SIEVWRIGHT_TABLE5[species][1] - lf * 2.0 / 3.0)
        if scale_with_table12:
            _, QV, _, QI = TABLE12_PURE[species]
            dinv = 1.0 / cond.T_K - 1.0 / SIEVWRIGHT_T_K
            vac = vac * np.exp(-QV * 1.0e3 / R_GAS * dinv)
            inter = inter * np.exp(-QI * 1.0e3 / R_GAS * dinv)
        return vac + inter
    return _f


for _sp in SIEVWRIGHT_TABLE5:
    _scaled = _sp in TABLE12_PURE
    _lv, _li, _lfmin, _ldmin = SIEVWRIGHT_TABLE5[_sp]
    _add(DiffusionCoefficient(
        key=f"mt_{_sp}_sievwright2020",
        mineral="magnetite", species=_sp,
        kind="effective" if _scaled else "chemical",
        transported_variable=f"{_sp} concentration",
        reference_state="1150 C, magnetite equilibrated with silicate melt, 1 bar",
        fixed_temperature_K=None if _scaled else SIEVWRIGHT_T_K,
        calibration_notes=(("HYPOTHESIS: temperature dependence borrowed from different tracer experiments. "
                            "Sievwright et al. measured only 1150 C; use the 1150 C-only entry for the original law.",)
                           if _scaled else ("Single-temperature law: only 1150 C is supported.",)),
        uncertainty_note="The 0.2 log10 D sampling width is an assumed representative error, not a published fit covariance.",
        label=(f"Magnetite {_sp}, Sievwright et al. (2020) at 1150 C"
               + (", T-scaled with Table 12" if _scaled else ", 1150 C ONLY")),
        citation="sievwright2020",
        equation_number="5 and Table 5",
        equation_text=(f"D = D_V1 a_O2^(2/3) + D_I1 a_O2^(-2/3) (eq. 5 with D_V1 = D_V,0 exp(-H_V/RT) "
                       f"and D_I1 = D_I,0 exp(-H_I/RT), Table 5 note d), at 1150 C, a_O2 taken as fO2 in bar. "
                       f"{_sp}: log D_V1 = {_lv:g}, log D_I1 = {_li:g} (m2/s). Minimum "
                       f"log D = {_ldmin:g} at log fO2 = {_lfmin:g}"
                       + (". Away from 1150 C each branch is scaled by "
                          "exp[-Q/R (1/T - 1/1423.15 K)] with the Q_V and Q_I of Van Orman "
                          "& Crispin (2010) Table 12 for pure magnetite" if _scaled else "")),
        func=_make_sievwright_func(_sp, _scaled),
        params={},
        sigma_logD=0.2,
        needs_fo2=True, fo2_unit="bar",
        T_range=Range(1423.15, 1423.15, "K (1150 C only; other temperatures are hypothetical)"),
        P_range=Range(1.0e5, 101325.0, "Pa (atmospheric; 1 bar to 1 atm)"),
        fo2_range=Range(*SIEVWRIGHT_FO2_RANGE,
                        "log10 bar (FMQ-1 to FMQ+4.89 at 1150 C, FMQ from O'Neill 1987 as in the paper)"),
        verified=True,
        verified_from=("read from the paper PDF: eq. 5 (p. 12) and Table 5 (p. 13). The "
                       "minimum of eq. 5 computed from the rounded Table 5 constants reproduces the "
                       "published minimum log D of every element to within 0.04 log units and the "
                       "published minimum log fO2 to within 0.09 (the printed constants are rounded, "
                       "so the two cannot agree exactly)"
                       + (". Diffusor added the temperature scaling. It is absent from the "
                          "paper" if _scaled else "")),
        secondary_citations=(("vanorman_crispin2010", "aggarwal_dieckmann2002")
                             if _scaled else ()),
        recommended=False,
        notes=(("Natural magnetite equilibrated with a silicate melt at 1 bar and FMQ-1 to "
                "FMQ+4.89, measured by LA-ICP-MS. Uncertainties on individual log D values are "
                "typically below 0.2 log units (1 sigma); the paper's Supplementary Table S1 lists "
                "the 1 sigma of each (median 0.12, 84 % of the 164 values that have one at or "
                "below 0.2, largest 0.51), and the fits of Table 5 were weighted by them but "
                "carry no published uncertainty. The paper calls its coefficients chemical "
                "diffusion coefficients. fO2 convention: eq. 5 is written with the activity a_O2 "
                "without a stated standard state; Fig. 6 labels its axis 'log fO2 (atm)', but the "
                "runs were at 1 bar and dFMQ comes from O'Neill (1987) (reference pressure 1 bar), so "
                "Diffusor takes fO2 in bar (the difference to atm is 0.0057 in log fO2, at most 0.004 "
                "in log D). ")
               + ("Sievwright et al. give no activation energy because all runs were at "
                  "1150 C, and do not scale their data to other temperatures. They compare their "
                  "Mn and Ti results with literature data (Dieckmann et al. 1978, Dieckmann & "
                  "Schmalzried 1986, Aggarwal & Dieckmann 2002) predicted for 1150 C by assuming "
                  "the temperature dependence of Fe; their Supplementary Fig. S3 plots the 1200 C "
                  "curves of Fe, Co, Ni, Mn, Ti, Al and Cr, the same curves shifted to 1150 C and "
                  "their Mn and Ti data, and the text calls the agreement very good. Rebuilt here "
                  "(the 1200 C Mn and Ti curves from the partial coefficients of Aggarwal & "
                  "Dieckmann's Table 5 for pure magnetite, shifted by the temperature dependence of "
                  "the Table 12 Fe law as a stand-in for the Fe law the authors used), the six 12 h "
                  "experiments lie within 0.43 log units of the Mn curve (mean +0.26) and 0.76 of "
                  "the Ti curve (mean -0.14), the largest misfit being at FMQ+4.89 and FMQ-1 for Ti; "
                  "that is the real size of the agreement. The scaling "
                  "used here is Diffusor's own: to use the data at other temperatures it borrows "
                  "the vacancy and interstitial activation energies of the same element from "
                  "Table 12 of Van Orman & Crispin (2010), for pure magnetite. The entry is "
                  "labelled 'effective' because it is applied beyond the temperature measured; the "
                  "1150 C-only entry is 'chemical', as in the paper. Over FMQ-1 to FMQ+4.89 at "
                  "1150 C (using the Table 12 laws outside their own data windows for Co, Cr and Al) "
                  "this entry and the pure-magnetite Table 12 entry differ by at most 0.50 log units "
                  "for Ti, 0.14 for Mn, 0.22 for Co, 0.43 for Cr (0.13 above FMQ+2) and 0.19 for Al. "
                  "Far from 1150 C the result "
                  "rests on the borrowed energies." if _scaled else
                  "No temperature dependence was measured; evaluating this entry away from "
                  "1150 C raises an error. V3+ and V4+ rows describe the same fitted V data, "
                  "not two independently calibrated transport fields.")),
    ))

# Preserve the old hypothesis keys for reproducibility, but also expose the
# measured fixed-temperature laws for the five species previously T-scaled.
from dataclasses import replace as _replace
for _c in list(COEFFICIENTS):
    if _c.key.endswith("_sievwright2020") and _c.species in TABLE12_PURE:
        _add(_replace(
            _c, key=_c.key + "_1150", label=f"Magnetite {_c.species}, Sievwright (2020), 1150 C ONLY",
            func=_make_sievwright_func(_c.species, False), fixed_temperature_K=SIEVWRIGHT_T_K,
            kind="chemical", secondary_citations=(),
            equation_text=_c.equation_text.split(". Away from")[0],
            verified_from="Primary PDF eq. 5 and Table 5, measured at 1150 C only",
            calibration_notes=("Single-temperature law; no activation energy was measured.",),
            notes=("Original fixed-temperature fit, without Diffusor's legacy temperature-scaling hypothesis. "
                   "Natural magnetite equilibrated with a silicate melt at 1 bar and FMQ-1 to FMQ+4.89. "
                   "Uncertainties on individual log D values are typically below 0.2 log units (1 sigma); "
                   "Supplementary Table S1 lists each (median 0.12, largest 0.51) and the Table 5 fits "
                   "were weighted by them but carry no published uncertainty. "
                   "fO2 convention: eq. 5 is written with the activity a_O2 without a stated standard "
                   "state; Fig. 6 labels its axis 'log fO2 (atm)', but the runs were at 1 bar and dFMQ "
                   "comes from O'Neill (1987) (reference pressure 1 bar), so Diffusor takes fO2 in bar "
                   "(the difference to atm is 0.0057 in log fO2, at most 0.004 in log D)."),
        ))


# --- Fe-Ti interdiffusion ---------------------------------------------------
def _feti_lnD(dc, cond: Conditions, p):
    xTi = np.asarray(cond.X.get("xTi", 0.0), dtype=float)
    return np.exp(p["c0"] + p["a"] * xTi - p["b"] / cond.T_K)


_add(DiffusionCoefficient(
    key="mt_FeTi_freer_hauptman1978",
    kind="interdiffusion",
    mineral="magnetite", species="Fe-Ti",
    label="Titanomagnetite Fe-Ti interdiffusion, Freer & Hauptman (1978)",
    citation="freer_hauptman1978",
    equation_number="Van Orman & Crispin (2010) Table 11",
    equation_text=("ln D = -15.17 + 13.3 x_Ti - 25870/T[K], D in m2/s, 'self-buffered' "
                   "(the combined expression as printed in Van Orman & Crispin 2010, Table 11; "
                   "Freer & Hauptman print D0 and Q per Ti content in their Table I and the "
                   "composition factor beta in their eq. 6)"),
    func=_feti_lnD,
    params={
        "c0": Parameter("c0", -15.17, 0.0, "ln(m2/s)", "1s", "intercept"),
        "a": Parameter("a", 13.3, 0.0, "-", "1s", "x_Ti coefficient"),
        "b": Parameter("b", 25870.0, 0.0, "K", "1s",
                       "Q/R = 215.1 kJ/mol, as printed by Van Orman & Crispin (2.23 eV/k = 25878 K)"),
    },
    sigma_logD=0.5,
    requires=("xTi",),
    needs_fo2=False,
    T_range=Range(1161.15, 1307.15, "K (888-1034 C: the runs that give beta and Fig. 5. The abstract "
                                    "says 600-1034 C, but the 600 and 800 C couples could not be "
                                    "analysed or were erratic. Runs at 1400 C in a controlled "
                                    "atmosphere (log pO2 = -5) are not covered by this self-buffered law)"),
    X_range=Range(0.0, 0.05, "x_Ti (mole fraction of Ti, i.e. Ti per cation site; the tracer-to-5 mol% Ti "
                             "interdiffusion data, TM20 = Fe2.8Ti0.2O4 is 0.067. Larger x_Ti is an "
                             "extrapolation, as in Van Orman & Crispin's Fig. 26 at x_Ti = 0.15)"),
    verified=True,
    verified_from=("Van Orman & Crispin (2010) Table 11 (p. 820), checked on 2 October 2026 against "
                   "Freer & Hauptman (1978): the abstract gives D = 3.85e-3 cm2/s exp(-2.23 eV/kT) at "
                   "3 mol% Ti (2.23 eV/k = 25,879 K, and ln 3.85e-7 m2/s = -15.17 + 13.3 x 0.03) and "
                   "p. 227 gives the composition factor beta = 13.3 +/- 4.8"),
    secondary_citations=("costa2008", "sievwright2020"),
    recommended=False,
    superseded_note=(
        "This is a 1978 calibration under self-buffered conditions, so its oxygen fugacity is "
        "poorly constrained, and it is the oldest entry in the registry. No direct replacement "
        "for Fe-Ti INTERDIFFUSION has been published since Aragon et al. (1984), but for Ti and "
        "Fe TRACER diffusion the Aggarwal & Dieckmann (2002) data tabulated by Van Orman & "
        "Crispin (2010) Table 12 are far better constrained, explicitly fO2-dependent and "
        "temperature-dependent, and are what Tomiya et al. (2013) used at Shinmoedake. Prefer "
        "the Table 12 entries unless you specifically need an interdiffusion coefficient. "
        "Sievwright et al. (2020) add modern magnetite diffusivities for Ti and many other "
        "elements as a function of fO2, but only at 1150 C. Diffusor lists them as the "
        "mt_*_sievwright2020 entries."),
    notes=("Interdiffusion between synthetic Fe3O4 and Fe2.8Ti0.2O4 under self-buffered "
           "conditions (sealed silica tubes), so the fO2 is only loosely constrained (about three "
           "orders of magnitude wide, according to Freer & Hauptman). Costa et al. (2008, Fig. 8) "
           "take their Fe-Ti oxide diffusion data from this paper. The composition variable of "
           "the source is the mole fraction of Ti; the data cover 0 to 5 mol% Ti. Compare with "
           "the Aragon et al. (1984) entry: the two expressions of Van Orman & Crispin (2010) "
           "Table 11 differ by about 1 log unit (this one faster) at 900-1200 C and x_Ti = 0.15, "
           "which is the order-of-magnitude agreement those authors describe."),
))

_add(DiffusionCoefficient(
    key="mt_FeTi_aragon1984",
    kind="interdiffusion",
    mineral="magnetite", species="Fe-Ti",
    label="Titanomagnetite Fe-Ti interdiffusion, Aragon et al. (1984)",
    citation="aragon1984",
    equation_number="Van Orman & Crispin (2010) Table 11",
    equation_text=("ln D = -22.71 + 15.09 x_Ti - 19630/T[K], D in m2/s, 'QFM buffer' "
                   "(as printed in Van Orman & Crispin 2010, Table 11)"),
    func=_feti_lnD,
    params={
        "c0": Parameter("c0", -22.71, 0.0, "ln(m2/s)", "1s", "intercept"),
        "a": Parameter("a", 15.09, 0.0, "-", "1s", "x_Ti coefficient"),
        "b": Parameter("b", 19630.0, 0.0, "K", "1s", "Q/R = 163.2 kJ/mol"),
    },
    sigma_logD=0.5,
    requires=("xTi",),
    needs_fo2=False,
    T_range=Range(1263.15, 1493.15, "K (990-1220 C, the anneals of Aragon et al. 1984, Table 1. "
                                    "Van Orman & Crispin print no temperature range for this row)"),
    X_range=Range(0.0, 0.2, "x_Ti (Aragon et al. couples span X_Usp 0-0.77, i.e. x_Ti up to 0.26 "
                            "with x_Ti = X_Usp/3. Van Orman & Crispin print no range for this row "
                            "and extrapolate to x_Ti = 0.15)"),
    verified=False,
    verified_from=("transcribed from Van Orman & Crispin (2010) RiMG 72, Table 11 (p. 820), checked "
                   "against the rendered table on 2 October 2026. The expression was not found in "
                   "Aragon et al. (1984) itself, so it is not checked against the primary paper"),
    secondary_citations=("sievwright2020",),
    notes=("Aragon et al. controlled fO2 with solid-state buffers (Ni-NiO, quartz-magnetite-fayalite, "
           "wustite-magnetite) and CO2/CO gas mixtures (their Table 1: 990-1220 C, log fO2 -12.6 to "
           "-7.4), so the redox state is better controlled than in the self-buffered experiments of "
           "Freer & Hauptman (1978). They report log D = c0 + c1 x for each couple and calculate no "
           "activation energy; the single Arrhenius expression with the label 'QFM buffer' is the "
           "one printed by Van Orman & Crispin (2010) in Table 11 and was not found in Aragon et al. "
           "Van Orman & Crispin (p. 796) and Aragon et al. (p. 184) both describe the Aragon and Freer & "
           "Hauptman results as being in order of magnitude agreement; at x_Ti = 0.15 the two "
           "Table 11 expressions differ by 0.05 log units at 600 C, 0.9 at 900 C and 1.3 at 1200 C "
           "(Freer & Hauptman faster)."),
))


# --- tracer Arrhenius entries along fixed buffers ------------------------------
def _buffer_arrhenius(dc, cond: Conditions, p):
    return arrhenius(p["D0"], p["Q"] * 1.0e3, cond.T_K)


# Temperature range (K, as printed) of each fit: Table 10 (pure magnetite) and Table 11 (x_Ti = 0.2).
for _key, _sp, _xti, _buf, _D0, _Q, _Trange in [
        ("mt_Fe_aggarwal2002_WM", "Fe", 0.0, "WM", 1.38e-5, 197.0, (1173.0, 1673.0)),
        ("mt_Fe_aggarwal2002_MH", "Fe", 0.0, "MH", 9.10e-6, 175.0, (1173.0, 1673.0)),
        ("mt_Ti_aggarwal2002_WM", "Ti", 0.0, "WM", 2.77e-5, 267.0, (1373.0, 1573.0)),
        ("mt_Ti_aggarwal2002_MH", "Ti", 0.0, "MH", 3.29e-5, 208.0, (1373.0, 1573.0)),
        ("mt_Ti_aggarwal2002_WM_xti02", "Ti", 0.2, "WM", 1.33e-2, 332.0, (1373.0, 1573.0)),
        ("mt_Ti_aggarwal2002_MH_xti02", "Ti", 0.2, "MH", 8.55e-4, 184.0, (1373.0, 1573.0)),
        ("mt_Fe_aggarwal2002_WM_xti02", "Fe", 0.2, "WM", 3.38e-7, 165.0, (1373.0, 1573.0)),
        ("mt_Fe_aggarwal2002_MH_xti02", "Fe", 0.2, "MH", 2.33e-4, 147.0, (1373.0, 1573.0))]:
    _add(DiffusionCoefficient(
        key=_key,
        kind="tracer",
        mineral="magnetite", species=_sp,
        label=(f"{'Titanomagnetite' if _xti else 'Magnetite'} {_sp} tracer along the {_buf} "
               f"buffer, Aggarwal & Dieckmann (2002b)"),
        citation="vanorman_crispin2010",
        equation_number=f"Table {11 if _xti else 10}",
        equation_text=f"D = {_D0:g} exp(-{_Q:g} kJ/mol / R T) m2/s, along the {_buf} buffer"
                      + (f", x_Ti = {_xti:g}" if _xti else ""),
        func=_buffer_arrhenius,
        params={
            "D0": Parameter("D0", _D0, 0.0, "m2/s", "1s", "pre-exponential factor"),
            "Q": Parameter("Q", _Q, 0.0, "kJ/mol", "1s", "activation energy"),
        },
        sigma_logD=0.3,
        needs_fo2=False,
        T_range=Range(_Trange[0], _Trange[1],
                      f"K (Table {11 if _xti else 10} of Van Orman & Crispin: "
                      f"{_Trange[0]:.0f}-{_Trange[1]:.0f} K, as printed)"),
        verified=True,
        verified_from="read from Van Orman & Crispin (2010) RiMG 72, Tables 10 and 11",
        notes=(f"Arrhenius parameters computed by Van Orman & Crispin (2010) *along the {_buf} "
               "buffer* from the fO2-dependent expressions of the source, using the buffer "
               "equations of Huebner (1971); the table footnote says they come from the data "
               "compilation in Table 5 of Aggarwal & Dieckmann (2002b). The paper says only that "
               "they were calculated along the buffer; Diffusor's inference is that they hold only "
               "on that buffer (a pure Arrhenius law cannot carry the fO2 dependence), so for any "
               "other fO2 use the Table 12 entry instead. Check: building D along the buffer from "
               "the partial coefficients of that Table 5 (Frost 1991 buffer here, Huebner 1971 "
               "there) at the temperatures of the data gives this law within 0.14 log units for "
               "x_Ti = 0.2 and within 0.27 for pure magnetite (the latter from only the 1200-1400 C "
               "values of Aggarwal & Dieckmann's own measurements). "
               + ("The range printed for pure magnetite Ti (1373-1573 K) is kept, but Table 5 "
                  "lists Ti for pure magnetite only at 1200, 1300 and 1400 C (1473-1673 K)."
                  if (_sp == "Ti" and not _xti) else
                  "The range printed for pure magnetite Fe (1173-1673 K) is kept, but the tracer "
                  "values of Table 5 of Aggarwal & Dieckmann start at 1100 C (1373 K)."
                  if (_sp == "Fe" and not _xti) else
                  "The range matches the 1100-1300 C measurements of Table 5.")),
    ))
