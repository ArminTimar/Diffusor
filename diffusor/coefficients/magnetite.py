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
eruption, and Diffusor reproduces their published values (see the test suite).
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
    "Al": (3.24e-13, -64.9, 6.92e1, 681.0),
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


def _table12_D(species: str, T_K: float, a_O2, table) -> float:
    DV0, QV, DI0, QI = table[species]
    vac = DV0 * np.exp(-QV * 1.0e3 / (R_GAS * T_K)) * a_O2 ** (2.0 / 3.0)
    inter = DI0 * np.exp(-QI * 1.0e3 / (R_GAS * T_K)) * a_O2 ** (-2.0 / 3.0)
    return vac + inter


def _make_table12_func(species: str):
    def _f(dc, cond: Conditions, p):
        """Vacancy + interstitial sum, log-interpolated in x_Ti."""
        a_O2 = 10.0 ** (cond.log_fo2_bar - np.log10(ATM_IN_BAR))   # fO2 in atm
        D_pure = _table12_D(species, cond.T_K, a_O2, TABLE12_PURE)
        if species not in TABLE12_XTI02:
            return D_pure
        xTi = np.asarray(cond.X.get("xTi", 0.0), dtype=float)
        D_02 = _table12_D(species, cond.T_K, a_O2, TABLE12_XTI02)
        # log-linear (geometric) interpolation between x_Ti = 0 and x_Ti = 0.2,
        # which is how Tomiya et al. (2013) obtained their Shinmoedake values.
        f = np.clip(xTi / 0.2, 0.0, 1.5)
        return 10.0 ** (np.log10(D_pure) + f * (np.log10(D_02) - np.log10(D_pure)))
    return _f


for _sp in ("Ti", "Fe", "Mn", "Co", "Cr", "Al"):
    _has02 = _sp in TABLE12_XTI02
    _add(DiffusionCoefficient(
        key=f"mt_{_sp}_vanorman_crispin2010",
        mineral="magnetite", species=_sp,
        label=f"Magnetite {_sp} tracer diffusion, Van Orman & Crispin (2010) Table 12",
        citation="vanorman_crispin2010",
        equation_number="Table 12",
        equation_text=("D* = D_V,0 exp(-Q_V/RT) a_O2^(2/3) + D_I,0 exp(-Q_I/RT) a_O2^(-2/3), "
                       "a_O2 = fO2 / fO2_0 with fO2_0 = 1 atm; "
                       f"{_sp} (x_Ti = 0): D_V,0 = {TABLE12_PURE[_sp][0]:g} m2/s, "
                       f"Q_V = {TABLE12_PURE[_sp][1]:g} kJ/mol, "
                       f"D_I,0 = {TABLE12_PURE[_sp][2]:g} m2/s, Q_I = {TABLE12_PURE[_sp][3]:g} kJ/mol"
                       + (f"; {_sp} (x_Ti = 0.2): D_V,0 = {TABLE12_XTI02[_sp][0]:g} m2/s, "
                          f"Q_V = {TABLE12_XTI02[_sp][1]:g} kJ/mol, "
                          f"D_I,0 = {TABLE12_XTI02[_sp][2]:g} m2/s, "
                          f"Q_I = {TABLE12_XTI02[_sp][3]:g} kJ/mol" if _has02 else "")),
        func=_make_table12_func(_sp),
        params={},
        sigma_logD=0.3,
        requires=("xTi",) if _has02 else (),
        needs_fo2=True, fo2_unit="atm",
        T_range=Range(1373.15, 1673.15, "K (1100-1400 C; extrapolated to magmatic T)"),
        X_range=Range(0.0, 0.2, "x_Ti (Ti per cation site; x_Ti = X_Usp / 3)"),
        verified=True,
        verified_from=("read from the PDF of Van Orman & Crispin (2010) RiMG 72, Table 12 "
                       "(p. 821) and its footnote giving the vacancy/interstitial sum"),
        secondary_citations=("tomiya2013",),
        recommended=(_sp in ("Ti", "Fe")),
        notes=("Underlying tracer data from Aggarwal & Dieckmann (2002b) (Fe, Co, Mn, Ti) and "
               "Dieckmann et al. (1987) (Cr, Al). Values for intermediate x_Ti are obtained "
               "by log-linear interpolation between the x_Ti = 0 and x_Ti = 0.2 entries, "
               "which reproduces the numbers published by Tomiya et al. (2013) for "
               "Shinmoedake: at 950 C, log fO2 = -11 and X_Usp = 0.3 (x_Ti = 0.1) this gives "
               "4.35e-16 m2/s for Ti against their published 4.3e-16 (and 6.84e-16 vs 6.9e-16 "
               "at 900 C), so the Ti/Al branch is reproduced to about 1%. The Fe/Mg branch "
               "agrees to 11% at 900 C (7.3e-15 vs 6.6e-15) but is a factor of 2.1 higher at "
               "950 C (9.4e-15 vs 4.4e-15): Tomiya et al. place the Fe minimum at about 950 C "
               "whereas this implementation puts it near 900-920 C, so they evidently treated "
               "the Fe composition dependence differently. Use the Fe entry with that caveat "
               "and prefer Ti for Shinmoedake timescales until the difference is resolved. "
               "Tomiya et al. took D_Al = D_Ti and D_Mg = D_Fe. Because D passes through a "
               "minimum near 950-980 C, the sensitivity of the retrieved timescale to "
               "temperature is unusually small here (less than a factor of 2)."),
    ))


# --- Fe-Ti interdiffusion ---------------------------------------------------
def _feti_lnD(dc, cond: Conditions, p):
    xTi = np.asarray(cond.X.get("xTi", 0.0), dtype=float)
    return np.exp(p["c0"] + p["a"] * xTi - p["b"] / cond.T_K)


_add(DiffusionCoefficient(
    key="mt_FeTi_freer_hauptman1978",
    mineral="magnetite", species="Fe-Ti",
    label="Titanomagnetite Fe-Ti interdiffusion, Freer & Hauptman (1978)",
    citation="freer_hauptman1978",
    equation_number="Van Orman & Crispin (2010) Table 11",
    equation_text="ln D = -15.17 + 13.3 x_Ti - 25870/T[K], D in m2/s, 'self-buffered'",
    func=_feti_lnD,
    params={
        "c0": Parameter("c0", -15.17, 0.0, "ln(m2/s)", "1s", "intercept"),
        "a": Parameter("a", 13.3, 0.0, "-", "1s", "x_Ti coefficient"),
        "b": Parameter("b", 25870.0, 0.0, "K", "1s", "Q/R = 215.1 kJ/mol"),
    },
    sigma_logD=0.5,
    requires=("xTi",),
    needs_fo2=False,
    T_range=Range(873.15, 1307.15, "K (600-1034 C)"),
    X_range=Range(0.0, 0.2, "x_Ti"),
    verified=True,
    verified_from=("transcribed from Van Orman & Crispin (2010) RiMG 72, Table 11 (p. 820); "
                   "not read from the 1978 primary paper"),
    secondary_citations=("saunders2012", "costa2008"),
    recommended=True,
    notes=("Interdiffusion between synthetic Fe3O4 and Fe2.8Ti0.2O4 under self-buffered "
           "conditions (sealed silica tubes), so the fO2 is only loosely constrained. This is "
           "the coefficient used for Fe-Ti oxide timescales by Costa et al. (2008, Fig. 8) "
           "and Saunders et al. (2012). Compare with the Aragon et al. (1984) entry, which "
           "differs by about an order of magnitude at the same x_Ti."),
))

_add(DiffusionCoefficient(
    key="mt_FeTi_aragon1984",
    mineral="magnetite", species="Fe-Ti",
    label="Titanomagnetite Fe-Ti interdiffusion, Aragon et al. (1984)",
    citation="aragon1984",
    equation_number="Van Orman & Crispin (2010) Table 11",
    equation_text="ln D = -22.71 + 15.09 x_Ti - 19630/T[K], D in m2/s, at the QFM buffer",
    func=_feti_lnD,
    params={
        "c0": Parameter("c0", -22.71, 0.0, "ln(m2/s)", "1s", "intercept"),
        "a": Parameter("a", 15.09, 0.0, "-", "1s", "x_Ti coefficient"),
        "b": Parameter("b", 19630.0, 0.0, "K", "1s", "Q/R = 163.2 kJ/mol"),
    },
    sigma_logD=0.5,
    requires=("xTi",),
    needs_fo2=False,
    T_range=Range(873.15, 1473.15, "K"),
    X_range=Range(0.0, 0.2, "x_Ti"),
    verified=True,
    verified_from=("transcribed from Van Orman & Crispin (2010) RiMG 72, Table 11 (p. 820); "
                   "not read from the 1984 primary paper"),
    notes=("Calibrated at the QFM buffer with solid-state buffering, so the redox state is "
           "better controlled than in Freer & Hauptman (1978). Van Orman & Crispin (2010) "
           "note the two data sets differ by about an order of magnitude when extrapolated "
           "to x_Ti = 0.15."),
))


# --- tracer Arrhenius entries along fixed buffers ------------------------------
def _buffer_arrhenius(dc, cond: Conditions, p):
    return arrhenius(p["D0"], p["Q"] * 1.0e3, cond.T_K)


for _key, _sp, _xti, _buf, _D0, _Q in [
        ("mt_Fe_aggarwal2002_WM", "Fe", 0.0, "WM", 1.38e-5, 197.0),
        ("mt_Fe_aggarwal2002_MH", "Fe", 0.0, "MH", 9.10e-6, 175.0),
        ("mt_Ti_aggarwal2002_WM", "Ti", 0.0, "WM", 2.77e-5, 267.0),
        ("mt_Ti_aggarwal2002_MH", "Ti", 0.0, "MH", 3.29e-5, 208.0),
        ("mt_Ti_aggarwal2002_WM_xti02", "Ti", 0.2, "WM", 1.33e-2, 332.0),
        ("mt_Ti_aggarwal2002_MH_xti02", "Ti", 0.2, "MH", 8.55e-4, 184.0),
        ("mt_Fe_aggarwal2002_WM_xti02", "Fe", 0.2, "WM", 3.38e-7, 165.0),
        ("mt_Fe_aggarwal2002_MH_xti02", "Fe", 0.2, "MH", 2.33e-4, 147.0)]:
    _add(DiffusionCoefficient(
        key=_key,
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
        T_range=Range(1373.15, 1673.15, "K"),
        verified=True,
        verified_from="read from Van Orman & Crispin (2010) RiMG 72, Tables 10 and 11",
        notes=(f"Arrhenius parameters computed by Van Orman & Crispin (2010) *along the {_buf} "
               "buffer* from the fO2-dependent expressions of the source, using the buffer "
               "equations of Huebner (1971). They are only valid on that buffer -- for any "
               "other fO2 use the Table 12 entry instead."),
    ))
