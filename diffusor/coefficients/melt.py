"""Diffusion of water in silicate melt.

Behrens, Zhang & Xu (2004) measured bulk water diffusion (H2Ot, the sum of H2O
molecules and OH groups) in dacitic and andesitic melts at 0.5-1.5 GPa and
1458-1858 K with diffusion couples of dry and hydrous glass. They treat it as
effective binary interdiffusion between H2O and the rest of the melt (abstract)
and find D_H2Ot proportional to the water content in dacite up to 6 wt% (their
eq. 4, D_H2Ot = C_H2Ot D_H2O(1 wt%)), with a weaker dependence in andesite.

The laws here therefore depend on the diffusing quantity itself: the profile
is the composition variable ``H2Ot`` in wt%, and D is zero where the melt is
dry. Behrens et al. (2004, p. 5147) note that eq. 4 may not apply as C_H2Ot
approaches 0 wt% (ppm levels), where OH diffusion takes over.

Read from the rendered pages 5146-5149 on 13 October 2026.
"""
from __future__ import annotations

import numpy as np

from .base import DiffusionCoefficient, Parameter, Range

COEFFICIENTS = []

_VERIFIED = ("Behrens, Zhang & Xu (2004) GCA 68:5139-5150, eqs 4, 7, 8 and 10 and Table 3 (pp. 5146-5149), "
             "read from the rendered pages on 13 October 2026")
_PROPORTIONAL = ("D is taken proportional to C_H2Ot (eq. 4), so it is zero in dry melt. The paper notes that this "
                 "may not hold as C_H2Ot approaches 0 wt% (ppm levels), where OH diffusion dominates (p. 5147).")


def _add(c):
    COEFFICIENTS.append(c)
    return c


def _with_latex(builder):
    def wrap(func):
        func.latex = builder
        return func
    return wrap


def _water(cond):
    return np.clip(np.asarray(cond.x("H2Ot"), dtype=float), 0.0, None)


@_with_latex(lambda c: (rf"\log_{{10}} D_{{H_2O_t}} = \log_{{10}} C_{{H_2O_t}} {c.params['a'].value:+g} "
                        rf"- \frac{{{abs(c.params['b'].value):g}}}{{T}}"
                        r"\quad (D\ \mathrm{in\ m^2/s},\ C\ \mathrm{in\ wt\%})"))
def _proportional(dc, cond, p):
    """log10 D = log10 C_H2Ot + a + b/T, D in m2/s (Behrens et al. 2004 eqs 7 and 8)."""
    return _water(cond) * 10.0 ** (p["a"] + p["b"] / cond.T_K)


_add(DiffusionCoefficient(
    key="melt_H2O_behrens2004_dacite",
    mineral="melt", species="H2O", kind="effective",
    transported_variable="C_H2Ot, total dissolved water in wt%",
    label="Melt H2O (dacite), Behrens et al. (2004) eq. 7",
    citation="behrens2004",
    equation_number="7, p. 5147",
    equation_text=("log D_H2Ot = log C_H2Ot - 5.93(+/-0.33) - 7266(+/-540)/T, D in m2/s, C_H2Ot in wt%, T in K, "
                   "1 sigma errors; 1458-1798 K and 0-6 wt% H2Ot"),
    func=_proportional,
    params={
        "a": Parameter("a", -5.93, 0.33, "log10(m2/s)", "1s", "intercept"),
        "b": Parameter("b", -7266.0, 540.0, "K", "1s", "1/T coefficient"),
    },
    sigma_logD=0.06,
    requires=("H2Ot",), needs_fo2=False,
    T_range=Range(1458.0, 1798.0, "K (eq. 7)"),
    P_range=Range(0.5e9, 1.5e9, "Pa (0.5-1.5 GPa, the experiments; no pressure effect was resolved)"),
    X_range=Range(0.0, 6.0, "wt% H2Ot"),
    verified=True, verified_from=_VERIFIED,
    notes=("Dacite with 67 wt% SiO2 (Table 3). " + _PROPORTIONAL + " The pressure effect was too small to "
           "resolve between 0.5 and 1.5 GPa (p. 5147). Extrapolated to 600 C the law lies 0.8 log units below "
           "the low-pressure dehydration data of Liu et al. (2004), and the paper advises caution with that "
           "extrapolation (p. 5147). The scatter about the fit is the standard error the paper gives, 0.06 log "
           "units."),
))

_add(DiffusionCoefficient(
    key="melt_H2O_behrens2004_andesite",
    mineral="melt", species="H2O", kind="effective",
    transported_variable="C_H2Ot, total dissolved water in wt%",
    label="Melt H2O (andesite), Behrens et al. (2004) eq. 8",
    citation="behrens2004",
    equation_number="8, p. 5147",
    equation_text=("log D_H2Ot = log C_H2Ot - 5.41(+/-0.54) - 7415(+/-926)/T, D in m2/s, C_H2Ot in wt%, T in K, "
                   "for D at 1 wt% H2Ot, 1573-1813 K"),
    func=_proportional,
    params={
        "a": Parameter("a", -5.41, 0.54, "log10(m2/s)", "1s", "intercept"),
        "b": Parameter("b", -7415.0, 926.0, "K", "1s", "1/T coefficient"),
    },
    sigma_logD=0.08,
    requires=("H2Ot",), needs_fo2=False,
    T_range=Range(1573.0, 1813.0, "K (eq. 8)"),
    P_range=Range(0.5e9, 1.5e9, "Pa (0.5-1.5 GPa, the experiments; no pressure effect was resolved)"),
    X_range=Range(0.0, 5.6, "wt% H2Ot (Table 3; the law is for D at 1 wt%)"),
    verified=True, verified_from=_VERIFIED,
    notes=("Andesite with 57 wt% SiO2 (Table 3). The paper obtained eq. 8 for D_H2Ot at 1 wt% H2Ot and writes "
           "it with log C_H2Ot. For andesite it finds the dependence of D on C_H2Ot weak and says that a rough "
           "proportionality at low water contents (0.1-0.4 wt%) cannot be ruled out (pp. 5147-5148). Diffusor evaluates eq. 8 as "
           "printed, so D is proportional to C_H2Ot at every content. Extrapolation to magmatic temperatures "
           "(1100-1300 K) may lead to a large uncertainty, and the dependence on C_H2Ot cannot be predicted "
           "there (p. 5147). The scatter about the fit is the standard error the paper gives, 0.08 log units."),
))


@_with_latex(lambda c: (r"\log_{10} D_{H_2O_t} = \log_{10} C_{H_2O_t} + (-0.757 - 0.0868\,S) + "
                        r"\frac{-14785 + 131.7\,S}{T} + (3.079 - 0.0490\,S)\,\frac{P}{T}"
                        rf"\quad (S = {c.params['S'].value:g}\ \mathrm{{wt\%\ SiO_2}},\ P\ \mathrm{{in\ MPa}})"))
def _behrens_sio2(dc, cond, p):
    """Behrens et al. (2004) eq. 10 at 1 wt% H2Ot, times C_H2Ot (eq. 4)."""
    S, T, P = p["S"], cond.T_K, cond.P_Pa / 1e6
    log_d1 = ((p["a0"] + p["a1"] * S) + (p["b0"] + p["b1"] * S) / T + (p["c0"] + p["c1"] * S) * P / T)
    return _water(cond) * 10.0 ** log_d1


_add(DiffusionCoefficient(
    key="melt_H2O_behrens2004_sio2",
    mineral="melt", species="H2O", kind="effective",
    transported_variable="C_H2Ot, total dissolved water in wt%",
    label="Melt H2O (rhyolite to andesite, SiO2 term), Behrens et al. (2004) eq. 10 with eq. 4",
    citation="behrens2004",
    equation_number="10 (p. 5148) with 4 (p. 5146)",
    equation_text=("log D_H2Ot(1 wt%) = (-0.757 - 0.0868 S) + (-14,785 + 131.7 S)/T + (3.079 - 0.0490 S) P/T, "
                   "D in m2/s, T in K, P in MPa, S in wt% SiO2; D_H2Ot = C_H2Ot D_H2Ot(1 wt%) (eq. 4)"),
    func=_behrens_sio2,
    params={
        "a0": Parameter("a0", -0.757, 0.0, "log10(m2/s)", "1s", "intercept"),
        "a1": Parameter("a1", -0.0868, 0.0, "1/wt%", "1s", "SiO2 term of the intercept"),
        "b0": Parameter("b0", -14785.0, 0.0, "K", "1s", "1/T coefficient"),
        "b1": Parameter("b1", 131.7, 0.0, "K/wt%", "1s", "SiO2 term of the 1/T coefficient"),
        "c0": Parameter("c0", 3.079, 0.0, "K/MPa", "1s", "P/T coefficient"),
        "c1": Parameter("c1", -0.0490, 0.0, "K/(MPa wt%)", "1s", "SiO2 term of the P/T coefficient"),
        "S": Parameter("S", 67.0, 0.0, "wt%", "1s",
                       "SiO2 content of the melt; 67 (the dacite of Table 3) is Diffusor's default, set the "
                       "value of the melt"),
    },
    sigma_logD=0.12,
    requires=("H2Ot",), needs_fo2=False,
    T_range=Range(803.0, 1848.0, "K (Table 3, all data)"),
    P_range=Range(0.1e6, 1.5e9, "Pa (0.1-1500 MPa, Table 3)"),
    X_range=Range(0.0, 2.0, "wt% H2Ot (proportionality stated for rhyolite to dacite up to 2 wt%)"),
    verified=True, verified_from=_VERIFIED,
    notes=("Fitted to 69 data at 1 wt% H2Ot for rhyolite (77 wt% SiO2), haplogranite (76), dacite (67) and "
           "andesite (57) (Table 3), with a standard error of 0.12 log units (p. 5148). The paper states that "
           "with the proportionality of eq. 4 it can also be used for rhyolite to dacite with up to 2 wt% H2Ot "
           "at magmatic temperatures, that for andesitic melts the relation between D and C_H2Ot is not known "
           "at magmatic temperatures, that no data constrain melts below 65 wt% SiO2 under 1500 K, and that "
           "applying it to other melt compositions is not recommended (trachyte is underestimated by 0.5 log "
           "units) (p. 5149). " + _PROPORTIONAL + " The SiO2 content S is a parameter of the entry."),
))
