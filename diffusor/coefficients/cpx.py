"""Fe-Mg and Ca-Mg diffusion coefficients for clinopyroxene."""
from __future__ import annotations

import numpy as np

from .base import Conditions, DiffusionCoefficient, Parameter, Range, arrhenius

COEFFICIENTS = []


def _add(c):
    COEFFICIENTS.append(c)
    return c


def _plain_arrhenius(dc, cond: Conditions, p):
    return arrhenius(p["D0"], p["Q"] * 1.0e3, cond.T_K)


# ---------------------------------------------------------------------------
# Mueller, Dohmen, Becker, ter Heege & Chakraborty (2013) CMP 166, 1563-1576
# ---------------------------------------------------------------------------
_add(DiffusionCoefficient(
    key="cpx_FeMg_muller2013",
    kind="interdiffusion",
    mineral="cpx", species="Fe-Mg",
    label="Cpx Fe-Mg // [001], Mueller et al. (2013)",
    citation="muller2013",
    equation_number="abstract / Fig. 5b",
    equation_text=("D_Fe-Mg = 2.77 (+/- 4.27) x 10^-7 exp(-320.7 +/- 16.0 kJ/mol / R T) m2/s"),
    func=_plain_arrhenius,
    params={
        "D0": Parameter("D0", 2.77e-7, 4.27e-7, "m2/s", "unstated", "pre-exponential factor"),
        "Q": Parameter("Q", 320.7, 16.0, "kJ/mol", "unstated", "activation energy"),
    },
    sigma_logD=0.5,
    reference_axis="c",
    needs_fo2=False,
    T_range=Range(1073.15, 1473.15, "K (800-1200 C)"),
    P_range=Range(1.0e5, 101325.0, "Pa (atmospheric; 1 bar to 1 atm)"),
    fo2_range=Range(-17.0, -11.0, "log10 bar (no dependence resolved)"),
    fo2_unit="bar",
    X_range=Range(0.07, 0.35, "XFe (Di93 to Di65, the compositional range the authors estimate; Di93Hd7 was the crystal measured)"),
    verified=True,
    verified_from="read from the paper PDF: abstract and Fig. 5b annotation",
    recommended=True,
    notes=("Measured along [001] on Di93Hd7 by RBS after pulsed-laser-deposition thin films. "
           "No fO2 dependence resolved between 1e-17 and 1e-11 bar, unlike Dimanov & "
           "Wiedenbeck (2006). The authors relate this to the high Al content of their "
           "crystals. D0 is printed as 2.77e-7 m2/s in the abstract and in the annotation of "
           "Fig. 5b (p. 1569) and as 2.77e-8 m2/s in the running text on p. 1569 (right column), "
           "which Diffusor reads as a typographical slip. It uses 2.77e-7, the value that "
           "appears twice; with Q = 320.7 kJ/mol it gives log D = -19.7 at 1000 C, on the "
           "data of Fig. 5b, whereas 2.77e-8 would give -20.7, a full log unit below them. "
           "The +/- values (D0 +/-4.27e-7, Q +/-16.0 kJ/mol) are printed without a stated "
           "confidence level, and the error on D0 is larger than D0 itself, so it cannot be a "
           "symmetric 1 sigma interval (it would allow negative D0); the parameters are therefore "
           "not sampled independently and the Monte Carlo uses the scatter of log D instead. "
           "The paper says the Arrhenius expression reproduces the measured coefficients "
           "within 1 log unit (envelope of Fig. 5b) and puts the reproducibility of individual "
           "coefficients at 0.4 log units. Reading the 1 log unit envelope as roughly two "
           "standard deviations is Diffusor's own conversion, hence sigma_logD = 0.5 "
           "(treated as 1 sigma). The composition range of the data is Di93 to Di65 (XFe about "
           "0.07 to 0.35); the authors find no composition dependence over it. XFe below 0.07 "
           "(towards pure diopside) is extrapolation."),
))


# ---------------------------------------------------------------------------
# Dimanov & Sautter (2000) -- the coefficient used by NIDIS
# ---------------------------------------------------------------------------
_add(DiffusionCoefficient(
    key="cpx_FeMg_dimanov_sautter2000",
    kind="interdiffusion",
    mineral="cpx", species="Fe-Mg",
    label="Cpx (Fe,Mn)-Mg, Dimanov & Sautter (2000) -- as used by NIDIS",
    citation="dimanov_sautter2000",
    equation_number="least-squares fit, p. 757 (abstract p. 749)",
    equation_text=("D = D0 exp(-dH / R T). dH = 406 +/- 64 kJ/mol and log D0 [cm2/s] = "
                   "-0.02 +/- 0.32 (1 sigma), i.e. D0 = 0.955 cm2/s = 9.55 x 10^-5 m2/s"),
    func=_plain_arrhenius,
    params={
        "D0": Parameter("D0", 9.55e-5, 0.0, "m2/s", "1s", "pre-exponential factor"),
        "Q": Parameter("Q", 406.0, 0.0, "kJ/mol", "1s", "apparent activation enthalpy"),
    },
    sigma_logD=0.5,
    reference_axis="c",
    needs_fo2=False,
    fo2_unit="atm",
    T_range=Range(1173.15, 1513.15, "K (900-1240 C)"),
    fo2_range=Range(-18.0, -13.0, "log10 atm (pO2 10^-18 to 10^-13 atm, abstract; it varied with T, Table 1: 4e-18 to 2.5e-13 atm)"),
    verified=True,
    verified_from=("read from the paper PDF: abstract (p. 749) and the fit on p. 757. Before "
                   "1 October 2026 this entry used 9.5e-5 m2/s from the Table 2 footnote of "
                   "Petrone et al. (2016), 0.5 per cent low. NIDIS createfit.m uses 9.55e-5"),
    secondary_citations=("petrone2016",),
    notes=("This is the coefficient behind the published NIDIS timescales, kept so that "
           "Diffusor can reproduce them. Reproduces the D values quoted by Petrone et al. "
           "(2016) to the digits printed: 3.26e-20 m2/s at 1098 C and 1.20e-19 m2/s at "
           "1150 C. These are (Fe,Mn)-Mg interdiffusion coefficients "
           "on natural diopside at pO2 between 1e-18 and 1e-13 atm (IW to QIF), not fixed, so "
           "406 kJ/mol is an apparent enthalpy with no fO2 correction. Leaving out the "
           "1240 C point gives 332 kJ/mol (p. 757). The published 1 sigma values (64 kJ/mol, "
           "0.32 log units) are strongly correlated and no covariance is given, so they are "
           "not sampled. Diffusion was measured along [001] only (film on surfaces cut "
           "perpendicular to [001]); the entry carries no anisotropy factors, so other "
           "directions are not calibrated. The paper reports pO2 in atm (the fO2 range shown is "
           "in log10 atm); the law itself has no fO2 term and the pO2 varied with temperature in "
           "the experiments. Mueller et al. (2013) is the more recent and better constrained "
           "calibration and gives markedly different timescales."),
))


# ---------------------------------------------------------------------------
# Brady & McCallister (1983) -- Ca-Mg interdiffusion
# ---------------------------------------------------------------------------
_add(DiffusionCoefficient(
    key="cpx_CaMg_brady1983",
    kind="interdiffusion",
    mineral="cpx", species="Ca-Mg",
    label="Cpx Ca-Mg effective binary interdiffusion, Brady & McCallister (1983)",
    citation="brady_mccallister1983",
    equation_number="eq. 5, p. 100 (abstract p. 95), homogenisation experiments at 25 kbar",
    equation_text=("D = 3.89 x 10^-7 exp(-360.87 kJ/mol / R T) m2/s "
                   "(= 3.89 x 10^-3 exp(-86.25 kcal/mol / R T) cm2/s), uncertainty a factor of 2"),
    func=_plain_arrhenius,
    params={
        "D0": Parameter("D0", 3.89e-7, 0.0, "m2/s", "1s", "pre-exponential factor"),
        "Q": Parameter("Q", 360.87, 0.0, "kJ/mol", "1s", "activation energy"),
    },
    sigma_logD=0.30,
    reference_axis="c",
    needs_fo2=False,
    T_range=Range(1423.15, 1523.15, "K (1150-1250 C)"),
    P_range=Range(2.5e9, 2.5e9, "Pa (25 kbar)"),
    verified=True,
    verified_from=("read from the paper PDF: abstract (p. 95) and eq. 5 (p. 100). "
                   "86.25 kcal/mol x 4.184 = 360.87 kJ/mol and 3.89e-3 cm2/s = 3.89e-7 m2/s"),
    notes=("From homogenisation of (001) pigeonite lamellae in sub-calcic diopside. The "
           "diffusion is one-dimensional and normal to the (001) lamellae, i.e. along c* "
           "(shown as the c direction; for monoclinic cpx c* is not exactly [001]), and no "
           "other direction was measured. The "
           "stated uncertainty is a factor of 2, hence sigma_logD = 0.30. The statement that "
           "Ca-Mg is slower than Fe-Mg in cpx is from the abstract of Mueller et al. (2013), "
           "not from Brady & McCallister; using Ca profiles to test whether a boundary formed "
           "by growth or by diffusion is Diffusor's suggestion, not a claim of the source. "
           "These are 'average' "
           "effective binary coefficients: the paper warns that near the solvus actual "
           "interdiffusion coefficients should vary by an order of magnitude or more with "
           "composition."),
))
