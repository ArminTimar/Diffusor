"""Diffusion coefficients for K-feldspar (sanidine, orthoclase).

The laws used in diffusion chronometry of silicic magmas are Sr and Ba in
sanidine (Chamberlain et al. 2014, Audetat, Grocolas & Mutch 2026, section 2.3).
Both were measured on the same Or61 sanidine by Rutherford backscattering, so
the pair can be compared directly. Neither shows a resolvable dependence on
crystal orientation, and the minerals are treated as isotropic.

The Arrhenius parameters come from the abstracts of the papers. Chamberlain et
al. (2014, Table 1) list the same D0 values for Sr and Ba. Neither paper
publishes a covariance of D0 and Q.
"""
from __future__ import annotations

from ..constants import R_GAS
from .base import Conditions, DiffusionCoefficient, Parameter, Range, arrhenius

COEFFICIENTS = []


def _add(c):
    COEFFICIENTS.append(c)
    return c


def _arrhenius_law(dc, cond: Conditions, p):
    """D = D0 exp(-Q / RT), D0 in m2/s and Q in kJ/mol."""
    return arrhenius(p["D0"], p["Q"] * 1.0e3, cond.T_K)


_add(DiffusionCoefficient(
    key="kfs_Sr_cherniak1996",
    mineral="kfeldspar", species="Sr",
    label="K-feldspar Sr, Cherniak (1996)",
    citation="cherniak1996",
    equation_number="abstract",
    equation_text="D = 8.4 exp(-450 +/- 13 kJ/mol / RT) m2/s, sanidine Or61, normal to (001)",
    func=_arrhenius_law,
    params={
        "D0": Parameter("D0", 8.4, 0.0, "m2/s", "1s", "pre-exponential factor"),
        "Q": Parameter("Q", 450.0, 13.0, "kJ/mol", "1s", "activation energy"),
    },
    sigma_logD=0.03,
    requires=(),
    needs_fo2=False,
    T_range=Range(998.15, 1348.15, "K (725-1075 C)"),
    P_range=Range(1.0e5, 1.0e5, "Pa (1 atm, dry)"),
    X_range=Range(0.61, 0.61, "X_Or (sanidine Or61)"),
    verified=True,
    verified_from=("the Arrhenius relation in the abstract. Chamberlain et al. (2014) Table 1 "
                   "lists the same D0 of 8.4 m2/s and Q of 450 kJ/mol"),
    secondary_citations=("chamberlain2014", "audetat2026"),
    recommended=True,
    notes=("Experiments in air and at the FMQ buffer showed a small fO2 effect, which the law "
           "does not include. The scatter of 0.03 log units is the uncertainty Chamberlain et al. "
           "(2014) derived for this law from D0 and Q. Audetat, Grocolas & Mutch (2026) report "
           "that the Or98 orthoclase experiments of Toussaint et al. (2025) give Sr diffusivities "
           "about 0.6 log units lower at 825 C. Those results are a conference abstract without "
           "an Arrhenius law, so Diffusor does not include them yet."),
))


_add(DiffusionCoefficient(
    key="kfs_Ba_cherniak2002",
    mineral="kfeldspar", species="Ba",
    label="K-feldspar Ba, Cherniak (2002)",
    citation="cherniak2002",
    equation_number="abstract",
    equation_text="D = 2.9e-1 exp(-455 +/- 20 kJ/mol / RT) m2/s, sanidine Or61, normal to (001)",
    func=_arrhenius_law,
    params={
        "D0": Parameter("D0", 0.29, 0.0, "m2/s", "1s", "pre-exponential factor"),
        "Q": Parameter("Q", 455.0, 20.0, "kJ/mol", "1s", "activation energy"),
    },
    sigma_logD=0.12,
    requires=(),
    needs_fo2=False,
    T_range=Range(1048.15, 1397.15, "K (775-1124 C)"),
    P_range=Range(1.0e5, 1.0e5, "Pa (1 atm, dry)"),
    X_range=Range(0.61, 0.61, "X_Or (sanidine Or61)"),
    verified=True,
    verified_from=("the abstract, which gives 455 +/- 20 kJ/mol (108.5 +/- 4.8 kcal/mol) and "
                   "2.9e-1 m2/s. Chamberlain et al. (2014) Table 1 lists the same D0"),
    secondary_citations=("chamberlain2014", "audetat2026"),
    recommended=True,
    notes=("Measured on the same Or61 sanidine as the Sr law of Cherniak (1996). Diffusion "
           "normal to (010) matched diffusion normal to (001). Ba is about 1.7 log units slower "
           "than Sr at 800 C, which Audetat, Grocolas & Mutch (2026) confirm. The scatter of 0.12 "
           "log units is the uncertainty Chamberlain et al. (2014) derived from D0 and Q. Ba "
           "partitions strongly into sanidine, so partial dissolution or fast growth can leave "
           "Ba steps that look like diffusion (Audetat, Grocolas & Mutch 2026, section 3.3)."),
))


_add(DiffusionCoefficient(
    key="kfs_Ti_cherniak_watson2020",
    mineral="kfeldspar", species="Ti",
    label="K-feldspar Ti, Cherniak & Watson (2020)",
    citation="cherniak_watson2020",
    equation_number="abstract",
    equation_text="D_Ksp = 3.01e-6 exp(-342 +/- 47 kJ/mol / RT) m2/s, normal to (001)",
    func=_arrhenius_law,
    params={
        "D0": Parameter("D0", 3.01e-6, 0.0, "m2/s", "1s", "pre-exponential factor"),
        "Q": Parameter("Q", 342.0, 47.0, "kJ/mol", "1s", "activation energy"),
    },
    sigma_logD=None,
    requires=(),
    needs_fo2=False,
    T_range=Range(1073.15, 1273.15, "K (800-1000 C)"),
    P_range=Range(1.0e5, 1.0e5, "Pa (1 atm)"),
    verified=True,
    verified_from="the Arrhenius relation in the open-access abstract",
    recommended=True,
    notes=("Experiments buffered at NNO matched those in air, water had little effect and Ti "
           "shows little anisotropy. The paper gives no uncertainty on D0 and no scatter about "
           "the fit, so the Monte Carlo can only vary Q. With D0 held fixed that overstates the "
           "spread of D, so switch off coefficient sampling or treat the result as an upper "
           "bound on the uncertainty."),
))
