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
    kind="chemical",
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
    kind="chemical",
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
           "normal to (010) appears to be slightly slower than normal to (001) (abstract). Ba is about 1.7 log units slower "
           "than Sr at 800 C, which Audetat, Grocolas & Mutch (2026) confirm. The scatter of 0.12 "
           "log units is the uncertainty Chamberlain et al. (2014) derived from D0 and Q. Ba "
           "partitions strongly into sanidine, so partial dissolution or fast growth can leave "
           "Ba steps that look like diffusion (Audetat, Grocolas & Mutch 2026, section 3.3)."),
))


_add(DiffusionCoefficient(
    key="kfs_Ti_cherniak_watson2020",
    kind="chemical",
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
           "the fit, so the coefficient is held fixed in the Monte Carlo by default. Varying Q "
           "alone with D0 fixed would overstate the spread of D."),
))


# Schaffer et al. (2014): Na-K interdiffusion measured directly by cation exchange.
# Abstract (p. 1300): D_NaK = D0 exp(-EA/RT) for two compositions and two directions.
# Their eq. 3 predicts D_NaK from Na and K tracer coefficients; the direct values are
# 5-10 times slower normal to (001) and almost 100 times slower normal to (010) than
# that prediction (p. 1313), which is why these exchange laws are listed separately.
_SCHAFFER = (
    ("001", 0.92, 5.18e-8, 179.0), ("001", 0.98, 1.82e-7, 182.0),
    ("010", 0.92, 5.81e-5, 272.0), ("010", 0.98, 1.34e-4, 269.0),
)
for _plane, _xor, _D0, _EA in _SCHAFFER:
    _add(DiffusionCoefficient(
        key=f"kfs_NaK_schaffer2014_{_plane}_or{round(_xor * 100)}",
        kind="interdiffusion",
        mineral="kfeldspar", species="Na-K",
        label=f"K-feldspar Na-K interdiffusion normal to ({_plane}), X_Or {_xor:g}, Schaffer et al. (2014)",
        citation="schaffer2014",
        equation_number="abstract (p. 1300) and p. 1310",
        equation_text=f"D_NaK = {_D0:g} exp(-{_EA:g} kJ/mol / R T) m2/s, normal to ({_plane}), X_Or = {_xor:g}",
        func=_arrhenius_law,
        params={
            "D0": Parameter("D0", _D0, 0.0, "m2/s", "1s", "pre-exponential factor"),
            "Q": Parameter("Q", _EA, 0.0, "kJ/mol", "1s", "activation energy"),
        },
        requires=(),
        needs_fo2=False,
        # (010) is normal to b in a monoclinic crystal; the normal to (001) is not a
        # crystal axis, so no axis or direction cosine is accepted for it.
        reference_axis="b" if _plane == "010" else "",
        allowed_axes=("b",) if _plane == "010" else ("normal_(001)",),
        T_range=Range(1073.15, 1273.15, "K (800-1000 C)"),
        P_range=Range(1.0e5, 1.0e5, "Pa (close to ambient pressure)"),
        X_range=(Range(0.65, 0.95, "X_Or (D_NaK nearly constant here)") if _xor < 0.95
                 else Range(0.95, 0.99, "X_Or (D_NaK rises steeply with X_Or)")),
        transported_variable="X_Or = K/(K+Na), mole fraction",
        reference_state="disordered Eifel sanidine, cation exchange with alkali-halide melt",
        verified=True,
        verified_from="read from the rendered abstract (p. 1300)",
        uncertainty_note="No uncertainties are given for D0 or EA; the coefficient is held fixed in the Monte Carlo.",
        calibration_notes=(
            ("For 0.65 < X_Or < 0.95 D_NaK is nearly independent of composition (abstract)."
             if _xor < 0.95 else
             "Above X_Or = 0.95 D_NaK rises steeply with X_Or; a profile that spans this range "
             "is not described by one constant D (p. 1310)."),
            "Directly measured exchange coefficient. Values calculated from Na and K tracer "
            "coefficients (eq. 3) are 5-10 times faster normal to (001) and almost 100 times "
            "faster normal to (010) (p. 1313).",
        ),
        notes=("Diffusion normal to (001) is about one order of magnitude faster and less "
               "temperature dependent than normal to (010) (abstract)."),
    ))
