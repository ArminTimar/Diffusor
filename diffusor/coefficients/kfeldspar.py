"""Diffusion coefficients for K-feldspar (sanidine, orthoclase).

The laws used in diffusion chronometry of silicic magmas are Sr and Ba in
sanidine (Chamberlain et al. 2014, Audetat, Grocolas & Mutch 2026, section 2.3).
Both were measured on the same Or61 sanidine by Rutherford backscattering, so
the pair can be compared directly. For Sr, Cherniak (1996) reports that the
(010) anneals gave results similar to the (001) anneals. For Ba, the abstract of Cherniak (2002) says diffusion normal to
(010) appears slightly slower than normal to (001), while its text calls the two
comparable (a single (010) run on sanidine). The minerals are treated as
isotropic.

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
        "Q": Parameter("Q", 450.0, 13.0, "kJ/mol", "unstated", "activation energy"),
    },
    sigma_logD=0.03,
    requires=(),
    needs_fo2=False,
    T_range=Range(998.15, 1348.15, "K (725-1075 C)"),
    P_range=Range(1.0e5, 101325.0, "Pa (atmospheric, dry; 1 bar to 1 atm)"),
    X_range=Range(0.61, 0.61, "X_Or (sanidine Or61)"),
    verified=True,
    verified_from=("the Arrhenius relation in the abstract. Chamberlain et al. (2014) Table 1 "
                   "lists the same D0 of 8.4 m2/s and Q of 450 kJ/mol"),
    secondary_citations=("chamberlain2014", "audetat2026"),
    recommended=True,
    notes=("Cherniak (1996) reports that anneals buffered at QFM (FMQ) and the (010) anneals gave "
           "results similar to the (001) anneals in air; the law is the fit to the (001) anneals in "
           "air and has no fO2 term. The paper prints log10 D0 = 0.9252 +/- 0.5893 for it (level "
           "not stated, no covariance with Q); D0 is held fixed here. "
           "Audetat, Grocolas & Mutch (2026, p. 13) summarise the paper as "
           "showing a small fO2 effect ascribed to defects from Fe oxidation. The scatter of "
           "0.03 log units is the uncertainty that Chamberlain et al. (2014, p. 14) state for D0 "
           "and E of this law; they derive it in their Electronic Appendix 7, which Diffusor does "
           "not have, and give no confidence level. It is much smaller than the +/-13 kJ/mol on Q "
           "alone would suggest (about 0.6 log units at 800 C, ignoring its correlation with D0). "
           "Audetat, Grocolas & Mutch (2026, p. 14) report that Cherniak's (1996) Or61 Sr "
           "diffusivity at 825 C is 0.6 log units higher than that of Toussaint et al. (2025) for "
           "Or98 orthoclase. Those results are a conference abstract without an Arrhenius law, so "
           "Diffusor does not include them yet."),
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
        "Q": Parameter("Q", 455.0, 20.0, "kJ/mol", "unstated", "activation energy"),
    },
    sigma_logD=0.12,
    requires=(),
    needs_fo2=False,
    T_range=Range(1101.15, 1348.15, "K (828-1075 C, the sanidine runs; the 775-1124 C of the abstract covers all three feldspars)"),
    P_range=Range(1.0e5, 101325.0, "Pa (atmospheric, dry; 1 bar to 1 atm)"),
    X_range=Range(0.61, 0.61, "X_Or (sanidine Or61)"),
    verified=True,
    verified_from=("the abstract, which gives 455 +/- 20 kJ/mol (108.5 +/- 4.8 kcal/mol) and "
                   "2.9e-1 m2/s. Chamberlain et al. (2014) Table 1 lists the same D0"),
    secondary_citations=("chamberlain2014", "audetat2026"),
    recommended=True,
    notes=("Measured on the same Or61 sanidine as the Sr law of Cherniak (1996). Diffusion "
           "normal to (010) appears to be slightly slower than normal to (001) (abstract); the "
           "text calls the two comparable, and only one (010) run was made on sanidine (976 C). "
           "The sanidine runs span 828-1075 C (Table 1); the 775-1124 C of the abstract covers "
           "all three feldspars. Ba is about 1.7 log units slower than Sr at 800 C (computed from "
           "the two laws), in line with the 'about 1.7 orders of magnitude' that Audetat, "
           "Grocolas & Mutch (2026, p. 15) give. The scatter of 0.12 log units is the uncertainty "
           "that Chamberlain et al. (2014, p. 14) state for D0 and E of this law; they derive it "
           "in their Electronic Appendix 7, which Diffusor does not have, and give no confidence "
           "level. It is much smaller than the +/-20 kJ/mol on Q alone would suggest (about 1 log "
           "unit at 800 C, ignoring its correlation with D0). Ba and Sr partition strongly into "
           "sanidine, so partial dissolution or fast growth can leave Ba-Sr steps that look like "
           "diffusion (Audetat, Grocolas & Mutch 2026, section 3.3)."),
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
        "D0": Parameter("D0", 3.01e-6, 0.0, "m2/s", "1s", "pre-exponential factor "
                        "(the paper prints log10 D0 = -5.52 +/- 2.06; held fixed, see notes)"),
        "Q": Parameter("Q", 342.0, 47.0, "kJ/mol", "unstated", "activation energy"),
    },
    sigma_logD=None,
    requires=(),
    needs_fo2=False,
    T_range=Range(1073.15, 1273.15, "K (800-1000 C)"),
    P_range=Range(1.0e5, 101325.0, "Pa (atmospheric; 1 bar to 1 atm)"),
    verified=True,
    verified_from="the Arrhenius relation in the open-access abstract",
    recommended=True,
    notes=("Experiments buffered at NNO gave diffusivities similar to those in air, hydrous species had little effect "
           "and Ti shows little anisotropy. The paper prints both uncertainties (log10 D0 = "
           "-5.52 +/- 2.06 and Q = 342 +/- 47 kJ/mol, Fig. 3 caption) but no covariance of the "
           "two and no summary scatter of log D about the fit, so the coefficient is held fixed "
           "in the Monte Carlo by default. Varying Q alone with D0 fixed would misrepresent the "
           "spread of D, and varying both independently would ignore their correlation."),
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
        P_range=Range(1.0e5, 101325.0, "Pa (close to ambient pressure; 1 bar to 1 atm)"),
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
             "is not described by one constant D (abstract and p. 1309)."),
            "Directly measured exchange coefficient. Values calculated from Na and K tracer "
            "coefficients (eq. 3) are 5-10 times faster normal to (001) and almost 100 times "
            "faster normal to (010) (p. 1313).",
            "The Arrhenius parameters printed in the abstract and on p. 1310 were extracted from "
            "the 800, 950 and 1000 C data only. A log-linear fit to Table 3 at those three "
            "temperatures does not reproduce them (it gives Q = 201, 208, 293 and 282 kJ/mol for "
            "the four laws, against the printed 179, 182, 272 and 269). This is an inconsistency "
            "inside the paper, so the printed parameters are used as published.",
            "The 850 C data (measured D higher than the law by factors of about 1.5 to 2.2) and "
            "the anomalous 920 C data (D higher than at 950 and 1000 C; normal to (001) the law "
            "is low by a factor of about 3.5 to 3.8) were not used for the fit and are not "
            "reproduced by it (p. 1310, Fig. 6, Table 3).",
            "The abstract says 'close to ambient pressure'; each plate was sealed with the salt "
            "in an evacuated silica glass tube and heated in a box furnace at atmospheric "
            "pressure (p. 1304). No numerical pressure is printed, so the numerical pressure "
            "bounds (1 bar to 1 atm) are Diffusor's rendering of that wording.",
        ),
        notes=("Diffusion normal to (001) is about one order of magnitude faster and less "
               "temperature dependent than normal to (010) (abstract)."),
    ))
