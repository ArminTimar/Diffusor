"""Primary-source scalar laws audited against the local September 2026 library.

Independent axis fits are projected only in an orthonormal principal frame.
No temperature dependence, covariance, or missing principal diffusivity is invented.
Each entry records the source location of its equation and its calibration limits.
"""
import numpy as np

from ..constants import R_GAS
from .base import DiffusionCoefficient, Parameter, Range, LN10

COEFFICIENTS = []
NO_COVARIANCE = (
    "No joint parameter covariance is published. Quoted marginal errors do not define a joint "
    "uncertainty distribution, so the coefficient is held fixed in the Monte Carlo by default."
)
# Diffusor's own cautions (none of the cited papers makes these statements)
GROWTH = ("Diffusor caution, not from the source: a fitted profile does not establish a diffusion origin; "
          "growth zoning must be evaluated independently.")
PRISTINE = "Crystalline, undamaged lattice only; radiation damage and fluid-assisted replacement are not modelled."


def _axis_law(axis):
    def evaluate(dc, cond, p):
        logD = p[f"logD0_{axis}"] - p[f"Q_{axis}"] * 1000 / (LN10 * R_GAS * cond.T_K)
        if "n_fo2" in p:
            logD += p["n_fo2"] * (cond.log_fo2_bar - p["log_fo2_ref_bar"])
        return 10.0 ** logD
    return evaluate


def _scalar(dc, cond, p):
    return 10.0 ** p["logD0"] * np.exp(-p["Q"] * 1000 / (R_GAS * cond.T_K))


# Pressure label of the calibration range. The laws are pressure independent in Diffusor (no activation
# volume), so the range only decides when the interface warns about extrapolation to high pressure.
P_STATED = "Pa (atmospheric; 1 atm stated by the authors; 1 bar to 1 atm)"
P_ASSUMED = "Pa (atmospheric; not stated by the authors, assumed ambient; 1 bar to 1 atm)"
STATE_DEFAULT = ("ambient pressure, no pressure correction; water content not stated by the authors "
                 "(anhydrous lattice assumed)")


def scalar(key, mineral, species, citation, logD0, Q, temps, *,
           source, notes=(), reference_axis="", allowed_axes=(), kind="chemical",
           uncertainty="", label="", reference_state=STATE_DEFAULT, p_stated=False,
           kind_stated=False, kind_note="", _registry=None):
    """One Arrhenius law ``D = 10**logD0 exp(-Q/RT)`` read from a primary paper.

    ``p_stated`` says the authors state the anneal pressure (otherwise ambient pressure is Diffusor's
    assumption); ``kind_stated`` says the authors call the coefficient by the name of ``kind``
    (otherwise the label is Diffusor's classification, with ``kind_note`` giving the reason). Neither
    changes any computed D: ``kind`` only adds the tracer advice to the warnings for ``tracer`` laws and
    the pressure range only triggers the out-of-range warning.
    """
    assumed = []
    if not kind_stated:
        assumed.append(kind_note or f"The transport kind ('{kind}') is Diffusor's classification; the "
                                    "authors do not label this coefficient as chemical, tracer or interdiffusion.")
    if not p_stated:
        assumed.append("The anneal pressure is not stated by the authors; ambient pressure is assumed for "
                       "the pressure range.")
    c = DiffusionCoefficient(
        key=key, mineral=mineral, species=species, citation=citation,
        label=label or f"{mineral.capitalize()} {species}, {citation}",
        equation_text=f"D = 10^({logD0:.12g}) exp(-{Q:g} kJ/mol / RT) m2/s",
        equation_number=source, func=_scalar,
        params={"logD0": Parameter("logD0", logD0, unit="log10(m2/s)"),
                "Q": Parameter("Q", Q, unit="kJ/mol")},
        T_range=Range(temps[0] + 273.15, temps[1] + 273.15, "K"),
        P_range=Range(1.0e5, 101325.0, P_STATED if p_stated else P_ASSUMED),
        verified=True, verified_from=f"Primary PDF: {citation}, {source}",
        kind=kind, transported_variable=f"{species} concentration",
        reference_state=reference_state, calibration_notes=notes, notes=" ".join(assumed),
        reference_axis=reference_axis, allowed_axes=allowed_axes,
        uncertainty_note=uncertainty or "Coefficient uncertainty is not quantified for sampling in this entry.",
    )
    (COEFFICIENTS if _registry is None else _registry).append(c)
    return c


def _principal(key, species, citation, fits, temps, notes, *, redox=False):
    funcs = {axis: _axis_law(axis) for axis in fits}
    params = {}
    for axis, (a, q) in fits.items():
        params[f"logD0_{axis}"] = Parameter(f"logD0_{axis}", a, unit="log10(m2/s)")
        params[f"Q_{axis}"] = Parameter(f"Q_{axis}", q, unit="kJ/mol")
    if redox:
        params.update(n_fo2=Parameter("n_fo2", .31),
                      log_fo2_ref_bar=Parameter("log_fo2_ref_bar", -12, unit="log10 bar"))
    c = DiffusionCoefficient(
        key=key, mineral="olivine", species=species, citation=citation,
        label=f"Olivine {species}, {citation}, independent a/b/c laws",
        equation_text="D_axis = 10^logD0_axis exp(-Q_axis/RT)" +
                      (" (fO2[bar]/1e-12)^0.31" if redox else "") + "; D = sum(D_axis cos²(angle_axis))",
        equation_number="abstract and Arrhenius fits", func=funcs["c"],
        principal_funcs=funcs, reference_axis="c", orientation_required=True,
        params=params, needs_fo2=redox, fo2_unit="bar",
        T_range=Range(temps[0]+273.15, temps[1]+273.15, "K"),
        P_range=Range(1.0e5, 101325.0, "Pa (atmospheric; assumed by Diffusor, no anneal pressure is stated in "
                                       "the audited text; 1 bar to 1 atm)"),
        kind="tracer" if redox else "chemical", transported_variable=f"dilute {species} concentration",
        verified=True, verified_from=f"Primary PDF {citation}, abstract; rendered equations inspected",
        reference_state="fO2 = 1e-12 bar; Fo83-92" if redox else "Fo100; natural Fo90 comparison",
        calibration_notes=notes, uncertainty_note=NO_COVARIANCE,
    )
    COEFFICIENTS.append(c)
    return c


ca = _principal("ol_Ca_coogan2005", "Ca", "coogan2005ca",
                {"a": (-10.78, 193), "b": (-10.46, 201), "c": (-10.02, 207)},
                (900, 1500), ("Calibrated for Fo83-Fo92; do not apply to pure forsterite.",), redox=True)
# Coogan p. 3691 explicitly defines these errors as standard errors.
for axis, sa, sq in (("a", .43, 11), ("b", .37, 10), ("c", .29, 8)):
    ca.params[f"logD0_{axis}"].sigma = sa
    ca.params[f"Q_{axis}"].sigma = sq

be = _principal("ol_Be_jollands2016", "Be", "jollands2016be",
                {"a": (-4.20, 326.1), "b": (-4.64, 285.8), "c": (-5.82, 227.6)},
                (950, 1475), ("Main fits are for synthetic forsterite; natural Fo90 comparison at 1160-1350 C. "
                             "One oxidising natural-olivine experiment was faster.",))
be.uncertainty_note = ("Quoted errors (logD0, Q kJ/mol): a (0.27, 7.9), b (0.38, 11.6), c (0.15, 4.1). "
                       "Covariance and error convention not yet verified; coefficient uncertainty is not sampled.")

scalar("ol_P_watson2015", "olivine", "P", "watson2015p", -10.06, 229, (650, 850),
       source="abstract p. 2053", p_stated=True,
       notes=(GROWTH, "San Carlos Fo90; evacuated silica ampoules at near-atmospheric pressure, generally "
              "Ni-NiO buffered (two wuestite-magnetite runs); the authors do not call the olivine runs dry."),
       kind_note=("The transport kind ('chemical') is Diffusor's classification: P diffuses in from a "
                  "powder source containing a non-native element; the authors do not name the kind "
                  "(they use 'chemical diffusion' only for P in basaltic melt)."),
       uncertainty="Quoted logD0 ±0.80 and Q ±16 kJ/mol (confidence level not stated); joint uncertainty "
                   "not transcribed and not sampled.")

ni = scalar("ol_Ni_petry2004", "olivine", "Ni", "petry2004", np.log10(3.84e-9), 216,
            (1000, 1400), source="p. 4184, Fig. 6 fixed-fO2 fit", kind="effective",
            reference_axis="c", reference_state="Fo90, fO2=1e-6 Pa, 1 atm",
            p_stated=True,
            kind_note="The transport kind ('effective') is Diffusor's classification; Petry et al. do not name it.",
            notes=("Restricted Fo90 fit at log10 fO2 = -6 Pa (-11 bar). The paper's global 220 kJ/mol slope "
                   "must not be combined with the 216 kJ/mol fit's intercept. No composition or redox extrapolation is implemented.",
                   "Sixfold anisotropy was measured at 1200 C only; fixed axis ratios away from that T are an approximation."))
ni.axis_factors = {"a": 1/6, "b": 1/6, "c": 1.0}
ni.requires = ("XFe",)
ni.X_range = Range(.10, .10, "XFe (Fo90)")
ni.fo2_unit = "Pa"
ni.fo2_range = Range(-6, -6, "log10 Pa")

scalar("qz_Ti_cherniak2007", "quartz", "Ti", "cherniak2007quartz", np.log10(7e-8), 273,
       (700, 1150), source="abstract p. 65 (D0 rounded to 7e-8; the results text gives 7.01e-8, log D0 = -7.154); "
                           "fitted parallel to c (Table 1, Fig. 3)",
       reference_axis="c", allowed_axes=("c",), p_stated=True,
       reference_state="dry, 1 atm (stated by the authors); no pressure correction",
       notes=(GROWTH,
           "Dry 1-atm calibration, fitted for diffusion parallel to c (the abstract writes 'parallel to (001)'; "
           "Table 1 and Fig. 3 say parallel to c). The abstract calls diffusion normal to c slightly slower and the "
           "results text says it does not differ significantly; no transverse law is inferred.",
           "Diffusor caution, not from the source: alternative Ti-in-quartz calibrations differ substantially; "
           "compare literature before interpreting a time."),
       uncertainty="Quoted Q ±12 kJ/mol (confidence level not stated); log D0 = -7.154 ±0.525 in the results "
                   "text; covariance not transcribed, coefficient uncertainty not sampled.")
scalar("rt_Zr_cherniak2007_c", "rutile", "Zr", "cherniak2007rutile", np.log10(9.8e-15), 170,
       (750, 1100), source="abstract p. 267", reference_axis="c", allowed_axes=("c",), p_stated=True,
       kind_stated=True, reference_state="anhydrous (stated by the authors); air, NNO or QFM; 1-atm furnaces",
       notes=("Anhydrous rutile; experiments in air, NNO and QFM. All Zr runs in Table 1 are on synthetic rutile "
              "(the natural-rutile runs are Hf only). Zr transverse law not measured here. Table 1 heads the Zr "
              "rows 'normal to (100)', but the abstract states parallel to c and the run identifiers repeat the Hf "
              "rows normal to (001); Diffusor follows the abstract.",),
       uncertainty="Quoted Q ±30 kJ/mol; joint uncertainty not transcribed and not sampled.")
# The abstract gives 800-1000 C for Hf parallel to c and 750-1050 C normal to c, but Table 1 lists the
# 'normal to (001)' (= parallel to c) block over 750-1050 C and the 'normal to (100)' block over 800-1000 C,
# the opposite assignment. Both laws therefore use 800-1000 C, the part common to both readings.
for ax, D0, Q in (("c", 9.1e-15, 169), ("a", 2.5e-12, 227)):
    scalar(f"rt_Hf_cherniak2007_{ax}", "rutile", "Hf", "cherniak2007rutile", np.log10(D0), Q,
           (800, 1000), source="abstract p. 267", reference_axis=ax, p_stated=True, kind_stated=True,
           allowed_axes=("c",) if ax == "c" else ("a", "b"),
           reference_state="anhydrous (stated by the authors); air, NNO or QFM; 1-atm furnaces",
           notes=("Direction-specific anhydrous rutile law; no fixed anisotropy multiplier applied.",
                  "Temperature range: the abstract gives 800-1000 C for Hf parallel to c and 750-1050 C normal to c, "
                  "while Table 1 lists the run block normal to (001) (parallel to c) over 750-1050 C and the block "
                  "normal to (100) over 800-1000 C. Diffusor uses 800-1000 C, the part common to both readings."))
scalar("ttn_Sr_cherniak1995", "titanite", "Sr", "cherniak1995titanite", -3.57, 415,
       (925, 1175), source="abstract p. 219 (visually checked)", p_stated=True, kind_stated=True,
       reference_state="anhydrous (stated by the authors); air, QFM or self-buffered; 1-atm furnaces",
       notes=("Anhydrous natural titanite. The law is the fit for transport parallel to (100) (abstract); a separate "
              "fit for transport normal to (100) is given in the paper but not implemented, and diffusion normal to "
              "(100) is reported as slightly slower (less than a factor of 2 on average). No resolved Sr redox "
              "effect in the experiments.",
              "The printed equation divides by 2.303 (a rounded ln 10); Diffusor uses the exact ln 10, which "
              "lowers D by about 0.7 % relative to the printed constant."),
       uncertainty="Quoted logD0 ±0.59 and Q ±27 kJ/mol (confidence level not stated); joint uncertainty "
                   "not transcribed and not sampled.")
scalar("ttn_Zr_cherniak2006_c", "titanite", "Zr", "cherniak2006titanite", np.log10(5.33e-7), 325,
       (753, 1100), source="abstract p. 639", reference_axis="c", allowed_axes=("c",), kind_stated=True,
       reference_state="anhydrous, NNO-buffered (stated by the authors); sealed under vacuum; no pressure correction",
       notes=("Anhydrous natural titanite buffered at NNO; QFM checks were similar. Only c-direction fitted. "
              "Capsules were sealed under vacuum in silica glass; the internal pressure is not given.",),
       uncertainty="Quoted Q ±30 kJ/mol and log D0 = -6.27 ±1.28 (confidence level not stated); joint "
                   "uncertainty not transcribed and not sampled.")
