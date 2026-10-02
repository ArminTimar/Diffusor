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
GROWTH = "A fitted profile does not establish a diffusion origin; growth zoning must be evaluated independently."
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


def scalar(key, mineral, species, citation, logD0, Q, temps, *,
           source, notes=(), reference_axis="", allowed_axes=(), kind="chemical",
           uncertainty="", label="", reference_state="dry lattice; no pressure correction", _registry=None):
    c = DiffusionCoefficient(
        key=key, mineral=mineral, species=species, citation=citation,
        label=label or f"{mineral.capitalize()} {species}, {citation}",
        equation_text=f"D = 10^({logD0:.12g}) exp(-{Q:g} kJ/mol / RT) m2/s",
        equation_number=source, func=_scalar,
        params={"logD0": Parameter("logD0", logD0, unit="log10(m2/s)"),
                "Q": Parameter("Q", Q, unit="kJ/mol")},
        T_range=Range(temps[0] + 273.15, temps[1] + 273.15, "K"),
        P_range=Range(1e5, 101325, "Pa (near atmospheric)"),
        verified=True, verified_from=f"Primary PDF: {citation}, {source}",
        kind=kind, transported_variable=f"{species} concentration",
        reference_state=reference_state, calibration_notes=notes,
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
        P_range=Range(1e5, 101325, "Pa (main calibration at atmospheric pressure)"),
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
       source="abstract p. 2053", notes=(GROWTH, "San Carlos Fo90; NNO-buffered dry experiments, with WM checks."),
       uncertainty="Quoted logD0 ±0.80 and Q ±16 kJ/mol; joint uncertainty not transcribed and not sampled.")

ni = scalar("ol_Ni_petry2004", "olivine", "Ni", "petry2004", np.log10(3.84e-9), 216,
            (1000, 1400), source="p. 4184, Fig. 6 fixed-fO2 fit", kind="effective",
            reference_axis="c", reference_state="Fo90, fO2=1e-6 Pa, 1 atm",
            notes=("Restricted Fo90 fit at log10 fO2 = -6 Pa (-11 bar). The paper's global 220 kJ/mol slope "
                   "must not be combined with the 216 kJ/mol fit's intercept. No composition or redox extrapolation is implemented.",
                   "Sixfold anisotropy was measured at 1200 C only; fixed axis ratios away from that T are an approximation."))
ni.axis_factors = {"a": 1/6, "b": 1/6, "c": 1.0}
ni.requires = ("XFe",)
ni.X_range = Range(.10, .10, "XFe (Fo90)")
ni.fo2_unit = "Pa"
ni.fo2_range = Range(-6, -6, "log10 Pa")

scalar("qz_Ti_cherniak2007", "quartz", "Ti", "cherniak2007quartz", np.log10(7e-8), 273,
       (700, 1150), source="abstract p. 65; fitted c-direction",
       reference_axis="c", allowed_axes=("c",), notes=(GROWTH,
           "Dry 1-atm calibration. Slightly slower transverse diffusion was observed; no transverse law inferred.",
           "Alternative Ti-in-quartz calibrations differ substantially; compare literature before interpreting a time."),
       uncertainty="Quoted Q ±12 kJ/mol; covariance not transcribed, coefficient uncertainty not sampled.")
scalar("rt_Zr_cherniak2007_c", "rutile", "Zr", "cherniak2007rutile", np.log10(9.8e-15), 170,
       (750, 1100), source="abstract p. 267", reference_axis="c", allowed_axes=("c",),
       notes=("Dry synthetic and natural rutile; experiments in air, NNO and QFM. Zr transverse law not measured here.",),
       uncertainty="Quoted Q ±30 kJ/mol; joint uncertainty not transcribed and not sampled.")
for ax, D0, Q, temps in (("c", 9.1e-15, 169, (800, 1000)), ("a", 2.5e-12, 227, (750, 1050))):
    scalar(f"rt_Hf_cherniak2007_{ax}", "rutile", "Hf", "cherniak2007rutile", np.log10(D0), Q,
           temps, source="abstract p. 267", reference_axis=ax,
           allowed_axes=("c",) if ax == "c" else ("a", "b"),
           notes=("Direction-specific dry rutile law; no fixed anisotropy multiplier applied.",))
scalar("ttn_Sr_cherniak1995", "titanite", "Sr", "cherniak1995titanite", -3.57, 415,
       (925, 1175), source="abstract p. 219 (visually checked)",
       notes=("Dry natural titanite; no resolved Sr orientation or redox effect in the experiments.",),
       uncertainty="Quoted logD0 ±0.59 and Q ±27 kJ/mol; joint uncertainty not transcribed and not sampled.")
scalar("ttn_Zr_cherniak2006_c", "titanite", "Zr", "cherniak2006titanite", np.log10(5.33e-7), 325,
       (753, 1100), source="abstract p. 639", reference_axis="c", allowed_axes=("c",),
       notes=("Dry natural titanite buffered at NNO; QFM checks were similar. Only c-direction fitted.",),
       uncertainty="Quoted Q ±30 kJ/mol; joint uncertainty not transcribed and not sampled.")
