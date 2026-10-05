"""Garnet: tracer families for Fe-Mg-Mn-Ca multicomponent diffusion, and the
Fe-Mg binary computed from the Borinski et al. (2012) tracer laws.

Families (multicomponent solver only)
-------------------------------------
``grt_carlson2006``
    Carlson (2006) eq. 3 with Table 4: ln D* = ln D*0,alm + k (a0 - a0,alm)
    - (Q + P dV)/(R T) + 1/6 ln(fO2 / fO2(graphite)), P in GPa and dV in
    cm3/mol so that P dV is in kJ/mol. a0 is the mole-fraction weighted mean
    of the end-member cell edges (p. 7).
``grt_chakraborty_ganguly1992``
    Chakraborty & Ganguly (1992) eq. 3 with the Arrhenius parameters on
    p. 79, ln D* = ln D0 - Q/(R T) - (P - 1) dV/(R T), D0 in cm2/s, Q in
    cal/mol, P in bar, dV in cm3/mol. D*Ca = 0.5 D*Fe (p. 81). fO2 is the
    graphite-capsule value implicit in Q.
``..._chen_chu2024_mn``
    The same families with the Mn law replaced by the recalibration of Chen
    & Chu (2024), Table 8, from a natural eclogite at about 510 C. The other
    parameters of each model are unchanged, as in that study.

The fO2 term of Carlson (2006) is relative to the graphite-oxygen
equilibrium. Diffusor takes the offset log10(fO2/fO2_graphite) as an option
(default 0, the condition the data were normalised to) instead of computing
the graphite buffer, which needs a CO2 equation of state.
"""
from __future__ import annotations

import numpy as np

from ..constants import R_GAS
from .base import Conditions, DiffusionCoefficient, Parameter, Range, LN10
from .families import TracerFamily, register
from .transport import interdiffusion_from_tracers

COEFFICIENTS = []
COMPONENTS = ("Fe", "Mg", "Mn", "Ca")
END_MEMBERS = {"Fe": "almandine", "Mg": "pyrope", "Mn": "spessartine", "Ca": "grossular"}
_ABBREVIATION = {"Fe": "alm", "Mg": "prp", "Mn": "sps", "Ca": "grs"}      # as in Carlson (2006)
LABELS = {c: f"X_{_ABBREVIATION[c]} = {c}/(Fe+Mg+Mn+Ca)" for c in COMPONENTS}

# Carlson (2006) p. 7: end-member unit-cell edges in nm (Ganguly et al. 1993; Geiger & Feenstra 1997)
A0_NM = {"Fe": 1.1525, "Mg": 1.1456, "Mn": 1.1614, "Ca": 1.1852}
A0_CALIBRATED_MAX = 1.1821       # p. 10: calibration restricted to a0 < 1.1821 nm

# Carlson (2006) Table 4: ln D*0,alm (D* in m2/s), k (1/nm), Q (kJ/mol), dV (cm3/mol)
CARLSON_TABLE4 = {
    "Fe": (-17.638, 465.7, 264.55, 13.568),
    "Mg": (-19.947, 419.9, 244.21, 8.565),
    "Mn": (-17.619, 496.9, 264.44, 9.631),
    "Ca": (-22.572, 511.1, 230.56, 9.795),
}
# Chakraborty & Ganguly (1992) p. 79: D0 (cm2/s), Q (cal/mol), dV (cm3/mol)
CG1992 = {"Fe": (6.4e-4, 65824.0, 5.6), "Mg": (1.1e-3, 67997.0, 5.3), "Mn": (5.1e-4, 60569.0, 6.0)}
CAL = 4.184
# Chen & Chu (2024) Table 8, "this study" columns: log10 D0 (m2/s), Q (kJ/mol)
CHEN_CHU_MN = {"carlson": (-7.75, 257.66), "cg": (-8.40, 225.27)}


def unit_cell_edge(X) -> np.ndarray:
    """Mole-fraction weighted unit-cell edge a0 in nm (Carlson 2006, p. 7)."""
    return sum(np.asarray(X[c], dtype=float) * A0_NM[c] for c in COMPONENTS)


def _carlson(lnD0, k, Q, dV):
    def law(T_K, P_Pa, X, opts):
        P_GPa = P_Pa / 1.0e9
        lnD = (lnD0 + k * (unit_cell_edge(X) - A0_NM["Fe"])
               - (Q + P_GPa * dV) * 1.0e3 / (R_GAS * T_K)
               + opts["dlogfo2_graphite"] * LN10 / 6.0)
        return np.exp(lnD)
    return law


def _cg(D0_cm2, Q_cal, dV):
    def law(T_K, P_Pa, X, opts):
        P_bar = P_Pa / 1.0e5
        return D0_cm2 * 1.0e-4 * np.exp(-(Q_cal * CAL + (P_bar - 1.0) * dV * 0.1) / (R_GAS * T_K))
    return law


def _half(law):
    def ca(T_K, P_Pa, X, opts):
        return 0.5 * law(T_K, P_Pa, X, opts)
    return ca


def _a0_inside(X) -> bool:
    return bool(np.all(unit_cell_edge(X) < A0_CALIBRATED_MAX))


_FO2_HELP = ("log10(fO2 / fO2 at the graphite-oxygen equilibrium). 0 is the condition Carlson's "
             "data were normalised to; D* changes by 1/6 log unit per log unit of fO2.")
_A0_CHECK = (("The garnet is more grossular-rich than the calibration: Carlson (2006, p. 10) restricts "
              "the model to a0 < 1.1821 nm.", _a0_inside),)

_carlson_laws = {c: _carlson(*CARLSON_TABLE4[c]) for c in COMPONENTS}
_carlson_text = ("ln D* = ln D*0,alm + k (a0 - a0,alm) - (Q + P dV)/(R T) + 1/6 ln(fO2/fO2_graphite), "
                 "a0 = sum X_i a0_i; " + "; ".join(
                     f"{c}: ln D*0,alm = {v[0]:g}, k = {v[1]:g} /nm, Q = {v[2]:g} kJ/mol, dV = {v[3]:g} cm3/mol"
                     for c, v in CARLSON_TABLE4.items()))
_CARLSON_COMMON = dict(
    mineral="garnet", components=COMPONENTS, citation="carlson2006",
    T_range=Range(853.15, 1753.15, "K (580-1480 C: natural profiles and experiments)"),
    P_range=Range(1.0e5, 8.5e9, "Pa (up to 8.5 GPa)"),
    options={"dlogfo2_graphite": 0.0}, option_help={"dlogfo2_graphite": _FO2_HELP},
    composition_labels=LABELS, default_dependent="Ca", composition_checks=_A0_CHECK,
    secondary_citations=("chakraborty_ganguly1992", "lasaga1979"),
)

register(TracerFamily(
    key="grt_carlson2006", label="Garnet Fe, Mg, Mn, Ca tracer laws, Carlson (2006)",
    equation_text=_carlson_text, equation_number="eq. 3 and Table 4",
    laws=_carlson_laws, recommended=True,
    verified=True, verified_from="eq. 3 read from the rendered p. 7 and Table 4 (p. 8)",
    sigma_logD=0.8 / 1.96,
    uncertainty_note=("Carlson (2006, p. 10) states that the model specifies each D* within +/-0.8 log "
                      "units at 95 % confidence; Diffusor samples log D* with sigma = 0.41, independently "
                      "for each component (their covariance is not published)."),
    calibration_notes=("The tracer coefficients were retrieved with the ideal multicomponent model, so they "
                       "are used with the same model here (Chakraborty & Ganguly 1992, p. 79).",),
    notes=("Combines high-temperature experiments with stranded diffusion profiles in natural garnet at "
           "580-839 C. The composition dependence enters through the unit-cell edge a0; diffusion is faster "
           "in Mn- and especially Ca-rich garnet."),
    **_CARLSON_COMMON))

# The Mn law of Chen & Chu (2024) was calibrated on a garnet that cooled from about 515 to 480 C
# (thermobarometry in the abstract; garnet mantle at 500-510 C and 1.9-2.0 GPa, p. 13), below the
# 580 C lower bound of Carlson's data. The family has one temperature range, so it is the hull of
# the Mn calibration and the data of the other components.


def _T_CHEN_CHU_HULL(label):
    return Range(753.15, 1753.15, label)


_carlson_chen = dict(_carlson_laws)
_carlson_chen["Mn"] = _carlson(CHEN_CHU_MN["carlson"][0] * LN10, CARLSON_TABLE4["Mn"][1],
                               CHEN_CHU_MN["carlson"][1], CARLSON_TABLE4["Mn"][3])
register(TracerFamily(
    key="grt_carlson2006_chen_chu2024_mn",
    label="Garnet tracer laws, Carlson (2006) with the Mn law of Chen & Chu (2024)",
    equation_text=(_carlson_text + ". Mn replaced by Chen & Chu (2024) Table 8: log10 D*0,alm = -7.75 +/- 0.66, "
                   "Q = 257.66 +/- 19.00 kJ/mol, k and dV unchanged"),
    equation_number="Carlson (2006) eq. 3 and Table 4; Chen & Chu (2024) Table 8",
    laws=_carlson_chen,
    verified=True, verified_from="Table 8 of Chen & Chu (2024); Carlson (2006) as for grt_carlson2006",
    uncertainty_note=("Chen & Chu (2024) give 1 SD on log D0 and Q of the Mn law without their correlation; "
                      "the coefficients are held fixed in the Monte Carlo."),
    calibration_notes=("The Mn law is recalibrated from one natural eclogite garnet at about 510 C and 2 GPa "
                       "with an assumed 5 +/- 3 Myr duration; the other components keep Carlson's laws.",
                       "The temperature range is the hull of the Mn calibration (about 480-515 C) and "
                       "Carlson's data (580-1480 C) for the other components. Fe, Mg and Ca are not "
                       "calibrated between 480 and 580 C."),
    **{**_CARLSON_COMMON,
       "T_range": _T_CHEN_CHU_HULL("K (480-1480 C: Mn law 480-515 C, other components 580-1480 C)")}))

_cg_laws = {c: _cg(*CG1992[c]) for c in ("Fe", "Mg", "Mn")}
_cg_laws["Ca"] = _half(_cg_laws["Fe"])
_cg_text = ("ln D* = ln D0 - Q/(R T) - (P - 1 bar) dV/(R T), D*Ca = 0.5 D*Fe; " + "; ".join(
    f"{c}: D0 = {v[0]:g} cm2/s, Q = {v[1]:g} cal/mol, dV = {v[2]:g} cm3/mol" for c, v in CG1992.items()))
_CG_COMMON = dict(
    mineral="garnet", components=COMPONENTS, citation="chakraborty_ganguly1992",
    T_range=Range(1373.15, 1753.15, "K (1100-1480 C)"),
    P_range=Range(1.4e9, 4.3e9, "Pa (14-43 kbar)"),
    composition_labels=LABELS, default_dependent="Ca", secondary_citations=("lasaga1979",),
)
_CG_NOTES = ("fO2 is that of the graphite capsule and is built into Q; no fO2 correction is applied "
             "(p. 79).",
             "Calibrated on almandine-spessartine couples; D*Ca = 0.5 D*Fe is an assumption based on "
             "Loomis et al. (1985) (p. 81).")
register(TracerFamily(
    key="grt_chakraborty_ganguly1992",
    label="Garnet Fe, Mg, Mn (Ca = 0.5 Fe) tracer laws, Chakraborty & Ganguly (1992)",
    equation_text=_cg_text, equation_number="eq. 3 and p. 79; D*Ca on p. 81",
    laws=_cg_laws,
    verified=True, verified_from=("Arrhenius parameters on p. 79 and eq. 2 read from the rendered PDF; "
                                  "the eq. 5 matrix is reproduced in the test suite"),
    uncertainty_note=("Q and dV carry +/-1 sigma from the data scatter (p. 79) but no covariance with D0; the "
                      "coefficients are held fixed in the Monte Carlo."),
    calibration_notes=_CG_NOTES, **_CG_COMMON))

_cg_chen = dict(_cg_laws)
_mn_log, _mn_Q = CHEN_CHU_MN["cg"]


def _cg_chen_mn(T_K, P_Pa, X, opts):
    return 10.0 ** _mn_log * np.exp(-(_mn_Q * 1.0e3 + (P_Pa / 1.0e5 - 1.0) * CG1992["Mn"][2] * 0.1)
                                    / (R_GAS * T_K))


_cg_chen["Mn"] = _cg_chen_mn
register(TracerFamily(
    key="grt_chakraborty_ganguly1992_chen_chu2024_mn",
    label="Garnet tracer laws, Chakraborty & Ganguly (1992) with the Mn law of Chen & Chu (2024)",
    equation_text=_cg_text + ". Mn replaced by Chen & Chu (2024) Table 8: log10 D0 = -8.40 +/- 0.67, "
                             "Q = 225.27 +/- 18.99 kJ/mol, dV unchanged",
    equation_number="Chakraborty & Ganguly (1992) p. 79; Chen & Chu (2024) Table 8",
    laws=_cg_chen,
    verified=True, verified_from="Table 8 of Chen & Chu (2024)",
    uncertainty_note="The coefficients are held fixed in the Monte Carlo.",
    calibration_notes=_CG_NOTES + ("The Mn law is recalibrated from one natural eclogite garnet at about "
                                   "510 C and 2 GPa (Chen & Chu 2024). Table 8 gives an activation volume "
                                   "only for the original models, so the dV = 6.0 cm3/mol of Chakraborty & "
                                   "Ganguly (1992) is kept for the new Mn law (Diffusor's assumption).",
                                   "The temperature range is the hull of the Mn calibration (about "
                                   "480-515 C) and the Chakraborty & Ganguly data (1100-1480 C) for Fe, "
                                   "Mg and Ca."),
    **{**_CG_COMMON,
       "T_range": _T_CHEN_CHU_HULL("K (480-1480 C: Mn law 480-515 C, other components 1100-1480 C)")}))


# ---------------------------------------------------------------------------
# Borinski et al. (2012): Fe and Mg tracer laws for almandine-pyrope garnet (eq. 11, p. 582)
# log D* = log D0 - (Q_1bar + dV (P - 1)) / (2.303 R T), P in bar, dV in J/bar/mol
# ---------------------------------------------------------------------------
def _borinski(dc, cond: Conditions, p):
    P_bar = cond.P_Pa / 1.0e5
    return p["D0"] * np.exp(-(p["Q"] * 1.0e3 + p["dV"] * (P_bar - 1.0)) / (R_GAS * cond.T_K))


for _el, _D0, _sD0, _Q, _sQ, _dV in (("Mg", 2.72e-10, 4.52e-10, 228.3, 20.3, 0.53),
                                     ("Fe", 1.64e-10, 2.54e-10, 226.9, 18.6, 0.56)):
    COEFFICIENTS.append(DiffusionCoefficient(
        key=f"grt_{_el}_borinski2012", kind="tracer",
        mineral="garnet", species=_el,
        label=f"Garnet {_el} tracer, Borinski et al. (2012)",
        citation="borinski2012", equation_number="eq. 11 and p. 582",
        equation_text=(f"log D*_{_el} = log D0 - (Q + dV (P - 1)) / (2.303 R T), D0 = ({_D0:g} +/- {_sD0:g}) m2/s, "
                       f"Q = {_Q:g} +/- {_sQ:g} kJ/mol, dV = {_dV:g} J/bar/mol (from Chakraborty & Ganguly 1992), P in bar"),
        func=_borinski,
        params={"D0": Parameter("D0", _D0, 0.0, "m2/s", "unstated",
                                f"pre-exponential factor (quoted +/- {_sD0:g}, level not stated, not sampled)"),
                "Q": Parameter("Q", _Q, 0.0, "kJ/mol", "unstated",
                               f"activation energy at 1 bar (quoted +/- {_sQ:g}, level not stated, not sampled)"),
                "dV": Parameter("dV", _dV, 0.0, "J/bar/mol", "1s", "activation volume")},
        T_range=Range(1330.15, 1705.15, "K (1057-1432 C, with the refitted Ganguly et al. 1998 runs)"),
        P_range=Range(2.0e9, 4.0e9, "Pa (20-40 kbar)"),
        transported_variable=f"{_el} tracer in almandine-pyrope garnet",
        reference_state="almandine-pyrope garnet; fO2 of graphite in the C-O-H system (built into dV and Q)",
        verified=True, verified_from="eq. 11 and the parameters on p. 582 of the PDF",
        uncertainty_note=("The quoted errors of D0 exceed D0 and have no stated covariance with Q; the "
                          "coefficient is held fixed in the Monte Carlo."),
        calibration_notes=("Tracer coefficients retrieved from Fe-Mg diffusion couples with a multicomponent "
                           "model. For a zoning profile use the Fe-Mg interdiffusion law computed from them.",),
        secondary_citations=("chakraborty_ganguly1992",),
    ))

COEFFICIENTS.append(interdiffusion_from_tracers(
    "grt_FeMg_borinski2012", tracer_a=COEFFICIENTS[1], tracer_b=COEFFICIENTS[0],
    species="Fe-Mg", comp_key="XFe",
    label="Garnet Fe-Mg interdiffusion from the Borinski et al. (2012) tracer laws",
    citation="borinski2012",
    equation_number=("eq. 8, p. 578 (interdiffusion in an ideal system of equally charged species, the binary "
                     "form of eq. 1) with the eq. 11 tracer laws"),
    X_range=Range(0.0, 1.0, "XFe = Fe/(Fe+Mg) (the whole binary join, not a calibration range)"),
    calibration_notes=(
        "Binary Fe-Mg exchange. Borinski et al. note that most metamorphic garnets are nearly binary "
        "almandine-pyrope solid solutions poor in spessartine and grossular (p. 573). Diffusor's own "
        "restriction: use it only where Mn and Ca are low and nearly uniform, since they are not "
        "modelled. Otherwise use the multicomponent garnet model (File > Multicomponent and isotope "
        "study).",
        "Ideal mixing; Borinski et al. (2012) found that non-ideality changes retrieved coefficients by "
        "less than a factor of 1.2 for most natural garnet compositions (abstract).",),
    notes="Computed from the two tracer laws for an ideal Fe-Mg binary of divalent ions.",
))
