"""Fe-Mg and Mg diffusion coefficients for orthopyroxene."""
from __future__ import annotations

import numpy as np

from ..constants import R_GAS
from .base import Conditions, DiffusionCoefficient, Parameter, Range, LN10, arrhenius

COEFFICIENTS = []


def _add(c: DiffusionCoefficient):
    COEFFICIENTS.append(c)
    return c


# ---------------------------------------------------------------------------
# Dohmen, ter Heege, Becker & Chakraborty (2016) Am Mineral 101, 2210-2221
# ---------------------------------------------------------------------------
def _dohmen2016_fs9(dc, cond: Conditions, p):
    """log D_FeMg//[001] = log D0 - Q/(ln10 R T) + n log fO2(Pa)   [eq. 1]

    with the compositional correction recommended on p. 2219:
        D(XFe, T, fO2) = D(T, fO2) * 10^(m (XFe - 0.09)),  m = 1
    """
    logD = (p["logD0"]
            - p["Q"] * 1.0e3 / (LN10 * R_GAS * cond.T_K)
            + p["n"] * cond.log_fo2_Pa)
    XFe = cond.X.get("XFe", None)
    if XFe is not None:
        logD = logD + p["m"] * (np.asarray(XFe, dtype=float) - 0.09)
    return 10.0 ** logD


_add(DiffusionCoefficient(
    key="opx_FeMg_dohmen2016",
    mineral="opx", species="Fe-Mg",
    label="Opx Fe-Mg, Dohmen et al. (2016) -- Fs9, fO2-dependent",
    citation="dohmen2016",
    equation_number="1 (+ compositional correction, p. 2219)",
    equation_text=("log D_Fe-Mg^c = log D0 - {Q / ln(10) R T} + n log fO2   with "
                   "Q = 308 +/- 23 kJ/mol, log(D0 [m2/s]) = -5.95 +/- 0.83, n = 0.053 +/- 0.027; "
                   "abstract form D_Fe-Mg [m2/s] = 1.12e-6 (fO2 [Pa])^0.053 exp[-308 kJ/mol /(R T)]; "
                   "D(XFe,T,fO2) = D(T,fO2) * 10^(m (XFe - 0.09)) with m = 1"),
    func=_dohmen2016_fs9,
    params={
        "logD0": Parameter("logD0", -5.95, 0.83, "log10(m2/s)", "1s", "pre-exponential factor"),
        "Q": Parameter("Q", 308.0, 23.0, "kJ/mol", "1s", "activation energy"),
        "n": Parameter("n", 0.053, 0.027, "-", "1s", "fO2 exponent (fO2 in Pa)"),
        "m": Parameter("m", 1.0, 0.0, "-", "1s", "composition exponent, 10^(m(XFe-0.09))"),
    },
    sigma_logD=0.10,
    axis_factors={"c": 1.0, "b": 1.0, "a": 1.0 / 3.5},
    reference_axis="c",
    requires=("XFe",),
    needs_fo2=True, fo2_unit="Pa",
    T_range=Range(1143.15, 1373.15, "K (870-1100 C)"),
    P_range=Range(1.0e5, 1.0e5, "Pa (1 atm)"),
    fo2_range=Range(-11.0, -7.0, "log10 Pa"),
    X_range=Range(0.09, 0.5, "XFe"),
    verified=True,
    verified_from="read from the paper PDF: abstract, eq. 1 (p. 2215), anisotropy p. 2215, composition p. 2219, run table p. 2213",
    secondary_citations=("sato2022", "polo_sanchez2023"),
    recommended=True,
    notes=("Reference axis is [001]. D//[100] = D//[001]/3.5; D//[010] is indistinguishable "
           "from D//[001]. The compositional term is only recommended for XFe = 0.09-0.5. "
           "Test case: at 950 C and log fO2 = -7 Pa the equation gives log D = -19.48, "
           "against the measured -19.49 +/- 0.07 (run OPXD_14, sample 7D_20)."),
))


def _dohmen2016_fs1(dc, cond: Conditions, p):
    """Opx8 (XFe = 0.01): single Arrhenius, no resolvable fO2 dependence."""
    return arrhenius(10.0 ** p["logD0"], p["Q"] * 1.0e3, cond.T_K)


_add(DiffusionCoefficient(
    key="opx_FeMg_dohmen2016_fs1",
    mineral="opx", species="Fe-Mg",
    label="Opx Fe-Mg, Dohmen et al. (2016) -- Fs1, fO2-independent",
    citation="dohmen2016",
    equation_number="1 with n = 0 (Opx8 fit, p. 2215)",
    equation_text=("D_Fe-Mg [m2/s] = 1.66e-4 exp[-377 +/- 30 kJ/mol /(R T)]; "
                   "log(D0 [m2/s]) = -3.78 +/- 1.26, Q = 377 +/- 30 kJ/mol"),
    func=_dohmen2016_fs1,
    params={
        "logD0": Parameter("logD0", -3.78, 1.26, "log10(m2/s)", "1s", "pre-exponential factor"),
        "Q": Parameter("Q", 377.0, 30.0, "kJ/mol", "1s", "activation energy"),
    },
    sigma_logD=0.10,
    axis_factors={"c": 1.0, "b": 1.0, "a": 1.0 / 3.5},
    reference_axis="c",
    needs_fo2=False,
    T_range=Range(1143.15, 1373.15, "K (870-1100 C)"),
    P_range=Range(1.0e5, 1.0e5, "Pa (1 atm)"),
    X_range=Range(0.0, 0.05, "XFe"),
    verified=True,
    verified_from="read from the paper PDF: abstract and p. 2215",
    notes="For near-end-member enstatite only (XFe ~ 0.01); no fO2 dependence resolved.",
))


# ---------------------------------------------------------------------------
# Ganguly & Tazzoli (1994), in the form used by Allan et al. (2013)
# ---------------------------------------------------------------------------
def _ganguly_tazzoli(dc, cond: Conditions, p):
    """log D [m2/s] = -9.54 + 2.6 XFe - 12530/T + (1/6) log(fO2_sample/fO2_IW)

    Quoted as Equation 3 of Dohmen et al. (2016), who attribute the general
    expression to Allan et al. (2013) derived from Ganguly & Tazzoli (1994).
    The fO2 term is optional: Ostorero et al. (2022) eq. 1 use the same
    expression without it (they write log D [cm2/s] = -5.54 + 2.6 XFe - 12530/T,
    which is the same law, 1 cm2/s = 1e-4 m2/s).
    """
    XFe = np.asarray(cond.x("XFe"), dtype=float)
    logD = p["c0"] + p["a"] * XFe - p["b"] / cond.T_K
    if p.get("use_fo2", 1.0) >= 0.5:
        from ..thermo.buffers import log_fo2_buffer
        d_iw = cond.log_fo2_bar - log_fo2_buffer("IW", cond.T_K, cond.P_Pa)
        logD = logD + p["n_fo2"] * d_iw
    return 10.0 ** logD


_add(DiffusionCoefficient(
    key="opx_FeMg_ganguly_tazzoli1994",
    mineral="opx", species="Fe-Mg",
    label="Opx Fe-Mg, Ganguly & Tazzoli (1994) / Allan et al. (2013) form",
    citation="ganguly_tazzoli1994",
    equation_number="Dohmen et al. (2016) eq. 3; Ostorero et al. (2022) eq. 1",
    equation_text=("log D_Fe-Mg [m2/s] = -9.54 + 2.6 XFe - 12530/T[K] "
                   "+ (1/6) log( fO2(sample,T) / fO2(IW,T) )"),
    func=_ganguly_tazzoli,
    params={
        "c0": Parameter("c0", -9.54, 0.0, "log10(m2/s)", "1s", "intercept"),
        "a": Parameter("a", 2.6, 0.0, "-", "1s", "XFe coefficient"),
        "b": Parameter("b", 12530.0, 0.0, "K", "1s", "Q/(ln10 R) = 239.9 kJ/mol"),
        "n_fo2": Parameter("n_fo2", 1.0 / 6.0, 0.0, "-", "1s", "fO2 exponent relative to IW"),
        "use_fo2": Parameter("use_fo2", 1.0, 0.0, "-", "1s", "1 = apply fO2 term, 0 = omit it"),
    },
    sigma_logD=0.5,
    axis_factors={"c": 1.0, "b": 1.0, "a": 1.0 / 3.5},
    reference_axis="c",
    requires=("XFe",),
    needs_fo2=True, fo2_unit="bar",
    T_range=Range(773.15, 1073.15, "K (500-800 C, extrapolated above)"),
    fo2_range=Range(None, None, "relative to IW"),
    X_range=Range(0.10, 0.50, "XFe"),
    verified=True,
    verified_from=("transcribed from Dohmen et al. (2016) eq. 3 (PDF p. 2216) and cross-checked "
                   "against Ostorero et al. (2022) eq. 1, which is the same law in cm2/s"),
    secondary_citations=("dohmen2016", "ostorero2022"),
    notes=("Retrieved indirectly from Fe-Mg order-disorder kinetics, not from a direct "
           "diffusion measurement. Dohmen et al. (2016) argue the fO2 exponent of 1/6 "
           "(taken by analogy with olivine) is not supported by their experiments, which "
           "give n = 0.053, and that it produces an artificially low activation energy "
           "(150 kJ/mol) at constant fO2. Set use_fo2 = 0 to drop the fO2 term, as "
           "Ostorero et al. (2022) did. Kept for comparison with older studies."),
))


# ---------------------------------------------------------------------------
# Schwandt, Cygan & Westrich (1998) -- Mg self-diffusion in orthoenstatite
# ---------------------------------------------------------------------------
def _schwandt(dc, cond: Conditions, p):
    return arrhenius(10.0 ** p["logD0"], p["Q"] * 1.0e3, cond.T_K)


for _axis, _label, _logD0, _slog, _Q, _sQ in [
        ("a", "(100)", -3.96, 2.48, 360.0, 52.0),
        ("b", "(010)", -5.16, 3.69, 339.0, 77.0),
        ("c", "(001)", -8.36, 3.27, 265.0, 66.0)]:
    _add(DiffusionCoefficient(
        key=f"opx_Mg_schwandt1998_{_axis}",
        mineral="opx", species="Mg",
        label=f"Opx Mg self-diffusion // {_axis}-axis {_label}, Schwandt et al. (1998)",
        citation="schwandt1998",
        equation_number="Table 3",
        equation_text=(f"D = D0 exp(-Ea/RT); {_label}: Ea = {_Q:g} +/- {_sQ:g} kJ/mol, "
                       f"log D0 [m2/s] = {_logD0:g} +/- {_slog:g}"),
        func=_schwandt,
        params={
            "logD0": Parameter("logD0", _logD0, _slog, "log10(m2/s)", "1s", "pre-exponential factor"),
            "Q": Parameter("Q", _Q, _sQ, "kJ/mol", "1s", "activation energy"),
        },
        reference_axis=_axis,
        needs_fo2=False,
        T_range=Range(1023.15, 1173.15, "K (750-900 C)"),
        P_range=Range(1.0e5, 1.0e5, "Pa (1 atm)"),
        X_range=Range(0.10, 0.10, "En90Fs10"),
        verified=True,
        verified_from="read from the paper PDF: abstract and Table 3",
        notes=("25Mg self-diffusion by SIMS depth profiling at the IW buffer on En90Fs10. "
               "The three axes agree within their (large) uncertainties, so the apparent "
               "anisotropy is not resolved. Dohmen et al. (2016) note these tracer data are "
               "consistent with their interdiffusion data."),
    ))
