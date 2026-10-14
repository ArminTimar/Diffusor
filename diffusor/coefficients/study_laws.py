"""Diffusion laws in the form published application studies printed and used.

Diffusor's other law modules hold the calibrations of the experimental papers.
A study that applies a calibration often writes it differently: another
composition reference, another anisotropy factor, a missing or added term, its
own regression of a compiled data set, or one element's law used for another.
To rerun such a study as it was done, this module holds each of those forms as
its own entry, cited to the study that printed it. Where a form differs from
the calibration it names, the entry says how, and ``variant_of`` names the
calibration's entry. None is recommended: they exist so that published
timescales can be reproduced and compared, including forms that later work has
superseded.

Two calibrations that a study used but the registry did not yet hold are here
too, cited to their own papers: Dimanov & Wiedenbeck (2006) for cpx, used by
Kilgour et al. (2014), and Giletti & Shanahan (1997) for Li in plagioclase,
used by Ellis et al. (2018).

Checked on 9 October 2026 against rendered pages of each study; the page of
every equation is in ``equation_number`` and ``verified_from``. Laws that a
study uses exactly as the calibration prints them are not repeated here.
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np

from ..constants import R_GAS
from .base import DiffusionCoefficient, Parameter, Range, LN10
from . import magnetite as _magnetite
from .olivine import _dohmen_chakraborty_tamed

COEFFICIENTS = []

STUDY_NOTE = ("The form a published study used, kept so that its timescales can be reproduced. "
              "It is not a recommendation.")


def _add(c):
    COEFFICIENTS.append(c)
    return c


def _with_latex(builder):
    """Attach a LaTeX builder (diffusor.coefficients.latex) to a law function."""
    def wrap(func):
        func.latex = builder
        return func
    return wrap


def _v(c, name):
    return c.params[name].value


def _n(v, digits=6):
    from .latex import num
    return num(v, digits)


_OL_ANISO6 = {"c": 1.0, "a": 1.0 / 6.0, "b": 1.0 / 6.0}
_OPX_ANISO35 = {"c": 1.0, "a": 1.0 / 3.5, "b": 1.0 / 3.5}
_FIXED = "No uncertainty of the law is sampled: the study gives none for it."


# ===========================================================================
# Olivine
# ===========================================================================
@_with_latex(lambda c: (rf"D_c = {_n(_v(c, 'D0'))}\,\left(\frac{{f_{{O_2}}}}{{10^{{-12}}}}\right)^{{1/6}}"
                        rf"\,10^{{3(X_{{Fe}} - 0.14)}}\,\exp\left(-\frac{{{_n(_v(c, 'Q') * 1e3)}}}"
                        rf"{{{_n(_v(c, 'R'))}\,T}}\right)\quad (f_{{O_2}}\ \mathrm{{in\ bar}})"))
def _costa_chakraborty2004(dc, cond, p):
    """Costa & Chakraborty (2004) eq. 4 in m2/s (printed in um2/s)."""
    XFe = np.asarray(cond.x("XFe"), dtype=float)
    return (p["D0"] * 10.0 ** (p["n_fo2"] * (cond.log_fo2_bar - p["log_fo2_ref_bar"]))
            * 10.0 ** (p["m"] * (XFe - p["XFe_ref"])) * np.exp(-p["Q"] * 1e3 / (p["R"] * cond.T_K)))


_add(DiffusionCoefficient(
    key="ol_FeMg_costa_chakraborty2004",
    mineral="olivine", species="Fe-Mg", kind="interdiffusion",
    transported_variable="XFe = Fe/(Fe+Mg), mole fraction",
    label="Olivine Fe-Mg // [001], Costa & Chakraborty (2004) eq. 4 (Chakraborty 1997 with fO2 and Fo terms)",
    citation="costa_chakraborty2004",
    equation_number="4, p. 523",
    equation_text=("D_c [um2/s] = 5.38e3 (fO2 / 1e-12)^(1/6) 10^(3 (86 - Fo)/100) exp(-226,000 / 8.314 T), "
                   "fO2 in bar, Fo in mol%; D_c ~ 6 D_a ~ 6 D_b (p. 523)"),
    func=_costa_chakraborty2004,
    params={
        "D0": Parameter("D0", 5.38e-9, 0.0, "m2/s", "1s", "5.38e3 um2/s as printed"),
        "n_fo2": Parameter("n_fo2", 1.0 / 6.0, 0.0, "-", "1s", "fO2 exponent"),
        "log_fo2_ref_bar": Parameter("log_fo2_ref_bar", -12.0, 0.0, "log10 bar", "1s", "reference fO2"),
        "m": Parameter("m", 3.0, 0.0, "-", "1s", "10^(3(86 - Fo)/100) = 10^(3 (XFe - 0.14))"),
        "XFe_ref": Parameter("XFe_ref", 0.14, 0.0, "-", "1s", "Fo86"),
        "Q": Parameter("Q", 226.0, 0.0, "kJ/mol", "1s", "activation energy"),
        "R": Parameter("R", 8.314, 0.0, "J/(mol K)", "1s", "gas constant as printed"),
    },
    axis_factors=_OL_ANISO6, reference_axis="c", requires=("XFe",),
    needs_fo2=True, fo2_unit="bar",
    T_range=Range(1253.15, 1573.15, "K (980-1300 C, the Chakraborty 1997 experiments; the study applies it "
                                    "at 840-990 C)"),
    X_range=Range(0.0, 0.5, "XFe (not stated)"),
    verified=True,
    verified_from=("Costa & Chakraborty (2004) EPSL 227:517-530, eq. 4 and the anisotropy sentence on p. 523, "
                   "read from the rendered page on 9 October 2026"),
    secondary_citations=("chakraborty1997",),
    variant_of="ol_FeMg_chakraborty1997",
    uncertainty_note=_FIXED,
    notes=(STUDY_NOTE + " The constants are those of Chakraborty (1997) at Fo86, // [001] and fO2 = 1e-12 "
           "bar. The study adds an fO2 term with exponent 1/6 (it says the exponent 'is not yet very well "
           "characterized'), the composition term 10^(3 (86 - Fo)/100), which is 3 (XFe - 0.14) in mole "
           "fractions, and D_c ~ 6 D_a ~ 6 D_b, independent of temperature. The traverse diffusivity is "
           "projected on the three axes with the study's eq. 5. The study applied it at 850 C and about "
           "NNO+1 for the Volcan San Pedro dacite."),
))


@_with_latex(lambda c: (rf"\log_{{10}} D_c = -9.21 - \frac{{201000 + (P - 10^5)\,7\times10^{{-6}}}}{{2.303\,RT}}"
                        rf" + \frac{{1}}{{6}}\log_{{10}}\frac{{f_{{O_2}}}}{{10^{{-7}}}} + 3\,(X_{{Mg}} - 0.9)"
                        r"\quad (P,\ f_{O_2}\ \mathrm{in\ Pa})"))
def _mourey2023(dc, cond, p):
    """Mourey et al. (2023) eq. 2: TaMED with the composition term 10^(3 (XMg - 0.9))."""
    return _dohmen_chakraborty_tamed(dc, cond, p)


_add(DiffusionCoefficient(
    key="ol_FeMg_mourey2023",
    mineral="olivine", species="Fe-Mg", kind="interdiffusion",
    transported_variable="XFe = Fe/(Fe+Mg), mole fraction",
    label="Olivine Fe-Mg // [001], Mourey et al. (2023) eq. 2 (TaMED with 3 (X_Mg - 0.9))",
    citation="mourey2023",
    equation_number="2, p. 4",
    equation_text=("D_c = 10^-9.21 (fO2 / 1e-7)^(1/6) 10^(3 (X_Mg - 0.9)) exp(-(201,000 + (P - 1e5) .7e-6) / RT), "
                   "D in m2/s, P and fO2 in Pa"),
    func=_mourey2023,
    params={
        "c0": Parameter("c0", -9.21, 0.0, "log10(m2/s)", "1s", "intercept"),
        "Q": Parameter("Q", 201.0, 0.0, "kJ/mol", "1s", "activation energy"),
        "dV": Parameter("dV", 7.0e-6, 0.0, "m3/mol", "1s", "activation volume"),
        "n_fo2": Parameter("n_fo2", 1.0 / 6.0, 0.0, "-", "1s", "fO2 exponent"),
        "log_fo2_ref": Parameter("log_fo2_ref", -7.0, 0.0, "log10 Pa", "1s", "reference fO2"),
        "m": Parameter("m", -3.0, 0.0, "-", "1s", "3 (X_Mg - 0.9) = -3 (X_Fe - 0.1)"),
        "XFe_ref": Parameter("XFe_ref", 0.1, 0.0, "-", "1s", "X_Mg = 0.9"),
    },
    axis_factors=_OL_ANISO6, reference_axis="c", requires=("XFe",),
    needs_fo2=True, fo2_unit="Pa",
    T_range=Range(973.15, 1473.15, "K (700-1200 C, as the Dohmen & Chakraborty 2007 entry)"),
    fo2_range=Range(-10.0, None, "log10 Pa (TaMED regime)"),
    X_range=Range(0.0, 0.5, "XFe"),
    verified=True,
    verified_from="Mourey et al. (2023) Bull. Volcanol. 85:18, eqs 1-3, p. 4, rendered page, 9 October 2026",
    secondary_citations=("dohmen_chakraborty2007",),
    variant_of="ol_FeMg_dohmen_chakraborty2007_tamed",
    uncertainty_note=_FIXED,
    notes=(STUDY_NOTE + " Eq. 2 prints the composition term as 10^(3 (X_Mg - 0.9)). Dohmen & Chakraborty "
           "(2007) as corrected by their erratum give 10^(3 (X_Fe - 0.1)), which has the opposite sign; "
           "this entry keeps the printed form. Eq. 3 of the study gives the fO2 used, "
           "10^(8.912 - 25160/T) 'in Pa', for the QFM buffer; read in bar it lies within 0.1 log units of "
           "Diffusor's QFM. The study ran its models at 60 MPa and along c with EBSD angles; the a and b "
           "factors here are the 1/6 of Dohmen & Chakraborty."),
))


# --- Mutch et al. (2019) multiple linear regressions -------------------------------------------------
_MUTCH_ORDER6 = ("A", "B", "C", "E", "D", "V")
_MUTCH_COV = {
    "global": [[0.497039964543917, 0.00363381377240425, -0.132319267271482, -378.453948506284, -2.76790679156379e-11, 2.69463943024172e-08],
               [0.00363381377240425, 0.000431162727069159, 0.00107609125043613, 10.2162808414763, -6.40646254133509e-13, -1.98771926260489e-10],
               [-0.132319267271482, 0.00107609125043613, 0.149124060624799, 51.0034229198428, -1.46324714807073e-13, -4.70997758933153e-09],
               [-378.453948506284, 10.2162808414763, 51.0034229198428, 839580.703219708, 1.33015759093729e-08, -3.94295445360981e-05],
               [-2.76790679156379e-11, -6.40646254133509e-13, -1.46324714807073e-13, 1.33015759093729e-08, 2.33230090554607e-19, -3.91323758633186e-16],
               [2.69463943024172e-08, -1.98771926260489e-10, -4.70997758933153e-09, -3.94295445360981e-05, -3.91323758633186e-16, 6.61290384810696e-13]],
    "tamed": [[0.719988227888066, 0.0136049442637915, -0.137206914066245, -316.542691119953, -5.10703480738867e-11, 3.56565852966258e-08],
              [0.0136049442637915, 0.000824692307620194, 0.00022482943325554, 11.8157816556034, -1.60520420649397e-12, 2.06662706021798e-10],
              [-0.137206914066245, 0.00022482943325554, 0.133506645898472, 44.4762047396756, 1.76263206439128e-12, -5.04771041057526e-09],
              [-316.542691119953, 11.8157816556034, 44.4762047396756, 819871.010459364, 8.11565428517953e-09, -3.61035571275858e-05],
              [-5.10703480738867e-11, -1.60520420649397e-12, 1.76263206439128e-12, 8.11565428517953e-09, 2.07869617518516e-19, -3.45519816705532e-16],
              [3.56565852966258e-08, 2.06662706021798e-10, -5.04771041057526e-09, -3.61035571275858e-05, -3.45519816705532e-16, 5.82629996012315e-13]],
    "Ni": [[3.32788082296725, 0.0108660968622702, -1.77375521776915, -2186.61589865981, -1.40203650629739e-10, 1.90307239581102e-07],
           [0.0108660968622702, 0.00216784683169149, -0.0152818744193823, 85.006525908944, -1.98169090924804e-12, -1.97548792434352e-09],
           [-1.77375521776915, -0.0152818744193823, 1.87806040681595, -339.86176883119, 2.67755733941985e-11, -2.60844977093618e-08],
           [-2186.61589865981, 85.006525908944, -339.86176883119, 6793292.28743474, 9.49540857947479e-08, -0.000320505347787139],
           [-1.40203650629739e-10, -1.98169090924804e-12, 2.67755733941985e-11, 9.49540857947479e-08, 2.22508365973983e-19, -3.68519013721069e-16],
           [1.90307239581102e-07, -1.97548792434352e-09, -2.60844977093618e-08, -0.000320505347787139, -3.68519013721069e-16, 6.24915157004132e-13]],
    "Mn": [[3.23909781874242, 0.00394323922982214, -0.678808849146587, -3680.51039678914, -1.94956354613448e-10, 2.68546039329982e-07],
           [0.00394323922982214, 0.00347810091754163, 0.002784587490056, 119.28822156565, -4.0328467747652e-12, -2.03960098719852e-09],
           [-0.678808849146587, 0.002784587490056, 0.323011631927977, 726.306288511926, 2.81952079772358e-11, -5.3690017575469e-08],
           [-3680.51039678914, 119.28822156565, 726.306288511926, 8787975.0418876, 9.60676117507748e-08, -0.00039878531957198],
           [-1.94956354613448e-10, -4.0328467747652e-12, 2.81952079772358e-11, 9.60676117507748e-08, 2.83349563026499e-19, -4.64960532776102e-16],
           [2.68546039329982e-07, -2.03960098719852e-09, -5.3690017575469e-08, -0.00039878531957198, -4.64960532776102e-16, 7.86972428160841e-13]],
    # Table S3, order A, B, F, E
    "Ni_aSiO2": [[21.5016970782961, 0.0452483455229286, 0.201837208485653, -34216.6589053779],
                 [0.0452483455229286, 0.00103766450588279, 0.00109265975861079, -58.0799716087498],
                 [0.201837208485653, 0.00109265975861079, 0.022595903294646, -223.229808610946],
                 [-34216.6589053779, -58.0799716087498, -223.229808610946, 55159024.0112905]],
    "Mn_aSiO2": [[6.09147555926196, 0.0046841103351481, 0.0500946657327513, -9812.57551228311],
                 [0.0046841103351481, 0.000133225885263336, 6.49692626413953e-05, -4.7258323253643],
                 [0.0500946657327513, 6.49692626413953e-05, 0.00776304656000876, -46.4574183809125],
                 [-9812.57551228311, -4.7258323253643, -46.4574183809125, 16061862.0217955]],
}


def _mutch_latex(c):
    s = rf"\ln D_{{[001]}} = {_n(_v(c, 'A'))}{' + ' if _v(c, 'B') >= 0 else ' - '}{_n(abs(_v(c, 'B')))}\,\ln f_{{O_2}}"
    if "C" in c.params:
        s += rf" - {_n(abs(_v(c, 'C')))}\,X_{{Fo}}" if _v(c, 'C') < 0 else rf" + {_n(_v(c, 'C'))}\,X_{{Fo}}"
    if "V" in c.params:
        s += (rf" + \frac{{{_n(_v(c, 'E'))} {'+' if _v(c, 'V') >= 0 else '-'} {_n(abs(_v(c, 'V')))}\,P}}{{T}}"
              rf" {'+' if _v(c, 'D') >= 0 else '-'} {_n(abs(_v(c, 'D')))}\,P")
    else:
        s += rf" {'+' if _v(c, 'E') >= 0 else '-'} \frac{{{_n(abs(_v(c, 'E')))}}}{{T}}"
        if "D" in c.params:
            s += rf" {'+' if _v(c, 'D') >= 0 else '-'} {_n(abs(_v(c, 'D')))}\,P"
    if "F" in c.params:
        s += rf" + {_n(_v(c, 'F'))}\,\ln a_{{SiO_2}}"
    return s + r"\quad (P\ \mathrm{in\ Pa},\ f_{O_2}\ \mathrm{in\ bar})"


@_with_latex(_mutch_latex)
def _mutch2019(dc, cond, p):
    """Mutch et al. (2019) eq. 1: ln D = a + b ln fO2 + c XFo + (q + h P)/T + j P + k ln aSiO2."""
    XFo = 1.0 - np.asarray(cond.x("XFe"), dtype=float)
    P = cond.P_Pa
    lnD = p["A"] + p["B"] * cond.log_fo2_bar * LN10 + p["E"] / cond.T_K
    if "C" in p:
        lnD = lnD + p["C"] * XFo
    if "V" in p:
        lnD = lnD + p["V"] * P / cond.T_K
    if "D" in p:
        lnD = lnD + p["D"] * P
    if "F" in p:
        lnD = lnD + p["F"] * np.log(p["aSiO2"])
    return np.exp(lnD)


_MUTCH_ROWS = {
    # key suffix: (species, Table S1 row, values in column order A B C E D V F, covariance key)
    "ol_FeMg_mutch2019_global": ("Fe-Mg", "FeMg (Global)", (-7.861, 0.187, -7.208, -26600, -4.15e-10, -1.54e-07, None), "global"),
    "ol_FeMg_mutch2019_tamed": ("Fe-Mg", "FeMg (TaMED)", (-6.755, 0.224, -7.181, -26700, -5.21e-10, -1.03e-07, None), "tamed"),
    "ol_Ni_mutch2019": ("Ni", "Ni", (-11.1, 0.277, -2.185, -25100, -1.25e-09, 9.97e-07, None), "Ni"),
    "ol_Mn_mutch2019": ("Mn", "Mn", (-7.55, 0.196, -7.153, -26700, -9.5e-10, 7.2e-07, None), "Mn"),
    "ol_Ni_mutch2019_asio2": ("Ni", "Ni (aSiO2)", (-14.444, -0.1073, None, -32980, None, None, 0.714), "Ni_aSiO2"),
    "ol_Mn_mutch2019_asio2": ("Mn", "Mn (aSiO2)", (-7.463, -0.09655, None, -44310, None, None, 0.7613), "Mn_aSiO2"),
}
_MUTCH_DESC = {"A": ("a, intercept", "ln(m2/s)"), "B": ("b, ln fO2 coefficient", "-"),
               "C": ("c, X_Fo coefficient", "-"), "E": ("q, Table S1 column E", "K"),
               "D": ("j, Table S1 column D, multiplies P", "1/Pa"),
               "V": ("h, Table S1 column V, multiplies P/T", "K/Pa"),
               "F": ("k, ln a_SiO2 coefficient", "-")}
for _key, (_species, _row, _vals, _cov) in _MUTCH_ROWS.items():
    _params = {name: Parameter(name, value, 0.0, _MUTCH_DESC[name][1], "1s", _MUTCH_DESC[name][0])
               for name, value in zip(("A", "B", "C", "E", "D", "V", "F"), _vals) if value is not None}
    _order = _MUTCH_ORDER6 if "V" in _params else ("A", "B", "F", "E")
    if "F" in _params:
        _params["aSiO2"] = Parameter("aSiO2", 0.62, 0.0, "-", "1s",
                                     "silica activity, held at the Borgarhraun value of the study (0.62 +/- 0.03, p. 7)")
    _kind = "interdiffusion" if _species == "Fe-Mg" else "chemical"
    _add(DiffusionCoefficient(
        key=_key, mineral="olivine", species=_species, kind=_kind,
        transported_variable=("XFe = Fe/(Fe+Mg), mole fraction" if _species == "Fe-Mg"
                              else f"{_species} concentration"),
        label=f"Olivine {_species} // [001], Mutch et al. (2019) regression, Table S1 row '{_row}'",
        citation="mutch2019",
        equation_number="1, p. 7; Supplementary Tables 1-3",
        equation_text=("ln D[001] = a + b ln fO2 + c X_Fo + (q + h P)/T + j P + k ln aSiO2, P in Pa, T in K, fO2 in "
                       f"bar; row '{_row}' of Supplementary Table 1 (columns A, B, C, E, D, V, F): "
                       + ", ".join(f"{n} = {v:g}" for n, v in zip(("A", "B", "C", "E", "D", "V", "F"), _vals)
                                   if v is not None)),
        func=_mutch2019,
        params=_params,
        covariance=np.array(_MUTCH_COV[_cov]), cov_order=_order,
        axis_factors=_OL_ANISO6, reference_axis="c", requires=("XFe",),
        needs_fo2=True, fo2_unit="bar",
        T_range=Range(None, None, "K (the compiled experiments; the paper prints no range)"),
        X_range=Range(0.0, 0.5, "XFe (not stated)"),
        verified=True,
        verified_from=("Mutch et al. (2019) Nature Geosci. 12:569-574, Methods eq. 1 (p. 7, rendered page) and "
                       "Supplementary Tables 1-3 (spreadsheets), read 9 October 2026"),
        secondary_citations=("mutch2021", "dohmen_chakraborty2007", "petry2004"),
        variant_of=("ol_FeMg_dohmen_chakraborty2007_tamed" if _species == "Fe-Mg"
                    else "ol_Ni_petry2004" if _species == "Ni" else ""),
        uncertainty_note=("Published covariance matrix (Supplementary Table 2 or 3) of the regression parameters, "
                          "in the parameters' own units; it is sampled in the Monte Carlo."),
        notes=(STUDY_NOTE + " A multiple linear regression through compiled olivine experiments, along [001]. "
               "Eq. 1 names the parameters a, b, c, q, h, j and k; Supplementary Table 1 heads its columns A, "
               "B, C, E, D, V and F. Diffusor reads E as q, V as h (it multiplies P/T) and D as j (it multiplies "
               "P): taking the columns in the order of the equation instead would put log D near -50 at 0.8 GPa, "
               "and the variances of Table 2 have the sizes of those readings. Anisotropy: 'six times faster "
               "along [001] than the [010] and [100] axes for each element' (p. 7). Mutch et al. (2021, DFENS, "
               "p. 6) refer to these regressions and use 10.7 for the Ni anisotropy instead (Spandler & O'Neill "
               "2010); set the traverse orientation accordingly to follow that choice."
               + (" The aSiO2 regression covers the experiments buffered for silica activity; a_SiO2 is held at "
                  "0.62, the Borgarhraun value. The covariance is the block of Supplementary Table 3. Mutch et "
                  "al. (2021, Table S4) print the same block as the one for fO2 in Pa and another block, with "
                  "a smaller variance of a, as the one for fO2 in bar, both for these parameters; Diffusor keeps "
                  "the 2019 block here." if "F" in _params else "")),
    ))


# --- Mutch et al. (2021) DFENS, Supporting Information Tables S1-S6 and data set S6 ----------------------
# Version 1 is the 2019 regression with more digits (Table S1 of the 2021 SI; the covariance of
# Table S2 is the 2019 one). Version 2 drops the P/T term. Table S1 of the SI rounds c to two
# decimals; the values here are those of data set S6 (sheet 'Olivine Regression Parameters').
_DFENS_COV_V2 = {
    "global": [[0.49159595420422, 0.00360999888813302, -0.130969497867373, -373.544911886026, -1.16305078479025e-11],
               [0.00360999888813302, 0.000427325177531545, 0.00106525801796523, 10.1150065124044, -7.51626284457795e-13],
               [-0.130969497867373, 0.00106525801796523, 0.147784017870059, 50.2781021859257, -2.9077843016064e-12],
               [-373.544911886026, 10.1150065124044, 50.2781021859257, 829892.975349277, -9.94326134150935e-09],
               [-1.16305078479025e-11, -7.51626284457795e-13, -2.9077843016064e-12, -9.94326134150935e-09, 1.64663773862295e-21]],
    "tamed": [[0.710758186702312, 0.0134588387422847, -0.135553844965028, -311.246851436954, -2.96309345291454e-11],
              [0.0134588387422847, 0.000816522364177152, 0.00022439478988668, 11.7124471558664, -1.46808847562356e-12],
              [-0.135553844965028, 0.00022439478988668, 0.132152489650859, 43.7297908884151, -1.21875033407287e-12],
              [-311.246851436954, 11.7124471558664, 43.7297908884151, 809605.744833811, -1.31644696860565e-08],
              [-2.96309345291454e-11, -1.46808847562356e-12, -1.21875033407287e-12, -1.31644696860565e-08, 2.93526358643008e-21]],
    "Ni": [[3.32503120827199, 0.0116609530670394, -1.79556930965568, -2124.21575282372, -2.84489580755306e-11],
           [0.0116609530670394, 0.00219802948765591, -0.0156232546662106, 85.4088076370308, -3.19968468802249e-12],
           [-1.79556930965568, -0.0156232546662106, 1.90860257974076, -359.192782501142, 1.15852727286911e-11],
           [-2124.21575282372, 85.4088076370308, -359.192782501142, 6740623.35062343, -9.56362388525762e-08],
           [-2.84489580755306e-11, -3.19968468802249e-12, 1.15852727286911e-11, -9.56362388525762e-08, 5.27628659851697e-21]],
    "Mn": [[3.1218164850397, 0.00460143502452992, -0.655106603108254, -3515.55210916565, -3.59977891835702e-11],
           [0.00460143502452992, 0.0034445213199113, 0.00262388596870162, 117.29124839146, -5.19521380386978e-12],
           [-0.655106603108254, 0.00262388596870162, 0.316746928121797, 693.404083925858, -3.49730226790862e-12],
           [-3515.55210916565, 117.29124839146, 693.404083925858, 8515946.53521086, -1.38406600114609e-07],
           [-3.59977891835702e-11, -5.19521380386978e-12, -3.49730226790862e-12, -1.38406600114609e-07, 8.57030785833625e-21]],
}
_DFENS_OL_ROWS = {
    # key: (species, row, version, values a b c q j h in Diffusor's names A B C E D V, covariance)
    "ol_FeMg_mutch2021_global": ("Fe-Mg", "FeMg (Global)", 1,
                                 (-7.861, 0.1868, -7.208, -26580, -4.148e-10, -1.54e-07), _MUTCH_COV["global"]),
    "ol_FeMg_mutch2021_tamed": ("Fe-Mg", "FeMg (TaMED)", 1,
                                (-6.755, 0.2244, -7.181, -26740, -5.213e-10, -1.028e-07), _MUTCH_COV["tamed"]),
    "ol_Ni_mutch2021": ("Ni", "Ni", 1, (-11.09, 0.2769, -2.185, -25080, -1.246e-09, 9.967e-07), _MUTCH_COV["Ni"]),
    "ol_Mn_mutch2021": ("Mn", "Mn", 1, (-7.548, 0.1963, -7.153, -26720, -9.504e-10, 7.195e-07), _MUTCH_COV["Mn"]),
    "ol_FeMg_mutch2021_global_v2": ("Fe-Mg", "FeMg (Global)", 2,
                                    (-7.855, 0.187, -7.209, -26590, -5.061e-10, None), _DFENS_COV_V2["global"]),
    "ol_FeMg_mutch2021_tamed_v2": ("Fe-Mg", "FeMg (TaMED)", 2,
                                   (-6.749, 0.2245, -7.182, -26740, -5.822e-10, None), _DFENS_COV_V2["tamed"]),
    "ol_Ni_mutch2021_v2": ("Ni", "Ni", 2, (-11.39, 0.28, -2.143, -24570, -6.584e-10, None), _DFENS_COV_V2["Ni"]),
    "ol_Mn_mutch2021_v2": ("Mn", "Mn", 2, (-7.794, 0.1982, -7.104, -26360, -5.253e-10, None), _DFENS_COV_V2["Mn"]),
}
_DFENS_NI_ANISO = {"c": 1.0, "a": 1.0 / 10.7, "b": 1.0 / 10.7}
_DFENS_VERIFIED = ("Mutch et al. (2021) G-Cubed 22:e2020GC009303, eq. 6 and the anisotropy sentence (p. 6), "
                   "Supporting Information Text S2 eqs 17-18 (pp. 6-7) and Tables S1-S3 (pp. 128-130), rendered "
                   "pages, and data set S6 (spreadsheet), read 13 October 2026")
for _key, (_species, _row, _version, _vals, _cov) in _DFENS_OL_ROWS.items():
    _params = {name: Parameter(name, value, 0.0, _MUTCH_DESC[name][1], "1s", _MUTCH_DESC[name][0]
                               .replace("Table S1 column E", "1/T coefficient")
                               .replace("Table S1 column D, multiplies P", "P coefficient")
                               .replace("Table S1 column V, multiplies P/T", "P/T coefficient"))
               for name, value in zip(("A", "B", "C", "E", "D", "V"), _vals) if value is not None}
    _v2 = _version == 2
    _form = ("ln D[001] = a + b ln fO2 + c X_Fo + q/T + j P" if _v2
             else "ln D[001] = a + b ln fO2 + c X_Fo + (q + h P)/T + j P")
    _add(DiffusionCoefficient(
        key=_key, mineral="olivine", species=_species,
        kind="interdiffusion" if _species == "Fe-Mg" else "chemical",
        transported_variable=("XFe = Fe/(Fe+Mg), mole fraction" if _species == "Fe-Mg"
                              else f"{_species} concentration"),
        label=(f"Olivine {_species} // [001], Mutch et al. (2021) DFENS regression version {_version}, "
               f"row '{_row}'"),
        citation="mutch2021",
        equation_number=("SI eq. 18, p. 7; SI Tables S1 and S3" if _v2
                         else "6, p. 6; SI eq. 17 and Tables S1-S2"),
        equation_text=(_form + ", P in Pa, T in K, fO2 in bar; data set S6 row '" + _row + f"' (version {_version}): "
                       + ", ".join(f"{n} = {v:g}" for n, v in zip("abcqjh", _vals) if v is not None)),
        func=_mutch2019,
        params=_params,
        covariance=np.array(_cov), cov_order=_MUTCH_ORDER6[:5] if _v2 else _MUTCH_ORDER6,
        axis_factors=_DFENS_NI_ANISO if _species == "Ni" else _OL_ANISO6, reference_axis="c",
        requires=("XFe",), needs_fo2=True, fo2_unit="bar",
        T_range=Range(None, None, "K (the compiled experiments; the paper prints no range)"),
        fo2_range=(Range(-15.0, None, "log10 bar (data set S6: TaMED for fO2 > 1e-15 bar)")
                   if "TaMED" in _row else Range()),
        X_range=Range(0.0, 0.5, "XFe (not stated)"),
        verified=True,
        verified_from=_DFENS_VERIFIED,
        secondary_citations=("mutch2019",),
        variant_of=_key.replace("mutch2021", "mutch2019").replace("_v2", ""),
        uncertainty_note=(f"Published covariance matrix of the regression parameters (SI Table "
                          f"{'S3' if _v2 else 'S2'}, full digits in data set S6), sampled in the Monte Carlo."),
        notes=(STUDY_NOTE + (" Version 2: the regression refitted without the P/T term (SI eq. 18); the paper "
                             "used version 1 for its models (Table S1 caption)." if _v2 else
                             " Version 1, the regression DFENS used (Table S1 caption). It is the Mutch et al. "
                             "(2019) regression with more digits: the 2019 table prints b, q, j and h rounded.")
               + " Data set S6 also gives each regression with fO2 in Pa, which is the same law with a shifted "
               "intercept; the authors recommend the bar form. "
               + ("Anisotropy: D[001] is 10.7 times D[100] and D[010] for Ni (p. 6, after Spandler & O'Neill "
                  "2010), which DFENS used instead of the factor 6 of Mutch et al. (2019)."
                  if _species == "Ni" else
                  "Anisotropy: D[001] is six times D[100] and D[010] (p. 6, after Chakraborty 2010).")
               + (" Data set S6 restricts the TaMED regression to fO2 > 1e-15 bar; the SI text says it covers "
                  "experiments at fO2 > 1e-10 Pa." if "TaMED" in _row else "")),
    ))


# DFENS plagioclase regressions, ln D = a + b X_An + c ln aSiO2 + q/T (eq. 7, SI eq. 20, Table S5)
def _dfens_plag_latex(c):
    s = rf"\ln D = {_n(_v(c, 'a'))}"
    if "b" in c.params:
        s += rf" {'+' if _v(c, 'b') >= 0 else '-'} {_n(abs(_v(c, 'b')))}\,X_{{An}}"
    if "c" in c.params:
        s += rf" + {_n(_v(c, 'c'))}\,\ln a_{{SiO_2}}"
    return s + rf" - \frac{{{_n(abs(_v(c, 'q')))}}}{{T}}\quad (D\ \mathrm{{in\ m^2/s}})"


@_with_latex(_dfens_plag_latex)
def _dfens_plag(dc, cond, p):
    """Mutch et al. (2021) eq. 7: ln D = a + b X_An + c ln aSiO2 + q/T, D in m2/s."""
    lnD = p["a"] + p["q"] / cond.T_K
    if "b" in p:
        lnD = lnD + p["b"] * np.asarray(cond.x("XAn"), dtype=float)
    if "c" in p:
        lnD = lnD + p["c"] * np.log(p["aSiO2"])
    return np.exp(lnD)


_DFENS_PLAG_ROWS = {
    # key: (species, row, (a, b, c, q), covariance over the present parameters in the order a, b, c, q,
    #       data set, calibration entry)
    "plag_Mg_mutch2021_f13": ("Mg", "Mg (F13)", (-11.77, None, 2.931, -34080),
                              [[25.3238145277826, -0.317911855537308, -36291.6915744149],
                               [-0.317911855537308, 0.199038284600889, 582.771314966343],
                               [-36291.6915744149, 582.771314966343, 52104130.1004335]],
                              "Faak et al. (2013)", ""),
    "plag_Mg_mutch2021_vo14": ("Mg", "Mg (VO14)", (-5.45, -7.983, None, -35360),
                               [[1.01310337218217, -0.206341842567178, -1105.15309042578],
                                [-0.206341842567178, 0.139790583067871, 143.213031775432],
                                [-1105.15309042578, 143.213031775432, 1281622.31678333]],
                               "Van Orman et al. (2014)", "plag_Mg_vanorman2014"),
    "plag_Mg_mutch2021_vof": ("Mg", "Mg (VO+F)", (-8.727, -6.125, 3.712, -32860),
                              [[1.06395533788672, -0.160078835464304, 0.19920338670931, -1236.87710775558],
                               [-0.160078835464304, 0.146591197055639, -0.0341384763646049, 83.337010235268],
                               [0.19920338670931, -0.0341384763646049, 0.0778839599467262, -211.998447693381],
                               [-1236.87710775558, 83.337010235268, -211.998447693381, 1536028.15994034]],
                              "Faak et al. (2013) and Van Orman et al. (2014), with aSiO2 = 1 assigned to the "
                              "Van Orman experiments (p. 7)", ""),
    "plag_Sr_mutch2021_combined": ("Sr", "Sr (Combined)", (-12.81, -5.712, None, -32390),
                                   [[0.947743129406916, -0.164931535260681, -1027.09737511478],
                                    [-0.164931535260681, 0.117454496730595, 123.581073990204],
                                    [-1027.09737511478, 123.581073990204, 1155796.40125077]],
                                   "Cherniak & Watson (1994) and Giletti & Casserly (1994)", ""),
    "plag_Sr_mutch2021_cw": ("Sr", "Sr (C+W)", (-13.42, -4.001, None, -32540),
                             [[1.60892175695837, -0.115484252667338, -1826.58831952283],
                              [-0.115484252667338, 0.251638910826196, 6.40860709556696],
                              [-1826.58831952283, 6.40860709556696, 2148695.00986041]],
                             "Cherniak & Watson (1994)", "plag_Sr_cherniak_watson1994"),
    "plag_Sr_mutch2021_gc": ("Sr", "Sr (G+C)", (-9.175, -8.021, None, -34920),
                             [[1.03781326741931, -0.270436704781618, -1048.06796268207],
                              [-0.270436704781618, 0.17237850457254, 183.867126777445],
                              [-1048.06796268207, 183.867126777445, 1148528.70913514]],
                             "Giletti & Casserly (1994)", "plag_Sr_giletti_casserly1994"),
    "plag_Ba_mutch2021": ("Ba", "Ba", (-12.32, -3.287, None, -39950),
                          [[2.54171087173606, -0.15103024911905, -2956.23933859547],
                           [-0.15103024911905, 0.304782062277467, -0.0512040884890124],
                           [-2956.23933859547, -0.0512040884890124, 3563137.39869231]],
                          "Cherniak (2002)", "plag_Ba_cherniak2002"),
    "plag_K_mutch2021": ("K", "K", (-9.079, -3.855, None, -33960),
                         [[0.620746461992886, -0.095302644457305, -635.156310256245],
                          [-0.095302644457305, 0.150926960891201, 66.207898069129],
                          [-635.156310256245, 66.207898069129, 668097.605313518]],
                         "Giletti & Shanahan (1997)", ""),
}
_DFENS_PLAG_REFS = {"Mg (F13)": ("faak2013",), "Mg (VO14)": ("vanorman2014",),
                    "Mg (VO+F)": ("faak2013", "vanorman2014"),
                    "Sr (Combined)": ("cherniak_watson1994", "giletti_casserly1994"),
                    "Sr (C+W)": ("cherniak_watson1994",), "Sr (G+C)": ("giletti_casserly1994",),
                    "Ba": ("cherniak2002",), "K": ("giletti_shanahan1997",)}
_DFENS_PLAG_DESC = {"a": ("a, intercept", "ln(m2/s)"), "b": ("b, X_An coefficient", "-"),
                    "c": ("c, ln aSiO2 coefficient", "-"), "q": ("q, 1/T coefficient", "K")}
for _key, (_species, _row, _vals, _cov, _data, _calib) in _DFENS_PLAG_ROWS.items():
    _params = {name: Parameter(name, value, 0.0, _DFENS_PLAG_DESC[name][1], "1s", _DFENS_PLAG_DESC[name][0])
               for name, value in zip("abcq", _vals) if value is not None}
    _order = tuple(_params)
    if "c" in _params:
        _params["aSiO2"] = Parameter("aSiO2", 1.0, 0.0, "-", "1s",
                                     "silica activity; 1 (quartz saturation) is Diffusor's default, set the "
                                     "value of the system")
    _add(DiffusionCoefficient(
        key=_key, mineral="plagioclase", species=_species,
        kind="effective" if _key == "plag_Sr_mutch2021_combined" else "chemical",
        transported_variable=f"{_species} concentration",
        label=f"Plagioclase {_species}, Mutch et al. (2021) DFENS regression, row '{_row}'",
        citation="mutch2021",
        equation_number="7, p. 6; SI eq. 20 and Tables S5-S6",
        equation_text=("ln D = a + b X_An + c ln aSiO2 + q/T, D in m2/s, T in K, X_An in mole fraction; row '"
                       + _row + "': " + ", ".join(f"{n} = {v:g}" for n, v in zip("abcq", _vals)
                                                  if v is not None)),
        func=_dfens_plag,
        params=_params,
        covariance=np.array(_cov), cov_order=_order,
        requires=("XAn",) if "b" in _params else (),
        needs_fo2=False,
        T_range=Range(None, None, f"K (the experiments of {_data.split(',')[0]})"),
        X_range=Range(0.0, 1.0, "X_An (not stated)"),
        verified=True,
        verified_from=("Mutch et al. (2021) G-Cubed 22:e2020GC009303, eq. 7 (p. 6, rendered page), SI eq. 20 (p. 7) "
                       "and Tables S5-S6 (pp. 131-132), rendered pages, and data set S6 (spreadsheet), read "
                       "13 October 2026" + ("; the K row is in Table S6 and data set S6, not in Table S5"
                                            if _species == "K" else "")),
        secondary_citations=_DFENS_PLAG_REFS[_row],
        variant_of=_calib,
        uncertainty_note="Published covariance matrix (SI Table S6, full digits in data set S6), sampled in the "
                         "Monte Carlo.",
        notes=(STUDY_NOTE + f" A regression through the data of {_data}, made for DFENS. "
               "It includes all data regardless of crystallographic direction (p. 7). "
               + ("DFENS models Mg with the activity term of eq. 15 of the SI and A_Mg = -34.099 kJ/mol, "
                  "the 'mutch2021' activity set. The study bases its interpretation on the F13 regression "
                  "(p. 18) and assigned aSiO2 a Gaussian prior from an independent estimate (p. 8); the "
                  "value is not printed in the text. "
                  if _species == "Mg" else "")
               + ("Combines a tracer and a chemical data set, so the transport kind is 'effective'. "
                  if _key == "plag_Sr_mutch2021_combined" else "")
               + ("Table S5 of the SI ends with the Ba row; the K row is in Table S6 and data set S6. "
                  if _species == "K" else "")),
    ))


# --- DIPRA (Girona & Costa 2013), Table 1 -------------------------------------------------------------
@_with_latex(lambda c: (rf"D_{{[001]}} = 3.84\times10^{{-9}}\,\left(\frac{{f_{{O_2}}}}{{10^{{-6}}}}\right)^{{1/4.25}}"
                        rf"\,10^{{1.5(0.9 - X_{{Fo}})}}\,\exp\left(-\frac{{{_n(_v(c, 'Q') * 1e3)} + 7\times10^{{-6}}"
                        r"\,(P - 10^5)}{RT}\right)\quad (P,\ f_{O_2}\ \mathrm{in\ Pa})"))
def _dipra_ni(dc, cond, p):
    XFo = 1.0 - np.asarray(cond.x("XFe"), dtype=float)
    return (p["D0"] * 10.0 ** ((cond.log_fo2_Pa - p["log_fo2_ref"]) / p["n_inv"])
            * 10.0 ** (p["m"] * (0.9 - XFo))
            * np.exp(-(p["Q"] * 1e3 + p["dV"] * (cond.P_Pa - 1e5)) / (R_GAS * cond.T_K)))


@_with_latex(lambda c: (rf"\log_{{10}} D_{{[001]}} = -9.21 + \frac{{1}}{{6}}\log_{{10}}\frac{{f_{{O_2}}}}{{10^{{-7}}}}"
                        rf" + 3\,(0.9 - X_{{Fo}}) - \frac{{{_n(_v(c, 'Q') * 1e3)} + 7\times10^{{-6}}\,(P - 10^5)}}"
                        r"{2.303\,RT}\quad (P,\ f_{O_2}\ \mathrm{in\ Pa})"))
def _dipra_mn(dc, cond, p):
    return _dohmen_chakraborty_tamed(dc, cond, p)


def _dipra_ca_axis(axis):
    def f(dc, cond, p):
        return (p[f"D0_{axis}"] * 10.0 ** (p["n_fo2"] * (cond.log_fo2_Pa - p["log_fo2_ref"]))
                * np.exp(-p[f"Q_{axis}"] * 1e3 / (R_GAS * cond.T_K)))
    f.latex = lambda c: "; ".join(
        rf"D_{{[{i}]}} = {_n(_v(c, f'D0_{a}'))}\,(f_{{O_2}}/10^{{-7}})^{{0.3}}\,"
        rf"\exp\left(-\frac{{{_n(_v(c, f'Q_{a}') * 1e3)}}}{{RT}}\right)"
        for a, i in (("a", "100"), ("b", "010"), ("c", "001"))) + r"\quad (f_{O_2}\ \mathrm{in\ Pa})"
    return f


_DIPRA_COMMON = dict(
    citation="girona_costa2013", equation_number="Table 1, p. 424",
    verified=True,
    verified_from=("Girona & Costa (2013) G3 14:422-431, Table 1 and its footnote, p. 424, rendered page, "
                   "9 October 2026"),
    uncertainty_note=_FIXED,
)
_DIPRA_USERS = ("DIPRA computes these laws, so they are the forms of the studies that model with it, for example "
                "Gordeychik et al. (2018), Ruth et al. (2018), Rasmussen et al. (2018), Sundermeyer et al. (2020) "
                "and Weller et al. (2026). fO2 is in Pa and DIPRA calculates it from the NNO buffer.")

_add(DiffusionCoefficient(
    key="ol_Mn_girona_costa2013", mineral="olivine", species="Mn", kind="effective",
    transported_variable="Mn concentration",
    label="Olivine Mn // [001], DIPRA (Girona & Costa 2013), the Fe-Mg form with E_Mn",
    equation_text=("D[001] = 10^-9.21 (fO2/1e-7)^(1/6) 10^(3 (0.9 - X_Fo)) exp(-(E_Mn + 7e-6 (p - 1e5))/RT), "
                   "E_Mn = 201,000 J/mol, D in m2/s, fO2 and p in Pa"),
    func=_dipra_mn,
    params={
        "c0": Parameter("c0", -9.21, 0.0, "log10(m2/s)", "1s", "intercept"),
        "Q": Parameter("Q", 201.0, 0.0, "kJ/mol", "1s", "E_Mn"),
        "dV": Parameter("dV", 7.0e-6, 0.0, "m3/mol", "1s", "pressure term"),
        "n_fo2": Parameter("n_fo2", 1.0 / 6.0, 0.0, "-", "1s", "fO2 exponent"),
        "log_fo2_ref": Parameter("log_fo2_ref", -7.0, 0.0, "log10 Pa", "1s", "reference fO2"),
        "m": Parameter("m", 3.0, 0.0, "-", "1s", "3 (0.9 - X_Fo) = 3 (X_Fe - 0.1)"),
        "XFe_ref": Parameter("XFe_ref", 0.1, 0.0, "-", "1s", "Fo90"),
    },
    axis_factors=_OL_ANISO6, reference_axis="c", requires=("XFe",), needs_fo2=True, fo2_unit="Pa",
    secondary_citations=("dohmen_chakraborty2007", "chakraborty2010"),
    variant_of="ol_FeMg_dohmen_chakraborty2007_tamed",
    notes=(STUDY_NOTE + " DIPRA takes the Fe-Mg interdiffusion law of Dohmen & Chakraborty (2007) for Mn, "
           "with the same activation energy, and lets it depend on the Fo content of the host. "
           + _DIPRA_USERS),
    **_DIPRA_COMMON,
))

_add(DiffusionCoefficient(
    key="ol_Ni_girona_costa2013", mineral="olivine", species="Ni", kind="chemical",
    transported_variable="Ni concentration",
    label="Olivine Ni // [001], DIPRA (Girona & Costa 2013), with fO2, Fo and pressure terms",
    equation_text=("D[001] = 3.84e-9 (fO2/1e-6)^(1/4.25) 10^(1.5 (0.9 - X_Fo)) exp(-(E_Ni + 7e-6 (p - 1e5))/RT), "
                   "E_Ni = 220,000 J/mol, D in m2/s, fO2 and p in Pa"),
    func=_dipra_ni,
    params={
        "D0": Parameter("D0", 3.84e-9, 0.0, "m2/s", "1s", "pre-exponential factor"),
        "log_fo2_ref": Parameter("log_fo2_ref", -6.0, 0.0, "log10 Pa", "1s", "reference fO2"),
        "n_inv": Parameter("n_inv", 4.25, 0.0, "-", "1s", "fO2 exponent 1/4.25"),
        "m": Parameter("m", 1.5, 0.0, "-", "1s", "Fo term"),
        "Q": Parameter("Q", 220.0, 0.0, "kJ/mol", "1s", "E_Ni"),
        "dV": Parameter("dV", 7.0e-6, 0.0, "m3/mol", "1s", "pressure term"),
    },
    axis_factors=_OL_ANISO6, reference_axis="c", requires=("XFe",), needs_fo2=True, fo2_unit="Pa",
    secondary_citations=("petry2004", "holzapfel2007"),
    variant_of="ol_Ni_petry2004",
    notes=(STUDY_NOTE + " The Petry et al. (2004) pre-exponential factor with an fO2 term referred to 1e-6 Pa, "
           "the exponent 1/4.25, a Fo term 10^(1.5 (0.9 - X_Fo)), the 7e-6 m3/mol pressure term and "
           "E_Ni = 220 kJ/mol; Diffusor's Petry entry is the fixed-fO2 fit with 216 kJ/mol. " + _DIPRA_USERS),
    **_DIPRA_COMMON,
))

_dipra_ca = {a: _dipra_ca_axis(a) for a in ("a", "b", "c")}
_add(DiffusionCoefficient(
    key="ol_Ca_girona_costa2013", mineral="olivine", species="Ca", kind="chemical",
    transported_variable="Ca concentration",
    label="Olivine Ca, DIPRA (Girona & Costa 2013), three axis laws with (fO2/1e-7 Pa)^0.3",
    equation_text=("D[100] = 16.59e-12 (fO2/1e-7)^0.3 exp(-193,000/RT); D[010] = 34.67e-12 (fO2/1e-7)^0.3 "
                   "exp(-201,000/RT); D[001] = 95.49e-12 (fO2/1e-7)^0.3 exp(-207,000/RT); D in m2/s, fO2 in Pa"),
    func=_dipra_ca["c"], principal_funcs=_dipra_ca, reference_axis="c", orientation_required=True,
    params={
        "D0_a": Parameter("D0_a", 16.59e-12, 0.0, "m2/s", "1s", "[100]"),
        "Q_a": Parameter("Q_a", 193.0, 0.0, "kJ/mol", "1s", "[100]"),
        "D0_b": Parameter("D0_b", 34.67e-12, 0.0, "m2/s", "1s", "[010]"),
        "Q_b": Parameter("Q_b", 201.0, 0.0, "kJ/mol", "1s", "[010]"),
        "D0_c": Parameter("D0_c", 95.49e-12, 0.0, "m2/s", "1s", "[001]"),
        "Q_c": Parameter("Q_c", 207.0, 0.0, "kJ/mol", "1s", "[001]"),
        "n_fo2": Parameter("n_fo2", 0.3, 0.0, "-", "1s", "fO2 exponent"),
        "log_fo2_ref": Parameter("log_fo2_ref", -7.0, 0.0, "log10 Pa", "1s", "reference fO2"),
    },
    needs_fo2=True, fo2_unit="Pa",
    secondary_citations=("coogan2005ca",),
    variant_of="ol_Ca_coogan2005",
    notes=(STUDY_NOTE + " The pre-exponential factors and activation energies are those of Coogan et al. "
           "(2005); the fO2 exponent is 0.3, where Diffusor's Coogan entry uses 0.31, and the reference "
           "1e-7 Pa equals Coogan's 1e-12 bar. " + _DIPRA_USERS),
    **_DIPRA_COMMON,
))


# ===========================================================================
# Orthopyroxene
# ===========================================================================
@_with_latex(lambda c: (rf"D = {_n(_v(c, 'D0'))}\,\exp\left(-\frac{{{_n(_v(c, 'Q'))}\ \mathrm{{kJ\,mol^{{-1}}}}}}{{RT}}\right)"
                        rf"\,10^{{m\,(X_{{Fe}} - {_n(_v(c, 'XFe_ref'))})}},\ m = {_n(_v(c, 'm_slope'))}/T + "
                        rf"{_n(_v(c, 'm_int'))}"))
def _winslow2026(dc, cond, p):
    XFe = np.asarray(cond.x("XFe"), dtype=float)
    m = p["m_slope"] / cond.T_K + p["m_int"]
    return p["D0"] * np.exp(-p["Q"] * 1e3 / (R_GAS * cond.T_K)) * 10.0 ** (m * (XFe - p["XFe_ref"]))


_add(DiffusionCoefficient(
    key="opx_FeMg_winslow2026",
    mineral="opx", species="Fe-Mg", kind="interdiffusion", transported_variable="XFe = Fe/(Fe+Mg)",
    label="Opx Fe-Mg, Winslow et al. (2026) eqs 2-3 (Dias & Dohmen 2024 referred to X_Fe = 0.09)",
    citation="winslow2026", equation_number="2-3, p. 16",
    equation_text=("D [m2/s] = 3.8e-9 exp[-261.07 +/- 24 kJ/mol / RT]; D(X_Fe) = D(X_Fe = 0.09) 10^(m (X_Fe - 0.09)), "
                   "m = -2.711e4/T + 23.5408"),
    func=_winslow2026,
    params={
        "D0": Parameter("D0", 3.8e-9, 0.0, "m2/s", "1s", "pre-exponential factor"),
        "Q": Parameter("Q", 261.07, 24.0, "kJ/mol", "unstated", "activation energy"),
        "m_slope": Parameter("m_slope", -27110.0, 0.0, "K", "1s", "slope of m on 1/T"),
        "m_int": Parameter("m_int", 23.5408, 0.0, "-", "1s", "intercept of m"),
        "XFe_ref": Parameter("XFe_ref", 0.09, 0.0, "-", "1s", "reference composition as printed"),
    },
    axis_factors={"c": 1.0, "b": 1.0, "a": 1.0 / 3.5}, reference_axis="c", requires=("XFe",),
    T_range=Range(1223.15, 1373.15, "K (950-1100 C, the Dias & Dohmen 2024 experiments; the study uses 1175 C)"),
    X_range=Range(0.09, 0.5, "XFe"),
    verified=True,
    verified_from="Winslow et al. (2026) Bull. Volcanol. 88:5, eqs 2-3 and the anisotropy sentence, p. 16, rendered page, 9 October 2026",
    secondary_citations=("dias_dohmen2024", "dohmen2016"),
    variant_of="opx_FeMg_dias_dohmen2024",
    uncertainty_note="The study propagates only the thermometer uncertainty (+/-13 C, 1 sigma); the law is held fixed.",
    notes=(STUDY_NOTE + " The Dias & Dohmen (2024) law with the composition correction referred to X_Fe = 0.09, "
           "as the study prints it; Diffusor's Dias & Dohmen entry refers it to 0.1 (their eq. 7). The study "
           "applies it at 1175 +/- 13 C to Mauna Loa 2022 enstatite, with diffusion along a 3.5 times slower "
           "than along b or c (Dohmen et al. 2016). Grain orientations were not measured, so the study gives a "
           "minimum (b/c) and a maximum (a) time for each traverse."),
))


@_with_latex(lambda c: (rf"D_c = {_n(_v(c, 'D0'))}\,f_{{O_2}}^{{{_n(_v(c, 'n'))}}}\,"
                        rf"\exp\left(-\frac{{{_n(_v(c, 'Q'))}\ \mathrm{{kJ\,mol^{{-1}}}}}}{{RT}}\right)"
                        r"\quad (f_{O_2}\ \mathrm{in\ Pa})"))
def _dohmen2016_no_composition(dc, cond, p):
    return p["D0"] * 10.0 ** (p["n"] * cond.log_fo2_Pa) * np.exp(-p["Q"] * 1e3 / (R_GAS * cond.T_K))


_add(DiffusionCoefficient(
    key="opx_FeMg_sato2022",
    mineral="opx", species="Fe-Mg", kind="interdiffusion", transported_variable="XFe = Fe/(Fe+Mg)",
    label="Opx Fe-Mg, Sato et al. (2022) form of Dohmen et al. (2016): no composition term, a = b = c/3.5",
    citation="sato2022", equation_number="p. 11",
    equation_text=("D(001)Fe-Mg [m2/s] = 1.12e-6 fO2^(0.053 +/- 0.027) exp(-(308 +/- 23)/RT), fO2 in Pa, "
                   "D along c 3.5 times D along a and b"),
    func=_dohmen2016_no_composition,
    params={
        "D0": Parameter("D0", 1.12e-6, 0.0, "m2/s", "1s", "pre-exponential factor"),
        "n": Parameter("n", 0.053, 0.027, "-", "unstated", "fO2 exponent"),
        "Q": Parameter("Q", 308.0, 23.0, "kJ/mol", "unstated", "activation energy"),
    },
    axis_factors=_OPX_ANISO35, reference_axis="c", needs_fo2=True, fo2_unit="Pa",
    T_range=Range(1143.15, 1373.15, "K (870-1100 C, the Dohmen et al. 2016 Fs9 experiments)"),
    verified=True,
    verified_from=("Sato et al. (2022) JVGR 432:107686, p. 11 (law and the 3.5 factor) and p. 8 (profiles parallel "
                   "to a and/or b); Ban et al. (2026) EPS 78:37, eq. 1 and the anisotropy sentence, p. 17; rendered "
                   "pages, 9 October 2026"),
    secondary_citations=("dohmen2016", "ban2026"),
    variant_of="opx_FeMg_dohmen2016",
    uncertainty_note=("The +/- on n and Q are those printed (confidence level not stated); the study's minimum and "
                      "maximum times vary T by +/-28 C and fO2 by +/-0.5 log units only, so the law is held fixed."),
    notes=(STUDY_NOTE + " The Dohmen et al. (2016) Fs9 law as Sato et al. (2022, Zao) and Ban et al. (2026, "
           "Azumayama) print it: without the composition correction 10^(X_Fe - 0.09) that Dohmen et al. "
           "recommend, and with diffusion along c 3.5 times faster than along both a and b. Diffusor's Dohmen "
           "entry keeps the composition term and sets b = c (Dohmen et al. give Dc >= Db > Da). Both studies "
           "measured profiles parallel to a and/or b, at 1025 C and NNO (Sato) and at the temperatures of "
           "their Table 3 (Ban); Ban et al. do not state the unit of fO2. Mangler et al. (2022) use the same "
           "law along a without the composition term."),
))


@_with_latex(lambda c: (rf"\log_{{10}} D_c\ [\mathrm{{cm^2\,s^{{-1}}}}] = {_n(_v(c, 'c0') + 4)} + 2.6\,X_{{Fe}} - "
                        r"\frac{12530}{T} + \frac{1}{6}\log_{10}\frac{f_{O_2}}{f_{O_2}^{IW}}"))
def _krimer_costa_eq1(dc, cond, p):
    from ..thermo.buffers import log_fo2_buffer
    XFe = np.asarray(cond.x("XFe"), dtype=float)
    log_iw = log_fo2_buffer("IW", cond.T_K, cond.P_Pa)
    return 10.0 ** (p["c0"] + p["a"] * XFe - p["b"] / cond.T_K + p["n_fo2"] * (cond.log_fo2_bar - log_iw))


_add(DiffusionCoefficient(
    key="opx_FeMg_krimer_costa2017_eq1",
    mineral="opx", species="Fe-Mg", kind="interdiffusion", transported_variable="XFe = Fe/(Fe+Mg)",
    label="Opx Fe-Mg // c, Krimer & Costa (2017) eq. 1 (Allan et al. 2013 form with -5.44)",
    citation="krimer_costa2017", equation_number="1, p. 5",
    equation_text=("D_c = exp(-5.44 + 2.6 (1 - Mg#) - 12530/T) (fO2(NNO,T) / fO2(IW,T))^(1/6), 'D is in m2/s' "
                   "(as printed); Diffusor reads it as log10 D [cm2/s] = -5.44 + 2.6 X_Fe - 12530/T + 1/6 log10(fO2/fO2_IW)"),
    func=_krimer_costa_eq1,
    params={
        "c0": Parameter("c0", -9.44, 0.0, "log10(m2/s)", "1s", "-5.44 read in cm2/s"),
        "a": Parameter("a", 2.6, 0.0, "-", "1s", "X_Fe coefficient"),
        "b": Parameter("b", 12530.0, 0.0, "K", "1s", "temperature term"),
        "n_fo2": Parameter("n_fo2", 1.0 / 6.0, 0.0, "-", "1s", "fO2 exponent relative to IW"),
    },
    axis_factors={"c": 1.0, "b": 1.0 / 1.7, "a": 1.0 / 3.5}, reference_axis="c", requires=("XFe",),
    needs_fo2=True, fo2_unit="bar",
    T_range=Range(None, None, "K (the study's models are at 1000 C)"),
    verified=True,
    verified_from="Krimer & Costa (2017) GCA 196:271-288, eqs 1-2 and the anisotropy paragraph, p. 5, rendered page, 9 October 2026",
    secondary_citations=("allan2013", "ganguly_tazzoli1994"),
    variant_of="opx_FeMg_ganguly_tazzoli1994",
    uncertainty_note=_FIXED,
    notes=(STUDY_NOTE + " The study writes the Allan et al. (2013) formulation 'corrected for self-consistent "
           "units' as exp(...) with D in m2/s and the constant -5.44. Read that way it gives about 3e-7 m2/s at "
           "1000 C, and read as log10 m2/s about 1e-15 m2/s; the study says eq. 1 overlaps eq. 2 (the Dohmen "
           "et al. 2016 law, about 1e-19 m2/s) at its conditions, which only the log10 cm2/s reading of -5.44 "
           "gives. Diffusor reads it so, which is Diffusor's reading. The anisotropy ratios a : b : c = "
           "1/3.5 : 1/1.7 : 1 are those of the study's anisotropy test (p. 5). The study quotes fO2 at 1000 C "
           "as NNO = 1e-10 Pa and IW = 1e-15 Pa after Frost (1991); Diffusor computes IW from its own buffer."),
))

_add(DiffusionCoefficient(
    key="opx_FeMg_krimer_costa2017_eq2",
    mineral="opx", species="Fe-Mg", kind="interdiffusion", transported_variable="XFe = Fe/(Fe+Mg)",
    label="Opx Fe-Mg // c, Krimer & Costa (2017) eq. 2 (Dohmen et al. 2016 at Mg# 0.91, fO2^0.05)",
    citation="krimer_costa2017", equation_number="2, p. 5",
    equation_text="D_c = 1.12e-6 exp(-308/RT) (fO2)^0.05; R, the unit of 308 and of fO2 not printed (Dohmen et al. 2016: kJ/mol, Pa)",
    func=_dohmen2016_no_composition,
    params={
        "D0": Parameter("D0", 1.12e-6, 0.0, "m2/s", "1s", "pre-exponential factor"),
        "n": Parameter("n", 0.05, 0.0, "-", "1s", "fO2 exponent as printed"),
        "Q": Parameter("Q", 308.0, 0.0, "kJ/mol", "1s", "activation energy"),
    },
    axis_factors={"c": 1.0, "b": 1.0 / 1.7, "a": 1.0 / 3.5}, reference_axis="c", needs_fo2=True, fo2_unit="Pa",
    T_range=Range(1143.15, 1373.15, "K (870-1100 C, the Dohmen et al. 2016 Fs9 experiments)"),
    verified=True,
    verified_from="Krimer & Costa (2017) GCA 196:271-288, eq. 2 and the anisotropy paragraph, p. 5, rendered page, 9 October 2026",
    secondary_citations=("dohmen2016",),
    variant_of="opx_FeMg_dohmen2016",
    uncertainty_note=_FIXED,
    notes=(STUDY_NOTE + " Dohmen et al. (2016) for Mg# 0.91 as the study prints it: fO2 exponent 0.05 (Dohmen "
           "et al.: 0.053) and no composition term. The units of the activation energy and of fO2 are not "
           "printed; Diffusor uses those of Dohmen et al., kJ/mol and Pa. Anisotropy as in the study's test, "
           "a : b : c = 1/3.5 : 1/1.7 : 1."),
))


# ===========================================================================
# Clinopyroxene
# ===========================================================================
@_with_latex(lambda c: (rf"D_c = {_n(_v(c, 'D0'))}\,\exp\left(-\frac{{{_n(_v(c, 'Q'))}\ \mathrm{{kJ\,mol^{{-1}}}}}}{{RT}}\right)"
                        r"\quad D_a = D_b = D_c/3.5"))
def _plain(dc, cond, p):
    return p["D0"] * np.exp(-p["Q"] * 1e3 / (R_GAS * cond.T_K))


_add(DiffusionCoefficient(
    key="cpx_FeMg_ban2026",
    mineral="cpx", species="Fe-Mg", kind="interdiffusion", transported_variable="XFe = Fe/(Fe+Mg)",
    label="Cpx Fe-Mg, Ban et al. (2026) eq. 2: Mueller et al. (2013) with a = b = c/3.5",
    citation="ban2026", equation_number="2, p. 17",
    equation_text="D_cpx [m2/s] = (2.77 +/- 4.27) x 10^-7 exp(-320.7 +/- 16.0 kJ/mol / RT), along c; a and b 3.5 times slower",
    func=_plain,
    params={
        "D0": Parameter("D0", 2.77e-7, 4.27e-7, "m2/s", "unstated", "pre-exponential factor"),
        "Q": Parameter("Q", 320.7, 16.0, "kJ/mol", "unstated", "activation energy"),
    },
    axis_factors=_OPX_ANISO35, reference_axis="c",
    T_range=Range(1073.15, 1373.15, "K (800-1100 C, the Mueller et al. 2013 experiments)"),
    verified=True,
    verified_from="Ban et al. (2026) Earth Planets Space 78:37, eqs 1-2 and the anisotropy sentence, p. 17, rendered page, 9 October 2026",
    secondary_citations=("muller2013",),
    variant_of="cpx_FeMg_muller2013",
    uncertainty_note=_FIXED,
    notes=(STUDY_NOTE + " The Mueller et al. (2013) law along c, with the study's statement that diffusion "
           "along c is faster than along a and b by a factor of 3.5, given for both pyroxenes. Diffusor's "
           "Mueller entry has no anisotropy. The study measured profiles parallel to a and/or b."),
))


@_with_latex(lambda c: (rf"D = {_n(_v(c, 'D0'))}\,\left(\frac{{p_{{O_2}}}}{{0.21}}\right)^{{{_n(_v(c, 'n'))}}}"
                        rf"\,\exp\left(-\frac{{{_n(_v(c, 'Q'))}\ \mathrm{{kJ\,mol^{{-1}}}}}}{{RT}}\right)"
                        r"\quad (p_{O_2}\ \mathrm{in\ atm})"))
def _dimanov_wiedenbeck(dc, cond, p):
    log_atm = cond.log_fo2_bar - np.log10(1.01325)
    return p["D0"] * 10.0 ** (p["n"] * (log_atm - np.log10(0.21))) * np.exp(-p["Q"] * 1e3 / (R_GAS * cond.T_K))


_add(DiffusionCoefficient(
    key="cpx_FeMg_dimanov_wiedenbeck2006",
    mineral="cpx", species="Fe-Mg", kind="interdiffusion", transported_variable="(Fe,Mn)/(Fe+Mn+Mg)",
    label="Cpx (Fe,Mn)-Mg, Dimanov & Wiedenbeck (2006), with pO2 dependence",
    citation="dimanov_wiedenbeck2006", equation_number="concluding remarks, p. 716; n from p. 712",
    equation_text=("D (cm2/s) = 7.68e-2 (pO2 / pO2_ref)^(-0.22 +/- 0.02) exp(-(297 +/- 31 kJ/mol)/RT), "
                   "pO2_ref = 0.21 atm, as printed on p. 716; p. 712: rates increase with pO2, n = 0.22 +/- 0.02"),
    func=_dimanov_wiedenbeck,
    params={
        "D0": Parameter("D0", 7.68e-6, 0.0, "m2/s", "1s", "7.68e-2 cm2/s"),
        "n": Parameter("n", 0.22, 0.02, "-", "unstated", "pO2 exponent (sign of p. 712)"),
        "Q": Parameter("Q", 297.0, 31.0, "kJ/mol", "unstated", "activation energy"),
    },
    needs_fo2=True, fo2_unit="atm",
    T_range=Range(1323.15, 1423.15, "K (1050-1150 C, concluding remarks)"),
    fo2_range=Range(-15.0, -7.0, "log10 atm (1e-7 to 1e-15 atm)"),
    verified=True,
    verified_from=("Dimanov & Wiedenbeck (2006) Eur. J. Mineral. 18:705-718, p. 712 (oxygen fugacity exponent) and "
                   "p. 716 (concluding remarks), rendered pages, 9 October 2026"),
    secondary_citations=("kilgour2014",),
    uncertainty_note=_FIXED,
    notes=("A calibration on natural diopside by thin-film exchange, SIMS profiles, normal to the film. The "
           "concluding remarks print the pO2 term as (pO2/pO2_ref)^(-0.22); p. 712 says (Fe,Mn)-Mg "
           "interdiffusion rates are enhanced at higher oxygen partial pressures and gives n = 0.22 for the "
           "EPM-free samples, so Diffusor uses +0.22. The concluding expression gives log D = -3.82 (cm2/s) at "
           "1e-13 atm before the exp term, where p. 712 gives log D0 = -4.47 for the Arrhenius fit at 1e-13 atm; "
           "that fit is the separate entry cpx_FeMg_dimanov_wiedenbeck2006_1e13. Kilgour et al. (2014, Ruapehu) "
           "use this paper at NNO; their Table 1 values are not reproduced by either form with Diffusor's NNO "
           "(ratios 0.1 to 2.8)."),
))

_add(DiffusionCoefficient(
    key="cpx_FeMg_dimanov_wiedenbeck2006_1e13",
    mineral="cpx", species="Fe-Mg", kind="interdiffusion", transported_variable="(Fe,Mn)/(Fe+Mn+Mg)",
    label="Cpx (Fe,Mn)-Mg, Dimanov & Wiedenbeck (2006), Arrhenius fit at pO2 = 1e-13 atm",
    citation="dimanov_wiedenbeck2006", equation_number="p. 712 (activation energy)",
    equation_text="log D0 (cm2/s) = -4.47 +/- 1.20 and E = 297 +/- 31 kJ/mol at pO2 = 1e-13 atm, samples without EPM",
    func=_plain,
    params={
        "D0": Parameter("D0", 10.0 ** -8.47, 0.0, "m2/s", "1s", "log D0 = -4.47 in cm2/s"),
        "Q": Parameter("Q", 297.0, 31.0, "kJ/mol", "unstated", "activation energy"),
    },
    T_range=Range(1323.15, 1423.15, "K (1050-1150 C)"),
    fo2_range=Range(-13.0, -13.0, "log10 atm (fit made at 1e-13 atm)"),
    verified=True,
    verified_from="Dimanov & Wiedenbeck (2006) Eur. J. Mineral. 18:705-718, p. 712, rendered page, 9 October 2026",
    uncertainty_note=_FIXED,
    notes=("The Arrhenius fit of the EPM-free samples at a fixed pO2 of 1e-13 atm (about IW at 1100 C). The fit "
           "to all samples gives E = 345 +/- 26 kJ/mol and log D0 = -2.57 +/- 0.99 (cm2/s). See "
           "cpx_FeMg_dimanov_wiedenbeck2006 for the pO2-dependent form."),
))


# ===========================================================================
# Plagioclase
# ===========================================================================
@_with_latex(lambda c: (rf"D_{{Sr}} = 2.92\times10^{{-4.1X_{{An}} - 4.08}}\,\exp\left(-\frac{{276000}}{{RT}}\right)"))
def _druitt_sr(dc, cond, p):
    XAn = np.asarray(cond.x("XAn"), dtype=float)
    return p["k"] * 10.0 ** (p["a"] * XAn + p["c0"]) * np.exp(-p["Q"] * 1e3 / (R_GAS * cond.T_K))


_add(DiffusionCoefficient(
    key="plag_Sr_druitt2012",
    mineral="plagioclase", species="Sr", kind="tracer", transported_variable="Sr concentration",
    label="Plagioclase Sr, Druitt et al. (2012) eq. 2 (Giletti & Casserly 1994 with the factor 2.92)",
    citation="druitt2012", equation_number="2, Methods p. 5",
    equation_text="D_Sr = 2.92 x 10^(-4.1 X_An - 4.08) exp(-276,000/RT), D in m2/s",
    func=_druitt_sr,
    params={
        "k": Parameter("k", 2.92, 0.0, "-", "1s", "factor as printed"),
        "a": Parameter("a", -4.1, 0.0, "-", "1s", "X_An coefficient"),
        "c0": Parameter("c0", -4.08, 0.0, "log10(m2/s)", "1s", "intercept"),
        "Q": Parameter("Q", 276.0, 0.0, "kJ/mol", "1s", "activation energy"),
    },
    requires=("XAn",),
    T_range=Range(None, None, "K (the study applies it at 855-900 C)"),
    X_range=Range(0.0, 1.0, "X_An"),
    verified=True,
    verified_from="Druitt et al. (2012) Nature 482:77-80, Methods eqs 1-2, p. 5, rendered page, 9 October 2026",
    secondary_citations=("giletti_casserly1994",),
    variant_of="plag_Sr_giletti_casserly1994",
    uncertainty_note=_FIXED,
    notes=(STUDY_NOTE + " Eq. 2 of the study is the Giletti & Casserly (1994) Sr law multiplied by 2.92, the "
           "factor of its Mg law (eq. 1, Costa et al. 2003, identical to Diffusor's plag_Mg_costa2003). The "
           "factor makes D 0.47 log units higher than Diffusor's Giletti & Casserly entry. The study states "
           "that by eqs 1 and 2 D_Mg > D_Sr at all temperatures for a given An content."),
))


def _li_entry(key, phase, logD0, s_logD0, Q, s_Q, T_lo, T_hi):
    return DiffusionCoefficient(
        key=key, mineral="plagioclase", species="Li", kind="tracer", transported_variable="Li (6Li tracer)",
        label=f"Plagioclase Li in {phase}, Giletti & Shanahan (1997) Table 3",
        citation="giletti_shanahan1997", equation_number="Table 3, p. 13",
        equation_text=f"log D0 (m2/s) = {logD0} +/- {s_logD0}, Q = {Q:g} +/- {s_Q:g} kJ/mol ({phase}, +/-2 sigma), normal to (010)",
        func=_plain,
        params={
            "D0": Parameter("D0", 10.0 ** logD0, 0.0, "m2/s", "2s", f"log D0 = {logD0} +/- {s_logD0} (2 sigma)"),
            "Q": Parameter("Q", Q, s_Q, "kJ/mol", "2s", "activation energy"),
        },
        T_range=Range(T_lo + 273.15, T_hi + 273.15, f"K ({T_lo}-{T_hi} C)"),
        verified=True,
        verified_from="Giletti & Shanahan (1997) Chem. Geol. 139:3-20, Table 3, p. 13, rendered page, 9 October 2026",
        secondary_citations=("ellis2018",),
        uncertainty_note=("Table 3 gives +/-2 sigma on log D0 and Q without their covariance; the law is held "
                          "fixed in the Monte Carlo."),
        notes=("Tracer diffusion of 6Li at 1 atm, dry, normal to (010). Ellis et al. (2018) take Li diffusion in "
               "plagioclase from this paper; they do not state the composition, so both the albite and the "
               "anorthite fits are entries. The paper finds Li diffusion in anorthite essentially the same as in "
               "albite."),
    )


_add(_li_entry("plag_Li_giletti_shanahan1997_ab", "albite", -3.8, 1.0, 146.0, 14.0, 200, 800))
_add(_li_entry("plag_Li_giletti_shanahan1997_an", "anorthite", -3.6, 3.6, 151.0, 62.0, 500, 850))


# ===========================================================================
# Titanomagnetite
# ===========================================================================
@_with_latex(lambda c: r"\log_{10} D_{Fe-Ti} = -7.14 - \frac{1.17\times10^4}{T} + 3.21\,X_{Usp},\ X_{Usp} = 3\,x_{Ti}")
def _araya_feti(dc, cond, p):
    X = 3.0 * np.asarray(cond.x("xTi"), dtype=float)
    return 10.0 ** (p["c0"] - p["b"] / cond.T_K + p["a"] * X)


def _araya_fe_branch(p, T, log_fo2_atm, which):
    v = p[f"DV_{which}"] * np.exp(p[f"QV_{which}"] * 1e3 / (R_GAS * T)) * p["Cv"] * 10.0 ** (2.0 / 3.0 * log_fo2_atm)
    i = p[f"DI_{which}"] * np.exp(-p[f"QI_{which}"] * 1e3 / (R_GAS * T)) * 10.0 ** (-2.0 / 3.0 * log_fo2_atm)
    return v + i


@_with_latex(lambda c: (r"\log_{10} D_{Fe} = \frac{(X_{Usp} - 0.3)\log_{10} D^{0.6}_{Fe} - (X_{Usp} - 0.6)"
                        r"\log_{10} D^{0.3}_{Fe}}{0.3};\ D^{X}_{Fe} = D_V e^{Q_V/RT} C_v f_{O_2}^{2/3} + "
                        r"D_I e^{-Q_I/RT} f_{O_2}^{-2/3}"))
def _araya_fe(dc, cond, p):
    X = 3.0 * np.asarray(cond.x("xTi"), dtype=float)
    lf = cond.log_fo2_bar - np.log10(1.01325)
    d3 = np.log10(_araya_fe_branch(p, cond.T_K, lf, "03"))
    d6 = np.log10(_araya_fe_branch(p, cond.T_K, lf, "06"))
    return 10.0 ** (((X - 0.3) * d6 - (X - 0.6) * d3) / 0.3)


_ARAYA_COMMON = dict(
    citation="araya2024",
    verified=True,
    verified_from="Araya et al. (2024) JGR Solid Earth, doi:10.1029/2023JB028558, eqs 3-6, pp. 16-17, rendered pages, 9 October 2026",
)

_add(DiffusionCoefficient(
    key="mt_FeTi_araya2024",
    mineral="magnetite", species="Fe-Ti", kind="interdiffusion", transported_variable="X_Usp",
    label="Titanomagnetite Fe-Ti, Araya et al. (2024) eq. 3 (fit to the NNO runs of Aragon et al. 1984)",
    equation_number="3, p. 16",
    equation_text="log D_Fe-Ti = -7.14 - 1.17e4/T + 3.21 X_Usp, D in m2/s, T in K, X_Usp the ulvospinel mole fraction",
    func=_araya_feti,
    params={
        "c0": Parameter("c0", -7.14, 0.0, "log10(m2/s)", "1s", "intercept"),
        "b": Parameter("b", 1.17e4, 0.0, "K", "1s", "temperature term"),
        "a": Parameter("a", 3.21, 0.0, "-", "1s", "X_Usp coefficient"),
    },
    sigma_logD=0.25,
    sigma_logD_basis=("published: Araya et al. (2024, p. 16) assume D_Fe-Ti had an uncertainty of 0.25 log units, "
                      "the agreement of eq. 3 with the original data"),
    requires=("xTi",),
    T_range=Range(1263.15, 1493.15, "K (990-1220 C, the Aragon et al. 1984 anneals)"),
    X_range=Range(0.0, 0.26, "x_Ti (X_Usp = 3 x_Ti)"),
    secondary_citations=("aragon1984",),
    variant_of="mt_FeTi_aragon1984",
    notes=(STUDY_NOTE + " The study's own formulation of the NNO-buffered interdiffusion experiments of Aragon "
           "et al. (1984); it reproduces those data within 0.25 log units and differs from Freer & Hauptman "
           "(1978) by less than 0.6 log units at the study's conditions (p. 16). It contains no fO2 term. "
           "Diffusor's profile variable is x_Ti = Ti/(Fe+Ti); X_Usp = 3 x_Ti. The study models a sphere "
           "(Crank eq. 1.8) for Sakurajima magnetite."),
    **_ARAYA_COMMON,
))

_ARAYA_FE_PARAMS = {
    "DV_03": Parameter("DV_03", 1.68e-13, 0.0, "m2/s", "1s", "vacancy term, X_Usp = 0.3 (eq. 4)"),
    "QV_03": Parameter("QV_03", 124.7, 0.0, "kJ/mol", "1s", "exp(+QV/RT), eq. 4"),
    "DI_03": Parameter("DI_03", 4.81e5, 0.0, "m2/s", "1s", "interstitial term, X_Usp = 0.3"),
    "QI_03": Parameter("QI_03", 667.7, 0.0, "kJ/mol", "1s", "eq. 4"),
    "DV_06": Parameter("DV_06", 1.36e-13, 0.0, "m2/s", "1s", "vacancy term, X_Usp = 0.6 (eq. 5)"),
    "QV_06": Parameter("QV_06", 167.5, 0.0, "kJ/mol", "1s", "eq. 5"),
    "DI_06": Parameter("DI_06", 2.34e1, 0.0, "m2/s", "1s", "interstitial term, X_Usp = 0.6"),
    "QI_06": Parameter("QI_06", 562.4, 0.0, "kJ/mol", "1s", "eq. 5"),
    "Cv": Parameter("Cv", 1.0, 0.0, "-", "1s", "correction factor; value not printed, assumed 1"),
}
_ARAYA_FE_NOTE = ("Eqs 4 and 5 are the study's formulations of the Aggarwal & Dieckmann (2002a) Fe tracer data at "
                  "X_Usp = 0.3 and 0.6 ('similar to that presented in Table 12 of Van Orman & Crispin 2010'); "
                  "eq. 5 equals the x_Ti = 0.2 row of that table to the printed digits, eq. 4 is the study's own. "
                  "Eq. 6 interpolates log D linearly in X_Usp between them. The study does not print the value of "
                  "the correction factor C_v (Aggarwal & Dieckmann 2002a, 2002b), the gas constant, the unit of "
                  "the energies or of fO2. Diffusor uses C_v = 1 (an assumption), kJ/mol and fO2 in atm, as in Van "
                  "Orman & Crispin Table 12. X_Usp = 3 x_Ti.")

for _key, _species, _label, _kind, _extra in (
        ("mt_Fe_araya2024", "Fe", "Titanomagnetite Fe, Araya et al. (2024) eqs 4-6", "tracer", ""),
        ("mt_Mg_araya2024", "Mg", "Titanomagnetite Mg taken equal to Fe, Araya et al. (2024) eqs 4-6", "effective",
         " The study takes D_Mg equal to D_Fe (p. 16-17) and holds it constant for the Mg zoning.")):
    _add(DiffusionCoefficient(
        key=_key, mineral="magnetite", species=_species, kind=_kind, transported_variable=f"{_species} concentration",
        label=_label, equation_number="4-6, p. 17",
        equation_text=("D_Fe^(X_Usp=0.3) = 1.68e-13 exp(124.7/RT) Cv fO2^(2/3) + 4.81e5 exp(-667.7/RT) fO2^(-2/3); "
                       "D_Fe^(X_Usp=0.6) = 1.36e-13 exp(167.5/RT) Cv fO2^(2/3) + 2.34e1 exp(-562.4/RT) fO2^(-2/3); "
                       "log D_Fe = (1/0.3)[(X_Usp - 0.3) log D^0.6 - (X_Usp - 0.6) log D^0.3]"),
        func=_araya_fe, params=dict(_ARAYA_FE_PARAMS),
        requires=("xTi",), needs_fo2=True, fo2_unit="atm",
        T_range=Range(1173.15, 1673.15, "K (900-1400 C, the Aggarwal & Dieckmann experiments; Diffusor's range)"),
        X_range=Range(0.1, 0.2, "x_Ti (X_Usp 0.3 to 0.6, the interpolation interval)"),
        secondary_citations=("aggarwal_dieckmann2002", "vanorman_crispin2010"),
        variant_of="mt_Fe_vanorman_crispin2010",
        uncertainty_note=_FIXED,
        notes=STUDY_NOTE + " " + _ARAYA_FE_NOTE + _extra,
        **_ARAYA_COMMON,
    ))


@_with_latex(lambda c: (r"D^*_{Ti} = D^0_0 e^{-E_0/RT} + (1 - X)\left(D^0_V e^{-E_V/RT}\left[\frac{1+X}{1-X}\right]^{8/3}"
                        r"f_{O_2}^{2/3} - D^0_I e^{-E_I/RT}\left[\frac{1-X}{1+X}\right]^{8/3}f_{O_2}^{-2/3}\right)"))
def _morgado2019(dc, cond, p):
    X = 3.0 * np.asarray(cond.x("xTi"), dtype=float)
    T = cond.T_K
    f = 10.0 ** (cond.log_fo2_bar - np.log10(1.01325))
    return (p["D0_0"] * np.exp(-p["E0"] * 1e3 / (R_GAS * T))
            + (1.0 - X) * (p["D0_V"] * np.exp(-p["EV"] * 1e3 / (R_GAS * T)) * ((1 + X) / (1 - X)) ** (8 / 3) * f ** (2 / 3)
                           - p["D0_I"] * np.exp(-p["EI"] * 1e3 / (R_GAS * T)) * ((1 - X) / (1 + X)) ** (8 / 3) * f ** (-2 / 3)))


_add(DiffusionCoefficient(
    key="mt_FeTi_morgado2019",
    mineral="magnetite", species="Fe-Ti", kind="interdiffusion", transported_variable="X_Ti",
    label="Titanomagnetite Fe-Ti, Morgado et al. (2019) eq. 5 (Aragon et al. 1984 point-defect form)",
    citation="morgado2019", equation_number="5, p. 4; constants SM2 'Errors short' K1:L7",
    equation_text=("D*_Ti = D0_0 exp(-E0/RT) + (1 - X) (D0_V exp(-E_V/RT) [(1+X)/(1-X)]^(8/3) fO2^(2/3) - "
                   "D0_I exp(-E_I/RT) [(1-X)/(1+X)]^(8/3) fO2^(-2/3)), fO2 in atm, X = X_Usp; D0_0 = 8.2e-7 m2/s, "
                   "E0 = 251.04 kJ/mol, D0_V = 1.75e-16, E_V = -117.15 kJ/mol, D0_I = 390, E_I = 665.256 kJ/mol"),
    func=_morgado2019,
    params={
        "D0_0": Parameter("D0_0", 8.2e-7, 0.0, "m2/s", "1s", "SM2 L1"),
        "E0": Parameter("E0", 251.04, 0.0, "kJ/mol", "1s", "SM2 L2 (60,000 cal/mol)"),
        "D0_V": Parameter("D0_V", 1.75e-16, 0.0, "m2/s atm^-2/3", "1s", "SM2 L4"),
        "EV": Parameter("EV", -117.15, 0.0, "kJ/mol", "1s", "SM2 L5"),
        "D0_I": Parameter("D0_I", 390.0, 0.0, "m2/s atm^2/3", "1s", "SM2 L7"),
        "EI": Parameter("EI", 665.256, 0.0, "kJ/mol", "1s", "SM2 L6"),
    },
    requires=("xTi",), needs_fo2=True, fo2_unit="atm",
    T_range=Range(1263.15, 1493.15, "K (990-1220 C, the Aragon et al. 1984 anneals)"),
    X_range=Range(0.0, 0.26, "x_Ti (X_Usp = 3 x_Ti)"),
    verified=True,
    verified_from=("Morgado et al. (2019) CMP 174:61, eqs 4-5 and their definitions, p. 4, rendered page; constants "
                   "from Supplementary Material 2, 'SM2 - Errors_Total_dependent.xlsx', sheet 'Errors short', "
                   "cells K1:L7; 9 October 2026"),
    secondary_citations=("aragon1984",),
    variant_of="mt_FeTi_aragon1984",
    uncertainty_note=("The study propagates T, fO2 and X_Usp uncertainties and adds 0.1 log units (Aragon et al. "
                      "1984); the law itself is held fixed here."),
    notes=(STUDY_NOTE + " The paper prints the form and the units (energies in cal/mol, R in cal/(mol K), fO2 "
           "in atm) and gives the constants in Supplementary Material 2, whose workbook holds them in kJ/mol "
           "(60,000 cal/mol = 251.04 kJ/mol). Eq. 5 subtracts the interstitial term; the workbook formulas add "
           "it, which changes D by about 0.1 % at the study's conditions. This entry follows the printed "
           "eq. 5. X is the ulvospinel fraction, X = 3 x_Ti."),
))


# --- one element's law used for another ---------------------------------------------------------------
def _alias(base_key, key, species, citation, label, note):
    base = next(c for c in _magnetite.COEFFICIENTS if c.key == base_key)

    def func(dc, cond, p):
        return base.func(dc, cond, p)

    def latex(c):
        from .latex import coefficient_latex
        return coefficient_latex(base)     # typeset the row of the element whose law is used

    func.latex = latex
    from .uncertainty_basis import SIGMA_LOGD_BASIS
    return replace(base, key=key, species=species, citation=citation, label=label, kind="effective", func=func,
                   sigma_logD_basis=SIGMA_LOGD_BASIS.get(base_key, base.sigma_logD_basis),
                   recommended=False, superseded_by="", superseded_note="", variant_of=base_key,
                   secondary_citations=(base.citation,),
                   verified_from=(base.verified_from + " The assignment to " + species + " is stated by "
                                  "Tomiya et al. (2013) Bull. Volcanol. 75:750, pp. 14-15 of 19, 'Diffusion "
                                  "coefficients in magnetite' (prose, read 9 October 2026)."),
                   notes=STUDY_NOTE + " " + note)


_add(_alias("mt_Ti_vanorman_crispin2010", "mt_Al_tomiya2013", "Al", "tomiya2013",
            "Titanomagnetite Al taken equal to Ti, Tomiya et al. (2013) (Van Orman & Crispin 2010 Table 12 Ti)",
            "Tomiya et al. (2013, Shinmoedake 2011) modelled Al zoning with the Ti law of Van Orman & Crispin "
            "(2010) Table 12. Diffusor's own Al entries give D_Al about 1.6 log units lower at their conditions."))
_add(_alias("mt_Fe_vanorman_crispin2010", "mt_Mg_tomiya2013", "Mg", "tomiya2013",
            "Titanomagnetite Mg taken equal to Fe, Tomiya et al. (2013) (Van Orman & Crispin 2010 Table 12 Fe)",
            "Tomiya et al. (2013, Shinmoedake 2011) took D_Mg equal to the Fe tracer coefficient of Van Orman & "
            "Crispin (2010) Table 12."))
