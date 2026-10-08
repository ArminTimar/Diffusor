"""Diffusion-driven isotope fractionation: beta factors, abundances and delta values.

Two isotopes e and f of one element diffuse at rates related by

    D_e / D_f = (M_f / M_e)^beta

(Richter et al. 1999, 2003, as written by Oeser et al. 2026, eq. 1). Mass
numbers are used, as in the papers' own formulas (Oeser et al. 2026, p. 51;
Richter et al. 2014, abstract). beta is a property of the element, the host,
the direction and of the diffusion model it was fitted with: the same Li
profiles need beta = 0.27 with a two-site model and 0.44 with a one-site model
(Richter et al. 2014, p. 363), and treating Fe and Mg isotopes as independent
binaries instead of one coupled system changes beta_Mg strongly (Oeser et al.
2026, Fig. 11). Each entry therefore names the model it belongs to.

delta values are per mil deviations of a heavy/light ratio R = N(heavy)/N(light)
from a reference ratio (the convention of IUPAC, Prohaska et al. 2022, eqs. 3
and 4; the isotope pairs 7Li/6Li, 26Mg/24Mg and 56Fe/54Fe are Diffusor's choice).
The reference ratio is NOT that of a delta reference material (such as L-SVEC,
DSM3 or IRMM-014): Diffusor splits an element into isotopes with the
representative abundances of Berglund & Wieser (2011), column 9 of their
Table 1, taken as the reference composition. Input and modelled delta values
are therefore on the same internal scale, which is offset from the published
scale of any reference material by an amount that is not quantified here.
Berglund & Wieser state that the representative values may match no real
sample and recommend samples of precisely known composition for precise work.
Only a modelled delta difference (not an absolute delta value) should be compared
with published data. A test raises the 6Li abundance by 1 % of its value (which lowers
the 7Li/6Li reference ratio by about 10.7 per mil) and the modelled Li profile changes by under 1 % of its
range; no such test was made for Mg and Fe.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional, Tuple

import numpy as np

from ..references import cite

# Berglund & Wieser (2011), Table 1 column 9, representative isotopic compositions (mole fractions).
# These are not delta reference materials. The authors print the Li values in square brackets because
# commercial laboratory reagents depleted in 6Li exist (their note c); Diffusor does not use those.
ABUNDANCES: Dict[str, Dict[int, float]] = {
    "Li": {6: 0.0759, 7: 0.9241},
    "Mg": {24: 0.7899, 25: 0.1000, 26: 0.1101},
    "Fe": {54: 0.05845, 56: 0.91754, 57: 0.02119, 58: 0.00282},
}
# the ratio each delta value refers to: (heavy, light)
DELTA_RATIO: Dict[str, Tuple[int, int]] = {"Li": (7, 6), "Mg": (26, 24), "Fe": (56, 54)}
DELTA_LABEL = {"Li": "δ7Li", "Mg": "δ26Mg", "Fe": "δ56Fe"}


@dataclass(frozen=True)
class IsotopeBeta:
    key: str
    element: str
    mineral: str
    beta: float
    sigma: float                 # one standard deviation of the published values (0 if none given)
    direction: str               # 'a', 'b', 'c' or a description
    citation: str
    where: str
    diffusion_model: str         # the model the beta belongs to
    conditions: str
    upper_bound: bool = False    # the published value is a maximum
    assumed: bool = False        # an assumption in the source, not a measurement
    notes: str = ""

    def describe(self) -> str:
        val = f"beta {'<= ' if self.upper_bound else '= '}{self.beta:g}"
        if self.sigma:
            val += f" +/- {self.sigma:g}"
        return (f"{self.element} in {self.mineral} // {self.direction}: {val} "
                f"({cite(self.citation)}, {self.where}; {self.diffusion_model}; {self.conditions})"
                + (" -- assumed in the source, not measured" if self.assumed else "")
                + (f". {self.notes}" if self.notes else ""))


_OL = "San Carlos olivine Fo91.5, 1100-1300 C, log fO2 -7.1 to -2.3 Pa"
_MC = "coupled seven-isotope Fe-Mg multicomponent model"
BETAS: Dict[str, IsotopeBeta] = {b.key: b for b in [
    IsotopeBeta("ol_Fe_oeser2026_a", "Fe", "olivine", 0.19, 0.07, "a", "oeser2026", "Table 5", _MC, _OL),
    IsotopeBeta("ol_Fe_oeser2026_b", "Fe", "olivine", 0.22, 0.06, "b", "oeser2026", "Table 5", _MC, _OL),
    IsotopeBeta("ol_Fe_oeser2026_c", "Fe", "olivine", 0.07, 0.02, "c", "oeser2026", "Table 5", _MC, _OL,
                upper_bound=True, notes="Fractionation parallel to c is hardly resolvable; the value is a maximum."),
    IsotopeBeta("ol_Mg_oeser2026_a", "Mg", "olivine", 0.18, 0.06, "a", "oeser2026", "Table 5", _MC, _OL),
    IsotopeBeta("ol_Mg_oeser2026_b", "Mg", "olivine", 0.17, 0.05, "b", "oeser2026", "Table 5", _MC, _OL),
    IsotopeBeta("ol_Mg_oeser2026_c", "Mg", "olivine", 0.10, 0.07, "c", "oeser2026", "Table 5", _MC, _OL,
                upper_bound=True, notes="Fractionation parallel to c is hardly resolvable; the value is a maximum."),
    IsotopeBeta("ol_Fe_teng2011", "Fe", "olivine", 0.05, 0.0, "unspecified", "teng2011", "p. 318-321",
                "binary Mg-Fe interdiffusion (one coefficient, Dohmen & Chakraborty 2007), beta_Mg = beta_Fe",
                "Hawaiian olivine (Kilauea Iki)", assumed=True,
                notes="Teng et al. take beta ~ 0.05 for Mg and Fe from silicate melts (Richter et al. 2009) because "
                      "it is not known for olivine; the crystal orientation was not taken into account. Used to "
                      "compare natural trends with theory; kept to reproduce that comparison."),
    IsotopeBeta("ol_Mg_teng2011", "Mg", "olivine", 0.05, 0.0, "unspecified", "teng2011", "p. 318-321",
                "binary Mg-Fe interdiffusion (one coefficient, Dohmen & Chakraborty 2007), beta_Mg = beta_Fe",
                "Hawaiian olivine (Kilauea Iki)", assumed=True,
                notes="Teng et al. take beta ~ 0.05 for Mg and Fe from silicate melts (Richter et al. 2009) because "
                      "it is not known for olivine; the crystal orientation was not taken into account. Used to "
                      "compare natural trends with theory; kept to reproduce that comparison."),
    IsotopeBeta("ol_Li_richter2017", "Li", "olivine", 0.4, 0.1, "a and c", "richter2017li", "abstract",
                "two-site Li model (interstitial and metal site)", "San Carlos olivine, 1 atm and piston cylinder",
                notes="Fractionation along b appears somewhat lower (abstract). A different model needs a "
                      "different beta."),
    IsotopeBeta("plag_Li_pohl2024_interstitial", "Li", "plagioclase", 0.49, 0.02, "unspecified", "pohl2024",
                "p. 997, section 4.7", "multispecies model, interstitial Li",
                "An61 labradorite, 606-1114 C, 200 MPa", notes="Pair with plag_Li_pohl2024_interstitial. The experiments went into a crystal surface perpendicular to [001]; the paper does not discuss an orientation dependence of beta (it only cites earlier Li studies that found none for the diffusion coefficient), so the direction is left unspecified."),
    IsotopeBeta("plag_Li_pohl2024_vacancy", "Li", "plagioclase", 0.24, 0.02, "unspecified", "pohl2024",
                "p. 997, section 4.7", "multispecies model, Li on A sites (vacancy mechanism)",
                "An61 labradorite, 606-1114 C, 200 MPa", notes="Pair with plag_Li_pohl2024_vacancy. The experiments went into a crystal surface perpendicular to [001]; the paper does not discuss an orientation dependence of beta (it only cites earlier Li studies that found none for the diffusion coefficient), so the direction is left unspecified."),
    IsotopeBeta("cpx_Li_richter2014", "Li", "cpx", 0.27, 0.0, "unspecified", "richter2014li", "abstract and p. 363",
                "two-site Li model (interstitial and metal site)", "Templeton augite, 900 C",
                notes="A single-species model needs beta = 0.44 and does not fit the shape of the isotope "
                      "profiles (p. 363). Natural pyroxenes were fitted with 0.25-0.30 (abstract)."),
]}


def betas_for(element: str, mineral: Optional[str] = None) -> List[IsotopeBeta]:
    return [b for b in BETAS.values() if b.element == element and (mineral is None or b.mineral == mineral)]


def mass_factor(m_ref: int, m: int, beta: float) -> float:
    """D(m) / D(m_ref) = (m_ref / m)^beta."""
    return (m_ref / m) ** beta


def _ratio_exponent(element: str, mass: int) -> float:
    heavy, light = DELTA_RATIO[element]
    return np.log(mass / light) / np.log(heavy / light)


def split_isotopes(element: str, total, delta) -> Dict[int, np.ndarray]:
    """Isotope concentrations of an element of concentration ``total`` whose
    principal delta value (per mil, see :data:`DELTA_RATIO`) is ``delta``.

    Isotopes other than the two in the ratio follow the exponential mass law,
    R_i/R_i,ref = (R/R_ref)^(ln(m_i/m_light)/ln(m_heavy/m_light)).
    """
    ab = ABUNDANCES[element]
    light = DELTA_RATIO[element][1]
    factor = 1.0 + np.asarray(delta, dtype=float) / 1000.0
    rel = {m: (ab[m] / ab[light]) * factor ** _ratio_exponent(element, m) for m in ab}
    s = sum(rel.values())
    total = np.asarray(total, dtype=float)
    return {m: total * r / s for m, r in rel.items()}


def delta_value(element: str, conc: Mapping[int, np.ndarray]) -> np.ndarray:
    """Principal delta value (per mil) from isotope concentrations."""
    heavy, light = DELTA_RATIO[element]
    ab = ABUNDANCES[element]
    with np.errstate(divide="ignore", invalid="ignore"):
        return (np.asarray(conc[heavy]) / np.asarray(conc[light]) / (ab[heavy] / ab[light]) - 1.0) * 1000.0
