"""Minerals available in the coefficient registry.

Sources for the structural/compositional conventions: Deer, Howie & Zussman,
*An Introduction to the Rock-Forming Minerals*, 3rd edition (2013), citekey
``deer2013``.  Each mineral's formula, crystal system and (where used) space
group is read from that mineral's data sheet in the 3rd edition; the printed
page is given in the mineral's ``notes`` (``DHZ 2013, p. N``).  The end-member
fractions come from Appendix 3, p. 488 (olivine, pyroxene) and p. 489
(feldspar), and the symbol Fe* from the list of abbreviations (p. x).  The
2nd edition is not used anywhere.  Xenotime has no data sheet in the 3rd
edition; its formula and the tetragonal system rest on the isostructural
relation to zircon stated by Cherniak (2006).  The diffusion statements in the
``notes`` fields come from the diffusion papers named in each note, not from
Deer, Howie & Zussman.  The diffusion papers themselves define the composition
variable used by each coefficient (see :mod:`diffusor.coefficients`).
"""
from __future__ import annotations

from .base import CompositionVariable, Mineral, Species

# --- composition variables ------------------------------------------------------
X_FE = CompositionVariable(
    "XFe", "X_Fe", "Fe/(Fe+Mg) molar, all Fe as FeO (Deer, Howie & Zussman 2013 use "
    "Fe* = Fe2+ + Fe3+, p. x. Binary X_Fe = Fe/(Fe+Mg) for orthopyroxene, p. 104)",
    default_mode="A/(A+B)", citation="deer2013")
# Appendix 3 (DHZ 2013, p. 488) gives the olivine end-member percentages as
# Fo = 100 Mg/(Mg + Fe*) with Fe* = Fe2+ + Fe3+ (Mn may be added to Fe*, giving a
# slightly lower Fo); p. 7 defines forsterite as Fo100-50 by the 50 per cent rule.
X_FO = CompositionVariable(
    "XFo", "X_Fo", "Mg/(Mg+Fe) molar (forsterite fraction, total Fe, DHZ 2013, "
    "Appendix 3, p. 488)",
    default_mode="B/(A+B)", citation="deer2013")
# Deer, Howie & Zussman (2013, Appendix 3, p. 489) give An = 100 Ca/(Ca+Na) for a
# plagioclase feldspar (binary) and An = 100 Ca/(Ca+Na+K) for a ternary feldspar.  p. 292
# says most plagioclases are ternary (An-Ab-Or), with Or usually at low concentration.
# The K-inclusive (ternary) form is used here, which is identical for K-free
# plagioclase and slightly lower for K-bearing plagioclase.
X_AN = CompositionVariable(
    "XAn", "X_An", "Ca/(Ca+Na+K) molar (anorthite fraction, ternary-feldspar form "
    "equal to Ca/(Ca+Na) without K, DHZ 2013, Appendix 3, p. 489)",
    default_mode="A/(A+B)", citation="deer2013")
# For ternary feldspar Appendix 3 (p. 489) gives Or = 100 K/(Ca+Na+K); for alkali
# feldspar 100 K/(Na+K).
X_OR = CompositionVariable(
    "XOr", "X_Or", "K/(K+Na+Ca) molar (orthoclase fraction, ternary-feldspar form "
    "equal to K/(Na+K) without Ca, DHZ 2013, Appendix 3, p. 489)",
    default_mode="A", citation="deer2013")
X_TI = CompositionVariable(
    "xTi", "x_Ti", "Ti per cation site in (Ti_x Fe_(1-x))3 O4, x_Ti = X_Usp/3",
    default_mode="A", citation="vanorman_crispin2010")

# --- minerals --------------------------------------------------------------------
OLIVINE = Mineral(
    key="olivine", name="Olivine", formula="(Mg,Fe)2SiO4", system="orthorhombic",
    axes=("a", "b", "c"),
    composition_variable=X_FO,
    species={
        "Fe-Mg": Species("Fe-Mg", "Fe-Mg interdiffusion", "interdiffusion",
                         "modelled as forsterite content X_Fo or X_Fe"),
        "Ni": Species("Ni", "Ni", "tracer"),
        "Mn": Species("Mn", "Mn", "tracer"),
        "Ca": Species("Ca", "Ca", "tracer"),
        "P": Species("P", "P", "trace"),
        "Be": Species("Be", "Be", "trace"),
    },
    notes=("Formula, orthorhombic system and the forsterite cell (a 4.75, b 10.20, c 5.98 A, "
           "space group Pbnm): DHZ 2013, p. 5. Fo is defined by the 50 per cent rule on p. 7. "
           "Each diffusion law states its own axis orientation. "
           "Diffusion is strongly anisotropic: D[001] is about 6x D[100] and D[010] "
           "(Dohmen & Chakraborty 2007). Hartley et al. (2016, p. 61) quote the same factor of "
           "about 6 from Nakamura & Schmalzried (1983) and Dohmen & Chakraborty (2007), and for "
           "some Laki crystals modelled profiles along two directions to confirm the expected "
           "anisotropy. They report no measured factor."),
)

ORTHOPYROXENE = Mineral(
    key="opx", name="Orthopyroxene", formula="(Mg,Fe)2Si2O6", system="orthorhombic",
    axes=("a", "b", "c"),
    composition_variable=X_FE,
    species={
        "Fe-Mg": Species("Fe-Mg", "Fe-Mg interdiffusion", "interdiffusion"),
        "Mg": Species("Mg", "Mg self-diffusion", "tracer"),
        "Lu": Species("Lu", "Lu", "trace", "rare earth element, Dias et al. (2025)"),
        "Ce": Species("Ce", "Ce", "trace", "rare earth element, Dias et al. (2025)"),
        "Eu": Species("Eu", "Eu", "trace", "rare earth element, assumed trivalent"),
    },
    notes=("Formula: DHZ 2013, p. 95. Orthorhombic, space group Pbca, a 18.22-18.43, "
           "b 8.81-9.08, c 5.17-5.24 A: DHZ 2013, p. 102. "
           "D//[001] >= D//[010] > D//[100] (diffusion papers, not DHZ). Dohmen et al. (2016) "
           "give D_a = D_c / 3.5."),
)

CLINOPYROXENE = Mineral(
    key="cpx", name="Clinopyroxene", formula="Ca(Mg,Fe)Si2O6", system="monoclinic",
    axes=("a", "b", "c"),
    composition_variable=X_FE,
    species={
        "Fe-Mg": Species("Fe-Mg", "Fe-Mg interdiffusion", "interdiffusion"),
        "Ca-Mg": Species("Ca-Mg", "Ca-(Mg,Fe) interdiffusion", "interdiffusion"),
    },
    notes=("Formula of the diopside-hedenbergite series: DHZ 2013, p. 95. Monoclinic, space "
           "group C2/c: DHZ 2013, p. 112 (augite has the more general formula "
           "(Ca,Mg,Fe2+,Al)2(Si,Al)2O6, p. 95). "
           "Mueller et al. (2013) measured along [001] only. Anisotropy of Fe-Mg in cpx "
           "is not well constrained (Cherniak & Dimanov 2010)."),
)

PLAGIOCLASE = Mineral(
    key="plagioclase", name="Plagioclase", formula="NaAlSi3O8-CaAl2Si2O8", system="triclinic",
    axes=("a", "b", "c"),
    composition_variable=X_AN,
    species={
        "Mg": Species("Mg", "Mg", "trace", "depends on An content"),
        "Sr": Species("Sr", "Sr", "trace", "depends on An content"),
        "Ba": Species("Ba", "Ba", "trace", "depends on An content"),
        "Li": Species("Li", "Li", "trace"),
        "NaSi-CaAl": Species("NaSi-CaAl", "coupled NaSi-CaAl interdiffusion", "interdiffusion",
                             "the anorthite profile itself, normally treated as frozen"),
    },
    isotropic=True,
    notes=("Formula (NaAlSi3O8-CaAl2Si2O8 with minor KAlSi3O8) and triclinic system: DHZ 2013, "
           "p. 292. Most plagioclases are ternary An-Ab-Or solid solutions with Or usually "
           "low (p. 292), hence the ternary X_An. Albite becomes monoclinic (C2/m) above "
           "about 950 C by a displacive transformation (p. 292, footnote). All other "
           "plagioclase is triclinic at all temperatures (p. 250). "
           "Trace-element diffusion is coupled to the anorthite gradient through the "
           "activity term of Costa et al. (2003). X_An is treated as frozen because "
           "NaSi-CaAl interdiffusion is orders of magnitude slower (Grove et al. 1984). "
           "Van Orman et al. (2014) found little anisotropy for Mg and recommend treating "
           "plagioclase as isotropic."),
)

KFELDSPAR = Mineral(
    key="kfeldspar", name="K-feldspar (sanidine, orthoclase)", formula="(K,Na)AlSi3O8",
    system="monoclinic", axes=("a", "b", "c"),
    composition_variable=X_OR,
    species={
        "Sr": Species("Sr", "Sr", "trace", "Cherniak (1996), sanidine Or61"),
        "Ba": Species("Ba", "Ba", "trace", "Cherniak (2002), sanidine Or61"),
        "Ti": Species("Ti", "Ti", "trace", "Cherniak & Watson (2020)"),
        "Na-K": Species("Na-K", "Na-K interdiffusion", "interdiffusion",
                        "the orthoclase fraction itself, Schaffer et al. (2014)"),
    },
    isotropic=True,
    notes=("Formula (K,Na)AlSi3O8 with minor CaAl2Si2O8: DHZ 2013, p. 253. High sanidine is "
           "monoclinic, C2/m (p. 253). Orthoclase is a highly ordered monoclinic K-feldspar "
           "(p. 250). Microcline, anorthoclase and albite are triclinic (p. 253), so the "
           "monoclinic system holds for sanidine and orthoclase only. "
           "Sr and Ba diffusion in sanidine showed no resolvable dependence on orientation "
           "(Cherniak 1996, 2002) and Ti little anisotropy (Cherniak & Watson 2020). Ba diffuses "
           "about 50 times slower than Sr (ratio of the Chamberlain et al. 2014 Table 1 laws at "
           "753-815 C), and Chamberlain et al. used the divergent behaviour of the two "
           "to constrain the initial width of the zone boundary."),
)

MAGNETITE = Mineral(
    key="magnetite", name="Magnetite / titanomagnetite", formula="(Ti_x Fe_(1-x))3 O4",
    system="cubic", axes=("a",),
    composition_variable=X_TI,
    species={
        "Fe-Ti": Species("Fe-Ti", "Fe-Ti interdiffusion", "interdiffusion"),
        "Ti": Species("Ti", "Ti", "tracer"),
        "Fe": Species("Fe", "Fe", "tracer"),
        "Mg": Species("Mg", "Mg", "tracer", "taken as equal to Fe by Tomiya et al. (2013)"),
        "Al": Species("Al", "Al", "tracer", "taken as equal to Ti by Tomiya et al. (2013)"),
        "Mn": Species("Mn", "Mn", "tracer"),
        "Co": Species("Co", "Co", "tracer"),
        "Cr": Species("Cr", "Cr", "tracer"),
    },
    isotropic=True,
    notes=("Cubic, so diffusion is isotropic. Spinel group, cubic, space group Fd3m. Magnetite "
           "Fe2+Fe3+2O4 and ulvospinel Fe2+2TiO4: DHZ 2013, p. 402. DHZ restrict the term "
           "titanomagnetite to specimens with a demonstrable ulvospinel component (p. 404). "
           "The (Ti_x Fe_(1-x))3 O4 formula and x_Ti = X_Usp/3 are the convention of the "
           "diffusion papers (x = 1/3 is Fe2TiO4), not a DHZ expression. "
           "Cation diffusion has a minimum with respect "
           "to both T and fO2 because a vacancy mechanism (D ~ fO2^(2/3)) competes with an "
           "interstitial mechanism (D ~ fO2^(-2/3)), see Van Orman & Crispin (2010) Table 12. "
           "Tomiya et al. (2013) used exactly this formulation for Shinmoedake 2011."),
)

GARNET = Mineral(
    key="garnet", name="Garnet", formula="(Fe,Mg,Mn,Ca)3Al2Si3O12", system="cubic",
    axes=("a",),
    composition_variable=X_FE,
    species={
        "Fe-Mg": Species("Fe-Mg", "Fe-Mg interdiffusion (binary)", "interdiffusion",
                         "almandine-pyrope garnet with little Mn and Ca"),
        "Fe": Species("Fe", "Fe tracer", "tracer"),
        "Mg": Species("Mg", "Mg tracer", "tracer"),
    },
    isotropic=True,
    notes=("Space group Ia3d, general formula X3Y2Si3O12: DHZ 2013, p. 18. The ideal "
           "formula (Mg,Fe2+,Mn,Ca)3(Al,Fe3+,Cr)2(Si,Al)3O12 is given in Appendix 3, p. 489 "
           "(the formula here is the Al-only almandine-pyrope-spessartine-grossular case). "
           "Cubic, so diffusion is isotropic. Fe, Mg, Mn and Ca zoning relaxes by coupled "
           "multicomponent exchange with large off-diagonal terms (Chakraborty & Ganguly 1992), "
           "which File > Multicomponent and isotope study models. The scalar Fe-Mg law here is "
           "for near-binary almandine-pyrope garnet."),
)

MINERALS = {m.key: m for m in (OLIVINE, ORTHOPYROXENE, CLINOPYROXENE, PLAGIOCLASE, KFELDSPAR,
                                   MAGNETITE, GARNET)}

_UNIAXIAL = ("quartz", "rutile", "apatite", "zircon", "xenotime")
# printed page of each mineral's data sheet in Deer, Howie & Zussman (2013)
_DHZ_PAGE = {
    "quartz": "DHZ 2013, p. 311 (alpha-quartz trigonal, P3121 or P3221, beta-quartz hexagonal)",
    "rutile": "DHZ 2013, p. 393 (tetragonal, P42/mnm)",
    "titanite": ("DHZ 2013, p. 15 (CaTi[SiO4](O,OH,F), monoclinic, P21/a. The end-member "
                 "CaTiSiO5 has O on the O(1) site that can also hold OH and F)"),
    "apatite": ("DHZ 2013, p. 473 (hexagonal, C63/m. The group formula is "
                "Ca5(PO4)3(OH,F,Cl), fluorapatite is Ca5(PO4)3F)"),
    "zircon": "DHZ 2013, p. 12 (Zr[SiO4], tetragonal, I41/amd)",
    "monazite": "DHZ 2013, p. 478 ((Ce,La,Th)PO4, monoclinic, P21/n)",
}
TRACE = CompositionVariable("none", "concentration", "measured trace-element concentration", default_mode="A")
for _key, _name, _formula, _system, _species in (
        ("quartz", "Quartz", "SiO2", "trigonal (alpha-quartz, beta-quartz is hexagonal)", ("Ti",)),
        ("rutile", "Rutile", "TiO2", "tetragonal", ("Zr", "Hf")),
        ("titanite", "Titanite", "CaTi(SiO4)(O,OH,F)", "monoclinic", ("Sr", "Zr")),
        ("apatite", "Fluorapatite", "Ca5(PO4)3F", "hexagonal", ("Sr", "Pb", "La", "Nd", "Sm", "Dy", "Yb")),
        ("zircon", "Zircon", "ZrSiO4", "tetragonal", ("Pb", "Sm", "Dy", "Yb", "Ti")),
        ("monazite", "Monazite", "(Ce,La,Th)PO4", "monoclinic", ("Pb",)),
        ("xenotime", "Xenotime", "YPO4", "tetragonal", ("Sm", "Dy", "Yb", "Pb"))):
    _note = "Direction restrictions are specific to the selected experimental law."
    if _key in _DHZ_PAGE:
        _note += f" Formula and crystal system: {_DHZ_PAGE[_key]}."
    if _key in _UNIAXIAL:
        _note += (" The crystal is uniaxial (two independent cell edges, a and c). The generic axis "
                  "labels a, b, c are Diffusor's convention and b is equivalent to a.")
    if _key == "quartz":
        _note += (" The trigonal system is that of alpha-quartz. DHZ 2013 (p. 311) give alpha-quartz as "
                  "stable up to 573 C and beta-quartz as stable from 573 to 870 C. The book states "
                  "no pressure for these limits. At magmatic temperatures quartz is therefore "
                  "the hexagonal beta form.")
    if _key == "monazite":
        _note += (" The formula is the general (Ce,La,Th)PO4 of Deer, Howie & Zussman (2013, p. 478). The laws "
                  "were measured on synthetic CePO4 or natural monazite as stated in each entry.")
    if _key == "xenotime":
        _note += (" Xenotime has no data sheet in Deer, Howie & Zussman (2013) and nothing here is taken from "
                  "that book (zircon, the isostructural mineral, is on p. 12). YPO4 is the composition of the synthetic "
                  "crystals of Cherniak (2006), who states that xenotime is isostructural with zircon, "
                  "hence the tetragonal system.")
    MINERALS[_key] = Mineral(_key, _name, _formula, _system, ("a", "b", "c"),
                            {s: Species(s, s, "trace") for s in _species}, TRACE, notes=_note)

for _sp in ("Ni", "Zn", "Sc", "Ga", "In", "Y", "Lu", "V3+", "V4+", "Zr", "Hf", "U", "Nb", "Ta", "Mo"):
    MAGNETITE.species[_sp] = Species(_sp, _sp, "trace")


def get_mineral(key: str) -> Mineral:
    try:
        return MINERALS[key]
    except KeyError as exc:
        raise KeyError(f"Unknown mineral '{key}'. Available: {list(MINERALS)}") from exc
