"""Minerals available in the coefficient registry.

Sources for the structural/compositional conventions: Deer, Howie & Zussman
(1992); the diffusion papers themselves define the composition variable used
by each coefficient (see :mod:`diffusor.coefficients`).
"""
from __future__ import annotations

from .base import CompositionVariable, Mineral, Species

# --- composition variables ------------------------------------------------------
X_FE = CompositionVariable(
    "XFe", "X_Fe", "Fe/(Fe+Mg) molar, all Fe as Fe2+",
    default_mode="A/(A+B)", citation="deer1992")
X_FO = CompositionVariable(
    "XFo", "X_Fo", "Mg/(Mg+Fe) molar (forsterite fraction)",
    default_mode="B/(A+B)", citation="deer1992")
X_AN = CompositionVariable(
    "XAn", "X_An", "Ca/(Ca+Na+K) molar (anorthite fraction)",
    default_mode="A/(A+B)", citation="deer1992")
X_OR = CompositionVariable(
    "XOr", "X_Or", "K/(K+Na+Ca) molar (orthoclase fraction)",
    default_mode="A", citation="deer1992")
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
    notes=("Diffusion is strongly anisotropic: D[001] is about 6x D[100] and D[010] "
           "(Dohmen & Chakraborty 2007). Hartley et al. (2016) confirm the factor of 6 in "
           "natural crystals)."),
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
    notes=("Space group Pbca. D//[001] >= D//[010] > D//[100]. Dohmen et al. (2016) "
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
    notes=("Mueller et al. (2013) measured along [001] only. Anisotropy of Fe-Mg in cpx "
           "is not well constrained (Cherniak & Dimanov 2010)."),
)

PLAGIOCLASE = Mineral(
    key="plagioclase", name="Plagioclase", formula="(Ca,Na)(Al,Si)4O8", system="triclinic",
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
    notes=("Trace-element diffusion is coupled to the anorthite gradient through the "
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
    notes=("Sr and Ba diffusion in sanidine showed no resolvable dependence on orientation "
           "(Cherniak 1996, 2002) and Ti little anisotropy (Cherniak & Watson 2020). Ba diffuses "
           "about 50 times slower than Sr, so paired Sr and Ba profiles across the same zone "
           "boundary test whether the boundary is diffusive (Chamberlain et al. 2014)."),
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
    notes=("Cubic, so diffusion is isotropic. Cation diffusion has a minimum with respect "
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
    notes=("Cubic, so diffusion is isotropic. Fe, Mg, Mn and Ca zoning relaxes by coupled "
           "multicomponent exchange with large off-diagonal terms (Chakraborty & Ganguly 1992), "
           "which File > Multicomponent and isotope study models. The scalar Fe-Mg law here is "
           "for near-binary almandine-pyrope garnet."),
)

MINERALS = {m.key: m for m in (OLIVINE, ORTHOPYROXENE, CLINOPYROXENE, PLAGIOCLASE, KFELDSPAR,
                                   MAGNETITE, GARNET)}

TRACE = CompositionVariable("none", "concentration", "measured trace-element concentration", default_mode="A")
for _key, _name, _formula, _system, _species in (
        ("quartz", "Quartz", "SiO2", "trigonal", ("Ti",)),
        ("rutile", "Rutile", "TiO2", "tetragonal", ("Zr", "Hf")),
        ("titanite", "Titanite", "CaTiSiO5", "monoclinic", ("Sr", "Zr")),
        ("apatite", "Fluorapatite", "Ca5(PO4)3F", "hexagonal", ("Sr", "Pb", "La", "Nd", "Sm", "Dy", "Yb")),
        ("zircon", "Zircon", "ZrSiO4", "tetragonal", ("Pb", "Sm", "Dy", "Yb", "Ti")),
        ("monazite", "Monazite", "CePO4", "monoclinic", ("Pb",)),
        ("xenotime", "Xenotime", "YPO4", "tetragonal", ("Sm", "Dy", "Yb", "Pb"))):
    MINERALS[_key] = Mineral(_key, _name, _formula, _system, ("a", "b", "c"),
                            {s: Species(s, s, "trace") for s in _species}, TRACE,
                            notes="Direction restrictions are specific to the selected experimental law.")

for _sp in ("Ni", "Zn", "Sc", "Ga", "In", "Y", "Lu", "V3+", "V4+", "Zr", "Hf", "U", "Nb", "Ta", "Mo"):
    MAGNETITE.species[_sp] = Species(_sp, _sp, "trace")


def get_mineral(key: str) -> Mineral:
    try:
        return MINERALS[key]
    except KeyError as exc:
        raise KeyError(f"Unknown mineral '{key}'. Available: {list(MINERALS)}") from exc
