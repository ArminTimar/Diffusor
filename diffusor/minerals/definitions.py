"""The five minerals supported in version 1.

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
X_TI = CompositionVariable(
    "xTi", "x_Ti", "Ti per cation site in (Ti_x Fe_(1-x))3 O4; x_Ti = X_Usp/3",
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
    },
    notes=("Diffusion is strongly anisotropic: D[001] is about 6x D[100] and D[010] "
           "(Dohmen & Chakraborty 2007; Hartley et al. 2016 confirm the factor ~6 in "
           "natural crystals)."),
)

ORTHOPYROXENE = Mineral(
    key="opx", name="Orthopyroxene", formula="(Mg,Fe)2Si2O6", system="orthorhombic",
    axes=("a", "b", "c"),
    composition_variable=X_FE,
    species={
        "Fe-Mg": Species("Fe-Mg", "Fe-Mg interdiffusion", "interdiffusion"),
        "Mg": Species("Mg", "Mg self-diffusion", "tracer"),
        "Lu": Species("Lu", "Lu", "trace", "rare earth element; Dias et al. (2025)"),
        "Ce": Species("Ce", "Ce", "trace", "rare earth element; Dias et al. (2025)"),
        "Eu": Species("Eu", "Eu", "trace", "rare earth element, assumed trivalent"),
    },
    notes=("Space group Pbca. D//[001] >= D//[010] > D//[100]; Dohmen et al. (2016) "
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
    notes=("Mueller et al. (2013) measured along [001] only; anisotropy of Fe-Mg in cpx "
           "is not well constrained (Cherniak & Dimanov 2010)."),
)

PLAGIOCLASE = Mineral(
    key="plagioclase", name="Plagioclase", formula="(Ca,Na)(Al,Si)4O8", system="triclinic",
    axes=("a", "b", "c"),
    composition_variable=X_AN,
    species={
        "Mg": Species("Mg", "Mg", "trace", "An-dependent; needs the activity term"),
        "Sr": Species("Sr", "Sr", "trace", "An-dependent; needs the activity term"),
        "Ba": Species("Ba", "Ba", "trace", "An-dependent; needs the activity term"),
        "Li": Species("Li", "Li", "trace"),
        "NaSi-CaAl": Species("NaSi-CaAl", "coupled NaSi-CaAl interdiffusion", "interdiffusion",
                             "the anorthite profile itself; normally treated as frozen"),
    },
    isotropic=True,
    notes=("Trace-element diffusion is coupled to the anorthite gradient through the "
           "activity term of Costa et al. (2003); X_An is treated as frozen because "
           "NaSi-CaAl interdiffusion is orders of magnitude slower (Grove et al. 1984). "
           "Van Orman et al. (2014) found little anisotropy for Mg and recommend treating "
           "plagioclase as isotropic."),
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
           "interstitial mechanism (D ~ fO2^(-2/3)); see Van Orman & Crispin (2010) Table 12. "
           "Tomiya et al. (2013) used exactly this formulation for Shinmoedake 2011."),
)

MINERALS = {m.key: m for m in (OLIVINE, ORTHOPYROXENE, CLINOPYROXENE, PLAGIOCLASE, MAGNETITE)}


def get_mineral(key: str) -> Mineral:
    try:
        return MINERALS[key]
    except KeyError as exc:
        raise KeyError(f"Unknown mineral '{key}'. Available: {list(MINERALS)}") from exc
