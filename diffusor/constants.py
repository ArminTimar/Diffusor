"""Physical constants and atomic masses used throughout Diffusor.

Every value is traceable:

* Fundamental constants: CODATA 2018 recommended values
  (Tiesinga et al. 2021, Rev. Mod. Phys. 93, 025010) -- citekey ``codata2018``.
* Atomic weights: IUPAC Commission on Isotopic Abundances and Atomic
  Weights, "Standard atomic weights 2021" (Prohaska et al. 2022, Pure Appl.
  Chem. 94, 573-600) -- citekey ``iupac2021``. Conventional values are used
  for elements with an interval notation.
"""

# --- fundamental constants (CODATA 2018) --------------------------------------
K_BOLTZMANN = 1.380649e-23   # J K^-1          (exact, 2019 SI redefinition)
N_AVOGADRO = 6.02214076e23   # mol^-1          (exact, 2019 SI redefinition)
# The molar gas constant is the exact product N_A k_B = 8.31446261815324... J mol^-1 K^-1.
# CODATA 2018 prints it as 8.314 462 618..., a truncation at nine decimals; the product is
# used here so that the value is exact to double precision (the difference, 1.84e-11 relative,
# has no numerical consequence).
R_GAS = N_AVOGADRO * K_BOLTZMANN   # J mol^-1 K^-1
EV_TO_J = 1.602176634e-19    # J per eV        (exact)

# --- unit conversions ----------------------------------------------------------
T_KELVIN_OFFSET = 273.15     # K  (t/degC = T/K - 273.15, the SI definition of the Celsius scale)
BAR_TO_PA = 1.0e5            # exact
GPA_TO_PA = 1.0e9            # exact
KBAR_TO_PA = 1.0e8           # exact
SEC_PER_MIN = 60.0
SEC_PER_HOUR = 3600.0
SEC_PER_DAY = 86400.0
SEC_PER_YEAR = 365.25 * SEC_PER_DAY   # Julian year a = 365.25 d = 31.5576 Ms, IAU Style Manual (citekey ``iau_julian_year``)
CM2_TO_M2 = 1.0e-4
UM_TO_M = 1.0e-6

# --- atomic weights, g mol^-1 (IUPAC 2021 Table 1: the standard atomic weight where it is
# a single value (Na, Al, Mn, Rb ...), the conventional (abridged) value where it is an
# interval (O, Mg, Si ...)) ----------------------------------------------------
ATOMIC_MASS = {
    "O": 15.999,
    "Na": 22.98976928,
    "Mg": 24.305,
    "Al": 26.9815384,
    "Si": 28.085,
    "K": 39.0983,
    "Ca": 40.078,
    "Ti": 47.867,
    "Cr": 51.9961,
    "Mn": 54.938043,
    "Fe": 55.845,
    "Ni": 58.6934,
    "Li": 6.94,
    "Sr": 87.62,
    "Ba": 137.327,
    "Rb": 85.4678,
}

# --- oxide molar masses derived from the atomic weights above ---------------
def _oxide(cations: dict, n_o: int) -> float:
    return sum(ATOMIC_MASS[el] * n for el, n in cations.items()) + n_o * ATOMIC_MASS["O"]

OXIDE_MOLAR_MASS = {
    "SiO2": _oxide({"Si": 1}, 2),
    "TiO2": _oxide({"Ti": 1}, 2),
    "Al2O3": _oxide({"Al": 2}, 3),
    "Cr2O3": _oxide({"Cr": 2}, 3),
    "FeO": _oxide({"Fe": 1}, 1),
    "Fe2O3": _oxide({"Fe": 2}, 3),
    "MnO": _oxide({"Mn": 1}, 1),
    "MgO": _oxide({"Mg": 1}, 1),
    "NiO": _oxide({"Ni": 1}, 1),
    "CaO": _oxide({"Ca": 1}, 1),
    "Na2O": _oxide({"Na": 2}, 1),
    "K2O": _oxide({"K": 2}, 1),
    "Li2O": _oxide({"Li": 2}, 1),
    "SrO": _oxide({"Sr": 1}, 1),
    "BaO": _oxide({"Ba": 1}, 1),
}

# number of cations per formula unit of each oxide (for wt% -> cation moles)
OXIDE_CATIONS = {
    "SiO2": 1, "TiO2": 1, "Al2O3": 2, "Cr2O3": 2, "FeO": 1, "Fe2O3": 2,
    "MnO": 1, "MgO": 1, "NiO": 1, "CaO": 1, "Na2O": 2, "K2O": 2, "Li2O": 2,
    "SrO": 1, "BaO": 1,
}
