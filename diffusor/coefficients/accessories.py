"""Accessory-mineral lattice laws; distinct experimental exchange mechanisms.

These scalar entries predict concentration relaxation, not radiometric ages,
damage evolution, fluid replacement, or multi-site Li transport.
"""
import numpy as np
from functools import partial

from .literature import scalar as _scalar_factory, PRISTINE

COEFFICIENTS = []
scalar = partial(_scalar_factory, _registry=COEFFICIENTS)

# Calories in the apatite papers are converted with the thermochemical calorie (4.184 J). The papers do not
# define the calorie; the international-table calorie (4.1868 J) would raise Q by 0.067 %.
CAL = 4.184
CAL_NOTE = ("Q was printed in cal/mol and is converted with the thermochemical calorie (4.184 J), which the "
            "paper does not define; the 4.1868 J calorie would raise Q by 0.07 %.")

scalar("ap_Sr_cherniak1993", "apatite", "Sr", "cherniak1993apatite", np.log10(2.7e-7),
       65000 * CAL / 1000, (700, 1050), source="p. 4657; abstract; D0 converted from cm2/s",
       p_stated=True,
       notes=("Dry Durango fluorapatite. The law is the fit for transport perpendicular to c (abstract); the paper "
              "gives a separate fit parallel to c, which is not implemented (parallel diffusivities are generally "
              "smaller, by up to 0.4 log units, with significant scatter). "
              "Hydrous calibrations differ and are not interchangeable.", CAL_NOTE),
       reference_state="dry Durango fluorapatite; Sr implantation and SrO reservoir",
       uncertainty="Quoted Q ±2200 cal/mol and asymmetric D0 error (confidence level not stated); "
                   "coefficient uncertainty not yet propagated.")
scalar("ap_Pb_cherniak1991", "apatite", "Pb", "cherniak1991apatite", np.log10(1.27e-8),
       54.6 * CAL, (600, 900), source="abstract p. 1663; cm2/s and kcal/mol converted to SI", kind="chemical",
       kind_note=("The transport kind ('chemical') is Diffusor's classification: Pb is ion-implanted and its "
                  "profile relaxes; the authors do not label the coefficient."),
       notes=("Pb introduced by ion implantation; rapid annealing was inferred for apatite. "
              "This is not the damaged-zircon law from the same paper or an (U-Th)/He model. The fit is the "
              "paper's own lower-temperature fit (abstract), not its 'global' fit with Watson et al. (1985). "
              "The diffusion direction is not stated: the apatite was cut parallel to the c-axis (some samples "
              "were prism faces), so the profile direction is probably perpendicular to c.", CAL_NOTE),
       reference_state="implanted Pb in natural apatite, annealed in air (dryness not stated by the authors); "
                       "no pressure correction",
       uncertainty="Quoted Q ±1.7 kcal/mol (confidence level not stated); joint uncertainty not "
                   "transcribed and not sampled.")

for sp, D0, Q in (("La", 2.6e-7, 324), ("Nd", 2.4e-6, 348),
                  ("Dy", 9.7e-7, 340), ("Yb", 1.3e-8, 292)):
    scalar(f"ap_{sp}_cherniak2000_in", "apatite", sp, "cherniak2000apatite", np.log10(D0), Q,
           (800, 1250), source="abstract p. 3871, REE silicate oxyapatite source", kind_stated=True,
           reference_state="anhydrous fluorapatite (stated by the author); REE silicate oxyapatite in-diffusion",
           notes=("Effective chemical diffusion for the stated charge-compensating exchange mechanism. "
                  "Do not substitute for isotope exchange or out-diffusion from Nd-doped apatite.",))
# The author attributes the implantation law to simple light-REE3+ <-> REE3+ exchange without charge
# compensation and contrasts it with "REE chemical diffusion" (coupled substitution), so it is classified
# as a tracer-type exchange coefficient (this adds the tracer advice to the warnings; D is unchanged).
scalar("ap_Sm_cherniak2000_implant", "apatite", "Sm", "cherniak2000apatite", np.log10(6.3e-7), 298,
       (750, 1100), source="abstract p. 3871, ion-implantation relaxation", kind="tracer",
       kind_note=("The transport kind ('tracer') is Diffusor's classification: Cherniak (2000) attributes the "
                  "implantation law to simple light-REE3+ <-> REE3+ exchange without charge-compensating species, "
                  "used for isotopic gradients, and contrasts it with 'REE chemical diffusion' by coupled "
                  "substitution. He does not use the word tracer, and the Durango host already holds native REE."),
       reference_state="anhydrous Durango fluorapatite (stated by the author); implanted Sm",
       notes=("Implanted Sm exchange (simple REE3+ <-> REE3+) is faster than coupled REE chemical diffusion; "
              "select the mechanism matching the study.",))
scalar("ap_Nd_cherniak2000_out", "apatite", "Nd", "cherniak2000apatite", np.log10(9.3e-6), 392,
       (950, 1400), source="abstract p. 3871, synthetic Nd-doped apatite out-diffusion",
       reference_state="synthetic Nd-doped fluorapatite (anhydrous, stated by the author); undoped apatite reservoir",
       notes=("Out-diffusion without abundant charge-compensating species. Not interchangeable with the in-diffusion law. "
              "The slabs were cut parallel to the c-axis; the diffusion direction is not stated and no axis is "
              "assigned.",))

scalar("zrn_Pb_cherniak2001", "zircon", "Pb", "cherniak2001zircon", np.log10(.11), 550,
       (1000, 1500), source="abstract p. 5 (D0 = 1.1e-1; the results text prints log D0 = -0.962)",
       p_stated=True,
       reference_state="1-atm runs sealed under vacuum in silica capsules; the authors find no resolvable effect of "
                       "water or pressure; no pressure correction",
       notes=(PRISTINE,
              "log D0 is taken from the rounded D0 = 1.1e-1 of the abstract (log10 = -0.959); the results text prints "
              "log D0 = -0.962, which lowers D by 0.8 %.",
              "The Arrhenius fit uses the in-diffusion data (1000-1302 C); the stated 1000-1500 C range "
              "also covers out-diffusion runs up to 1500 C that agree with it."),
       uncertainty="Quoted Q ±30 kJ/mol and log D0 = -0.962 ±1.074 (confidence level not stated); joint "
                   "uncertainty not transcribed and not sampled.")
# Per-element data ranges (Tables 2 and 3): the abstract's 1150-1400 C is the union for the three REE.
for sp, logD0, Q, temps in (("Sm", 8.46, 841, (1200, 1400)), ("Dy", 5.36, 734, (1150, 1350)),
                            ("Yb", 7.40, 769, (1150, 1350))):
    scalar(f"zrn_{sp}_cherniak1997", "zircon", sp, "cherniak1997zircon", logD0, Q,
           temps, source="abstract p. 289, low-temperature RBS fits",
           reference_state="annealed in air in Pt capsules (dryness not stated by the authors); "
                           "no pressure correction",
           notes=(PRISTINE, "No significant anisotropy resolved for REE; REE-phosphate source, synthetic and natural zircon.",
                  "Temperature range: the abstract gives 1150-1400 C for the three REE together; the data of this "
                  f"law span {temps[0]}-{temps[1]} C" +
                  (" (the 1650 C electron-microprobe point is excluded; a fit including it gives 691 ±47 kJ/mol)."
                   if sp == "Yb" else ".")))

zti = scalar("zrn_Ti_bloch2022_c", "zircon", "Ti", "bloch2022zircon", 1.34, 555.425,
             (1100, 1540), source="abstract p. 1", reference_axis="c", allowed_axes=("c",), p_stated=True,
             reference_state="1 atm, Ni-NiO buffered, parallel to c; water content not stated; no pressure correction",
             notes=(PRISTINE, "Weak concentration dependence was observed; this is the paper's effective Arrhenius fit. "
                    "Diffusion transverse to c is many orders of magnitude slower.",
                    "The printed equation divides by 2.303 (a rounded ln 10); Diffusor uses the exact ln 10, which "
                    "lowers log D by 0.003-0.004 (0.7-0.9 % in D)."), kind="effective",
             kind_note=("The transport kind ('effective') is Diffusor's classification (the fit lumps a weak "
                        "concentration dependence); the authors do not name the kind."))
zti.orientation_required = True
zti = scalar("zrn_Ti_cherniak2007_perp_c", "zircon", "Ti", "cherniak2007zircon", np.log10(333), 754,
             (1350, 1550), source="abstract p. 470, 1-atm fit", reference_axis="a", allowed_axes=("a", "b"),
             p_stated=True, kind_stated=True,
             reference_state="anhydrous, 1 atm (stated by the authors); no pressure correction",
             notes=(PRISTINE, "Perpendicular to c only. The faster c-direction law is a separate entry; no fixed axis ratio is valid."))
zti.orientation_required = True

scalar("mnz_Pb_cherniak2004", "monazite", "Pb", "cherniak2004monazite", np.log10(.94), 592,
       (1100, 1350), source="abstract p. 829", p_stated=True,
       reference_state="dry, 1 atm (stated by the authors); no pressure correction",
       notes=(PRISTINE,
              "Dry synthetic CePO4 and natural monazite. Fluid-mediated recrystallization can reset ages without lattice diffusion.",
              "The law is the fit to the synthetic CePO4 in-diffusion data (RBS) only; a fit to all synthetic and "
              "natural in-diffusion data gives 563 ±35 kJ/mol (D0 = 0.089 m2/s). The results text and the Fig. 3 "
              "caption print 594 kJ/mol for the synthetic fit; the abstract and the conclusions print 592 ±39, which "
              "Diffusor uses."))

# Xenotime experiments measured normals to (101), not a crystallographic axis.
# No crystallographic axis is assigned and the limitation stays explicit.
for sp, D0, Q, temps in (("Sm", 1.5e-4, 441, (1052, 1400)), ("Dy", 9e-8, 349, (1000, 1400)),
                         ("Yb", 3.9e-7, 362, (1002, 1400)), ("Pb", 3e-9, 382, (1200, 1500))):
    extra = ()
    if sp == "Pb":
        extra = ("Temperature range: Table 1 lists Pb experiments from 1200 to 1400 C; the experimental section and "
                 "Fig. 2 extend to 1500 C, which Diffusor uses. D0 = 3.0e-9 is the abstract value; the results "
                 "text and the Fig. 2 caption print 2.95e-9 (log D0 = -8.529), which lowers D by 1.7 %.",)
    c = scalar(f"xtm_{sp}_cherniak2006", "xenotime", sp, "cherniak2006xenotime", np.log10(D0), Q,
               temps, source="abstract pp. 1-2; Table 1 for per-species temperature ranges", p_stated=True,
               kind_note=("The transport kind ('chemical') is Diffusor's classification: diffusion from a "
                          "single-REE phosphate (or PbTiO3) source by exchange for Y3+; the author does not label "
                          "the coefficient as chemical or tracer."),
               notes=(PRISTINE, "Measured normal to (101) only. No angular dependence or tensor is available.", *extra),
               reference_state="dry synthetic xenotime; transport normal to (101)")
    # Sentinel prevents any a/b/c or direction-cosine request being accepted.
    c.allowed_axes = ("normal_(101)",)
