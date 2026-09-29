"""Accessory-mineral lattice laws; distinct experimental exchange mechanisms.

These scalar entries predict concentration relaxation, not radiometric ages,
damage evolution, fluid replacement, or multi-site Li transport.
"""
import numpy as np
from functools import partial

from .literature import scalar as _scalar_factory, PRISTINE

COEFFICIENTS = []
scalar = partial(_scalar_factory, _registry=COEFFICIENTS)

scalar("ap_Sr_cherniak1993", "apatite", "Sr", "cherniak1993apatite", np.log10(2.7e-7),
       65000 * 4.184 / 1000, (700, 1050), source="p. 4657; abstract; D0 converted from cm2/s",
       notes=("Dry Durango fluorapatite, fit perpendicular to c; parallel measurements show additional scatter. "
              "Hydrous calibrations differ and are not interchangeable.",),
       reference_state="dry Durango fluorapatite; Sr implantation and SrO reservoir",
       uncertainty="Quoted Q ±2200 cal/mol and asymmetric D0 error; coefficient uncertainty not yet propagated.")
scalar("ap_Pb_cherniak1991", "apatite", "Pb", "cherniak1991apatite", np.log10(1.27e-8),
       54.6 * 4.184, (600, 900), source="abstract p. 1663; cm2/s and kcal/mol converted to SI", kind="tracer",
       notes=("Pb introduced by ion implantation; rapid annealing was inferred for apatite. "
              "This is not the damaged-zircon law from the same paper or an (U-Th)/He model.",),
       reference_state="dry apatite, implanted Pb",
       uncertainty="Quoted Q ±1.7 kcal/mol; joint uncertainty not transcribed and not sampled.")

for sp, D0, Q in (("La", 2.6e-7, 324), ("Nd", 2.4e-6, 348),
                  ("Dy", 9.7e-7, 340), ("Yb", 1.3e-8, 292)):
    scalar(f"ap_{sp}_cherniak2000_in", "apatite", sp, "cherniak2000apatite", np.log10(D0), Q,
           (800, 1250), source="abstract p. 3871, REE silicate oxyapatite source",
           reference_state="dry fluorapatite; REE silicate oxyapatite in-diffusion",
           notes=("Effective chemical diffusion for the stated charge-compensating exchange mechanism. "
                  "Do not substitute for isotope exchange or out-diffusion from Nd-doped apatite.",))
scalar("ap_Sm_cherniak2000_implant", "apatite", "Sm", "cherniak2000apatite", np.log10(6.3e-7), 298,
       (750, 1100), source="abstract p. 3871, ion-implantation relaxation", kind="tracer",
       reference_state="dry fluorapatite; implanted Sm",
       notes=("Implanted Sm exchange is faster than coupled REE chemical diffusion; select the mechanism matching the study.",))
scalar("ap_Nd_cherniak2000_out", "apatite", "Nd", "cherniak2000apatite", np.log10(9.3e-6), 392,
       (950, 1400), source="abstract p. 3871, synthetic Nd-doped apatite out-diffusion",
       reference_state="synthetic Nd-doped fluorapatite; undoped apatite reservoir",
       notes=("Out-diffusion without abundant charge-compensating species. Not interchangeable with the in-diffusion law.",))

scalar("zrn_Pb_cherniak2001", "zircon", "Pb", "cherniak2001zircon", np.log10(.11), 550,
       (1000, 1500), source="abstract p. 5", notes=(PRISTINE,),
       uncertainty="Quoted Q ±30 kJ/mol; joint uncertainty not transcribed and not sampled.")
for sp, logD0, Q in (("Sm", 8.46, 841), ("Dy", 5.36, 734), ("Yb", 7.40, 769)):
    scalar(f"zrn_{sp}_cherniak1997", "zircon", sp, "cherniak1997zircon", logD0, Q,
           (1150, 1400), source="abstract p. 289, low-temperature RBS fits",
           notes=(PRISTINE, "No significant anisotropy resolved for REE; REE-phosphate source, synthetic and natural zircon."))

zti = scalar("zrn_Ti_bloch2022_c", "zircon", "Ti", "bloch2022zircon", 1.34, 555.425,
             (1100, 1540), source="abstract p. 1", reference_axis="c", allowed_axes=("c",),
             notes=(PRISTINE, "Weak concentration dependence was observed; this is the paper's effective Arrhenius fit. "
                    "Diffusion transverse to c is many orders of magnitude slower."), kind="effective")
zti.orientation_required = True
zti = scalar("zrn_Ti_cherniak2007_perp_c", "zircon", "Ti", "cherniak2007zircon", np.log10(333), 754,
             (1350, 1550), source="abstract p. 470, 1-atm fit", reference_axis="a", allowed_axes=("a", "b"),
             notes=(PRISTINE, "Perpendicular to c only. The faster c-direction law is a separate entry; no fixed axis ratio is valid."))
zti.orientation_required = True

scalar("mnz_Pb_cherniak2004", "monazite", "Pb", "cherniak2004monazite", np.log10(.94), 592,
       (1100, 1350), source="abstract p. 829", notes=(PRISTINE,
           "Dry synthetic CePO4 and natural monazite. Fluid-mediated recrystallization can reset ages without lattice diffusion."))

# Xenotime experiments measured normals to (101), not a crystallographic axis.
# No crystallographic axis is assigned and the limitation stays explicit.
for sp, D0, Q, temps in (("Sm", 1.5e-4, 441, (1052, 1400)), ("Dy", 9e-8, 349, (1000, 1400)),
                         ("Yb", 3.9e-7, 362, (1002, 1400)), ("Pb", 3e-9, 382, (1200, 1500))):
    c = scalar(f"xtm_{sp}_cherniak2006", "xenotime", sp, "cherniak2006xenotime", np.log10(D0), Q,
               temps, source="abstract pp. 1-2; Table 1 for per-species temperature ranges",
               notes=(PRISTINE, "Measured normal to (101) only. No angular dependence or tensor is available."),
               reference_state="dry synthetic xenotime; transport normal to (101)")
    # Sentinel prevents any a/b/c or direction-cosine request being accepted.
    c.allowed_axes = ("normal_(101)",)
