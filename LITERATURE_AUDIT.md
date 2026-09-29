# Literature audit and implementation status

Updated 29 September 2026. This is the current implementation record; the supplied
20 September gap analysis is a research roadmap, not an authoritative specification.
The local library contains 106 PDFs. The inventory records file hashes and extraction
status; indexing is not equivalent to reading, validating or implementing a paper.
No copyrighted full text is committed to this repository.

## Corrections established from the sources

- **Olivine Fe–Mg:** Dohmen & Chakraborty (2007), pp. 424–425, equations 27–28,
  contain `+3 XFe`, not `+3 (XFe-0.1)`. The latter lowered D by 0.3 dex and
  doubled inferred times under otherwise equivalent conditions. Both branches
  are corrected. The erratum has not been checked in full; these legacy entries
  remain unverified and TaMED is no longer recommended by default.
- **Ni in olivine:** Petry et al. (2004), p. 4184, give the fixed-fO2 Fo90 fit
  `D0=3.84e-9 m²/s, Q=216 kJ/mol`. Its intercept must not be paired with the
  separate global activation energy of 220 kJ/mol. The implemented fit is restricted
  to the stated host composition and redox state, with extrapolation warnings.
  The sixfold anisotropy measurement is at 1200 °C; extending that ratio is explicitly
  an approximation. A complete composition/redox Ni model remains to be added.
- **Mn/pressure in olivine:** Holzapfel et al. (2007) do not justify a generic Mn
  Arrhenius entry assembled from Fe–Mg parameters. Activation volumes at fixed fO2
  differ from apparent values along a buffer. Mn remains unimplemented.
- **Magnetite:** Sievwright et al. (2020), Table 5, measures one temperature,
  1150 °C. All table rows now have fixed-temperature entries. The five existing
  temperature-scaled entries remain for reproducibility, explicitly labelled as
  hypotheses. V3+ and V4+ rows are speciation interpretations of the same total-V
  fit, not two independently calibrated populations for coupled transport.
- **Anisotropy:** independent Ca and Be Arrhenius functions must be evaluated
  before projection. Zircon Ti c and transverse laws remain separate; no constant
  ratio is imposed. An unspecified direction is not a powder average. A plane normal
  in a non-orthogonal crystal cannot simply be renamed a crystal axis.
- **Mineral metadata:** quartz is trigonal. Xenotime data normal to (101) are labelled
  as such; requests to use their unknown full tensor are rejected.
- **Garnet attribution:** the 2024 *Bridging the Gap in Garnet Diffusion Models at
  Low Temperatures* paper is by Junxing Chen and Xu Chu, not Fan et al.
- **Statistics and geometry:** ln D0 and positive Q are normally positively
  correlated. Independent sampling omits covariance; missing covariance is not zero
  covariance. A 1-D inferred time is not a universal upper bound on a 3-D history.

## Implemented source families

Each entry carries its mechanism, transported variable, calibration ranges,
reference state, equation location, uncertainty status and caveats into the chooser,
JSON, methods text and workbook. `verified=True` denotes source transcription checks;
it does not certify a natural sample, a mechanism, or complete uncertainty propagation.

| Source | Implementation | Principal restrictions |
| --- | --- | --- |
| Coogan et al. 2005 | Olivine Ca, a/b/c equations and fO2^0.31 | Fo83–92; independent marginal-error sampling approximates missing covariance |
| Jollands et al. 2016 | Olivine Be, a/b/c equations | Fo100 principal fits; explicit direction required |
| Watson et al. 2015 | Olivine P | 650–850 °C; growth zoning is not proof of diffusion |
| Petry et al. 2004 | Olivine Ni restricted fit | Fo90, fO2=1e-6 Pa; incomplete global model |
| Sievwright et al. 2020 | All 21 table rows at 1150 °C | Temperature changes rejected for fixed-temperature laws |
| Cherniak et al. 2007 | Quartz Ti | c-direction only; competing calibrations require comparison |
| Cherniak et al. 2007 | Rutile Zr/Hf | Direction-specific fits; no unmeasured transverse Zr law |
| Cherniak 1995, 2006 | Titanite Sr and Zr | Dry experiments; Zr c-direction only |
| Cherniak & Ryerson 1993; Cherniak et al. 1991 | Apatite Sr and Pb | cm²/s and cal/mol converted explicitly; dry experiments |
| Cherniak 2000 | Apatite La, Nd, Dy, Yb in-diffusion; implanted Sm; Nd out-diffusion | Distinct charge-compensation and experimental mechanisms |
| Cherniak & Watson 2001; Cherniak et al. 1997 | Zircon Pb and Sm/Dy/Yb | Undamaged crystalline lattice; no radiogenic production |
| Bloch et al. 2022; Cherniak & Watson 2007 | Zircon Ti c and transverse fits | Explicit orientation; no universal anisotropy ratio |
| Cherniak et al. 2004 | Monazite Pb | Dry lattice; fluid-assisted replacement excluded |
| Cherniak 2006 | Xenotime Pb/Sm/Dy/Yb | Normal to (101); no inferred tensor |

Published uncertainties are not silently converted into a joint probability model.
Entries without verified uncertainty conventions remain fixed during coefficient
sampling and state that omission explicitly. A fixed-temperature coefficient cannot
be used for a Monte Carlo that varies temperature.

## Calculation and application changes

The olivine coordinate selector distinguishes XFe, XFo and Fo mol%. For an Fo-percent
profile the coefficient receives `XFe=1-C/100`, while the solver and plots retain the
measured units. Trace-element concentrations do not become host Fe fractions.

Numerical transport now assembles one flux per face, including the activity term,
and applies its equal and opposite contributions to neighbouring control volumes.
Exact shell volumes are used for cylinder/sphere operators and mass diagnostics.
Closed boundaries impose zero total external flux. Invalid grids/diffusivities and
exhausted step budgets raise errors. A time-varying boundary forces the numerical
route. A supplied isothermal history uses its own temperature in the analytical
integral instead of silently using the nominal conditions.

**File → Shared-duration study** stores independent profiles with their own models
and fits a common duration by minimising the sum of squared residuals divided by
measurement variance. Every profile needs finite positive uncertainties in its own
units. The joint degrees of freedom are total point count minus one; per-profile
contributions are also reported. Conditions, boundaries, orientation, initial states
and coefficients stay fixed. This is not a posterior sampler or a coupled diffusion
matrix. Residual panels expose inconsistencies between the common-event assumption
and individual profiles. Common systematic errors and cross-profile measurement
covariance are not yet represented.

Exports include numeric Excel results, source/caveat sheets, profile and residual
charts, CSV, JSON, methods text, and PNG/SVG plots. The workbook is a completed-run
snapshot, not an independently recalculating diffusion solver. Joint studies also
export each constituent profile's complete provenance. The 68% profile envelope is
an ensemble band, not an observation prediction interval.

## Remaining roadmap, with implementation gates

| Earlier item | Current disposition | Required evidence / implementation |
| --- | --- | --- |
| 0A Mn and general Ni | Partial | Full redox/composition model and pressure conventions, independent numerical benchmarks |
| 0B magnetite table | Implemented with restrictions | Do not manufacture thermal dependence for new rows |
| 1A Cr and Li in olivine | Pending | Valence/site populations and exchange rates, not an effective Fe–Mg law |
| 1B quartz | Cherniak c law implemented | Audit alternative calibrations and review's natural-profile criteria |
| 1C spinel Fe–Mg | Pending; Liermann/Ganguly and Vogt PDFs indexed and inspected | Distinguish self/tracer D from chemical interdiffusion; reproduce thermodynamic/defect model |
| 1D rutile | Principal laws implemented | Validate cooling/residence examples against Blackburn et al. |
| 1E apatite | Dry Cherniak families implemented | Watson et al. hydrous alternatives require separate verification |
| 1F titanite | Sr/Zr implemented | Nd substitution/redox, Pb damage; Nb/Ta primary PDF absent from indexed library |
| 1G feldspar alkalis | Pending; Schäffer PDF inspected | Transcribe composition dependence and plane-normal directions; benchmark exchange couples |
| 2A zircon | Pb/REE/Ti implemented | No damage evolution or radiogenic ages |
| 2B zircon Li | Pending | Coupled modes, REE/Y charge balance and isotope transport |
| 2C monazite/xenotime | Lattice laws implemented | Dissolution/reprecipitation is a separate reaction problem |
| 3A garnet | Pending | Vector state, site/charge constraints, thermodynamic matrix, nonlinear positivity, calibration-family benchmarks |
| 3B Li/H multisite | Pending | Mobile/immobile populations, reaction stoichiometry and validated exchange kinetics |
| 3C real 2-D/3-D geometry | Pending | Tensor rotation, mesh/surface boundaries and section extraction; 2026 PDF text extraction hit a decompression limit |
| 3D growth/resorption | Pending | Moving material coordinates, accreted-shell ages and melt mass balance; existing callable fixed-domain boundaries are supported |
| 3E stable isotopes | Pending | Isotope-specific fields and mass dependence; constraints on beta from the actual experiments |
| 3F joint/Bayesian | Shared-time likelihood implemented | Correlated measurement errors, common systematic uncertainty, nuisance parameters, priors/posterior diagnostics |
| 4 thermochronology | Pending | Radiogenic production, decay/recoil/ejection, damage/annealing and benchmarked age calculation |

Books and review papers guide derivations; they do not substitute for experimental
calibration tables. The full roadmap remains unfinished. New model classes must pass
mass/site balance, limiting-case, grid/time convergence and primary-paper reproduction
tests before being offered as named mineral methods.

## Reproducibility and verification

`docs/literature_inventory.json` records the PDF inventory and SHA256 checksums.
`scripts/index_literature.py` regenerates extracted text under ignored `tmp/`.
The supplementary Ostorero spreadsheets were located; existing Kizimen extraction
continues to provide the measured example. They are not treated as additional
independent validations of unrelated laws.

Regression tests cover printed equations and units, changing directional ratios,
unsupported orientation/temperature rejection, composition-coordinate equivalence,
thermal-history integration, radial/activity mass conservation, known-duration joint
fits, concentration-unit invariance, incompatible-profile residuals, cancellation,
and workbook/JSON provenance. Existing Crank analytical-solution and measured-example
checks remain part of the full suite. Synthetic recovery validates software behaviour;
it does not independently validate geological interpretation.

## Garnet matrix implementation note (next stage)

Borinski et al. (2012), p. 573, equation 1 was checked against the rendered PDF.
For equal-valence Fe–Mg–Mn–Ca components in the ideal-solution limit, with component
n dependent and `sum(X)=1`, it reduces to

`D_ij = D_i* delta_ij - (D_i* X_i / sum_k(D_k* X_k)) (D_j* - D_n*)`,

for the independent components i,j=1,…,n−1. Thus a published tracer D_i* is not
itself a diagonal chemical-diffusion matrix entry. The equal-tracer limit gives
`D_ij=D* delta_ij`; the binary limit gives
`D_inter = D_1* D_2* / (X_1 D_1* + X_2 D_2*)`. These are necessary limiting-case
checks for a future implementation of this particular constrained ionic model.
They are not universal substitutions for other diffusion mechanisms.

Equation 2 adds activity-coefficient derivatives; equations 3–5 and Table 1 define
one thermodynamic model. These derivatives must respect the dependent component.
The reduced Fick matrix need not be symmetric and may contain negative off-diagonal
entries. Matrix validation must follow the mobility/thermodynamic formulation rather
than rejecting physically meaningful cross terms. The paper's discussion distinguishes
Fe–Mg-rich and Mn-rich compositional vectors, and identifies larger non-ideality effects
for unusual Mg–Mn–Ca-rich compositions. No generic scalar "garnet D" has been registered.

Implementation still requires conservative vector fluxes, site-balance and positivity
checks, thermal/pressure conventions, separately audited tracer calibration families,
and reproduction of the paper's diffusion couples. Shared-time fitting above does
not supply these missing equations.

## Validated working checkpoint

Final full suite: **187 passed** (29 September 2026), plus successful module compilation
and a clean whitespace check with the repository's CRLF convention. The remaining
warnings are the existing Qt high-DPI deprecation emitted by Matplotlib.

A synthetic paired apatite Sr–Nd study at 1050 °C recovered **97.66 years** from
profiles generated at **100 years** with independent Gaussian noise and seed 20260929.
The exported study workbook was reopened and checked against the JSON duration and
121 joint degrees of freedom; individual workbooks contain both charts, all 61
measurements and source metadata. Report plots and the joint-study interface were
rendered and visually inspected. This is a software-recovery demonstration, not an
experimental validation of the apatite calibrations.
