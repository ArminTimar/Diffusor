# Diffusor

A local desktop application for diffusion chronometry. It solves 1-D diffusion
problems analytically and numerically, carries a registry of literature
diffusion coefficients for olivine, orthopyroxene, clinopyroxene, plagioclase
and magnetite, fits the diffusion time to a measured profile, and propagates
uncertainties by Monte Carlo without pretending that temperature, oxygen
fugacity and the Arrhenius parameters are independent.

Everything is traceable. Every equation, coefficient, constant and convention
carries a citation key, the test suite fails if a key does not resolve, and
every run exports a methods block listing the sources it actually used.

## Install

```bash
python -m pip install -e .
```

Needs Python 3.10 or later. Dependencies are numpy, scipy, pandas, matplotlib,
openpyxl and PySide6.

## Run

```bash
python -m diffusor
```

The window walks through six steps and then opens a results view.

1. **Data** loads your file, or one of the bundled example datasets. This step
   also states the input format in full, so you never have to guess.
2. **Mineral** picks the phase, the diffusing species and the traverse
   orientation.
3. **Conditions** takes temperature, pressure and oxygen fugacity with their
   uncertainties, plus the analytical resolution.
4. **Model** sets the geometry, the boundaries and the initial condition.
5. **Coefficient** lists everything published for that mineral and species,
   tagged *recommended*, *unverified* or *superseded*.
6. **Uncertainty** chooses the Monte Carlo draws, the seed and what to sample.

The results view puts the fitted time and its interval at the top of a narrow
summary of every setting, with an *edit* link beside each group that jumps back
to the relevant step. The plot takes the rest of the width. Methods, the Monte
Carlo histogram and the log open on demand rather than occupying the window.

## Input format

One row per measurement point, one column of distance, one or two columns of
composition, and optionally an uncertainty column for each. Column names are
free: you map them after loading. CSV, TSV, plain text and Excel are accepted.

```
Distance_um,FeO_wt,MgO_wt,FeO_err,MgO_err
0.0,20.9,22.4,0.15,0.20
1.7,20.8,22.5,0.15,0.20
```

With two composition columns Diffusor forms the molar ratio A/(A+B), which is
what the diffusion coefficients are calibrated against. Name the oxides when
mapping and it converts weight per cent to cation moles first. With one column
the values are modelled as they stand, which suits a trace element in ppm, a
forsterite content in mol% or a calibrated grey value. Distances may run either
way and need not be evenly spaced. Extra columns are ignored, so a plagioclase
anorthite column can be carried along and used to drive the activity term.

## Example datasets

One per mineral, listed on the data step with their provenance. The catalogue
lives in `diffusor/datasets.py`.

| dataset | kind | source |
| --- | --- | --- |
| Plagioclase, Santorini Minoan | **measured** | Crystal S82-30A 12 from Supplementary Table 1 of Druitt et al. (2012), Nature 482:77-80 |
| Orthopyroxene Fe-Mg | synthetic | forward model, conditions after Tomiya et al. (2013) and Sato et al. (2022) |
| Clinopyroxene Fe-Mg | synthetic | forward model, conditions after Petrone et al. (2016, 2018) |
| Olivine Fe-Mg | synthetic | forward model, conditions after Hartley et al. (2016) |
| Titanomagnetite Ti | synthetic | forward model, exactly the conditions of Tomiya et al. (2013) at Shinmoedake |
| Clinopyroxene BSE greyscale | synthetic | forward model plus a linear grey response, with microprobe anchors |

Only the Santorini set is real measured data. The synthetic ones were generated
by Diffusor's own forward model with a known time and are labelled as such in
the catalogue, in the interface and in every exported methods block. The
Santorini set is shipped as a worked example of loading real data and of how
much the initial condition matters, **not** as a validation: Diffusor does not
reproduce the published 47-year timescale, because Druitt et al. reconstructed
their initial profile from the Sr-An correlation and a two-melt history, which
Diffusor's built-in initial conditions cannot express. The catalogue says so.

## What it does

**Analytical solutions** come from Crank (1975), *The Mathematics of
Diffusion*, 2nd edition: the semi-infinite step (eq. 2.14), the finite band
(eq. 2.15), the plane sheet (eq. 4.17), the cylinder (eq. 5.22), the sphere
(eq. 6.18) and the corresponding fractional-uptake series. They are used
automatically whenever they are exactly valid, and the application says which
route it took and why.

**Numerical solutions** use a Crank-Nicolson scheme on a uniform grid with the
diffusion coefficient evaluated at half-nodes, following Crank (1975)
section 8.4 and the finite-difference appendix of Dohmen, Faak & Blundy (2017).
Plane, cylindrical and spherical geometry are handled by a single conservative
operator; the centre of a cylinder or sphere uses the L'Hopital limit. The time
step starts at the explicit stability limit and grows geometrically, which
keeps a 100 kyr run to a few hundred steps while still resolving the sharp
initial transient. Against the closed forms the solver is accurate to about
3 parts in a million and converges at second order, and a closed system
conserves mass to 1 part in 10^8.

**Composition-dependent diffusivity** is supported throughout, so the fitted
profile itself sets the local D. For plagioclase trace elements the flux
carries the activity term of Costa et al. (2003), driven by the frozen
anorthite gradient.

**Anisotropy** uses the direction-cosine relation of Costa & Chakraborty
(2004), with the per-axis ratios published for each mineral.

**Analytical resolution** is handled by convolving the model with a Gaussian
before comparing it with the data (Ganguly et al. 1988; Bradshaw & Kent 2017),
so the fitted time is corrected for beam smearing rather than inflated by it.

## Error propagation

This is the part most existing tools get wrong, and the reason the project
exists.

The time retrieved from a profile is a strongly non-linear function of
parameters that are not independent:

- Oxygen fugacity is normally known as an offset from a mineral buffer, and the
  buffer is itself a function of temperature. Diffusor re-evaluates the buffer
  at every sampled temperature, so the correlation is exact by construction.
- The Arrhenius parameters ln D0 and Q are strongly anti-correlated by the
  regression that produced them. Sampling them independently inflates the
  spread of D at the temperature of interest by a large factor. Diffusor
  samples from the published covariance where it exists, and otherwise samples
  ln D directly at the working temperature using the scatter the source
  actually reports. The `independent` mode is offered only so the older
  assumption can be reproduced for comparison.
- Composition enters both the diffusion coefficient and the profile being
  fitted, so measurement noise propagates by two routes at once.

Every draw re-runs the whole fit, so all of this is honoured automatically.
Times are log-normally distributed, so the median and the 16th, 84th, 2.5th and
97.5th percentiles are reported rather than a symmetric standard deviation.
A variance decomposition re-runs the Monte Carlo with one source active at a
time, which ranks what is worth measuring better.

## Traceability

- `diffusor/references.py` is the single citation registry. `REFERENCES.md` and
  `references.bib` are generated from it by `scripts/build_references.py`.
- Every `DiffusionCoefficient` records the equation as printed in its source,
  the equation number, the calibration ranges, the published uncertainties and
  a `verified_from` field saying exactly where the numbers were read.
- Coefficients whose primary publication was not available offline are flagged
  **unverified** in the registry, shown in red in the application and called
  out in every exported methods block. Check them against the paper before
  publishing.
- Exporting writes the results as JSON, the profile and model as CSV, the
  figure as PNG and SVG, the Monte Carlo times as CSV, and a methods block with
  the full reference list for that run.

## Verification against published values

The test suite checks the registry against numbers printed in the sources:

| check | published | Diffusor |
| --- | --- | --- |
| Opx Fe-Mg at 950 C, log fO2 = -7 Pa, //c (Dohmen et al. 2016, run OPXD_14) | log D = -19.49 ± 0.07 | -19.47 |
| Opx anisotropy D//[001] / D//[100] (Dohmen et al. 2016) | 3.5 | 3.5 |
| Cpx Fe-Mg at 1098 C and 1150 C (Petrone et al. 2016, Table 2 footnote) | 3.26e-20 and 1.20e-19 m²/s | 3.25e-20 and 1.19e-19 |
| Titanomagnetite Ti at 950 C, log fO2 = -11, X_Usp = 0.3 (Tomiya et al. 2013, Shinmoedake) | 4.3e-16 m²/s | 4.35e-16 |
| Titanomagnetite Ti at 900 C, same conditions | 6.9e-16 m²/s | 6.84e-16 |

The magnetite Fe branch reproduces Tomiya et al. (2013) to 11 per cent at
900 C but is a factor of 2.1 higher at 950 C, because they place the diffusion
minimum at a different temperature. That discrepancy is documented in the
coefficient's notes rather than hidden.

Run the tests with:

```bash
python -m pytest
```

## Examples

`examples/` holds synthetic profiles generated by the forward model with known
times, so the whole chain can be checked against an answer that is known
exactly. Regenerate them with `python scripts/make_examples.py`. The round-trip
tests in `tests/test_examples_roundtrip.py` load each one, fit it and check the
recovered time.

## Coefficient registry

Thirty-four entries across five minerals. The best-constrained are:

| mineral | species | recommended source |
| --- | --- | --- |
| Orthopyroxene | Fe-Mg | Dohmen, ter Heege, Becker & Chakraborty (2016) |
| Clinopyroxene | Fe-Mg | Mueller et al. (2013) |
| Plagioclase | Mg | Van Orman, Cherniak & Kita (2014) |
| Plagioclase | Sr | none recommended, see below |
| Magnetite | Ti, Fe | Van Orman & Crispin (2010), Table 12 |
| Magnetite | Fe-Ti | Aragon et al. (1984), or the Table 12 tracer entries |
| Olivine | Fe-Mg | Dohmen & Chakraborty (2007), **unverified transcription** |

Older calibrations are kept alongside them so published timescales can be
reproduced and compared: Ganguly & Tazzoli (1994) for orthopyroxene, Dimanov &
Sautter (2000) for clinopyroxene (the coefficient behind the published NIDIS
results), Costa et al. (2003) for Mg in plagioclase, Freer & Hauptman (1978)
for Fe-Ti in titanomagnetite. `REFERENCES.md` lists all of them with their
verification status.

### Superseded coefficients

Entries known to be out of date are flagged in the registry, demoted from
*recommended*, marked in the chooser and called out in every fit.

**Sr in plagioclase.** Grocolas, Bloch, Bouvier & Müntener (2025), EPSL
651:119141, measured Sr diffusion in oligoclase and labradorite between 900 and
1200 °C with silica activity buffered, and found it **1.5 to 2 orders of
magnitude slower** than Giletti & Casserly (1994) and Cherniak & Watson (1994),
which they attribute to feldspar stability not having been controlled in the
earlier experiments. Timescales from the 1990s calibrations are therefore
likely too short by a factor of roughly 30 to 100. Their Arrhenius parameters
are **not implemented**: the paper is open access but could not be retrieved
offline here, so the numbers would have had to be invented. Add them to
`diffusor/coefficients/plagioclase.py` from the PDF. Ba diffusion, by contrast,
they found similar to the earlier work, so the Ba entry is not flagged.

**Fe-Ti in magnetite.** No direct replacement for Freer & Hauptman (1978) or
Aragon et al. (1984) has been published, but for Ti and Fe *tracer* diffusion
the Aggarwal & Dieckmann (2002) data tabulated by Van Orman & Crispin (2010)
Table 12 are far better constrained, explicitly dependent on both temperature
and oxygen fugacity, and are what Tomiya et al. (2013) used at Shinmoedake.
Those are now the recommended entries. Sievwright et al. (2020) add modern
magnetite diffusivities for Ti and many other elements against fO2, but only at
1150 °C, so they give no activation energy and cannot be extrapolated to
magmatic temperatures on their own.

## Known limitations

- Version 1 is 1-D only. Modelling a 3-D crystal in 1-D returns a maximum
  estimate, and sectioning adds further bias (Shea et al. 2015; Krimer & Costa
  2017). The application says so on every fit.
- The olivine entries were transcribed from secondary sources and are flagged
  unverified. Check them against Dohmen & Chakraborty (2007) and its erratum.
- Two recent calibrations are recorded in the reference list but not
  implemented, because their numbers could not be retrieved without the PDFs:
  the orthopyroxene recalibration of Dias, Dohmen & Behrens (2025) and the
  plagioclase Sr and Ba diffusivities of Grocolas et al. (2025).
- Multi-component and isotopic diffusion are not implemented.

## Licence

MIT.
