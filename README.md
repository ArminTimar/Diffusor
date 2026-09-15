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

The window has six control groups on the left, the plot in the middle and
results, methods and a log on the right.

1. **Data** loads a CSV, TSV or Excel table and asks which columns are the
   distance, the two elements and their uncertainties.
2. **Mineral and species** picks the phase, the diffusing species and the
   traverse orientation, either a named axis or angles to *a*, *b* and *c*.
3. **Conditions** takes temperature, pressure and oxygen fugacity with their
   uncertainties. Oxygen fugacity is normally given as an offset from a buffer.
4. **Model** sets the geometry, the boundary conditions, the initial condition
   and the grid, and can switch to a linear cooling path.
5. **Diffusion coefficients** lists everything published for that
   mineral-species pair. Tick one to fit, several to compare.
6. **Monte Carlo** chooses the number of draws, the seed and which sources of
   uncertainty to sample.

Fitting is fast. A Monte Carlo of 500 draws takes a few minutes and runs on a
background thread, so the window stays responsive and can be stopped.

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
| Plagioclase | Sr | Giletti & Casserly (1994) |
| Magnetite | Ti, Fe | Van Orman & Crispin (2010), Table 12 |
| Magnetite | Fe-Ti | Freer & Hauptman (1978) |
| Olivine | Fe-Mg | Dohmen & Chakraborty (2007), **unverified transcription** |

Older calibrations are kept alongside them so published timescales can be
reproduced and compared: Ganguly & Tazzoli (1994) for orthopyroxene, Dimanov &
Sautter (2000) for clinopyroxene (the coefficient behind the published NIDIS
results), Costa et al. (2003) for Mg in plagioclase, Aragon et al. (1984) for
Fe-Ti in titanomagnetite. `REFERENCES.md` lists all of them with their
verification status.

## Known limitations

- Version 1 is 1-D only. Modelling a 3-D crystal in 1-D returns a maximum
  estimate, and sectioning adds further bias (Shea et al. 2015; Krimer & Costa
  2017). The application says so on every fit.
- The olivine entries were transcribed from secondary sources and are flagged
  unverified. Check them against Dohmen & Chakraborty (2007) and its erratum.
- The 2025 orthopyroxene recalibration of Dias, Dohmen & Behrens is recorded in
  the reference list but not implemented; add it to
  `diffusor/coefficients/opx.py` once the paper is to hand.
- Multi-component and isotopic diffusion are not implemented.

## Licence

MIT.
