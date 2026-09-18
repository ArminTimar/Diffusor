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

1. **Data** loads your file or a bundled example. Diffusor guesses the distance
   column, its unit and the two composition columns from the names and asks you
   to confirm them. Uncertainties, oxides and a fit window sit under Advanced.
   An example also fills in the later steps, and each step shows what it set.
2. **Mineral** picks the phase, the diffusing species and the traverse
   orientation.
3. **Conditions** takes temperature, pressure and oxygen fugacity with their
   uncertainties, plus the analytical resolution. Pick the instrument and
   Diffusor fills in a beam sigma with the published convention it came from.
4. **Model** sets the geometry, the initial condition and what each end of the
   profile is: a plateau that continues, a crystal rim held by the melt, a
   closed rim or the crystal centre. The card states which solver the run will
   use and why.
5. **Coefficient** lists everything published for that mineral and species,
   tagged *recommended*, *unverified* or *superseded*.
6. **Uncertainty** chooses the Monte Carlo draws, the seed and what to sample.

The results view puts the fitted time and its interval at the top of a narrow
summary of every setting, with an *edit* link beside each group that jumps back
to the relevant step. The plot takes the rest of the width. Methods, the
coefficient details and the log open as formatted reading panes.

A Monte Carlo run draws itself while it works. Every draw adds a faint fitted
curve and, when measurement noise is being sampled, the perturbed data points,
so the cloud of points builds up into what an error bar really means here. Below
the profile the histogram of times fills in with the running median and the 68
per cent interval, and beside it the sampled temperatures, the temperature
against oxygen fugacity (which shows the buffer correlation directly) and the
sampled log D. A panel whose input is held fixed says so instead.

No page scrolls. Only lists and reading panes do, and scrolling never changes a
number or a drop-down. Values change only when you type or click.

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
anorthite column can be carried along and used to drive the activity term. An
optional fit window leaves points out of the fit without deleting them from the
file, for example a later overgrowth at the very rim. Any extra column such as
anorthite follows the same window.

## Example datasets

One per mineral, listed on the data step with their provenance. The catalogue
lives in `diffusor/datasets.py`.

| dataset | kind | source |
| --- | --- | --- |
| Plagioclase, Santorini Minoan | **measured** | Crystal S82-30A 12 from Supplementary Table 1 of Druitt et al. (2012), Nature 482:77-80 |
| Orthopyroxene Fe-Mg, Kizimen 2010 | **measured** | Crystal K9_L10C4 from Supplementary Data 2 of Ostorero et al. (2022), Commun. Earth Environ. 3:290 |
| Orthopyroxene Fe-Mg | synthetic | forward model (Dias et al. 2025) at conditions typical of an arc andesite |
| Clinopyroxene Fe-Mg | synthetic | forward model, conditions after Petrone et al. (2016, 2018) |
| Olivine Fe-Mg | synthetic | forward model, conditions after Hartley et al. (2016) |
| Titanomagnetite Ti | synthetic | forward model, exactly the conditions of Tomiya et al. (2013) at Shinmoedake |
| Clinopyroxene BSE greyscale | synthetic | forward model plus a linear grey response, with microprobe anchors |
| Sanidine Ba | synthetic | forward model (Cherniak 2002) at 790 °C, inside the Bishop Tuff range of Chamberlain et al. (2014) |

Two sets are real measured data. The synthetic ones were generated by
Diffusor's own forward model with a known time and are labelled as such in the
catalogue, in the interface and in every exported methods block.

Every example records where its temperature, pressure and oxygen fugacity come
from, and says plainly when a value is a placeholder that the chosen
coefficient never uses. The Kizimen conditions are the mean and standard
deviation of 21 magnetite-ilmenite pairs from the andesites, Supplementary Data
3 of Ostorero et al. (2022): 850 ± 57 °C and NNO +1.3 ± 0.35. The main text of
that paper swaps the andesite and dacite temperatures. The table does not.

The **Kizimen** set is a validation against a published timescale. Ostorero et
al. (2022) modelled this reverse zone at 850 °C with the Ganguly & Tazzoli
(1994) law and no fO2 term and obtained 2.32 years (+7.16/−1.75). With the same
law and temperature Diffusor fits about 3 years to the microprobe traverse,
inside their uncertainty, and the test suite checks it. It is not exact
because Ostorero et al. fitted high-resolution BSE grey-scale profiles, which
the supplement does not include, whereas the 2 µm microprobe traverse resolves
the 4 µm zone with only a handful of points. Switching the same crystal to the
newest opx calibration (Dias et al. 2025) gives roughly ten times longer, which
is worth knowing before choosing a law.

The **Santorini** set is shipped as a worked example of loading real data and of
how much the initial condition matters, **not** as a validation: Diffusor does
not reproduce the published 47-year timescale, because Druitt et al.
reconstructed their initial profile from the Sr-An correlation and a two-melt
history, which Diffusor's built-in initial conditions cannot express. The
catalogue says so.

## What it does

**Boundaries.** Each end of the profile is described in words rather than as a
mathematical condition. A plateau that continues stands in for a plateau of
infinite length, which is what the error-function solution assumes and what
Ostorero et al. (2022) and Chamberlain et al. (2014) used. A crystal rim held
by the melt is a fixed concentration, the open boundary of Costa et al. (2008,
p. 555). A closed rim and the crystal centre are both zero flux, the centre by
symmetry (Crank 1975, section 4.3), and a closed system equilibrates far faster
than an open one (Costa et al. 2008, Fig. 6). A crystal core is not truly held
at a fixed composition. While both plateaus survive, the far ends make no
difference to the fitted time, and Diffusor warns as soon as the diffusion front
reaches the end of the traverse.

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
operator. The centre of a cylinder or sphere uses the L'Hopital limit. The time
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
before comparing it with the data (Ganguly et al. 1988, Bradshaw & Kent 2017),
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
  spread of D at the temperature of interest by a large factor. Only the Sr and
  Ba laws of Grocolas et al. (2025) come with anything like a covariance: the
  authors state that log D0 and Q are strongly correlated and treat them as
  perfectly correlated in their own Monte Carlo, and Diffusor copies that
  assumption. No other paper in the registry publishes one, so those laws
  instead sample ln D directly at the working temperature using the scatter the
  source reports. The chooser greys out an option a law cannot support and names
  the one it will use. The `independent` mode is offered only so the older
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
| Opx Fe-Mg for Fs9 at 950, 1050 and 1100 C, log fO2 = -7 Pa (Dias & Dohmen 2024, Table 1) | fitted D0 and m of each run | within 0.2 log units, for both the 2024 and 2025 laws |
| Plagioclase Sr, gap to Giletti & Casserly (1994) at An36 and 750 C (Audétat et al. 2026) | 2.8 log units | 2.77 |
| Magnetite, minimum of D against fO2 for 21 elements at 1150 C (Sievwright et al. 2020, Table 5) | log fO2 and log D at the minimum | within 0.1 and 0.05 log units |
| Opx K9_L10C4, Kizimen (Ostorero et al. 2022, Supplementary Data 4) | 2.32 yr (+7.16/−1.75) | about 3 yr |

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
recovered time. The two measured files come from elsewhere. The
Kizimen traverse is extracted from the published spreadsheet by
`python scripts/extract_kizimen.py <folder>`.

## Coefficient registry

Fifty-four entries across six minerals, 47 of them read from the primary
publication or from a named re-fit. The recommended ones are:

| mineral | species | recommended source |
| --- | --- | --- |
| Orthopyroxene | Fe-Mg | Dias, Dohmen & Behrens (2025), two fO2 regimes |
| Clinopyroxene | Fe-Mg | Mueller et al. (2013) |
| Plagioclase | Mg | Van Orman, Cherniak & Kita (2014), and Audétat, Grocolas & Mutch (2026) with silica activity |
| Plagioclase | Sr, Ba | Grocolas, Bloch, Bouvier & Müntener (2025) |
| Plagioclase | Li | Pohl et al. (2024), interstitial mechanism |
| K-feldspar | Sr | Cherniak (1996), sanidine Or61 |
| K-feldspar | Ba | Cherniak (2002), same crystal as the Sr law |
| K-feldspar | Ti | Cherniak & Watson (2020) |
| Magnetite | Ti, Fe | Van Orman & Crispin (2010), Table 12 |
| Olivine | Fe-Mg | Dohmen & Chakraborty (2007), **unverified transcription** |

Sr and Ba in sanidine were measured on the same Or61 crystal by the same
method, so the pair can be compared directly: Ba is about 1.7 log units slower
at 800 °C, which is why paired profiles across one zone boundary test whether
that boundary is diffusive at all (Chamberlain et al. 2014, Audétat, Grocolas &
Mutch 2026). Both laws carry the log D uncertainties Chamberlain et al. (2014)
derived from them, 0.03 log units for Sr and 0.12 for Ba. Audétat, Grocolas &
Mutch (2026) report that Toussaint et al. (2025) find Sr in Or98 orthoclase
about 0.6 log units slower again, but that work is a conference abstract with no
Arrhenius law, so it is not in the registry.

Also available: the rare earth elements Lu, Ce and Eu in orthopyroxene from
Dias, Dohmen & Hartmann (2025), and magnetite Ti, Mn, Co, Cr, Al and Mg from
Sievwright et al. (2020). Sievwright et al. measured only at 1150 °C, so
Diffusor scales those entries to other temperatures with the Table 12
activation energies of the same element. Mg has no such energies and is valid
at 1150 °C only. The entry notes say plainly that the scaling is Diffusor's
construction, not the authors'.

Older calibrations are kept alongside them so published timescales can be
reproduced and compared: Dohmen et al. (2016) and Ganguly & Tazzoli (1994) for
orthopyroxene, the latter also in the no-fO2 form used by Ostorero et al.
(2022). Dimanov & Sautter (2000) for clinopyroxene, the coefficient behind the
published NIDIS results. Costa et al. (2003) for Mg in plagioclase. Giletti &
Casserly (1994), Cherniak & Watson (1994) and Cherniak (2002) for Sr and Ba in
plagioclase, the last two in the re-fitted form of Grocolas et al. (2025). And
Freer & Hauptman (1978) for Fe-Ti in titanomagnetite. `REFERENCES.md` lists all
of them with their verification status.

### Superseded coefficients

Entries known to be out of date are flagged in the registry, demoted from
*recommended*, marked in the chooser and called out in every fit.

**Sr in plagioclase.** Grocolas, Bloch, Bouvier & Müntener (2025), EPSL
651:119141, measured Sr diffusion in oligoclase and labradorite between 900 and
1200 °C with silica activity buffered, and found it one to two orders of
magnitude slower than Giletti & Casserly (1994) and Cherniak & Watson (1994).
They attribute the older, faster values to reaction fronts, because Sr-feldspar
was not stable in those source materials. The activation energy is higher, so
the gap grows as temperature falls: about 1.6 log units at 1100 °C, 2.2 at
900 °C and 2.8 at 750 °C for An36. Timescales from the 1990s calibrations are
therefore too short by roughly 100 to 600 times at 750-900 °C. Their equations
7 and 8 are implemented and recommended, and their Monte Carlo re-fits of the
older data (equations 12-14) replace two entries that Diffusor previously held
only as unverified transcriptions. Ba, by contrast, agrees with Cherniak (2002)
to within about half a log unit, so that entry is kept as a legitimate
alternative rather than flagged.

**Fe-Mg in orthopyroxene.** Dias & Dohmen (2024) and Dias, Dohmen & Behrens
(2025) re-fitted the Dohmen et al. (2016) experiments with a temperature-
dependent composition effect, added Fe-rich crystals, and found two regimes
separated at log fO2 = −10 Pa. For Fs10 at log fO2 = −7 Pa the new law agrees
with Dohmen et al. (2016) to within 0.3 log units, but the two differ by up to
one log unit for Fe-rich opx at 1100 °C or under reducing conditions, so Dohmen
et al. (2016) is marked superseded. Arc magmas usually sit above the calibrated
fO2 range, where the fO2^(1/4) term is extrapolated.

**Fe-Ti in magnetite.** No direct replacement for Freer & Hauptman (1978) or
Aragon et al. (1984) has been published, but for Ti and Fe *tracer* diffusion
the Aggarwal & Dieckmann (2002) data tabulated by Van Orman & Crispin (2010)
Table 12 are far better constrained, explicitly dependent on both temperature
and oxygen fugacity, and are what Tomiya et al. (2013) used at Shinmoedake.
Those are the recommended entries. At 1150 °C they agree with Sievwright et al.
(2020) to within half a log unit for Ti, Mn and Co.

## Known limitations

- A Monte Carlo on a composition-dependent law with the numerical solver takes
  a few seconds per draw, so 500 draws is a coffee break rather than a moment.
  The live plot shows how far along it is, and Stop ends it cleanly.
- Version 1 is 1-D only. Modelling a 3-D crystal in 1-D returns a maximum
  estimate, and sectioning adds further bias (Shea et al. 2015, Krimer & Costa
  2017). The application says so on every fit.
- The olivine entries were transcribed from secondary sources and are flagged
  unverified. Check them against Dohmen & Chakraborty (2007) and its erratum.
- Li in plagioclase is modelled with one effective coefficient per mechanism.
  Pohl et al. (2024) fitted a multispecies model with interstitial and
  lattice-site Li exchanging, which Diffusor does not implement.
- The Lu law of Dias et al. (2025) quotes its fO2 term relative to an fO2_0 that
  the paper does not define. Diffusor assumes the IW buffer and flags the entry
  unverified until that is checked.
- The Audétat et al. (2026) Mg law is taken from the accepted manuscript. Its
  silica activity is fixed at 1 (quartz saturation) unless overridden from
  Python. The interface does not expose it yet.
- Multi-component and isotopic diffusion are not implemented.

## Licence

MIT.
