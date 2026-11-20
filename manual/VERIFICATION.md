# Diffusor verification

This document lists how the diffusion laws, the solvers and the example datasets
are checked against published values. The models themselves are described in
[METHODS.md](METHODS.md).

## Test suite

The test suite contains 549 tests and runs on Windows, macOS and Linux with
Python 3.10 and 3.13 on every change to the repository. It does not access the
network. The tests are run with:

```bash
python -m pip install -e .[dev,images]
python -m pytest
```

The tests cover citation integrity, comparison with published values, solver
accuracy, fitting and Monte Carlo, and input, output and the interface. Each is
described below.

**Citation integrity.** Every coefficient must cite an existing reference key
for its primary and secondary sources, contain the transcribed equation, and
state where in the source its values were read. Every measured example dataset
and every superseded entry must cite an existing key, every reference must
render as text and as BibTeX, and an unknown key raises an error.

**Solver accuracy.** The finite-difference solver is compared with the
closed-form solutions of Crank (1975) for the semi-infinite step, the plane
sheet, the cylinder and the sphere. On an 801-node grid the largest deviation
is 3 to 14 parts per million of the concentration step (14 for the plane sheet,
3 to 4 for the others), and the tests require less than 1e-4. The error
decreases by a factor of about four with each halving of the node spacing
(second-order convergence). A closed plane system conserves mass to 1e-8 and a
closed sphere to 1e-6, and a closed system relaxes to its mean composition. The
explicit and implicit schemes agree with each other, the fractional-uptake
series agree with the integrated profiles, and a constant temperature history
gives the same result as the isothermal solution. The multicomponent solver
agrees with the eigen-component solution for a constant matrix (Toor 1964, as
used by Chakraborty & Ganguly 1992) within 1e-3 in mole fraction on a 401-node
grid, conserves each component to 1e-9, reduces to the scalar solver for two
components, and gives the same profiles whichever component is taken as
dependent. With β = 0 the isotope models give zero delta values and the element
profiles of the scalar model.

**Fitting and Monte Carlo.** Fits recover a known time from noise-free and noisy
synthetic profiles, with and without composition-dependent D. The bundled
synthetic examples are fitted from file and the recovered time is compared with
the known one (see [Example datasets](#example-datasets)). Further tests check
that ignoring beam convolution lengthens the apparent time, that a zero
uncertainty budget gives no spread, that the spread grows with the temperature
uncertainty, that fO2 follows the buffer when temperature is sampled and along
a cooling path, that independent sampling of ln D<sub>0</sub> and Q gives a
wider interval than sampling log D at temperature, and that parallel runs
reproduce the single-core draws exactly.

**Input, output and interface.** Unit and oxide conversions, buffer equations,
grey-value calibration, image reading and profile extraction, the contents of
the exported JSON and Excel files, and the main interface workflows (built off
screen) are also tested.

## Source transcription

160 of the 161 registry entries have been compared with the rendered pages of
their primary publications. The Fe-Ti law of Aragon et al. (1984) is taken from
Table 11 of Van Orman & Crispin (2010) and is marked as not checked against the
original paper.

In October 2026 every coefficient, equation, unit, range, uncertainty and
attributed statement in the code (2,574 items) was compared with its source
publication. The items were inventoried with OpenAI ChatGPT and compared by
Anthropic Claude models (Claude Sonnet 5.5 for the item checks, Claude Opus 5.5
for second checks of disagreements), which read numbers and equations from
rendered pages of the source PDFs and recomputed every unit conversion and
derived value in Python. The disagreements found were corrected in version
0.3.0. Values that a source does not state are labelled in the code as
Diffusor's assumptions. The comparison was not repeated item by item by a
person, so transcription errors may remain. Please report any you find.

## Comparison with published values

The table below mixes measured runs, worked evaluations and comparisons between
independently fitted models. A global fit need not reproduce an individual
measured D exactly. For OPXD_14 the printed log D is −19.49 ± 0.07, whereas the
implemented fit gives −19.473888. Passing a tolerance test establishes
approximate agreement, not exact source transcription. Source coefficients,
equation forms and units require separate checks.

| check | published | Diffusor |
| --- | --- | --- |
| Opx Fe-Mg at 950 °C, log fO2 = −7 (Pa), parallel to c (Dohmen et al. 2016, run OPXD_14) | log D = −19.49 ± 0.07 | −19.47 (0.016 above the run) |
| Opx anisotropy D[001]/D[100] (Dohmen et al. 2016) | 3.5 | 3.5 |
| Cpx Fe-Mg at 1098 and 1150 °C (Petrone et al. 2016, Table 2) | 3.26e-20 and 1.20e-19 m²/s | 3.262374e-20 and 1.198555e-19 m²/s, with D<sub>0</sub> = 9.55e-5 m²/s. The Table 2 footnote prints 9.5e-5, which gives 3.25e-20 and 1.19e-19 |
| Titanomagnetite Ti at 950 °C, log fO2 = −11, X<sub>Usp</sub> = 0.3 (Tomiya et al. 2013) | 4.3e-16 m²/s | 4.15e-16 m²/s (3 % lower) |
| Titanomagnetite Ti at 900 °C, same conditions | 6.9e-16 m²/s | 6.83e-16 m²/s (1 % lower) |
| Titanomagnetite Fe at 950 °C, same conditions (printed for Mg, taken as Fe, Tomiya et al. 2013) | 4.4e-15 m²/s | 4.17e-15 m²/s (5 % lower) |
| Titanomagnetite Fe at 900 °C, same conditions | 6.6e-15 m²/s | 6.54e-15 m²/s (1 % lower) |
| Opx Fe-Mg for Fs9 at 950 to 1102 °C, nominal log fO2 = −7 (Pa) (Dias & Dohmen 2024, Table 1) | log D at X<sub>Fe</sub> = 0.09 from the fitted D<sub>0</sub> and m of each of five runs. Measured fO2 differs slightly from nominal | four runs (950, 1050, 1050 and 1102 °C) within 0.10 log units of the 2024 and 2025 laws. The 1000 °C run lies 0.37 (2024 law) and 0.34 (2025 law) log units below them |
| Olivine Ni at 1005 °C, Fo90 (Petry et al. 2004, experiment Ni10) | log D = −17.12 | −17.24 (0.12 lower) |
| Plagioclase Sr, difference from Giletti & Casserly (1994) at An36 and 750 °C (Audétat et al. 2026) | 2.8 log units | 2.77 |
| Magnetite, minimum of D with fO2 for 21 elements at 1150 °C (Sievwright et al. 2020, Table 5) | log fO2 and log D at the minimum, as printed | minimum computed from the two printed branches differs by at most 0.09 (log fO2) and 0.04 (log D) |
| Opx K9_L10C4, Kizimen (Ostorero et al. 2022, Supplementary Data 4) | 2.32 yr (+7.16/−1.75) | 2.5 yr |
| Olivine Fe-Mg, 18 Kilauea 2020 profiles refitted with the archived inputs (Lynn et al. 2024) | 14 to 521 d, ± 30 % | 0.75 to 1.06 times the published time, largest RMS difference from the archived model curve 0.10 mol% Fo |
| Olivine Fo ages from the fitted widths and Dt products of Tables SM4-A, SM5-A and SM6-A (Gordeychik et al. 2018) | 32 ages | relative difference below 4e-16 |
| Opx Fe-Mg perpendicular to c at 966 °C, NNO, after Dohmen et al. (2016) (Araya et al. 2024, section 5.2) | log D = −19.78 | −19.782 (a axis) |
| Olivine 3311_2_ol6, Fo only, interface position fitted (Mourey et al. 2023) | 395 d (+149/−106) | 406 d, with equation 3 of the paper read in bar |
| Olivine Fe-Mg, 59 Mauna Loa 2022 profiles, the paper's law, sheet conditions and angles, step position and plateaus fitted (Lynn et al. 2024, Nature Commun.) | MCMC medians, 1.1 to 232 d | 58 of 59 inside the published 95 % interval, median ratio 1.07 |
| Sanidine Ba, 23 greyscale and X-ray traverses at 930 °C, Cherniak (2002) law, sharp step with free plateaus (Iovine et al. 2017, Tables 1 and 2) | 1 to 151 yr | 0.90 to 1.28 times the published time, 18 of 23 equal after rounding to whole years |
| Sr, Ba and Ti times at T ± 30 °C (Chamberlain et al. 2014, Electronic Appendix 7), from the Table 1 laws | 222 rows | relative difference below 1e-4 |
| Cpx Fe-Mg times from √(4Dt) and T with the Dimanov & Sautter (2000) law, D0 = 9.55e-5 m²/s (Petrone et al. 2018, Table S4, semi-infinite rows) | 71 times | 63 within 5 % or rounding, median ratio 1.0005 |
| Garnet D matrix at X<sub>Fe</sub> 0.61, X<sub>Mn</sub> 0.20, X<sub>Mg</sub> 0.18, X<sub>Ca</sub> 0.01 (Chakraborty & Ganguly 1992, eq. 5) | nine printed elements | within 0.45 %, with tracer coefficients chosen to reproduce them (see below) |
| Olivine D<sub>Fe-Mg</sub> from the Fe and Mg tracer fits, 1100 to 1250 °C, a, b and c (Oeser et al. 2026, eq. 6 and Table 4) | the paper's own D<sub>Fe-Mg</sub> fit | within 0.2 log units (largest 0.18) |
| Garnet Mn tracer at 510 °C, 2 GPa, eclogite garnet (Chen & Chu 2024) | about 1e-24 m²/s (an order of magnitude, not a fitted value) | 2.2e-24 m²/s |

The Kizimen, Kilauea, Mauna Loa, Gordeychik, Araya, Mourey, Iovine, Chamberlain
and Petrone rows come from the published-study validation in
[examples/validation/report.md](../examples/validation/report.md), which also
lists what each downloaded supplement does not contain. Two reconstructions of
Sundermeyer et al. (2020) do not reproduce the published ages (49 and 579 d
against 277 and 358 d at Piton de la Fournaise), because the authors' fit
windows, initial profiles and, for the Eifel case, fO2 are not published. They
are reported as such and are not offered as examples. The extractions and runs
were prepared with OpenAI ChatGPT and rechecked with Anthropic Claude (Claude
Opus 5.5), which reran every script and compared the inputs with the source
workbooks and with rendered pages of the papers. The Iovine, Chamberlain and
Petrone checks and the Lynn Mauna Loa, Sato and Morgado traverses were added
later. Claude Haiku subagents read the papers and supplements, and Claude Opus
5.5 checked the laws, conditions and published times on rendered pages and ran
the refits.

The garnet check tests the form of the ideal-solution matrix (eq. 2), not the
tracer coefficients. The paper does not print the tracer coefficients behind
eq. 5. Its text refers to Table 1, which holds microprobe analyses, and Table 2
covers only 1200 °C and 20 to 35 kb. The values used here, D*<sub>Mn</sub> =
2.495e-12, D*<sub>Mg</sub> = D*<sub>Fe</sub> = 6.0e-13 and D*<sub>Ca</sub> =
0.5 D*<sub>Fe</sub> cm²/s, were chosen to reproduce the printed matrix. Only
D*<sub>Ca</sub> = 0.5 D*<sub>Fe</sub> is stated by the authors. Their eq. 3 at
the 41 kb and 1430 °C of that experiment gives 1.52e-12, 4.46e-13 and 4.52e-13
cm²/s for Mn, Mg and Fe, and with these every element of the matrix is 17 to
36 % below the printed one. Two printed elements differ from the computed ones
in the last digit (D[Fe,Mn] −8.231e-13 against −8.24e-13, D[Fe,Mg] −1.125e-13
against −1.12e-13), which is rounding.

The titanomagnetite rows compare the Ti and Fe entries with worked values that
Tomiya et al. (2013) printed. The two sets agree within 5.2 % (largest at
950 °C for Fe). Tomiya et al. do not state the unit of fO2 or how they obtained
the x<sub>Ti</sub> dependence. Diffusor takes log fO2 in bar and interpolates
each branch of Table 12 of Van Orman & Crispin (2010) log-linearly in
x<sub>Ti</sub> between the printed values for x<sub>Ti</sub> = 0 and 0.2, with
x<sub>Ti</sub> = X<sub>Usp</sub>/3 = 0.1. Tomiya et al. place the Fe minimum at
about 950 °C, and the entry has it at 955 °C.

Further tests reproduce the Arrhenius parameters, equation forms and figures of
the other sources, for example the olivine composition term of the Dohmen &
Chakraborty (2007) erratum and the Fig. 3 line of Grove et al. (1984).

## Published validation library

The **File > Published validation library** browser includes 492 measured
profiles from twelve studies of olivine, sanidine, orthopyroxene,
clinopyroxene and titanomagnetite. The Lynn (2024) Kilauea single-event
examples also load the published conditions, orientations and tabulated
initial states, the Lynn (2024) Mauna Loa examples the conditions, orientations
and step model of the paper, and the Iovine et al. (2017) examples the published
law, temperature and initial step. Data-only entries identify the setup that
still needs to be supplied. The
[validation report](../examples/validation/report.md) compares published times,
reconstructed fits and alternative diffusion laws, and explains which checks
are exact age calculations, approximate profile reconstructions or coefficient
comparisons. Extracted measurements retain their source rows and workbook
hashes.

## Example datasets

One example per mineral is bundled and listed on the Data step with its
provenance. The list is grouped by mineral, then by study, then by example. It
starts closed, and a study with many profiles takes one row until it is opened.
Four studies supply published profiles from their supplementary tables: one
orthopyroxene traverse from Ostorero et al. (2022), 18 olivine traverses from
Lynn et al. (2024, Kilauea), 59 olivine traverses from Lynn et al. (2024, Mauna
Loa) and 20 sanidine Ba traverses from Iovine et al. (2017). The others were
generated with Diffusor's forward model at a set time and given Gaussian noise.
They are marked as synthetic in the catalogue, in the interface and in every
exported methods paragraph. The conditions of each example are taken from a
published study, and a value that the chosen law does not use is labelled as a
placeholder.

| dataset | kind | source | known time | fitted by Diffusor |
| --- | --- | --- | --- | --- |
| Orthopyroxene Fe-Mg, Kizimen 2010 | published | crystal K9_L10C4, Supplementary Data 2 of Ostorero et al. (2022) | 2.32 yr (+7.16/−1.75), published | 2.5 yr |
| Olivine Fe-Mg, Kilauea 2020, 18 crystals | published | Lynn et al. (2024), EarthChem doi:10.60520/IEDA/113087 | 14 to 521 d (± 30 %), published | 0.75 to 1.06 times the published time, all 18 within ± 30 % |
| Olivine Fe-Mg, Mauna Loa 2022, 59 crystals | published | Lynn et al. (2024, Nature Commun.), Supplementary Data 1 | 1.1 to 232 d, posterior medians with 95 % intervals | 0.81 to 1.34 times the published median, 58 of 59 inside the 95 % interval |
| Sanidine Ba, Agnano-Monte Spina, 20 greyscale and X-ray traverses | published | Iovine et al. (2017), Electronic Supplementary Material 1 to 3 | 1 to 151 yr at 930 °C, published (Tables 1 and 2) | 0.90 to 1.28 times the published time |
| Orthopyroxene Fe-Mg | synthetic | Dias et al. (2025), 950 °C, NNO+1, arc andesite conditions | 1.5 yr | 1.41 yr |
| Clinopyroxene Fe-Mg | synthetic | Mueller et al. (2013), 1100 °C, conditions after Petrone et al. (2016, 2018) | 45 d | 42 d |
| Olivine Fe-Mg | synthetic | Dohmen & Chakraborty (2007), 1150 °C, FMQ−1, conditions after Hartley et al. (2016) | 120 d | 116 d |
| Titanomagnetite Ti | synthetic | Van Orman & Crispin (2010), 950 °C, log fO2 = −11, the Shinmoedake conditions of Tomiya et al. (2013) | 8 d | 7.4 d |
| Clinopyroxene BSE grey values | synthetic | Mueller et al. (2013), 1000 °C, linear grey-value response with five microprobe anchor points | 3 yr | 3.0 yr |
| Sanidine Ba | synthetic | Cherniak (2002), 790 °C, within the Bishop Tuff range of Chamberlain et al. (2014) | 5 kyr | not tested |

The fitted values of the synthetic sets are those of the round-trip tests
(`tests/test_examples_roundtrip.py`). Each file is fitted for the time and the
interface position, with the plateaus taken from the outer points of the
traverse. The titanomagnetite traverse ends before the profile has flattened,
so plateaus taken from its outer points lie inside the true ones and shorten
the time to 5.5 d. Its plateaus are therefore fitted as well, and the example
sets this option when loaded. The remaining deviations lie within the scatter
produced by the added noise, one standard deviation of 3 to 13 % in the
fitted time.

The conditions of a synthetic set are not always those of the study named. The
titanomagnetite example evaluates D at the Shinmoedake conditions of Tomiya et
al. (2013), where Diffusor gives 4.15e-16 m²/s against their 4.3e-16 m²/s, but
its TiO<sub>2</sub> values (6.5 to 4.2 wt%) are illustrative. They are lower
than the 9.2 to 10.4 wt% of the type A magnetite of Shinmoedake
(X<sub>Usp</sub> 0.27 to 0.31) and correspond in proportion to X<sub>Usp</sub>
of about 0.12 to 0.19. D uses a fixed x<sub>Ti</sub> of 0.1, so this does not
change the recovered time. Tomiya et al. modelled a sphere, whereas the example
is a plane sheet. The orthopyroxene example at NNO+1 and 950 °C has
log fO2 = −5.0 (Pa), above the calibrated range of −11 to −7 of its law, so the
fO2 term is extrapolated in the file and in the fit. The sanidine example at
790 °C lies below the 828 to 1075 °C of the Ba measurements. The olivine example
takes its temperature (1150 ± 30 °C) and oxygen fugacity (FMQ −1 ± 0.5) from
Hartley et al. (2016, p. 61). The pressure, the traverse along [001], the
composition range and the 120 d duration are Diffusor choices.

The Kizimen conditions, 850 ± 57 °C and NNO +1.3 ± 0.35, are the mean and
standard deviation of 21 magnetite-ilmenite pairs from the andesites in
Supplementary Data 3 of Ostorero et al. (2022). The temperature statistics are
those printed in the table. The ΔNNO statistics (1.28 ± 0.35 before rounding)
are computed from its ΔNNO column, which has no mean. The main text of that
paper lists 832 ± 37 °C for the andesite and 850 ± 57 °C for the dacite, while
the table and the caption of its Fig. 3 give 850 ± 57 °C for the andesite.
Diffusor uses the table and caption values. Ostorero et al. modelled this
reverse zone at 850 °C with the Ganguly & Tazzoli (1994) law without an fO2
term. With the same law and temperature, Diffusor fits 2.5 years to the
microprobe traverse, within the published interval. The remaining difference is
expected, because Ostorero et al. fitted high-resolution BSE grey-value profiles
that are not included in the supplement, whereas the 2 µm microprobe traverse
resolves the 4 µm zone with few points. With the Dias et al. (2025) law instead,
the same traverse gives 12 years.

The Lynn et al. (2024) examples load the inputs archived with each profile in
the EarthChem workbook: 1200 ± 10 °C, 42 MPa, log fO2 = −8.2 (bar), the three
EBSD angles to a, b and c, and the tabulated initial profile, which Diffusor
interpolates linearly between the analysed positions. They use the paper's
law, Dohmen & Chakraborty (2007) with its pressure term (eq. 2 of the paper),
with D<sub>a</sub> = D<sub>b</sub> = D<sub>c</sub>/6 along the measured
orientation, fixed compositions at both ends and only the time free. The paper's
Methods give 45 MPa, and the workbook value is used. The sub-grid interface and
the numerical grid of the authors' MATLAB model are not published, so the refits
agree with the published times (ratios 0.75 to 1.06, and Ol 20 is 25 % shorter)
rather than reproduce them digit for digit. Ol 8 records two diffusion events
and is not offered as an example.
