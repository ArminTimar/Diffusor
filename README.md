# Diffusor

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23074989.svg)](https://doi.org/10.5281/zenodo.23074989)
[![tests](https://github.com/ArminTimar/Diffusor/actions/workflows/tests.yml/badge.svg)](https://github.com/ArminTimar/Diffusor/actions/workflows/tests.yml)

Diffusor is a desktop application and Python package for diffusion chronometry
of zoned crystals. It models one-dimensional diffusion profiles with the
closed-form solutions of Crank (1975) or with a Crank-Nicolson finite-difference
solver, fits the diffusion time to a measured compositional profile, and
propagates the uncertainties of temperature, oxygen fugacity, pressure, the
diffusion coefficient and the measurements to the fitted time by Monte Carlo
sampling. Oxygen fugacity is sampled as an offset from a mineral buffer and
recalculated at each sampled temperature, and the Arrhenius parameters are
sampled from a published covariance where one exists, so correlated inputs are
not sampled as if they were independent.

The program contains 113 published diffusion laws for 14 minerals: olivine,
orthopyroxene, clinopyroxene, plagioclase, K-feldspar, magnetite, garnet,
quartz, rutile, titanite, apatite, zircon, monazite and xenotime. Each entry
states whether it is a tracer, an interdiffusion (exchange), a chemical
(trace-element) or an effective coefficient. Each entry stores the
citation, the equation as printed in the source with its equation number, the
calibration ranges in temperature, pressure, oxygen fugacity and composition,
the published parameter uncertainties, and the table, page or figure from
which the values were transcribed. 112 of the 113 entries have been compared
with the rendered pages of their primary publications; the Fe-Ti law of Aragon
et al. (1984) is taken from Table 11 of Van Orman & Crispin (2010) and is marked
as not checked against the original paper. Each equation is shown typeset on
the coefficient page and exported as LaTeX source with the results. A run outside the calibration range of the selected
law produces a warning, and every export contains a methods paragraph that
names the diffusion law, solver, boundary conditions and uncertainty model of
that run, followed by the references it cites. The checks that are run against
published values are listed under [Verification](#verification).

## Installation on Windows

1. Install Python 3.10 or newer from [python.org](https://www.python.org/downloads/).
   Tick *Add python.exe to PATH* in the installer.
2. On this page, click **Code > Download ZIP**, then right-click the ZIP and
   choose **Extract All**. Diffusor does not run from inside the ZIP.
3. In the extracted folder, double-click **Start Diffusor.bat**. If Windows
   warns that the file came from the internet, choose **More info > Run
   anyway**.

The first start creates a separate Python environment in a `.venv` folder
inside the Diffusor folder and downloads the required libraries into it, so no
other Python installation on the computer is changed. This takes a few minutes
and needs an internet connection. Later starts open the application directly.
The first start also places a **Diffusor** shortcut with the application icon
next to the batch file; it can be moved to the desktop or the taskbar. If the
folder is moved, start once with **Start Diffusor.bat** to recreate the
shortcut. If the environment is damaged, delete the `.venv` folder and run
**Start Diffusor.bat** again.

## Installation on macOS and Linux

1. Install Python 3.10 or newer with the macOS installer from
   [python.org](https://www.python.org/downloads/). The `python3` supplied
   with macOS is too old.
2. On this page, click **Code > Download ZIP** and unpack it.
3. In the Diffusor folder, double-click **Start Diffusor.command**. It opens in
   Terminal.

macOS blocks unsigned files downloaded from the internet, so the first attempt
may report that the file cannot be opened or that Apple could not verify it.
It has to be allowed once:

* macOS 15 (Sequoia) and later: click **Done**, open **System Settings >
  Privacy & Security**, find the message about *Start Diffusor.command*, click
  **Open Anyway** and confirm.
* Earlier versions: right-click (or Control-click) the file, choose **Open**,
  then **Open** again.
* Any version: in Terminal, type `bash ` (with a trailing space), drag
  **Start Diffusor.command** into the window and press Return.

The first start sets up the `.venv` environment as on Windows. Later starts
open the application directly, and the Terminal window can be closed once the
application is running. To keep Diffusor in the Dock, right-click its icon
while it runs and choose **Options > Keep in Dock**.

On Linux, run `./Start\ Diffusor.command` from a terminal in the Diffusor
folder.

## Installation as a Python package

```bash
python -m pip install -e .
python -m diffusor
```

Dependencies are numpy, scipy, pandas, matplotlib, Pillow, openpyxl and
PySide6. `python -m pip install -e .[images]` adds `tifffile`, which reads
BigTIFF and compressed TIFF files that Pillow cannot open; the image extractor
works without it. The installation also creates a `diffusor` launcher
(`Scripts\diffusor.exe` on Windows) that opens the window without a console.

## Updates

At start-up, at most once a day, Diffusor asks `api.github.com` for the newest
release. Nothing else is sent, and nothing happens if the computer is offline.
When a newer version exists, a bar under the title offers **Download**, which
opens the release page in the browser, and **Later**, which hides the bar until
the next check. Diffusor does not download or replace any files itself. To
update, extract the new ZIP into a new folder and start it as described above;
the first start sets up the environment again. The automatic check is switched
off under **Help > Check for updates at startup**, and **Help > Check for
updates...** runs it manually.

## Workflow

The main window leads through six steps and then opens a results view.

1. **Data.** Load a profile file or one of the bundled examples. Diffusor
   proposes the distance column, its unit and the composition columns from the
   column names; the assignment is confirmed by the user. Uncertainty columns,
   oxide names and a fit window are set under *Advanced*. Files loaded before
   are listed under *Recent profiles*, which shares a panel with the profile
   plot; a file from that list whose columns have not changed opens with the
   mapping chosen last time, without the column dialog, and the file dialog
   opens in the folder used last. The spatial resolution of the analysis is
   set here as well: selecting an instrument fills in a beam width together
   with the published convention it is based on. Loading an example also
   fills in the following steps.
2. **Mineral.** Phase, diffusing species and orientation of the traverse.
3. **Coefficient.** All diffusion laws in the registry for the chosen mineral
   and species, marked as recommended or superseded. Several laws can be
   selected and fitted side by side. The law is chosen before the conditions
   because it determines which of them enter the calculation.
4. **Conditions.** Temperature with an optional linear cooling path, pressure
   and oxygen fugacity with their uncertainties, and the host composition
   (for example x<sub>Ti</sub> or X<sub>An</sub>) where the chosen law depends
   on it. Inputs that the law does not use are disabled; pressure remains
   active when oxygen fugacity is given relative to a buffer, because the
   buffer depends on pressure.
5. **Model.** Geometry (plane, cylinder or sphere), initial condition, and the
   boundary condition at each end of the profile (see
   [Numerical methods](#numerical-methods)). The panel states which solver
   will be used and why.
6. **Uncertainty.** Number of Monte Carlo draws, random seed, the sources to
   sample and the number of processor cores.

The results view shows the fitted time and its interval above a summary of all
settings. Each group in the summary links back to the step where it is set.
During a Monte Carlo run the plot is updated after every draw: the fitted
profile of each draw, the perturbed data when measurement noise is sampled,
the histogram of fitted times with the running median and 68 % interval, and
the sampled temperature, oxygen fugacity and log D. The methods text, the
details of the diffusion coefficient and the log are under the **View** menu.

## Input data

One row per analytical point, with a distance column, one or two composition
columns and optional uncertainty columns. Column names are free and are
assigned after loading. CSV, TSV, plain text and Excel files are read.

```
Distance_um,FeO_wt,MgO_wt,FeO_err,MgO_err
0.0,20.9,22.4,0.15,0.20
1.7,20.8,22.5,0.15,0.20
```

With two composition columns Diffusor forms the molar ratio A/(A+B), which is
the variable the Fe-Mg diffusion laws are calibrated against; if the oxides are
named, weight per cent is first converted to cation moles. A single column is
modelled as given, for example a trace element in ppm, forsterite in mol% or a
calibrated grey value. Olivine Fe-Mg profiles have an explicit selector for
X<sub>Fe</sub>, X<sub>Fo</sub> or Fo per cent. Distances may increase in either
direction and need not be evenly spaced. Additional columns are carried along,
so that an anorthite column can supply the activity term for plagioclase trace
elements. A fit window excludes points from the fit without removing them from
the file, for example a later overgrowth at the rim.

## Profiles from images

**File > Extract profile from image** reads profiles from back-scattered
electron images and element maps, following the `greyvalues.m` script of NIDIS
(Petrone et al. 2016). The user draws a guideline along the zone boundary.
Values are read along lines perpendicular to it, one line per pixel of
guideline and by default 50 px to either side, and averaged at each position.
As in NIDIS, the first value of each line lies to the right of the guideline as
seen from its first point. The guideline may have any number of vertices, so
curved boundaries can be followed.

* Formats: TIFF (8, 16 and 32 bit, floating point, multi-page stacks as
  channels), PNG, JPEG, BMP, GIF, WebP, ENVI `.hdr` with its data file, NumPy
  `.npy`, text grids of counts as exported by microprobe software, and
  headerless binary files of known width, height and data type. Values are not
  rescaled.
* Scale: read from Zeiss, Thermo Fisher (FEI), Tescan and ImageJ TIFF tags and
  from JEOL and Hitachi text sidecar files; otherwise entered by hand or
  measured on the scale bar.
* Values: grey value, luminance or mean of a colour image, a single channel or
  page, or a colour scale. For a false-colour element map the legend is
  sampled by clicking along it and entering its end values; the legend may be
  in the map or in a separate file. Known colour maps can be selected by name.
  Pixels whose colour is not on the legend (cracks, labels, epoxy) are left
  without a value.
* Filtering: lower and upper value limits for cracks, holes and bright
  inclusions, with optional growth of the rejected areas by a few pixels;
  polygons drawn around inclusions or lamellae; outlier rejection across the
  lines at each position (median ± k MAD by default, or mean ± 1 SD as in
  NIDIS); and rejection of lines with too many excluded values. Rejected values
  are marked on the image, and the averaged profile is recalculated as settings
  change.
* Output: an Excel workbook with the averaged profile (N, minimum, maximum,
  mean, median, SD, relative SD and SE, before and after filtering), every
  individual line, a code for the reason each value was rejected, the line
  geometry, the colour legend if used, an image of the lines on the picture,
  and all settings, so that the extraction can be repeated with *Reuse
  settings from a workbook*.
* Calibration: the workbook is loaded with **File > Load profile** or passed on
  with *Use in Diffusor*. Diffusor asks for the pixel size if the workbook has
  none and for the conversion from value to composition: a line through two
  reference points, a linear or quadratic fit to microprobe anchor points
  (entered by hand or taken from a microprobe traverse along the same line), or
  none if the values are already compositions. The uncertainty combines the
  scatter between lines with the uncertainty of the calibration.

The standard error of the mean over several hundred lines is small, and
neighbouring lines share pixels, so it underestimates the uncertainty of the
profile. The standard deviation across lines is offered as the alternative.
The folder `examples/images/` contains two synthetic test images: a 16-bit BSE
image of a clinopyroxene with a curved Fe-rich rim, a crack, an inclusion and a
lamella (`cpx_bse_zoned.tif`), and an MgO map in the jet colour scale with its
legend and a scale bar (`opx_mg_map_jet.png`).

## Shared-duration fits

**File > Shared-duration study** collects profiles configured in the main
window and fits one common duration to all of them, for example several
elements across the same zone boundary. Each profile keeps its own diffusion
law, geometry and boundary conditions, and each requires measurement
uncertainties, because the misfits of profiles in different units are weighted
by those uncertainties. Only the duration is a free parameter. The individual
residuals should be inspected before the profiles are interpreted as recording
a single event. The export contains a study workbook and a full report for each
profile.

## Multicomponent and isotope profiles

**File > Multicomponent and isotope study** has two tabs.

*Garnet (Fe-Mg-Mn-Ca).* A table with one column per cation is normalised to
Fe + Mg + Mn + Ca = 1 per row, and the four components are modelled together.
The diffusion matrix is that of an ideal ionic solution (Lasaga 1979, in the
form of Chakraborty & Ganguly 1992, eq. 2), recalculated at every node and
time step from the tracer coefficients of one of four published sets: Carlson
(2006), Chakraborty & Ganguly (1992), and each of them with the Mn law
recalibrated by Chen & Chu (2024). The off-diagonal terms can move a component
that has no initial gradient (uphill diffusion of Ca, Carlson 2006, Fig. 4).
One duration is fitted to all components by weighted least squares; the step
position can be fitted as well, and the Monte Carlo varies temperature and
each tracer coefficient (Carlson 2006 gives ±0.8 log units at 95 %).

*Isotopes.* The profile and model configured in the main window are combined
with the delta-value columns of the same table. For a dilute element (Li),
each isotope diffuses with D<sub>m</sub> = D (m<sub>ref</sub>/m)<sup>β</sup>.
For Fe-Mg in olivine, the three Mg and four Fe isotopes are seven components
of one exchange (Oeser et al. 2026, eqs 3-4), with the Fe and Mg tracer
coefficients of Oeser et al. (2026) or with the main-window Fe-Mg law and a
ratio D*<sub>Fe</sub>/D*<sub>Mg</sub>. The β values offered are those of
Oeser et al. (2026, Table 5; olivine Fe and Mg per crystal axis) and Richter
et al. (2014, 2017; Li in augite and olivine). β depends on the diffusion
model it was fitted with: the Li data of Richter et al. (2014) need β = 0.27
with a two-site model and 0.44 with a one-site model.

Tracer and exchange coefficients can differ by orders of magnitude. Na-K
interdiffusion in K-feldspar measured directly (Schäffer et al. 2014) is 5 to
10 times slower normal to (001), and almost 100 times slower normal to (010),
than the value calculated from Na and K tracer coefficients. A fit with a
tracer coefficient therefore carries a warning, and interdiffusion laws
computed from two tracer laws (olivine, Oeser et al. 2026; garnet, Borinski et
al. 2012) name the laws they were computed from.

## Exported results

**File > Export results** writes, for each run:

* the results as JSON, including the coefficient provenance, the uncertainty
  model and the Monte Carlo statistics;
* the measured profile, the initial condition, the best-fit model, the
  residuals and the 16th and 84th percentile envelope of the Monte Carlo fits
  as CSV;
* an Excel workbook with sheets for the summary, the profile with charts of the
  fit and the residuals, the diffusion coefficient, the methods text and the
  Monte Carlo draws;
* the fitted times of all Monte Carlo draws as CSV;
* the figure as PNG (300 dpi) and SVG;
* the methods paragraph as plain text, with the reference list of the run.

**Help > All references** lists every source used anywhere in the program, with
a search field, and saves the list as BibTeX.

## Example datasets

One example per mineral is bundled and listed on the Data step with its
provenance. One is measured data from a published supplementary table. The
others were generated with Diffusor's forward model at a set time and given
Gaussian noise; they are marked as synthetic in the catalogue, in the
interface and in every exported methods paragraph. The conditions of each
example are taken from a published study, and a value that the chosen law does
not use is labelled as a placeholder.

| dataset | kind | source | known time | fitted by Diffusor |
| --- | --- | --- | --- | --- |
| Orthopyroxene Fe-Mg, Kizimen 2010 | measured | crystal K9_L10C4, Supplementary Data 2 of Ostorero et al. (2022) | 2.32 yr (+7.16/−1.75), published | 2.6 yr |
| Orthopyroxene Fe-Mg | synthetic | Dias et al. (2025), 950 °C, NNO+1, arc andesite conditions | 1.5 yr | 1.39 yr |
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
the time to 5.3 d; its plateaus are therefore fitted as well, and the example
sets this option when loaded. The remaining deviations lie within the scatter
produced by the added noise, one standard deviation of 3 to 13 % in the
fitted time.

The Kizimen conditions, 850 ± 57 °C and NNO +1.3 ± 0.35, are the mean and
standard deviation of 21 magnetite-ilmenite pairs from the andesites in
Supplementary Data 3 of Ostorero et al. (2022). The main text of that paper
gives the andesite and dacite temperatures in reverse order; the table is used
here. Ostorero et al. modelled this reverse zone at 850 °C with the Ganguly &
Tazzoli (1994) law without an fO2 term. With the same law and temperature,
Diffusor fits 2.6 years to the microprobe traverse, within the published
interval. The remaining difference is expected, because Ostorero et al. fitted
high-resolution BSE grey-value profiles that are not included in the
supplement, whereas the 2 µm microprobe traverse resolves the 4 µm zone with
few points. With the Dias et al. (2025) law instead, the same traverse gives
13 years.

Earlier versions included a plagioclase traverse from Druitt et al. (2012)
(Santorini). It was removed because the example modelled it with a different
Mg diffusion law (Van Orman et al. 2014 instead of Costa et al. 2003), initial
profile and anorthite activity factor than the paper, so it neither reproduced
nor tested the published result.

## Numerical methods

Boundary conditions are selected by their physical meaning. A plateau that
continues corresponds to a semi-infinite medium, the assumption of the
error-function solution used by Ostorero et al. (2022) and Chamberlain et al.
(2014). A crystal rim in contact with the melt is a fixed concentration, the
open boundary of Costa et al. (2008, p. 555). A closed rim and the crystal
centre are zero-flux boundaries, the centre by symmetry (Crank 1975, section
4.3). As long as both plateaus are preserved, the far boundaries do not affect
the fitted time; Diffusor issues a warning when the diffusion front reaches the
end of the traverse.

Closed-form solutions from Crank (1975) are used where they are exact: the
semi-infinite step (eq. 2.14), the finite band (eq. 2.15), the plane sheet
(eq. 4.17), the cylinder (eq. 5.22), the sphere (eq. 6.18) and the
corresponding fractional-uptake series. Otherwise the diffusion equation is
solved with a Crank-Nicolson scheme on a uniform grid with the diffusion
coefficient evaluated at half-nodes (Crank 1975, section 8.4; Dohmen, Faak &
Blundy 2017). Plane, cylindrical and spherical geometries share one
conservative operator, with the L'Hôpital limit at the centre of a cylinder or
sphere. The time step starts at the explicit stability limit and increases
geometrically, which resolves the initial transient and keeps a 100 kyr run to
a few hundred steps.

The diffusion coefficient may depend on composition, in which case it is
evaluated at each node from the evolving profile. For plagioclase trace
elements the flux includes the anorthite activity term of Costa et al. (2003),
-D C (A/RT) dX<sub>An</sub>/dx. The factor A is taken from one of the two sets
in Table 1 of Dohmen, Faak & Blundy (2017), chosen on the Model step: Dohmen &
Blundy (2014), using their 900 or 1200 °C column, whichever is nearer the run
temperature, or Bindeman et al. (1998), the set used by Costa et al. (2003) and
Druitt et al. (2012). The two differ in sign for Mg. A/RT is evaluated at the
temperature of each time step, and every export states the set.
Anisotropy is handled with the direction-cosine relation of Costa &
Chakraborty (2004) and the axial ratios published for each mineral. A linear
cooling path replaces D t by the time integral of D(T(t)) (Crank 1975,
section 7.2; Lasaga 1983), with the buffer-referenced fO2 following the
temperature. The model profile is convolved with a Gaussian of the analytical
beam width before comparison with the data (Ganguly et al. 1988; Bradshaw &
Kent 2017), so that beam broadening is not interpreted as diffusion.

## Uncertainty propagation

The fitted time depends non-linearly on inputs that are correlated with each
other. Three correlations are treated explicitly.

1. Oxygen fugacity is usually known as an offset from a mineral buffer, and the
   buffer depends on temperature. Diffusor recalculates the buffer at every
   sampled temperature, so temperature and fO2 remain correlated.
2. ln D<sub>0</sub> and the activation energy Q of an Arrhenius law are
   positively correlated through the regression that produced them, and
   sampling them independently overestimates the spread of D at the working
   temperature. Where a source gives the correlation, the parameters are
   sampled jointly. Only the Sr and Ba laws of Grocolas et al. (2025) provide
   this: the authors treat log D<sub>0</sub> and Q as perfectly correlated in
   their own Monte Carlo, and Diffusor uses the same assumption. For 68 other
   laws log D is sampled at the working temperature with a stated scatter. For
   8 of them the scatter is published or derived from a statement in the source;
   for the other 60 it is Diffusor's assumption, because the source gives only
   errors of D<sub>0</sub> and Q without their covariance. The coefficient page
   and the methods paragraph state which applies. Laws with neither a covariance
   nor a scatter are held fixed by default (43 laws). Independent sampling of
   D<sub>0</sub> and Q is never the default and is available only to reproduce
   results obtained that way.
3. Composition enters both the diffusion coefficient and the fitted profile,
   so measurement noise affects the result through both.

Each draw perturbs the inputs and refits the profile. The reported values are
the median and the 2.5th, 16th, 84th and 97.5th percentiles of the fitted
times, not a symmetric standard deviation. The band on the profile plot is the
spread of the profiles fitted in the draws and is about as wide as the
measurement uncertainty. Temperature and the diffusion coefficient change the
time needed to reach the fitted diffusion length √(Dt), not the shape of the
fitted profile, so their contribution appears in the time histogram and the
interval rather than in the band. Correlations that are not supplied and
systematic errors that are not modelled are not included in the interval.
From Python, `diffusor.fitting.montecarlo.contributions` repeats the Monte
Carlo with one source active at a time, to rank the sources of uncertainty.

Draws are fitted in parallel, and a given seed gives the same draws on any
number of cores.

## Diffusion coefficient registry

| mineral | entries | species |
| --- | --- | --- |
| Olivine | 10 | Fe-Mg, Fe and Mg tracers, Ca, Be, P, Ni |
| Orthopyroxene | 12 | Fe-Mg, Mg, Lu, Ce, Eu |
| Clinopyroxene | 3 | Fe-Mg, Ca-Mg |
| Plagioclase | 11 | Mg, Sr, Ba, Li, NaSi-CaAl |
| K-feldspar | 7 | Sr, Ba, Ti, Na-K |
| Magnetite | 42 | Ti, Fe, Fe-Ti and 21 trace elements |
| Garnet | 3 | Fe-Mg, Fe and Mg tracers (multicomponent sets listed above) |
| Quartz | 1 | Ti |
| Rutile | 3 | Zr, Hf |
| Titanite | 2 | Sr, Zr |
| Apatite | 8 | Sr, Pb, rare earth elements |
| Zircon | 6 | Ti, Pb, rare earth elements |
| Monazite | 1 | Pb |
| Xenotime | 4 | Pb, rare earth elements |

The recommended law for each mineral and species, where one is designated:

| mineral | species | recommended source |
| --- | --- | --- |
| Olivine | Fe-Mg | Dohmen & Chakraborty (2007), TaMED mechanism, with the composition term of the erratum |
| Orthopyroxene | Fe-Mg | Dias, Dohmen & Behrens (2025), two fO2 regimes |
| Clinopyroxene | Fe-Mg | Mueller et al. (2013) |
| Plagioclase | Mg | Van Orman, Cherniak & Kita (2014); Audétat, Grocolas & Mutch (2026) with silica activity |
| Plagioclase | Sr, Ba | Grocolas, Bloch, Bouvier & Müntener (2025) |
| Plagioclase | Li | Pohl et al. (2024), interstitial mechanism |
| K-feldspar | Sr | Cherniak (1996), sanidine Or61 |
| K-feldspar | Ba | Cherniak (2002), same crystal as the Sr law |
| K-feldspar | Ti | Cherniak & Watson (2020) |
| Magnetite | Ti, Fe | Van Orman & Crispin (2010), Table 12 |

Older calibrations are kept so that published timescales can be reproduced and
compared, for example Ganguly & Tazzoli (1994) for orthopyroxene, also in the
form without an fO2 term used by Ostorero et al. (2022), Dimanov & Sautter
(2000) for clinopyroxene, which is the law behind published NIDIS results,
Costa et al. (2003) for Mg in plagioclase, and Freer & Hauptman (1978) and
Aragon et al. (1984) for Fe-Ti in titanomagnetite. **Help > All references**
lists every source.

Sievwright et al. (2020) measured 21 elements in magnetite at 1150 °C only.
Each row of their Table 5 is therefore an entry valid at that temperature, and
a Monte Carlo run that samples temperature is refused for these entries. Five
further entries (Ti, Mn, Co, Cr and Al) extrapolate these values to other
temperatures with activation energies from Van Orman & Crispin (2010); they
are labelled as Diffusor's own extrapolation, since Sievwright et al. did not
determine activation energies.

Sr and Ba in sanidine were measured on the same Or61 crystal with the same
method (Cherniak 1996, 2002). At 800 °C Ba diffuses about 1.7 log units more
slowly than Sr, so paired Sr and Ba profiles across one zone boundary indicate
whether the boundary has been modified by diffusion (Chamberlain et al. 2014;
Audétat, Grocolas & Mutch 2026). Both laws carry the log D uncertainties that
Chamberlain et al. (2014) derived for them, 0.03 log units for Sr and 0.12 for
Ba.

### Superseded laws

Four entries are marked as superseded. They are not recommended, are labelled
in the coefficient list, and produce a note in every fit that uses them.

Sr in plagioclase: Grocolas, Bloch, Bouvier & Müntener (2025) measured Sr
diffusion in oligoclase and labradorite between 900 and 1200 °C with buffered
silica activity and found it one to two orders of magnitude slower than
Giletti & Casserly (1994) and Cherniak & Watson (1994). They attribute the
faster older values to reaction fronts, because Sr-feldspar was not stable in
the earlier source materials. Because the new activation energy is higher, the
difference increases with decreasing temperature, from 1.6 log units at
1100 °C to 2.2 at 900 °C and 2.8 at 750 °C for An36. Timescales calculated with
the older laws at 750 to 900 °C are therefore about 100 to 600 times too short.
Their equations 7 and 8 are the recommended entries, and their refits of the
older data (equations 12 to 14) replace both older Sr laws. For Ba the new law
agrees with Cherniak (2002) within about half a log unit, and the older law is
kept as an alternative.

Fe-Mg in orthopyroxene: Dias & Dohmen (2024) and Dias, Dohmen & Behrens (2025)
refitted the experiments of Dohmen et al. (2016) with a temperature-dependent
composition term, added Fe-rich crystals, and found two regimes separated at
log fO2 = −10 (fO2 in Pa). For Fs10 at log fO2 = −7 the 2025 law agrees with
Dohmen et al. (2016) within 0.3 log units, but the two differ by up to one log
unit for Fe-rich orthopyroxene at 1100 °C or under reducing conditions.
The general Dohmen et al. (2016) law and the intermediate law of Dias &
Dohmen (2024) are therefore marked as superseded. Arc magmas generally lie above the calibrated
fO2 range, where the fO2<sup>1/4</sup> term is extrapolated.

Fe-Ti in magnetite has no direct successor to Freer & Hauptman (1978) or Aragon
et al. (1984), and these entries are not marked as superseded. For Ti and Fe
tracer diffusion, however, the data of Aggarwal & Dieckmann (2002) as
tabulated by Van Orman & Crispin (2010, Table 12) depend explicitly on both
temperature and fO2, are better constrained, and were used by Tomiya et al.
(2013) at Shinmoedake. They are the recommended entries. At 1150 °C they agree
with Sievwright et al. (2020) within half a log unit for Ti, Mn and Co.

## Verification

The test suite contains 293 tests and runs on Windows, macOS and Linux
with Python 3.10 and 3.13 on every change to the repository (see the badge at
the top). It does not access the network. The tests cover the following.

Citation integrity. Every coefficient must cite an existing reference key for
its primary and secondary sources, contain the transcribed equation, and state
where in the source its values were read. Every measured example dataset and
every superseded entry must cite an existing key, every reference must render
as text and as BibTeX, and an unknown key raises an error.

Agreement with published values. Coefficients are evaluated at the conditions
of published experiments or worked examples and compared with the printed
values:

| check | published | Diffusor |
| --- | --- | --- |
| Opx Fe-Mg at 950 °C, log fO2 = −7 (Pa), parallel to c (Dohmen et al. 2016, run OPXD_14) | log D = −19.49 ± 0.07 | −19.47 |
| Opx anisotropy D[001]/D[100] (Dohmen et al. 2016) | 3.5 | 3.5 |
| Cpx Fe-Mg at 1098 and 1150 °C (Petrone et al. 2016, Table 2) | 3.26e-20 and 1.20e-19 m²/s | 3.26e-20 and 1.20e-19 m²/s |
| Titanomagnetite Ti at 950 °C, log fO2 = −11, X<sub>Usp</sub> = 0.3 (Tomiya et al. 2013) | 4.3e-16 m²/s | 4.35e-16 m²/s |
| Titanomagnetite Ti at 900 °C, same conditions | 6.9e-16 m²/s | 6.84e-16 m²/s |
| Opx Fe-Mg for Fs9 at 950 to 1100 °C, log fO2 = −7 (Pa) (Dias & Dohmen 2024, Table 1) | fitted D<sub>0</sub> and m of each run | within 0.2 log units, for the 2024 and 2025 laws |
| Olivine Ni at 1005 °C, Fo90 (Petry et al. 2004, experiment Ni10) | log D = −17.12 | within 0.16 log units |
| Plagioclase Sr, difference from Giletti & Casserly (1994) at An36 and 750 °C (Audétat et al. 2026) | 2.8 log units | 2.77 |
| Magnetite, minimum of D with fO2 for 21 elements at 1150 °C (Sievwright et al. 2020, Table 5) | log fO2 and log D at the minimum | within 0.1 and 0.05 log units |
| Opx K9_L10C4, Kizimen (Ostorero et al. 2022, Supplementary Data 4) | 2.32 yr (+7.16/−1.75) | 2.6 yr |
| Garnet D matrix at X<sub>Fe</sub> 0.61, X<sub>Mn</sub> 0.20, X<sub>Mg</sub> 0.18, X<sub>Ca</sub> 0.01 (Chakraborty & Ganguly 1992, eq. 5) | nine printed elements | within 0.4 % |
| Olivine D<sub>Fe-Mg</sub> from the Fe and Mg tracer fits, 1100 to 1250 °C, a, b and c (Oeser et al. 2026, eq. 6 and Table 4) | the paper's own D<sub>Fe-Mg</sub> fit | within 0.2 log units |
| Garnet Mn tracer at 510 °C, 2 GPa, eclogite garnet (Chen & Chu 2024) | about 1e-24 m²/s | 2.2e-24 m²/s |

Further tests reproduce the Arrhenius parameters, equation forms and figures of
the other sources, for example the olivine composition term of the Dohmen &
Chakraborty (2007) erratum and the Fig. 3 line of Grove et al. (1984). The
magnetite Fe entry agrees with Tomiya et al. (2013) within 11 % at 900 °C but
is 2.1 times higher at 950 °C, because Tomiya et al. place the diffusion minimum
at a different temperature. This difference is documented in the entry.

Solver accuracy. The finite-difference solver is compared with the closed-form
solutions of Crank (1975) for the semi-infinite step, the plane sheet, the
cylinder and the sphere. On an 801-node grid the largest deviation is 3 to 4
parts per million of the concentration step; the tests require less than
1e-4. The error decreases by a factor of about four with each halving of the
node spacing (second-order convergence). A closed plane system conserves mass
to 1e-8 and a closed sphere to 1e-6, and a closed system relaxes to its mean
composition. The explicit and implicit schemes agree with each other, the
fractional-uptake series agree with the integrated profiles, and a constant
temperature history gives the same result as the isothermal solution. The
multicomponent solver agrees with the eigen-component solution for a constant
matrix (Toor 1964, as used by Chakraborty & Ganguly 1992) within 1e-3 in mole
fraction on a 401-node grid, conserves each component to 1e-9, reduces to the
scalar solver for two components, and gives the same profiles whichever
component is taken as dependent. With β = 0 the isotope models give zero delta
values and the element profiles of the scalar model.

Fitting and Monte Carlo. Fits recover a known time from noise-free and noisy
synthetic profiles, with and without composition-dependent D. The bundled
synthetic examples are fitted from file and the recovered time is compared
with the known one (table above). Further tests check that ignoring beam
convolution lengthens the apparent time, that a zero uncertainty budget gives
no spread, that the spread grows with the temperature uncertainty, that fO2
follows the buffer when temperature is sampled and along a cooling path, that
independent sampling of ln D<sub>0</sub> and Q gives a wider interval than
sampling log D at temperature, and that parallel runs reproduce the
single-core draws exactly.

Input, output and interface. Unit and oxide conversions, buffer equations,
grey-value calibration, image reading and profile extraction, the contents of
the exported JSON and Excel files, and the main interface workflows (built off
screen) are also tested.

The tests are run with:

```bash
python -m pip install -e .[dev,images]
python -m pytest
```

## Limitations

* The model is one-dimensional. A time obtained from a 1-D profile through a
  3-D crystal is subject to geometric and sectioning effects (Shea et al. 2015;
  Krimer & Costa 2017); every fit states this.
* A Monte Carlo run with a composition-dependent law and the numerical solver
  takes a few seconds per draw and core, so several hundred draws take
  minutes. The plot shows the progress, and **Stop** ends the run.
* The garnet D matrix assumes an ideal solution. Borinski et al. (2012) found
  that non-ideality changes retrieved coefficients by less than a factor of
  1.2 for most natural garnets; their non-ideal matrix is not implemented.
* The fO2 term of Carlson (2006) is relative to the graphite-oxygen
  equilibrium, which Diffusor does not calculate; the offset from it is an
  input (0 by default).
* Isotopes of a dilute element are modelled with one diffusing species. The
  two-site Li models of Richter et al. (2014, 2017) are not implemented, so a
  β taken from them belongs to a different model.
* Reactions between sites, growth and resorption, and Bayesian inference are
  not implemented.
* The Chakraborty (1997) olivine law is a single fit for Fo86 along [001] at
  fO2 = 1e-12 bar. Other directions are refused, and other fO2 values produce
  a warning.
* Li in plagioclase is modelled with one effective coefficient per mechanism.
  The multispecies model of Pohl et al. (2024), with exchange between
  interstitial and lattice-site Li, is not implemented.
* Ni in olivine is restricted to the Fo90, fixed-fO2 fit of Petry et al.
  (2004); Mn in olivine is not included.
* The Lu law of Dias et al. (2025) gives its fO2 dependence relative to a
  reference fO2<sub>0</sub> that the paper does not define. Diffusor uses the
  IW buffer, because the paper describes the laws as valid along it; the entry
  states that this is an interpretation.
* The plagioclase Mg law of Audétat et al. (2026) is taken from the accepted
  manuscript. Its silica activity is 1 (quartz saturation) unless changed from
  Python; the interface does not offer this setting.

Corrections to results of earlier versions are listed in
[CHANGELOG.md](CHANGELOG.md).

## Citation

If Diffusor contributes to published work, cite it with the version used; the
version is stated in every exported methods paragraph. Each release is archived
on Zenodo. [doi:10.5281/zenodo.23074989](https://doi.org/10.5281/zenodo.23074989)
resolves to the newest release, and each version's Zenodo page gives its own
DOI. GitHub's **Cite this repository** button, generated from
[CITATION.cff](CITATION.cff), gives the citation in APA and BibTeX format. The
diffusion coefficients should be cited as well; the methods paragraph lists the
sources of each run, and **Help > All references** exports them as BibTeX.

## Licence

MIT, see [LICENSE](LICENSE).
