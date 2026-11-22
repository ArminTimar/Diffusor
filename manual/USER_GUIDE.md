# Diffusor user guide

This guide describes the windows and files of Diffusor in detail. The
[README](../README.md) covers installation and the first run. The models and
the diffusion laws are described in [METHODS.md](METHODS.md), and the checks
against published values in [VERIFICATION.md](VERIFICATION.md).

## Installation in detail

The installers contain Diffusor together with its own copy of Python and of
every library it uses. No Python installation, terminal or internet connection
is needed to install or run it, and no other software on the computer is
changed.

### Windows

Double-click the downloaded `Setup.exe`. The setup wizard asks for the licence,
the folder and whether to put a shortcut on the desktop. By default it installs
for the current user only, into `%LOCALAPPDATA%\Programs\Diffusor`, and does not
need administrator rights. The first page also offers an installation for all
users, which does. Diffusor is then started from the Start Menu, or from the
last page of the wizard.

The installer is not yet signed with a code-signing certificate, so Windows
SmartScreen shows *Windows protected your PC* the first time: click **More
info**, then **Run anyway**. On computers where Smart App Control is switched
on, Windows does not run unsigned programs at all.

Diffusor is removed under **Settings > Apps > Installed apps**.

### macOS

Double-click the downloaded `.dmg` file. A window opens with the Diffusor icon
and the Applications folder. Drag **Diffusor** onto **Applications**, then start
it from Launchpad or the Applications folder. To keep it in the Dock,
right-click its icon while it runs and choose **Options > Keep in Dock**.

The **Apple menu > About This Mac** window shows which processor a Mac has:
*Chip: Apple M…* for Apple silicon, *Processor: … Intel …* for Intel.

The application is not yet notarised by Apple, so macOS refuses the first start
with a message that Apple could not verify it. It has to be allowed once:

* macOS 15 (Sequoia) and later: click **Done**, open **System Settings >
  Privacy & Security**, find the message about *Diffusor*, click **Open
  Anyway** and confirm.
* Earlier versions: right-click (or Control-click) Diffusor in Applications,
  choose **Open**, then **Open** again.

Diffusor is removed by dragging it from Applications to the Bin.

### Linux

* **AppImage** (any distribution): right-click the downloaded file, choose
  **Properties > Permissions** and tick *Allow executing file as program*,
  then double-click it. The AppImage is the program itself and can be kept in
  any folder.
* **.deb** (Ubuntu, Debian, Linux Mint): double-click the file. The software
  centre installs it after asking for the user's password, and Diffusor then
  appears among the applications.

## Updates

At start-up, at most once a day, Diffusor asks `api.github.com` for the newest
release. Nothing else is sent, and nothing happens if the computer is offline.
When a newer version exists, a bar under the title offers **Download**, which
opens the release page in the browser, and **Later**, which hides the bar until
the next check. Diffusor does not download or replace any files itself. To
update, download the new installer and install it as described above: on
Windows the setup wizard replaces the installed version, on macOS the new
Diffusor is dragged onto Applications and replaces the old one, and on Linux
the new AppImage or `.deb` takes the place of the old one. Settings and the
list of recent files are kept. The automatic check is switched off under
**Help > Check for updates at startup**, and **Help > Check for updates...**
runs it manually.

## Running from the source code

This is for development, or for systems without an installer. It needs Python
3.10 or newer from [python.org](https://www.python.org/downloads/). Download the
source with **Code > Download ZIP** and unpack it, then double-click
**Start Diffusor.bat** (Windows) or **Start Diffusor.command** (macOS, and on Linux
run `./Start\ Diffusor.command` in a terminal). The first start creates a
separate Python environment in a `.venv` folder inside the Diffusor folder and
downloads the required libraries into it, which takes a few minutes and needs
an internet connection. Later starts open the application directly. If the
environment is damaged, delete the `.venv` folder and start again.

### As a Python package

```bash
python -m pip install -e .
python -m diffusor
```

Dependencies are numpy, scipy, pandas, matplotlib, Pillow, openpyxl and
PySide6. `python -m pip install -e .[images]` adds `tifffile`, which reads
BigTIFF and compressed TIFF files that Pillow cannot open. The image extractor
works without it. The installation also creates a `diffusor` launcher
(`Scripts\diffusor.exe` on Windows) that opens the window without a console.

### Building the installers

```bash
python -m pip install -e .[build,images]
python packaging/build.py
```

`packaging/build.py` makes the installer for the system it runs on, in `dist/`.
PyInstaller first collects Python, the libraries and Diffusor into one folder
(recipe in `packaging/diffusor.spec`). The bundled program is then started with
`--self-check`, which builds the window off screen and runs a Monte Carlo fit on
two worker processes. Finally the installer is made: Inno Setup 6 on Windows
(`packaging/windows/diffusor.iss`), a disk image on macOS, and an AppImage and
a `.deb` package on Linux (`appimagetool` and `dpkg-deb`). When a release is
published, the `installers` workflow runs these steps on Windows, macOS (Apple
silicon and Intel) and Ubuntu 22.04 and attaches the five files to the release.

## Code signing policy

The Windows installer is to be signed through the free code-signing programme
of the [SignPath Foundation](https://signpath.org) for open-source projects. The
application is pending, and until it is accepted the installer is unsigned.

* **Signed files.** `Diffusor.exe` and the setup program
  `Diffusor-<version>-Windows-Setup.exe`, built by the `installers` workflow
  from the source code in this repository. The libraries in the installation
  (Python, Qt, numpy, scipy and the others) are included as their authors
  publish them and are not signed again.
* **Roles.** Committer and reviewer: [Armin Timar](https://github.com/ArminTimar).
  Every release is approved for signing by
  [Armin Timar](https://github.com/ArminTimar).
* **Privacy.** Diffusor sends one request on its own: at start-up, at most once
  a day, it asks `api.github.com` for the version number of the newest release
  (see [Updates](#updates)). The request carries Diffusor's own version number
  and no data about the user, the computer or the files opened. The check is
  switched off under **Help > Check for updates at startup**. Diffusor transfers
  nothing else to other networked systems unless the user requests it.

## The steps

The main window leads through six steps to a seventh, the results view. The
numbered steps in the header can be clicked to move between them.

**Data.** Load a profile file or one of the bundled examples. Diffusor proposes
the distance column, its unit and the composition columns from the column
names, and the user confirms the assignment. Uncertainty columns, oxide names
and a fit window are set under *Advanced*. Files loaded before are listed under
*Recent profiles*, which shares a panel with the profile plot. A file from that
list whose columns have not changed opens with the mapping chosen last time,
without the column dialog, and the file dialog opens in the folder used last.
The spatial resolution of the analysis is set here as well: selecting an
instrument fills in a beam width together with the published convention it is
based on. Loading an example also fills in the following steps.

**Mineral.** Phase, diffusing species and orientation of the traverse. The line
under the two boxes says what the profile of the chosen species is: the mineral's
composition variable itself for an exchange such as Fe-Mg in orthopyroxene
(X<sub>Fe</sub>), or the measured content of a trace or tracer element.

**Coefficient.** All diffusion laws in the registry for the chosen mineral and
species, marked as recommended or superseded. Several laws can be selected and
fitted side by side. The selected law is shown with its equation and its
calibration range in °C. The law is chosen before the conditions because it
determines which of them enter the calculation.

**Conditions.** Temperature with an optional linear cooling path, pressure and
oxygen fugacity with their uncertainties, and the host composition (for example
x<sub>Ti</sub> or X<sub>An</sub>) where the chosen law depends on it. Inputs
that the law does not use are disabled. Pressure remains active when oxygen
fugacity is given relative to a buffer, because the buffer depends on pressure.

**Model.** Geometry (plane, cylinder or sphere), initial condition, and the
boundary condition at each end of the profile (see
[METHODS.md](METHODS.md#numerical-methods)). The panel states which solver will
be used and why.

**Uncertainty.** Number of Monte Carlo draws, random seed, the sources to sample
and the number of processor cores.

The results view shows the fitted time and its interval above a summary of all
settings. Each group in the summary links back to the step where it is set.
During a Monte Carlo run the plot is updated after every draw: the fitted
profile of each draw, the perturbed data when measurement noise is sampled,
the histogram of fitted times with the running median and 68 % interval, and
the sampled temperature, oxygen fugacity and log D. While the run goes on,
**Details** shows the statistics of the draws finished so far and **Methods**
the methods text for the best fit, which the run makes first. The methods text,
the details of the diffusion coefficient and the log are also under the
**View** menu.

**Edit graph** in the toolbar above the plot opens a column on the left of the
plot, in place of the summary, with everything that can be changed: the title,
the legend and its position, the text size and grid lines; every axis with its
label, range (automatic or typed) and linear or logarithmic scale; and every
curve with its legend text, colour, opacity, line style, width, marker and
marker size. *Done*, or the toolbar button again, closes it and brings the
summary back. Changes show on the plot at once and are
kept when the plot is drawn again after the next fit. Typed axis ranges are
dropped when another profile is loaded, and *Reset all* returns to the plot as
Diffusor draws it. Figures for publication are saved with the toolbar's save
button or **File > Export figure for Inkscape or CorelDRAW**.

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
the variable the Fe-Mg diffusion laws are calibrated against. If the oxides are
named, weight per cent is first converted to cation moles. A single column is
modelled as given, for example a trace element in ppm, forsterite in mol% or a
calibrated grey value. Olivine Fe-Mg profiles have an explicit selector for
X<sub>Fe</sub>, X<sub>Fo</sub> or Fo per cent. Distances may increase in either
direction and need not be evenly spaced. Additional columns are carried along,
so that an anorthite column can supply the activity term for plagioclase trace
elements. A fit window excludes points from the fit without removing them from
the file, for example a later overgrowth at the rim. Single points, or a group
inside a box, can also be cut on the plot of the Results step (**Cut points**
in its toolbar). A spurious analysis or the stretch across a second zone
boundary then stays out of the fit, and out of the Monte Carlo, while it
remains visible as a grey cross. At least three points must stay.

## Profiles from images

**File > Extract profile from image** reads profiles from back-scattered
electron images and element maps, following the `greyvalues.m` script of NIDIS
(Petrone et al. 2016). With the **Boundary** tool the user draws a straight
line along the zone boundary: one click where it starts and one where it ends.
After the first click a dashed line follows the mouse; a further click starts a
new line and a right-click removes it. Values are read along lines
perpendicular to it, one line per pixel of boundary line and by default 50 px
to either side, and averaged at each position. As in NIDIS, the first value of
each line lies to the right of the boundary line as seen from its first point
(*Start on the left* reverses this). For a curved boundary, draw the line along
the stretch being profiled. A bent guideline of several points saved in a
workbook by an earlier version is still read with *Reuse settings from a
workbook*.

A line under the drawing tools says what the next click does and, when no
profile appears, why: a second end still to click, a colour legend not yet set,
or an error. Choosing a drawing tool switches off the toolbar's pan and zoom,
which otherwise take the mouse clicks.

* Formats: TIFF (8, 16 and 32 bit, floating point, multi-page stacks as
  channels), PNG, JPEG, BMP, GIF, WebP, ENVI `.hdr` with its data file, NumPy
  `.npy`, text grids of counts as exported by microprobe software, and
  headerless binary files of known width, height and data type. Values are not
  rescaled.
* Scale: read from Zeiss, Thermo Fisher (FEI), Tescan and ImageJ TIFF tags and
  from JEOL and Hitachi text sidecar files, otherwise entered by hand or
  measured on the scale bar.
* Values: grey value, luminance or mean of a colour image, a single channel or
  page, or a colour scale. For a false-colour element map the legend is
  sampled by clicking along it and entering its end values. The legend may be
  in the map or in a separate file. Known colour maps can be selected by name.
  Pixels whose colour is not on the legend (cracks, labels, epoxy) are left
  without a value.
* Filtering: lower and upper value limits for cracks, holes and bright
  inclusions, with optional growth of the rejected areas by a few pixels,
  polygons drawn around inclusions or lamellae (*Exclude area*: click the
  corners, then the first corner again, or right-click, to close the area),
  outlier rejection across the
  lines at each position (median ± k MAD by default, or mean ± 1 SD as in
  NIDIS), and rejection of lines with too many excluded values. Rejected values
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
by those uncertainties. Only the duration is a free parameter. Check that the
model follows each profile (the export lists the standardised residuals of
each) before the profiles are interpreted as recording a single event. *Fit shared duration* becomes available once two profiles have
been added. The export contains a study workbook and a full report for each
profile.

## Crystal setups and batches

**File > Save crystal setup** writes one crystal to a `.diffusor` file: the
profile table as it was loaded, the column mapping, the points cut from the fit
on the plot, every setting of the steps and, when it is current, the fitted time
and the Monte Carlo times. **File > Open crystal setup** restores the data and
all settings and fits again, so the plot and the numbers on screen belong to the
restored settings. The file is self-contained and does not need the original
data file.

**File > Batch of crystals** fits many profiles that share most of their
settings. The window is a table with one row per crystal and one column per
setting (diffusion law, temperature, pressure, oxygen fugacity, host
composition, traverse direction and angles, boundaries, step position, plateaus,
fitted parameters, beam width and others, and *Columns* chooses which are
shown). The first row, *All crystals*, holds the settings of every crystal. A
value entered there applies to all crystals, and Diffusor asks whether crystals
that have their own value should keep it. A value entered in a crystal's row
belongs to that crystal only and is shown in bold. The context menu returns it
to the shared value or copies it to all crystals.

A batch starts from the settings of the main window. Profiles added with *Add
profiles* share one column mapping (the mapping dialog opens again for a file
whose columns differ), and each takes its step position and plateaus from its
own data. Crystals can also be added from saved setups, from another batch, or
from the main window. *Open in main window* loads one crystal into the steps to
inspect its fit or move its step, and *Take back from main window* returns the
changed settings to its row. *Fit* and *Fit and Monte Carlo* run the selected
crystals, or all of them, one after another, using the same model as the main
window. The table shows the status, best-fit time, Monte Carlo 68 % interval and
reduced χ² of each, and a crystal whose fit produced warnings is marked for
checking. *Send to population* passes the results to the timescale population
window with their groups. A batch is saved as one `.diffusor` file holding the
shared settings and every crystal's data, own values and result.

## Timescale populations

**File > Timescale population** compares the timescales of many crystals in one
figure. Timescales are added from the run in the main window, from a batch of
crystals, from saved crystal setups and batches (`.diffusor`), from results
files saved with *Export results* (`*_results.json`, with the Monte Carlo times
read from the `*_montecarlo_times.csv` beside it when that file was exported),
or from a CSV or Excel table with one row per crystal. Each timescale belongs to
a group, such as an eruption, a sample or a mineral, and up to eight groups are
drawn with their own colour and marker. The figure has three panels on one
logarithmic time axis, each of which can be switched off:

* **Ranked timescales**: the best-fit time of each crystal with its 68 % and
  95 % intervals, sorted within each group.
* **Summed probability density**: for each group, the mean of its crystals'
  probability densities in log<sub>10</sub> time, each normalised to one, so
  every crystal has the same weight. A crystal with Monte Carlo times
  contributes a Gaussian kernel density estimate of log<sub>10</sub> t with the
  rule-of-thumb bandwidth of Silverman (1986). A crystal with only a 68 %
  interval contributes two half-normals joined at the median, each holding
  half of the probability, so the 16th and 84th percentiles fall at the given
  bounds. A crystal without an uncertainty appears only in the other two
  panels.
* **Cumulative distribution**: the empirical distribution of the best-fit
  times of each group.

A table is read by its column names, ignoring case: a label (`label`,
`crystal`, `name` or `id`), an optional group (`group`, `eruption`,
`population` or `event`), the time (`time`, `t`, `best` or `median`), and
optionally a 68 % interval (`p16` and `p84`, or `lower` and `upper`), a 95 %
interval (`p2.5` and `p97.5`) or a standard deviation in log<sub>10</sub> units
(`sigma_log10`). The unit of a column is given in its name, for example
`time (days)` or `p16_yr`, or in a `unit` column. Otherwise it is asked for when
the table is loaded. The export contains the timescales and group statistics
as CSV, the density curves, a methods text and the figure as PNG, SVG and PDF.

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
One duration is fitted to all components by weighted least squares. The step
position can be fitted as well, and the Monte Carlo varies temperature and
each tracer coefficient (Carlson 2006 gives ±0.8 log units at 95 %).

*Isotopes.* The profile and model configured in the main window are combined
with the delta-value columns of the same table. For a dilute element (Li),
each isotope diffuses with D<sub>m</sub> = D (m<sub>ref</sub>/m)<sup>β</sup>.
For Fe-Mg in olivine, the three Mg and four Fe isotopes are seven components
of one exchange (Oeser et al. 2026, eqs 3-4), with the Fe and Mg tracer
coefficients of Oeser et al. (2026) or with the main-window Fe-Mg law and a
ratio D*<sub>Fe</sub>/D*<sub>Mg</sub>. The β values offered are those of
Oeser et al. (2026, Table 5, olivine Fe and Mg per crystal axis) and Richter
et al. (2014, 2017, Li in augite and olivine). β depends on the diffusion
model it was fitted with. The Li data of Richter et al. (2014) need β = 0.27
with a two-site model and 0.44 with a one-site model.

Tracer and exchange coefficients can differ by orders of magnitude. Na-K
interdiffusion in K-feldspar measured directly (Schäffer et al. 2014) is 5 to
10 times slower normal to (001), and almost 100 times slower normal to (010),
than the value calculated from Na and K tracer coefficients. A fit with a
tracer coefficient therefore carries a warning, and interdiffusion laws
computed from two tracer laws (olivine, Oeser et al. 2026, and garnet, Borinski
et al. 2012) name the laws they were computed from.

## Exported results

**File > Export results** first asks which outputs you want, with a tick box for
each, and then for the folder. What you switch off stays off the next time. It
can write, for each run:

* the results as JSON, including the coefficient provenance, the uncertainty
  model and the Monte Carlo statistics
* the measured profile, the initial condition, the best-fit model, the
  residuals and the 16th and 84th percentile envelope of the Monte Carlo fits
  as CSV
* an Excel workbook with sheets for the summary, the profile with charts of the
  fit and the residuals, the diffusion coefficient, the methods text and the
  Monte Carlo draws
* the fitted times of all Monte Carlo draws as CSV
* the figure as PNG (300 dpi) and, for editing, as SVG and PDF (see below)
* the methods paragraph as plain text, with the reference list of the run

**File > Export figure for Inkscape or CorelDRAW** saves the plot on screen as
SVG, PDF or EPS, and the save button of the plot toolbar does the same for these
three types. The files are made for editing afterwards. Text stays text, in Arial
where it is installed. Each data point is its own group, named for its series
and number (for example `measured_point_007`), holding the marker and its error
bar, so one point can be selected, recoloured, enlarged or removed. A marker is
a plain path, not a copy of a shared shape. The fit, the Monte Carlo bands, the
legend and each axis have names as well. In Inkscape, **Edit > Find** with the
name of a point selects it. CorelDRAW 2018 and later opens the SVG directly.
The PDF is the alternative where fonts matter, because it embeds them.

**Help > All references** lists every source used anywhere in the program, with
a search field, and saves the list as BibTeX.
