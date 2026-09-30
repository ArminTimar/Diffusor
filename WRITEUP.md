# Diffusor — code and sources writeup

> **Current implementation:** section 5 and [LITERATURE_AUDIT.md](LITERATURE_AUDIT.md) supersede conflicting baseline descriptions below. Sections 1–4 retain the earlier audit as historical context.

A single reference that answers three questions for every part of Diffusor:

1. **What does this script do, and when does the application call it?**
2. **Which variables does it use, and where do they come from?**
3. **Which paper is behind every equation, coefficient, constant, convention
   and standard in it?**

and, from the other direction, for every paper in the registry:

4. **What is this paper about?**
5. **Exactly what does Diffusor take from it?**

The point is to make the app improvable without re-reading the whole tree, and
to make a disagreement between the code and a source something you can find by
looking rather than by accident. Section 4 collects the disagreements that are
already known.

| | |
| --- | --- |
| Covers | September 29 working tree; section 5 updates the September 18 baseline |
| Package version | 0.1.0 |
| Modules documented | 70 Python modules |
| Coefficients | 103 entries, 13 minerals, 98 with source-transcription checks |
| References | 95 keys in `diffusor/references.py` |

---

## 0. Keeping this document current

This file is **not** generated. `REFERENCES.md` and `references.bib` are
generated from `diffusor/references.py` by `scripts/build_references.py`; this
file is written by hand because the parts that matter — what a module is for,
when it runs, what is doubtful about a number — are not in the code as data.

`tests/test_writeup.py` fails when the document falls behind the code. It
checks that every Python module, every coefficient key, every citation key,
every example dataset key and every oxygen buffer named in the code also
appears here. It cannot check that the prose is still true.

**When you change Diffusor, update this file in the same commit:**

| change | what to update here |
| --- | --- |
| new or renamed module | section 2, and the package map in section 1 |
| new diffusion coefficient | section 2.5 table, section 3.5, section 4.1 and 4.4 |
| new citation key | section 3 (overview + what is used from it) |
| changed equation or constant | section 4.1 or 4.2, and the module entry |
| a coefficient verified against its paper | section 4.4, and drop it from 4.3 |
| a discrepancy found or resolved | section 4.3 |
| new GUI step or worker | section 2.10 |
| new example dataset | section 2.9 |

Run `python -m pytest tests/test_writeup.py` to see what is missing.

---

## 1. Map of the package

Diffusor is a layered library with a desktop front end. Nothing in the lower
layers imports from a higher one.

```
                        gui/            PySide6 window, six steps, workers, plots
                          |
   datasets.py  <---------+------------>  dataio/     load, calibrate, export
                          |
                       fitting/          model -> fit -> Monte Carlo
                       /      \
                solvers/       coefficients/     Crank solutions, Crank-Nicolson
                    |               |            registry of published laws
                minerals/       thermo/          phases, axes, X variables
                    \               /            buffers, units, conversions
                     constants.py, references.py
```

**A run, from click to number.** The order below is the order of calls for a
single fit, which is the spine of everything else:

1. `gui/main_window.MainWindow.load_file` → `dataio.profiles.read_table` →
   `suggest_spec` → `ColumnDialog` → `build_profile` → a `Profile`.
2. `MainWindow._conditions` reads the Conditions step and, if fO2 was given as
   a buffer offset, calls `thermo.buffers.log_fo2_from_delta` → a `Conditions`.
3. `MainWindow._model` assembles the `Conditions`, the chosen
   `DiffusionCoefficient` from `coefficients.registry`, an `InitialCondition`,
   a `Geometry`, two `BoundaryCondition`s, an optional `ThermalHistory` and the
   beam sigma → a `fitting.model.DiffusionModel`.
4. `gui/workers.FitWorker` runs `fitting.fit.fit_time` on a `QThread`.
5. `fit_time` scans chi² over a logarithmic grid of times, then refines with
   `scipy.optimize.least_squares`. Each trial time calls
   `DiffusionModel.profile`, which asks `can_use_analytical()` and then either
   `solvers.analytical.step_infinite` or `solvers.numerical.solve_1d`, and
   convolves with `solvers.convolution.gaussian_convolve` if a beam sigma is set.
6. `MainWindow._fit_done` draws the result with `gui/plot_widget.ProfilePlot`
   and rebuilds the summary sidebar.
7. **Uncertainty** (optional) repeats steps 3–5 once per draw inside
   `fitting.montecarlo.run`, called from `gui/workers.MonteCarloWorker`.
8. **Export** calls `dataio.export.save_results`, which calls
   `methods_paragraph` and `collect_citations`, which walk back through the
   model and pull citation keys out of `references.REFERENCES`.

**A profile from an image.** File > Extract profile from image (or "From an
image..." on the Data step) opens `gui/image_extractor.ImageExtractorDialog`:
`dataio.images.load_image` → `value_map` (grey, channel or inverted colour
legend) → the user draws a guideline → `dataio.image_profiles.extract_profiles`
→ `write_workbook`. Loading that workbook, or pressing "Use in Diffusor",
goes through `MainWindow.load_image_table` →
`gui/image_calibration.ImageCalibrationDialog` (pixel size, value-to-composition
map from `dataio.greyscale.calibrate`) → `image_profiles.composition_table` →
the same `build_profile` as step 1 above.

**Internal units.** The library works in **micrometres and seconds** inside
the model layer: `fitting.model.M2_PER_S_TO_UM2_PER_S = 1e12` converts D from
m²/s once, so the grid stays O(1–1000) instead of pairing a 1e-5 m grid with
D ≈ 1e-20 m²/s. Coefficients always return **m²/s**. Distances in a loaded
file are converted to micrometres by `dataio.profiles.build_profile`.
Temperature is K inside, °C in the interface. Pressure is Pa inside, MPa in
the interface. fO2 is log10 bar inside, and each coefficient declares the unit
its own law expects (`fo2_unit`: `Pa`, `bar` or `atm`).

---

## 2. Scripts

Each entry gives the role, when it runs, the variables that matter and the
source of every equation or number in it. "Called by" means the direct caller
in a normal run.

### 2.1 Entry points

#### `diffusor/__init__.py`
Package docstring and `__version__ = "0.1.0"`. The version string is written
into every exported JSON and methods block by `dataio/export.py`, shown in
Help > About, and compared with the newest GitHub release by `updates.py`.
It is the only place the version is written: `pyproject.toml` declares it
`dynamic` and reads it with `[tool.setuptools.dynamic] version = { attr =
"diffusor.__version__" }`. To release, change it here, commit, tag `vX.Y.Z` and
publish a GitHub release with that tag.

#### `diffusor/__main__.py`
`from .gui.app import main`, then `raise SystemExit(main())` under an
`if __name__ == "__main__"` guard. Makes `python -m diffusor` start the window.
The guard is required by the parallel Monte Carlo: on Windows and macOS each
worker process re-imports the main module under another name, and without the
guard every worker would open its own window.

#### `Start Diffusor.bat`
The double-click start for someone who downloaded the ZIP. On the first run
it finds Python 3.10 or newer (the `py -3` launcher, then `python`), makes a
virtual environment in `.venv` beside it, runs `pip install -e .`, checks that
`diffusor.gui.main_window` imports, and writes `.venv/diffusor-ready.txt`.
Later runs skip straight to the start: `pythonw -m diffusor` from the venv,
so no console stays open. Each run also remakes a `Diffusor.lnk` shortcut
with `icons/diffusor.ico` beside it (ignored by git), so the shortcut follows
the folder when it is moved. `.gitattributes` keeps the file's Windows line
endings, which cmd needs, in every checkout and GitHub ZIP.

#### `Start Diffusor.command`
The same for macOS (double-click in Finder, which runs it in Terminal) and
Linux (run from a terminal). It looks for Python 3.10 or newer as
`python3.14` down to `python3.10`, then `python3` and `python`, because the
`python3` that ships with macOS is often 3.9; makes `.venv`, installs with
`pip install -e .`, checks the import and writes the same ready marker. It
then starts `python -m diffusor` with `nohup` in the background, so closing
the Terminal window leaves Diffusor running. It also accepts a Windows-layout
`.venv` (`Scripts/python.exe`), which lets it be tested in Git Bash. Git keeps
it executable (mode 100755) and `.gitattributes` forces Unix line endings; the
GitHub ZIP preserves both. It is not signed, so macOS Gatekeeper asks the user
to allow it once (the README gives the steps).

#### `diffusor/gui/app.py`
Creates the `QApplication`, forces the matplotlib `QtAgg` backend before the
window is imported, sets the application icon (`icons/diffusor.png`, which
every window inherits) and shows a `MainWindow`. On Windows it first sets an
explicit AppUserModelID, so the taskbar shows Diffusor's icon instead of
grouping the window under `pythonw.exe`. Also the `diffusor` launcher, declared
under `[project.gui-scripts]` in `pyproject.toml` so it opens no console window;
a desktop shortcut to it can use `icons/diffusor.ico`. 1.5 s after the window
is shown it calls `MainWindow.startup_update_check`, so a slow network never
delays the start. `MainWindow` itself never checks on construction, which keeps
the tests off the network.

#### `diffusor/updates.py`
Finds out whether a newer Diffusor has been released. Standard library only,
no Qt. `fetch_latest` GETs `https://api.github.com/repos/ArminTimar/Diffusor/releases/latest`
(the newest published release; drafts and pre-releases are not returned) with a
6 s timeout, a `User-Agent` (GitHub refuses requests without one) and a 1 MB
cap on the reply, and returns a `ReleaseInfo(version, tag, url, notes)`. Every
failure (no release yet, rate limit, offline, timeout, unreadable or oversized
reply) raises `UpdateError` with a sentence fit to show the user.
`parse_version` reads `v0.2.0` as `(0, 2, 0)` and rejects anything else,
including `0.2.0rc1`; `is_newer` compares as numbers, padded with zeros, so
`0.10.0 > 0.9.0` and `1.0 == 1.0.0`; a tag it cannot read gives no notice.
The `html_url` from the reply is opened in the user's browser, so it is used
only if it starts with `https://github.com/ArminTimar/Diffusor/releases/`;
otherwise the page is rebuilt from the tag. `due` is true when no check has
succeeded in the last 20 hours (`CHECK_EVERY`). The module downloads and
changes nothing: the user reads the release page and downloads the ZIP as at
the first install.

---

### 2.2 Foundations

#### `diffusor/constants.py`
Physical constants, atomic weights, oxide molar masses and unit factors.
Nothing computes here except the oxide masses.

* **Called by** everything. `thermo/units.py`, `coefficients/base.py`
  (`R_GAS`), every coefficient module, `fitting/fit.py` (`SEC_PER_YEAR`).
* **Variables.** `R_GAS = 8.314462618` J/mol/K, `K_BOLTZMANN`, `N_AVOGADRO`,
  `EV_TO_J`; `T_KELVIN_OFFSET = 273.15`; `BAR_TO_PA`, `GPA_TO_PA`,
  `KBAR_TO_PA`; `SEC_PER_YEAR`; `ATOMIC_MASS` (16 elements);
  `OXIDE_MOLAR_MASS` (15 oxides, computed from `ATOMIC_MASS` by `_oxide()`);
  `OXIDE_CATIONS` (cations per formula unit, for wt% to cation moles).
* **Sources.**
  * Fundamental constants: CODATA 2018 (`codata2018`). R, k_B, N_A and the
    eV conversion are exact under the 2019 SI redefinition.
  * `T_KELVIN_OFFSET`: ITS-90 definition of the Celsius scale.
  * `SEC_PER_YEAR = 365.25 x 86400`, the **Julian year**, IAU convention
    (`iau_julian_year`). Note that `dataio/export.py` writes the Monte Carlo
    CSV with a hard-coded `3.15576e7` s/yr, which is the same number.
  * Atomic weights: IUPAC 2021 standard atomic weights, conventional values
    for elements quoted as an interval (`iupac2021`).

#### `diffusor/references.py`
The single citation registry, and the reason the rest of the traceability
machinery can work.

* **Called by** `coefficients/base.py` (`cite`, `get`), `dataio/export.py`
  (`format_reference`), `gui/richtext.py`, `scripts/build_references.py`, and
  the test suite.
* **Contents.** `Reference` is a frozen dataclass (`key, authors, year, title,
  journal, volume, pages, doi, note, kind`) with `short()` giving
  "Author (year)" and `full()` a bibliographic line. `REFERENCES` holds 80
  entries. `cite(key, where)` builds "Crank (1975), eq. 2.14". `to_bibtex`
  renders a BibTeX record, switching field names on `kind`
  (`article` / `book` / `incollection` / `misc`).
* **Rule.** A key that does not resolve raises `KeyError` with the instruction
  to add it here, and the test suite fails. Nothing else in Diffusor is
  allowed to name a source in a machine-readable field without a key.

---

### 2.3 `diffusor/thermo/` — buffers and conversions

#### `thermo/units.py`
Temperature, pressure, time and length conversions, plus the composition
variables the diffusion laws are written against.

* **Called by** `dataio/profiles.build_profile` (`length_to_m`,
  `composition_variable`), `fitting/fit.py` and the GUI (`human_time`),
  `scripts/make_examples.py`.
* **Functions.** `c_to_k` / `k_to_c`; `pressure_to_pa` (Pa, bar, kbar, MPa,
  GPa); `time_to_seconds` / `seconds_to` / `human_time` (s, min, h, d, yr,
  kyr, Myr); `length_to_m` (m, mm, um, nm, cm);
  `oxide_wt_to_cation_moles`; `ratio_A_over_AplusB`; `x_fe_from_oxides`;
  `x_an_from_oxides`; `composition_variable(a, b, mode, oxide_a, oxide_b)`
  with modes `A`, `A/(A+B)`, `B/(A+B)`, `A-B`.
* **Sources and conventions.**
  * `X_Fe = Fe/(Fe+Mg)` molar, all Fe as Fe2+: the convention of
    Dohmen & Chakraborty (2007), Dohmen et al. (2016) and Mueller et al.
    (2013), that is, of the laws that consume it.
  * `X_An = Ca/(Ca+Na+K)` molar: Deer, Howie & Zussman (1992) (`deer1992`).
  * `Fo = 100 Mg/(Mg+Fe)` and `Mg# = 100 Mg/(Mg+Fe)`: same source.
  * `X_Usp` from Ti apfu on a 3-cation basis, `x_Ti = X_Usp/3` for
    (Ti_x Fe_(1-x))3 O4: the convention of Van Orman & Crispin (2010).
    **A Stormer-style recalculation is deliberately not applied**; the user
    supplies Ti apfu or X_Usp.
  * Molar masses from `constants.py`, that is IUPAC 2021.

#### `thermo/buffers.py`
Solid oxygen buffers as `log10 fO2` in bar, for T in K and P in Pa.

* **Called by** `gui/main_window._conditions` and `_update_fo2_label` (every
  time the temperature, pressure, buffer or offset spinner changes),
  `fitting/montecarlo._draw_conditions` (**once per draw**, at the sampled
  temperature), and the two coefficient laws that quote fO2 relative to IW
  (`opx_FeMg_ganguly_tazzoli1994`, `opx_Lu_dias2025`).
* **Variables.** `_FROST` holds the (A, B, C) triples for IW, WM, FMQ, NNO and
  HM. `_ALIASES` maps QFM to FMQ and IQF to IW. `BUFFER_CITATIONS` is what
  `scripts/build_references.py` renders into `REFERENCES.md`.
* **Equations and sources.**
  * `frost1991` parameterisation (the default):
    `log fO2 = A/T + B + C (P_bar - 1)/T`, with the coefficients of
    Frost (1991) Reviews in Mineralogy 25, Table 1 (`frost1991`):

    | buffer | A | B | C |
    | --- | --- | --- | --- |
    | IW | -27489 | 6.702 | 0.055 |
    | WM | -32807 | 13.012 | 0.083 |
    | FMQ | -25096.3 | 8.735 | 0.110 |
    | NNO | -24930 | 9.36 | 0.046 |
    | HM | -25700.6 | 14.558 | 0.019 |

  * `oneill` parameterisation, 1 bar free energies converted by
    `log fO2 = dG/(R T ln10)`:
    * NNO, 2 Ni + O2 = 2 NiO:
      `dG = -478967 + 248.514 T - 9.7961 T lnT` J/mol, from
      O'Neill & Pownceby (1993) (`oneill_pownceby1993`).
    * FMQ, 3 Fe2SiO4 + O2 = 2 Fe3O4 + 3 SiO2:
      `dG = -587474 + 1584.427 T - 203.3164 T lnT + 0.092710 T^2` J/mol, from
      O'Neill (1987) (`oneill1987`).
    * The Frost (1991) pressure term is added to both so the two
      parameterisations share a pressure dependence.
* **Why it matters for uncertainty.** The buffer is a function of T. The Monte
  Carlo re-evaluates it at every sampled temperature, so fO2 and T stay
  correlated by construction instead of being drawn independently.

---

### 2.4 `diffusor/minerals/` — phases, axes, composition variables

#### `minerals/base.py`
The `CompositionVariable`, `Species` and `Mineral` dataclasses, and the
anisotropy relation.

* **Called by** `gui/main_window` (mineral and species drop-downs, and
  `mineral.composition_variable.key` when building `Conditions.X`), and by
  `coefficients/base.DiffusionCoefficient.D` when a traverse is given as
  angles.
* **Equation.** `direction_factor`:
  `D_V = D_a cos^2(alpha) + D_b cos^2(beta) + D_c cos^2(gamma)`,
  Costa & Chakraborty (2004) (`costa_chakraborty2004`); also DMG Short Course
  2025 Lecture 6 and the course script `MCdiff_OlFo_MO.m` (`dmg2025`). The
  function refuses angles whose direction cosines do not close to 1 within
  0.05, which catches the common case of three angles measured independently.

#### `minerals/definitions.py`
The six supported minerals and their conventions.

* **Composition variables.** `X_FE` (Fe/(Fe+Mg)), `X_FO` (Mg/(Mg+Fe)),
  `X_AN` (Ca/(Ca+Na+K)), `X_OR` (K/(K+Na+Ca)), `X_TI` (Ti per cation site).
* **Minerals.** olivine, opx, cpx, plagioclase, kfeldspar, magnetite, each
  with formula, crystal system, axes, a species dict and an `isotropic` flag
  (plagioclase, K-feldspar and magnetite are flagged isotropic).
* **Sources behind the notes.**
  * Olivine, about 6x anisotropy of D//[001]: Dohmen & Chakraborty (2007),
    confirmed in natural crystals by Hartley et al. (2016).
  * Opx, `D_a = D_c/3.5`, space group Pbca: Dohmen et al. (2016).
  * Cpx anisotropy of Fe-Mg is unconstrained; Mueller et al. (2013) measured
    [001] only: Cherniak & Dimanov (2010).
  * Plagioclase X_An treated as frozen because NaSi-CaAl interdiffusion is
    orders of magnitude slower: Grove et al. (1984). Trace elements coupled to
    the An gradient: Costa et al. (2003). Isotropy for Mg:
    Van Orman et al. (2014).
  * K-feldspar: no resolvable orientation dependence for Sr and Ba
    (Cherniak 1996, 2002); Ba about 50x slower than Sr, so paired profiles
    test whether a boundary is diffusive at all (Chamberlain et al. 2014).
  * Magnetite is cubic and so isotropic; competing vacancy (fO2^(2/3)) and
    interstitial (fO2^(-2/3)) mechanisms give D a minimum in both T and fO2
    (Van Orman & Crispin 2010, Table 12). This is the formulation
    Tomiya et al. (2013) used for Shinmoedake.

---

### 2.5 `diffusor/coefficients/` — the registry of published laws

This is where most of the citation weight sits. Every entry is one published
diffusion law, wrapped so that the equation, its number, the calibration
ranges, the published uncertainties and the provenance of the numbers travel
with the function that evaluates it.

#### `coefficients/base.py`
The framework. No diffusion law lives here.

* **Called by** every coefficient module, `fitting/model.py`
  (`Conditions`, `D_sampled`), `fitting/montecarlo.py` (`sample`),
  `gui/main_window` (`Conditions`), `gui/richtext.py` (`describe`).
* **Classes.**
  * `Conditions(T_K, P_Pa, log_fo2_bar, X, axis, angles_deg)`. `X` is a dict
    of composition variables and each value may be a scalar or an array, so
    the numerical solver can evaluate D node by node. `log_fo2_Pa` adds 5 to
    the bar value (1 bar = 1e5 Pa).
  * `Parameter(name, value, sigma, unit, sigma_level, description)` with
    `sigma_1s` halving a published 2-sigma value.
  * `Range(lo, hi, unit)` with `contains`, used for the extrapolation
    warnings.
  * `DiffusionCoefficient`: the law itself. Key fields: `func` (the evaluator),
    `params`, `covariance` + `cov_order`, `sigma_logD`, `axis_factors` +
    `reference_axis`, `requires`, `needs_fo2` + `fo2_unit`, the four `Range`
    fields, `verified` + `verified_from`, `superseded_by` +
    `superseded_note`, `secondary_citations`, `notes`, `recommended`.
* **Methods that matter.**
  * `D(cond)` evaluates `func` along the reference axis, then applies either
    the direction-cosine relation (if `angles_deg` is set) or a single axis
    factor (if `axis` is set).
  * `sample(rng, mode)` draws one Monte Carlo perturbation, in one of four
    modes described below.
  * `D_sampled(cond, overrides)` applies those overrides, including the
    special `__log10_D_offset__` key.
  * `check_conditions(cond)` returns the warning strings shown on every fit:
    out-of-range T, fO2, P or composition, the unverified flag, and the
    superseded flag.
  * `describe()` and `methods_sentence()` feed the reading panes and the
    exported methods block.
* **Sampling modes, and the argument behind them.** Arrhenius fits give
  typically positively correlated `ln D0` and `Q`; drawing them independently
  inflates the spread of D at the temperature of interest by a large factor.
  * `covariance`: multivariate normal over `cov_order`. The approach of
    Mutch et al. (2021), DFENS. Only the two Grocolas et al. (2025) laws
    carry one, and theirs is constructed (section 4.3).
  * `logD_at_T`: draw `log10 D` directly at the working temperature from the
    scatter the paper reports (e.g. "reproduces the data within 1 log unit",
    Mueller et al. 2013). The honest fallback, and the default whenever
    `sigma_logD` is set.
  * `independent`: each parameter separately. Offered only to reproduce the
    older assumption, for instance NIDIS (Petrone et al. 2016).
  * `none`.
  `default_sampling_mode()` picks covariance, then logD_at_T, then
  independent.
* **Helper.** `arrhenius(D0, Q_J, T_K, dV, P_Pa, P0_Pa)`:
  `D = D0 exp(-(Q + (P - P0) dV)/(R T))`. The standard form with an activation
  volume, Crank (1975) section 11 and Costa et al. (2008) eq. 20.

#### `coefficients/registry.py`
Builds `REGISTRY` by importing the six mineral modules in a fixed order
(olivine, opx, cpx, plagioclase, kfeldspar, magnetite) and refusing duplicate
keys at import time.

* **Called by** `gui/main_window._refresh_coefficients` (step 5 list, filtered
  by mineral and species), `scripts/build_references.py`, the tests.
* `list_coefficients()` sorts recommended first, then verified, then
  alphabetically. That sort is what the chooser shows.

#### The six law modules

Totals: **54 entries, 47 verified, 13 recommended.**

| module | entries | species covered |
| --- | --- | --- |
| `olivine.py` | 3 | Fe-Mg (TaMED, PED, Chakraborty 1997) |
| `opx.py` | 12 | Fe-Mg (5), Mg (3 axes), Lu, Ce, Eu |
| `cpx.py` | 3 | Fe-Mg (2), Ca-Mg |
| `plagioclase.py` | 11 | Mg (3), Sr (3), Ba (2), Li (2), NaSi-CaAl |
| `kfeldspar.py` | 3 | Sr, Ba, Ti |
| `magnetite.py` | 22 | Ti, Fe, Mn, Co, Cr, Al, Mg, Fe-Ti |

Each law's equation and its source are listed in section 4.1 and discussed
paper by paper in section 3.5. What follows is what is specific to each
module as code.

##### `coefficients/olivine.py`
Two regimes of Dohmen & Chakraborty (2007) plus one older comparison law.
`_dohmen_chakraborty_tamed` and `_dohmen_chakraborty_ped` differ only in the
intercept, the activation energy and whether an fO2 term is present. The
anisotropy dict `_OLIVINE_ANISO = {"c": 1, "b": 1/6, "a": 1/6}` is shared by
all three. The two Dohmen & Chakraborty entries are **verified** against the
printed equations 27 and 28 (pp. 424-425) and the erratum (PCM 34:597-598),
which corrects their composition term from `3 X_Fe` to `3 (X_Fe - 0.1)`; that
reference is the `XFe_ref = 0.1` parameter. TaMED's `fo2_range` starts at the
paper's 1e-10 Pa regime boundary and has no upper limit, as the paper gives
none. Only `ol_FeMg_chakraborty1997` remains unverified.
Diffusor does not switch between TaMED and PED automatically; the user picks.

##### `coefficients/opx.py`
The busiest module.
* `LOG_FO2_PA_SWITCH = -10.0` and `dias2025_regime()` implement the two-regime
  structure of Dias et al. (2025): eq. 22 above the switch, eq. 23 at or below.
  The two do not join smoothly there, by design, because the authors infer a
  change of mechanism.
* `_dias2025_m` is the temperature-dependent composition exponent, eqs 24 and
  25, whose slopes have opposite signs in the two regimes.
* `_ganguly_tazzoli` carries a `use_fo2` pseudo-parameter (1 or 0) so the same
  law can be run with or without the hypothesised fO2^(1/6) term. The
  no-fO2 form is also registered separately as
  `opx_FeMg_ganguly_tazzoli1994_nofo2`, because that is the entry that
  reproduces Ostorero et al. (2022).
* The three REE entries and the three Schwandt axes are generated in loops.
  The Lu entry alone takes `verified=False`.

##### `coefficients/cpx.py`
Three plain Arrhenius laws through the shared `_plain_arrhenius`. The
Mueller et al. (2013) entry documents an internal inconsistency in its own
source (section 4.3). The Dimanov & Sautter (2000) entry exists to reproduce
published NIDIS timescales and is marked unverified because the numbers come
from a table footnote in Petrone et al. (2016).

##### `coefficients/plagioclase.py`
Also holds the plagioclase-specific physics that the solver needs.
* `ACTIVITY_A` (kJ/mol): Mg +15.8, Sr -17.4, Ba -35.1, Li -1.7, K -8.0,
  Rb -15.9. From Dohmen, Faak & Blundy (2017) RiMG 83, Appendix Fig. A1
  caption and Table 1, computed at 1200 C from the lattice-strain model of
  Dohmen & Blundy (2014), through `-R T ln(gamma_i) = A_i X_An + B_i`
  (Appendix eq. A6). `ACTIVITY_A_NOTE` records that the same appendix also
  tabulates 900 C values and that the DMG 2025 course script uses
  `A_Sr = -15.1` kJ/mol instead.
* `activity_theta(species, T_K) = A_i x 1000 / (R T)` is the `theta` the
  numerical solver multiplies the `dX_An/dx` flux term by.
* `equilibrium_profile(X_An, T_K, species, C_ref, X_An_ref)`:
  `C_eq(x) = C_ref exp(theta (X_An(x) - X_An_ref))`, Dohmen et al. (2017)
  Appendix eq. A13 with the rim condition of eq. A14. This is the
  `equilibrium_plag` initial condition offered for plagioclase.
* `_compensated_covariance(s_a, s_b, s_Q)` builds the 3x3 matrix over
  (a, b, Q) in which `b` and `Q` are **perfectly** correlated and `a` is
  independent. This encodes what Grocolas et al. (2025) section 4.3 state they
  did in their own Monte Carlo. It is Diffusor's construction of their stated
  assumption, not a matrix printed in the paper.

##### `coefficients/kfeldspar.py`
Three abstract-level Arrhenius laws (Cherniak 1996 Sr, Cherniak 2002 Ba,
Cherniak & Watson 2020 Ti) on `_arrhenius_law`. The Sr and Ba `sigma_logD`
values (0.03 and 0.12) are **not** from the diffusion papers: they are the
uncertainties Chamberlain et al. (2014) derived from those D0 and Q. The Ti
entry has `sigma_logD=None` and a zero sigma on D0, so Monte Carlo coefficient
sampling there can only vary Q, which the notes warn overstates the spread.

##### `coefficients/magnetite.py`
The most structured module.
* `TABLE12_PURE` and `TABLE12_XTI02`: the (D_V0, Q_V, D_I0, Q_I) quadruples of
  Van Orman & Crispin (2010) Table 12, for x_Ti = 0 (Cr, Al, Fe, Co, Mn, Ti)
  and x_Ti = 0.2 (Fe, Co, Mn, Ti). The vacancy activation energies are
  **negative**, which is what produces the diffusion minimum.
* `_table12_D` evaluates
  `D* = D_V0 exp(-Q_V/RT) a_O2^(2/3) + D_I0 exp(-Q_I/RT) a_O2^(-2/3)` with
  `a_O2 = fO2/1 atm`, and `ATM_IN_BAR = 1.01325` does that conversion.
* Intermediate x_Ti is obtained by **log-linear interpolation** between the
  two tabulated compositions, clipped to `x_Ti/0.2 <= 1.5`. That interpolation
  is Diffusor's reconstruction of how Tomiya et al. (2013) reached their
  published Shinmoedake numbers, and it reproduces their Ti values to about
  1 per cent.
* `SIEVWRIGHT_TABLE5` holds 21 elements as
  (log D_V1, log D_I1, log fO2 at the minimum, log D at the minimum) at
  1150 C, with the last two kept only so the tests can verify the
  transcription. `sievwright_D_1150` evaluates
  `D = D_V1 fO2^(2/3) + D_I1 fO2^(-2/3)` with fO2 in bar.
  **Only 6 of the 21 elements are registered** as coefficients (Ti, Mn, Co,
  Cr, Al, Mg); the rest are transcribed but unexposed.
* `_make_sievwright_func(species, scale_with_table12)`: away from 1150 C each
  branch is scaled by `exp(-Q/R (1/T - 1/1423.15))` using the Table 12
  energies of the same element. **This scaling is Diffusor's construction and
  is absent from Sievwright et al. (2020)**, which the entry says plainly. Mg
  has no Table 12 energies, so that entry is valid at 1150 C only and warns
  at any other temperature.
* The eight `mt_*_aggarwal2002_*` entries are plain Arrhenius fits computed by
  Van Orman & Crispin (2010) **along** the WM or MH buffer (their Tables 10
  and 11, using the buffer equations of Huebner 1971). They are valid on that
  buffer only.

---

### 2.6 `diffusor/solvers/` — the diffusion equation

#### `solvers/geometry.py`
`GEOMETRIES = {"plane": 0, "cylinder": 1, "sphere": 2}` is the index `m` in
`dC/dt = (1/x^m) d/dx [x^m D dC/dx]`. Crank (1975) sections 2 and 4 (plane),
5 (cylinder), 6 (sphere). `make_grid` and `suggest_grid` build the uniform
grid the numerical solver requires. `stability_dt(dx, D_max, courant)` returns
`courant dx^2 / D`, the explicit stability limit of Crank (1975) eq. 8.33; the
docstring records that the DMG 2025 course scripts use a Courant number of
0.2 to 0.4. `GEOMETRY_NOTES` is the text shown on the Model step, and it says
that choosing the geometry is a petrological decision.

#### `solvers/boundary.py`
`BoundaryCondition(kind, value)` with kinds `dirichlet`, `neumann` and
`symmetry`. The Dirichlet value may be a constant or a callable `f(t)`.
Crank (1975) section 1.3 and section 8.4; Costa et al. (2008), "Boundary
conditions". The interface never shows these words: `gui/main_window` maps the
four plain-language choices (`far`, `rim_melt`, `rim_closed`, `centre`) onto
them in `_boundary()`, and `far`/`rim_melt` both become Dirichlet at the
plateau value while `rim_closed`/`centre` become zero flux.

#### `solvers/initial.py`
The initial condition, kept explicit because it is the largest systematic
error in the method (Costa et al. 2008, "Initial conditions";
Shea et al. 2015).

* **Forms.** `step`, `multi_step`, `plateau_rim`, `table`, and
  `equilibrium_plag` (which calls back into
  `coefficients.plagioclase.equilibrium_profile`, Dohmen et al. 2017 App.
  eq. A13).
* `step` optionally smooths with an error function. A node sitting exactly on
  the interface is given the cell average of the two plateaus; without that
  the discretised step sits half a cell off and the numerical solution is
  shifted by dx/2 relative to the Crank solution. The docstring notes that the
  DMG 2025 script `Diff_Model_Sr_in_Plag_implicit.m` applies 10 explicit
  smoothing sweeps for the same reason, and warns that smoothing must stay
  small compared with the diffusion length.
* `guess_step_from_data(x, C, plateau_fraction=0.15)` is what the **Guess**
  button on the Model step calls. Plateaus are the medians of the outer 15 per
  cent of points at each end; the interface is placed at the **midpoint
  crossing** inside the interior window, not at the steepest gradient, because
  on a broad coarsely sampled profile the noisiest single interval can
  out-gradient the real step. It falls back to the steepest gradient of a
  lightly smoothed profile if there is no crossing.

#### `solvers/analytical.py`
Closed forms for constant D. Every function takes the product `Dt`, so a
non-isothermal problem can be handled by substituting the effective integral
from `history.py`.

| function | problem | source |
| --- | --- | --- |
| `step_infinite` | infinite medium, step at x0 | Crank (1975) eq. 2.14; Costa et al. (2008) eq. 9 |
| `semi_infinite_fixed_surface` | x >= 0, surface held at C_s | Crank eq. 2.45 |
| `band_infinite` | band of half-width h in an infinite medium | Crank eq. 2.15 |
| `plane_sheet` | -l < x < l, surfaces at C_1 | Crank eq. 4.17 |
| `cylinder` | 0 < r < a, surface at C_1 | Crank eq. 5.22 |
| `sphere` | 0 < r < a, surface at C_1 | Crank eq. 6.18 |
| `fraction_plane_sheet` | M_t/M_inf | Crank eq. 4.18 |
| `fraction_cylinder` | M_t/M_inf | Crank eq. 5.23 |
| `fraction_sphere` | M_t/M_inf | Crank eq. 6.20 |

`step_infinite` is the only one the fitting path uses today; it is the form
NIDIS (Petrone et al. 2016, Methods) fits with free x0 and 2 sqrt(Dt). The
series solutions use 100 to 200 terms; the cylinder uses `scipy.special.jn_zeros`
for the roots of J0, the sphere takes the L'Hopital limit at r = 0.

#### `solvers/convolution.py`
`gaussian_convolve(x, C, sigma)` convolves the **model** with a Gaussian
before it is compared with the data, so the fitted time is corrected for beam
smearing rather than inflated by it. Edges are padded with the end values.
Source: Ganguly, Bhattacharya & Chakraborty (1988) (`ganguly1988`);
Bradshaw & Kent (2017) (`bradshaw_kent2017`); DMG 2025 Practical 5.
`resolution_warning(Dt, sigma, factor=3)` fires when `2 sqrt(Dt) < 3 sigma`,
the Bradshaw & Kent rule of thumb for "report this as an upper bound".

#### `solvers/history.py`
`ThermalHistory` is a piecewise-linear T(t) path with `isothermal`, `linear`
and `piecewise` constructors and `shifted_to_end(t_total)`, which rescales the
time axis so a cooling path of fixed shape can have its duration fitted.
`effective_Dt(D_of_T, history, t_total)` evaluates
`integral_0^t D(T(t')) dt'` by the trapezoidal rule on 2001 points.
Source: Crank (1975) section 7.2, eq. 7.7; used in geospeedometry by
Lasaga (1983); reviewed by Costa et al. (2008). The numerical solver instead
re-evaluates D at every step, so the two routes agree.

#### `solvers/numerical.py`
The Crank-Nicolson solver. The single most load-bearing module in the package.

* **Called by** `fitting.model.DiffusionModel._profile_numerical`, once per
  trial time, so hundreds of times per fit and hundreds of thousands of times
  in a Monte Carlo.
* **Equation solved.**
  `dC/dt = (1/x^m) d/dx [ x^m ( D(C,x,T) dC/dx - theta D C dX_An/dx ) ]`.
  The first part is Crank (1975) eq. 1.7 (plane), eq. 5.4 (cylinder), eq. 6.3
  (sphere). The second is the activity term of Costa et al. (2003) eq. 7 with
  `theta = A_i/(R T)`, in the form of Dohmen, Faak & Blundy (2017) Appendix
  eqs A7-A8; it is active only when an `an_profile` is supplied.
* **Discretisation.** Conservative finite volumes on a uniform grid with D at
  half nodes, `D_{i+1/2} = (D_i + D_{i+1})/2` (Dohmen et al. 2017 App. eqs
  A17-A19). Time integration by the theta-scheme
  `(I - theta dt L) C^{j+1} = (I + (1-theta) dt L) C^j`, with theta = 1/2
  giving Crank-Nicolson (Crank section 8.4 eq. 8.35; Dohmen et al. 2017 App.
  eq. A21), 0 explicit (Crank eq. 8.31) and 1 fully implicit. The operator is
  lagged at C^j. The activity term follows Dohmen et al. (2017) App. eq. A20.
  The banded system is solved by `scipy.linalg.solve_banded`.
* **Boundary rows.** Dirichlet becomes an identity row carrying the value at
  t^{j+1}. Zero flux at x > 0 uses a mirror ghost node. Symmetry at x = 0 for
  m > 0 uses the L'Hopital limit `dC/dt = (m+1) D d2C/dx2`, Crank section 8.5
  eq. 8.45.
* **Time-step plan.** dt starts at the explicit stability limit
  `courant dx^2 / D_max` (Crank eq. 8.33) and grows geometrically by
  `dt_growth = 1.07` up to a cap of `t_total / min_steps` (default 400
  steps). Crank-Nicolson is unconditionally stable, so dt is limited by
  accuracy, not stability; the small early steps resolve the sharp initial
  transient and the large later ones keep a 100 kyr run to a few hundred steps
  instead of millions. If even the cap would need more than `max_steps`
  (200000) the step is enlarged and a warning is attached to the result.
* **Diagnostics.** `NumericalResult` carries `mass_initial` and `mass_final`
  (computed with the `x^m` weight), the step count, the nominal dt and any
  warnings. A `progress` callback can abort a run.
* **Accuracy, from the test suite.** Second-order convergence against the
  Crank solutions, agreement to about 3 parts in 1e6, and mass conservation to
  1 part in 1e8 in a closed system.

---

### 2.7 `diffusor/fitting/` — forward model, fit, Monte Carlo

#### `fitting/model.py`
`DiffusionModel` bundles everything needed to turn a time into a predicted
profile, and is used identically by the least-squares fit and by every Monte
Carlo draw.

* **Called by** `gui/main_window._model` (builds it), `fitting/fit.py`
  (`profile`), `fitting/montecarlo.py` (shallow-copied per draw),
  `dataio/export.py` (reads its fields for the methods block),
  `scripts/make_examples.py`.
* **Fields.** `coefficient`, `conditions`, `initial`, `geometry`, `bc_left`,
  `bc_right`, `history`, `beam_sigma_um`, `n_nodes` (default 401), `x_grid`,
  `composition_dependent`, `comp_key`, `an_profile`, `activity_theta`,
  `force_numerical`, `boundaries_far`, `fo2_buffer`.
* **`fo2_buffer`** is `(buffer name, offset)` when fO2 was given relative to a
  buffer (the interface sets it). `conditions_at(T)` then recomputes fO2 from
  the buffer at every temperature, so along a cooling path fO2 falls with the
  buffer instead of staying at its starting value. The numerical solver
  (`_D_um2s`) and the effective-Dt integral both evaluate D through it.
* **`reference_C()`** is the midpoint of the two plateaus. Wherever one D
  stands for the whole profile (the closed form, the diffusion lengths in the
  warnings) a composition-dependent law is evaluated there, never at the
  single placeholder composition in `conditions.X`, which the numerical solver
  replaces node by node.
* **`M2_PER_S_TO_UM2_PER_S = 1e12`.** D is converted from m^2/s to um^2/s once
  in `_D_um2s`, which is what keeps the linear algebra well conditioned.
* **`can_use_analytical()`** returns a decision **and a reason**, and the
  reason is shown on the Model step and written into the methods block. The
  closed form is refused when: the user forced numerical; an end of the
  profile is a real rim or centre; D depends on composition (Crank 1975
  section 7.2); the anorthite activity term is on (no closed form,
  Costa et al. 2003); the initial profile is not a sharp step; the geometry is
  not plane; or the step is smoothed. Otherwise the route is
  "sharp step, constant D, plane geometry, plateaus continue"
  (Crank 1975 eq. 2.14).
* **`profile(t, x_out, overrides)`** is the single entry point every fit goes
  through. On the analytical route it evaluates `step_infinite` directly, and
  if a beam sigma is set it re-evaluates on the model grid, convolves, and
  interpolates back. On the numerical route it evaluates the initial condition
  on the grid, runs `solve_1d`, convolves, and interpolates.
* **`warnings(t)`** collects the coefficient's range and provenance warnings
  (ranges tested at the profile's own plateau compositions through
  `_range_conditions`, and at every temperature of a cooling path with its
  own fO2),
  the beam-resolution warning, a **far-field warning** when
  `2 sqrt(Dt) > 0.6 x` (Dt at `reference_C()`) the distance from the interface
  to the end of the
  traverse (the semi-infinite assumption is breaking down), and, for plane
  geometry, the standing caveat that 1-D modelling of a 3-D crystal gives a
  maximum estimate and that sectioning biases it further (Shea et al. 2015;
  Krimer & Costa 2017).

#### `fitting/objective.py`
Weighted least squares bookkeeping. `weights_from_sigma` returns `1/sigma`, or
uniform weights with no uncertainty column (in which case the reported chi2 is
only a relative measure). `residuals` returns `(obs - model)/sigma`.
`statistics` returns `FitStatistics(chi2, reduced_chi2, dof, rmse, r_squared,
n_points, n_params)`. Conventional weighted least squares; no external source.

#### `fitting/fit.py`
Finds the time.

* **Called by** `gui/workers.FitWorker`, `gui/workers.CompareWorker`,
  `fitting/montecarlo.run` (once for the base fit and once per draw),
  the round-trip tests.
* **Strategy.** A coarse logarithmic scan of chi2 over time locates the global
  minimum, then `scipy.optimize.least_squares` refines in log10 t with bounds.
  The scan matters because chi2(t) is very flat at long times, so a local
  optimiser started in the wrong decade converges to the edge.
* **Search window.** `T_MIN_DEFAULT = 1e2` s, `T_MAX_DEFAULT = 3.2e14` s
  (about 10 Myr). The upper end is deliberately far above any plausible answer
  because the coldest Monte Carlo draws can need hundreds of times the best-fit
  time; a fit that stops at a bound is reported, never silently clipped. The
  scan uses about six points per decade, 60 to 150 points.
* **Free parameters.** `t` always; optionally `x0`, `C_left`, `C_right` and
  `beam_sigma`. Bounds: x0 inside the data range, plateaus within half the
  data range of the observed extremes, beam sigma from 0 to 10 per cent of the
  traverse.
* **`t_guess`.** A Monte Carlo draw passes the base fit's time, which limits
  the scan to a factor of 30 either side and cuts forward solves per draw by
  roughly five, because each draw only perturbs the conditions slightly.
* **`FitResult`** carries the time, the model actually used (with nuisance
  parameters applied), the data, the model profile, the statistics, the scan
  bounds, the solver `route` string and the warnings.
* **Source note.** The scan range is justified in the docstring by the span of
  timescales reviewed in Costa, Shea & Ubide (2020).

#### `fitting/montecarlo.py`
Error propagation, and the reason the project exists.

* **Called by** `gui/workers.MonteCarloWorker`; `contributions` is called from
  the same worker when the variance decomposition box is ticked.
* **`UncertaintyBudget`** holds `sigma_T_K`; `fo2_mode` (`buffer` or
  `absolute`), `buffer`, `delta_buffer`, `sigma_delta_buffer`,
  `sigma_log_fo2`; `sigma_P_Pa`; `sample_coefficient` + `coefficient_mode`;
  `sample_measurement_noise`; `sigma_distance_scale` (relative, e.g. 0.02 for
  a 2 per cent image calibration); `sigma_boundary`; `refit_each_draw`.
  `active_sources()` names the seven possible sources and is what the
  interface, the plot panels and the methods block all read.
* **What one draw does.** `_draw_conditions` samples T and P, then either
  re-evaluates the buffer at the **sampled** temperature (the correlation that
  independent sampling misses) or samples an absolute log fO2. It returns the
  sampled buffer offset too, which becomes the draw's `fo2_buffer`. With a
  cooling path, the **whole path is shifted** by the sampled temperature minus
  the nominal one, keeping its shape; before 29 September 2026 the path stayed
  fixed and sampled temperatures acted only through fO2, so temperature
  uncertainty was largely lost with cooling on.
  `_draw_data` adds Gaussian measurement noise to C and a multiplicative
  scale error to x. `coefficient.sample()` supplies parameter overrides.
  Plateau compositions are perturbed if asked. Then the **whole fit is re-run**
  on the perturbed data with the perturbed conditions, seeded from the base
  time.
* **Parallel draws.** `run(..., workers=N)` fits the draws in `N` processes
  (`default_workers()`: all cores but one). Everything random is drawn first,
  in draw order, in the calling process (`_sample_draw` into `_Draw` records),
  so a given seed gives identical times on any number of workers; only
  `_evaluate` (the fit) runs in the workers. Each worker gets the base model
  once through `_worker_init`, with the coefficient sent as its registry key
  (`_portable`), because many laws are closures that cannot be pickled, then
  batches of draws sized to about half a second (`_worker_batch`). Processes
  are started with `spawn` everywhere, each with one BLAS thread
  (`BLAS_THREAD_VARIABLES`) so the workers do not oversubscribe the cores.
  On a 4-core, 8-thread laptop (i7-10510U) 7 workers fit olivine draws
  about 2.6 times faster than one. A run expected to take less than
  `PARALLEL_THRESHOLD_S = 6` s stays in the calling process, since each worker
  costs about a second to start. If the workers cannot start or die, the
  remaining draws are finished on one core and a note says so. Stopping cancels
  the batches not yet started. `MonteCarloResult.workers` records how many
  processes were used.
* **Reported statistics.** Times are log-normal, so the median and the 16th,
  84th, 2.5th and 97.5th percentiles are reported rather than a symmetric
  sigma; `sigma_log10` is the standard deviation of log10 t. Mutch et al.
  (2021) find the same log-normality.
* **`profiles`** keeps up to `keep_profiles` (default 200) fitted model
  profiles. `envelope(16, 84)` is the band drawn on the profile plot: it is
  the spread of the profiles the Monte Carlo actually fitted, **not** the
  profile at the two ends of the time interval. The distinction is the subject
  of the most recent commit.
* **`on_draw` callback.** Every successful draw reports its time, sampled T,
  sampled log fO2, log10 D, the perturbed data and the fitted profile. This is
  what lets the interface draw the Monte Carlo while it runs.
* **Warnings.** More than 10 per cent of draws failing, and any draw landing
  at the edge of the search range, both attach a warning.
* **`contributions()`** re-runs the Monte Carlo with one source active at a
  time and returns `sigma(log10 t)` for each. It is a one-at-a-time
  sensitivity analysis, so the parts do not add in quadrature to the total; it
  is meant for ranking what to measure better.
* **The argument, with sources.** Linear propagation as used by NIDIS
  (Petrone et al. 2016, Methods) treats `sqrt(4Dt)` and T as independent,
  which captures neither the buffer-temperature correlation nor the
  `ln D0`-`Q` correlation. Sampling from a covariance matrix is the DFENS
  approach (Mutch et al. 2021). Both citations are in the module docstring.

---

### 2.8 `diffusor/dataio/` — in and out

#### `dataio/profiles.py`
Loading and mapping measured profiles.

* **Called by** `gui/main_window.load_file` and `_load`, the dataset loader,
  the tests.
* **`read_table(path)`** dispatches on the suffix: Excel via `pandas.read_excel`,
  `.tsv`/`.tab` with a tab separator, everything else through
  `pandas.read_csv(sep=None, engine="python")`, which sniffs the delimiter.
* **`ProfileSpec`** records the mapping so it can be written into the methods
  block: distance column and unit, columns A and B, their uncertainty columns,
  the sigma level, the composition mode, the oxide names, and an optional
  `x_min`/`x_max` fit window.
* **`suggest_spec(df)`** is the guess the column dialog opens with. Distance is
  the first numeric column whose name tokenises to a distance word, else the
  first numeric column. A known oxide pair (FeO+MgO, CaO+Na2O, FeO+MnO) is
  preferred and modelled as a molar ratio with the oxides named; otherwise a
  cation-name pair; otherwise two bare columns are taken as A and B; otherwise
  the first composition column alone. Uncertainty columns are matched to their
  value column by name stem, and a `2s` token in the name switches the sigma
  level.
* **`build_profile(df, spec)`** converts the distance to micrometres, builds
  the composition variable through `thermo.units.composition_variable`,
  propagates the uncertainty, applies the fit window, drops non-finite rows,
  sorts, and records `row_index` so that `Profile.column(name)` can pull any
  other column of the source table **aligned point for point** with the
  profile. That alignment is what keeps an anorthite column in step after
  windowing and sorting.
* **Uncertainty propagation on a ratio.** For `A/(A+B)`,
  `sigma = sqrt((B sigma_A)^2 + (A sigma_B)^2)/(A+B)^2`, the usual
  first-order quadrature. If only one uncertainty column is given, the other
  element is assumed to carry the same **relative** uncertainty, and a note
  says so. Non-positive or missing sigmas are replaced by the median and
  counted in the notes.

#### `dataio/greyscale.py`
Calibrating BSE grey values to composition.

* **Called by** the greyscale example workflow and the tests.
* **Why it works.** BSE intensity varies with mean atomic number, so for a
  binary substitution such as Fe-Mg it is a linear proxy for composition over
  a narrow range. This is the basis of the grey-value approach of NIDIS
  (Petrone et al. 2016, Methods, "Rationale of working with greyscale values
  of BSE images") and of Morgan et al. (2004). It buys spatial resolution far
  better than a microprobe traverse at the cost of needing a calibration.
* **Rule.** Diffusor never fits raw grey values. `calibrate(grey, comp,
  degree)` fits a polynomial (degree 1 by default) through microprobe anchors
  with `numpy.polyfit(..., cov=True)` when there are degrees of freedom, and
  reports residuals, RMSE and R^2. A degree-1 fit with at least 4 anchors and
  R^2 < 0.95 attaches a note suggesting a non-linear BSE response or a second
  element varying along the traverse.
* **`apply_calibration`** combines the calibration prediction uncertainty with
  the scatter of grey values across averaged raster lines, propagated through
  the local slope of the calibration. NIDIS exports that scatter as the
  standard error of the mean of 200 to 600 lines.
* **`anchors_from_microprobe`** averages the grey profile within `window_um`
  of each probe spot, which should be comparable to the interaction volume so
  the two measurements sample the same material.

#### `dataio/images.py`
Reads any micrograph or element map into one float array, finds the pixel
size if the instrument wrote it, and turns a pixel into one number.

* **Called by** `gui/image_extractor.py` and the tests.
* **`load_image(path, raw=None)`** dispatches on the suffix. PNG, JPEG, BMP,
  GIF, WebP, PNM and TIFF go through Pillow; TIFF goes through `tifffile`
  first **when it is installed** (optional; it reads BigTIFF and compressed
  scientific TIFFs Pillow cannot). Multi-page TIFFs of equal size become
  channels (`page 1`, `page 2`, ...), which is how multi-element map stacks
  are usually saved. `load_text_grid` reads the matrices of counts that
  microprobe software exports (tab, space, comma or semicolon separated,
  decimal commas allowed), dropping header lines, label columns and a leading
  row-number column. `load_envi` reads an ENVI `.hdr` and its binary file.
  `load_raw` reads a headerless binary dump of given width, height, type,
  byte order, header length and band interleave; the extractor offers it for
  any file no other reader recognises. Values are **never rescaled**: a 16-bit
  image keeps 0-65535 and a float map keeps its units; `value_range` records
  the storage type's nominal full scale. Alpha channels are dropped with a note.
* **Pixel size** (`pixel_size_from_tags`, `pixel_size_from_sidecar`): Zeiss
  SmartSEM tag 34118 (`Image Pixel Size = 24.39 nm`), Thermo Fisher/FEI tag
  34682 (`PixelWidth` in m), Tescan tag 50431 (`PixelSizeX` in m), ImageJ
  calibration (description `unit=micron` plus XResolution in pixels per
  unit), a JEOL `.txt` sidecar (`$$SM_MICRON_BAR` px labelled
  `$$SM_MICRON_MARKER`), a Hitachi sidecar (`PixelSize` in nm). With
  `tifffile` the vendor tags are still read through Pillow, because
  `tifffile` pre-parses them into dictionaries. The value is shown to the user
  as a suggestion with its source; plain TIFF resolution tags are ignored
  because most software writes a meaningless 72 dpi.
* **`value_map(image, mode, channel, color_scale)`** returns values, a mask of
  pixels without a value, and a label. `luminance` is Rec. 601
  (0.299 R + 0.587 G + 0.114 B, Pillow's greyscale conversion), `mean` the
  unweighted mean ImageJ uses by default, `channel` one channel or page
  (NIDIS reads only the first channel), `colour_scale` a legend inversion.
* **`ColorScale`** holds an ordered run of legend colours and their values,
  built `from_legend` (sampled along a line the user draws through the legend
  bar, averaged over 3 px across it, values linear or logarithmic between the
  two end values), `from_colormap` (a named matplotlib map with vmin and vmax,
  for when the scheme is known but no legend is in the image) or `from_table`
  (value, R, G, B). `to_values` converts colours to CIELAB (`srgb_to_lab`,
  D65) so "nearest" means "looks most alike", finds the nearest legend colour
  with a k-d tree, refines the value by projecting onto the segment to the
  closer neighbouring legend colour, and gives **no value** to colours farther
  than `max_distance` (Delta E, default 20) from the whole legend: cracks,
  epoxy, labels and the scale bar. Unique colours are looked up once, so a
  full map is quick. A legend with discrete colour steps returns stepped
  values; a cyclic map such as `hsv` is ambiguous at its ends.

#### `dataio/image_profiles.py`
Profiles perpendicular to a guideline drawn on an image, cleaned and
averaged, written to and read back from Excel. The method is that of
`greyvalues.m` in NIDIS (Petrone et al. 2016, code at github.com/cpetrone/NIDIS).

* **Called by** `gui/image_extractor.py` (extraction, workbook),
  `gui/image_calibration.py` (`composition_table`), `gui/main_window.py`
  (`is_extraction_workbook`, `read_extraction`) and the tests.
* **Geometry** (`stations`, `line_geometry`). The guideline is a polyline of
  (x, y) pixel coordinates, x right and y down. Lines are placed every
  `line_spacing_px` (1 px, as NIDIS) along its arc length, each perpendicular
  to the local direction, which is taken over +/- 2 px of guideline so a
  vertex turns the lines gradually. The right-hand normal of direction
  (dx, dy) in image coordinates is (-dy, dx). Each line runs
  `length_before_px` on the right and `length_after_px` on the left (both 50,
  NIDIS's `hp`), sampled every `sample_step_px`. **Distance 0 is the
  right-hand end, looking from the first guideline point to the last**, as
  NIDIS; `flip` starts on the left. `Offset` is distance from the guideline.
* **`extract_profiles(image, settings, color_scale, values)`** samples the
  value map with `scipy.ndimage.map_coordinates`, bilinear (default) or
  nearest (MATLAB `improfile`'s default). Samples off the image are NaN.
* **Cleaning, in the image** (`impurity_masks`): values below `low` (cracks,
  holes) or above `high` (bright inclusions), pixels inside drawn
  `exclusions` polygons, and pixels with no value. Each mask can be grown by
  `grow_px` with a disk, because the edge of a crack is a blend of crack and
  crystal. A sample is flagged if **any** flagged pixel contributes to it
  (the mask is interpolated like the values), so a flagged pixel never leaks
  into a bilinear value. **Across the lines** (`_clip`), at each position:
  `nidis` rejects values outside mean +/- k SD once (k = 1 in NIDIS);
  `mad` (default) rejects outside median +/- k 1.4826 MAD, k = 3, repeated
  up to 5 times; `none`. **Whole lines** with more than
  `max_rejected_fraction` (50 %) of their values rejected are dropped, since a
  lamella lying along a line spoils all of it.
* **Reason codes**: every value keeps why it was not used (`KEPT`, `OUTSIDE`,
  `THRESHOLD`, `EXCLUDED`, `OFF_SCALE`, `OUTLIER`, `LINE_DROPPED`); `REASONS`
  describes them and the workbook records them.
* **`ProfileExtraction.statistics("raw"|"clean")`** gives NIDIS's columns at
  each position: N, min, max, mean, median, SD (n - 1), relative SD and
  SE = SD / sqrt(N). "Raw" is every value inside the image, as NIDIS's
  first CSV; "clean" is what survived. `profile_table` joins both with the
  distances in px and, if a pixel size is set, um.
* **`write_workbook`** writes `Profile` (first, so any reader sees the
  averaged profile; with a chart of raw and cleaned means), `Summary`
  (source, methods text, rejection counts, notes, how to load it, and the
  picture of the lines drawn over the image), `Raw lines` and `Clean lines`
  (every value, NIDIS's `Line_1 ... Line_N` columns), `Rejection codes`,
  `Reason key`, `Geometry` (each line's station and end points, values used,
  dropped or not, plus the guideline vertices), `Colour scale` when a legend
  was used, and `Settings` (a `format` marker and `settings_json`, the whole
  `ExtractionSettings`, so the extraction can be repeated exactly; NIDIS saves
  its guideline coordinates in a `.mat` file for the same purpose). More than
  16380 lines do not fit an Excel sheet and are refused.
* **Back into Diffusor.** `is_extraction_workbook` checks the marker,
  `read_extraction` returns an `ExtractionTable`, `table_from_extraction` the
  same without a file. `composition_table(table, pixel_size_um, statistic,
  uncertainty, calibration, name)` multiplies pixel distances by the pixel
  size, maps the chosen statistic (cleaned mean by default) to composition
  through a `GreyscaleCalibration` with `apply_calibration` (so the
  uncertainty combines the chosen scatter, SE by default, with the
  calibration's own), or passes it through unchanged when `calibration` is
  None (a quantitative map, or a legend already in composition units).
  **Without a pixel size it refuses**, rather than modelling pixels as
  micrometres.

#### `dataio/export.py`
The point of the whole citation machinery.

* **Called by** `gui/main_window.export_results`; `methods_paragraph` is also
  rendered as HTML by `gui/richtext.methods_html`.
* **`collect_citations(fit, mc, extra)`** walks the model and assembles the
  keys this run actually used, in a stable order: always `crank1975` and
  `costa2008`; the coefficient's own citation and its secondary citations;
  `costa_chakraborty2004` if anisotropy was applied; `costa2003`,
  `dohmen2017`, `dohmen_blundy2014` and `grove1984` if the activity term was
  on; `ganguly1988` and `bradshaw_kent2017` if a beam sigma was set;
  `lasaga1983` for a non-isothermal history; `dohmen2017` if the numerical
  solver ran; `frost1991` if fO2 came from a buffer; then always `shea2015`,
  `krimer_costa2017`, `codata2018` and `iupac2021`.
* **`methods_paragraph`** writes the sections Data, Diffusion coefficient
  (the full `describe()`), Conditions, Model, Result, Caveats and References.
  It states the solver route and why, the convolution if used, the Monte Carlo
  budget and sampling mode with a plain-language explanation of each, the
  percentile intervals with the note that times are log-normal, and the
  variance contributions.
* **`result_dict`** is the JSON: version, timestamp, the coefficient with its
  equation text and every parameter, the conditions, the model set-up, the fit
  and its statistics, the warnings, the full citation list, the data spec and
  the Monte Carlo block including the whole budget.
* **`save_results`** writes `*_results.json`, `*_profile.csv` (measured,
  model, residual, sigma, initial condition and the p16/p84 envelope),
  `*_methods.txt`, `*_montecarlo_times.csv` and, if a figure is passed, PNG at
  300 dpi and SVG.

---

### 2.9 `diffusor/datasets.py` — the bundled examples

The catalogue of the eight example datasets, each with an `ExampleDataset`
record: key, name, mineral, species, filename, **kind** (`measured` or
`synthetic`), a full `provenance` string, a citation key for measured data,
the `ProfileSpec` keyword arguments, a `settings` dict that pre-fills the
later steps, an `expected` answer, notes, and a `sources` dict saying where
the T, P and fO2 of the example come from.

* **Called by** `gui/main_window._page_data` (the list), `load_example` and
  `_apply_dataset_settings` (which fills in steps 2 to 6 and logs what it
  set), `gui/richtext.example_html`, and the tests.
* **The `sources` dict is the honest part.** Every entry states whether a
  value is published or a placeholder the chosen coefficient never uses, so a
  placeholder can never be mistaken for a measurement. The tests enforce that
  every dataset has all three.

| key | kind | source of the data |
| --- | --- | --- |
| `plag_santorini` | measured | Druitt et al. (2012) Nature 482:77-80, Supplementary Table 1, crystal S82-30A 12 |
| `opx_kizimen` | measured | Ostorero et al. (2022) Commun. Earth Environ. 3:290, Supplementary Data 2, crystal K9_L10C4 |
| `opx_shinmoedake` | synthetic | forward model, Dias et al. (2025), 950 C, NNO+1, 150 MPa, true time 1.5 yr |
| `cpx_stromboli` | synthetic | forward model, Mueller et al. (2013), 1100 C, conditions after Petrone et al. (2016, 2018), true time 45 d |
| `olivine_laki` | synthetic | forward model, Dohmen & Chakraborty (2007) TaMED, conditions after Hartley et al. (2016), true time 120 d |
| `magnetite_shinmoedake` | synthetic | forward model, Van Orman & Crispin (2010) Table 12, exactly the conditions of Tomiya et al. (2013), true time 8 d |
| `cpx_greyscale` | synthetic | forward model plus a linear grey response, workflow after Petrone et al. (2016), true time 3 yr |
| `sanidine_ba` | synthetic | forward model, Cherniak (2002), 790 C inside the Bishop Tuff range of Chamberlain et al. (2014), true time 5 kyr |

* **Kizimen is the validation.** Ostorero et al. (2022) modelled the reverse
  zone at 850 C with Ganguly & Tazzoli (1994) and no fO2 term and obtained
  2.32 yr (+7.16/-1.75). Diffusor fits about 3 yr with the same law and
  temperature, and `tests/test_datasets_and_ui.py` checks it. Its conditions
  are the mean and standard deviation of 21 magnetite-ilmenite pairs from the
  andesites, Supplementary Data 3: 850 +/- 57 C and NNO +1.28 +/- 0.35. The
  record notes that the main text of that paper swaps the andesite and dacite
  temperatures while the table does not.
* **Santorini is not a validation** and says so: Druitt et al. reconstructed
  their initial profile from the Sr-An correlation and a two-melt history,
  which Diffusor's built-in initial conditions cannot express.

---

### 2.10 `diffusor/gui/` — the desktop application

Six steps and a results view. Pages never scroll; only lists and reading panes
do, and scrolling never changes a number (`widgets.NoWheelOnInputs`).

#### `gui/main_window.py` (2003 lines)
The whole flow. Everything the user chooses ends up in `_model()`,
`_free_parameters()` and `_budget()`, which are the three functions the rest
of the library sees.

**Module-level tables.**
* `STEPS = ["Data", "Mineral", "Conditions", "Model", "Coefficient",
  "Uncertainty", "Results"]`.
* `RESOLUTION_PRESETS`: label, width kind, default width, fixed sigma, hint.
  The conversions are `sigma = d/4` for an evenly lit round spot and
  `sigma = w/sqrt(12)` for an evenly lit slit. The presets and the studies
  behind them: BSE grey-value profiles resolve better than 0.5 um
  (Petrone et al. 2016); a focused microprobe beam rarely exceeds
  sigma = 0.6 um (Ganguly et al. 1988); a 5 um defocused beam is common on
  feldspar (Chamberlain et al. 2014; Grocolas et al. 2025); a 10 um laser spot
  and 10 to 15 um ion beams (Druitt et al. 2012); a 7.5 um line-scan slit
  (Grocolas et al. 2025).
* `BOUNDARIES`: the four plain-language ends (`far`, `rim_melt`,
  `rim_closed`, `centre`).
* `SAMPLING_MODES`: the four coefficient-sampling options shown on step 6.
* `SOURCE_NAMES`: display names for the seven Monte Carlo sources.

**Step by step.**
1. **Data.** `load_file` reads the table, `suggest_spec` guesses, a
   `ColumnDialog` confirms, `build_profile` applies the mapping. `load_example`
   does the same for a bundled dataset and then `_apply_dataset_settings`
   fills steps 2 to 6 and shows a prefill bar on each saying what it set.
2. **Mineral.** Mineral, diffusing species, traverse direction (the
   experiments' own axis, a named axis, or three angles which are passed to
   the Costa & Chakraborty 2004 relation), and the representative composition
   plus the "D follows the composition along the profile" switch.
3. **Conditions.** T and P with 1-sigma values, fO2 as a buffer offset or an
   absolute log fO2, and the analytical resolution preset. `_update_fo2_label`
   shows the resulting absolute log fO2 in both bar and Pa, naming Frost
   (1991) as the buffer source, and re-runs whenever T, P, the buffer or the
   offset changes.
4. **Model.** Geometry, initial condition (with the **Guess** button calling
   `guess_step_from_data`), what each end of the profile is, the solver
   override and the node count. `_update_solver_note` calls
   `can_use_analytical()` on a trial model and states which solver will run
   and why; the node spinner is disabled when the closed form will be used.
5. **Coefficient.** `_refresh_coefficients` lists everything registered for
   that mineral and species, tagged recommended, unverified or superseded.
   Ticking two or more enables Compare.
6. **Uncertainty.** Draws, seed, processor cores (`sp_cores`, default all but
   one), which sources to sample, the coefficient sampling mode and whether to
   run the variance decomposition.
   `_on_dmode_changed` greys out a mode the chosen law cannot support and
   names the one that will actually be used.
7. **Results.** A narrow summary of every setting with an *edit* link beside
   each group, the plot, and the Methods, coefficient and log reading panes.

**The three functions that matter.**
* `_conditions(coef)` builds `Conditions`: T in K, P in Pa, absolute
  `log_fo2_bar` (through `log_fo2_from_delta` when the buffer mode is on), an
  `X` dict keyed by the mineral's composition variable plus anything the
  coefficient `requires`, and the axis or angles.
* `_model(coef_key)` assembles the `DiffusionModel`. It chooses between a step
  and the `equilibrium_plag` initial condition, maps the two boundary
  drop-downs onto Dirichlet or zero flux, builds a linear `ThermalHistory` if
  the cooling box is ticked, sets `comp_key` only when the coefficient
  actually requires that variable, switches on the plagioclase activity term
  when the mineral is plagioclase, the species has an `ACTIVITY_A` value and
  an anorthite column was loaded, and lays a uniform `x_grid` of `n_nodes`
  across the data span with the anorthite profile interpolated onto it.
* `_budget()` reads the six Monte Carlo check boxes and spinners into an
  `UncertaintyBudget`.

**Running.** `run_fit`, `run_compare` and `run_mc` each build the model,
create a worker and hand it to `_launch`, which refuses to start a second job
while one is running. Worker signals are connected only to `@Slot` methods of
the window, never to lambdas, so Qt delivers them on the interface thread.

**Updates.** `_build_update_banner` adds a hidden strip under the header:
"Diffusor X is available. This is Y." with **Download** (opens the release page
in the browser, then hides the strip) and **Later** (hides it; nothing is
remembered, so the next check shows it again). `check_for_updates(manual)`
starts an `UpdateWorker` kept in `_update_job`, separate from `_jobs` so an
update check never blocks or disables Fit, Compare and Monte Carlo.
`_update_found` records the time in `QSettings("Diffusor", "Diffusor")` under
`updates/last_check` (only on success, so an offline start is retried at the
next launch) and shows the strip when `updates.is_newer`; the release notes are
its tooltip. A check the user asked for (Help > Check for updates...) always
answers, with a message box saying "newest version" or why the check failed; a
startup check is silent apart from a line in the log. `startup_update_check`
runs the silent check unless Help > Check for updates at startup is unticked
(`updates/check_at_startup`, default on) or the last success was under 20 hours
ago. `closeEvent` waits for a running check, which ends within its own 6 s
timeout.

#### `gui/image_extractor.py`
`ImageExtractorDialog`, the nonmodal window that draws profiles on images;
opened by `MainWindow.show_image_extractor` from File > Extract profile from
image, from "From an image..." on the Data step, and when an image file is
picked in File > Load profile (`IMAGE_ONLY_SUFFIXES`). The image fills the
left with a matplotlib toolbar for zoom and pan; the averaged profile (raw
mean, cleaned mean +/- 1 SD, the guideline as a dashed line) sits under it;
the settings are on the right in cards: Image (file, value mode, channel),
Draw (tools), Scale, Colour legend (only in colour-scale mode), Profile lines,
Cleaning and Result.

* **Tools** (mouse, only while the plot toolbar's pan and zoom are off):
  Boundary adds guideline points, right-click or Backspace removes the last;
  Scale bar takes two clicks and asks for the bar's length, setting um per
  px; Legend takes the low-value end then the high end of a colour legend;
  Exclude area adds polygon vertices, right-click or Enter closes it.
* **Legends.** The Colour legend card offers three sources: drawn on this
  image (the Legend tool), on a separate image, or a named matplotlib map.
  Microprobe software often saves the legend as its own file, so "On a
  separate image" opens that file in `LegendPickerDialog`, where the low and
  high ends are clicked (right-click undoes, the toolbar zooms on a small
  legend; the button is enabled only once both ends are set);
  `set_legend_image` keeps that image and its two points, and "Pick the ends
  again" reopens it. The legend must be a colour image.
* **Every change re-extracts** after a 120 ms debounce (`_schedule` →
  `recompute`), so what is saved is what is shown. The value map is computed
  once per image, mode and legend (`_values_changed`) and passed to
  `extract_profiles`, so moving a slider does not redo a colour inversion.
  Up to 40 of the lines, the start-side and far-side ends, the rejected
  samples (red, at most 20000 drawn) and the exclusion areas are drawn over
  the image. A colour map is shown in colour so the legend can be found;
  anything else is shown as the grey values the profile reads.
* **Opening a second image of the same size keeps the guideline and
  exclusion areas**, so a Mg map and a Fe map of the same area give profiles
  along the same line. "Reuse settings from a workbook" loads every setting
  from an earlier extraction (NIDIS's "use existing coordinates").
* An unrecognised file offers `RawDialog` (width, height, type, byte order,
  header bytes, channels) and reads it with `load_raw`.
* All three windows are sized with `widgets.fit_to_screen`, and the settings
  column scrolls, so nothing is off screen on a small laptop display.
* `save_workbook` proposes `<image>_profile<n>.xlsx` beside the image, as
  NIDIS names its output folders, and embeds `overlay_png()`.
  `use_in_diffusor` hands `table_from_extraction` to the main window.

#### `gui/image_calibration.py`
`ImageCalibrationDialog`: what the image cannot tell Diffusor. Pixel size
(prefilled from the workbook); which statistic (cleaned mean, cleaned median,
raw mean) and which uncertainty (SE as NIDIS, SD across lines, none); and the
value-to-composition map: linear through two reference points, fitted to
anchor points (linear or quadratic, typed in or filled from a microprobe
traverse with `anchors_from_microprobe` and an averaging window), or none.
The calibration's `describe()` and notes show live. A profile read through a
colour legend starts on "none", since legends are usually drawn in
composition units. The form scrolls and the buttons stay under it, so the
dialog fits any screen. On accept it builds the
table with `composition_table` and a `ProfileSpec` (mode A, distance in um),
and `MainWindow.load_image_table` passes both to `_use_table`, the part of
loading shared with ordinary files, after logging `summary()`.

#### `gui/workers.py`
`FitWorker`, `CompareWorker`, `MonteCarloWorker` and `UpdateWorker`, each a `QObject` moved
onto a fresh `QThread` by `start()`. `CompareWorker` fits the same profile with
several coefficients and reports progress per coefficient.
`MonteCarloWorker` buffers the `on_draw` callbacks and emits them in batches
at most once a second (`BATCH_SECONDS = 1.0`), so the live plot can grow
without flooding the event loop, and runs the variance decomposition
afterwards if asked. `UpdateWorker` calls `updates.fetch_latest` and emits the
`ReleaseInfo` or the failure sentence; it has no `abort()`, because a web
request cannot be interrupted, and relies on the 6 s timeout. Both long
workers support `abort()`. `MonteCarloWorker`
passes the core count from the Uncertainty step to `run` and `contributions`,
which fit the draws in worker processes; the `QThread` only waits for them.

#### `gui/plot_widget.py`
The matplotlib canvas, with `SaveToolbar`: matplotlib's toolbar whose save
dialog opens in `start_dir` with `default_name` (the main window sets both to
the loaded profile's folder and `<profile>_diffusor.png`; matplotlib on its own
starts wherever the last figure went). `show_data`, `show_fit` (profile plus residual panel),
`show_comparison`, `show_histogram`, and the live Monte Carlo view:
`start_monte_carlo` lays out the panels, `add_monte_carlo_draws` adds each
batch as a faint fitted curve plus the perturbed points (capped at 250
curves), `finish_monte_carlo` draws the final band and the best curve.
`_redraw_mc_panels` fills the histogram of times with the running median and
68 per cent interval, and the sampled temperature, T against log fO2 (which
shows the buffer correlation directly) and sampled log D. A panel whose input
was held fixed says so instead of showing an empty axis.

#### `gui/richtext.py`
Formatted HTML for the reading panes, using the subset of HTML 4 and CSS that
Qt renders reliably. `Doc` is a small builder (`h`, `p`, `ul`, `kv`, `table`,
`box`, `equation`). `coefficient_html` renders a law with its equation,
parameters, ranges, verification status and notes; `example_html` renders a
dataset with its provenance; `methods_html` mirrors
`dataio.export.methods_paragraph`; `compare_html`, `mc_html`, `log_html`,
`format_html` and `boundaries_html` cover the rest. `reference_items` turns
citation keys into formatted entries.

#### `gui/widgets.py`, `gui/theme.py`, `gui/format_help.py`, `gui/icons/`
Layout helpers (`card`, `field`, `row`, `pair`, `callout`, `collapsible`,
`page_columns`, `WrapLabel`, which reserves the height a wrapped label
actually needs (0 while it is empty: Qt reports -1, and a -1 minimum height
prints a "Negative sizes" warning), and `fit_to_screen`, which opens a window at its preferred
size or smaller so it fits the available screen with room for the title bar);
the palette and the Qt stylesheet, with the plot colours kept
in step with the interface; and the one place that describes what an input
file has to look like, shared by the hint text and the Format dialog.

---

### 2.11 `scripts/` — maintenance, not runtime

Nothing in `diffusor/` imports these. They are run by hand from the project
root.

#### `scripts/build_references.py`
Regenerates `references.bib` and `REFERENCES.md` from `diffusor/references.py`
plus `coefficients.list_coefficients()` and `thermo.buffers.BUFFER_CITATIONS`.
The Markdown file gets a table of every coefficient with its source and
verification status, the buffer citations, and the full bibliography annotated
with which coefficients use each key (primary or secondary).
**Run it after any change to the registry**, or the two generated files drift.

#### `scripts/make_examples.py`
Regenerates the seven synthetic files in `examples/` and their README table.
Each is a forward model run at a known time with `numpy` Gaussian noise added
and a fixed seed (`default_rng(20260915)`), so the whole chain can be checked
against an exactly known answer. `x_fe_to_oxides` converts X_Fe back to FeO
and MgO wt% for a pyroxene stoichiometry so the files look like real
microprobe output. The true times are 1.5 yr (opx), 45 d (cpx), 20 yr (plag
Mg), 120 d (olivine), 8 d (magnetite), 3 yr (greyscale) and 5 kyr (sanidine
Ba). **Re-running it changes the files, so the round-trip test tolerances
should be re-checked afterwards.**

#### `scripts/make_example_images.py`
Writes the two SYNTHETIC images in `examples/images/` for trying the image
extractor (`default_rng(20260929)`). `cpx_bse_zoned.tif` is 16-bit with an
ImageJ calibration of 0.05 um/px, a curved Fe-rich rim with X_Fe 0.16 to
0.26 as grey = 60 + 400 X_Fe (8-bit scale, x 257), an error-function
boundary of half-width 0.8 um 6 um inside the crystal face, plus a crack, an
oxide inclusion and a darker lamella. `opx_mg_map_jet.png` is a 'jet' MgO
map with a drawn legend (16 to 30 wt%) and a 20 um scale bar of 40 px, but no
pixel size in its metadata, as a microprobe export would be.

#### `scripts/make_icon.py`
Draws the application icon into `diffusor/gui/icons/`: a teal tile in the
interface accent colour, a crystal shaded across a gold zone boundary by an
error function, and the error-function profile across it in white. Drawn at
1024 px (the curve stamped with a round brush so it stays smooth) and reduced
to `diffusor.png` (256 px, the window icon) and `diffusor.ico` (16 to 256 px).

#### `scripts/extract_kizimen.py`
Rebuilds `examples/opx_kizimen_ostorero2022.csv` from the published
spreadsheet. Reads Supplementary Data 2, sheet "Opx profiles - compositions
And", takes crystal `K9_L10C4`, starts at the first row the authors coloured
green (`FF92D050`, their good-quality flag) and keeps the analytical standard
deviations from row 3. Nothing is smoothed or re-scaled. Usage:
`python scripts/extract_kizimen.py <folder with the supplementary files>`.
The files are archived at doi:10.5281/zenodo.7307563.

---

### 2.12 `tests/` — what is actually guaranteed

The suite is part of the traceability argument: it fails when a transcription
stops reproducing the number printed in its source.

| file | what it holds the code to |
| --- | --- |
| `test_coefficients.py` | Every coefficient reproduces a number printed in its source: Dohmen et al. (2016) run OPXD_14 and the 3.5 anisotropy, the Ganguly & Tazzoli form of Ostorero et al. (2022) eq. 1, the Mueller et al. (2013) abstract, the Dimanov & Sautter values quoted in Petrone et al. (2016) Table 2, Schwandt et al. (1998) Table 3, Van Orman et al. (2014) eq. 4, the course-script form of Giletti & Casserly (1994), the Tomiya et al. (2013) Shinmoedake diffusivity, the magnetite minimum and its 2/3 exponents, the olivine 1/6 fO2 exponent and 6x anisotropy. Also that unverified and out-of-range entries produce warnings, and that direction cosines must close. |
| `test_recent_literature.py` | The 2020-2026 calibrations: Grocolas et al. (2025) eqs 7, 8 and 12-14 and their covariance assumption, the Sr gap against the 1990s laws and how it widens as T falls, Audetat et al. (2026) eq. 1 and its silica-activity term, Pohl et al. (2024) Li, Dias & Dohmen (2024) m(T) against their Table 1, Dias et al. (2025) eq. 22 by hand and the regime switch, the opx REE laws, and the Sievwright et al. (2020) Table 5 transcription against the published minima. |
| `test_solvers.py` | The numerical solver against Crank (1975) eqs 2.14, 4.17, 5.22 and 6.18; second-order convergence; mass conservation in a closed system, including with the radial weight; explicit against implicit; the fractional-uptake series; the non-isothermal integral; the convolution; the resolution warning. |
| `test_fitting.py` | A known time is recovered with and without noise and with composition-dependent D; free x0 is recovered; ignoring the beam convolution lengthens the apparent time; the Monte Carlo spread grows with the temperature uncertainty and vanishes with a zero budget; **T and fO2 stay correlated through the buffer**; independent sampling overestimates relative to logD_at_T; the variance decomposition ranks temperature first for opx. |
| `test_thermo_and_io.py` | Buffer ordering and the two parameterisations agreeing; oxide molar masses; molar (not weight) ratios; time conversions; every reference rendering and producing BibTeX; the column guesser; sigma propagation and 2-sigma halving; the greyscale calibration; and that the methods paragraph lists the sources actually used. |
| `test_examples_roundtrip.py` | Each synthetic example is loaded, fitted and checked against its known time, plus the greyscale calibrate-then-fit path and the far-field warning. |
| `test_datasets_and_ui.py` | Every dataset declares its provenance and never calls synthetic data measured; superseded entries are flagged and demoted; the window builds with one page per step; loading an example fills in the later steps; the Kizimen fit lands inside the Ostorero uncertainty; no label on any step is clipped and no page scrolls. |
| `test_kfeldspar_and_interface.py` | The K-feldspar laws against their abstracts and the Sr-Ba gap; only the Grocolas entries carry a covariance; every example says where T, P and fO2 come from; boundary choices reach the model and closed ends hold the mass in; the Monte Carlo reports every draw and the worker batches them; the band on the profile is the spread of the refitted draws. |
| `test_image_profiles.py` | Lines start on the right of the guideline as in NIDIS and recover a known error function; lines stay perpendicular to a curved guideline; samples off the image are not used; value limits with grown masks remove a crack and an inclusion; NIDIS's 1 SD test rejects about a third of clean Gaussian data and the MAD test almost none; exclusion polygons and whole-line dropping; PNG, BMP, JPEG, 8/16-bit and float TIFF, multi-page TIFF, NumPy, text grids with headers, ENVI and raw binary are read; ImageJ, FEI and JEOL pixel sizes are found; a jet map is inverted to under 1 % of its range with off-scale pixels flagged, from a named map and from a legend drawn in the image; the workbook round trip and the calibrated profile; the extractor window saves and loads into Diffusor through both paths; a legend saved as a separate image is picked and used; the image windows fit the screen and their content can shrink. |
| `test_updates.py` | Versions compare as numbers (`0.10.0 > 0.9.0`, `1.0 == 1.0.0`) and a pre-release or unreadable tag never triggers a notice; the version is written only in `diffusor/__init__.py`; a release link that does not point at this project's own release pages is replaced; every failure (404, rate limit, offline, timeout, bad JSON, oversized or non-UTF-8 reply) becomes a readable `UpdateError`; the request carries a User-Agent; the worker reports a release or a reason; the banner appears only for a newer release, **Later** hides it and **Download** opens the release page; a silent check never opens a dialog and a manual one always does; only a successful check is remembered; the startup switch and the daily limit are respected and the switch survives a restart. The network is faked and the settings live in a temporary file, so the suite never goes online or touches the real settings. |
| `test_writeup.py` | This document lists every module, coefficient key, citation key, dataset key and buffer that exists in the code. |

---

## 3. Papers

### 3.1 Scope of this section

All 80 keys in `diffusor/references.py`, each with a short overview and an
exact list of what Diffusor takes from it. Grouped by the role the source
plays in the code rather than by the order they sit in the file; every key
appears exactly once below.

"**Used for**" lists only what is actually taken. Where a paper is in the
registry but nothing in the code points at its key, it is listed in
section 3.7 instead, with a note on where it is cited in prose.

---

### 3.2 Textbooks and reviews

**`crank1975`** — Crank, J. (1975) *The Mathematics of Diffusion*, 2nd ed.,
Oxford University Press, 414 pp.
The standard reference for solutions of the diffusion equation: closed forms
for every simple geometry and initial condition, plus the finite-difference
chapter. It is the mathematical backbone of Diffusor and the most-cited key in
the package.
*Used for:* the governing equation in plane, cylindrical and spherical form
(eq. 1.7, 5.4, 6.3); the step in an infinite medium (eq. 2.14, the only
analytical route the fit uses); the semi-infinite fixed surface (eq. 2.45);
the band (eq. 2.15); plane sheet, cylinder and sphere (eqs 4.17, 5.22, 6.18)
with their fractional-uptake series (eqs 4.18, 5.23, 6.20); the time-dependent
D substitution (section 7.2, eq. 7.7); the explicit scheme (eq. 8.31), the
stability limit (eq. 8.33), Crank-Nicolson (section 8.4, eq. 8.35) and the
L'Hopital limit at r = 0 (section 8.5, eq. 8.45); the Arrhenius pressure term
(section 11); boundary condition types (section 1.3).

**`costa2008`** — Costa, Dohmen & Chakraborty (2008) Time scales of magmatic
processes from modelling the zoning patterns of crystals. *RiMG* 69:545-594.
The review that made diffusion chronometry a standard tool: it sets out the
workflow, the choice of initial and boundary conditions, and the error sources.
*Used for:* the open-boundary (fixed-concentration) description, p. 555; the
observation that a closed system equilibrates far faster than an open one
(Fig. 6); the step solution as their eq. 9; the Arrhenius form with an
activation volume, eq. 20; the olivine activation volume of 7e-6 m3/mol
quoted on p. 571; the warning that 1-D modelling of a 3-D crystal
over-estimates the time. It is added to **every** exported citation list.

**`dohmen2017`** — Dohmen, Faak & Blundy (2017) Chronometry and speedometry of
magmatic processes using chemical diffusion in olivine, plagioclase and
pyroxenes. *RiMG* 83:535-575.
The successor review, with an appendix that gives the finite-difference
scheme and the plagioclase activity formalism explicitly enough to implement.
*Used for:* the non-ideality flux term (App. eqs A7-A8); the activity slopes
A_i for Mg, Sr, Ba, Li, K and Rb at 1200 C (App. Fig. A1 caption and Table 1,
with 900 C values also tabulated there) through `-R T ln(gamma_i) = A_i X_An +
B_i` (App. eq. A6); the quasi-steady-state profile `C = C0 exp(A X_An/RT)`
(App. eq. A13) and its rim condition (eq. A14); the half-node discretisation
(App. eqs A17-A19), the activity-term coefficients (eq. A20) and the
theta-scheme (eq. A21).

**`costa2020`** — Costa, Shea & Ubide (2020) Diffusion chronometry and the
timescales of magmatic processes. *Nat. Rev. Earth Environ.* 1:201-214.
A shorter, more recent review aimed at what the method can and cannot resolve.
*Used for:* justifying the default time search window, which has to bracket
everything from syn-eruptive ascent to long crustal residence
(`fitting/fit.py` docstring). Prose only, no key reference in code.

**`chakraborty2010`** — Chakraborty (2010) Diffusion coefficients in olivine,
wadsleyite and ringwoodite. *RiMG* 72:603-639.
The compilation of olivine diffusion data, including the point-defect
framework behind the TaMED and PED regimes.
*Used for:* a secondary citation on the olivine TaMED entry, and one of the
secondary sources the unverified olivine transcription was checked against.

**`cherniak2010`** — Cherniak (2010) Cation diffusion in feldspars. *RiMG*
72:691-733.
The feldspar compilation, covering Sr, Ba, Pb, alkalis and the coupled
NaSi-CaAl exchange.
*Used for:* a secondary citation on the Giletti & Casserly (1994) Sr entry.

**`cherniak_dimanov2010`** — Cherniak & Dimanov (2010) Diffusion in pyroxene,
mica and amphibole. *RiMG* 72:641-690.
The pyroxene compilation.
*Used for:* the statement in `minerals/definitions.py` that Fe-Mg anisotropy
in clinopyroxene is not well constrained. Prose only.

**`vanorman_crispin2010`** — Van Orman & Crispin (2010) Diffusion in oxides.
*RiMG* 72:757-825.
The oxide compilation, and the source Diffusor leans on hardest for magnetite:
Tables 10 and 11 give Arrhenius fits along fixed buffers, Table 11 also lists
the Fe-Ti interdiffusion laws, and Table 12 gives the vacancy plus
interstitial parameterisation as a function of T and fO2.
*Used for:* **14 registry entries**. Table 12 (D_V0, Q_V, D_I0, Q_I for Cr,
Al, Fe, Co, Mn, Ti at x_Ti = 0 and Fe, Co, Mn, Ti at x_Ti = 0.2) and the
vacancy/interstitial sum in its footnote; Tables 10 and 11 for the eight
buffer-specific Arrhenius entries; Table 11 for the Freer & Hauptman (1978)
and Aragon et al. (1984) Fe-Ti laws; the `x_Ti = X_Usp/3` convention; the
activation energies borrowed to give the Sievwright entries a temperature
dependence.

**`lasaga1983`** — Lasaga (1983) Geospeedometry: an extension of
geothermometry. In *Kinetics and Equilibrium in Mineral Reactions*, 81-114.
The paper that established reading cooling rates from diffusion profiles, and
the use of the time-integrated diffusion coefficient.
*Used for:* the effective-Dt integral in `solvers/history.py`; added to the
citation list of any run with a non-isothermal history.

**`deer1992`** — Deer, Howie & Zussman (1992) *An Introduction to the
Rock-Forming Minerals*, 2nd ed., Longman.
The standard mineralogy text.
*Used for:* the composition conventions X_Fe, X_Fo, X_An and X_Or, and the
`citation` field of every `CompositionVariable` except x_Ti.

**`dmg2025`** — DMG Short Course: Diffusion in minerals, Ruhr-Universität
Bochum, October 2025. Unpublished course material (lectures, practicals and
the scripts `Diffusion_equation_for_plag.pdf`,
`Diff_Model_Sr_in_Plag_implicit.m`, `MCdiff_OlFo_MO.m`).
Course material from Dohmen, Chakraborty and colleagues, with working MATLAB
implementations of the plagioclase and olivine problems.
*Used for:* the Courant numbers of 0.2 to 0.4 noted in `stability_dt`; the
smoothing of a sharp initial step; the An-dependent form of the Giletti &
Casserly (1994) Sr law, transcribed verbatim from
`Diff_Model_Sr_in_Plag_implicit.m` and the source of the
`plag_Sr_giletti_casserly1994` entry; the direction-cosine relation as taught
in Lecture 6; the convolution correction in Practical 5; the note that the
course script uses A_Sr = -15.1 kJ/mol rather than the -17.4 in `ACTIVITY_A`.
**This is the only unpublished source in the registry**, and one registry
entry rests on it.

---

### 3.3 Constants, units and oxygen buffers

**`codata2018`** — Tiesinga, Mohr, Newell & Taylor (2021) CODATA recommended
values of the fundamental physical constants: 2018. *Rev. Mod. Phys.*
93:025010.
*Used for:* R = 8.314462618 J/mol/K, k_B, N_A and the eV conversion in
`constants.py`. Added to every exported citation list.

**`iupac2021`** — Prohaska, Irrgeher, Benefield et al. (2022) Standard atomic
weights of the elements 2021. *Pure Appl. Chem.* 94:573-600.
*Used for:* all 16 atomic masses and hence the 15 oxide molar masses, which
are what converts wt% oxide to cation moles. Added to every exported citation
list.

**`iau_julian_year`** — IAU (1976) Resolution on the Julian year of 365.25
days.
*Used for:* `SEC_PER_YEAR = 365.25 x 86400 = 3.15576e7` s, the unit every
reported timescale in years is divided by.

**`frost1991`** — Frost (1991) Introduction to oxygen fugacity and its
petrologic importance. *Reviews in Mineralogy* 25:1-9.
The standard short introduction, whose Table 1 compiles the solid buffers in
the compact `A/T + B + C(P-1)/T` form.
*Used for:* the default buffer parameterisation and all five buffers (IW, WM,
FMQ, NNO, HM) with their pressure terms; the pressure term is also added to
the O'Neill parameterisations. Added to the citation list of any run whose
fO2 came from a buffer.

**`oneill_pownceby1993`** — O'Neill & Pownceby (1993) Thermodynamic data from
redox reactions at high temperatures I. *CMP* 114:296-314.
Electrochemical measurements with revised free energies for the Fe-FeO,
Co-CoO, Ni-NiO and Cu-Cu2O buffers.
*Used for:* the NNO free-energy expression in the `oneill` parameterisation.

**`oneill1987`** — O'Neill (1987) Quartz-fayalite-iron and
quartz-fayalite-magnetite equilibria and the free energy of formation of
fayalite and magnetite. *Am. Mineral.* 72:67-75.
*Used for:* the FMQ free-energy expression in the `oneill` parameterisation.

**`huebner1971`** — Huebner (1971) Buffering techniques for hydrostatic
systems at elevated pressures. In *Research Techniques for High Pressure and
High Temperature*, 123-177.
The classic buffer-technique chapter, including the buffer equations used to
convert fO2-dependent diffusion data onto a fixed buffer.
*Used for:* named in `magnetite.py` as the buffer equations Van Orman &
Crispin (2010) used to compute their Tables 10 and 11. Prose only; no key
reference in code.

---

### 3.4 Method papers

**`costa2003`** — Costa, Chakraborty & Dohmen (2003) Diffusion coupling
between trace and major elements and a model for calculation of magma
residence times using plagioclase. *GCA* 67:2189-2200.
Showed that trace-element diffusion in plagioclase is driven by the
chemical-potential gradient, so the anorthite gradient enters the flux
directly, and gave the first Mg-in-plagioclase parameterisation.
*Used for:* the activity flux term, their eq. 7, implemented in the numerical
solver; the `plag_Mg_costa2003` coefficient entry (in the re-written form of
Van Orman et al. 2014); the reason the analytical route is refused when the
activity term is on. Added to the citation list of any run using it.

**`costa_chakraborty2004`** — Costa & Chakraborty (2004) Decadal time gaps
between mafic intrusion and silicic eruption obtained from chemical zoning
patterns in olivine. *EPSL* 227:517-530.
*Used for:* the direction-cosine relation
`D_V = D_a cos^2 a + D_b cos^2 b + D_c cos^2 g`, applied whenever a traverse
is given as three angles. Added to the citation list when anisotropy is used.

**`dohmen_blundy2014`** — Dohmen & Blundy (2014) A predictive thermodynamic
model for element partitioning between plagioclase and melt. *Am. J. Sci.*
314:1319-1372.
A lattice-strain model for plagioclase-melt partitioning as a function of P,
T and composition.
*Used for:* the model from which Dohmen et al. (2017) derived the activity
slopes A_i that Diffusor uses. Added to the citation list when the activity
term is on.

**`ganguly1988`** — Ganguly, Bhattacharya & Chakraborty (1988) Convolution
effect in the determination of compositional profiles and diffusion
coefficients by microprobe step scans. *Am. Mineral.* 73:901-909.
Established that a microbeam measures a concentration averaged over its
interaction volume, and that the correct treatment is to convolve the model
rather than deconvolve the data.
*Used for:* the Gaussian convolution in `solvers/convolution.py`; the
"sigma rarely exceeds 0.6 um on a modern microprobe" preset; the definition of
the custom sigma as the standard deviation of the beam profile. Added to the
citation list whenever a beam sigma is set.

**`bradshaw_kent2017`** — Bradshaw & Kent (2017) The analytical limits of
modeling short diffusion timescales. *Chem. Geol.* 466:667-677.
Quantified how short a timescale can be retrieved before the analytical
resolution dominates.
*Used for:* the rule that a profile with `2 sqrt(Dt) < 3 sigma` should be
reported as an upper bound only, which is the `resolution_warning` text.

**`shea2015`** — Shea, Lynn & Garcia (2015) Cracking the olivine zoning code:
distinguishing between crystal growth and diffusion. *Geology* 43:935-938.
Showed with 3-D models how much of an apparent diffusion profile can be
growth zoning, and how sectioning biases a 1-D interpretation.
*Used for:* the standing caveat attached to every plane-geometry fit, and the
warning in `solvers/initial.py` that a mixed growth-and-diffusion profile
gives a spuriously long time. Added to every exported citation list.

**`krimer_costa2017`** — Krimer & Costa (2017) Evaluation of the effects of 3D
diffusion, crystal geometry, and initial conditions on retrieved time-scales
from Fe-Mg zoning in natural oriented orthopyroxene crystals. *GCA*
196:271-288.
The quantitative study of how much a 1-D treatment of a 3-D orthopyroxene
biases the answer.
*Used for:* the same standing caveat. Added to every exported citation list.

**`mutch2021`** — Mutch, Maclennan & Madden-Nadeau (2021) DFENS: diffusion
chronometry using finite elements and nested sampling. *G-cubed*
22:e2020GC009303.
A Bayesian diffusion-chronometry code that samples the diffusion-law
parameters from their covariance matrix and finds times to be log-normally
distributed.
*Used for:* the justification for the `covariance` sampling mode; the
log-normal reporting convention; the statement that the DFENS source uses an
olivine anisotropy of 6.0. Prose only; no key reference in code.

**`girona_costa2013`** — Girona & Costa (2013) DIPRA: a user-friendly program
to model multi-element diffusion in olivine. *G-cubed* 14:422-431.
An earlier desktop tool for olivine diffusion chronometry.
*Used for:* comparison context only. Not referenced by key or prose.

**`morgan2004`** — Morgan, Blake, Rogers et al. (2004) Time scales of crystal
residence and magma chamber volume from modelling of diffusion profiles in
phenocrysts: Vesuvius 1944. *EPSL* 222:933-946.
One of the first studies to use calibrated BSE grey-value profiles for
diffusion chronometry.
*Used for:* named in `dataio/greyscale.py` as a source of the grey-value
approach. Prose only.

**`petrone2016`** — Petrone, Bugatti, Braschi & Tommasini (2016) Pre-eruptive
magmatic processes re-timed using a non-isothermal approach to magma chamber
dynamics. *Nat. Commun.* 7:12946.
Introduced NIDIS, a non-isothermal diffusion approach applied to Stromboli
clinopyroxene, working from BSE grey-value profiles.
*Used for:* the grey-value calibration workflow; the statement that NIDIS fits
the error-function step with free x0 and 2 sqrt(Dt); the Dimanov & Sautter
(2000) D0 and dH quoted in its Table 2 footnote, which is the whole
`cpx_FeMg_dimanov_sautter2000` entry and the test that checks it; the "BSE
profiles resolve better than 0.5 um" resolution preset; and, throughout the
uncertainty documentation, as the example of propagating errors as if T and
sqrt(4Dt) were independent. Its `greyvalues.m` (github.com/cpetrone/NIDIS) is
the model for `dataio/image_profiles.py`: a guideline along the feature, one
perpendicular line per guideline pixel, 50 px either side, the first value to
the right of the guideline looking from its first point, averaging position
by position, the mean +/- 1 SD outlier test offered as the `nidis` preset, raw
and cleaned outputs with min, max, mean, SD, relative SD and SE, and saving
the guideline so the extraction can be repeated.

**`bindeman1998`** — Bindeman, Davis & Drake (1998) Ion microprobe study of
plagioclase-basalt partition experiments at natural concentration levels of
trace elements. *GCA* 62:1175-1193.
Experimental plagioclase-melt partition coefficients for trace elements as a
function of anorthite content.
*Used for:* nothing in the code yet. See section 3.7.

**`kress_carmichael1991`** — Kress & Carmichael (1991) The compressibility of
silicate liquids containing Fe2O3 and the effect of composition, temperature,
oxygen fugacity and pressure on their redox states. *CMP* 108:82-92.
The standard model relating melt Fe3+/Fe2+ to fO2, T, P and composition.
*Used for:* nothing in the code yet. See section 3.7.

---

### 3.5 Diffusion-coefficient sources, by mineral

#### Olivine

**`dohmen_chakraborty2007`** — Dohmen & Chakraborty (2007) Fe-Mg diffusion in
olivine II: point defect chemistry, change of diffusion mechanisms and a model
for calculation of diffusion coefficients in natural olivine. *Phys. Chem.
Minerals* 34:409-430, with erratum 34:597-598.
Part II of the pair that established the modern olivine Fe-Mg law. It
separates a transition-metal extrinsic (TaMED) regime, in which D depends on
fO2 through an exponent of 1/6, from a pure extrinsic (PED) regime at reducing
conditions in which it does not, and gives a model equation in T, P, fO2 and
X_Fe. The erratum (published 31 August 2007, free to read at
doi:10.1007/s00269-007-0185-3) corrects the composition term of eqs 27 and 28
from `3 X_Fe` to `3 (X_Fe - 0.1)`, and of the global eq. 29 from `3 X_Fe` to
`3 (X_Fe - 0.14)`; nothing else changes.
*Used for:* both olivine Fe-Mg entries —
`log D = -9.21 - (201000 + (P-1e5) 7e-6)/(2.303 R T) + (1/6) log(fO2/1e-7) +
3 (X_Fe - 0.1)` for TaMED and the -8.91 / 220 kJ/mol form for PED, D//[001],
fO2 and P in Pa; the ~6x anisotropy of [001] over [100] and [010]; and the
TaMED/PED mechanism change itself.
**Both entries are verified** (29 September 2026) against the printed pages
and the erratum. From 20 to 29 September 2026 the code used the misprinted
`3 X_Fe`, which makes D twice as large and olivine times half as long.

**`dohmen2007`** — Dohmen, Becker & Chakraborty (2007) Fe-Mg diffusion in
olivine I: experimental determination between 700 and 1200 C as a function of
composition, crystal orientation and oxygen fugacity. *Phys. Chem. Minerals*
34:389-407.
The experimental half of the pair.
*Used for:* a secondary citation on the olivine TaMED entry.

**`chakraborty1997`** — Chakraborty (1997) Rates and mechanisms of Fe-Mg
interdiffusion in olivine at 980-1300 C. *JGR* 102:12317-12331.
The earlier calibration, a single Arrhenius law at fixed fO2 around Fo86.
*Used for:* the `ol_FeMg_chakraborty1997` comparison entry,
D0 = 1.0e-9 m2/s and Q = 226 kJ/mol. Unverified, and kept only so older
published timescales can be reproduced.

**`holzapfel2007`** — Holzapfel, Chakraborty, Rubie & Frost (2007) Effect of
pressure on Fe-Mg, Ni and Mn diffusion in olivine. *PEPI* 162:186-198.
The source of the olivine activation volume.
*Used for:* dV = 7e-6 m3/mol in both olivine entries (reached through
Costa et al. 2008 p. 571); a secondary citation on the TaMED entry.

**`petry2004`** — Petry, Chakraborty & Palme (2004) Experimental determination
of Ni diffusion coefficients in olivine. *GCA* 68:4179-4188.
Ni diffusion in olivine as a function of T, composition, fO2 and orientation.
*Used for:* nothing yet; olivine declares a Ni species with no coefficient
behind it. See section 3.7.

**`coogan2005ca`** — Coogan, Hain, Stahl & Chakraborty (2005) Experimental
determination of the diffusion coefficient for calcium in olivine between
900 C and 1500 C. *GCA* 69:3683-3694.
*Used for:* nothing yet; olivine declares a Ca species with no coefficient
behind it. See section 3.7.

#### Orthopyroxene

**`dias2025`** — Dias, Dohmen & Behrens (2025) Fe-Mg interdiffusion in
orthopyroxene: complex interdependencies of temperature, composition and
oxygen fugacity. *GCA* 395:195-211.
The current recommended opx law. It refits the Dohmen et al. (2016)
experiments together with new Fe-rich runs and finds two regimes separated at
log fO2 = -10 Pa, with a composition exponent that depends on temperature.
*Used for:* the `opx_FeMg_dias2025` entry, eqs 22-25:
above the switch, `D = 3.085e-8 (fO2[Pa]/1e-7)^0.25 exp(-284 kJ/mol / RT)
10^(m1 (X_Fe - 0.1))` with `m1 = -2.37 (1e4/T) + 21.09`; at or below it,
`D = 1.93e-10 exp(-246 kJ/mol / RT) 10^(m2 (X_Fe - 0.1))` with
`m2 = 2.96 (1e4/T) - 21.08`. Also the `superseded_by` target for both
Dohmen et al. (2016) and Dias & Dohmen (2024), and the statement that the
authors advise against extrapolating below 900 C.

**`dias_dohmen2024`** — Dias & Dohmen (2024) Experimental determination of
Fe-Mg interdiffusion in orthopyroxene as a function of Fe content. *CMP*
179:36.
The first calibration of the temperature-dependent Fe effect, at a single
oxygen fugacity.
*Used for:* the `opx_FeMg_dias_dohmen2024` entry, eqs 12-13:
`D = 3.8e-9 exp(-261.07 kJ/mol / RT) 10^(m (X_Fe - 0.09))` with
`m = -2.711e4/T + 23.5408`, valid at log fO2 = -7 Pa only; the fitted m values
of its Table 1 (3.7, 3.0, 2.4, 1.1 at 1102, 1050, 1000 and 950 C), which are
what disambiguated the typeset m(T) expression and are checked by the tests;
and the finding that the Ganguly & Tazzoli (1994) composition dependence is
too strong below about 1000 C, which lengthened the Mount St Helens timescales
of Saunders et al. (2012) from about 5 to 11 weeks.

**`dias2025ree`** — Dias, Dohmen & Hartmann (2025) Diffusion of Eu, Ce and Lu
in orthopyroxene. *GCA* 410:85-100.
TOF-SIMS measurements on thin-film diffusion couples along [001] in natural
orthopyroxene, giving REE diffusivities faster and with lower activation
energies than earlier opx and diopside data.
*Used for:* the three REE entries — Lu (eq. 6,
D0 = 1.51e-9 m2/s, Q = 263 +/- 52 kJ/mol, with an fO2 exponent of 1/7), Ce
(eq. 7, 5.75e-14, 166 +/- 40) and Eu (eq. 8, 1.90e-14, 147 +/- 22, assumed
trivalent). **The Lu entry is unverified** because the reference fugacity in
its fO2 term is not defined in the paper; Diffusor assumes the IW buffer and
says so.

**`dohmen2016`** — Dohmen, ter Heege, Becker & Chakraborty (2016) Fe-Mg
interdiffusion in orthopyroxene. *Am. Mineral.* 101:2210-2221.
The first modern experimental opx Fe-Mg calibration, with an explicit fO2
dependence far weaker than the olivine analogy had assumed, and a measured
anisotropy.
*Used for:* two entries — the Fs9 fO2-dependent law (eq. 1,
`log D = -5.95 - 308 kJ/mol/(ln10 R T) + 0.053 log fO2[Pa]`, with the
composition correction `10^(m (X_Fe - 0.09))`, m = 1, from p. 2219) and the
Fs1 fO2-independent fit (log D0 = -3.78, Q = 377 kJ/mol); the anisotropy
`D_a = D_c/3.5` used by **every** opx entry including the Dias ones; the
Allan et al. (2013) form of the Ganguly & Tazzoli law as their eq. 3; the
critique of the 1/6 fO2 exponent; the test case at 950 C and log fO2 = -7 Pa
(run OPXD_14, published log D = -19.49 +/- 0.07). Marked superseded by
Dias et al. (2025).

**`ganguly_tazzoli1994`** — Ganguly & Tazzoli (1994) Fe2+-Mg interdiffusion in
orthopyroxene: retrieval from the data on intracrystalline exchange reaction.
*Am. Mineral.* 79:930-937.
Derived an interdiffusion coefficient indirectly from Fe-Mg order-disorder
kinetics rather than from any direct diffusion measurement. It was the
standard opx law for two decades, including in Saunders et al. (2012).
*Used for:* two entries, `log D [m2/s] = -9.54 + 2.6 X_Fe - 12530/T` with and
without the `(1/6) log(fO2/fO2_IW)` term. The version without it is what
Ostorero et al. (2022) used and is the entry that reproduces their Kizimen
timescale.

**`schwandt1998`** — Schwandt, Cygan & Westrich (1998) Magnesium
self-diffusion in orthoenstatite. *CMP* 130:390-396.
25Mg SIMS depth profiling on En90Fs10 at the IW buffer, along three axes.
*Used for:* the three `opx_Mg_schwandt1998_*` entries from its Table 3:
(100) log D0 = -3.96 +/- 2.48, Ea = 360 +/- 52; (010) -5.16 +/- 3.69,
339 +/- 77; (001) -8.36 +/- 3.27, 265 +/- 66. The entries note that the three
axes agree within their large uncertainties, so the apparent anisotropy is not
resolved.

**`allan2013`** — Allan, Morgan, Wilson & Millet (2013) From mush to eruption
in centuries: assembly of the super-sized Oruanui magma body. *CMP*
166:143-164.
*Used for:* named in `opx.py` as the origin of the fO2-dependent form of the
Ganguly & Tazzoli law that Dohmen et al. (2016) quote as their eq. 3. Prose
only.

#### Clinopyroxene

**`muller2013`** — Mueller, Dohmen, Becker, ter Heege & Chakraborty (2013)
Fe-Mg interdiffusion rates in clinopyroxene. *CMP* 166:1563-1576.
Thin-film diffusion couples on Di93Hd7 measured by RBS along [001], with no
fO2 dependence resolved between 1e-17 and 1e-11 bar, which the authors relate
to the high Al content of their crystals.
*Used for:* the recommended cpx entry,
`D = 2.77(+/-4.27)e-7 exp(-320.7 +/- 16.0 kJ/mol / RT)` m2/s; `sigma_logD =
0.5` because the fit reproduces the measurements within about 1 log unit; and
the statement that Ca-Mg is slower than Fe-Mg in cpx.
**The paper is internally inconsistent about D0** (see section 4.3).

**`dimanov_sautter2000`** — Dimanov & Sautter (2000) "Average" interdiffusion
of (Fe,Mn)-Mg in natural diopside. *Eur. J. Mineral.* 12:749-760.
The older diopside calibration, and the coefficient behind the published NIDIS
timescales.
*Used for:* the `cpx_FeMg_dimanov_sautter2000` entry, D0 = 9.5e-5 m2/s and
dH = 406 kJ/mol. **Unverified**: the numbers come from the Table 2 footnote of
Petrone et al. (2016), not the primary paper, and the NIDIS script itself uses
9.55e-5.

**`dimanov_wiedenbeck2006`** — Dimanov & Wiedenbeck (2006) (Fe,Mn)-Mg
interdiffusion in natural diopside: effect of pO2. *Eur. J. Mineral.*
18:705-718.
Found an oxygen-fugacity dependence of Fe-Mg interdiffusion in diopside.
*Used for:* the contrast drawn in the Mueller et al. (2013) notes, which
resolved no such dependence. Prose only.

**`brady_mccallister1983`** — Brady & McCallister (1983) Diffusion data for
clinopyroxenes from homogenization and self-diffusion experiments.
*Am. Mineral.* 68:95-105.
Ca-Mg interdiffusion from the homogenisation of (001) pigeonite lamellae in
sub-calcic diopside at 25 kbar, with a stated uncertainty of a factor of 2.
*Used for:* the `cpx_CaMg_brady1983` entry,
`D = 3.89e-7 exp(-360.87 kJ/mol / RT)` m2/s, `sigma_logD = 0.30`.
**Unverified**: from a secondary summary.

**`vanorman2001`**, **`sneeringer1984`**, **`coogan2005li`** — REE diffusion
in diopside with an elastic model (*CMP* 141:687-703); Sr and Sm diffusion in
diopside (*GCA* 48:1589-1608); Li geospeedometry of cooling oceanic crust
(*EPSL* 240:415-424).
Context for cpx trace-element diffusion. None is implemented; see section 3.7.

#### Plagioclase

**`vanorman2014`** — Van Orman, Cherniak & Kita (2014) Magnesium diffusion in
plagioclase: dependence on composition, and implications for thermal resetting
of the 26Al-26Mg early solar system chronometer. *EPSL* 385:79-88.
A joint fit across An23 to An93 plus the An95 data of LaTourrette & Wasserburg
(1998), finding little anisotropy between the b and c directions and
recommending that plagioclase be treated as isotropic for Mg.
*Used for:* the recommended Mg entry, their eq. 4,
`ln D = -6.06(+/-1.10) - 7.96(+/-0.42) x_An - 287(+/-10) kJ/mol/(R T)` with
2-sigma uncertainties; the re-written form of the Costa et al. (2003)
expression used for the `plag_Mg_costa2003` entry; and the finding that
Costa et al. (2003) over-predicts D at low T and low An, so timescales built
on it, including Druitt et al. (2012) for Santorini, are too short by factors
of 2.25 to 5.5.

**`faak2013`** — Faak, Chakraborty & Coogan (2013) Mg in plagioclase:
experimental calibration of a new geothermometer and diffusion coefficients.
*GCA* 123:195-217.
Measured Mg diffusion with silica activity controlled, and found D faster at
higher aSiO2.
*Used for:* a secondary citation on the Audetat et al. (2026) entry, whose
eq. 1 is a joint fit of the Faak and Van Orman data; the silica-activity
effect it established.

**`audetat2026`** — Audétat, Grocolas & Mutch (2026) Ti-in-quartz and
Sr-Ba-Mg-in-feldspars diffusion chronometry: a review of available diffusion
data, and a critical evaluation of applications to natural samples.
*J. Petrology* egag078 (accepted manuscript).
A critical review of the feldspar and quartz diffusion datasets, which
recommends particular calibrations and warns where zoning that looks diffusive
is not.
*Used for:* the recommended silica-activity-dependent Mg entry, their eq. 1,
`log10 D_Mg = -2.99(+/-0.35) X_An - 4.03(+/-0.63) - 262914(+/-15529)/(2.303 R
T) - 1.87(+/-0.45)(1 - aSiO2)`, evaluated at aSiO2 = 1 by default; the
2.8-log-unit Sr gap at An36 and 750 C that the tests check; the confirmation
that Ba is about 1.7 log units slower than Sr in sanidine at 800 C; the
warning that Ba partitions strongly into sanidine so partial dissolution can
mimic diffusion (their section 3.3); and the report of the unpublished
Toussaint et al. (2025) Or98 result. **Taken from the accepted manuscript**;
the confidence level of its uncertainties is not stated and Diffusor treats
them as 1 sigma.

**`grocolas2025`** — Grocolas, Bloch, Bouvier & Müntener (2025) Diffusion of
Sr and Ba in plagioclase: composition and silica activity dependencies, and
application to volcanic rocks. *EPSL* 651:119141. Open access.
Sr and Ba in-diffusion into oriented oligoclase and labradorite at 900-1200 C
with silica activity buffered. Sr comes out 1.5 to 2 orders of magnitude
slower than the 1990s calibrations, which the authors attribute to
SrO-plagioclase reaction fronts in the older source materials; Ba agrees with
Cherniak (2002) to within about half a log unit. No resolvable dependence on
aSiO2 or crystal orientation.
*Used for:* four registry entries — their eq. 7 (Sr) and eq. 8 (Ba) as the
recommended laws, and their Monte Carlo re-fits of the older data as eq. 12
(Giletti & Casserly), eq. 13 (Cherniak & Watson) and eq. 14 (Cherniak 2002);
the covariance assumption in `_compensated_covariance`, taken from their
section 4.3 statement that log10 D0 and Ea follow a linear trend with no
uncertainty envelope; the 7.5 um line-scan resolution preset; and the
`superseded_by` flag on both 1990s Sr entries.

**`giletti_casserly1994`** — Giletti & Casserly (1994) Strontium diffusion
kinetics in plagioclase feldspars. *GCA* 58:3785-3793.
The classic Sr calibration, with a strong anorthite dependence, used for
decades of plagioclase timescales.
*Used for:* the `plag_Sr_giletti_casserly1994` entry,
`D = 8.3176e-5 exp(-276000/RT) 10^(-4.1 X_An)` m2/s, transcribed from the
DMG 2025 course script and cross-checked against the Grocolas eq. 12 re-fit.
Marked **superseded** by Grocolas et al. (2025).

**`cherniak_watson1994`** — Cherniak & Watson (1994) A study of strontium
diffusion in plagioclase using Rutherford backscattering spectroscopy. *GCA*
58:5179-5190.
RBS measurements of Sr at 1 atm, reporting diffusion parallel to b about 0.7
log units slower than parallel to c in some compositions.
*Used for:* the `plag_Sr_cherniak_watson1994` entry **in the re-fitted form of
Grocolas et al. (2025) eq. 13**, which treats the data as isotropic. Marked
superseded.

**`cherniak2002`** — Cherniak (2002) Ba diffusion in feldspar. *GCA*
66:1641-1650.
RBS measurements of Ba in both plagioclase and Or61 sanidine.
*Used for:* two entries in two minerals — the plagioclase entry in the
Grocolas et al. (2025) eq. 14 re-fitted form, and the K-feldspar entry
`D = 2.9e-1 exp(-455 +/- 20 kJ/mol / RT)` m2/s taken from its abstract, which
is the same Or61 crystal as the Cherniak (1996) Sr law and is the coefficient
behind the `sanidine_ba` example.

**`pohl2024`** — Pohl, Behrens, Oeser, Marxer & Dohmen (2024) Li diffusion in
plagioclase crystals and glasses. *Eur. J. Mineral.* 36:985-1003. Open access.
Fitted a multispecies model in which interstitial Li and A1-site Li exchange,
on An61 between 606 and 1114 C, and found chemical Li diffusion charge
balanced by Na to be 1.5 to 2 orders of magnitude slower than the tracer
diffusion of Giletti & Shanahan (1997).
*Used for:* two entries, eq. 21 (interstitial, log D0 = -3.76 +/- 0.58,
Q = 180 +/- 12 kJ/mol, recommended) and eq. 22 (vacancy, -5.53 +/- 0.16,
151.7 +/- 3.2). **Diffusor applies each mechanism as one effective
coefficient, which the multispecies model is not**; the entries say so.

**`giletti_shanahan1997`** — Giletti & Shanahan (1997) Alkali diffusion in
plagioclase feldspar. *Chem. Geol.* 139:3-20.
Tracer diffusion of the alkalis, including Li.
*Used for:* a secondary citation on both Pohl entries, as the faster tracer
comparison.

**`grove1984`** — Grove, Baker & Kinzler (1984) Coupled CaAl-NaSi diffusion in
plagioclase feldspar. *GCA* 48:2113-2121.
The coupled NaSi-CaAl exchange, which is what would relax an anorthite
profile.
*Used for:* the `plag_NaSiCaAl_grove1984` entry,
`D = 1.1e-4 exp(-520 kJ/mol / RT)` m2/s, included **only** so the "is X_An
frozen?" assumption can be checked quantitatively; it is unverified and an
order-of-magnitude placeholder. Also the justification for treating X_An as
frozen. Added to the citation list when the activity term is on.

**`liu_yund1992`** — Liu & Yund (1992) NaSi-CaAl interdiffusion in
plagioclase. *Am. Mineral.* 77:275-283.
*Used for:* named in the Grove entry's notes as the other source to check
before reporting any NaSi-CaAl timescale. Prose only.

**`latourrette_wasserburg1998`** — LaTourrette & Wasserburg (1998) Mg
diffusion in anorthite. *EPSL* 158:91-108.
*Used for:* a secondary citation on the `plag_Mg_costa2003` entry; its An95
data are part of the Van Orman et al. (2014) fit.

#### K-feldspar

**`cherniak1996`** — Cherniak (1996) Strontium diffusion in sanidine and
albite, and general comments on strontium diffusion in alkali feldspars. *GCA*
60:5037-5043.
RBS measurements normal to (001) on Or61 sanidine between 725 and 1075 C, with
a small fO2 effect between air and FMQ that the published law does not
include.
*Used for:* the recommended K-feldspar Sr entry,
`D = 8.4 exp(-450 +/- 13 kJ/mol / RT)` m2/s, with `sigma_logD = 0.03` taken
from Chamberlain et al. (2014).

**`cherniak_watson2020`** — Cherniak & Watson (2020) Ti diffusion in feldspar.
*Am. Mineral.* 105:1040-1051. Open access.
Ti diffusion normal to (001) between 800 and 1000 C, with NNO-buffered
experiments matching those in air, little effect of water and little
anisotropy.
*Used for:* the K-feldspar Ti entry,
`D = 3.01e-6 exp(-342 +/- 47 kJ/mol / RT)` m2/s. The paper publishes no D0
uncertainty and no scatter about the fit, which is why that entry's Monte
Carlo behaviour is limited (section 4.3).

#### Magnetite and titanomagnetite

**`aggarwal_dieckmann2002`** — Aggarwal & Dieckmann (2002) Point defects and
cation tracer diffusion in (Ti_x Fe_(1-x))_(3-d) O4. II. Cation tracer
diffusion. *Phys. Chem. Minerals* 29:707-718.
The tracer-diffusion measurements that underlie the modern magnetite
parameterisation, covering Fe, Co, Mn and Ti as a function of temperature,
oxygen fugacity and Ti content, and resolving the competing vacancy and
interstitial mechanisms.
*Used for:* the primary data behind every `mt_*_vanorman_crispin2010` and
`mt_*_aggarwal2002_*` entry; a secondary citation on the temperature-scaled
Sievwright entries. Diffusor reads the numbers from Van Orman & Crispin
(2010), not from this paper.

**`dieckmann1987`** — Dieckmann, Mason, Hodge & Schmalzried, *Ber.
Bunsenges. Phys. Chem.* 82:778-783. Defects and cation diffusion in magnetite
III: tracer diffusion of foreign tracer cations as a function of temperature
and oxygen potential.
The source of the Cr and Al entries in Van Orman & Crispin (2010) Table 12.
*Used for:* named in `magnetite.py` as the origin of the Cr and Al Table 12
rows. Prose only. **Note the key says 1987 while the record's year field says
1978** (section 4.3).

**`freer_hauptman1978`** — Freer & Hauptman (1978) An experimental study of
magnetite-titanomagnetite interdiffusion. *PEPI* 16:223-231.
Interdiffusion between synthetic Fe3O4 and Fe2.8Ti0.2O4 in sealed silica
tubes, so self-buffered and with poorly constrained fO2. It is the coefficient
behind the Fe-Ti oxide timescales of Costa et al. (2008, Fig. 8) and
Saunders et al. (2012).
*Used for:* the `mt_FeTi_freer_hauptman1978` entry,
`ln D = -15.17 + 13.3 x_Ti - 25870/T` m2/s, transcribed from Van Orman &
Crispin (2010) Table 11. Demoted in favour of the Table 12 tracer entries,
with a long `superseded_note` explaining why no direct interdiffusion
replacement exists.

**`aragon1984`** — Aragon, McCallister & Harrison (1984) Cation diffusion in
titanomagnetites. *CMP* 85:174-185.
Fe-Ti interdiffusion calibrated at the QFM buffer with solid-state buffering,
so better redox control than Freer & Hauptman, but about an order of magnitude
different when extrapolated to x_Ti = 0.15.
*Used for:* the `mt_FeTi_aragon1984` entry,
`ln D = -22.71 + 15.09 x_Ti - 19630/T` m2/s, also from Van Orman & Crispin
Table 11.

**`sievwright2020`** — Sievwright, O'Neill, Tolley, Wilkinson & Berry (2020)
Diffusion and partition coefficients of minor and trace elements in magnetite
as a function of oxygen fugacity at 1150 C. *CMP* 175:40.
Modern LA-ICP-MS measurements on natural magnetite equilibrated with a
silicate melt at 1 bar over FMQ-1 to FMQ+4.89, for 21 elements, fitted to the
same vacancy-plus-interstitial form. All runs are at 1150 C, so the dataset
carries no activation energy.
*Used for:* six registry entries (Ti, Mn, Co, Cr, Al, Mg) from their eq. 5 and
Table 5, `D = D_V1 fO2^(2/3) + D_I1 fO2^(-2/3)` with fO2 in bar; the full
21-element Table 5 is transcribed in `SIEVWRIGHT_TABLE5` and checked against
the published minima by the tests; and as the modern comparison in the
Freer & Hauptman `superseded_note`. **The temperature scaling applied to five
of those entries is Diffusor's own construction**, using Van Orman & Crispin
Table 12 activation energies.

---

### 3.6 Application, validation and data-source papers

**`tomiya2013`** — Tomiya, Miyagi, Saito & Geshi (2013) Short time scales of
magma-mixing processes prior to the 2011 eruption of Shinmoedake volcano,
Kirishima volcanic group, Japan. *Bull. Volcanol.* 75:750.
Titanomagnetite diffusion chronometry of the 2011 Shinmoedake sub-Plinian
eruption, giving mixing-to-eruption timescales of days. It uses the Van Orman
& Crispin (2010) Table 12 formulation at 950 C, log fO2 = -11 and X_Usp = 0.3,
and takes D_Al = D_Ti and D_Mg = D_Fe.
*Used for:* the numerical check the magnetite implementation is validated
against (4.3e-16 m2/s for Ti at 950 C, 6.9e-16 at 900 C, reproduced to about
1 per cent); the log-linear x_Ti interpolation, which is how those numbers are
reachable; the conditions of the `magnetite_shinmoedake` example; and the
documented Fe-branch disagreement (section 4.3). **This is the most directly
relevant paper in the registry to the Shinmoedake work the app is being built
for.**

**`saunders2012`** — Saunders, Blundy, Dohmen & Cashman (2012) Linking
petrology and seismology at an active volcano. *Science* 336:1023-1027.
The Mount St Helens study that tied orthopyroxene and Fe-Ti oxide diffusion
timescales to the seismic record — the template for calibrating diffusion
times against earthquake data.
*Used for:* a secondary citation on the Freer & Hauptman Fe-Ti entry; the
example in the Ganguly & Tazzoli notes of a published timescale that changed
(from about 5 to 11 weeks) when Dias & Dohmen (2024) corrected the composition
dependence.

**`ostorero2022`** — Ostorero, Balcone-Boissard, Boudon et al. (2022)
Correlated petrology and seismicity indicate rapid magma accumulation prior to
eruption of Kizimen volcano, Kamchatka. *Commun. Earth Environ.* 3:290.
Orthopyroxene Fe-Mg diffusion chronometry of the 2010-2013 Kizimen eruption,
modelled at 850 C with the Ganguly & Tazzoli (1994) law and no fO2 term, and
compared with the seismicity.
*Used for:* the only published-timescale validation in the test suite
(2.32 yr +7.16/-1.75 against Diffusor's ~3 yr); their eq. 1, which is the
`opx_FeMg_ganguly_tazzoli1994_nofo2` entry; the measured traverse shipped as
the `opx_kizimen` example, extracted from Supplementary Data 2; the conditions
from Supplementary Data 3 (850 +/- 57 C, NNO +1.28 +/- 0.35 from 21
magnetite-ilmenite pairs); and the plateau-that-continues boundary
description.

**`druitt2012`** — Druitt, Costa, Deloule, Dungan & Scaillet (2012) Decadal to
monthly timescales of magma transfer and reservoir growth at a caldera
volcano. *Nature* 482:77-80.
Plagioclase diffusion chronometry of the Minoan eruption of Santorini,
reconstructing the initial Mg profile from the Sr-An correlation and a
two-melt history.
*Used for:* the measured `plag_santorini` example, crystal S82-30A 12 from
Supplementary Table 1; the 900 C and +/-25 C the example uses; the 10 um laser
and 10-15 um ion-beam resolution presets; and, explicitly, as an example
Diffusor **cannot** reproduce, because its built-in initial conditions cannot
express that reconstruction.

**`chamberlain2014`** — Chamberlain, Morgan & Wilson (2014) Timescales of
mixing and mobilisation in the Bishop Tuff magma body. *CMP* 168:1034.
Feldspar Sr and Ba diffusion chronometry of the Bishop Tuff, measured with a
5 um defocused microprobe beam at 753 to 815 C.
*Used for:* the `sigma_logD` of both sanidine entries (0.03 for Sr, 0.12 for
Ba), which are their derived uncertainties rather than the diffusion papers';
the confirmation that Cherniak's D0 values are as transcribed; the 5 um
defocused-beam resolution preset; the temperature range of the `sanidine_ba`
example; the plateau-that-continues boundary description; and the argument
that paired Sr and Ba profiles test whether a boundary is diffusive.

**`hartley2016`** — Hartley, Morgan, Maclennan, Edmonds & Thordarson (2016)
Tracking timescales of short-term precursors to large basaltic fissure
eruptions through Fe-Mg diffusion in olivine. *EPSL* 439:58-70.
Olivine Fe-Mg chronometry of the Laki 1783-84 fissure eruption.
*Used for:* independent confirmation of the ~6x olivine anisotropy in natural
crystals (p. 60), which is part of the unverified olivine entry's
`verified_from`; and the conditions of the `olivine_laki` example
(1150 +/- 30 C, FMQ-1 +/- 0.5).

**`sato2022`** — Sato, Ban, Yoshida & Andrews (2022) Magma plumbing system and
eruption processes of the Okama pyroclastics, Zao volcano, revealed by
orthopyroxene Fe-Mg diffusion chronometry. *JVGR* 429:107607.
An application of the Dohmen et al. (2016) opx law to a Japanese arc volcano.
*Used for:* a secondary citation on the `opx_FeMg_dohmen2016` entry.

**`polo_sanchez2023`** — Polo-Sanchez, Druitt, Cluzel & Devidal (2023)
Pyroxene diffusion chronometry of the magmatic plumbing system of Santorini
volcano. *Front. Earth Sci.* 11:1149446.
*Used for:* a secondary citation on the `opx_FeMg_dohmen2016` entry.

**`zellmer1999`** — Zellmer, Blake, Vance, Hawkesworth & Turner (1999)
Plagioclase residence times at two island arc volcanoes determined by Sr
diffusion systematics. *CMP* 136:345-357.
An early application of Sr-in-plagioclase chronometry, and one of the sources
for the inverse correlation between equilibrium Sr and anorthite content.
*Used for:* named in the Giletti & Casserly entry's notes alongside
Dohmen et al. (2017) and Costa et al. (2003) for the sign of the activity
term. Prose only.

**`grocolas2025cmp`** — Grocolas, Müntener, Bloch, Escrig, Ulyanov & Bouvier
(2025) Cooling rates and melt extraction timescales determined by diffusion
chronometry on shallow crustal plutonic rocks. *CMP* 180:45. Open access.
Applies the new Sr and Ba diffusivities to the Adamello batholith.
*Used for:* nothing in the code; an application of the recommended laws. See
section 3.7.

---

### 3.7 Registry entries not wired into the code

Nineteen keys resolve but are never named by key anywhere outside
`references.py`, so they cannot reach an exported methods block unless passed
through `collect_citations(extra=...)` by hand. Eight of them are cited in
user-visible prose, where the citation is a string rather than a key.

| key | cited in prose? | why it is in the registry |
| --- | --- | --- |
| `allan2013` | yes, `opx.py` | origin of the fO2-dependent Ganguly & Tazzoli form |
| `bindeman1998` | no | plagioclase-melt partitioning, for future Mg work |
| `cherniak_dimanov2010` | yes, `minerals/definitions.py` | pyroxene diffusion compilation |
| `coogan2005ca` | no | Ca in olivine, a declared species with no law |
| `coogan2005li` | no | Li geospeedometry in cpx |
| `costa2020` | yes, `fitting/fit.py` | justifies the default time search window |
| `dieckmann1987` | yes, `magnetite.py` | primary source of the Cr and Al Table 12 rows |
| `dimanov_wiedenbeck2006` | yes, `cpx.py` | the fO2 dependence Mueller et al. did not resolve |
| `girona_costa2013` | no | DIPRA, comparison software |
| `grocolas2025cmp` | no | application of the recommended Sr and Ba laws |
| `huebner1971` | yes, `magnetite.py` | buffer equations behind Tables 10 and 11 |
| `kress_carmichael1991` | no | melt Fe3+/Fe2+ redox model |
| `liu_yund1992` | yes, `plagioclase.py` | second check on NaSi-CaAl |
| `morgan2004` | yes, `dataio/greyscale.py` | origin of the grey-value approach |
| `mutch2021` | yes, four modules | DFENS, the covariance-sampling argument |
| `petry2004` | no | Ni in olivine, a declared species with no law |
| `sneeringer1984` | no | Sr and Sm in diopside |
| `vanorman2001` | yes, `opx.py` | REE in diopside, the slower comparison |
| `zellmer1999` | yes, `plagioclase.py` | sign of the Sr activity term |

The straightforward fix for the eight prose cases is to add them to the
relevant coefficient's `secondary_citations`, which is the mechanism that
already exists for exactly this.

---

## 4. Cross-checks

### 4.1 Equation index

Every equation, formula and convention in the code, with its module and its
source. Rows marked **Diffusor** have no literature source: they are the
package's own constructions, and they are the first places to look when a
result disagrees with a published one.

#### Solvers, fitting and corrections

| what | where | source |
| --- | --- | --- |
| `dC/dt = (1/x^m) d/dx [x^m D dC/dx]` | `solvers/numerical.py` | Crank (1975) eq. 1.7 (plane), 5.4 (cylinder), 6.3 (sphere) |
| step in an infinite medium | `analytical.step_infinite` | Crank eq. 2.14; Costa et al. (2008) eq. 9 |
| semi-infinite, fixed surface | `analytical.semi_infinite_fixed_surface` | Crank eq. 2.45 |
| band of half-width h | `analytical.band_infinite` | Crank eq. 2.15 |
| plane sheet | `analytical.plane_sheet` | Crank eq. 4.17 |
| cylinder | `analytical.cylinder` | Crank eq. 5.22 |
| sphere | `analytical.sphere` | Crank eq. 6.18 |
| fractional uptake, plane / cylinder / sphere | `analytical.fraction_*` | Crank eqs 4.18, 5.23, 6.20 |
| explicit scheme | `numerical.solve_1d` (theta = 0) | Crank eq. 8.31 |
| stability limit `dt <= dx^2/(2D)` | `geometry.stability_dt`, `numerical.solve_1d` | Crank eq. 8.33; Courant 0.2-0.4 from DMG (2025) |
| theta-scheme / Crank-Nicolson | `numerical.solve_1d` | Crank section 8.4 eq. 8.35; Dohmen et al. (2017) App. eq. A21 |
| half-node D, conservative form | `numerical._operator` | Dohmen et al. (2017) App. eqs A17-A19 |
| activity flux term `-theta D C dX_An/dx` | `numerical._operator` | Costa et al. (2003) eq. 7; Dohmen et al. (2017) App. eqs A7-A8, coefficients A20 |
| `theta = A_i/(R T)` | `plagioclase.activity_theta` | Dohmen et al. (2017) App. eq. A6 |
| equilibrium profile `C0 exp(A X_An/RT)` | `plagioclase.equilibrium_profile` | Dohmen et al. (2017) App. eqs A13, A14 |
| mirror ghost node at a zero-flux end | `numerical._apply_bc_rows` | Crank section 8.4 |
| `dC/dt = (m+1) D d2C/dx2` at r = 0 | `numerical._apply_bc_rows` | Crank section 8.5 eq. 8.45 (L'Hopital) |
| boundary condition types | `solvers/boundary.py` | Crank section 1.3; Costa et al. (2008) |
| `Dt_eff = integral D(T(t)) dt` | `history.effective_Dt` | Crank section 7.2 eq. 7.7; Lasaga (1983) |
| Gaussian convolution of the model | `convolution.gaussian_convolve` | Ganguly et al. (1988); Bradshaw & Kent (2017); DMG (2025) Practical 5 |
| upper-bound rule `2 sqrt(Dt) < 3 sigma` | `convolution.resolution_warning` | Bradshaw & Kent (2017) |
| `sigma = d/4` (round spot), `sigma = w/sqrt(12)` (slit) | `gui/main_window.RESOLUTION_PRESETS` | **Diffusor** — standard beam-geometry conversions, no citation in code |
| `D_V = D_a cos^2 a + D_b cos^2 b + D_c cos^2 g` | `minerals.base.direction_factor` | Costa & Chakraborty (2004); DMG (2025) Lecture 6 |
| `D = D0 exp(-(Q + (P-P0) dV)/(R T))` | `coefficients.base.arrhenius` | Crank section 11; Costa et al. (2008) eq. 20 |
| weighted residuals `(obs - model)/sigma`, chi2, reduced chi2, RMSE, R2 | `fitting/objective.py` | standard weighted least squares, no citation |
| log scan then bounded `least_squares` in log10 t | `fitting/fit.py` | **Diffusor** |
| search window 100 s to 3.2e14 s | `fitting/fit.py` | **Diffusor**, range justified by Costa et al. (2020) |
| median and 16/84, 2.5/97.5 percentiles | `fitting/montecarlo.py` | log-normality of times, Mutch et al. (2021) |
| buffer re-evaluated at each sampled T | `montecarlo._draw_conditions` | **Diffusor**, the correction to independent sampling |
| one-at-a-time variance decomposition | `montecarlo.contributions` | **Diffusor** |
| ratio uncertainty `sqrt((B sA)^2 + (A sB)^2)/(A+B)^2` | `dataio/profiles.build_profile` | first-order quadrature, no citation |
| grey value to composition by polynomial fit | `dataio/greyscale.py` | Petrone et al. (2016) Methods; Morgan et al. (2004) |

#### Conventions, buffers and constants

| what | where | source |
| --- | --- | --- |
| `X_Fe = Fe/(Fe+Mg)` molar, all Fe as Fe2+ | `thermo/units.py`, `minerals/definitions.py` | Deer et al. (1992), and the usage of Dohmen & Chakraborty (2007), Dohmen et al. (2016), Mueller et al. (2013) |
| `X_Fo`, `Mg#`, `X_An = Ca/(Ca+Na+K)`, `X_Or` | same | Deer et al. (1992) |
| `x_Ti = X_Usp/3` for (Ti_x Fe_(1-x))3 O4 | `thermo/units.py`, `minerals/definitions.py` | Van Orman & Crispin (2010) |
| Stormer-style recalculation **not** applied | `thermo/units.py` | deliberate, stated in the docstring |
| wt% oxide to cation moles | `thermo/units.oxide_wt_to_cation_moles` | IUPAC (2021) masses |
| `log fO2 = A/T + B + C (P_bar - 1)/T` | `thermo/buffers.py` | Frost (1991) Table 1 |
| NNO free energy | `thermo/buffers.py` | O'Neill & Pownceby (1993) |
| FMQ free energy | `thermo/buffers.py` | O'Neill (1987) |
| Frost pressure term added to the O'Neill forms | `thermo/buffers.py` | **Diffusor**, so the two agree on P |
| R, k_B, N_A, eV | `constants.py` | CODATA 2018 |
| atomic and oxide masses | `constants.py` | IUPAC 2021 |
| `SEC_PER_YEAR = 3.15576e7` s | `constants.py` | IAU Julian year |
| 1 atm = 1.01325 bar | `coefficients/magnetite.py` | SI definition |
| 1 bar = 1e5 Pa (`log_fo2_Pa = log_fo2_bar + 5`) | `coefficients/base.py` | SI definition |
| D converted m2/s to um2/s by 1e12 | `fitting/model.py` | **Diffusor**, conditioning |

#### Constructions that are Diffusor's, not a paper's

These five are the ones most likely to explain a disagreement with a published
number, and every one is stated in the relevant entry's notes:

1. **Log-linear interpolation in x_Ti** between the Van Orman & Crispin (2010)
   Table 12 rows at x_Ti = 0 and 0.2, clipped at 1.5. Reproduces Tomiya et al.
   (2013) for Ti to about 1 per cent, which is the evidence that it is what
   they did, but they do not say so.
2. **Temperature scaling of the Sievwright et al. (2020) entries** by the Van
   Orman & Crispin Table 12 activation energies. Absent from the paper.
3. **The Grocolas covariance matrix**, built by `_compensated_covariance` from
   their prose statement that log10 D0 and Ea are perfectly compensated.
4. **The Lu reference fugacity**, assumed to be the IW buffer because the
   paper does not define it.
5. **The beam-sigma conversions** `d/4` and `w/sqrt(12)`.

### 4.2 The coefficient registry, one row per entry

| key | species | source and equation | MC scatter | flags |
| --- | --- | --- | --- | --- |
| `cpx_CaMg_brady1983` | Ca-Mg | Brady & McCallister (1983), (homogenisation experiments, 25 kbar) | 0.3 | **unver** |
| `cpx_FeMg_dimanov_sautter2000` | Fe-Mg | Dimanov & Sautter (2000), as tabulated by Petrone et al. (2016), Table 2 fo... | 0.5 | **unver** |
| `cpx_FeMg_muller2013` | Fe-Mg | Mueller et al. (2013), abstract / Fig. 5b | 0.5 | rec |
| `kfs_Ba_cherniak2002` | Ba | Cherniak (2002), abstract | 0.12 | rec |
| `kfs_Sr_cherniak1996` | Sr | Cherniak (1996), abstract | 0.03 | rec |
| `kfs_Ti_cherniak_watson2020` | Ti | Cherniak & Watson (2020), abstract | none | rec |
| `mt_Al_sievwright2020` | Al | Sievwright et al. (2020), 5 and Table 5 | 0.2 | - |
| `mt_Al_vanorman_crispin2010` | Al | Van Orman & Crispin (2010), Table 12 | 0.3 | - |
| `mt_Co_sievwright2020` | Co | Sievwright et al. (2020), 5 and Table 5 | 0.2 | - |
| `mt_Co_vanorman_crispin2010` | Co | Van Orman & Crispin (2010), Table 12 | 0.3 | - |
| `mt_Cr_sievwright2020` | Cr | Sievwright et al. (2020), 5 and Table 5 | 0.2 | - |
| `mt_Cr_vanorman_crispin2010` | Cr | Van Orman & Crispin (2010), Table 12 | 0.3 | - |
| `mt_Fe_aggarwal2002_MH` | Fe | Van Orman & Crispin (2010), Table 10 | 0.3 | - |
| `mt_Fe_aggarwal2002_MH_xti02` | Fe | Van Orman & Crispin (2010), Table 11 | 0.3 | - |
| `mt_Fe_aggarwal2002_WM` | Fe | Van Orman & Crispin (2010), Table 10 | 0.3 | - |
| `mt_Fe_aggarwal2002_WM_xti02` | Fe | Van Orman & Crispin (2010), Table 11 | 0.3 | - |
| `mt_Fe_vanorman_crispin2010` | Fe | Van Orman & Crispin (2010), Table 12 | 0.3 | rec |
| `mt_FeTi_aragon1984` | Fe-Ti | Aragon et al. (1984), Van Orman & Crispin (2010) Table 11 | 0.5 | - |
| `mt_FeTi_freer_hauptman1978` | Fe-Ti | Freer & Hauptman (1978), Van Orman & Crispin (2010) Table 11 | 0.5 | - |
| `mt_Mg_sievwright2020` | Mg | Sievwright et al. (2020), 5 and Table 5 | 0.2 | - |
| `mt_Mn_sievwright2020` | Mn | Sievwright et al. (2020), 5 and Table 5 | 0.2 | - |
| `mt_Mn_vanorman_crispin2010` | Mn | Van Orman & Crispin (2010), Table 12 | 0.3 | - |
| `mt_Ti_aggarwal2002_MH` | Ti | Van Orman & Crispin (2010), Table 10 | 0.3 | - |
| `mt_Ti_aggarwal2002_MH_xti02` | Ti | Van Orman & Crispin (2010), Table 11 | 0.3 | - |
| `mt_Ti_aggarwal2002_WM` | Ti | Van Orman & Crispin (2010), Table 10 | 0.3 | - |
| `mt_Ti_aggarwal2002_WM_xti02` | Ti | Van Orman & Crispin (2010), Table 11 | 0.3 | - |
| `mt_Ti_sievwright2020` | Ti | Sievwright et al. (2020), 5 and Table 5 | 0.2 | - |
| `mt_Ti_vanorman_crispin2010` | Ti | Van Orman & Crispin (2010), Table 12 | 0.3 | rec |
| `ol_FeMg_chakraborty1997` | Fe-Mg | Chakraborty (1997), (Arrhenius fit at fO2 = 1e-7 Pa, Fo86) | 0.3 | **unver** |
| `ol_FeMg_dohmen_chakraborty2007_ped` | Fe-Mg | Dohmen & Chakraborty (2007), eq. 28 as corrected by the erratum | 0.21 | - |
| `ol_FeMg_dohmen_chakraborty2007_tamed` | Fe-Mg | Dohmen & Chakraborty (2007), eq. 27 as corrected by the erratum | 0.21 | rec |
| `opx_Ce_dias2025` | Ce | Dias et al. (2025), 7 | 0.3 | - |
| `opx_Eu_dias2025` | Eu | Dias et al. (2025), 8 | 0.3 | - |
| `opx_FeMg_dias2025` | Fe-Mg | Dias et al. (2025), 22-25 | 0.2 | rec |
| `opx_FeMg_dias_dohmen2024` | Fe-Mg | Dias & Dohmen (2024), 12-13 | 0.2 | superseded |
| `opx_FeMg_dohmen2016` | Fe-Mg | Dohmen et al. (2016), 1 (+ compositional correction, p. 2219) | 0.1 | superseded |
| `opx_FeMg_dohmen2016_fs1` | Fe-Mg | Dohmen et al. (2016), 1 with n = 0 (Opx8 fit, p. 2215) | 0.1 | - |
| `opx_FeMg_ganguly_tazzoli1994` | Fe-Mg | Ganguly & Tazzoli (1994), Dohmen et al. (2016) eq. 3. Ostorero et al. (2022... | 0.5 | - |
| `opx_FeMg_ganguly_tazzoli1994_nofo2` | Fe-Mg | Ganguly & Tazzoli (1994), Ostorero et al. (2022) eq. 1 | 0.5 | - |
| `opx_Lu_dias2025` | Lu | Dias et al. (2025), 6 | 0.3 | **unver** |
| `opx_Mg_schwandt1998_a` | Mg | Schwandt et al. (1998), Table 3 | none | - |
| `opx_Mg_schwandt1998_b` | Mg | Schwandt et al. (1998), Table 3 | none | - |
| `opx_Mg_schwandt1998_c` | Mg | Schwandt et al. (1998), Table 3 | none | - |
| `plag_Ba_cherniak2002` | Ba | Cherniak (2002), Grocolas et al. (2025) eq. 14 | 0.3 | - |
| `plag_Ba_grocolas2025` | Ba | Grocolas et al. (2025), 8 | cov | rec |
| `plag_Li_pohl2024_interstitial` | Li | Pohl et al. (2024), 21 | 0.3 | rec |
| `plag_Li_pohl2024_vacancy` | Li | Pohl et al. (2024), 22 | 0.3 | - |
| `plag_Mg_audetat2026` | Mg | Audetat et al. (2026), 1 | 0.25 | rec |
| `plag_Mg_costa2003` | Mg | Costa et al. (2003), as re-written by Van Orman et al. (2014), p. 84 | 0.3 | - |
| `plag_Mg_vanorman2014` | Mg | Van Orman et al. (2014), 4 | 0.25 | rec |
| `plag_NaSiCaAl_grove1984` | NaSi-CaAl | Grove et al. (1984), (cooling-rate speedometry calibration) | 0.7 | **unver** |
| `plag_Sr_cherniak_watson1994` | Sr | Cherniak & Watson (1994), Grocolas et al. (2025) eq. 13 | 0.3 | superseded |
| `plag_Sr_giletti_casserly1994` | Sr | Giletti & Casserly (1994), (An-dependent form as implemented in the DMG Shor... | 0.3 | superseded |
| `plag_Sr_grocolas2025` | Sr | Grocolas et al. (2025), 7 | cov | rec |

"MC scatter" is what a Monte Carlo draw uses: `cov` for the two entries with a
covariance matrix, otherwise the 1-sigma scatter of log10 D about the fit, or
`none` where the source publishes neither and only Q can be varied.

### 4.3 Known discrepancies and open questions

The list the document exists for. Each item says what is wrong or unresolved,
what Diffusor currently does, and what would settle it.

#### Verification gaps

1. **Settled: the Dohmen & Chakraborty (2007) olivine entries.** Checked on
   29 September 2026 against pp. 424-425 and the erratum (Phys. Chem.
   Minerals 34:597-598): -9.21, 201 kJ/mol, -8.91, 220 kJ/mol, 7e-6 m3/mol,
   the 1/6 exponent on fO2/1e-7 Pa, the 1e-10 Pa regime boundary, log 6 for
   [100] and [010], and the erratum's `3 (X_Fe - 0.1)`. The printed pages
   alone say `3 X_Fe`, which is why reading them without the erratum (as the
   20 September audit did) gives D twice too large.

2. **Five entries are unverified.** Besides `ol_FeMg_chakraborty1997`:
   `cpx_FeMg_dimanov_sautter2000` (numbers from a table footnote in
   Petrone et al. 2016, whose own NIDIS script uses 9.55e-5 rather than
   9.5e-5 m2/s, a 0.5 per cent difference); `cpx_CaMg_brady1983` (secondary
   summary); `plag_NaSiCaAl_grove1984` (order-of-magnitude placeholder, and
   the notes say not to report a timescale from it); `opx_Lu_dias2025` (see
   item 4).

3. **`sigma_logD` is often Diffusor's judgement.** Where a paper states no
   scatter, the value was chosen from what the paper says qualitatively (for
   example 0.5 for Mueller et al. 2013 from "within about 1 log unit", 0.21
   for olivine from the scatter of the experimental database, 0.2 for
   Dias et al. 2025 from how well the regressions reproduce the experiments).
   These propagate straight into every reported uncertainty. They are
   documented per entry but they are not published numbers.

#### Places where the code interprets a source

4. **The Lu fO2 term has an undefined reference.** Dias, Dohmen & Hartmann
   (2025) eq. 6 writes `(fO2/fO2_0)^(1/7)` without defining fO2_0. Diffusor
   assumes the IW buffer, on the strength of the paper stating the law for
   "fO2 close to the IW buffer", and flags the entry unverified.

5. **The Grocolas covariance is constructed.** The paper publishes no matrix;
   `_compensated_covariance` encodes their stated assumption that log10 D0 and
   Ea are perfectly correlated. If they later publish a real matrix, replace
   it. Note also that they say this assumption *slightly underestimates* the
   uncertainty, so the intervals Diffusor reports for Sr and Ba are, by the
   authors' own account, slightly narrow.

6. **The Sievwright temperature scaling is an addition.** Entries for Ti, Mn,
   Co, Cr and Al are anchored at 1150 C and scaled with Table 12 activation
   energies. At 1150 C they agree with the Table 12 entries within 0.5 log
   units for Ti, Mn and Co, and within 0.4 for Cr above FMQ+2, but **Al
   differs by up to 2 log units at FMQ-1**. Far from 1150 C the answer rests
   entirely on the borrowed energies.

7. **Magnetite Mg is a single-temperature entry.** `mt_Mg_sievwright2020`
   returns the 1150 C value at any temperature and warns. Tomiya et al. (2013)
   instead set D_Mg = D_Fe, which is a different and probably better choice
   for magmatic work.

8. **The x_Ti interpolation is a reconstruction** of what Tomiya et al. (2013)
   must have done, supported by reproducing their Ti numbers to about 1 per
   cent.

9. **Audétat et al. (2026) is an accepted manuscript.** Page numbers are not
   assigned, the confidence level of the +/- values is not stated and Diffusor
   treats them as 1 sigma, and the silica activity is **fixed at aSiO2 = 1**
   and cannot be changed from the interface, only from Python. A
   silica-undersaturated melt would give D up to 1.87 log units lower.

10. **Pohl et al. (2024) Li is approximated.** The paper fits a multispecies
    model with interstitial and A1-site Li in exchange; Diffusor applies each
    mechanism as one effective coefficient.

**Image profiles depart from NIDIS `greyvalues.m`** (`dataio/image_profiles.py`):

* Sampling is **bilinear** by default; MATLAB `improfile` defaults to nearest
  neighbour, which is offered as an option.
* NIDIS's guideline is straight; Diffusor's may be a **polyline**, with each
  line perpendicular to the smoothed local direction. On the concave side of
  a tight curve, neighbouring lines converge and share pixels.
* NIDIS reads `2 hp` values over `2 hp` px (a step slightly over 1 px) and
  `int32(len)` lines; Diffusor reads `2 hp + 1` values including both ends,
  and `round(len / spacing) + 1` lines including both guideline ends.
* The default outlier test is **median +/- 3 MAD, repeated**, not NIDIS's
  single mean +/- 1 SD, which removes about 32 % of perfectly good Gaussian
  values and, without value limits, lets a crack inflate the SD so much that
  the crack itself survives. NIDIS's condition
  `v > m + s || v < m - s && v ~= 0` exempts zeros only on the low side (in
  MATLAB `&&` binds before `||`); Diffusor has no zero exemption and relies on
  the value limits for zero-valued pixels.
* **The standard error understates the uncertainty.** Lines 1 px apart share
  pixels through bilinear sampling and the electron interaction volume, so
  they are not independent, and SE = SD / sqrt(N) with N of several hundred
  is too small. NIDIS exports it and Diffusor uses it by default for
  continuity, but the calibration dialog offers the SD across lines, and the
  Monte Carlo measurement-noise term is the safer place for the real scatter.

#### Disagreements with published numbers

11. **Mueller et al. (2013) contradicts itself about D0.** The abstract and
    the Fig. 5b annotation give 2.77e-7 m2/s; the running text on p. 1570
    gives 2.77e-8. Diffusor uses 2.77e-7, the value that appears twice. **This
    is a factor of ten in every clinopyroxene timescale** and should be
    checked against the authors before publishing.

12. **The magnetite Fe branch disagrees with Tomiya et al. (2013).** Diffusor
    reproduces their Ti values to about 1 per cent, and Fe to 11 per cent at
    900 C, but is a **factor of 2.1 high at 950 C**: they place the Fe
    diffusion minimum near 950 C while this implementation puts it near
    900-920 C, so they evidently treated the Fe composition dependence
    differently. Prefer Ti for Shinmoedake timescales until this is resolved.
    Resolving it matters directly for the intended Shinmoedake work.

13. **Santorini is not reproduced, by design.** Diffusor returns about 90
    years against the published 47, with a reduced chi-squared in the
    hundreds, because Druitt et al. (2012) reconstructed the initial Mg
    profile from the Sr-An correlation and a two-melt history. Implementing
    that class of initial condition is the obvious improvement.

14. **Kizimen is reproduced only approximately**, about 3 yr against 2.32
    (+7.16/-1.75) yr, because the supplement contains only the 2 um microprobe
    traverse while the authors fitted high-resolution BSE grey-value profiles.
    The 4 um zone is crossed by three or four points.

15. **Ostorero et al. (2022) swap the andesite and dacite temperatures**
    between their main text and their Supplementary Data 3. Diffusor follows
    the table (850 +/- 57 C).

#### Extrapolation beyond calibration

16. **Arc fO2 is above the Dias et al. (2025) calibrated range.** Natural arc
    magmas usually sit above log fO2 = -7 Pa, so the fO2^(1/4) term is
    extrapolated; at 950 C and NNO+1 that makes D about 0.4 log units faster
    than Dohmen et al. (2016). The authors also advise against extrapolating
    below 900 C, which the Kizimen example at 850 C does.

17. **Ganguly & Tazzoli (1994) is calibrated for 500-800 C** and the Kizimen
    example runs it at 850 C. Dias & Dohmen (2024) find it overstates the
    composition dependence below about 1000 C.

18. **The two Dias et al. (2025) regimes do not join at log fO2 = -10 Pa.**
    This is intentional — the authors infer a change of mechanism — but a
    profile modelled near the switch will jump.

#### Registry hygiene

19. **`dieckmann1987` carries the year 1978.** The key says 1987, the record's
    `year` field and the volume (Ber. Bunsenges. 82) say 1978, and
    `magnetite.py` cites it in prose as "Dieckmann et al. (1987)". The key
    should be renamed or the prose corrected. `Reference.short()` renders it
    as "Dieckmann et al. (1978)", so `REFERENCES.md` and the key disagree.

20. **Seven sources are named in user-visible prose with no registry key**, so
    they cannot appear in a methods block and are invisible to the test that
    checks citation keys resolve: Toussaint et al. (2025) and Mutch et al.
    (2022) (feldspar notes), Kroll et al. (1997) and Cherniak & Liang (2007)
    (opx notes), Sauerzapf et al. (2008) and Petrone et al. (2018) (dataset
    records), Stormer (1983) (`thermo/units.py`).

21. **Nineteen registry keys are never referenced by key** (section 3.7), so
    eight sources that the prose actually relies on never reach an export.

22. **Three declared species have no coefficient**: olivine Ni, Mn and Ca.
    They appear in the species drop-down and then the coefficient list is
    empty. Petry et al. (2004) and Coogan et al. (2005) are already in the
    registry for Ni and Ca.

23. **Only 6 of the 21 elements in `SIEVWRIGHT_TABLE5` are registered.** The
    other 15 (Ni, Zn, Sc, Ga, In, Y, Lu, V3+, V4+, Zr, Hf, U, Nb, Ta, Mo) are
    transcribed and tested but not reachable from the application.

24. **`3.15576e7` is hard-coded** in `dataio/export.py` for the Monte Carlo
    CSV instead of importing `SEC_PER_YEAR`. Same value today; two places to
    change tomorrow.

### 4.4 Verification status

| mineral | entries | verified | recommended |
| --- | --- | --- | --- |
| olivine | 3 | 2 | 1 |
| orthopyroxene | 12 | 11 | 1 |
| clinopyroxene | 3 | 1 | 1 |
| plagioclase | 11 | 10 | 5 |
| K-feldspar | 3 | 3 | 3 |
| magnetite | 22 | 22 | 2 |
| **total** | **54** | **49** | **13** |

Checks the test suite makes against numbers printed in the sources:

| check | published | Diffusor |
| --- | --- | --- |
| Opx Fe-Mg, 950 C, log fO2 = -7 Pa, //c (Dohmen et al. 2016, run OPXD_14) | log D = -19.49 +/- 0.07 | -19.47 |
| Opx anisotropy D//[001] / D//[100] (Dohmen et al. 2016) | 3.5 | 3.5 |
| Cpx Fe-Mg at 1098 and 1150 C (Petrone et al. 2016, Table 2 footnote) | 3.26e-20 and 1.20e-19 m2/s | 3.25e-20 and 1.19e-19 |
| Titanomagnetite Ti, 950 C, log fO2 = -11, X_Usp = 0.3 (Tomiya et al. 2013) | 4.3e-16 m2/s | 4.35e-16 |
| Titanomagnetite Ti, 900 C, same conditions | 6.9e-16 m2/s | 6.84e-16 |
| Opx Fs9 at 950, 1050, 1100 C (Dias & Dohmen 2024, Table 1) | fitted D0 and m per run | within 0.2 log units, both laws |
| Plagioclase Sr, gap to Giletti & Casserly at An36 and 750 C (Audétat et al. 2026) | 2.8 log units | 2.77 |
| Magnetite, minimum of D against fO2 for 21 elements at 1150 C (Sievwright et al. 2020, Table 5) | log fO2 and log D at the minimum | within 0.1 and 0.05 log units |
| Opx K9_L10C4, Kizimen (Ostorero et al. 2022, Supplementary Data 4) | 2.32 yr (+7.16/-1.75) | about 3 yr |
| Numerical solver against Crank (1975) closed forms | exact | ~3 parts in 1e6, second-order convergence |
| Mass conservation, closed system | exact | 1 part in 1e8 |
| Olivine Fe-Mg TaMED and PED, 1100 C, Fo90, 1e-7 Pa (Dohmen & Chakraborty 2007 eqs 27-28 with the erratum) | log D = -16.856 and -17.279 | -16.856 and -17.279 |

### 4.5 What to improve first

Ordered by how much it would change a published answer, not by effort:

1. Done: the olivine transcription is verified against the primary paper and
   its erratum (item 1).
2. Settle the Mueller et al. (2013) D0 (item 11). A factor of ten sits under
   every clinopyroxene timescale.
3. Resolve the magnetite Fe branch against Tomiya et al. (2013) (item 12).
   This is the Shinmoedake case directly.
4. Implement a reconstructed-initial-profile class so the Santorini kind of
   problem can be expressed (item 13).
5. Expose the silica activity of the Audétat law in the interface (item 9).
6. Add the seven missing citation keys and wire the eight prose-only ones into
   `secondary_citations` (items 20, 21), which makes every claim the app shows
   traceable in an export.
7. Register the remaining Sievwright elements and the olivine Ni and Ca laws
   (items 22, 23).

## 5. September 29 implementation update

The current registry has 103 entries across 13 minerals. The earlier counts, recommended olivine flag, olivine reference subtraction, partial magnetite registration, missing workbook export and independent-only fitting statements are historical. Source-transcription verification does not establish natural-sample suitability or complete uncertainty validation.

### New modules

| Module | Responsibility |
| --- | --- |
| `diffusor/coefficients/literature.py` | Scalar and independent principal-axis laws with calibration/uncertainty metadata. |
| `diffusor/coefficients/accessories.py` | Accessory lattice laws separated by experimental mechanism and direction. |
| `diffusor/dataio/workbook.py` | Numeric Excel reports, profile/residual charts and provenance sheets. |
| `diffusor/dataio/figures.py` | Headless profile/residual and shared-duration PNG/SVG layouts. |
| `diffusor/dataio/joint.py` | Aggregate shared-duration export plus full per-profile provenance. |
| `diffusor/fitting/joint.py` | Shared-duration likelihood, strict uncertainties and joint degrees of freedom. |
| `diffusor/gui/joint_study.py` | Nonmodal study collection, background fitting, cancellation, plot and export. |
| `scripts/index_literature.py` | Hashes and extracts local PDFs into an ignored working index. |
| `tests/test_literature_expansion.py` | Equations, units, directions, coordinate mapping, conservative fluxes and exports. |
| `tests/test_joint_study.py` | Duration recovery, unit invariance, inconsistent profiles, cancellation and exports. |
| `diffusor/dataio/images.py` | Image and map readers (Pillow, optional tifffile, text grids, ENVI, raw), pixel size from metadata, grey/channel/colour-legend values. |
| `diffusor/dataio/image_profiles.py` | Perpendicular profiles after NIDIS `greyvalues.m`, cleaning with reason codes, statistics, Excel workbook out and back, composition table. |
| `diffusor/gui/image_extractor.py` | The image window: drawing tools, live extraction, overlays, workbook, hand-over to the main window. |
| `diffusor/gui/image_calibration.py` | Pixel size and value-to-composition map for an image profile. |
| `scripts/make_example_images.py` | Synthetic BSE TIFF and jet element map in examples/images/. |
| `scripts/make_icon.py` | Draws the application icon (PNG for the window, ICO for shortcuts). |
| `tests/test_image_profiles.py` | Readers, geometry, cleaning, colour legends, workbook round trip and the window. |

### New coefficient keys

| Key | Citation | Source location |
| --- | --- | --- |
| `ap_Dy_cherniak2000_in` | `cherniak2000apatite` | abstract p. 3871, REE silicate oxyapatite source |
| `ap_La_cherniak2000_in` | `cherniak2000apatite` | abstract p. 3871, REE silicate oxyapatite source |
| `ap_Nd_cherniak2000_in` | `cherniak2000apatite` | abstract p. 3871, REE silicate oxyapatite source |
| `ap_Nd_cherniak2000_out` | `cherniak2000apatite` | abstract p. 3871, synthetic Nd-doped apatite out-diffusion |
| `ap_Pb_cherniak1991` | `cherniak1991apatite` | abstract p. 1663; cm2/s and kcal/mol converted to SI |
| `ap_Sm_cherniak2000_implant` | `cherniak2000apatite` | abstract p. 3871, ion-implantation relaxation |
| `ap_Sr_cherniak1993` | `cherniak1993apatite` | p. 4657; abstract; D0 converted from cm2/s |
| `ap_Yb_cherniak2000_in` | `cherniak2000apatite` | abstract p. 3871, REE silicate oxyapatite source |
| `mnz_Pb_cherniak2004` | `cherniak2004monazite` | abstract p. 829 |
| `mt_Al_sievwright2020_1150` | `sievwright2020` | 5 and Table 5 |
| `mt_Co_sievwright2020_1150` | `sievwright2020` | 5 and Table 5 |
| `mt_Cr_sievwright2020_1150` | `sievwright2020` | 5 and Table 5 |
| `mt_Ga_sievwright2020` | `sievwright2020` | 5 and Table 5 |
| `mt_Hf_sievwright2020` | `sievwright2020` | 5 and Table 5 |
| `mt_In_sievwright2020` | `sievwright2020` | 5 and Table 5 |
| `mt_Lu_sievwright2020` | `sievwright2020` | 5 and Table 5 |
| `mt_Mn_sievwright2020_1150` | `sievwright2020` | 5 and Table 5 |
| `mt_Mo_sievwright2020` | `sievwright2020` | 5 and Table 5 |
| `mt_Nb_sievwright2020` | `sievwright2020` | 5 and Table 5 |
| `mt_Ni_sievwright2020` | `sievwright2020` | 5 and Table 5 |
| `mt_Sc_sievwright2020` | `sievwright2020` | 5 and Table 5 |
| `mt_Ta_sievwright2020` | `sievwright2020` | 5 and Table 5 |
| `mt_Ti_sievwright2020_1150` | `sievwright2020` | 5 and Table 5 |
| `mt_U_sievwright2020` | `sievwright2020` | 5 and Table 5 |
| `mt_V3+_sievwright2020` | `sievwright2020` | 5 and Table 5 |
| `mt_V4+_sievwright2020` | `sievwright2020` | 5 and Table 5 |
| `mt_Y_sievwright2020` | `sievwright2020` | 5 and Table 5 |
| `mt_Zn_sievwright2020` | `sievwright2020` | 5 and Table 5 |
| `mt_Zr_sievwright2020` | `sievwright2020` | 5 and Table 5 |
| `ol_Be_jollands2016` | `jollands2016be` | abstract and Arrhenius fits |
| `ol_Ca_coogan2005` | `coogan2005ca` | abstract and Arrhenius fits |
| `ol_Ni_petry2004` | `petry2004` | p. 4184, Fig. 6 fixed-fO2 fit |
| `ol_P_watson2015` | `watson2015p` | abstract p. 2053 |
| `qz_Ti_cherniak2007` | `cherniak2007quartz` | abstract p. 65; fitted c-direction |
| `rt_Hf_cherniak2007_a` | `cherniak2007rutile` | abstract p. 267 |
| `rt_Hf_cherniak2007_c` | `cherniak2007rutile` | abstract p. 267 |
| `rt_Zr_cherniak2007_c` | `cherniak2007rutile` | abstract p. 267 |
| `ttn_Sr_cherniak1995` | `cherniak1995titanite` | abstract p. 219 (visually checked) |
| `ttn_Zr_cherniak2006_c` | `cherniak2006titanite` | abstract p. 639 |
| `xtm_Dy_cherniak2006` | `cherniak2006xenotime` | abstract pp. 1-2; Table 1 for per-species temperature ranges |
| `xtm_Pb_cherniak2006` | `cherniak2006xenotime` | abstract pp. 1-2; Table 1 for per-species temperature ranges |
| `xtm_Sm_cherniak2006` | `cherniak2006xenotime` | abstract pp. 1-2; Table 1 for per-species temperature ranges |
| `xtm_Yb_cherniak2006` | `cherniak2006xenotime` | abstract pp. 1-2; Table 1 for per-species temperature ranges |
| `zrn_Dy_cherniak1997` | `cherniak1997zircon` | abstract p. 289, low-temperature RBS fits |
| `zrn_Pb_cherniak2001` | `cherniak2001zircon` | abstract p. 5 |
| `zrn_Sm_cherniak1997` | `cherniak1997zircon` | abstract p. 289, low-temperature RBS fits |
| `zrn_Ti_bloch2022_c` | `bloch2022zircon` | abstract p. 1 |
| `zrn_Ti_cherniak2007_perp_c` | `cherniak2007zircon` | abstract p. 470, 1-atm fit |
| `zrn_Yb_cherniak1997` | `cherniak1997zircon` | abstract p. 289, low-temperature RBS fits |

### Additional reference keys

- `cherniak1993apatite`: Cherniak, D. J. and Ryerson, F. J. (1993) A study of strontium diffusion in apatite using Rutherford backscattering spectroscopy and ion implantation. Geochimica et Cosmochimica Acta 57:4653-4662. https://doi.org/10.1016/0016-7037(93)90190-8 Used for the named calibration above.
- `cherniak1991apatite`: Cherniak, D. J. and Lanford, W. A. and Ryerson, F. J. (1991) Lead diffusion in apatite and zircon using ion implantation and Rutherford Backscattering techniques. Geochimica et Cosmochimica Acta 55:1663-1673. https://doi.org/10.1016/0016-7037(91)90137-T Used for the named calibration above.
- `cherniak2000apatite`: Cherniak, D. J. (2000) Rare earth element diffusion in apatite. Geochimica et Cosmochimica Acta 64:3871-3885. https://doi.org/10.1016/S0016-7037(00)00467-1 Used for the named calibration above.
- `cherniak2001zircon`: Cherniak, D. J. and Watson, E. B. (2001) Pb diffusion in zircon. Chemical Geology 172:5-24. https://doi.org/10.1016/S0009-2541(00)00233-3 Used for the named calibration above.
- `cherniak1997zircon`: Cherniak, D. J. and Hanchar, J. M. and Watson, E. B. (1997) Rare-earth diffusion in zircon. Chemical Geology 134:289-301. https://doi.org/10.1016/S0009-2541(96)00098-8 Used for the named calibration above.
- `cherniak2007zircon`: Cherniak, D. J. and Watson, E. B. (2007) Ti diffusion in zircon. Chemical Geology 242:470-483. https://doi.org/10.1016/j.chemgeo.2007.05.005 Used for the named calibration above.
- `bloch2022zircon`: Bloch, E. M. and Jollands, M. C. and Tollan, P. and others (2022) Diffusion anisotropy of Ti in zircon and implications for Ti-in-zircon thermometry. Earth and Planetary Science Letters 578:117317. https://doi.org/10.1016/j.epsl.2021.117317 Used for the named calibration above.
- `cherniak2004monazite`: Cherniak, D. J. and Watson, E. B. and Grove, M. and Harrison, T. M. (2004) Pb diffusion in monazite: A combined RBS/SIMS study. Geochimica et Cosmochimica Acta 68:829-840. https://doi.org/10.1016/j.gca.2003.07.012 Used for the named calibration above.
- `cherniak2006xenotime`: Cherniak, D. J. (2006) Pb and rare earth element diffusion in xenotime. Lithos 88:1-14. https://doi.org/10.1016/j.lithos.2005.08.002 Used for the named calibration above.
- `watson2015p`: Watson, E. B. and Cherniak, D. J. and Holycross, M. E. (2015) Diffusion of phosphorus in olivine and molten basalt. American Mineralogist 100:2053-2065. https://doi.org/10.2138/am-2015-5416 Used for the named calibration above.
- `jollands2016be`: Jollands, M. C. and Burnham, A. D. and O'Neill, H. St. C. and Hermann, J. and Qian, Q. (2016) Beryllium diffusion in olivine: A new tool to investigate timescales of magmatic processes. Earth and Planetary Science Letters 450:71-82. https://doi.org/10.1016/j.epsl.2016.06.028 Used for the named calibration above.
- `cherniak2007quartz`: Cherniak, D. J. and Watson, E. B. and Wark, D. A. (2007) Ti diffusion in quartz. Chemical Geology 236:65-74. https://doi.org/10.1016/j.chemgeo.2006.09.001 Used for the named calibration above.
- `cherniak2007rutile`: Cherniak, D. J. and Manchester, J. and Watson, E. B. (2007) Zr and Hf diffusion in rutile. Earth and Planetary Science Letters 261:267-279. https://doi.org/10.1016/j.epsl.2007.06.027 Used for the named calibration above.
- `cherniak1995titanite`: Cherniak, D. J. (1995) Sr and Nd diffusion in titanite. Chemical Geology 125:219-232. https://doi.org/10.1016/0009-2541(95)00074-V Used for the named calibration above.
- `cherniak2006titanite`: Cherniak, D. J. (2006) Zr diffusion in titanite. Contributions to Mineralogy and Petrology 152:639-647. https://doi.org/10.1007/s00410-006-0133-0 Used for the named calibration above.

### Mathematical and reporting changes

Olivine Fo-percent profiles map to coefficient XFe=1-C/100. Dohmen & Chakraborty equations 27–28 use the erratum's +3 (XFe − 0.1) (checked 29 September 2026; the 20 September reading of the printed +3 XFe was wrong and is reverted), so both entries are verified and TaMED is recommended again. Ca/Be principal functions are evaluated before direction-cosine projection. Fixed-temperature tables reject other temperatures and thermal Monte Carlo. No unmeasured tensor or parameter covariance is invented.

The numerical operator uses exact node-centred control volumes in plane/cylinder/sphere geometry and includes activity fluxes at boundary-adjacent faces. Mass uses the same weights. Zero external total flux conserves mass. Invalid grids, diffusivity arrays and exhausted step budgets raise errors. Callable boundaries force numerical integration; an isothermal supplied history uses its own temperature.

The joint objective is sum over profiles/points of ((observed-model)/sigma)^2, with one shared duration and N-1 joint degrees of freedom. All other settings remain fixed. Independent Gaussian errors are assumed. This is independent scalar transport with a common parameter, not coupled transport or posterior sampling. Per-profile chi-squared values are contributions to the joint objective.

Excel workbooks are numeric snapshots with charts and source/caveat sheets. JSON retains mechanism, ranges and coordinate transform. PNG/SVG reports show data, model, initial state and residuals. Joint studies export aggregate and complete individual reports. The corrected olivine synthetic example is regenerated with a documented fixed seed. Remaining roadmap items and library evidence are in LITERATURE_AUDIT.md.
