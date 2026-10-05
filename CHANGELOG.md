# Changelog

What changed between releases, newest first. The text under a version heading
is what goes into the GitHub release, so it is written for someone deciding
whether to update. Commit messages hold the details.

## Unreleased

### Fixed

- **One best-fit time on screen.** After a Fit, a change of settings and a Monte
  Carlo run, the sidebar could show the earlier fit (for example 2.8 yr) while
  the plot legend showed the best fit of the Monte Carlo run (338 d). The run
  refits the data with the current settings, so the window now takes that fit
  as its best fit, and the sidebar names the median of the draws separately.
  A setting changed after a result was computed now says so in the sidebar
  until the next Fit.
- **Equations are sharp on scaled displays.** The typeset laws were drawn at a
  fixed resolution and stretched by the window, so they looked soft, and small
  exponents could disappear. They are now drawn at the screen's own pixel
  density and a little smaller, and long laws are broken over several lines at
  plus and minus signs instead of being shrunk.
- **Beam preset without a mineral in its hint.** "Microprobe, defocused beam"
  explained itself with feldspar studies, which was confusing for the
  orthopyroxene example. It is now "Microprobe, stated beam diameter" and asks
  for the diameter given in the study's methods.
- **Source notes say what the papers state.** Notes no longer call a paper
  inconsistent. Where two places in a paper give different values, the note says
  what each gives and which one Diffusor uses. The Dias et al. (2025) opx law
  uses 900-1100 C, the interval of its Section 4.3. Reference notes no longer
  name local files.

## 0.3.0 (5 October 2026)

### Fixed

- **Examples set the analytical resolution.** Loading an example left the
  resolution on whatever the previous profile used. Each example now selects
  its own. The Kizimen microprobe traverse uses the 2 um focused beam stated
  by Ostorero et al. (2022), a beam sigma of 0.5 um. The synthetic examples
  were made without beam broadening and load with no correction, or with the
  BSE preset for the two that mimic BSE profiles. The "Where these values come
  from" page lists the resolution and its source.
- **Taskbar icon on Windows.** The running window showed a default icon beside
  the pinned Diffusor icon. The window now sets its own icon and the shell's
  relaunch properties, so its taskbar button and any pin made from it show
  Diffusor's icon.
- **Notes, citations and labels checked against the papers.** The laws,
  constants and methods were compared with the rendered pages of their
  sources. Values were kept unless listed under "Results that change"; the
  text around them was corrected:
  - Equation numbers from Crank (1975) were wrong in the code and in the
    exported Methods text (plane flow is eq. 1.5, Crank-Nicolson is section
    8.5, the time-dependent D substitution is eqs 7.2-7.3). The exported text
    no longer quotes equation numbers from the appendix of Dohmen et al.
    (2017), which could not be checked, and calls the finite-volume scheme
    Diffusor's own implementation.
  - Parameters whose source prints a plus-minus without saying what it is now
    read "level not stated, sampled as 1σ" instead of "1σ" (49 parameters in
    23 laws). Sampling is unchanged.
  - The fO2 term of the Ganguly & Tazzoli (1994) orthopyroxene law is
    attributed to Allan et al. (2013), who added it; Ganguly & Tazzoli only
    speculate about it. Counting all Fe as Fe2+ in X_Fe is stated to be
    Diffusor's convention, not that of Dohmen & Chakraborty or Mueller et al.
  - The covariance used for Sr and Ba in plagioclase (Grocolas et al. 2025) is
    described as Diffusor's reading of one sentence of the paper, which prints
    no matrix. The Ti law for K-feldspar is held fixed because the paper
    prints no covariance of D0 and Q, not because it gives no uncertainty
    for D0. The Petry et al. (2004) Ni note now shows the errors the paper
    prints.
  - Aluminium in magnetite (Van Orman & Crispin 2010, Table 12) is labelled
    chemical, not tracer, since it comes from interdiffusion data; implanted
    Sm in apatite is labelled tracer. Pressure and transport kind are marked as
    assumed wherever the authors do not state them.
  - Composition definitions and the mineral notes cite the third edition of
    Deer, Howie & Zussman (2013) by printed page, and the second edition is no
    longer cited. The plagioclase formula is NaAlSi3O8-CaAl2Si2O8 and the
    titanite formula CaTi(SiO4)(O,OH,F), as printed there. The quartz note
    gives the book's limits (alpha-quartz up to 573 C, beta-quartz from 573 to
    870 C) without a pressure, and the xenotime note says that the book has no
    data sheet for it. The Julian year cites the IAU Style Manual, and the gas
    constant is the exact product of the Avogadro and Boltzmann constants (a
    relative change of 1.8e-11).
  - The Laki olivine example: its temperature (1150 ± 30 C) and oxygen
    fugacity (FMQ-1 ± 0.5) were checked against Hartley et al. (2016, p. 61),
    and the pressure, the traverse direction, the composition range and the
    120 d duration are now stated to be Diffusor's choices, not the paper's.
    The olivine note no longer says that Hartley et al. confirmed the
    anisotropy factor of 6 in natural crystals. They quote it from Nakamura &
    Schmalzried (1983) and Dohmen & Chakraborty (2007).
  - The Sievwright et al. (2020) supplement gives the 1σ of the individual log D
    values (median 0.12) and none for the fitted constants, so the 0.2 log
    units sampled for their magnetite laws stay a stated assumption.

### Results that change

- **Monte Carlo defaults.** Five laws (Mg tracer in orthopyroxene, Schwandt et
  al. 1998, three axes; Ti in K-feldspar; Ca in olivine) sampled log D0 and Q
  independently when the coefficient was left to the default, which spreads
  log D by orders of magnitude. They are now held fixed by default.
  Ganguly & Tazzoli (1994) is sampled with 1.0 log unit (was 0.5), the
  standard error the paper states, and Dias et al. (2025) with 0.34 (was 0.2),
  the average misfit the paper states.
- **Dias & Dohmen (2024) orthopyroxene law.** The composition term used
  X_Fe - 0.09 instead of the paper's X_Fe - 0.1 (eqs 7 and 12). D changes by
  0.02 to 0.04 log units.
- **Plagioclase activity term.** The anorthite activity factor A was the
  Dohmen & Blundy (2014) 1200 C value at every temperature. The 900 C column is
  now used when the run temperature is nearer 900 C, A/RT follows the
  temperature along a cooling path, and the Bindeman et al. (1998) set used by
  Costa et al. (2003) and Druitt et al. (2012) can be chosen on the Model step.
  Refit plagioclase trace-element results that used the anorthite column.
- **Santorini example removed.** It modelled the Druitt et al. (2012) traverse
  with a different Mg law, initial profile and activity factor from the paper.
- **Aragon et al. (1984)** is now marked as not checked against the original
  paper; its expression comes from Van Orman & Crispin (2010) Table 11.
- **Titanomagnetite between x_Ti 0 and 0.2.** Van Orman & Crispin (2010)
  print Table 12 for x_Ti 0 and 0.2 only. Diffusor interpolated the total D
  between the two; it now interpolates the vacancy and the interstitial branch
  separately, each log-linearly in x_Ti, the form that Aggarwal & Dieckmann
  (2002) find for the concentrations of the defects that carry the diffusion.
  No published source gives D between the two rows, so this is Diffusor's own
  construction and the coefficient notes say so. Against the tracer
  coefficients that Aggarwal & Dieckmann measured at x_Ti 0.1 (1200 and
  1300 C) it has an rms error of 0.30 log units, and interpolating the total D
  is as close (0.29), so those data do not decide between the two. Using their
  x_Ti 0.1 and 0.3 values as extra interpolation nodes was tested and not
  adopted: each node rests on two temperatures and on activation energies the
  authors call unreliable, and at 900-1000 C the nodes would change D by up
  to a factor of 40. No measurement at intermediate x_Ti exists below 1200 C,
  so at magmatic temperatures the real uncertainty is larger than 0.3 log
  units. The end members are unchanged. For Ti, Fe, Mn and Co at x_Ti in
  between, D is lower by up to
  0.67 log units (1173-1573 K, log fO2 -14 to -6). At 950 C, log fO2 -11 and
  x_Ti 0.1: Ti 4.35e-16 to 4.15e-16 m2/s, Fe 9.40e-15 to 4.17e-15. The
  Shinmoedake values of Tomiya et al. (2013) for Ti and Fe at 950 and 900 C are
  now reproduced to within 5 %; Fe was 2.1 times too high at 950 C before.
  Refit titanomagnetite results. The bundled titanomagnetite example was made
  with the old D and has been regenerated with the new one; the old file would
  now give a time about 5 % longer than its true 8 days.
- **Orthopyroxene, Dohmen et al. (2016).** The scatter used by the default
  Monte Carlo is 0.2 log units, the overall uncertainty the paper states, not
  0.1. The spread of the log timescale for the two laws (Fs9 and Fs1) about
  doubles.
- **Uncertainty of ratio profiles built from oxides.** The uncertainty of
  Fe/(Fe+Mg) and similar profiles was propagated through the wt% columns while
  the profile itself is calculated from cation moles. It now follows the cation
  moles. The change is up to 11 % where the oxides differ in cations per
  formula unit or in relative error, and close to zero for FeO and MgO with
  equal relative errors. A difference profile (A-B) now has a propagated
  uncertainty instead of that of the first column. Fit weights, chi-squared and
  slightly the best-fit time follow.
- **Oxygen buffers.** All rows of Frost (1991) Table 1 are used, chosen by
  temperature; one high-temperature row per buffer was used before. FMQ below
  the alpha-beta quartz transition (573 C at 1 bar) and HM below 682 C change:
  FMQ at 500 C is -23.874, was -23.725; HM at 600 C is -14.841, was -14.876.
  IW, WM, NNO and the other temperatures are unchanged. `IQF` now means
  quartz-iron-fayalite, which lies 0.7 to 1.1 log units below the iron-wustite
  value it used to return between 800 and 1100 C, and QIF is in the buffer
  list. A temperature outside the range the source prints now raises a warning
  and still returns the extrapolated value.
- **Sharp initial step on the numerical route.** A fitted interface position
  between two grid nodes left the discretised step up to half a grid spacing
  off. Each node now takes the average over its cell, so fits with a free
  interface position and the numerical solver are more accurate. Fits where the
  interface falls exactly on a node are unchanged.
- **Fit statistics.** A fit with as many fitted parameters as points reports a
  reduced chi-squared of NaN and a warning; it divided by 1 before. A NaN,
  infinite or zero uncertainty is replaced by the median of the valid ones and
  a warning says so. A greyscale calibration with as many anchors as
  coefficients reports its own uncertainty as unknown, with a warning, instead
  of a silent zero.
- **Validity ranges and warnings.** D does not change, but warnings do:
  - Magnetite Table 12 laws warn outside the temperatures of the data behind
    them. The Fe, Co, Mn and Ti rows reproduce when refitted to the tracer
    coefficients of Aggarwal & Dieckmann (2002), so their windows are those
    of that data, 1373-1573 K for Fe, Co and Mn and 1473-1573 K for Ti. Cr
    (1483-1683 K) and Al (1553-1773 K) keep the windows printed by Van Orman &
    Crispin (all were 1373-1673 K). 900-1000 C always warns.
    The Aggarwal and Dieckmann entries, Aragon et al. (990-1220 C) and Freer &
    Hauptman (888-1034 C, x_Ti 0-0.05) follow their data.
  - Sr and Ba in plagioclase (Grocolas et al. 2025) warn above 1150 C, the
    hottest run used; Giletti & Casserly (1994) covers 550-1300 C and An0.6 to
    An95.6; Cherniak & Watson (1994) warns above An67; Ba in K-feldspar is
    828-1075 C; zircon rare earths and Hf in rutile have per-element ranges.
  - The fO2 range of the Mueller et al. (2013) clinopyroxene law was compared in
    the wrong unit and was off by 5 log units. The olivine PED law now warns
    above 1e-10 Pa, and the Dimanov & Sautter (2000) clinopyroxene law has an
    fO2 range. The Sievwright et al. (2020) range moved by 0.12 log units to the
    O'Neill (1987) FMQ their paper uses.
  - A pressure of exactly 1 atm no longer warns for the laws calibrated at
    1 atm.

- **Titanomagnetite example in the interface.** Loading the example replaced
  its x_Ti of 0.1 by the mean of the TiO2 profile, clamped to 1, which made D
  about 250 times too large: the fit gave about 0.02 days instead of 8. A
  bundled example now keeps its own representative composition.
- **Representative composition from your own profile.** Loading a file or
  pressing *Guess* wrote the profile mean into the representative composition,
  whatever the profile was. A TiO2 wt% traverse set x_Ti to 1 (D about 250
  times too large), and a ppm profile did the same to X_An or x_Ti. The mean is
  now used only when the profile is itself the host composition (Fe-Mg and the
  other exchange pairs, olivine through the chosen coordinate), and it is
  re-checked when you change the mineral or species. Refit any magnetite or
  trace-element result whose law needs a composition, and check the
  *Representative value* on the Conditions step.

### New

- **Typeset equations.** Every law is shown as a typeset equation on the
  coefficient pages and the Coefficient step, and its LaTeX source is written
  to the JSON, workbook and methods text.
- **Source of every uncertainty.** Each sampled log D scatter states whether
  it is published, derived from a statement in the source, or assumed.
- **Multicomponent garnet profiles.** File > Multicomponent and isotope study
  fits one duration to Fe, Mg, Mn and Ca profiles together. The diffusion
  matrix is that of an ideal ionic solution and is recalculated from the
  tracer coefficients at every node and step, so the cross terms (uphill
  diffusion) are included. Tracer sets: Carlson (2006), Chakraborty & Ganguly
  (1992), and both with the Mn law of Chen & Chu (2024).
- **Isotope profiles.** The same window fits a duration to a concentration
  profile together with its delta values: Li isotopes diffusing separately
  with D proportional to m^-beta, or the seven Fe and Mg isotopes of olivine
  as one coupled exchange (Oeser et al. 2026). beta values from Oeser et al.
  (2026) and Richter et al. (2014, 2017) are offered with the model each was
  fitted with.
- **New laws.** Fe and Mg tracer diffusion in olivine and the Fe-Mg
  interdiffusion they imply (Oeser et al. 2026, a, b and c axes, 1100-1250
  C); Na-K interdiffusion in K-feldspar normal to (001) and (010) at X_Or 0.92
  and 0.98 (Schaffer et al. 2014); garnet as a mineral, with the Fe and Mg
  tracer laws of Borinski et al. (2012) and the Fe-Mg interdiffusion computed
  from them. 113 laws for 14 minerals.
- **Tracer or exchange.** Every law now states whether it is a tracer, an
  interdiffusion (exchange), a chemical (trace-element) or an effective
  coefficient; 40 laws had no stated kind before. A fit with a tracer
  coefficient warns that major-element zoning relaxes by interdiffusion, which
  can differ from either tracer coefficient.
- **Recent profiles on the Data step.** The last ten files you loaded are
  listed in the panel that shows the profile; one button switches between the
  two, and loading a file brings the profile to the front. A recent file whose
  columns have not changed opens with the columns you chose last time, without
  the column dialog. If the columns changed, the dialog opens with your
  previous choices filled in. Files that have moved are shown as not found.
- **The file dialog opens in the folder you used last.**

### Changed

- **The coefficient is chosen before the conditions.** The steps are now Data,
  Mineral, Coefficient, Conditions, Model, Uncertainty. The law decides which
  host composition is needed and whether pressure and fO2 enter at all, so
  those are now asked for after it, and inputs the law does not use are
  disabled. Pressure stays active when fO2 is given relative to a buffer,
  because the buffer moves with pressure.
- **Settings moved to where they belong.** The analytical resolution (beam
  width, distance scale error) is on the Data step with the profile it
  describes. The linear cooling path is on the Conditions step next to the
  temperature. The host composition and *D follows the composition* moved
  from the Mineral step to the Conditions step.
- **Titanomagnetite example, plateaus.** The traverse ends before the profile
  has flattened, so the example now fits the plateaus. Held at the values of
  the outer points they sit inside the true plateaus and the fit gave 5.3 days;
  with the plateaus fitted it gives 7.4 days for a true 8. The same applies to
  your own data: tick *Fit the plateaus* when a profile still slopes at the
  ends of the traverse.

## 0.2.0 (1 October 2026)

The first published release.

### Results that change

Refit anything computed with these before this release:

- **Olivine Fe-Mg, Dohmen & Chakraborty (2007).** Between 20 and 29 September
  2026 the composition term was the misprinted 3 XFe instead of the erratum's
  3 (XFe - 0.1), which made olivine Fe-Mg times about half as long as they
  should be.
- **Olivine Fe-Mg, Chakraborty (1997).** D0 was 1.0e-9 m²/s instead of the
  paper's 5.38e-9 m²/s, so D was 5.4 times too low and times 5.4 times too
  long. The law is now refused along axes other than [001] and warns at fO2
  values other than 1e-12 bar, the only conditions it was fitted for.
- **Plagioclase NaSi-CaAl, Grove et al. (1984).** The entry was a placeholder
  about 15 times too slow. It now reproduces the paper's Fig. 3 line.
- **Monte Carlo with a cooling path.** Sampled temperatures used to change only
  fO2 and left the path itself at its nominal temperatures, so most of the
  temperature uncertainty was lost. Each draw now shifts the whole path.
- **Profile band.** The band under the fit is now the spread of the profiles
  the Monte Carlo draws fitted, not the model at the ends of the time
  interval.

### New

- A six-step interface (Data, Mineral, Conditions, Model, Coefficient,
  Uncertainty) with a results view, a live Monte Carlo view and boundary
  conditions described in words.
- 103 diffusion coefficients for 13 minerals, every one checked against its
  primary source, including the 2024-2026 calibrations for opx Fe-Mg (Dias et
  al. 2025), plagioclase Sr and Ba (Grocolas et al. 2025), Mg (Audétat et al.
  2026) and Li (Pohl et al. 2024), and K-feldspar, quartz and accessory phases.
  Superseded laws are kept, flagged and demoted.
- Profiles read from BSE images and element maps, perpendicular to a drawn
  boundary, after the NIDIS `greyvalues.m` script.
- Shared-duration fits of several profiles.
- Monte Carlo draws run in parallel. The same seed gives the same answer on
  any number of cores.
- Measured examples from Kizimen (Ostorero et al. 2022), reproduced within the
  published uncertainty, and Santorini (Druitt et al. 2012).
- Excel workbook exports with charts, calibration metadata and methods.
- Help > All references, with search and BibTeX export.
- Double-click starts for Windows and macOS, an application icon, and a notice
  in the app when a newer release is out.

## 0.1.0 (15 September 2026, not released)

The first working version: Crank analytical solutions, a Crank-Nicolson
solver, 34 coefficients for olivine, pyroxenes, plagioclase and magnetite, and
Monte Carlo error propagation with buffer-correlated fO2.
