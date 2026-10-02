# Changelog

What changed between releases, newest first. The text under a version heading
is what goes into the GitHub release, so it is written for someone deciding
whether to update. Commit messages hold the details.

## Unreleased

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
