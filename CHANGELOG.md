# Changelog

What changed between releases, newest first. The text under a version heading
is what goes into the GitHub release, so it is written for someone deciding
whether to update. Commit messages hold the details.

## Unreleased

### Results that change

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
