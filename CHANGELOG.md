# Changelog

What changed between releases, newest first. The text under a version heading
is what goes into the GitHub release, so it is written for someone deciding
whether to update. Commit messages hold the details.

## 0.2.0 (unreleased)

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
  the Monte Carlo actually fitted, not the model at the ends of the time
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
