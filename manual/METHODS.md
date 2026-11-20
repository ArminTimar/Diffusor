# Diffusor methods

This document describes the solvers, the uncertainty model and the registry of
diffusion laws. The checks against published values are in
[VERIFICATION.md](VERIFICATION.md), and the windows of the program in the
[user guide](USER_GUIDE.md).

## Numerical methods

Boundary conditions are selected by their physical meaning. A plateau that
continues corresponds to a semi-infinite medium, the assumption of the
error-function solution used by Ostorero et al. (2022) and Chamberlain et al.
(2014). A crystal rim in contact with the melt is a fixed concentration, the
open boundary of Costa et al. (2008, p. 555). A closed rim and the crystal
centre are zero-flux boundaries, the centre by symmetry (Crank 1975, section
4.3). As long as both plateaus are preserved, the far boundaries do not affect
the fitted time. Diffusor issues a warning when the diffusion front reaches the
end of the traverse.

Closed-form solutions from Crank (1975) are used where they are exact: the
semi-infinite step (eq. 2.14), the finite band (eq. 2.15), the plane sheet
(eq. 4.17), the cylinder (eq. 5.22), the sphere (eq. 6.18) and the
corresponding fractional-uptake series. Otherwise the diffusion equation is
solved with a Crank-Nicolson scheme on a uniform grid with the diffusion
coefficient evaluated at half-nodes (Crank 1975, section 8.4, and Dohmen, Faak
& Blundy 2017). Plane, cylindrical and spherical geometries share one
conservative operator, with the L'Hôpital limit at the centre of a cylinder or
sphere. The time step starts at the explicit stability limit and increases
geometrically, which resolves the initial transient and keeps a 100 kyr run to
a few hundred steps.

The diffusion coefficient may depend on composition, in which case it is
evaluated at each node from the evolving profile. For plagioclase trace
elements the flux includes the anorthite activity term of Costa et al. (2003),
-D C (A/RT) dX<sub>An</sub>/dx. The factor A is taken from one of the two sets
in Table 1 of Dohmen, Faak & Blundy (2017), chosen on the Model step. One is
Dohmen & Blundy (2014), using their 900 or 1200 °C column, whichever is nearer
the run temperature. The other is Bindeman et al. (1998), the set used by Costa
et al. (2003) and Druitt et al. (2012). The two differ in sign for Mg. A third
choice is the A<sub>Mg</sub> of Mutch et al. (2021), -34.099 kJ/mol, used by
DFENS for Mg. A/RT is evaluated at the temperature of each time step, and every
export states the set. For water in melt the profile itself is the variable D
depends on. The laws of Behrens et al. (2004) are proportional to the water
content in wt%, so D is evaluated node by node from C<sub>H2O</sub>.

Anisotropy is handled with the direction-cosine relation of Costa &
Chakraborty (2004) and the axial ratios published for each mineral. A linear
cooling path replaces D t by the time integral of D(T(t)) (Crank 1975,
section 7.2, and Lasaga 1983), with the buffer-referenced fO2 following the
temperature. The model profile is convolved with a Gaussian of the analytical
beam width before comparison with the data (Ganguly et al. 1988, Bradshaw &
Kent 2017), so that beam broadening is not interpreted as diffusion.

## Uncertainty propagation

The fitted time depends non-linearly on inputs that are correlated with each
other. Three correlations are treated explicitly.

**Temperature and oxygen fugacity.** Oxygen fugacity is usually known as an
offset from a mineral buffer, and the buffer depends on temperature. Diffusor
recalculates the buffer at every sampled temperature, so temperature and fO2
remain correlated.

**Arrhenius parameters.** ln D<sub>0</sub> and the activation energy Q of an
Arrhenius law are positively correlated through the regression that produced
them, and sampling them independently overestimates the spread of D at the
working temperature. Where a source gives the correlation, the parameters are
sampled jointly. The Sr and Ba laws of Grocolas et al. (2025) provide this. The
authors treat log D<sub>0</sub> and Q as perfectly correlated in their own Monte
Carlo, and Diffusor uses the same assumption. The six olivine regressions of
Mutch et al. (2019) and the eight olivine and eight plagioclase regressions of
Mutch et al. (2021) are sampled from their published covariance matrices. For
74 other laws log D is sampled at the working temperature with a stated
scatter. For 11 of them the scatter is published or derived from a statement in
the source. For the other 63 it is Diffusor's assumption, because the source
gives only errors of D<sub>0</sub> and Q without their covariance. The
coefficient page and the methods paragraph state which applies. Laws with
neither a covariance nor a scatter are held fixed by default (63 laws).
Independent sampling of D<sub>0</sub> and Q is never the default and is
available only to reproduce results obtained that way.

**Composition.** Composition enters both the diffusion coefficient and the
fitted profile, so measurement noise affects the result through both.

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

Diffusor contains 161 published diffusion laws. Each entry states whether it is
a tracer, an interdiffusion (exchange), a chemical (trace-element) or an
effective coefficient. Each entry also stores the citation, the equation as
printed in the source with its equation number, the calibration ranges in
temperature, pressure, oxygen fugacity and composition, the published parameter
uncertainties, and the table, page or figure from which the values were
transcribed. A run outside the calibration range of the selected law produces a
warning. Each equation is shown typeset on the coefficient page and exported as
LaTeX source with the results.

| mineral | entries | species |
| --- | --- | --- |
| Olivine | 29 | Fe-Mg, Fe and Mg tracers, Ca, Be, P, Ni, Mn |
| Orthopyroxene | 16 | Fe-Mg, Mg, Lu, Ce, Eu |
| Clinopyroxene | 6 | Fe-Mg, Ca-Mg |
| Plagioclase | 22 | Mg, Sr, Ba, Li, K, NaSi-CaAl |
| K-feldspar | 7 | Sr, Ba, Ti, Na-K |
| Magnetite | 48 | Ti, Fe, Fe-Ti and 21 trace elements |
| Garnet | 4 | Fe-Mg, Fe and Mg tracers, Ca-(Fe,Mg) (multicomponent sets in the [user guide](USER_GUIDE.md#multicomponent-and-isotope-profiles)) |
| Quartz | 2 | Ti (Cherniak et al. 2007, Jollands et al. 2020) |
| Rutile | 3 | Zr, Hf |
| Titanite | 2 | Sr, Zr |
| Apatite | 8 | Sr, Pb, rare earth elements |
| Zircon | 6 | Ti, Pb, rare earth elements |
| Monazite | 1 | Pb |
| Xenotime | 4 | Pb, rare earth elements |
| Silicate melt | 3 | H<sub>2</sub>O (Behrens et al. 2004: dacite, andesite, and rhyolite to andesite with an SiO<sub>2</sub> term) |

The recommended law for each mineral and species, where one is designated:

| mineral | species | recommended source |
| --- | --- | --- |
| Olivine | Fe-Mg | Dohmen & Chakraborty (2007), TaMED mechanism, with the composition term of the erratum |
| Orthopyroxene | Fe-Mg | Dias, Dohmen & Behrens (2025), two fO2 regimes |
| Clinopyroxene | Fe-Mg | Mueller et al. (2013) |
| Plagioclase | Mg | Van Orman, Cherniak & Kita (2014), and Audétat, Grocolas & Mutch (2026) with silica activity |
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

### Laws as used by published studies

43 entries are laws in the form in which a published study applied them. These
are the study's own regressions (Mutch et al. 2019, and the olivine and
plagioclase regressions of the DFENS code, Mutch et al. 2021), the laws of the
DIPRA code (Girona & Costa 2013), or a calibration with another composition
term, anisotropy factor or reference value (for example Mourey et al. 2023, Sato
et al. 2022, Winslow et al. 2026, Morgado et al. 2019). Each is cited to the
study, names the calibration it was taken from and how it differs, and is never
recommended. They allow published timescales to be recomputed with the law
their authors used, including forms that later calibrations replace.

### Notes on individual laws

Sievwright et al. (2020) measured 21 elements in magnetite at 1150 °C only.
Each row of their Table 5 is therefore an entry valid at that temperature, and
a Monte Carlo run that samples temperature is refused for these entries. Five
further entries (Ti, Mn, Co, Cr and Al) extrapolate these values to other
temperatures with activation energies from Van Orman & Crispin (2010). They
are labelled as Diffusor's own extrapolation, since Sievwright et al. did not
determine activation energies.

Sr and Ba in sanidine were measured on the same Or61 crystal with the same
method (Cherniak 1996, 2002). At 800 °C Ba diffuses about 1.7 log units more
slowly than Sr, so paired Sr and Ba profiles across one zone boundary indicate
whether the boundary has been modified by diffusion (Chamberlain et al. 2014,
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
agrees with Cherniak (2002) within 0.6 log units, and the older law is kept as
an alternative.

Fe-Mg in orthopyroxene: Dias & Dohmen (2024) and Dias, Dohmen & Behrens (2025)
refitted the experiments of Dohmen et al. (2016) with a temperature-dependent
composition term, added Fe-rich crystals, and found two regimes separated at
log fO2 = −10 (fO2 in Pa). For Fs10 at log fO2 = −7 the 2025 law agrees with
Dohmen et al. (2016) within 0.3 log units, but the two differ by up to one log
unit for Fe-rich orthopyroxene at 1100 °C or under reducing conditions. The
general Dohmen et al. (2016) law and the intermediate law of Dias & Dohmen
(2024) are therefore marked as superseded. Arc magmas generally lie above the
calibrated fO2 range, where the fO2<sup>1/4</sup> term is extrapolated.

Fe-Ti in magnetite has no direct successor to Freer & Hauptman (1978) or Aragon
et al. (1984), and these entries are not marked as superseded. For Ti and Fe
tracer diffusion, however, the data of Aggarwal & Dieckmann (2002) as
tabulated by Van Orman & Crispin (2010, Table 12) depend explicitly on both
temperature and fO2, are better constrained, and were used by Tomiya et al.
(2013) at Shinmoedake. They are the recommended entries. At 1150 °C they agree
with Sievwright et al. (2020) within half a log unit for Ti, Mn and Co.

## Limitations

* The model is one-dimensional. A time obtained from a 1-D profile through a
  3-D crystal is subject to geometric and sectioning effects (Shea et al. 2015,
  Krimer & Costa 2017), and every fit states this.
* A Monte Carlo run with a composition-dependent law and the numerical solver
  takes a few seconds per draw and core, so several hundred draws take
  minutes. The plot shows the progress, and **Stop** ends the run.
* The garnet D matrix assumes an ideal solution. Borinski et al. (2012) found
  that non-ideality changes retrieved coefficients by less than a factor of
  1.2 for most natural garnets. Their non-ideal matrix is not implemented.
* The fO2 term of Carlson (2006) is relative to the graphite-oxygen
  equilibrium, which Diffusor does not calculate. The offset from it is an
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
* Ni and Mn in olivine are available only as the laws of published studies
  (Mutch et al. 2019, 2021, and Girona & Costa 2013) and, for Ni, the Fo90,
  fixed-fO2 fit of Petry et al. (2004).
* The Lu law of Dias et al. (2025) gives its fO2 dependence relative to a
  reference fO2<sub>0</sub> that the paper does not define. Diffusor uses the
  IW buffer, because the paper describes the laws as valid along it, and the
  entry states that this is an interpretation.
* The plagioclase Mg law of Audétat et al. (2026) is taken from the accepted
  manuscript. Its silica activity is 1 (quartz saturation) unless changed from
  Python, and the interface does not offer this setting.

Corrections to results of earlier versions are listed in
[CHANGELOG.md](../CHANGELOG.md).
