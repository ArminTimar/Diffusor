# Example datasets

Two files here are **real measurements** and are not produced by this script:

* `plagioclase_santorini_druitt2012.csv` -- Druitt et al. (2012) Nature 482:77-80,
  Supplementary Table 1, plagioclase S82-30A 12.
* `opx_kizimen_ostorero2022.csv` -- Ostorero et al. (2022) Commun. Earth Environ.
  3:290, Supplementary Data 2, orthopyroxene K9_L10C4; regenerate it with
  `python scripts/extract_kizimen.py <folder with the supplementary files>`.

The files below are **synthetic**: each was produced by running Diffusor's forward
model with a known time, then adding Gaussian noise. They exist so the fitting,
the Monte Carlo and the export can be checked against a known answer.
Regenerate them with `python scripts/make_examples.py`.

| file | true time | conditions |
| --- | --- | --- |
| `opx_femg_step.csv` | **1.5 yr** | Opx Fe-Mg, Dias, Dohmen & Behrens (2025), 950 C, NNO+1, 150 MPa, //c, X_Fe about 0.22-0.32. Map FeO_wt and MgO_wt with oxides FeO/MgO and mode A/(A+B). |
| `cpx_femg_step.csv` | **45 d** | Cpx Fe-Mg, Mueller et al. (2013), 1100 C, 200 MPa, [001]. No fO2 dependence. |
| `plag_mg_an.csv` | **20 yr** | Plagioclase Mg, Van Orman et al. (2014), 900 C, An45-62. Map Mg_ppm with mode 'A'. The anorthite column is supplied so the activity term of Costa et al. (2003) can be switched on. |
| `olivine_fo.csv` | **120 d** | Olivine Fe-Mg, Dohmen & Chakraborty (2007) TaMED, 1150 C, FMQ-1, //[001]. Map Fo_mol with mode 'A' and remember that the coefficients of this entry are flagged as UNVERIFIED. |
| `magnetite_ti.csv` | **8 d** | Titanomagnetite Ti, Van Orman & Crispin (2010) Table 12, 950 C, log fO2 = -11, X_Usp = 0.3 (x_Ti = 0.1): the conditions Tomiya et al. (2013) used for the 2011 Shinmoedake eruption. Map TiO2_wt with mode 'A'. |
| `cpx_greyscale.csv` | **3 yr** | BSE grey values across a cpx zone boundary, with five microprobe anchor points in cpx_greyscale_anchors.csv. Calibrate first (the true response is X_Fe = (grey - 90)/620), then fit with Mueller et al. (2013) at 1000 C. |

A fit will not return the true time exactly: the noise, and for the
composition-dependent cases the difference between the fitting grid and the
generating grid, shift it by a few per cent. Recovering the true value to
within the Monte Carlo interval is the test that matters.
