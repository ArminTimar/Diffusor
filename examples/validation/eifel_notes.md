# Sundermeyer Eifel Fe-Mg comparison

Run the saved calculation from the repository root with:

```powershell
.\.venv\Scripts\python.exe examples/validation/check_eifel.py
```

It recreates `eifel_results.json` using only the tracked E41-4-1 profile CSV and input values transcribed into the script. It does not need the downloaded PDF or workbooks.

## Case and comparison

E41-4-1 is the clearest Eifel case for a one-boundary comparison. The paper’s Table 1 reports an Mg-Fe time of **291 days** with asymmetric uncertainty **−213/+103 days**. Online Resource 4 gives measured T=1159±23 °C, effective T=1101 °C, line orientation angles (81.9°, 9.2°, 85.6°) to a/b/c, and Fo core/rim entries 0.838/0.859. The profile CSV includes both an outer low-Fo rim and the inner core–high-Fo mantle boundary. The reconstruction keeps points at 30–135.3 µm to fit that inner boundary and exclude the separate outer-rim gradient.

| Model | Fitted time | RMSE | Relative to published 291 d |
|---|---:|---:|---:|
| Dohmen & Chakraborty (2007) TaMED law, 241 nodes | 457.07 d | 0.1893 mol% | 1.57× |
| Oeser et al. (2026) interdiffusion, 241 nodes | 172.98 d | 0.1894 mol% | 0.59× |

A 121-to-241 node refinement changes fitted time by less than 0.02% for both laws. The very similar residuals do not decide between the laws; this is a coefficient sensitivity comparison.

## What is reconstructed and what is assumed

The paper describes a finite 1-D DIPRA model with composition-dependent diffusivity and RMS curve matching. It describes the core-to-high-Fo-mantle step, and says the outermost low-Fo rim is a separate zone. The article does not publish E41-4-1’s exact fit window, fO2, DIPRA initial/boundary setup, or original source code. Online Resource 4’s `Dx` and `Dt` columns are spatial/time increments for the numerical search, not the published diffusion age; the age is in Table 1.

The script therefore labels this a **reconstruction, not an exact replication**. It uses a one-dimensional plane, the OR4 core/rim Fo values as the step endpoints, lets the interface position vary, and fixes Dirichlet plateau values at the cropped window ends. It uses the paper’s effective temperature and its stated similar 2 kbar pressure estimate for E41. The paper gives no fO2, so the comparison uses log10 fO2 = −5 in Pa only as a shared reference near the Oeser law’s stated 1e−5 Pa calibration condition. This is not an inferred Eifel fugacity. The original law’s time changes with the unknown fugacity.

The Oeser result is a sensitivity only: E41’s Fo composition (XFe about 0.14–0.16) is outside the law’s San Carlos XFe=0.085 calibration, and 2 kbar is outside its atmospheric-pressure calibration. The conditions/warnings are preserved in the JSON. Neither this case nor its assumed fO2 should be shipped as a GUI preset. A defensible exact preset would need the source’s original DIPRA setup or a documented fO2 and fit-window specification.

## Source trail

- Sundermeyer et al. (2020), *Timescales from magma mixing to eruption in alkaline volcanism in the Eifel volcanic fields, western Germany*, DOI [10.1007/s00410-020-01715-y](https://doi.org/10.1007/s00410-020-01715-y). Methods: pp. 12–14 (pressure and DIPRA approach); results and Table 1: pp. 16–17.
- Online Resource 2: measured Eifel line profiles, with Fo in mol% and distance in µm.
- Online Resource 4, row E41-4-1: temperatures, numerical increments, orientation, and compositional endpoints. These values are embedded in `check_eifel.py` so the check is runnable without the downloaded supplement.
