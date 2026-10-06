# Mutch et al. (2019): limits on a single-profile comparison

For `BORG14_OL_C2_P1`, `Data_S1.xlsx` (`Olivine_EPMA_data_Initial_Condi`) provides the measured `XFo (mol frac)` and `Inital_XFo (mol frac)` arrays. The extracted CSV retains these as `Fo_mol` and `Initial_Fo_mol`, multiplied by 100. `Table_S4.xlsx` (`Olivine_EBSD_angles`) gives the traverse angles `angle100P = 102.6415°`, `angle010P = 16.6905°`, and `angle001P = 100.7185°`.

`Data_S3.xlsx` (`Inversion_medians`, Al-based initial-condition row) reports median `T (°C) = 1233.0893507764`, `P (kbar) = 7.8278509629`, and `Fe3+/Fetotal = 0.1456491767` for this profile. This table does not report median oxygen fugacity. Deriving fO2 from the oxidation ratio requires a thermodynamic mapping and melt-composition inputs that are not supplied with this profile record, so a conditional Fe-Mg fit at the published medians cannot be run without adding an unsupported assumption.

The paper's published result is a joint Bayesian inversion of Fe-Mg, Ni, and Mn with study-specific diffusion regressions and parameter distributions. A single-species least-squares duration using a registry coefficient would be a different conditional sensitivity calculation, not a reproduction of that posterior. Keep this profile data-only until a source-supported fO2 mapping and the intended coefficient law are established.
