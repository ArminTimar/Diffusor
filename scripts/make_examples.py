"""Generate the synthetic example datasets in examples/.

Each file is a forward model run with a known time, so the answer is known and
the examples double as an end-to-end check of the whole pipeline.
Run from the project root:  python scripts/make_examples.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from diffusor.coefficients import Conditions, get                     # noqa: E402
from diffusor.constants import OXIDE_MOLAR_MASS, SEC_PER_DAY, SEC_PER_YEAR   # noqa: E402
from diffusor.fitting import DiffusionModel                            # noqa: E402
from diffusor.solvers import Geometry, InitialCondition, dirichlet, zero_flux  # noqa: E402
from diffusor.thermo import log_fo2_from_delta                         # noqa: E402

OUT = ROOT / "examples"
OUT.mkdir(exist_ok=True)
rng = np.random.default_rng(20260915)
ANSWERS = {}


def _model(coef_key, cond, ic, comp_key=None, bc=None, n_nodes=401, geometry="plane"):
    bc_l, bc_r = bc or (dirichlet(ic.params["C_left"]), dirichlet(ic.params["C_right"]))
    return DiffusionModel(coefficient=get(coef_key), conditions=cond, initial=ic,
                          geometry=Geometry(geometry), bc_left=bc_l, bc_right=bc_r,
                          comp_key=comp_key, composition_dependent=bool(comp_key),
                          n_nodes=n_nodes)


def x_fe_to_oxides(xfe, total_cations=2.0):
    """Turn X_Fe back into FeO and MgO wt% for a pyroxene-like stoichiometry."""
    fe = xfe * total_cations
    mg = (1 - xfe) * total_cations
    feo = fe * OXIDE_MOLAR_MASS["FeO"]
    mgo = mg * OXIDE_MOLAR_MASS["MgO"]
    sio2 = 2.0 * OXIDE_MOLAR_MASS["SiO2"]
    tot = feo + mgo + sio2
    return 100 * feo / tot, 100 * mgo / tot


# --- 1. orthopyroxene Fe-Mg, the Shinmoedake-style case ------------------------
def opx_example():
    T = 950 + 273.15
    cond = Conditions(T_K=T, P_Pa=1.5e8, log_fo2_bar=log_fo2_from_delta("NNO", 1.0, T),
                      X={"XFe": 0.28}, axis="c")
    ic = InitialCondition("step", {"x0": 0.0, "C_left": 0.32, "C_right": 0.22})
    m = _model("opx_FeMg_dias2025", cond, ic, comp_key="XFe")
    t = 1.5 * SEC_PER_YEAR
    x = np.linspace(-35, 35, 43)
    xfe = m.profile(t, x) + rng.normal(0, 0.003, x.size)
    feo, mgo = x_fe_to_oxides(xfe)
    pd.DataFrame({"Distance_um": x - x[0], "FeO_wt": np.round(feo, 3),
                  "MgO_wt": np.round(mgo, 3), "FeO_err": 0.15, "MgO_err": 0.20}
                 ).to_csv(OUT / "opx_femg_step.csv", index=False)
    ANSWERS["opx_femg_step.csv"] = (
        t, "Opx Fe-Mg, Dias, Dohmen & Behrens (2025), 950 C, NNO+1, 150 MPa, //c, X_Fe about 0.22-0.32. "
           "Map FeO_wt and MgO_wt with oxides FeO/MgO and mode A/(A+B).")


# --- 2. clinopyroxene Fe-Mg -----------------------------------------------------
def cpx_example():
    T = 1100 + 273.15
    cond = Conditions(T_K=T, P_Pa=2.0e8, X={"XFe": 0.15})
    ic = InitialCondition("step", {"x0": 0.0, "C_left": 0.20, "C_right": 0.11})
    m = _model("cpx_FeMg_muller2013", cond, ic)
    t = 45 * SEC_PER_DAY
    x = np.linspace(-12, 12, 41)
    xfe = m.profile(t, x) + rng.normal(0, 0.002, x.size)
    feo, mgo = x_fe_to_oxides(xfe)
    pd.DataFrame({"Distance_um": x - x[0], "FeO_wt": np.round(feo, 3),
                  "MgO_wt": np.round(mgo, 3), "FeO_err": 0.12, "MgO_err": 0.15}
                 ).to_csv(OUT / "cpx_femg_step.csv", index=False)
    ANSWERS["cpx_femg_step.csv"] = (
        t, "Cpx Fe-Mg, Mueller et al. (2013), 1100 C, 200 MPa, [001]. No fO2 dependence.")


# --- 3. plagioclase Mg with an anorthite gradient -------------------------------
def plag_example():
    T = 900 + 273.15
    xan = np.linspace(0.62, 0.45, 81)
    x = np.linspace(0, 400, 81)
    cond = Conditions(T_K=T, P_Pa=1.0e8, X={"XAn": 0.55})
    ic = InitialCondition("step", {"x0": 200.0, "C_left": 180.0, "C_right": 95.0})
    m = _model("plag_Mg_vanorman2014", cond, ic, n_nodes=301)
    t = 20 * SEC_PER_YEAR
    mg = m.profile(t, x) + rng.normal(0, 2.0, x.size)
    pd.DataFrame({"Distance_um": x, "Mg_ppm": np.round(mg, 1), "XAn": np.round(xan, 4),
                  "Mg_err": 4.0}).to_csv(OUT / "plag_mg_an.csv", index=False)
    ANSWERS["plag_mg_an.csv"] = (
        t, "Plagioclase Mg, Van Orman et al. (2014), 900 C, An45-62. Map Mg_ppm with mode 'A'. "
           "The anorthite column is supplied so the activity term of Costa et al. (2003) can be "
           "switched on.")


# --- 4. olivine Fe-Mg -----------------------------------------------------------
def olivine_example():
    T = 1150 + 273.15
    cond = Conditions(T_K=T, P_Pa=1.0e5, log_fo2_bar=log_fo2_from_delta("FMQ", -1.0, T),
                      X={"XFe": 0.15}, axis="c")
    ic = InitialCondition("step", {"x0": 0.0, "C_left": 0.20, "C_right": 0.12})
    m = _model("ol_FeMg_dohmen_chakraborty2007_tamed", cond, ic, comp_key="XFe")
    t = 120 * SEC_PER_DAY
    x = np.linspace(-150, 150, 61)
    xfe = m.profile(t, x) + rng.normal(0, 0.0025, x.size)
    fo = 100 * (1 - xfe)
    pd.DataFrame({"Distance_um": x - x[0], "Fo_mol": np.round(fo, 3), "Fo_err": 0.25}
                 ).to_csv(OUT / "olivine_fo.csv", index=False)
    ANSWERS["olivine_fo.csv"] = (
        t, "Olivine Fe-Mg, Dohmen & Chakraborty (2007) TaMED, 1150 C, FMQ-1, //[001]. "
           "Map Fo_mol with mode 'A' and remember that the coefficients of this entry are "
           "flagged as UNVERIFIED.")


# --- 5. titanomagnetite Ti, the Shinmoedake case of Tomiya et al. (2013) --------
def magnetite_example():
    T = 950 + 273.15
    cond = Conditions(T_K=T, P_Pa=1.0e5, log_fo2_bar=-11.0, X={"xTi": 0.1})
    ic = InitialCondition("step", {"x0": 0.0, "C_left": 6.5, "C_right": 4.2})
    m = _model("mt_Ti_vanorman_crispin2010", cond, ic, n_nodes=301)
    t = 8 * SEC_PER_DAY
    x = np.linspace(-50, 50, 41)
    ti = m.profile(t, x) + rng.normal(0, 0.05, x.size)
    pd.DataFrame({"Distance_um": np.round(x - x[0], 4), "TiO2_wt": np.round(ti, 3),
                  "TiO2_err": 0.08}).to_csv(OUT / "magnetite_ti.csv", index=False)
    ANSWERS["magnetite_ti.csv"] = (
        t, "Titanomagnetite Ti, Van Orman & Crispin (2010) Table 12, 950 C, log fO2 = -11, "
           "X_Usp = 0.3 (x_Ti = 0.1): the conditions Tomiya et al. (2013) used for the 2011 "
           "Shinmoedake eruption. Map TiO2_wt with mode 'A'.")


# --- 6. greyscale profile plus microprobe anchors ------------------------------
def greyscale_example():
    T = 1000 + 273.15
    cond = Conditions(T_K=T, P_Pa=2.0e8, X={"XFe": 0.15})
    ic = InitialCondition("step", {"x0": 0.0, "C_left": 0.19, "C_right": 0.12})
    m = _model("cpx_FeMg_muller2013", cond, ic)
    t = 3.0 * SEC_PER_YEAR
    x = np.linspace(-15, 15, 301)
    xfe = m.profile(t, x)
    grey = 90.0 + 620.0 * xfe + rng.normal(0, 0.8, x.size)     # linear BSE response
    pd.DataFrame({"Distance_um": np.round(x - x[0], 4), "GreyValue": np.round(grey, 2),
                  "Grey_sd": 0.8}).to_csv(OUT / "cpx_greyscale.csv", index=False)
    xi = np.array([-13.0, -7.0, 0.0, 7.0, 13.0])
    probe_xfe = np.interp(xi, x, xfe) + rng.normal(0, 0.002, xi.size)
    probe_grey = np.interp(xi, x, grey)
    pd.DataFrame({"Distance_um": np.round(xi - x[0], 4),
                  "GreyValue": np.round(probe_grey, 2),
                  "XFe_probe": np.round(probe_xfe, 5)}
                 ).to_csv(OUT / "cpx_greyscale_anchors.csv", index=False)
    ANSWERS["cpx_greyscale.csv"] = (
        t, "BSE grey values across a cpx zone boundary, with five microprobe anchor points in "
           "cpx_greyscale_anchors.csv. Calibrate first (the true response is "
           "X_Fe = (grey - 90)/620), then fit with Mueller et al. (2013) at 1000 C.")


# --- 7. sanidine Ba, at Bishop Tuff temperatures --------------------------------
def sanidine_example():
    T = 790 + 273.15
    cond = Conditions(T_K=T, P_Pa=1.0e5)
    ic = InitialCondition("step", {"x0": 0.0, "C_left": 4800.0, "C_right": 1600.0})
    m = _model("kfs_Ba_cherniak2002", cond, ic, n_nodes=401)
    t = 5000 * SEC_PER_YEAR
    x = np.linspace(-10, 10, 81)
    ba = m.profile(t, x) + rng.normal(0, 90.0, x.size)
    pd.DataFrame({"Distance_um": np.round(x - x[0], 4), "Ba_ppm": np.round(ba, 0),
                  "Ba_err": 90.0}).to_csv(OUT / "sanidine_ba.csv", index=False)
    ANSWERS["sanidine_ba.csv"] = (
        t, "Sanidine Ba, Cherniak (2002), 790 C, a bright Ba-rich rim against a darker core. "
           "790 C lies inside the 753-815 C range Chamberlain et al. (2014) used for Bishop Tuff "
           "sanidine. Map Ba_ppm with mode 'A'.")


def main():
    opx_example()
    cpx_example()
    plag_example()
    olivine_example()
    magnetite_example()
    greyscale_example()
    sanidine_example()
    lines = ["# Example datasets", "",
             "Two files here are **real measurements** and are not produced by this script:",
             "",
             "* `plagioclase_santorini_druitt2012.csv` -- Druitt et al. (2012) Nature 482:77-80,",
             "  Supplementary Table 1, plagioclase S82-30A 12.",
             "* `opx_kizimen_ostorero2022.csv` -- Ostorero et al. (2022) Commun. Earth Environ.",
             "  3:290, Supplementary Data 2, orthopyroxene K9_L10C4. Regenerate it with",
             "  `python scripts/extract_kizimen.py <folder with the supplementary files>`.",
             "",
             "The files below are **synthetic**: each was produced by running Diffusor's forward",
             "model with a known time, then adding Gaussian noise. They exist so the fitting,",
             "the Monte Carlo and the export can be checked against a known answer.",
             "Regenerate them with `python scripts/make_examples.py`.", "",
             "| file | true time | conditions |", "| --- | --- | --- |"]
    from diffusor.thermo.units import human_time
    for f, (t, note) in ANSWERS.items():
        lines.append(f"| `{f}` | **{human_time(t)}** | {note} |")
    lines += ["", "A fit will not return the true time exactly: the noise, and for the",
              "composition-dependent cases the difference between the fitting grid and the",
              "generating grid, shift it by a few per cent. Recovering the true value to",
              "within the Monte Carlo interval is the test that matters.", ""]
    (OUT / "README.md").write_text("\n".join(lines), encoding="utf-8")
    print("wrote", len(ANSWERS) + 1, "example files to", OUT)
    for f, (t, _) in ANSWERS.items():
        print(f"  {f:<28s} true time {human_time(t)}")


if __name__ == "__main__":
    main()
