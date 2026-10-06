"""Re-fit selected Sundermeyer et al. Reunion olivine Fo traverses.

The source workflow only fits the diffusion-gradient section between the
unmodified core/outer-core plateau and the rim in contact with melt. This
script keeps those masks explicit, applies the archived fixed rim boundary,
and compares the source-era Chakraborty/Dohmen law with Diffusor's newer
Oeser et al. tracer-derived interdiffusion law. The latter is a sensitivity
case outside its published XFe/high-silica calibration for these Reunion
olivines; it is not a preferred calibration.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from diffusor.coefficients import Conditions, get
from diffusor.constants import SEC_PER_DAY
from diffusor.fitting import DiffusionModel, fit_time
from diffusor.solvers import InitialCondition, dirichlet, zero_flux
from diffusor.thermo import log_fo2_from_delta

ROOT = Path(__file__).resolve().parent
PROFILE_ROOT = ROOT / "profiles"
OUTPUT = ROOT / "reunion_results.json"

# Values transcribed from Sundermeyer et al. (2020), Online Resource 3,
# sheet REU 150915-1. Profile endpoints/windows are explicitly selected below
# from Online Resource 1 using the paper's Fig. 3 section definitions.
CASES = [
    {
        "sample": "150915-1-3",
        "published_days": 49.0,
        "published_minus_days": 12.0,
        "published_plus_days": 8.0,
        "initial_fo": 0.837,
        "boundary_fo": 0.805,
        "angles_to_axes_deg": [21.7, 87.0, 68.5],
        "dx_um": 3.1,
        "temperature_c": 1139.839,
        "temperature_basis": "mean of the three listed 150915-1 olivine-hosted MI temperatures; this crystal has no same-ID MI analysis",
        "profile_file": "sundermeyer2020_reunion_150915_1_3.csv",
        "mask_min_source_distance_um": 211.18172,
        "mask_max_source_distance_um": 264.18172,
        "mask_note": "retain the right-hand outer-core-to-margin event only; start when Fo first departs the tabulated outer-core Fo by >0.2; omit final low-Fo late-rim point at 268.18 um. This is an explicit reconstruction of the selected gradient, not a source-tabulated mask.",
        "mask_classification": "approximate event-window reconstruction",
    },
    {
        "sample": "150915-1-9",
        "published_days": 579.0,
        "published_minus_days": 72.0,
        "published_plus_days": 104.0,
        "initial_fo": 0.8475,
        "boundary_fo": 0.819,
        "angles_to_axes_deg": [55.1, 37.3, 78.6],
        "dx_um": 1.98288,
        "temperature_c": 1139.576,
        "temperature_basis": "matching 150915-1-9 MI analysis in Online Resource 2",
        "profile_file": "sundermeyer2020_reunion_150915_1_9.csv",
        "mask_min_source_distance_um": 227.882614,
        "mask_max_source_distance_um": 319.372774,
        "mask_note": "right-hand core-to-margin gradient only; start at the first point >0.2 Fo from the tabulated core composition when approaching from the core plateau; end at the archived boundary-composition point. The left-hand side of this transect has a separate rim/gradient and is excluded.",
        "mask_classification": "approximate event-window reconstruction",
    },
]


def model_for(case: dict, law_key: str, domain_extension_um: float = 20.0):
    t_k = case["temperature_c"] + 273.15
    p_pa = 1.0e8  # 1 kbar, explicitly stated in the article's modelling section.
    log_fo2_bar = log_fo2_from_delta("NNO", -0.5, t_k, p_pa)
    source = pd.read_csv(PROFILE_ROOT / case["profile_file"])
    mask = source["Distance_um"].between(
        case["mask_min_source_distance_um"], case["mask_max_source_distance_um"]
    )
    measured = source.loc[mask, ["Distance_um", "Fo_mol"]].copy()
    if measured.empty:
        raise ValueError(f"Empty selected source window for {case['sample']}")
    # The source spreadsheet's distance increases from core towards the
    # tabulated melt-equilibrated margin. Diffusor's left boundary is set to
    # that margin; reverse distances so x=0 is the fixed rim and x increases
    # into the crystal.
    x_boundary = float(measured["Distance_um"].max())
    measured["x_um"] = x_boundary - measured["Distance_um"]
    measured = measured.sort_values("x_um")
    x = measured["x_um"].to_numpy(float)
    y = measured["Fo_mol"].to_numpy(float)
    if x[0] != 0.0 or len(x) < 4:
        raise ValueError(f"Expected a rim point and >=4 points for {case['sample']}")

    # Initial crystal interior is an unzoned core/outer-core plateau; the melt
    # imposes a fixed rim composition from t=0. The far interior is closed.
    core = case["initial_fo"] * 100.0
    rim = case["boundary_fo"] * 100.0
    domain_end = max(float(x.max()) + domain_extension_um, 60.0)
    # Hold the base spatial resolution when extending the closed far-field.
    base_domain_end = max(float(x.max()) + 20.0, 60.0)
    nodes = int(np.ceil(400 * domain_end / base_domain_end)) + 1
    grid = np.linspace(0.0, domain_end, nodes)
    initial = InitialCondition("table", {"x_table": [0.0, domain_end], "C_table": [core, core]},
                               "homogeneous initial outer-core plateau; fixed melt-equilibrium rim")
    conditions = Conditions(
        T_K=t_k,
        P_Pa=p_pa,
        log_fo2_bar=log_fo2_bar,
        X={"XFe": 1.0 - case["initial_fo"]},
        angles_deg=tuple(case["angles_to_axes_deg"]),
    )
    model = DiffusionModel(
        coefficient=get(law_key),
        conditions=conditions,
        initial=initial,
        bc_left=dirichlet(rim),
        bc_right=zero_flux(),
        boundaries_far=False,
        n_nodes=grid.size,
        x_grid=grid,
        composition_dependent=True,
        comp_key="XFe",
        comp_scale=-0.01,
        comp_offset=1.0,
    )
    return model, x, y, log_fo2_bar, measured


def run_case(case: dict) -> dict:
    result = {
        **case,
        "source": "Sundermeyer et al. (2020), doi:10.1007/s00410-019-1642-y",
        "source_setup": "Online Resource 3, REU 150915-1; measured Fo profile Online Resource 1, sample sheet REU 150915-1",
        "conditions": {
            "pressure_kbar": 1.0,
            "fo2_buffer": "NNO-0.5",
            "fo2_log10_bar": None,
            "fo2_log10_pa": None,
        },
        "mask_points": None,
        "fits": {},
    }
    for label, law_key in (
        ("source_law", "ol_FeMg_dohmen_chakraborty2007_tamed"),
        ("alternative_law", "ol_FeMg_oeser2026"),
    ):
        model, x, y, log_fo2_bar, measured = model_for(case, law_key)
        fit = fit_time(
            model,
            x,
            y,
            sigma=np.full(y.size, 0.2),  # source states +/-0.2 mol% Fo
            free_parameters=("t",),
            t_min=0.01 * SEC_PER_DAY,
            t_max=1.0e5 * SEC_PER_DAY,
            scan_points=42,
        )
        result["conditions"]["fo2_log10_bar"] = float(log_fo2_bar)
        result["conditions"]["fo2_log10_pa"] = float(log_fo2_bar + 5.0)
        if result["mask_points"] is None:
            result["mask_points"] = int(len(y))
        model_at_published = model.profile(case["published_days"] * SEC_PER_DAY, x)
        result["fits"][label] = {
            "coefficient": law_key,
            "fit_days": float(fit.t_seconds / SEC_PER_DAY),
            "fit_success": bool(fit.success),
            "fit_message": fit.message,
            "fit_rmse_fo_mol_pct": float(np.sqrt(np.mean((fit.C_model - y) ** 2))),
            "published_time_rmse_fo_mol_pct": float(np.sqrt(np.mean((model_at_published - y) ** 2))),
            "D_m2_s_at_initial_fo": model.D_bulk(C_ref=case["initial_fo"] * 100.0),
            "warnings": model.warnings(fit.t_seconds),
        }
    result["fits"]["source_law"]["fit_to_published_ratio"] = (
        result["fits"]["source_law"]["fit_days"] / case["published_days"]
    )
    result["fits"]["alternative_law"]["alternative_to_source_time_ratio"] = (
        result["fits"]["alternative_law"]["fit_days"] / result["fits"]["source_law"]["fit_days"]
    )
    if case["sample"] == "150915-1-9":
        model, x, y, _, _ = model_for(case, "ol_FeMg_dohmen_chakraborty2007_tamed", domain_extension_um=120.0)
        extended_fit = fit_time(
            model, x, y, sigma=np.full(y.size, 0.2), free_parameters=("t",),
            t_min=0.01 * SEC_PER_DAY, t_max=1.0e5 * SEC_PER_DAY, scan_points=42,
        )
        result["domain_sensitivity"] = {
            "law": "ol_FeMg_dohmen_chakraborty2007_tamed",
            "base_extension_um": 20.0,
            "extended_extension_um": 120.0,
            "base_nodes": 401,
            "extended_nodes": int(model.n_nodes),
            "base_fit_days": result["fits"]["source_law"]["fit_days"],
            "extended_fit_days": float(extended_fit.t_seconds / SEC_PER_DAY),
            "relative_change": float(extended_fit.t_seconds / SEC_PER_DAY / result["fits"]["source_law"]["fit_days"] - 1.0),
        }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    results = [run_case(case) for case in CASES]
    args.output.write_text(json.dumps(results, indent=2, ensure_ascii=False) + "\n", encoding="utf8")
    for r in results:
        print(
            f"{r['sample']}: published {r['published_days']:g} d; "
            f"source law {r['fits']['source_law']['fit_days']:.4g} d; "
            f"alternative {r['fits']['alternative_law']['fit_days']:.4g} d"
        )


if __name__ == "__main__":
    main()
