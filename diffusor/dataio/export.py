"""Exporting results: data, figures and a fully referenced methods paragraph.

The methods block is the point of the whole citation machinery: every run
writes out the equations, coefficients, conventions and software versions it
actually used, with a reference list that is generated from the citation keys
touched by that run, so a result can always be traced back to its sources.
"""
from __future__ import annotations

import json
import platform
from dataclasses import asdict, is_dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

from .. import __version__
from ..references import format_reference, get as get_reference
from ..thermo.units import human_time


def _fmt(v):
    if isinstance(v, (np.floating, float)):
        return float(v)
    if isinstance(v, (np.integer, int)):
        return int(v)
    if isinstance(v, np.ndarray):
        return v.tolist()
    return v


def collect_citations(fit_result, mc_result=None, extra: Sequence[str] = ()) -> List[str]:
    """Every citation key used by this run, in a stable order."""
    keys: List[str] = ["crank1975", "costa2008"]
    model = fit_result.model
    c = model.coefficient
    keys.append(c.citation)
    keys.extend(c.secondary_citations)
    if c.axis_factors and (model.conditions.angles_deg or model.conditions.axis):
        keys.append("costa_chakraborty2004")
    if model.activity_theta:
        keys.extend(["costa2003", "dohmen2017", "dohmen_blundy2014", "grove1984"])
    if model.beam_sigma_um:
        keys.extend(["ganguly1988", "bradshaw_kent2017"])
    if model.history is not None and not model.history.is_isothermal:
        keys.append("lasaga1983")
    ok, _ = model.can_use_analytical()
    if not ok:
        keys.append("dohmen2017")
    if mc_result is not None and mc_result.budget.fo2_mode == "buffer":
        keys.append("frost1991")
    keys.extend(["shea2015", "krimer_costa2017", "codata2018", "iupac2021"])
    keys.extend(extra)
    seen, out = set(), []
    for k in keys:
        if k not in seen:
            seen.add(k)
            out.append(k)
    return out


def methods_paragraph(fit_result, mc_result=None, profile=None) -> str:
    """A methods block describing exactly what this run did, with references."""
    model = fit_result.model
    c = model.coefficient
    cond = model.conditions
    ok, why = model.can_use_analytical()
    L: List[str] = []

    L.append("METHODS")
    L.append("=" * 70)
    L.append(f"Diffusion chronometry performed with Diffusor v{__version__} "
             f"(Python {platform.python_version()}) on "
             f"{datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}.")
    L.append("")

    L.append("Data")
    L.append("-" * 70)
    if profile is not None:
        L.append(f"Source file: {profile.source}")
        if profile.spec:
            L.append(f"Column mapping: {profile.spec.describe()}")
        L.append(f"{len(profile)} points spanning "
                 f"{profile.x.min():.2f} to {profile.x.max():.2f} um.")
        for n in profile.notes:
            L.append(f"  Note: {n}")
    else:
        L.append(f"{fit_result.x_data.size} points spanning "
                 f"{fit_result.x_data.min():.2f} to {fit_result.x_data.max():.2f} um.")
    L.append("")

    L.append("Diffusion coefficient")
    L.append("-" * 70)
    L.append(c.describe())
    L.append("")

    L.append("Conditions")
    L.append("-" * 70)
    L.append(f"T = {cond.T_K - 273.15:.1f} C ({cond.T_K:.2f} K), P = {cond.P_Pa/1e6:.1f} MPa")
    if cond.log_fo2_bar is not None:
        L.append(f"log10 fO2 = {cond.log_fo2_bar:.3f} (bar) = {cond.log_fo2_Pa:.3f} (Pa)")
        if mc_result is not None and mc_result.budget.fo2_mode == "buffer":
            b = mc_result.budget
            L.append(f"  set as {b.buffer} {b.delta_buffer:+.2f} log units, with the buffer "
                     f"evaluated from Frost (1991) Table 1. The Monte Carlo re-evaluates the "
                     f"buffer at every sampled temperature, so fO2 and T stay correlated.")
    for k, v in cond.X.items():
        if np.ndim(v) == 0:
            L.append(f"{k} = {float(v):.4f}")
    if cond.axis:
        L.append(f"Traverse along the {cond.axis}-axis.")
    if cond.angles_deg:
        L.append(f"Traverse at {cond.angles_deg} degrees to a, b, c. The direction-cosine "
                 f"relation of Costa & Chakraborty (2004) was applied.")
    if model.history is not None and not model.history.is_isothermal:
        L.append(f"Non-isothermal history ({model.history.label}). The diffusion integral "
                 f"int D(T(t)) dt was evaluated numerically (Crank 1975 eq. 7.7, Lasaga 1983).")
    L.append("")

    L.append("Model")
    L.append("-" * 70)
    L.append(f"Geometry: {model.geometry.kind}. {model.geometry.describe()}")
    L.append(f"Initial condition: {model.initial.describe()}")
    L.append(f"Boundaries: left {model.bc_left.describe()}, right {model.bc_right.describe()}")
    L.append(f"Solver: {'analytical, ' + why if ok else 'numerical Crank-Nicolson (theta = 1/2), ' + why}")
    if not ok:
        L.append(f"  Finite-difference scheme after Crank (1975) section 8.4 with D evaluated at "
                 f"half-nodes (Dohmen et al. 2017, Appendix eqs A17-A21), {model.n_nodes} grid nodes.")
    if model.beam_sigma_um:
        L.append(f"Model profiles were convolved with a Gaussian of sigma = "
                 f"{model.beam_sigma_um:.2f} um to account for the analytical spatial "
                 f"resolution (Ganguly et al. 1988, Bradshaw & Kent 2017).")
    L.append("")

    L.append("Result")
    L.append("-" * 70)
    L.append(f"Best-fit time: {human_time(fit_result.t_seconds)} ({fit_result.t_seconds:.6g} s)")
    L.append(f"Fit quality: {fit_result.stats.describe()}")
    for k, v in fit_result.free.items():
        if k != "t":
            L.append(f"Fitted {k} = {v:.5g}")
    if mc_result is not None:
        b = mc_result.budget
        L.append("")
        L.append(f"Uncertainties were propagated by Monte Carlo with {mc_result.n_draws} draws "
                 f"(seed {mc_result.seed}), re-fitting the time for every draw. Sampled sources: "
                 f"{', '.join(b.active_sources())}.")
        if "temperature" in b.active_sources():
            L.append(f"  T: 1 sigma = {b.sigma_T_K:.1f} K")
        if "fo2" in b.active_sources():
            L.append(f"  fO2: 1 sigma = {b.sigma_delta_buffer:.2f} log units on the buffer offset")
        if "pressure" in b.active_sources():
            L.append(f"  P: 1 sigma = {b.sigma_P_Pa/1e6:.1f} MPa")
        if "diffusion_coefficient" in b.active_sources():
            mode = c.default_sampling_mode() if b.coefficient_mode == "auto" else b.coefficient_mode
            expl = {"covariance": "sampled from the published parameter covariance matrix",
                    "logD_at_T": "ln D sampled at the working temperature from the scatter "
                                 "reported by the source, so ln D0 and Q stay correlated",
                    "independent": "each Arrhenius parameter sampled independently. This "
                                   "ignores the strong ln D0 - Q correlation and overstates "
                                   "the uncertainty"}[mode]
            L.append(f"  Diffusion coefficient: {expl}.")
        L.append(f"Median time {human_time(mc_result.median)}. 68% interval "
                 f"{human_time(mc_result.p16)} to {human_time(mc_result.p84)}, 95% interval "
                 f"{human_time(mc_result.p2_5)} to {human_time(mc_result.p97_5)}. Times are "
                 f"log-normally distributed, so percentiles are reported.")
        if mc_result.contributions:
            L.append("Contribution of each source to sigma(log10 t), one at a time:")
            for k, v in sorted(mc_result.contributions.items(), key=lambda kv: -kv[1]):
                L.append(f"  {k:<24s} {v:.3f}")
    L.append("")

    warns = list(fit_result.warnings) + (list(mc_result.warnings) if mc_result else [])
    if warns:
        L.append("Caveats")
        L.append("-" * 70)
        for w in dict.fromkeys(warns):
            L.append(f"- {w}")
        L.append("")

    L.append("References")
    L.append("-" * 70)
    for k in collect_citations(fit_result, mc_result):
        L.append(f"- {format_reference(k)}")
    return "\n".join(L)


def result_dict(fit_result, mc_result=None, profile=None) -> Dict:
    model = fit_result.model
    c = model.coefficient
    cond = model.conditions
    d = {
        "diffusor_version": __version__,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "coefficient": {
            "key": c.key, "label": c.label, "citation": c.citation,
            "equation_number": c.equation_number, "equation": c.equation_text,
            "verified": c.verified, "verified_from": c.verified_from,
            "parameters": {k: {"value": p.value, "sigma": p.sigma, "unit": p.unit,
                               "sigma_level": p.sigma_level}
                           for k, p in c.params.items()},
        },
        "conditions": {
            "T_K": cond.T_K, "T_C": cond.T_K - 273.15, "P_Pa": cond.P_Pa,
            "log_fo2_bar": cond.log_fo2_bar,
            "X": {k: _fmt(v) for k, v in cond.X.items() if np.ndim(v) == 0},
            "axis": cond.axis, "angles_deg": cond.angles_deg,
        },
        "model": {
            "geometry": model.geometry.kind,
            "initial_condition": {"kind": model.initial.kind,
                                  "params": {k: _fmt(v) for k, v in model.initial.params.items()
                                             if np.ndim(v) == 0}},
            "bc_left": model.bc_left.kind, "bc_right": model.bc_right.kind,
            "beam_sigma_um": model.beam_sigma_um, "n_nodes": model.n_nodes,
            "route": fit_result.route,
        },
        "fit": {
            "t_seconds": fit_result.t_seconds,
            "t_years": fit_result.t_years,
            "t_human": human_time(fit_result.t_seconds),
            "chi2": fit_result.stats.chi2,
            "reduced_chi2": fit_result.stats.reduced_chi2,
            "r_squared": fit_result.stats.r_squared,
            "rmse": fit_result.stats.rmse,
            "free_parameters": {k: _fmt(v) for k, v in fit_result.free.items()},
        },
        "warnings": list(fit_result.warnings),
        "citations": {k: format_reference(k) for k in collect_citations(fit_result, mc_result)},
    }
    if profile is not None and profile.spec is not None:
        d["data"] = {"source": profile.source, "spec": asdict(profile.spec),
                     "n_points": len(profile), "notes": profile.notes}
    if mc_result is not None:
        d["montecarlo"] = {
            "n_draws": mc_result.n_draws, "n_failed": mc_result.n_failed,
            "seed": mc_result.seed,
            "median_s": mc_result.median, "p16_s": mc_result.p16, "p84_s": mc_result.p84,
            "p2.5_s": mc_result.p2_5, "p97.5_s": mc_result.p97_5,
            "sigma_log10": mc_result.sigma_log10,
            "budget": asdict(mc_result.budget),
            "contributions": mc_result.contributions,
        }
    return d


def save_results(directory, fit_result, mc_result=None, profile=None,
                 figure=None, basename: str = "diffusor_run") -> Dict[str, str]:
    """Write JSON, the profile CSV, the methods block and (optionally) the figure."""
    out = Path(directory)
    out.mkdir(parents=True, exist_ok=True)
    written: Dict[str, str] = {}

    j = out / f"{basename}_results.json"
    j.write_text(json.dumps(result_dict(fit_result, mc_result, profile), indent=2,
                            default=_fmt), encoding="utf-8")
    written["json"] = str(j)

    import pandas as pd
    df = pd.DataFrame({"x_um": fit_result.x_data,
                       "C_measured": fit_result.C_data,
                       "C_model": fit_result.C_model,
                       "residual": fit_result.C_data - fit_result.C_model})
    if fit_result.sigma is not None:
        df["sigma"] = fit_result.sigma
    try:
        df["C_initial"] = fit_result.model.initial.evaluate(fit_result.x_data)
    except Exception:
        pass
    if mc_result is not None and mc_result.profiles is not None:
        lo, hi = mc_result.envelope()
        df["C_model_p16"], df["C_model_p84"] = lo, hi
    p = out / f"{basename}_profile.csv"
    df.to_csv(p, index=False)
    written["profile"] = str(p)

    m = out / f"{basename}_methods.txt"
    m.write_text(methods_paragraph(fit_result, mc_result, profile), encoding="utf-8")
    written["methods"] = str(m)

    if mc_result is not None:
        t = out / f"{basename}_montecarlo_times.csv"
        pd.DataFrame({"t_seconds": mc_result.times,
                      "t_years": mc_result.times / 3.15576e7}).to_csv(t, index=False)
        written["montecarlo"] = str(t)

    if figure is not None:
        f = out / f"{basename}_figure.png"
        figure.savefig(f, dpi=300, bbox_inches="tight")
        written["figure"] = str(f)
        fs = out / f"{basename}_figure.svg"
        figure.savefig(fs, bbox_inches="tight")
        written["figure_svg"] = str(fs)
    return written
