"""Buffers, unit conversions, profile loading, greyscale calibration and export."""
import numpy as np
import pandas as pd
import pytest

from diffusor.constants import OXIDE_MOLAR_MASS, R_GAS
from diffusor.dataio import (ProfileSpec, build_profile, calibrate, apply_calibration,
                             anchors_from_microprobe, methods_paragraph, save_results,
                             suggest_spec)
from diffusor.references import REFERENCES, cite, format_reference, to_bibtex
from diffusor.thermo import (available_buffers, delta_from_log_fo2, log_fo2_buffer,
                             log_fo2_from_delta, x_an_from_oxides, x_fe_from_oxides)
from diffusor.thermo.units import human_time, pressure_to_pa, time_to_seconds


# --- buffers -------------------------------------------------------------------
def test_buffer_ordering_is_physically_correct():
    """IW is the most reducing, HM the most oxidising (Frost 1991 Fig. 1)."""
    T = 1273.15
    vals = {b: log_fo2_buffer(b, T) for b in ("IW", "WM", "FMQ", "NNO", "HM")}
    assert vals["IW"] < vals["WM"] < vals["FMQ"] < vals["NNO"] < vals["HM"]


def test_nno_and_fmq_agree_between_the_two_parameterisations():
    for T in (1073.15, 1273.15, 1473.15):
        for b in ("NNO", "FMQ"):
            a = log_fo2_buffer(b, T, parameterisation="frost1991")
            c = log_fo2_buffer(b, T, parameterisation="oneill")
            assert abs(a - c) < 0.25, f"{b} at {T} K: {a:.3f} vs {c:.3f}"


def test_nno_is_about_minus_ten_at_1000C():
    """Widely quoted value: log fO2 at NNO, 1000 C, 1 bar is about -10.2."""
    assert log_fo2_buffer("NNO", 1273.15) == pytest.approx(-10.2, abs=0.3)


def test_buffer_offsets_round_trip():
    T = 1223.15
    lf = log_fo2_from_delta("NNO", 1.3, T)
    assert delta_from_log_fo2("NNO", lf, T) == pytest.approx(1.3)


def test_qfm_is_an_alias_of_fmq():
    assert log_fo2_buffer("QFM", 1273.15) == log_fo2_buffer("FMQ", 1273.15)


def test_pressure_raises_the_buffer_fo2():
    assert log_fo2_buffer("NNO", 1273.15, 5.0e8) > log_fo2_buffer("NNO", 1273.15, 1.0e5)


def test_unknown_buffer_is_rejected():
    with pytest.raises(ValueError):
        log_fo2_buffer("XYZ", 1273.15)


# --- units and compositions ------------------------------------------------------
def test_oxide_molar_masses():
    assert OXIDE_MOLAR_MASS["SiO2"] == pytest.approx(60.08, abs=0.02)
    assert OXIDE_MOLAR_MASS["FeO"] == pytest.approx(71.84, abs=0.02)
    assert OXIDE_MOLAR_MASS["MgO"] == pytest.approx(40.30, abs=0.02)


def test_x_fe_from_oxides_is_molar_not_by_weight():
    xfe = float(x_fe_from_oxides(10.0, 30.0))
    fe = 10.0 / OXIDE_MOLAR_MASS["FeO"]
    mg = 30.0 / OXIDE_MOLAR_MASS["MgO"]
    assert xfe == pytest.approx(fe / (fe + mg))
    assert xfe != pytest.approx(10.0 / 40.0)


def test_x_an_from_oxides():
    assert float(x_an_from_oxides(20.0, 0.0)) == pytest.approx(1.0)
    assert 0.0 < float(x_an_from_oxides(12.0, 4.0)) < 1.0


def test_time_conversions_and_human_time():
    assert time_to_seconds(1.0, "yr") == pytest.approx(31557600.0)
    assert pressure_to_pa(1.0, "GPa") == 1e9
    assert "yr" in human_time(3.2e7)
    assert "d" in human_time(2.0e5)


# --- references ------------------------------------------------------------------
def test_every_reference_renders_and_has_a_bibtex_entry():
    for key, ref in REFERENCES.items():
        assert ref.full()
        assert to_bibtex(ref).startswith("@")
        assert cite(key)


def test_cite_formats_author_counts():
    assert cite("crank1975") == "Crank (1975)"
    assert "et al." in cite("costa2008")
    assert "&" in cite("freer_hauptman1978")


def test_unknown_citation_key_raises():
    with pytest.raises(KeyError):
        format_reference("not_a_real_key")


# --- profile loading ---------------------------------------------------------------
def _frame():
    return pd.DataFrame({
        "Distance_um": np.linspace(0, 40, 21),
        "FeO_wt": np.linspace(12.0, 18.0, 21),
        "MgO_wt": np.linspace(28.0, 22.0, 21),
        "FeO_err": np.full(21, 0.15),
        "MgO_err": np.full(21, 0.20),
    })


def test_suggest_spec_finds_the_obvious_columns():
    spec = suggest_spec(_frame())
    assert spec.distance_column == "Distance_um"
    assert spec.column_a == "FeO_wt" and spec.column_b == "MgO_wt"
    assert spec.sigma_a_column == "FeO_err"


def test_build_profile_makes_a_molar_ratio_and_propagates_sigma():
    df = _frame()
    spec = ProfileSpec("Distance_um", "FeO_wt", "MgO_wt", "FeO_err", "MgO_err",
                       mode="A/(A+B)", oxide_a="FeO", oxide_b="MgO")
    p = build_profile(df, spec, "memory")
    assert len(p) == 21
    assert np.all((p.C > 0) & (p.C < 1))
    assert p.C[0] < p.C[-1]                       # Fe increases along the traverse
    assert p.sigma is not None and np.all(p.sigma > 0)
    assert np.all(p.sigma < 0.02)


def test_build_profile_drops_missing_rows():
    df = _frame()
    df.loc[3, "FeO_wt"] = np.nan
    spec = ProfileSpec("Distance_um", "FeO_wt", "MgO_wt", mode="A/(A+B)")
    p = build_profile(df, spec)
    assert len(p) == 20
    assert any("dropped" in n for n in p.notes)


def test_two_sigma_column_is_halved():
    df = _frame()
    spec1 = ProfileSpec("Distance_um", "FeO_wt", "MgO_wt", "FeO_err", "MgO_err",
                        mode="A/(A+B)", sigma_level="1s")
    spec2 = ProfileSpec("Distance_um", "FeO_wt", "MgO_wt", "FeO_err", "MgO_err",
                        mode="A/(A+B)", sigma_level="2s")
    assert np.allclose(build_profile(df, spec2).sigma, build_profile(df, spec1).sigma / 2.0)


def test_single_column_mode_passes_values_through():
    df = _frame()
    p = build_profile(df, ProfileSpec("Distance_um", "FeO_wt", mode="A"))
    assert np.allclose(p.C, df["FeO_wt"].to_numpy())


# --- greyscale ----------------------------------------------------------------------
def test_linear_greyscale_calibration_is_exact_for_linear_data():
    g = np.array([100.0, 120.0, 140.0, 160.0])
    c = 0.5 - 0.002 * g
    cal = calibrate(g, c, degree=1)
    assert cal.rmse < 1e-12
    assert cal.r_squared == pytest.approx(1.0)
    assert np.allclose(cal.apply(g), c)


def test_greyscale_calibration_flags_a_poor_linear_fit():
    g = np.linspace(100, 200, 8)
    c = 0.3 + 1e-5 * (g - 150) ** 2
    cal = calibrate(g, c, degree=1)
    assert any("non-linear" in n for n in cal.notes)


def test_greyscale_needs_enough_anchors():
    with pytest.raises(ValueError):
        calibrate([100.0], [0.3], degree=1)


def test_apply_calibration_propagates_grey_scatter():
    g = np.linspace(100, 180, 6)
    c = 0.5 - 0.002 * g
    cal = calibrate(g, c, degree=1)
    x = np.linspace(0, 50, 6)
    _, C, sig = apply_calibration(x, g, cal, grey_scatter=np.full(6, 2.0))
    assert np.all(sig >= 0.002 * 2.0 * 0.99)


def test_anchors_from_microprobe_averages_the_window():
    xg = np.linspace(0, 100, 501)
    grey = 100 + xg
    g, c = anchors_from_microprobe(xg, grey, [10.0, 50.0], [0.2, 0.4], window_um=1.0)
    assert g[0] == pytest.approx(110.0, abs=0.5)
    assert list(c) == [0.2, 0.4]


# --- export -------------------------------------------------------------------------
def test_methods_paragraph_lists_the_sources_actually_used(tmp_path):
    from diffusor.coefficients import Conditions, get
    from diffusor.constants import SEC_PER_YEAR
    from diffusor.fitting import DiffusionModel, fit_time, run_montecarlo, UncertaintyBudget
    from diffusor.solvers import Geometry, InitialCondition, dirichlet

    coef = get("opx_FeMg_dohmen2016")
    cond = Conditions(T_K=1223.15, P_Pa=2e8, log_fo2_bar=-10.0, X={"XFe": 0.2}, axis="c")
    ic = InitialCondition("step", {"x0": 0.0, "C_left": 0.30, "C_right": 0.18})
    m = DiffusionModel(coefficient=coef, conditions=cond, initial=ic, geometry=Geometry("plane"),
                       bc_left=dirichlet(0.30), bc_right=dirichlet(0.18), n_nodes=151)
    x = np.linspace(-30, 30, 31)
    C = m.profile(SEC_PER_YEAR, x)
    r = fit_time(m, x, C)
    mc = run_montecarlo(m, x, C, np.full_like(x, 0.003), n_draws=20, seed=1,
                        budget=UncertaintyBudget(sigma_T_K=20.0, buffer="NNO",
                                                 delta_buffer=1.0, sigma_delta_buffer=0.3))
    text = methods_paragraph(r, mc)
    assert "Dohmen" in text and "Crank" in text
    assert "re-evaluated at every sampled temperature" in text
    assert "References" in text
    assert "https://doi.org/10.2138/am-2016-5815" in text

    written = save_results(tmp_path, r, mc, basename="t")
    assert set(written) >= {"json", "profile", "methods", "montecarlo"}
    import json
    d = json.loads((tmp_path / "t_results.json").read_text(encoding="utf-8"))
    assert d["coefficient"]["key"] == "opx_FeMg_dohmen2016"
    assert d["fit"]["t_seconds"] > 0
    assert "citations" in d and len(d["citations"]) > 5
