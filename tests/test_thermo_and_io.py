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
    import warnings
    from diffusor.thermo import BufferRangeWarning
    for T in (1073.15, 1273.15, 1473.15):
        for b in ("NNO", "FMQ"):
            a = log_fo2_buffer(b, T, parameterisation="frost1991")
            with warnings.catch_warnings():     # 1473.15 K is above O'Neill's printed 1420 K for FMQ
                warnings.simplefilter("ignore", BufferRangeWarning)
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


# Frost (1991) Table 1 rows (A, B, C), read from the rendered page
_FROST_TABLE_1 = {
    "aQIF": (-29435.7, 7.391, 0.044), "bQIF": (-29520.8, 7.492, 0.050),
    "IW": (-27489.0, 6.702, 0.055), "WM": (-32807.0, 13.012, 0.083),
    "FMaQ": (-26455.3, 10.344, 0.092), "FMbQ": (-25096.3, 8.735, 0.110),
    "NiNiO": (-24930.0, 9.36, 0.046),
    "MH300": (-25497.5, 14.330, 0.019), "MH573": (-26452.6, 15.455, 0.019),
    "MH682": (-25700.6, 14.558, 0.019),
}


def _table(row, T, P_bar):
    A, B, C = _FROST_TABLE_1[row]
    return A / T + B + C * (P_bar - 1.0) / T


@pytest.mark.parametrize("buffer,row,T_C", [
    ("IW", "IW", 900.0), ("WM", "WM", 900.0), ("NNO", "NiNiO", 900.0),
    ("FMQ", "FMbQ", 900.0), ("FMQ", "FMaQ", 500.0),
    ("QIF", "bQIF", 900.0), ("QIF", "aQIF", 400.0),
    ("HM", "MH682", 900.0), ("HM", "MH573", 600.0), ("HM", "MH300", 400.0)])
def test_frost_rows_are_selected_by_temperature(buffer, row, T_C):
    T = T_C + 273.15
    assert log_fo2_buffer(buffer, T) == pytest.approx(_table(row, T, 1.0), abs=1e-12)


def test_quartz_row_follows_the_pressure_shifted_transition():
    """T(C) = 573 + 0.025 P(bar): 1 kbar moves the transition to 598 C."""
    P = 1000.0 * 1.0e5
    assert log_fo2_buffer("FMQ", 590.0 + 273.15, P) == pytest.approx(_table("FMaQ", 863.15, 1000.0), abs=1e-12)
    assert log_fo2_buffer("FMQ", 605.0 + 273.15, P) == pytest.approx(_table("FMbQ", 878.15, 1000.0), abs=1e-12)


def test_iqf_is_the_quartz_iron_fayalite_buffer_not_iw():
    T = 1273.15
    assert log_fo2_buffer("IQF", T) == log_fo2_buffer("QIF", T)
    assert log_fo2_buffer("IQF", T) == pytest.approx(_table("bQIF", T, 1.0), abs=1e-12)
    assert log_fo2_buffer("IW", T) - log_fo2_buffer("IQF", T) == pytest.approx(0.806, abs=0.001)


def test_buffers_inside_the_printed_range_do_not_warn():
    import warnings
    from diffusor.thermo import BufferRangeWarning
    with warnings.catch_warnings():
        warnings.simplefilter("error", BufferRangeWarning)
        for b, (lo, hi) in {"IW": (565, 1200), "WM": (565, 1200), "FMQ": (400, 1200), "NNO": (600, 1200),
                            "HM": (300, 1100), "QIF": (150, 1200)}.items():
            for T_C in (lo, 0.5 * (lo + hi), hi):
                log_fo2_buffer(b, T_C + 273.15)
        log_fo2_buffer("FMQ", 1200.0, 1.0e5, "oneill")
        log_fo2_buffer("NNO", 1000.0, 1.0e5, "oneill")


@pytest.mark.parametrize("buffer,T_C,parameterisation", [
    ("IW", 500.0, "frost1991"), ("WM", 1250.0, "frost1991"), ("NNO", 550.0, "frost1991"),
    ("FMQ", 350.0, "frost1991"), ("HM", 1150.0, "frost1991"), ("QIF", 100.0, "frost1991"),
    ("FMQ", 1150.0, "oneill"), ("FMQ", 600.0, "oneill"), ("NNO", 1500.0, "oneill")])
def test_extrapolation_beyond_the_printed_range_warns_but_still_returns(buffer, T_C, parameterisation):
    from diffusor.thermo import BufferRangeWarning, buffer_range_status
    T = T_C + 273.15
    inside, message = buffer_range_status(buffer, T, parameterisation)
    assert not inside and "printed for" in message
    with pytest.warns(BufferRangeWarning, match="printed for"):
        value = log_fo2_buffer(buffer, T, 1.0e5, parameterisation)
    assert np.isfinite(value)


def test_buffer_accepts_temperature_arrays():
    T = np.array([873.15, 1073.15, 1273.15])
    vec = log_fo2_buffer("FMQ", T)
    assert vec.shape == T.shape
    assert [float(v) for v in vec] == pytest.approx([log_fo2_buffer("FMQ", float(t)) for t in T])


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
    assert "buffer at every sampled temperature" in text
    assert "References" in text
    assert "https://doi.org/10.2138/am-2016-5815" in text

    written = save_results(tmp_path, r, mc, basename="t")
    assert set(written) >= {"json", "profile", "methods", "montecarlo"}
    import json
    d = json.loads((tmp_path / "t_results.json").read_text(encoding="utf-8"))
    assert d["coefficient"]["key"] == "opx_FeMg_dohmen2016"
    assert d["fit"]["t_seconds"] > 0
    assert "citations" in d and len(d["citations"]) > 5


# --- propagation of the oxide uncertainties through the cation-mole conversion ------
def _fd_sigma(a, b, sa, sb, mode, oa, ob, h=1e-6):
    from diffusor.thermo.units import composition_variable

    def f(x, y):
        return composition_variable(np.array([x]), np.array([y]), mode, oa, ob)[0]
    da = (f(a + h, b) - f(a - h, b)) / (2 * h)
    db = (f(a, b + h) - f(a, b - h)) / (2 * h)
    return float(np.hypot(da * sa, db * sb))


@pytest.mark.parametrize("a,b,sa,sb,oa,ob", [
    (10.0, 10.0, 0.2, 0.2, "FeO", "MgO"),           # one cation each
    (10.0, 5.0, 0.1, 0.05, "CaO", "Na2O"),          # Na2O carries two cations
    (8.0, 3.0, 0.05, 0.2, "Na2O", "K2O"),           # two cations each
    (20.0, 15.0, 0.2, 0.15, None, None),            # raw columns, no oxide conversion
])
@pytest.mark.parametrize("mode", ["A/(A+B)", "B/(A+B)", "A-B"])
def test_composition_sigma_matches_finite_differences(a, b, sa, sb, oa, ob, mode):
    from diffusor.dataio.profiles import composition_sigma
    s = composition_sigma(np.array([a]), np.array([b]), np.array([sa]), np.array([sb]), mode, oa, ob)[0]
    assert s == pytest.approx(_fd_sigma(a, b, sa, sb, mode, oa, ob), rel=1e-6)


def test_composition_sigma_matches_monte_carlo():
    from diffusor.dataio.profiles import composition_sigma
    from diffusor.thermo.units import composition_variable
    rng = np.random.default_rng(7)
    for a, b, sa, sb, oa, ob in [(10.0, 10.0, 0.2, 0.2, "FeO", "MgO"),
                                 (10.0, 5.0, 0.1, 0.05, "CaO", "Na2O")]:
        draws = composition_variable(rng.normal(a, sa, 400000), rng.normal(b, sb, 400000),
                                     "A/(A+B)", oa, ob)
        s = composition_sigma(np.array([a]), np.array([b]), np.array([sa]), np.array([sb]),
                              "A/(A+B)", oa, ob)[0]
        assert s == pytest.approx(np.std(draws), rel=0.01)


def test_oxide_profile_sigma_is_propagated_through_the_cation_moles():
    """FeO = MgO = 10 wt%, 0.2 wt% each: the molar ratio has sigma 6.51e-3. The old
    propagation from the raw wt% ratio gave 7.07e-3 (8.6 % too high)."""
    df = pd.DataFrame({"x": [0.0, 1.0, 2.0], "FeO": [10.0] * 3, "MgO": [10.0] * 3,
                       "eA": [0.2] * 3, "eB": [0.2] * 3})
    spec = ProfileSpec("x", "FeO", "MgO", "eA", "eB", mode="A/(A+B)", oxide_a="FeO", oxide_b="MgO")
    p = build_profile(df, spec)
    assert p.sigma[0] == pytest.approx(6.51e-3, rel=2e-3)
    assert p.sigma[0] < 0.93 * 7.07e-3
    # without oxides the raw-column formula is unchanged
    spec_raw = ProfileSpec("x", "FeO", "MgO", "eA", "eB", mode="A/(A+B)")
    assert build_profile(df, spec_raw).sigma[0] == pytest.approx(np.sqrt(2) * 10 * 0.2 / 400, rel=1e-9)


def test_exactly_determined_greyscale_calibration_does_not_claim_zero_uncertainty():
    cal = calibrate([50.0, 200.0], [0.2, 0.5], degree=1)
    assert cal.exactly_determined and np.all(np.isnan(cal.sigma_from_calibration([100.0, 150.0])))
    x = np.array([0.0, 1.0])
    with pytest.warns(UserWarning, match="exactly determined"):
        _, C, sig = apply_calibration(x, np.array([100.0, 150.0]), cal)
    assert np.all(np.isnan(sig))                          # unknown, not 0
    with pytest.warns(UserWarning, match="exactly determined"):
        _, C, sig = apply_calibration(x, np.array([100.0, 150.0]), cal, grey_scatter=np.array([3.0, 3.0]))
    assert np.allclose(sig, (0.3 / 150.0) * 3.0)          # only the grey-value scatter, with the slope
    over = calibrate([50.0, 100.0, 200.0], [0.2, 0.31, 0.5], degree=1)
    assert not over.exactly_determined and np.all(over.sigma_from_calibration([120.0]) > 0)
