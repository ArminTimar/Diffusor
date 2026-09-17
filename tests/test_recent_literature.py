"""The 2020-2026 calibrations: each entry must reproduce numbers printed in its paper."""
import numpy as np
import pytest

from diffusor.coefficients import Conditions, get, list_coefficients
from diffusor.coefficients.magnetite import SIEVWRIGHT_TABLE5, TABLE12_PURE, sievwright_D_1150
from diffusor.coefficients.opx import dias2025_regime
from diffusor.constants import R_GAS
from diffusor.dataio import ProfileSpec, build_profile, read_table
from diffusor.thermo.buffers import log_fo2_buffer

LN10 = np.log(10.0)


def logD(key, T_C, **kw):
    X = kw.pop("X", {})
    return float(np.log10(get(key).D(Conditions(T_K=T_C + 273.15, X=X, **kw))))


# --- plagioclase Sr and Ba: Grocolas et al. (2025) ---------------------------------
def test_grocolas_sr_equation_7_by_hand():
    T = 1000 + 273.15
    expected = -1.65 * 0.5 - 3.03 - 368142.0 / (LN10 * R_GAS * T)
    assert logD("plag_Sr_grocolas2025", 1000, X={"XAn": 0.5}) == pytest.approx(expected, abs=1e-9)


def test_grocolas_sr_matches_its_own_time_series_means():
    # section 3.2: oligoclase An28 at 1000 C gave log D = -18.23 +/- 0.24 and labradorite
    # An67 at 1100 C gave -18.61 +/- 0.26; the global fit (eq. 7) must sit within about
    # two standard deviations of both
    assert logD("plag_Sr_grocolas2025", 1000, X={"XAn": 0.28}) == pytest.approx(-18.23, abs=0.5)
    assert logD("plag_Sr_grocolas2025", 1100, X={"XAn": 0.67}) == pytest.approx(-18.61, abs=0.55)


def test_new_sr_is_one_to_two_orders_slower_and_the_gap_widens_as_t_falls():
    # Grocolas et al. (2025) section 4.1: "differ by ~1-2 log units"; the higher activation
    # energy makes the gap grow towards magmatic temperatures (2.8 log units at 750 C and
    # An36 according to Audetat et al. 2026)
    gaps = []
    for T_C in (1100, 1000, 900, 750):
        gaps.append(logD("plag_Sr_giletti_casserly1994", T_C, X={"XAn": 0.36})
                    - logD("plag_Sr_grocolas2025", T_C, X={"XAn": 0.36}))
    assert 1.3 < gaps[0] < 1.9 and 1.9 < gaps[2] < 2.5
    assert np.all(np.diff(gaps) > 0)
    assert gaps[3] == pytest.approx(2.8, abs=0.3)


def test_short_course_giletti_casserly_agrees_with_the_grocolas_refit_eq_12():
    for T_C in (800, 950, 1100):
        for XAn in (0.2, 0.5, 0.9):
            T = T_C + 273.15
            eq12 = -3.76 * XAn - 4.52 - 270607.0 / (LN10 * R_GAS * T)
            assert logD("plag_Sr_giletti_casserly1994", T_C, X={"XAn": XAn}) == pytest.approx(eq12, abs=0.25)


def test_grocolas_ba_is_close_to_cherniak_2002():
    # Grocolas et al.: Ba "similar to slightly faster (~0.5 log unit)" than Cherniak (2002)
    for T_C, XAn in [(950, 0.28), (1100, 0.67)]:
        diff = (logD("plag_Ba_grocolas2025", T_C, X={"XAn": XAn})
                - logD("plag_Ba_cherniak2002", T_C, X={"XAn": XAn}))
        assert -0.3 < diff < 0.9


def test_grocolas_covariance_follows_their_perfect_compensation_assumption():
    c = get("plag_Sr_grocolas2025")
    assert c.default_sampling_mode() == "covariance"
    rng = np.random.default_rng(1)
    draws = np.array([[d["a"], d["b"], d["Q"]] for d in (c.sample(rng) for _ in range(4000))])
    assert np.corrcoef(draws[:, 1], draws[:, 2])[0, 1] > 0.99
    assert abs(np.corrcoef(draws[:, 0], draws[:, 1])[0, 1]) < 0.05
    assert np.std(draws[:, 2]) == pytest.approx(27.141, rel=0.05)


def test_recommended_plagioclase_sr_and_ba_are_the_2025_laws():
    rec = {c.species: c.key for c in list_coefficients("plagioclase") if c.recommended}
    assert rec["Sr"] == "plag_Sr_grocolas2025"
    assert rec["Ba"] == "plag_Ba_grocolas2025"


# --- plagioclase Mg and Li ----------------------------------------------------------
def test_audetat_mg_equation_1_and_silica_activity():
    T = 1100 + 273.15
    expected = -2.99 * 0.6 - 4.03 - 262914.0 / (LN10 * R_GAS * T)
    c = get("plag_Mg_audetat2026")
    cond = Conditions(T_K=T, X={"XAn": 0.6})
    assert np.log10(c.D(cond)) == pytest.approx(expected, abs=1e-9)
    assert np.log10(c.D(cond, {"aSiO2": 0.0})) == pytest.approx(expected - 1.87, abs=1e-9)


def test_audetat_mg_is_consistent_with_van_orman_at_silica_saturation():
    for T_C, XAn in [(900, 0.4), (1100, 0.7)]:
        assert logD("plag_Mg_audetat2026", T_C, X={"XAn": XAn}) == pytest.approx(
            logD("plag_Mg_vanorman2014", T_C, X={"XAn": XAn}), abs=0.6)


def test_pohl_li_interstitial_is_faster_by_point_two_to_one_order():
    for T_C in (700, 900, 1100):
        diff = logD("plag_Li_pohl2024_interstitial", T_C) - logD("plag_Li_pohl2024_vacancy", T_C)
        assert 0.1 < diff < 1.1, (T_C, diff)
    assert logD("plag_Li_pohl2024_interstitial", 1000) == pytest.approx(
        -3.76 - 180000.0 / (LN10 * R_GAS * 1273.15), abs=1e-9)


# --- orthopyroxene: Dias & Dohmen (2024), Dias et al. (2025) --------------------------
DIAS2024_TABLE1 = [   # T (C), fitted m, log D0 at XFe = 0 (m2/s), all at log fO2 = -7 Pa
    (1102, 3.7, -18.69), (1050, 3.0, -18.95), (1050, 2.9, -19.11), (950, 1.1, -19.72)]


def test_dias2024_m_of_T_reproduces_the_fitted_exponents():
    c = get("opx_FeMg_dias_dohmen2024")
    for T_C, m_tab, _ in DIAS2024_TABLE1:
        m = c.params["m_slope"].value / (T_C + 273.15) + c.params["m_int"].value
        assert m == pytest.approx(m_tab, abs=0.3)


@pytest.mark.parametrize("key,kw", [
    ("opx_FeMg_dias_dohmen2024", {}),
    ("opx_FeMg_dias2025", {"log_fo2_bar": -12.0}),       # -7 Pa
])
def test_dias_laws_reproduce_the_2024_experiments_at_fs9(key, kw):
    for T_C, m_tab, logD0 in DIAS2024_TABLE1:
        measured = logD0 + m_tab * 0.09
        assert logD(key, T_C, X={"XFe": 0.09}, **kw) == pytest.approx(measured, abs=0.2)


def test_dias2025_equation_22_by_hand_and_regime_switch():
    T = 1050 + 273.15
    m1 = -2.37e4 / T + 21.09
    lf_pa = -8.0
    expected = (np.log10(3.085e-8) + 0.25 * (lf_pa + 7.0)
                - 284000.0 / (LN10 * R_GAS * T) + m1 * (0.3 - 0.1))
    assert logD("opx_FeMg_dias2025", 1050, X={"XFe": 0.3}, log_fo2_bar=lf_pa - 5) == pytest.approx(expected, abs=1e-9)
    m2 = 2.96e4 / T - 21.08
    expected2 = np.log10(1.93e-10) - 246000.0 / (LN10 * R_GAS * T) + m2 * (0.3 - 0.1)
    assert logD("opx_FeMg_dias2025", 1050, X={"XFe": 0.3}, log_fo2_bar=-16.0) == pytest.approx(expected2, abs=1e-9)
    assert "eq. 22" in dias2025_regime(-9.9) and "eq. 23" in dias2025_regime(-10.0)


def test_dohmen2016_is_now_superseded_by_dias2025():
    old = get("opx_FeMg_dohmen2016")
    assert old.superseded_by == "dias2025" and not old.recommended
    assert get("opx_FeMg_dias2025").recommended


def test_ostorero_form_is_ganguly_tazzoli_without_the_fo2_term():
    with_fo2 = get("opx_FeMg_ganguly_tazzoli1994")
    cond = Conditions(T_K=1123.15, X={"XFe": 0.3}, log_fo2_bar=-12.0)
    no = get("opx_FeMg_ganguly_tazzoli1994_nofo2").D(cond)
    assert no == pytest.approx(with_fo2.D(cond, {"use_fo2": 0.0}), rel=1e-12)
    # Ostorero et al. (2022) eq. 1 in cm2/s
    assert np.log10(no * 1e4) == pytest.approx(-5.54 + 2.6 * 0.3 - 12530.0 / 1123.15, abs=1e-9)


def test_opx_ree_laws():
    T = 1000 + 273.15
    assert logD("opx_Ce_dias2025", 1000) == pytest.approx(
        np.log10(5.75e-14) - 166000.0 / (LN10 * R_GAS * T), abs=1e-9)
    assert logD("opx_Eu_dias2025", 1000) == pytest.approx(
        np.log10(1.90e-14) - 147000.0 / (LN10 * R_GAS * T), abs=1e-9)
    iw = log_fo2_buffer("IW", T, 1e5)
    at_iw = logD("opx_Lu_dias2025", 1000, log_fo2_bar=iw)
    assert at_iw == pytest.approx(np.log10(1.51e-9) - 263000.0 / (LN10 * R_GAS * T), abs=1e-6)
    assert logD("opx_Lu_dias2025", 1000, log_fo2_bar=iw + 7) == pytest.approx(at_iw + 1.0, abs=1e-6)
    assert not get("opx_Lu_dias2025").verified, "the Lu reference fO2 is an interpretation"


# --- magnetite: Sievwright et al. (2020) ------------------------------------------
def test_sievwright_table5_transcription_reproduces_the_published_minima():
    for sp, (lv, li, lf_min, ld_min) in SIEVWRIGHT_TABLE5.items():
        lf = 0.75 * (li - lv)                      # where the two branches are equal
        assert lf == pytest.approx(lf_min, abs=0.1), sp
        assert np.log10(sievwright_D_1150(sp, lf)) == pytest.approx(ld_min, abs=0.05), sp


def test_sievwright_entries_equal_table5_at_1150_and_agree_with_table12_for_ti():
    for sp in ("Ti", "Mn", "Co", "Cr", "Al", "Mg"):
        c = get(f"mt_{sp}_sievwright2020")
        for lf in (-9.0, -7.0, -5.0):
            D = c.D(Conditions(T_K=1423.15, log_fo2_bar=lf))
            assert D == pytest.approx(sievwright_D_1150(sp, lf), rel=1e-9)
    fmq = log_fo2_buffer("FMQ", 1423.15, 1e5)
    for d in (0.0, 2.0, 4.0):
        s20 = logD("mt_Ti_sievwright2020", 1150, log_fo2_bar=fmq + d)
        t12 = logD("mt_Ti_vanorman_crispin2010", 1150, log_fo2_bar=fmq + d, X={"xTi": 0.0})
        assert s20 == pytest.approx(t12, abs=0.5)


def test_sievwright_temperature_scaling_uses_the_table12_energies():
    c = get("mt_Ti_sievwright2020")
    lf = -10.0
    T = 1223.15
    lv, li = SIEVWRIGHT_TABLE5["Ti"][:2]
    _, QV, _, QI = TABLE12_PURE["Ti"]
    dinv = 1 / T - 1 / 1423.15
    expected = (10 ** (lv + lf * 2 / 3) * np.exp(-QV * 1e3 / R_GAS * dinv)
                + 10 ** (li - lf * 2 / 3) * np.exp(-QI * 1e3 / R_GAS * dinv))
    assert c.D(Conditions(T_K=T, log_fo2_bar=lf)) == pytest.approx(expected, rel=1e-12)
    # the Mg entry has no temperature dependence and must say so when used elsewhere
    w = get("mt_Mg_sievwright2020").check_conditions(Conditions(T_K=T, log_fo2_bar=lf))
    assert any("outside the calibration range" in s for s in w)


def test_the_old_magnetite_note_no_longer_says_sievwright_is_missing():
    assert "mt_*_sievwright2020" in get("mt_FeTi_freer_hauptman1978").superseded_note


# --- profile window ------------------------------------------------------------------
def test_fit_window_keeps_other_columns_aligned(tmp_path):
    f = tmp_path / "p.csv"
    f.write_text("d,a,an\n30,3,0.3\n0,0,0.0\n10,1,0.1\n20,2,0.2\n40,4,0.4\n", encoding="utf-8")
    df = read_table(f)
    prof = build_profile(df, ProfileSpec(distance_column="d", column_a="a", mode="A",
                                         x_min=5, x_max=35))
    assert prof.x == pytest.approx([10, 20, 30])
    assert prof.C == pytest.approx([1, 2, 3])
    assert prof.column("an") == pytest.approx([0.1, 0.2, 0.3])
    assert any("outside the fitted window" in n for n in prof.notes)
    assert "fitted window 5 to 35 um" in prof.spec.describe()
    with pytest.raises(ValueError):
        build_profile(df, ProfileSpec(distance_column="d", column_a="a", mode="A", x_min=100))
