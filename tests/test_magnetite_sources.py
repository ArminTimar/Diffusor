"""Magnetite laws against what their sources support (ranges, units, interpolation)."""
import numpy as np
import pytest

from diffusor.coefficients import get
from diffusor.coefficients.base import Conditions
from diffusor.coefficients.magnetite import (
    SIEVWRIGHT_TABLE5, TABLE12_PURE, TABLE12_XTI02, sievwright_D_1150)
from diffusor.thermo.buffers import log_fo2_from_delta


def _D(sp, T_C, lf=-11.0, x=0.1):
    c = get(f"mt_{sp}_vanorman_crispin2010")
    return float(c.D(Conditions(T_K=T_C + 273.15, log_fo2_bar=lf, X={"xTi": x})))


# --- Van Orman & Crispin (2010) Tables 10-12: data windows ----------------------------
@pytest.mark.parametrize("sp,lo,hi", [("Ti", 1473, 1573), ("Fe", 1373, 1573), ("Mn", 1373, 1573),
                                      ("Co", 1373, 1573), ("Cr", 1483, 1683), ("Al", 1553, 1773)])
def test_table12_temperature_range_is_the_data_window(sp, lo, hi):
    r = get(f"mt_{sp}_vanorman_crispin2010").T_range
    assert (r.lo, r.hi) == (lo, hi)
    assert "Table 12 prints no range" in r.unit


@pytest.mark.parametrize("key,lo,hi", [
    ("mt_Fe_aggarwal2002_WM", 1173, 1673), ("mt_Fe_aggarwal2002_MH", 1173, 1673),
    ("mt_Ti_aggarwal2002_WM", 1373, 1573), ("mt_Ti_aggarwal2002_MH", 1373, 1573),
    ("mt_Ti_aggarwal2002_WM_xti02", 1373, 1573), ("mt_Ti_aggarwal2002_MH_xti02", 1373, 1573),
    ("mt_Fe_aggarwal2002_WM_xti02", 1373, 1573), ("mt_Fe_aggarwal2002_MH_xti02", 1373, 1573)])
def test_aggarwal_buffer_path_ranges_are_those_of_tables_10_and_11(key, lo, hi):
    r = get(key).T_range
    assert (r.lo, r.hi) == (lo, hi)


def test_cr_and_al_are_pure_magnetite_laws_and_al_is_not_a_tracer_law():
    for sp in ("Cr", "Al"):
        c = get(f"mt_{sp}_vanorman_crispin2010")
        assert (c.X_range.lo, c.X_range.hi) == (0.0, 0.0)
        assert not c.requires
        # no x_Ti dependence in the law
        assert _D(sp, 1300.0, x=0.0) == _D(sp, 1300.0, x=0.2)
    al, cr = get("mt_Al_vanorman_crispin2010"), get("mt_Cr_vanorman_crispin2010")
    assert al.kind == "chemical" and "interdiffusion" in al.transported_variable
    assert cr.kind == "tracer" and "tracer" in cr.transported_variable
    assert "Petuskey" in al.notes


def test_cr_and_al_rows_are_the_1987_summarizing_expressions():
    # Dieckmann, Hilton & Mason (1987) p. 61: D0 [cm2/s], exponent [K]; Table 12 prints the
    # Al interstitial D0 as "6.92x10" without a visible exponent.
    R = 8.314462618
    expr = {"Cr": (5.12e-9, -880.0, 3.84e9, 90500.0), "Al": (3.24e-9, -7800.0, 6.92e7, 81900.0)}
    for sp, (dv, kv, di, ki) in expr.items():
        DV0, QV, DI0, QI = TABLE12_PURE[sp]
        assert DV0 == pytest.approx(dv * 1e-4) and DI0 == pytest.approx(di * 1e-4)
        assert QV == pytest.approx(kv * R / 1e3, abs=0.05) and QI == pytest.approx(ki * R / 1e3, abs=0.05)
    # and their Table 1 values for Al (log10 cm2/s) at 1280, 1380 and 1500 C
    for T_C, lv, li in [(1280, -6.32, -15.03), (1380, -6.42, -13.74), (1500, -6.59, -12.19)]:
        T = T_C + 273.15
        DV0, QV, DI0, QI = TABLE12_PURE["Al"]
        assert np.log10(DV0 * np.exp(-QV * 1e3 / (R * T)) * 1e4) == pytest.approx(lv, abs=0.06)
        assert np.log10(DI0 * np.exp(-QI * 1e3 / (R * T)) * 1e4) == pytest.approx(li, abs=0.07)


# --- x_Ti interpolation of Table 12 --------------------------------------------------
@pytest.mark.parametrize("sp", ["Ti", "Fe", "Co", "Mn"])
def test_table12_end_members_are_the_printed_rows(sp):
    R = 8.314462618
    for x, table in ((0.0, TABLE12_PURE), (0.2, TABLE12_XTI02)):
        DV, QV, DI, QI = table[sp]
        T, lf = 1450.0, -9.0
        a = 10.0 ** (lf - np.log10(1.01325))
        expected = (DV * np.exp(-QV * 1e3 / (R * T)) * a ** (2 / 3)
                    + DI * np.exp(-QI * 1e3 / (R * T)) * a ** (-2 / 3))
        got = float(get(f"mt_{sp}_vanorman_crispin2010").D(
            Conditions(T_K=T, log_fo2_bar=lf, X={"xTi": x})))
        assert got == pytest.approx(expected, rel=1e-3)


def test_table12_interpolates_each_branch_log_linearly_in_xti():
    """ln D0 and Q of each branch are linear in x_Ti, so at x_Ti = 0.1 each branch is the
    geometric mean of the two printed rows."""
    R = 8.314462618
    T, lf = 1250.0, -10.0
    a = 10.0 ** (lf - np.log10(1.01325))
    def branches(row):
        DV, QV, DI, QI = row
        return (DV * np.exp(-QV * 1e3 / (R * T)) * a ** (2 / 3),
                DI * np.exp(-QI * 1e3 / (R * T)) * a ** (-2 / 3))
    v0, i0 = branches(TABLE12_PURE["Fe"])
    v1, i1 = branches(TABLE12_XTI02["Fe"])
    expected = np.sqrt(v0 * v1) + np.sqrt(i0 * i1)
    assert _D("Fe", T - 273.15, lf, 0.1) == pytest.approx(expected, rel=1e-3)


def test_table12_reproduces_all_four_published_tomiya_values():
    """Tomiya et al. (2013), log fO2 = -11, X_Usp = 0.3 (x_Ti = 0.1): Ti 4.3e-16 (950 C) and
    6.9e-16 (900 C); Fe (their Mg) 4.4e-15 and 6.6e-15 m2/s."""
    for sp, d950, d900 in (("Ti", 4.3e-16, 6.9e-16), ("Fe", 4.4e-15, 6.6e-15)):
        assert _D(sp, 950.0) == pytest.approx(d950, rel=0.06)
        assert _D(sp, 900.0) == pytest.approx(d900, rel=0.06)


def test_minimum_positions_quoted_in_the_notes():
    T = np.arange(700.0, 1300.0, 1.0)
    for sp, tmin, band in (("Ti", 1023.0, (930.0, 1070.0)), ("Fe", 955.0, (880.0, 1000.0))):
        d = np.array([_D(sp, t) for t in T])
        assert T[d.argmin()] == pytest.approx(tmin, abs=3.0), sp
        inside = T[d <= 2.0 * d.min()]
        assert inside.min() == pytest.approx(band[0], abs=5.0), sp
        assert inside.max() == pytest.approx(band[1], abs=5.0), sp
    # Tomiya et al. place the Fe minimum at about 950 C
    assert 930.0 < T[np.array([_D("Fe", t) for t in T]).argmin()] < 980.0


def test_table12_notes_do_not_claim_to_know_how_tomiya_interpolated():
    for sp in ("Ti", "Fe", "Co", "Mn", "Cr", "Al"):
        n = get(f"mt_{sp}_vanorman_crispin2010").notes
        assert "which is how Tomiya" not in n
        assert "evidently" not in n
    assert "do not state the fO2 unit or how they obtained the x_Ti dependence" in \
        get("mt_Ti_vanorman_crispin2010").notes
    assert "D_Al = D_Ti" in get("mt_Ti_vanorman_crispin2010").notes


# --- Sievwright et al. (2020) --------------------------------------------------------
def test_sievwright_fo2_range_uses_oneill_fmq_as_in_the_paper():
    lo = log_fo2_from_delta("FMQ", -1.0, 1423.15, 1e5, "oneill")
    hi = log_fo2_from_delta("FMQ", 4.89, 1423.15, 1e5, "oneill")
    for sp in SIEVWRIGHT_TABLE5:
        c = get(f"mt_{sp}_sievwright2020")
        assert c.fo2_unit == "bar"
        assert c.fo2_range.lo == pytest.approx(lo, abs=0.01) and c.fo2_range.lo <= lo
        assert c.fo2_range.hi == pytest.approx(hi, abs=0.01) and c.fo2_range.hi >= hi
        assert "O'Neill" in c.fo2_range.unit
        assert "fO2 in bar" in c.notes or "takes fO2 in bar" in c.notes


def test_sievwright_table5_minima_within_the_stated_tolerances():
    worst_logD = worst_logf = 0.0
    for sp, (lv, li, lf_min, ld_min) in SIEVWRIGHT_TABLE5.items():
        lf = 0.75 * (li - lv)
        worst_logf = max(worst_logf, abs(lf - lf_min))
        worst_logD = max(worst_logD, abs(np.log10(sievwright_D_1150(sp, lf)) - ld_min))
    assert worst_logD < 0.04 and worst_logf < 0.09


def test_sievwright_vs_table12_differences_quoted_in_the_notes():
    lo = log_fo2_from_delta("FMQ", -1.0, 1423.15, 1e5, "oneill")
    hi = log_fo2_from_delta("FMQ", 4.89, 1423.15, 1e5, "oneill")
    lfs = np.linspace(lo, hi, 300)
    limits = {"Ti": 0.51, "Mn": 0.15, "Co": 0.23, "Cr": 0.44, "Al": 0.2}
    for sp, lim in limits.items():
        s, v = get(f"mt_{sp}_sievwright2020"), get(f"mt_{sp}_vanorman_crispin2010")
        d = [abs(np.log10(s.D(Conditions(1423.15, log_fo2_bar=l)))
                 - np.log10(v.D(Conditions(1423.15, log_fo2_bar=l, X={"xTi": 0.0})))) for l in lfs]
        assert max(d) < lim, sp


# --- Fe-Ti interdiffusion ------------------------------------------------------------
def test_freer_hauptman_composition_range_is_the_data_range():
    c = get("mt_FeTi_freer_hauptman1978")
    assert c.X_range.hi == 0.05            # tracer to 5 mol% Ti; TM20 = 0.2/3 = 0.067
    assert 0.2 / 3 > c.X_range.hi
    assert (c.T_range.lo, c.T_range.hi) == pytest.approx((888 + 273.15, 1034 + 273.15))
    warned = c.check_conditions(Conditions(T_K=1200.0, X={"xTi": 0.15}))
    assert any("x_Ti" in w and "outside the calibration range" in w for w in warned)
    assert "saunders2012" not in c.secondary_citations
    assert "order-of-magnitude agreement" in c.notes


def test_aragon_range_is_the_range_of_the_anneals():
    c = get("mt_FeTi_aragon1984")
    assert (c.T_range.lo, c.T_range.hi) == pytest.approx((990 + 273.15, 1220 + 273.15))
    assert "differ by about an order of magnitude" not in c.notes
    assert "order of magnitude agreement" in c.notes
    assert not c.verified


# --- Aggarwal & Dieckmann (2002a,b): the data behind Table 12 and the x_Ti interpolation ---
_AD_TABLE5 = {   # (T/C, x) -> species -> (log10 D_V, log10 D_I); Aggarwal & Dieckmann (2002b) Table 5, own values
    (1100, 0.2): {'Fe': (-6.57, -19.96), 'Co': (-6.699, -20.018), 'Mn': (-6.629, -20.16), 'Ti': (-7.431, -22.394)},
    (1200, 0.0): {'Fe': (-9.311, -17.975), 'Ti': (-10.021, -20.096)},
    (1200, 0.1): {'Fe': (-8.353, -17.994), 'Co': (-8.419, -17.855), 'Mn': (-8.363, -18.122), 'Ti': (-9.207, -19.418)},
    (1200, 0.2): {'Fe': (-6.768, -18.707), 'Co': (-6.845, -18.708), 'Mn': (-6.757, -18.908), 'Ti': (-7.763, -20.445)},
    (1200, 0.3): {'Fe': (-5.136, -19.308), 'Co': (-5.14, -19.31), 'Mn': (-5.177, -19.59), 'Ti': (-5.873, -20.867)},
    (1300, 0.0): {'Fe': (-9.672, -16.453), 'Ti': (-10.339, -18.089)},
    (1300, 0.1): {'Fe': (-8.634, -16.489), 'Co': (-8.657, -16.431), 'Mn': (-8.635, -16.636), 'Ti': (-9.15, -17.462)},
    (1300, 0.2): {'Fe': (-7.391, -17.231), 'Co': (-7.426, -17.18), 'Mn': (-7.366, -17.401), 'Ti': (-8.062, -18.319)},
    (1300, 0.3): {'Fe': (-6.249, -17.797), 'Co': (-6.32, -17.819), 'Mn': (-6.254, -18.001), 'Ti': (-6.951, -19.034)},
    (1400, 0.0): {'Fe': (-9.905, -15.317), 'Ti': (-10.474, -17.196)},
}

_AD_RAW = {   # (T/C, x) -> [(log10 aO2, {species: log10 D*})]; Aggarwal & Dieckmann (2002b) Tables 1-4
    (1100, 0.2): [
        (-8.65, {'Fe': -12.338, 'Co': -12.485, 'Mn': -12.441, 'Ti': -13.163}),
        (-9.03, {'Fe': -12.67, 'Co': -12.789, 'Mn': -12.695}),
        (-9.55, {'Ti': -13.83}),
        (-10.2, {'Fe': -12.938, 'Co': -12.978, 'Mn': -13.075, 'Ti': -14.249}),
        (-10.74, {'Fe': -12.721, 'Co': -12.812, 'Mn': -12.919, 'Ti': -14.571}),
        (-11.5, {'Fe': -12.32, 'Co': -12.39, 'Mn': -12.511, 'Ti': -14.553}),
    ],
    (1200, 0.0): [
        (-3.21, {'Fe': -11.536, 'Ti': -12.31}),
        (-3.37, {'Fe': -11.677, 'Ti': -12.494}),
        (-4.25, {'Fe': -12.171, 'Ti': -12.844}),
        (-5.25, {'Fe': -12.86, 'Ti': -13.523}),
        (-7.19, {'Fe': -13.147, 'Ti': -14.603}),
        (-8.29, {'Fe': -12.462, 'Ti': -14.383}),
        (-8.75, {'Fe': -12.113, 'Ti': -14.407}),
    ],
    (1200, 0.1): [
        (-5.64, {'Ti': -12.915}),
        (-5.95, {'Fe': -12.315, 'Co': -12.402, 'Mn': -12.33}),
        (-6.48, {'Ti': -13.633}),
        (-6.5, {'Fe': -12.661, 'Co': -12.703, 'Mn': -12.71}),
        (-6.61, {'Fe': -12.753, 'Co': -12.755, 'Mn': -12.74}),
        (-7.51, {'Fe': -12.837, 'Co': -12.76, 'Mn': -12.946}),
        (-7.53, {'Ti': -14.037}),
        (-8.54, {'Fe': -12.301, 'Co': -12.163, 'Mn': -12.424}),
        (-8.56, {'Ti': -13.73}),
        (-9.04, {'Ti': -13.327}),
        (-9.18, {'Fe': -11.858, 'Co': -11.715, 'Mn': -11.981}),
    ],
    (1200, 0.2): [
        (-7.63, {'Fe': -11.903, 'Co': -12.023, 'Mn': -11.932, 'Ti': -12.949}),
        (-7.91, {'Ti': -13.108}),
        (-7.95, {'Fe': -12.096, 'Co': -12.153, 'Mn': -12.086}),
        (-8.57, {'Fe': -12.392, 'Co': -12.397, 'Mn': -12.379, 'Ti': -13.404}),
        (-9.51, {'Fe': -12.296}),
        (-9.55, {'Fe': -12.257, 'Co': -12.28, 'Mn': -12.414}),
        (-9.92, {'Ti': -13.684}),
        (-9.98, {'Fe': -12.043, 'Co': -12.047, 'Mn': -12.257}),
        (-10.3, {'Fe': -11.816, 'Co': -11.845, 'Mn': -12.026}),
        (-10.31, {'Fe': -11.848, 'Ti': -13.573}),
    ],
    (1200, 0.3): [
        (-10.98, {'Fe': -11.894, 'Ti': -13.024}),
        (-11.02, {'Fe': -11.798, 'Co': -11.848, 'Mn': -12.024}),
        (-11.07, {'Fe': -11.882, 'Co': -11.931, 'Mn': -12.13}),
        (-11.21, {'Fe': -11.78, 'Co': -11.817, 'Mn': -12.026, 'Ti': -13.122}),
        (-11.37, {'Fe': -11.624, 'Co': -11.648, 'Mn': -11.877}),
        (-11.42, {'Fe': -11.623, 'Co': -11.717, 'Mn': -11.878}),
        (-11.45, {'Fe': -11.623, 'Ti': -13.029}),
        (-11.82, {'Fe': -11.486, 'Co': -11.418, 'Mn': -11.743}),
    ],
    (1300, 0.0): [
        (-1.64, {'Ti': -11.745}),
        (-2.63, {'Fe': -11.535, 'Ti': -12.182}),
        (-6.04, {'Fe': -12.299, 'Ti': -13.928}),
        (-7.34, {'Fe': -11.659, 'Ti': -13.164}),
    ],
    (1300, 0.1): [
        (-4.59, {'Fe': -11.818, 'Co': -11.851, 'Mn': -11.84, 'Ti': -12.255}),
        (-5.35, {'Fe': -12.1, 'Co': -12.085, 'Mn': -12.096, 'Ti': -12.795}),
        (-6.51, {'Fe': -12.066, 'Co': -12.057, 'Mn': -12.242, 'Ti': -12.847}),
        (-7.32, {'Fe': -11.631, 'Co': -11.545, 'Mn': -11.735, 'Ti': -12.651}),
    ],
    (1300, 0.2): [
        (-6.32, {'Fe': -11.722, 'Co': -11.765, 'Mn': -11.699, 'Ti': -12.357}),
        (-6.97, {'Ti': -12.725}),
        (-7.07, {'Fe': -11.958, 'Co': -11.947, 'Mn': -11.994}),
        (-8.22, {'Fe': -11.72, 'Co': -11.678, 'Mn': -11.877, 'Ti': -12.834}),
        (-8.53, {'Fe': -11.59, 'Co': -11.528, 'Mn': -11.725}),
        (-8.86, {'Fe': -11.273, 'Co': -11.239, 'Mn': -11.46, 'Ti': -12.337}),
    ],
    (1300, 0.3): [
        (-8.18, {'Fe': -11.782, 'Co': -11.839, 'Mn': -11.833, 'Ti': -12.563}),
        (-9.42, {'Fe': -11.475, 'Co': -11.502, 'Mn': -11.663, 'Ti': -12.709}),
        (-10.0, {'Fe': -11.128, 'Co': -11.15, 'Mn': -11.326, 'Ti': -12.292}),
    ],
    (1400, 0.0): [
        (-1.4, {'Fe': -11.01, 'Ti': -11.668}),
        (-3.6, {'Fe': -12.245, 'Ti': -12.81}),
        (-4.83, {'Fe': -12.066}),
        (-5.44, {'Fe': -11.672, 'Ti': -13.459}),
    ],
}

# literature entries (superscripts a, b, c) of Table 5 at x = 0, (T/C) -> species -> [(log D_V, log D_I)]
_AD_LIT_X0 = {
    1100: {"Fe": [(-9.097, -19.458)], "Co": [(-9.148, -19.373)], "Mn": [(-8.917, -19.564)]},
    1200: {"Fe": [(-9.464, -17.866), (-9.403, -18.082), (-9.44, -18.000)],
           "Co": [(-9.522, -17.795), (-9.53, -17.76)], "Mn": [(-9.479, -18.140), (-9.49, -18.03)]},
    1300: {"Fe": [(-9.739, -16.444)], "Co": [(-9.848, -16.414)], "Mn": [(-9.764, -16.634)]},
    1400: {"Fe": [(-9.943, -15.357)]},
}
# Part I (2002a) Table 4: (T/C) -> x -> (log10 [V]0, log10 [I]0)
_AD_PART1_T4 = {1100: {0.0: (1.301, -9.444), 0.1: (2.779, -10.891), 0.2: (3.646, -11.086), 0.25: (4.777, -11.728)}, 1200: {0.0: (0.699, -8.495), 0.1: (2.076, -9.607), 0.2: (3.19, -10.017), 0.25: (3.962, -10.2)}, 1300: {0.0: (0.197, -7.357), 0.1: (1.553, -7.941), 0.2: (2.554, -8.612), 0.25: (3.335, -9.079)}}

_R = 8.314462618


def _refit(sp, x):
    """Fit ln D = ln D0 - Q/RT to the Table 5 values (own measurements plus, at x = 0, the
    tabulated literature values); returns ((D_V0, Q_V), (D_I0, Q_I)) in m2/s and kJ/mol."""
    pts = [(T + 273.15, d[sp]) for (T, xx), d in _AD_TABLE5.items() if xx == x and sp in d]
    if x == 0.0:
        pts += [(T + 273.15, v) for T, d in _AD_LIT_X0.items() for v in d.get(sp, [])]
    T = np.array([p[0] for p in pts])
    out = []
    for k in (0, 1):
        y = np.array([p[1][k] for p in pts]) * np.log(10.0)
        A = np.vstack([np.ones_like(T), -1.0 / (_R * T)]).T
        c = np.linalg.lstsq(A, y, rcond=None)[0]
        out.append((float(np.exp(c[0])), float(c[1]) / 1e3))
    return out


@pytest.mark.parametrize("sp", ["Fe", "Co", "Mn", "Ti"])
@pytest.mark.parametrize("x,table", [(0.0, TABLE12_PURE), (0.2, TABLE12_XTI02)])
def test_table12_rows_are_refits_of_aggarwal_dieckmann_table5(sp, x, table):
    """The Table 12 footnote says the rows were calculated from all data in Table 5 of
    Aggarwal & Dieckmann (2002b); refitting that table reproduces them to the printed digits."""
    (dv, qv), (di, qi) = _refit(sp, x)
    DV, QV, DI, QI = table[sp]
    assert dv == pytest.approx(DV, rel=0.01) and di == pytest.approx(DI, rel=0.01)
    assert qv == pytest.approx(QV, abs=0.2) and qi == pytest.approx(QI, abs=0.2)


def test_table12_data_windows_match_the_temperatures_of_table5():
    """Pure magnetite: Ti 1200-1400 C only (not the 1373-1573 K printed in Table 10); Fe 1100-1400 C;
    Co and Mn 1100-1300 C.  x_Ti = 0.2: 1100-1300 C."""
    for sp, (lo, hi) in {"Fe": (1100, 1400), "Co": (1100, 1300), "Mn": (1100, 1300), "Ti": (1200, 1400)}.items():
        temps = {T for T, d in _AD_LIT_X0.items() if sp in d}
        temps |= {T for (T, x), d in _AD_TABLE5.items() if x == 0.0 and sp in d}
        assert (min(temps), max(temps)) == (lo, hi), sp
        assert {T for (T, x), d in _AD_TABLE5.items() if x == 0.2 and sp in d} == {1100, 1200, 1300}
    # the registered windows are the overlap of those data windows
    for sp, lo, hi in (("Fe", 1373, 1573), ("Co", 1373, 1573), ("Mn", 1373, 1573), ("Ti", 1473, 1573)):
        r = get(f"mt_{sp}_vanorman_crispin2010").T_range
        assert (r.lo, r.hi) == (lo, hi), sp
        assert "Aggarwal & Dieckmann" in r.unit


def _cv(T_C, x, loga):
    """Correction factor C_V of Aggarwal & Dieckmann (2002a) eqs. 5, 12, 15 (not part of Table 12)."""
    tab = _AD_PART1_T4[T_C]
    xs = sorted(tab)

    def at(k):
        ys = [tab[q][k] for q in xs]
        if x <= xs[-1]:
            return np.interp(x, xs, ys)
        return ys[-1] + (ys[-1] - ys[-2]) / (xs[-1] - xs[-2]) * (x - xs[-1])

    V, I = 10.0 ** at(0), 10.0 ** at(1)
    a23 = 10.0 ** (loga * 2.0 / 3.0)

    def C(d):
        return (((1 + 3 * x - (3 + x) * d) ** 3 / (2 - 6 * x + (2 + 2 * x) * d) ** 2)
                * ((2 - 6 * x) ** 2 / (1 + 3 * x) ** 3))

    lo, hi = -0.1, 0.2
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        glo = V * C(lo) * a23 - I / a23 - lo
        gmid = V * C(mid) * a23 - I / a23 - mid
        lo, hi = (lo, mid) if glo * gmid <= 0 else (mid, hi)
    return C(0.5 * (lo + hi))


def _branch_logs(sp, T_K, x):
    """log10 of the (vacancy, interstitial) terms at a_O2 = 1, each branch log-interpolated in x."""
    def br(tab):
        DV, QV, DI, QI = tab[sp]
        return (np.log(DV) - QV * 1e3 / (_R * T_K), np.log(DI) - QI * 1e3 / (_R * T_K))
    b0, b1 = br(TABLE12_PURE), br(TABLE12_XTI02)
    f = x / 0.2
    return [((1 - f) * b0[k] + f * b1[k]) / np.log(10.0) for k in (0, 1)]


def _scheme_logD(sp, T_C, x, loga, scheme, cv=False):
    """log10 D (m2/s) at a_O2 = 10**loga: 'branch' = Diffusor's scheme, 'total' = log-linear in the total D."""
    T = T_C + 273.15
    a23 = 10.0 ** (loga * 2.0 / 3.0)
    c = _cv(T_C, x, loga) if cv else 1.0
    if scheme == "branch":
        lv, li = _branch_logs(sp, T, x)
        return np.log10(10.0 ** lv * a23 * c + 10.0 ** li / a23)
    ends = []
    for tab in (TABLE12_PURE, TABLE12_XTI02):
        DV, QV, DI, QI = tab[sp]
        ends.append(np.log10(DV * np.exp(-QV * 1e3 / (_R * T)) * a23 * c
                             + DI * np.exp(-QI * 1e3 / (_R * T)) / a23))
    f = x / 0.2
    return (1 - f) * ends[0] + f * ends[1]


def _residuals(xsel, scheme, cv=False):
    r = []
    for (T, x), rows in _AD_RAW.items():
        if x not in xsel or T not in (1200, 1300):
            continue
        for loga, d in rows:
            for sp, obs in d.items():
                if sp in TABLE12_XTI02:
                    r.append(_scheme_logD(sp, T, x, loga, scheme, cv) - obs)
    return np.array(r)


def test_part1_correction_and_table5_reproduce_the_raw_tracer_data():
    """Validates the transcription: Table 5 with the C_V of Part I reproduces Tables 2-3 to 0.06."""
    r = []
    for (T, x), rows in _AD_RAW.items():
        if T not in (1200, 1300):
            continue
        for loga, d in rows:
            for sp, obs in d.items():
                if sp in _AD_TABLE5[(T, x)]:
                    lv, li = _AD_TABLE5[(T, x)][sp]
                    a23 = 10.0 ** (loga * 2.0 / 3.0)
                    r.append(np.log10(10.0 ** lv * a23 * _cv(T, x, loga) + 10.0 ** li / a23) - obs)
    r = np.array(r)
    assert len(r) == 139
    assert np.sqrt(np.mean(r ** 2)) < 0.06


def test_interpolation_scheme_against_aggarwal_dieckmann_at_x01_and_x03():
    """The numbers quoted in the notes and in the comment above _make_table12_func."""
    r = _residuals((0.0, 0.2), "branch")                  # the end members themselves
    assert len(r) == 65 and np.sqrt(np.mean(r ** 2)) == pytest.approx(0.125, abs=0.01)
    assert np.sqrt(np.mean(_residuals((0.0, 0.2), "branch", cv=True) ** 2)) == pytest.approx(0.11, abs=0.01)
    rb, rt = _residuals((0.1,), "branch"), _residuals((0.1,), "total")
    assert len(rb) == len(rt) == 39
    assert (rb.mean(), np.sqrt(np.mean(rb ** 2)), np.abs(rb).max()) == pytest.approx((-0.09, 0.30, 0.76), abs=0.01)
    assert (rt.mean(), np.sqrt(np.mean(rt ** 2)), np.abs(rt).max()) == pytest.approx((0.12, 0.29, 0.52), abs=0.01)
    assert abs(np.sqrt(np.mean(rb ** 2)) - np.sqrt(np.mean(rt ** 2))) < 0.03     # the data do not decide
    rb3, rt3 = _residuals((0.3,), "branch"), _residuals((0.3,), "total")           # extrapolation
    assert len(rb3) == len(rt3) == 35
    assert (rb3.mean(), np.sqrt(np.mean(rb3 ** 2))) == pytest.approx((0.31, 0.34), abs=0.01)
    assert (rt3.mean(), np.sqrt(np.mean(rt3 ** 2))) == pytest.approx((0.22, 0.33), abs=0.01)


def test_branch_wise_interpolation_against_the_partial_coefficients_at_x01():
    """Vacancy branch within 0.3 log units, interstitial branch 0.26-0.94 too low (notes)."""
    dv, di = [], []
    for T in (1200, 1300):
        for sp, (lv, li) in _AD_TABLE5[(T, 0.1)].items():
            pv, pi = _branch_logs(sp, T + 273.15, 0.1)
            dv.append(pv - lv)
            di.append(pi - li)
    dv, di = np.array(dv), np.array(di)
    assert np.abs(dv).max() < 0.31
    assert np.all(di < 0)
    assert (-di).min() == pytest.approx(0.26, abs=0.01) and (-di).max() == pytest.approx(0.94, abs=0.01)
    # the partial coefficients are not log-linear in x: the Fe interstitial branch is flat to 0.1, then falls
    li = [_AD_TABLE5[(1200, x)]["Fe"][1] for x in (0.0, 0.1, 0.2, 0.3)]
    assert abs(li[1] - li[0]) < 0.05 and li[2] - li[1] < -0.6


def test_minimum_depth_and_position_against_aggarwal_dieckmann():
    g = np.arange(-16.0, -2.0, 0.01)

    def minimum(sp, T_C, x, scheme):
        d = np.array([_scheme_logD(sp, T_C, x, l, scheme) for l in g])
        k = d.argmin()
        return g[k], d[k]

    for T, tol_branch in ((1200, 0.06), (1300, 0.16)):          # x = 0.1, Fe: depth of the minimum, eq. 8
        lv, li = _AD_TABLE5[(T, 0.1)]["Fe"]
        measured = np.log10(2.0) + 0.5 * (lv + li)
        assert abs(minimum("Fe", T, 0.1, "branch")[1] - measured) < tol_branch
        assert abs(minimum("Fe", T, 0.1, "total")[1] - measured) > 0.3
    for T, sp in ((1200, "Fe"), (1200, "Ti"), (1300, "Fe"), (1300, "Ti")):    # x = 0.3, position, eq. 6
        lv, li = _AD_TABLE5[(T, 0.3)][sp]
        measured = 0.75 * (li - lv)
        e_branch = abs(minimum(sp, T, 0.3, "branch")[0] - measured)
        e_total = abs(minimum(sp, T, 0.3, "total")[0] - measured)
        assert e_branch < 0.95 and e_total > 0.9 and e_branch < e_total


def test_the_two_interpolation_forms_differ_by_the_quoted_amounts():
    def worst(sp, x):
        return max(abs(_scheme_logD(sp, T, x, l, "branch") - _scheme_logD(sp, T, x, l, "total"))
                   for T in np.arange(850.0, 1301.0, 25.0) for l in np.arange(-16.0, -4.9, 0.25))
    assert worst("Fe", 0.1) == pytest.approx(0.50, abs=0.03)
    assert worst("Ti", 0.1) == pytest.approx(0.84, abs=0.03)


def test_tomiya_values_with_both_interpolation_forms():
    a = -11.0 - np.log10(1.01325)
    tomiya = {("Ti", 950): 4.3e-16, ("Ti", 900): 6.9e-16, ("Fe", 950): 4.4e-15, ("Fe", 900): 6.6e-15}
    dev_b = {k: _scheme_logD(k[0], k[1], 0.1, a, "branch") - np.log10(v) for k, v in tomiya.items()}
    dev_t = {k: _scheme_logD(k[0], k[1], 0.1, a, "total") - np.log10(v) for k, v in tomiya.items()}
    assert all(abs(v) < np.log10(1.06) for v in dev_b.values())
    assert abs(dev_t[("Ti", 950)]) < np.log10(1.02) and abs(dev_t[("Fe", 900)]) < np.log10(1.11)
    assert 10 ** dev_t[("Fe", 950)] == pytest.approx(2.13, abs=0.05)
    # linear interpolation of the printed parameters (D0 and Q themselves) is 0.2-0.5 log units off
    for (sp, TC), v in tomiya.items():
        T = TC + 273.15
        p = [0.5 * (u + w) for u, w in zip(TABLE12_PURE[sp], TABLE12_XTI02[sp])]
        a23 = 10.0 ** (a * 2.0 / 3.0)
        D = p[0] * np.exp(-p[1] * 1e3 / (_R * T)) * a23 + p[2] * np.exp(-p[3] * 1e3 / (_R * T)) / a23
        assert 0.2 < np.log10(D / v) < 0.5


def test_composition_notes_state_that_the_scheme_is_diffusors_own():
    for sp in ("Ti", "Fe", "Co", "Mn"):
        c = get(f"mt_{sp}_vanorman_crispin2010")
        assert "Diffusor's own construction" in c.notes
        assert "no source gives D between them" in c.notes
        assert "do not decide" in c.notes
        assert "Refitting ln D = ln D0 - Q/RT" in c.notes
        assert "aggarwal_dieckmann2002a" in c.secondary_citations
        assert "aggarwal_dieckmann2002" in c.secondary_citations
        assert "unavailable" not in c.notes.lower()


# --- node scheme: the measured x_Ti = 0.1 and 0.3 coefficients as Arrhenius nodes (evaluated, not adopted) ---
# Q_V, Q_I (kJ/mol) that Aggarwal & Dieckmann (2002b) print in their Table 6 (p. 716) for x = 0.1 and 0.3
_AD_TABLE6_Q = {("Fe", 0.1): (-125.0, 667.7), ("Fe", 0.3): (-493.8, 670.4),
                ("Co", 0.1): (-105.4, 631.7), ("Co", 0.3): (-409.7, 661.7),
                ("Mn", 0.1): (-120.7, 659.5), ("Mn", 0.3): (-477.7, 705.1),
                ("Ti", 0.1): (25.3, 867.8), ("Ti", 0.3): (-478.0, 813.3)}


def _two_point(sp, x):
    """(D_V0, Q_V, D_I0, Q_I) from the 1200 and 1300 C values of Table 5, the only two temperatures there."""
    T = np.array([1473.15, 1573.15])
    out = []
    for k in (0, 1):
        y = np.array([_AD_TABLE5[(1200, x)][sp][k], _AD_TABLE5[(1300, x)][sp][k]]) * np.log(10.0)
        c = np.linalg.solve(np.vstack([np.ones(2), -1.0 / (_R * T)]).T, y)
        out += [float(np.exp(c[0])), float(c[1]) / 1e3]
    return tuple(out)


def _nodes():
    n = {}
    for sp in ("Fe", "Co", "Mn", "Ti"):
        n[(sp, 0.0)], n[(sp, 0.2)] = TABLE12_PURE[sp], TABLE12_XTI02[sp]
        n[(sp, 0.1)], n[(sp, 0.3)] = _two_point(sp, 0.1), _two_point(sp, 0.3)
    return n


_NODES = _nodes()


def _node_logD(sp, T_C, x, loga):
    xs = [0.0, 0.1, 0.2, 0.3]
    i = min(max(int(np.searchsorted(xs, x, side="right")) - 1, 0), 2)
    f = (x - xs[i]) / 0.1
    T = T_C + 273.15
    a23 = 10.0 ** (loga * 2.0 / 3.0)

    def br(par):
        DV, QV, DI, QI = par
        return np.log(DV) - QV * 1e3 / (_R * T), np.log(DI) - QI * 1e3 / (_R * T)

    b0, b1 = br(_NODES[(sp, xs[i])]), br(_NODES[(sp, xs[i + 1])])
    return np.log10(np.exp((1 - f) * b0[0] + f * b1[0]) * a23 + np.exp((1 - f) * b0[1] + f * b1[1]) / a23)


def _anchored_logD(sp, T_C, x, loga, anchor_C):
    """Branch-wise 0/0.2 interpolation with each branch rescaled so that it equals the measured value at anchor_C."""
    f = x / 0.2
    a23 = 10.0 ** (loga * 2.0 / 3.0)

    def interp(T):
        b0, b1 = [], []
        for tab in (TABLE12_PURE, TABLE12_XTI02):
            DV, QV, DI, QI = tab[sp]
            (b0 if tab is TABLE12_PURE else b1).extend([np.log(DV) - QV * 1e3 / (_R * T), np.log(DI) - QI * 1e3 / (_R * T)])
        return [(1 - f) * b0[k] + f * b1[k] for k in (0, 1)]

    here, ref = interp(T_C + 273.15), interp(anchor_C + 273.15)
    meas = [v * np.log(10.0) for v in _AD_TABLE5[(anchor_C, x)][sp]]
    return np.log10(np.exp(here[0] + meas[0] - ref[0]) * a23 + np.exp(here[1] + meas[1] - ref[1]) / a23)


def test_node_parameters_reproduce_table6_activation_energies():
    for (sp, x), (qv, qi) in _AD_TABLE6_Q.items():
        _, QV, _, QI = _NODES[(sp, x)]
        if (sp, x) == ("Co", 0.3):
            assert QV == pytest.approx(-523.5, abs=0.5) and abs(QV - qv) > 100     # A&D print -409.7 (error unknown)
        else:
            assert QV == pytest.approx(qv, abs=0.5), (sp, x)
        assert QI == pytest.approx(qi, abs=0.5), (sp, x)
    # the physically odd energies that make the nodes unsafe to extrapolate
    assert _NODES[("Ti", 0.1)][1] > 0
    assert all(_NODES[(sp, 0.3)][1] < -400 for sp in ("Fe", "Co", "Mn", "Ti"))


def test_node_scheme_fits_its_own_data_but_not_the_low_temperature_law():
    def rms(xsel):
        r = []
        for (T, x), rows in _AD_RAW.items():
            if x != xsel or T not in (1200, 1300):
                continue
            for loga, d in rows:
                r += [_node_logD(sp, T, x, loga) - obs for sp, obs in d.items() if sp in TABLE12_XTI02]
        return np.sqrt(np.mean(np.square(r)))
    assert rms(0.1) == pytest.approx(0.06, abs=0.01) and rms(0.3) == pytest.approx(0.07, abs=0.01)
    # ... but at 900-1000 C it departs from the 0/0.2 laws
    a = [-11.0 - np.log10(1.01325), -13.0 - np.log10(1.01325)]
    d01 = {sp: [_node_logD(sp, T, 0.1, l) - _scheme_logD(sp, T, 0.1, l, "branch")
                for T in (900.0, 950.0, 1000.0) for l in a] for sp in ("Fe", "Co", "Mn", "Ti")}
    for sp in ("Fe", "Co", "Mn"):
        assert -0.58 < min(d01[sp]) and max(d01[sp]) < 0.15
    assert min(d01["Ti"]) == pytest.approx(-1.61, abs=0.03) and max(d01["Ti"]) < 0.15
    d03 = [_node_logD(sp, T, 0.3, l) - _scheme_logD(sp, T, 0.3, l, "branch")
           for sp in ("Fe", "Co", "Mn", "Ti") for T in (900.0, 1000.0) for l in a]
    assert min(d03) > 1.8 and max(d03) < 4.1


def test_node_scheme_does_not_reproduce_tomiya():
    a = -11.0 - np.log10(1.01325)
    tomiya = {("Ti", 950): 4.3e-16, ("Ti", 900): 6.9e-16, ("Fe", 950): 4.4e-15, ("Fe", 900): 6.6e-15}
    dev = {k: _node_logD(k[0], k[1], 0.1, a) - np.log10(v) for k, v in tomiya.items()}
    assert dev[("Ti", 950)] == pytest.approx(-1.35, abs=0.02) and dev[("Ti", 900)] == pytest.approx(-1.61, abs=0.02)
    assert dev[("Fe", 950)] == pytest.approx(-0.35, abs=0.02) and dev[("Fe", 900)] == pytest.approx(-0.38, abs=0.02)


def test_measured_amplitude_helps_at_x01_out_of_sample_but_not_at_x03():
    """Rescale the 0/0.2 branches to the measured coefficients at one temperature and predict the other
    (1200 <-> 1300 C): better than the unanchored interpolations at x_Ti = 0.1, not at 0.3."""
    def rms(x, anchor, target, fn):
        r = []
        for loga, d in _AD_RAW[(target, x)]:
            r += [fn(sp, target, x, loga) - obs for sp, obs in d.items() if sp in TABLE12_XTI02]
        return np.sqrt(np.mean(np.square(r)))
    for anchor, target in ((1200, 1300), (1300, 1200)):
        h = rms(0.1, anchor, target, lambda sp, T, x, l: _anchored_logD(sp, T, x, l, anchor))
        b = rms(0.1, anchor, target, lambda sp, T, x, l: _scheme_logD(sp, T, x, l, "branch"))
        t = rms(0.1, anchor, target, lambda sp, T, x, l: _scheme_logD(sp, T, x, l, "total"))
        assert h < 0.19 and min(b, t) > 0.25
    h3 = rms(0.3, 1200, 1300, lambda sp, T, x, l: _anchored_logD(sp, T, x, l, 1200))
    b3 = rms(0.3, 1200, 1300, lambda sp, T, x, l: _scheme_logD(sp, T, x, l, "branch"))
    assert h3 > b3


def test_sievwright_notes_use_the_supplement():
    for sp in ("Mn", "Ti", "Co"):
        n = get(f"mt_{sp}_sievwright2020").notes
        assert "Supplementary Table S1" in n and "carry no published uncertainty" in n
    n = get("mt_Mn_sievwright2020").notes
    assert "Supplementary Fig. S3" in n and "within 0.43 log units of the Mn curve" in n
    assert "Supplementary Table S1" in get("mt_Mn_sievwright2020_1150").notes
    assert "Their Supplementary Table S1" in get("mt_Mn_sievwright2020").sigma_logD_basis


def test_node_scheme_is_documented_as_evaluated_and_not_adopted():
    n = get("mt_Fe_vanorman_crispin2010").notes
    assert "as extra interpolation nodes was tested and not adopted" in n
