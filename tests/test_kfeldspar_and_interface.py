"""K-feldspar laws, example condition sources, boundaries, the live Monte Carlo and the layout."""
import os
from pathlib import Path

import numpy as np
import pytest

from diffusor import datasets as ds
from diffusor.coefficients import Conditions, get, list_coefficients
from diffusor.fitting import DiffusionModel, UncertaintyBudget, run_montecarlo
from diffusor.minerals import get_mineral
from diffusor.references import REFERENCES
from diffusor.solvers import Geometry, InitialCondition, dirichlet, zero_flux


# --- K-feldspar ------------------------------------------------------------------
def test_kfeldspar_is_a_mineral_with_sr_ba_ti_and_na_k():
    m = get_mineral("kfeldspar")
    assert set(m.species_keys()) == {"Sr", "Ba", "Ti", "Na-K"}
    assert m.isotropic
    for species in ("Sr", "Ba", "Ti", "Na-K"):
        assert list_coefficients("kfeldspar", species), species


def test_kfeldspar_laws_match_their_abstracts():
    T = 1073.15   # 800 C
    c = Conditions(T_K=T)
    R = 8.314462618
    assert get("kfs_Sr_cherniak1996").D(c) == pytest.approx(8.4 * np.exp(-450e3 / (R * T)), rel=1e-9)
    assert get("kfs_Ba_cherniak2002").D(c) == pytest.approx(0.29 * np.exp(-455e3 / (R * T)), rel=1e-9)
    assert get("kfs_Ti_cherniak_watson2020").D(c) == pytest.approx(
        3.01e-6 * np.exp(-342e3 / (R * T)), rel=1e-9)


def test_ba_is_about_1_7_log_units_slower_than_sr_in_sanidine():
    """Audetat, Grocolas & Mutch (2026) section 2.3 quote this gap for the two Cherniak laws."""
    c = Conditions(T_K=1073.15)
    gap = np.log10(get("kfs_Sr_cherniak1996").D(c) / get("kfs_Ba_cherniak2002").D(c))
    assert 1.6 < gap < 1.8, gap


def test_kfeldspar_uncertainties_come_from_chamberlain():
    assert get("kfs_Sr_cherniak1996").sigma_logD == pytest.approx(0.03)
    assert get("kfs_Ba_cherniak2002").sigma_logD == pytest.approx(0.12)
    for key in ("cherniak1996", "cherniak2002", "cherniak_watson2020", "chamberlain2014"):
        assert key in REFERENCES


def test_only_grocolas_entries_carry_a_covariance():
    with_cov = sorted(c.key for c in list_coefficients() if c.covariance is not None)
    assert with_cov == ["plag_Ba_grocolas2025", "plag_Sr_grocolas2025"]


# --- example conditions ------------------------------------------------------------
def test_every_example_says_where_t_p_and_fo2_come_from():
    for d in ds.DATASETS:
        for key in ("T", "P", "fO2"):
            assert d.sources.get(key), (d.key, key)


def test_kizimen_conditions_follow_supplementary_data_3():
    s = ds.get("opx_kizimen").settings
    assert (s["T_C"], s["sigma_T_K"]) == (850.0, 57.0)
    assert (s["buffer"], s["delta_buffer"], s["sigma_delta"]) == ("NNO", 1.3, 0.35)
    assert (s["P_MPa"], s["sigma_P_MPa"]) == (200.0, 50.0)
    assert "Supplementary Data 3" in ds.get("opx_kizimen").sources["T"]


def test_sanidine_example_recovers_its_time():
    from diffusor.fitting import fit_time
    d = ds.get("sanidine_ba")
    from diffusor.dataio import ProfileSpec, build_profile, read_table
    p = build_profile(read_table(d.path), ProfileSpec(**d.spec))
    ic = InitialCondition("step", {"x0": 10.0, "C_left": 4800.0, "C_right": 1600.0})
    m = DiffusionModel(coefficient=get("kfs_Ba_cherniak2002"), conditions=Conditions(T_K=1063.15),
                       initial=ic, bc_left=dirichlet(4800.0), bc_right=dirichlet(1600.0))
    r = fit_time(m, p.x, p.C, p.sigma, ["t", "x0"])
    years = r.t_seconds / (365.25 * 86400)
    assert 3500 < years < 7000, years


# --- boundaries and solver ------------------------------------------------------------
def _step_model(bc_left=None, bc_right=None, **kw):
    ic = InitialCondition("step", {"x0": 0.0, "C_left": 1.0, "C_right": 0.0})
    return DiffusionModel(coefficient=get("kfs_Sr_cherniak1996"), conditions=Conditions(T_K=1100.0),
                          initial=ic, bc_left=bc_left or dirichlet(1.0),
                          bc_right=bc_right or dirichlet(0.0),
                          x_grid=np.linspace(-50, 50, 201), **kw)


def test_a_crystal_rim_or_centre_switches_to_the_numerical_solver():
    assert _step_model().can_use_analytical()[0]
    ok, why = _step_model(boundaries_far=False).can_use_analytical()
    assert not ok and "rim or centre" in why


def test_far_ends_give_the_same_profile_with_either_solver():
    """While the plateaus survive, fixed far ends reproduce the closed form."""
    x = np.linspace(-50, 50, 101)
    t = 5.0e2 * 365.25 * 86400        # 2 sqrt(Dt) is about 15 um, well inside the traverse
    a = _step_model().profile(t, x)
    n = _step_model(force_numerical=True).profile(t, x)
    assert np.max(np.abs(a - n)) < 5e-3


def test_closed_ends_hold_the_mass_in_and_change_the_profile():
    """Zero flux at both ends is the closed system of Costa et al. (2008) Fig. 6."""
    x = np.linspace(-50, 50, 201)
    t = 2.0e5 * 365.25 * 86400
    closed = _step_model(boundaries_far=False, bc_left=zero_flux(), bc_right=zero_flux())
    C = closed.profile(t, x)
    assert np.trapezoid(C, x) == pytest.approx(50.0, rel=0.02)
    fixed = _step_model(boundaries_far=False).profile(t, x)
    assert np.max(np.abs(C - fixed)) > 0.01


# --- live Monte Carlo --------------------------------------------------------------------
def test_monte_carlo_reports_every_draw():
    x = np.linspace(-20, 20, 41)
    m = _step_model()
    C = m.profile(1.0e4 * 365.25 * 86400, x) + np.random.default_rng(1).normal(0, 0.01, x.size)
    seen = []
    res = run_montecarlo(m, x, C, np.full(x.size, 0.01),
                         budget=UncertaintyBudget(sigma_T_K=20.0), n_draws=12, seed=3,
                         on_draw=seen.append)
    assert len(seen) == res.n_draws == 12
    first = seen[0]
    for key in ("t", "T_K", "log10_D", "x", "C", "C_model"):
        assert key in first
    assert len({round(s["T_K"], 6) for s in seen}) > 1, "temperature was not sampled"
    assert first["C"].shape == x.shape


def test_worker_sends_all_draws_in_batches():
    from diffusor.gui.workers import MonteCarloWorker
    x = np.linspace(-20, 20, 41)
    m = _step_model()
    C = m.profile(1.0e4 * 365.25 * 86400, x)
    w = MonteCarloWorker(m, x, C, np.full(x.size, 0.01), UncertaintyBudget(sigma_T_K=20.0), 15, 1,
                         ["t"], 1e2, 3.2e12)
    batches, done = [], []
    w.draws.connect(batches.append)
    w.finished.connect(done.append)
    w.run()
    assert done and sum(len(b) for b in batches) == 15


def test_the_band_on_the_profile_is_the_spread_of_the_refitted_draws(app):
    """The band must show what the data constrain, not the ends of the time interval.

    Every draw is re-fitted, so all of them pass through the points and the band is
    about the size of the measurement uncertainty. Temperature moves the time, not
    the curve, so drawing the profile at the 16th and 84th percentile times with D
    held fixed would be many times too wide.
    """
    from matplotlib.collections import PolyCollection
    from diffusor.fitting import fit_time
    from diffusor.gui.plot_widget import ProfilePlot
    d = ds.get("sanidine_ba")
    from diffusor.dataio import ProfileSpec, build_profile, read_table
    p = build_profile(read_table(d.path), ProfileSpec(**d.spec))
    ic = InitialCondition("step", {"x0": 10.0, "C_left": 4800.0, "C_right": 1600.0})
    m = DiffusionModel(coefficient=get("kfs_Ba_cherniak2002"), conditions=Conditions(T_K=1063.15),
                       initial=ic, bc_left=dirichlet(4800.0), bc_right=dirichlet(1600.0))
    free = ["t", "x0"]
    r = fit_time(m, p.x, p.C, p.sigma, free)
    mc = run_montecarlo(m, p.x, p.C, p.sigma, budget=UncertaintyBudget(sigma_T_K=30.0),
                        n_draws=40, seed=2, free_parameters=free)
    plot = ProfilePlot()
    plot.show_fit(r, mc)

    def widest(ax):
        out = 0.0
        for coll in ax.collections:
            if not isinstance(coll, PolyCollection):
                continue
            v = coll.get_paths()[0].vertices
            for xv in np.unique(np.round(v[:, 0], 6)):
                ys = v[np.isclose(v[:, 0], xv), 1]
                if ys.size >= 2:
                    out = max(out, float(ys.max() - ys.min()))
        return out

    sigma = float(np.mean(p.sigma))
    drawn = widest(plot.ax)
    assert 0 < drawn < 3 * sigma, (drawn, sigma)
    time_band = float(np.max(np.abs(m.profile(mc.p84, p.x) - m.profile(mc.p16, p.x))))
    assert time_band > 5 * drawn, "the old time band was not this narrow, so this test is blind"
    plot.deleteLater()


# --- interface ----------------------------------------------------------------------------
@pytest.fixture(scope="module")
def app():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    import matplotlib
    matplotlib.use("QtAgg")
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _settle(app, n=20):
    import time
    for _ in range(n):
        app.processEvents()
        time.sleep(0.005)


def test_no_step_page_scrolls_and_every_page_fits_a_laptop(app):
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import QAbstractScrollArea
    from diffusor.gui.main_window import MainWindow
    w = MainWindow()
    w.resize(1366, 690)
    w.show()
    for page in range(6):
        assert not isinstance(w.pages.widget(page), QAbstractScrollArea), page
    for i in range(w.lst_examples.count()):
        w._go(0)
        w.lst_examples.setCurrentRow(i)
        w.load_example()
        for step in range(7):
            w._go(step)
            _settle(app)
        assert w.minimumSizeHint().height() <= 700, (w.lst_examples.item(i).data(Qt.UserRole),
                                                    w.minimumSizeHint().height())
    w.close()


def test_boundary_choices_reach_the_model(app):
    from PySide6.QtCore import Qt
    from diffusor.gui.main_window import MainWindow
    w = MainWindow()
    for i in range(w.lst_examples.count()):
        if w.lst_examples.item(i).data(Qt.UserRole) == "sanidine_ba":
            w.lst_examples.setCurrentRow(i)
    w.load_example()
    m = w._model(w._checked_keys()[0])
    assert m.boundaries_far and m.can_use_analytical()[0]
    w.cmb_bcr.setCurrentIndex(w.cmb_bcr.findData("centre"))
    m = w._model(w._checked_keys()[0])
    assert m.bc_right.kind == "neumann" and not m.can_use_analytical()[0]
    assert "Numerical" in w.lbl_solver.text()
    w.cmb_bcr.setCurrentIndex(w.cmb_bcr.findData("rim_melt"))
    assert w._model(w._checked_keys()[0]).bc_right.kind == "dirichlet"
    w.close()


def test_covariance_sampling_is_offered_only_where_it_exists(app):
    from diffusor.gui.main_window import MainWindow
    w = MainWindow()
    idx = w.cmb_dmode.findData("covariance")
    keys = [w.cmb_mineral.itemData(i) for i in range(w.cmb_mineral.count())]
    w.cmb_mineral.setCurrentIndex(keys.index("opx"))
    w._on_dmode_changed()
    assert not w.cmb_dmode.model().item(idx).isEnabled()
    w.cmb_mineral.setCurrentIndex(keys.index("plagioclase"))
    w.cmb_species.setCurrentText("Sr")
    for i in range(w.lst_coef.count()):
        it = w.lst_coef.item(i)
        from PySide6.QtCore import Qt
        it.setCheckState(Qt.Checked if it.data(Qt.UserRole) == "plag_Sr_grocolas2025" else Qt.Unchecked)
    w._on_dmode_changed()
    assert w.cmb_dmode.model().item(idx).isEnabled()
    w.close()


def test_reading_panes_render_for_every_coefficient_and_example(app):
    from diffusor.gui import richtext
    for c in list_coefficients():
        html = richtext.coefficient_html(c)
        assert "<h2>" in html and c.label.split(",")[0] in html
    for d in ds.DATASETS:
        coef = get(d.settings["coefficient"])
        assert "Conditions it sets" in richtext.example_html(d, coef)
    assert "<table" in richtext.boundaries_html()
    assert "<table" in richtext.format_html()


def test_unused_conditions_are_named():
    from diffusor.gui.richtext import unused_conditions
    assert "pressure or oxygen fugacity" in unused_conditions(get("opx_FeMg_ganguly_tazzoli1994_nofo2"))
    assert unused_conditions(get("opx_FeMg_dias2025")).find("oxygen") == -1


def test_theme_keeps_selection_on_hover_ticks_boxes_and_uses_dark_hover_text():
    from diffusor.gui import theme
    css = theme.stylesheet()
    assert "QListWidget::item:selected:hover" in css
    assert "check.svg" in css and (theme.ICONS / "check.svg").exists()
    primary_hover = css.split("QPushButton#Primary:hover")[1].split("}")[0]
    assert f"color: {theme.TEXT}" in primary_hover
    tooltip = css.split("QToolTip")[1].split("}")[0]
    assert f"color: {theme.TEXT}" in tooltip
    assert theme.ACCENT_SOFT in css.split("QListWidget::item:hover")[1].split("}")[0]
