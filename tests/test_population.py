"""Timescale populations: records, densities, loaders, the figure and the window."""
import os

import numpy as np
import pytest

from diffusor.constants import SEC_PER_DAY, SEC_PER_YEAR
from diffusor.population import (Timescale, describe, extend_unique, group_density, load_result_json,
                                 load_table, log_grid, save_population, summary)


def _integral(d, grid):
    return float(np.sum((d[1:] + d[:-1]) / 2 * np.diff(grid)))


@pytest.fixture(scope="module")
def run():
    from diffusor.coefficients import Conditions, get
    from diffusor.fitting import DiffusionModel, UncertaintyBudget, fit_time, run_montecarlo
    from diffusor.solvers import Geometry, InitialCondition, dirichlet
    cond = Conditions(T_K=1223.15, P_Pa=2e8, log_fo2_bar=-10.0, X={"XFe": 0.2}, axis="c")
    ic = InitialCondition("step", {"x0": 0.0, "C_left": 0.30, "C_right": 0.18})
    m = DiffusionModel(coefficient=get("opx_FeMg_dohmen2016"), conditions=cond, initial=ic,
                       geometry=Geometry("plane"), bc_left=dirichlet(0.30),
                       bc_right=dirichlet(0.18), n_nodes=151)
    x = np.linspace(-30, 30, 31)
    r = fit_time(m, x, m.profile(SEC_PER_YEAR, x))
    mc = run_montecarlo(m, x, m.profile(SEC_PER_YEAR, x), np.full_like(x, 0.003), n_draws=30,
                        seed=1, budget=UncertaintyBudget(sigma_T_K=20.0, buffer="NNO",
                                                         delta_buffer=1.0))
    return r, mc


# --- records and densities ---------------------------------------------------------
def test_every_kind_of_density_integrates_to_one_per_log10_unit():
    rng = np.random.default_rng(0)
    kde = Timescale("a", 30 * SEC_PER_DAY, samples=30 * SEC_PER_DAY * 10 ** rng.normal(0, .2, 500))
    split = Timescale("b", 30 * SEC_PER_DAY, p16=10 * SEC_PER_DAY, p84=200 * SEC_PER_DAY)
    grid = log_grid([kde, split], n=20000)       # fine enough for the step at the median
    for r in (kde, split):
        assert _integral(r.log_density(grid), grid) == pytest.approx(1, abs=1e-3)
    assert Timescale("c", 5.0).log_density(grid) is None


def test_half_normals_reach_the_16th_and_84th_percentiles():
    r = Timescale("b", 30.0, median=30.0, p16=10.0, p84=200.0)
    grid = np.linspace(-2, 5, 20001)
    cdf = np.cumsum(r.log_density(grid)) * (grid[1] - grid[0])
    assert np.interp(np.log10(10.0), grid, cdf) == pytest.approx(.1587, abs=2e-3)
    assert np.interp(np.log10(200.0), grid, cdf) == pytest.approx(.8413, abs=2e-3)


def test_monte_carlo_times_set_the_percentiles():
    s = np.geomspace(1, 1000, 999)
    r = Timescale("a", 40.0, p16=1.0, p84=2.0, samples=s)
    assert r.p16 == pytest.approx(np.percentile(s, 16))
    assert r.median == pytest.approx(np.percentile(s, 50))
    assert "kernel density" in r.density_basis


def test_bad_times_are_refused():
    for t in (0, -1, np.nan):
        with pytest.raises(ValueError):
            Timescale("x", t)


def test_each_crystal_weighs_the_same_in_its_group_density():
    narrow = Timescale("n", 10.0, p16=9.0, p84=11.0)
    wide = Timescale("w", 1e4, p16=1e3, p84=1e5)
    grid = log_grid([narrow, wide], n=4000)
    d, n = group_density([narrow, wide], grid)
    assert n == 2
    # all of the narrow crystal lies below 10^2.5 s, and 6.7% of the wide one (1.5 sigma)
    half = grid < 2.5
    assert _integral(d[half], grid[half]) == pytest.approx(.5 + .5 * .0668, abs=2e-3)


def test_summary_and_methods_count_groups_and_density_kinds():
    recs = [Timescale(f"a{i}", 10.0 * (i + 1), group="A", p16=5.0, p84=40.0) for i in range(3)]
    recs.append(Timescale("b", 1000.0, group="B"))
    rows = {r["group"]: r for r in summary(recs)}
    assert rows["A"]["n"] == 3 and rows["A"]["median_s"] == 20.0
    assert rows["B"]["n_with_uncertainty"] == 0
    text = describe(recs)
    assert "0 densities are Gaussian kernel" in text and "3 are two half-normals" in text
    assert "1 records without uncertainty" in text


def test_clashing_labels_are_renamed():
    recs = [Timescale("ol1", 1.0)]
    added = extend_unique(recs, [Timescale("ol1", 2.0), Timescale("ol1", 3.0)])
    assert added == ["ol1 (2)", "ol1 (3)"]


# --- loaders -----------------------------------------------------------------------
def test_a_table_with_units_in_the_headers(tmp_path):
    p = tmp_path / "t.csv"
    p.write_text("Crystal,Eruption,Time (days),p16_d,p84_d\n"
                 "ol1,2011,30,20,45\nol2,2018,2,,\nbad,2018,n/a,,\n", encoding="utf-8")
    recs = load_table(p, default_unit="yr")
    assert [r.label for r in recs] == ["ol1", "ol2"]
    assert [r.group for r in recs] == ["2011", "2018"]
    assert recs[0].t_best == pytest.approx(30 * SEC_PER_DAY)
    assert recs[0].p84 == pytest.approx(45 * SEC_PER_DAY)
    assert not recs[1].has_interval


def test_a_table_with_a_unit_column_and_a_log_sigma(tmp_path):
    p = tmp_path / "t.csv"
    p.write_text("name;t;unit;sigma_log10\nopx1;1.5;yr;0.2\nopx2;40;days;0.1\n", encoding="utf-8")
    recs = load_table(p, group="Kirishima")
    assert recs[0].t_best == pytest.approx(1.5 * SEC_PER_YEAR)
    assert recs[1].t_best == pytest.approx(40 * SEC_PER_DAY)
    assert recs[0].p84 / recs[0].t_best == pytest.approx(10 ** 0.2)
    assert {r.group for r in recs} == {"Kirishima"}


def test_a_table_without_a_time_column_is_refused(tmp_path):
    p = tmp_path / "t.csv"
    p.write_text("label,value\na,1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="no time column"):
        load_table(p)


def test_exported_results_load_with_their_monte_carlo_times(run, tmp_path):
    from diffusor.dataio import save_results
    r, mc = run
    save_results(tmp_path, r, mc, include=["json", "montecarlo"], basename="ol7")
    rec = load_result_json(tmp_path / "ol7_results.json")
    assert rec.label == "ol7"
    assert rec.t_best == pytest.approx(r.t_seconds)
    assert rec.samples.size == mc.times.size
    assert rec.median == pytest.approx(mc.median)
    assert rec.coefficient == "opx_FeMg_dohmen2016"
    # without the times file the percentiles of the JSON are used
    (tmp_path / "ol7_montecarlo_times.csv").unlink()
    rec = load_result_json(tmp_path / "ol7_results.json")
    assert rec.samples is None and rec.p16 == pytest.approx(mc.p16)


def test_a_shared_duration_file_is_not_taken_for_one_crystal(tmp_path):
    p = tmp_path / "joint_results.json"
    p.write_text('{"model_family": "x", "t_seconds": 5, "profiles": []}', encoding="utf-8")
    with pytest.raises(ValueError, match="shared-duration"):
        load_result_json(p)


# --- figure and export -------------------------------------------------------------
def _population():
    rng = np.random.default_rng(3)
    recs = [Timescale(f"a{i}", 10 ** rng.normal(6.5, .3), group="A",
                      samples=10 ** rng.normal(6.5, .3, 50)) for i in range(5)]
    recs += [Timescale(f"b{i}", 10 ** rng.normal(7.5, .3), group="B", p16=1e7, p84=1e8) for i in range(4)]
    recs.append(Timescale("c", 2e7, group="B"))
    return recs


def test_the_figure_draws_the_chosen_panels_on_one_time_axis():
    from diffusor.dataio.figures import population_figure
    fig = population_figure(_population())
    axes = fig.get_axes()
    assert len(axes) == 3
    assert all(ax.get_xscale() == "log" for ax in axes)
    assert axes[-1].get_xlabel() == "Time (days)"
    assert len(population_figure(_population(), panels=["density"], unit="yr").get_axes()) == 1


def test_more_groups_than_colours_are_refused():
    from diffusor.dataio.figures import population_figure
    recs = [Timescale(str(i), 10.0 + i, group=str(i)) for i in range(9)]
    with pytest.raises(ValueError, match="at most 8"):
        population_figure(recs)


def test_export_writes_tables_methods_and_figures(tmp_path):
    import pandas as pd
    from diffusor.dataio.figures import population_figure
    recs = _population()
    written = save_population(tmp_path, recs, population_figure(recs))
    assert set(written) == {"timescales", "groups", "density", "methods", "figure",
                            "figure_svg", "figure_pdf"}
    df = pd.read_csv(written["timescales"])
    assert len(df) == len(recs)
    dens = pd.read_csv(written["density"])
    assert {"density_A", "density_B"} <= set(dens.columns)
    assert "Silverman" in (tmp_path / "population_methods.txt").read_text(encoding="utf-8")


# --- the window --------------------------------------------------------------------
@pytest.fixture(scope="module")
def app():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def test_the_window_adds_the_current_run_and_follows_edits(app, run):
    from types import SimpleNamespace
    from PySide6.QtWidgets import QWidget
    from diffusor.gui.population_study import PopulationStudyDialog
    r, mc = run
    owner = QWidget()
    owner.fit_result, owner.mc_result, owner._stale = r, mc, False
    owner.profile = SimpleNamespace(source="C:/data/opx_12.csv")
    owner.output_dir = lambda: ""
    dlg = PopulationStudyDialog(owner)
    dlg.group_edit.setText("2011")
    dlg.add_current()
    dlg.add_current()
    assert [x.label for x in dlg.records] == ["opx_12", "opx_12 (2)"]
    assert dlg.records[0].group == "2011" and dlg.records[0].samples is not None
    assert dlg.figure is not None and dlg.export.isEnabled()
    dlg.table.item(1, 1).setText("2018")
    assert dlg.records[1].group == "2018"
    dlg.table.item(1, 0).setText("opx_12")          # a repeated label is put back
    assert dlg.records[1].label == "opx_12 (2)"
    dlg.panel_boxes["ranked"].setChecked(False)
    assert len(dlg.figure.get_axes()) == 2
    dlg.table.selectRow(0)
    dlg.remove_selected()
    assert [x.label for x in dlg.records] == ["opx_12 (2)"]
    dlg.clear_all()
    assert dlg.figure is None and not dlg.export.isEnabled()


def test_the_file_menu_opens_the_population_window(app):
    from diffusor.gui.main_window import MainWindow
    w = MainWindow()
    names = [a.text() for m in w.menuBar().actions() for a in m.menu().actions()]
    assert "Timescale population..." in names
    w.show_population_study()
    assert w.population_study.isVisible()
