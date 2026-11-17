"""Regressions from the November 2026 audit: table loading, the Monte Carlo distance
scale, the temperature of a time step along a cooling path, and export names."""
import os
import warnings

import numpy as np
import pandas as pd
import pytest

from diffusor.dataio.profiles import read_table, suggest_spec

ROWS = [(0.5, 10.2, 30.1), (1.5, 10.4, 30.0), (2.5, 11.0, 29.5), (3.5, 11.6, 29.1)]


def _check(df, distance="Distance (µm)"):
    assert list(df.columns)[0] == distance
    assert all(pd.api.types.is_numeric_dtype(df[c]) for c in list(df.columns)[:3])
    assert np.allclose(df.iloc[:, :3].to_numpy(dtype=float), ROWS)


# --- tables as spreadsheet programs write them -------------------------------------
def test_semicolons_and_decimal_commas(tmp_path):
    p = tmp_path / "eu.csv"
    p.write_text("Distance (µm);FeO;MgO;Comment\n" + "".join(
        f"{x};{a};{b};rim\n".replace(".", ",") for x, a, b in ROWS), encoding="utf-8")
    df = read_table(p)
    _check(df)
    assert not pd.api.types.is_numeric_dtype(df["Comment"])
    spec = suggest_spec(df)
    assert (spec.column_a, spec.column_b, spec.distance_unit) == ("FeO", "MgO", "um")


@pytest.mark.parametrize("encoding", ["cp1252", "utf-8-sig", "utf-16"])
def test_encodings_excel_writes(tmp_path, encoding):
    p = tmp_path / "t.csv"
    text = "Distance (µm),FeO,MgO\n" + "".join(f"{x},{a},{b}\n" for x, a, b in ROWS)
    p.write_bytes(text.encode(encoding))
    _check(read_table(p))


def test_a_text_column_with_commas_stays_text(tmp_path):
    p = tmp_path / "t.csv"
    p.write_text('x,name,val\n1,"a,b",3.5\n2,"1,5 mm",3.1\n', encoding="utf-8")
    df = read_table(p)
    assert not pd.api.types.is_numeric_dtype(df["name"])
    assert df["val"].tolist() == [3.5, 3.1]


def test_a_few_text_cells_do_not_make_a_column_text(tmp_path):
    """b.d. in one cell: the column is numbers, the cell is missing, names stay text."""
    cells = [30.1, "b.d.", 29.5, 29.0, 28.4, 28.0]
    frame = pd.DataFrame({"x": np.arange(6.0), "MgO": cells, "Sample": list("abcdef")})
    frame.to_csv(tmp_path / "t.csv", index=False)
    frame.to_excel(tmp_path / "t.xlsx", index=False)
    for name in ("t.csv", "t.xlsx"):
        df = read_table(tmp_path / name)
        assert pd.api.types.is_numeric_dtype(df["MgO"])
        assert np.isnan(df["MgO"][1]) and df["MgO"][2] == 29.5
        assert not pd.api.types.is_numeric_dtype(df["Sample"])


def test_an_uncertainty_column_without_numbers_gives_no_uncertainties():
    """The spread of the whole profile used to be invented as every point's sigma."""
    from diffusor.dataio.profiles import ProfileSpec, build_profile
    df = pd.DataFrame({"x": np.arange(6.0), "C": np.linspace(0.3, 0.2, 6), "err": [np.nan] * 6})
    p = build_profile(df, ProfileSpec("x", "C", sigma_a_column="err", mode="A"))
    assert p.sigma is None and any("weighted equally" in n for n in p.notes)
    df["err"] = [0.01, np.nan, 0.03, 0.0, 0.01, 0.02]
    p = build_profile(df, ProfileSpec("x", "C", sigma_a_column="err", mode="A"))
    assert np.allclose(p.sigma, [0.01, 0.015, 0.03, 0.015, 0.01, 0.02])


def test_population_labels_that_are_mostly_numbers_stay_text(tmp_path):
    from diffusor.population import load_table
    p = tmp_path / "times.csv"
    p.write_text("crystal;time (days)\n12;3,5\n13;4,0\n14;2,5\n15;6,0\n16b;5,0\n", encoding="utf-8")
    recs = load_table(p)
    assert [r.label for r in recs] == ["12", "13", "14", "15", "16b"]
    assert recs[0].t_best == pytest.approx(3.5 * 86400)


def test_a_tab_separated_map_with_decimal_commas(tmp_path):
    from diffusor.dataio.images import load_text_grid
    p = tmp_path / "map.txt"
    p.write_text("1,5\t2,5\t3,5\n4,5\t5,5\t6,5\n7,5\t8,5\t9,5\n", encoding="utf-8")
    img = load_text_grid(p)
    assert img.data.shape == (3, 3) and img.data[1, 1] == 5.5


def test_calibration_scatter_uses_the_degrees_of_freedom_of_the_fit():
    from diffusor.dataio.greyscale import calibrate
    g = np.array([100.0, 150.0, 200.0])
    cal = calibrate(g, [0.10, 0.21, 0.30], 1)
    s2 = float(np.sum(cal.residuals ** 2)) / (3 - 2)
    powers = np.vstack([g, np.ones(3)])
    var = np.einsum("ij,ik,jk->k", cal.covariance, powers, powers)
    assert np.allclose(cal.sigma_from_calibration(g), np.sqrt(var + s2))
    assert np.all(cal.sigma_from_calibration(g) > cal.rmse)


def test_multicomponent_cooling_takes_the_temperature_at_mid_step():
    """A scalar problem written as a 1x1 matrix must match the diffusion integral."""
    from scipy.optimize import minimize_scalar
    from diffusor.solvers import analytical
    from diffusor.solvers.history import ThermalHistory, effective_Dt
    from diffusor.solvers.initial import step
    from diffusor.solvers.multicomponent import solve_multicomponent

    def D_of_T(T):
        return 1e-6 * np.exp(-300e3 / (8.314462618 * np.asarray(T, dtype=float))) * 1e12

    x = np.linspace(0, 200, 401)
    C0 = step(x, 100, 1, 0)[None, :]
    t = 3e7
    h = ThermalHistory.linear(1273.15, 973.15, t)
    Dt = effective_Dt(D_of_T, h, t, n=200001)
    C = solve_multicomponent(x, C0, lambda C, T: np.full((1, 1, x.size), float(D_of_T(T))), t,
                             history=h).C_final[0]
    best = minimize_scalar(
        lambda lg: np.sum((analytical.step_infinite(x, 100, 1, 0, 10 ** lg) - C) ** 2),
        bounds=(np.log10(Dt) - 1, np.log10(Dt) + 1), method="bounded")
    assert 10 ** best.x / Dt == pytest.approx(1.0, abs=0.005)


def test_the_step_may_start_outside_the_measured_range():
    """A start outside the bounds stopped least_squares; the scan minimum was returned."""
    from diffusor import crystal_setup as cs
    from diffusor.dataio.profiles import ProfileSpec
    from diffusor.fitting import fit_time
    from diffusor.solvers import analytical
    x = np.linspace(5, 200, 40)
    s = cs.complete({"mineral": "opx", "species": "Fe-Mg", "coefficients": ["opx_FeMg_dohmen2016"],
                     "T_C": 950, "x0_um": 0.0, "C_left": 0.3, "C_right": 0.2,
                     "composition_dependent": False, "x_composition": 0.25})
    data = cs.CrystalData(pd.DataFrame({"x": x, "C": 0 * x}), ProfileSpec("x", "C", mode="A"))
    model = cs.build_model(s, data)
    t_true = 400.0 / (model.D_bulk() * 1e12)
    C = analytical.step_infinite(x, 0.0, 0.3, 0.2, 400.0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = fit_time(model, x, C, None, ["t", "x0"])
    assert r.success and r.t_seconds == pytest.approx(t_true, rel=1e-3)
    assert not any("refinement" in w for w in r.warnings)


# --- interface ---------------------------------------------------------------------
@pytest.fixture(scope="module")
def app():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    import matplotlib
    matplotlib.use("QtAgg")
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def test_column_dialog_takes_a_table_with_a_text_column(app):
    """numpy's issubdtype raises on the string dtype of pandas 3."""
    from diffusor.gui.main_window import ColumnDialog
    df = pd.DataFrame({"x": [0.0, 1.0, 2.0], "FeO": [1.0, 2.0, 3.0], "MgO": [3.0, 2.0, 1.0],
                       "Sample": ["a", "b", "c"]})
    dlg = ColumnDialog(df, suggest_spec(df))
    names = [dlg.a.itemText(i) for i in range(dlg.a.count())]
    assert names == ["x", "FeO", "MgO"]
    dlg.xmin.setText("0,5")
    dlg.xmax.setText("not a number")
    spec = dlg.spec()
    assert spec.x_min == 0.5 and spec.x_max is None
    dlg.close()


def test_too_few_rows_keep_the_loaded_profile(app, tmp_path, monkeypatch):
    from PySide6.QtWidgets import QMessageBox
    from diffusor.dataio.profiles import ProfileSpec
    from diffusor.gui.main_window import MainWindow
    monkeypatch.setattr(QMessageBox, "warning", lambda *a, **k: None)
    w = MainWindow()
    good = pd.DataFrame({"x": np.arange(10.0), "C": np.linspace(0.3, 0.2, 10)})
    spec = ProfileSpec("x", "C", mode="A")
    assert w._use_table(good, spec, tmp_path / "good.csv", None, (None, True))
    before = w.profile
    bad = pd.DataFrame({"x": [0.0, 1.0, 2.0], "C": [0.3, np.nan, np.nan]})
    assert not w._use_table(bad, spec, tmp_path / "bad.csv", None, (None, True))
    assert w.profile is before
    w.close()


def test_export_is_named_after_the_profile(app, tmp_path, monkeypatch):
    from PySide6.QtWidgets import QDialog, QFileDialog, QMessageBox
    from diffusor.fitting import fit_time
    from diffusor.gui import main_window as mw
    w = mw.MainWindow()
    x = np.linspace(0, 100, 30)
    df = pd.DataFrame({"x": x, "C": 0.25 + 0.05 * np.tanh((x - 50) / 8)})
    df.to_csv(tmp_path / "ol7.csv", index=False)
    from diffusor.dataio.profiles import ProfileSpec
    assert w._use_table(df, ProfileSpec("x", "C", mode="A"), tmp_path / "ol7.csv", None, (None, True))
    p = w._fit_profile()
    w._fit_done(fit_time(w._model(w._checked_keys()[0]), p.x, p.C, p.sigma, w._free_parameters()))
    monkeypatch.setattr(mw.ExportDialog, "exec", lambda self: QDialog.Accepted)
    monkeypatch.setattr(mw.ExportDialog, "selected", lambda self: ["json", "methods"])
    monkeypatch.setattr(QFileDialog, "getExistingDirectory", lambda *a, **k: str(tmp_path / "out"))
    monkeypatch.setattr(QMessageBox, "information", lambda *a, **k: None)
    monkeypatch.setattr(QMessageBox, "warning", lambda *a, **k: None)
    w.export_results()
    assert sorted(f.name for f in (tmp_path / "out").iterdir()) == [
        "ol7_diffusor_methods.txt", "ol7_diffusor_results.json"]
    w.close()


# --- numerics ----------------------------------------------------------------------
def test_cooling_path_matches_the_diffusion_integral():
    """T at the middle of each step: the summed D dt equals the integral to 0.5 per cent."""
    from scipy.optimize import minimize_scalar
    from diffusor.solvers import analytical
    from diffusor.solvers.history import ThermalHistory, effective_Dt
    from diffusor.solvers.initial import step
    from diffusor.solvers.numerical import solve_1d

    def D_of_T(T):
        return 1e-6 * np.exp(-300e3 / (8.314462618 * np.asarray(T, dtype=float))) * 1e12

    x = np.linspace(0, 200, 401)
    C0 = step(x, 100, 1, 0)
    for T0, T1, t in ((1273.15, 1173.15, 3e7), (1273.15, 973.15, 3e7), (1373.15, 773.15, 3e7)):
        h = ThermalHistory.linear(T0, T1, t)
        Dt = effective_Dt(D_of_T, h, t, n=200001)
        C = solve_1d(x, C0, lambda C, T: D_of_T(T), t, history=h).C_final
        best = minimize_scalar(
            lambda lg: np.sum((analytical.step_infinite(x, 100, 1, 0, 10 ** lg) - C) ** 2),
            bounds=(np.log10(Dt) - 1, np.log10(Dt) + 1), method="bounded")
        assert 10 ** best.x / Dt == pytest.approx(1.0, abs=0.005)


def _opx_case(fit_x0):
    from diffusor import crystal_setup as cs
    from diffusor.dataio.profiles import ProfileSpec
    from diffusor.solvers import analytical
    x = np.linspace(0, 200, 41)
    s = cs.complete({"mineral": "opx", "species": "Fe-Mg", "coefficients": ["opx_FeMg_dohmen2016"],
                     "T_C": 950, "x0_um": 100.0, "C_left": 0.3, "C_right": 0.2,
                     "composition_dependent": False, "x_composition": 0.25, "fit_x0": fit_x0})
    C = analytical.step_infinite(x, 100.0, 0.3, 0.2, 25.0)
    data = cs.CrystalData(pd.DataFrame({"x": x, "C": C}), ProfileSpec("x", "C", mode="A"))
    return s, data, cs.build_model(s, data), x, C


@pytest.mark.parametrize("fit_x0", [False, True])
def test_distance_scale_moves_the_model_with_the_data(fit_x0):
    """t scales with length squared, so sigma(log10 t) = 2 sigma_scale / ln 10, whether or
    not the step position is fitted."""
    from diffusor import crystal_setup as cs
    from diffusor.fitting import UncertaintyBudget, run_montecarlo
    s, data, model, x, C = _opx_case(fit_x0)
    budget = UncertaintyBudget(buffer=s["buffer"], delta_buffer=s["fo2_value"],
                               sample_coefficient=False, sample_measurement_noise=False,
                               sigma_distance_scale=0.03)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = run_montecarlo(model, x, C, None, budget=budget, n_draws=150, seed=4,
                           free_parameters=cs.free_parameters(s, data))
    assert r.sigma_log10 == pytest.approx(2 * 0.03 / np.log(10), rel=0.2)
    assert r.median / r.t_best == pytest.approx(1.0, abs=0.02)


def test_no_budget_keeps_the_fo2_of_the_model():
    from diffusor.fitting import run_montecarlo
    s, data, model, x, C = _opx_case(True)
    assert model.fo2_buffer == ("NNO", 1.0)
    seen = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run_montecarlo(model, x, C, None, n_draws=3, seed=1, free_parameters=("t",),
                       on_draw=lambda d: seen.append(d["log_fo2_bar"]))
    assert seen and all(v == pytest.approx(model.conditions.log_fo2_bar) for v in seen)
