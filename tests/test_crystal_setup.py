"""Crystal setups and batches: settings in and out of the window, files, shared values, fitting."""
import json
import os
import time

import numpy as np
import pandas as pd
import pytest

from diffusor import crystal_setup as cs
from diffusor.constants import SEC_PER_YEAR
from diffusor.dataio.profiles import ProfileSpec


def _synthetic(t_years=1.0, x0=0.0, left=0.30, right=0.18, seed=0, extra=None):
    """An opx Fe-Mg step profile made by the forward model, as a table."""
    from diffusor.coefficients import Conditions, get
    from diffusor.fitting import DiffusionModel
    from diffusor.solvers import Geometry, InitialCondition, dirichlet
    cond = Conditions(T_K=1223.15, P_Pa=2e8, log_fo2_bar=-10.0, X={"XFe": 0.2}, axis="c")
    ic = InitialCondition("step", {"x0": x0, "C_left": left, "C_right": right})
    m = DiffusionModel(coefficient=get("opx_FeMg_dohmen2016"), conditions=cond, initial=ic,
                       geometry=Geometry("plane"), bc_left=dirichlet(left),
                       bc_right=dirichlet(right), n_nodes=151)
    x = np.linspace(-30, 30, 31)
    C = m.profile(t_years * SEC_PER_YEAR, x) + np.random.default_rng(seed).normal(0, 0.002, x.size)
    df = pd.DataFrame({"x_um": x, "XFe": C, "sXFe": np.full(x.size, 0.002)})
    if extra:
        for k, v in extra.items():
            df[k] = v
    return df


SPEC = ProfileSpec(distance_column="x_um", column_a="XFe", sigma_a_column="sXFe", mode="A")


def _settings(**kw):
    s = cs.complete({"mineral": "opx", "species": "Fe-Mg", "coefficients": ["opx_FeMg_dohmen2016"],
                     "T_C": 950.0, "P_MPa": 200.0, "fo2_mode": "absolute", "fo2_value": -10.0,
                     "x_composition": 0.2, "x_composition_source": "explicit", "axis": "c",
                     "composition_dependent": False, "fit_x0": True, "fit_plateaus": True})
    s.update(kw)
    return s


# --- the pure module ---------------------------------------------------------------
def test_tables_survive_json_with_missing_values_and_text():
    df = pd.DataFrame({"x": [0.0, 1.5, 2.0], "C": [0.1, np.nan, 0.3], "note": ["a", None, "c"]})
    back = cs.table_from_json(json.loads(json.dumps(cs.table_to_json(df))))
    assert back["x"].tolist() == [0.0, 1.5, 2.0]
    assert np.isnan(back["C"][1]) and back["C"][2] == 0.3
    assert back["note"].tolist()[0] == "a" and pd.isna(back["note"].tolist()[1])


def test_complete_fills_defaults_and_drops_unknown_keys():
    s = cs.complete({"T_C": 800.0, "not_a_setting": 1})
    assert s["T_C"] == 800.0 and "not_a_setting" not in s
    assert set(s) == set(cs.DEFAULTS)


def test_differences_names_only_what_differs():
    base = _settings()
    own = dict(base, T_C=900.0, angles_deg=[10.0, 80.0, 90.0])
    assert cs.differences(own, base) == {"T_C": 900.0, "angles_deg": [10.0, 80.0, 90.0]}
    assert cs.differences(dict(base, T_C=950.0 + 1e-14), base) == {}


def test_a_crystal_is_fitted_with_its_own_step_and_plateaus():
    data = cs.CrystalData(_synthetic(t_years=2.0, x0=4.0, left=0.32, right=0.16), SPEC, "ol1.csv")
    s = cs.crystal_settings(_settings(), data)
    assert s["x0_um"] == pytest.approx(4.0, abs=3.0)
    assert s["C_left"] == pytest.approx(0.32, abs=0.01)
    fit, mc = cs.run_crystal(s, data)
    assert mc is None
    assert fit.t_seconds / SEC_PER_YEAR == pytest.approx(2.0, rel=0.1)


def test_points_cut_by_hand_are_kept_out_of_the_fit():
    df = _synthetic()
    df.loc[3, "XFe"] = 0.9                          # an outlier
    data = cs.CrystalData(df, SPEC, "ol.csv")
    mask = [False] * len(df)
    mask[3] = True                                  # the profile is sorted by x, as is the table
    data.excluded = mask
    assert len(data.profile().kept()) == len(df) - 1
    fit, _ = cs.run_crystal(_settings(fit_plateaus=False, C_left=0.30, C_right=0.18), data)
    assert fit.t_seconds / SEC_PER_YEAR == pytest.approx(1.0, rel=0.1)


def test_setup_files_hold_data_settings_and_result(tmp_path):
    data = cs.CrystalData(_synthetic(), SPEC, str(tmp_path / "opx7.csv"))
    s = cs.crystal_settings(_settings(T_C=940.0), data)
    fit, mc = cs.run_crystal(dict(s, mc_draws=20, mc_seed=3), data, monte_carlo=True)
    path = cs.save_setup(tmp_path / "opx7", s, data, cs.result_to_json(fit, mc), group="2011")
    assert path.suffix == ".diffusor"
    c = cs.load_setup(path)
    assert c.label == "opx7" and c.group == "2011"
    assert c.settings == cs.complete(s)
    assert c.data.table.equals(data.table) and c.data.spec == data.spec
    assert c.result["t_seconds"] == pytest.approx(fit.t_seconds)
    assert len(c.result["montecarlo"]["times_s"]) == mc.times.size


def test_files_from_a_newer_version_or_of_another_kind_are_refused(tmp_path):
    p = tmp_path / "x.diffusor"
    p.write_text(json.dumps({"format": cs.FORMAT_SETUP, "version": 99}), encoding="utf-8")
    with pytest.raises(ValueError, match="newer Diffusor"):
        cs.read_file(p)
    p.write_text("{}", encoding="utf-8")
    with pytest.raises(ValueError, match="not a Diffusor"):
        cs.read_file(p)


# --- a batch -----------------------------------------------------------------------
def _batch():
    b = cs.Batch(_settings())
    for k, (x0, left) in enumerate(((0.0, 0.30), (5.0, 0.33), (-4.0, 0.28))):
        b.add(f"opx{k}", cs.CrystalData(_synthetic(x0=x0, left=left, seed=k), SPEC, f"opx{k}.csv"))
    return b


def test_each_crystal_keeps_its_own_plateaus_and_shares_the_rest():
    b = _batch()
    assert {"x0_um", "C_left", "C_right"} <= set(b.crystals[1].overrides)
    assert "T_C" not in b.crystals[1].overrides
    assert b.settings(1)["C_left"] == pytest.approx(0.33, abs=0.01)
    b.set_shared("T_C", 900.0)
    assert all(b.settings(i)["T_C"] == 900.0 for i in range(3))


def test_an_own_value_survives_a_shared_change_only_when_asked():
    b = _batch()
    b.set_own(0, "T_C", 980.0)
    assert b.is_own(0, "T_C") and b.own_count("T_C") == 1
    b.set_shared("T_C", 900.0, keep_own=True)
    assert b.settings(0)["T_C"] == 980.0 and b.settings(1)["T_C"] == 900.0
    b.set_shared("T_C", 910.0)
    assert b.settings(0)["T_C"] == 910.0 and not b.is_own(0, "T_C")
    b.set_own(2, "T_C", 910.0)                       # equal to the shared value: not its own
    assert not b.is_own(2, "T_C")


def test_a_changed_setting_clears_the_result():
    b = _batch()
    b.crystals[0].result = {"t_seconds": 1.0}
    b.set_own(0, "P_MPa", 300.0)
    assert b.crystals[0].result is None and "fit again" in b.crystals[0].status


def test_a_batch_file_round_trips(tmp_path):
    b = _batch()
    b.set_own(2, "angles_deg", [30.0, 70.0, 70.0])
    b.crystals[1].group = "2018"
    back = cs.Batch.from_json(cs.read_file(b.save(tmp_path / "batch")))
    assert [c.label for c in back.crystals] == ["opx0", "opx1", "opx2"]
    assert back.crystals[1].group == "2018"
    for i in range(3):
        assert back.settings(i) == b.settings(i)
    with pytest.raises(ValueError, match="batch"):
        cs.load_setup(tmp_path / "batch.diffusor")


def test_population_reads_fitted_crystals_from_setup_files(tmp_path):
    from diffusor.population import load_setup_file
    b = _batch()
    for i in range(2):
        fit, _ = cs.run_crystal(b.settings(i), b.crystals[i].data)
        b.crystals[i].result = cs.result_to_json(fit)
    b.crystals[0].group = "2011"
    recs = load_setup_file(b.save(tmp_path / "b"), group="other")
    assert [(r.label, r.group) for r in recs] == [("opx0", "2011"), ("opx1", "other")]


# --- the window --------------------------------------------------------------------
@pytest.fixture(scope="module")
def app():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _window(app):
    from diffusor.gui.main_window import MainWindow
    w = MainWindow()
    w.show()
    app.processEvents()
    return w


def _load(w, df, path="opx_test.csv"):
    from pathlib import Path
    assert w._use_table(df, SPEC, Path(path), None, (None, True))


@pytest.mark.parametrize("example", ["opx_shinmoedake", "plag_mg_synthetic"])
def test_the_window_gives_back_the_settings_it_was_given(app, example):
    from diffusor import datasets as ds
    if example not in ds.BY_KEY:
        pytest.skip(f"{example} is not bundled")
    d = ds.get(example)
    w = _window(app)
    w._load(d.path, dataset=d)
    w.sp_T.setValue(912.0)
    w.cmb_axis.setCurrentIndex(4)
    w.sp_beta.setValue(60.0)
    w.sp_gamma.setValue(30.0)
    s = w.current_settings()
    w2 = _window(app)
    w2._load(d.path, dataset=d)
    w2.sp_T.setValue(700.0)
    assert w2.apply_settings(s) == []
    assert w2.current_settings() == s
    t = 3.0e6
    x = w.profile.x
    m1, m2 = w._model(s["coefficients"][0]), w2._model(s["coefficients"][0])
    assert np.allclose(m1.profile(t, x), m2.profile(t, x))
    w.close(); w2.close()


def test_a_saved_setup_reopens_with_its_data_cut_points_and_settings(app, tmp_path, monkeypatch):
    from PySide6.QtWidgets import QFileDialog, QMessageBox
    df = _synthetic(t_years=1.5)
    w = _window(app)
    _load(w, df)
    w.apply_settings(cs.crystal_settings(_settings(T_C=930.0), w.current_data()))
    mask = np.zeros(len(w.profile), dtype=bool)
    mask[[0, 5]] = True
    w._on_points_cut(mask)
    target = tmp_path / "crystal.diffusor"
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: (str(target), ""))
    w.save_setup()
    assert target.exists()
    saved = w.current_settings()
    w.close()

    w2 = _window(app)
    monkeypatch.setattr(QMessageBox, "warning", lambda *a, **k: None)
    w2.open_setup(str(target))
    app.processEvents()
    assert w2.current_settings() == saved
    assert w2.profile.excluded is not None and w2.profile.excluded.tolist() == mask.tolist()
    assert len(w2._fit_profile()) == len(df) - 2
    w2.close()


def test_the_main_window_has_setup_and_batch_menus(app):
    w = _window(app)
    names = [a.text() for m in w.menuBar().actions() for a in m.menu().actions()]
    for name in ("Open crystal setup...", "Save crystal setup...", "Batch of crystals..."):
        assert name in names
    w.close()


def _wait(app, cond, seconds=120):
    end = time.time() + seconds
    while not cond() and time.time() < end:
        app.processEvents()
        time.sleep(0.02)
    assert cond()


def test_the_batch_window_shares_values_fits_and_feeds_the_population(app, monkeypatch):
    from PySide6.QtWidgets import QMessageBox
    from diffusor.gui.batch_study import KEYS, SHARED_ROW, BatchStudyDialog
    w = _window(app)
    _load(w, _synthetic(seed=1))
    w.apply_settings(cs.crystal_settings(_settings(), w.current_data()))
    dlg = BatchStudyDialog(w)
    dlg.group_edit.setText("2011")
    dlg.add_from_main()
    for k, x0 in enumerate((3.0, -3.0)):
        dlg.batch.add(f"opx{k}", cs.CrystalData(_synthetic(x0=x0, seed=k + 5), SPEC, f"opx{k}.csv"),
                      group="2018")
    dlg.refresh()
    assert dlg.table.rowCount() == 4
    assert dlg.table.item(SHARED_ROW, KEYS.index("_label")).text() == "All crystals"
    # an own value in one row, then a shared change that keeps it
    dlg.set_value(2, "T_C", 980.0)
    assert dlg.table.item(2, KEYS.index("T_C")).font().bold()
    monkeypatch.setattr(QMessageBox, "question", lambda *a, **k: QMessageBox.No)
    dlg.set_value(SHARED_ROW, "T_C", 940.0)
    assert [dlg.batch.settings(i)["T_C"] for i in range(3)] == [940.0, 980.0, 940.0]
    # an angle typed in the shared row reaches every crystal
    dlg.set_value(SHARED_ROW, "angle_b", 45.0)
    assert all(dlg.batch.settings(i)["angles_deg"][1] == 45.0 for i in range(3))
    dlg.set_value(SHARED_ROW, "angle_b", 90.0)
    dlg.set_value(SHARED_ROW, "T_C", 950.0)
    dlg.batch.clear_own(1, "T_C")
    dlg.run(False)
    _wait(app, lambda: dlg.thread is not None and dlg.thread.isFinished())
    app.processEvents()
    assert all(c.result is not None for c in dlg.batch.crystals), [c.status for c in dlg.batch.crystals]
    for c in dlg.batch.crystals:
        assert c.result["t_seconds"] / SEC_PER_YEAR == pytest.approx(1.0, rel=0.15)
    dlg.send_to_population()
    recs = w.population_study.records
    assert [r.group for r in recs] == ["2011", "2018", "2018"]
    w.close()
