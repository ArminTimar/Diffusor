"""The bundled dataset catalogue, the superseded-coefficient flags, and the interface."""
import os

import numpy as np
import pytest

from diffusor import datasets as ds
from diffusor.coefficients import Conditions, get, list_coefficients
from diffusor.dataio import ProfileSpec, build_profile, read_table
from diffusor.references import REFERENCES


# --- catalogue -----------------------------------------------------------------
def test_every_dataset_declares_its_provenance():
    for d in ds.DATASETS:
        assert d.kind in ("measured", "synthetic"), d.key
        assert d.provenance.strip(), d.key
        if d.kind == "measured":
            assert d.citation in REFERENCES, f"{d.key} cites unknown key {d.citation}"
            assert "Real measured data" in d.provenance
        else:
            assert "SYNTHETIC" in d.provenance, (
                f"{d.key} must say in its provenance that it is not a measurement")


def test_provenance_banner_never_calls_synthetic_data_measured():
    for d in ds.DATASETS:
        banner = d.provenance_banner()
        if d.kind == "synthetic":
            assert "not a measurement" in banner
        else:
            assert banner.startswith("MEASURED DATA")


def test_there_is_a_dataset_for_every_mineral():
    covered = {d.mineral for d in ds.DATASETS}
    assert {"olivine", "opx", "cpx", "plagioclase", "magnetite"} <= covered


def test_every_dataset_file_exists_and_loads_with_its_own_spec():
    for d in ds.DATASETS:
        assert d.exists, f"{d.filename} is missing; run scripts/make_examples.py"
        prof = build_profile(read_table(d.path), ProfileSpec(**d.spec), str(d.path))
        assert len(prof) >= 10, d.key
        assert np.all(np.isfinite(prof.C))
        assert np.all(np.diff(prof.x) >= 0), "loader must return a sorted traverse"


def test_every_dataset_names_a_coefficient_that_exists():
    for d in ds.DATASETS:
        key = d.settings.get("coefficient")
        assert key, d.key
        c = get(key)
        assert c.mineral == d.mineral, f"{d.key} points at a {c.mineral} coefficient"


def test_santorini_dataset_is_the_real_published_traverse():
    d = ds.get("plag_santorini")
    assert d.kind == "measured" and d.citation == "druitt2012"
    df = read_table(d.path)
    assert len(df) == 14
    an = df["An_mol_percent"].to_numpy()
    # normally zoned: anorthite-rich core, sodic rim
    assert an.min() == pytest.approx(36.6)
    assert an.max() == pytest.approx(80.3)
    assert df["Distance_from_rim_um"].iloc[0] == pytest.approx(21.0)
    assert "does NOT reproduce" in d.expected, (
        "the catalogue must be honest that this is not a validation")


# --- superseded flags -----------------------------------------------------------
def test_superseded_coefficients_are_flagged_and_not_recommended():
    for c in list_coefficients():
        if c.superseded_by or c.superseded_note:
            assert not c.recommended, f"{c.key} is superseded but still recommended"
            assert c.superseded_note.strip(), c.key
            if c.superseded_by:
                assert c.superseded_by in REFERENCES


def test_plagioclase_sr_points_at_the_2025_recalibration():
    c = get("plag_Sr_giletti_casserly1994")
    assert c.superseded_by == "grocolas2025"
    assert "orders of magnitude" in c.superseded_note.lower()
    w = c.check_conditions(Conditions(T_K=1173.15, X={"XAn": 0.5}))
    assert any("SUPERSEDED" in s for s in w)


def test_magnetite_feti_1978_is_demoted_in_favour_of_the_table12_entries():
    old = get("mt_FeTi_freer_hauptman1978")
    assert not old.recommended
    assert "Table 12" in old.superseded_note
    assert get("mt_Ti_vanorman_crispin2010").recommended


def test_recommended_entries_are_never_superseded_or_unverified():
    for c in list_coefficients():
        if c.recommended:
            assert not (c.superseded_by or c.superseded_note), c.key


# --- interface ------------------------------------------------------------------
@pytest.fixture(scope="module")
def app():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    import matplotlib
    matplotlib.use("QtAgg")
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _window(app):
    from diffusor.gui.main_window import MainWindow
    w = MainWindow()
    w.resize(1440, 900)
    w.show()
    app.processEvents()
    return w


def test_window_builds_with_one_page_per_step(app):
    from diffusor.gui.main_window import STEPS
    w = _window(app)
    assert w.pages.count() == len(STEPS)
    assert len(w.step_labels) == len(STEPS)
    w.close()


def test_the_flow_starts_on_the_data_step_and_will_not_advance_without_data(app):
    w = _window(app)
    assert w.step == 0
    w._go(3)
    assert w.step == 0, "must not skip ahead before a profile is loaded"
    w.close()


def test_loading_an_example_fills_in_the_later_steps(app):
    from PySide6.QtCore import Qt
    w = _window(app)
    for i in range(w.lst_examples.count()):
        if w.lst_examples.item(i).data(Qt.UserRole) == "opx_shinmoedake":
            w.lst_examples.setCurrentRow(i)
            break
    w.load_example()
    app.processEvents()
    assert w.profile is not None and len(w.profile) > 10
    assert w.dataset.key == "opx_shinmoedake"
    assert w.cmb_mineral.currentData() == "opx"
    assert w.cmb_species.currentText() == "Fe-Mg"
    assert w.sp_T.value() == pytest.approx(950.0)
    assert w._checked_keys() == ["opx_FeMg_dias2025"]
    w._go(6)
    assert w.step == 6
    w.close()


def test_santorini_example_switches_on_the_activity_term(app):
    from PySide6.QtCore import Qt
    w = _window(app)
    for i in range(w.lst_examples.count()):
        if w.lst_examples.item(i).data(Qt.UserRole) == "plag_santorini":
            w.lst_examples.setCurrentRow(i)
            break
    w.load_example()
    app.processEvents()
    assert w.an_values is not None
    assert 0.3 < float(np.min(w.an_values)) < 0.9
    assert w.cmb_ic.currentIndex() == 1, "should use the equilibrium initial condition"
    m = w._model(w._checked_keys()[0])
    assert m.activity_theta != 0.0
    assert m.an_profile is not None
    assert m.initial.kind == "equilibrium_plag"
    w.close()


def test_summary_sidebar_keeps_a_fixed_width_and_does_not_clip(app):
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import QLabel, QSplitter
    w = _window(app)
    for i in range(w.lst_examples.count()):
        if w.lst_examples.item(i).data(Qt.UserRole) == "opx_shinmoedake":
            w.lst_examples.setCurrentRow(i)
            break
    w.load_example()
    w._go(6)
    app.processEvents()
    app.processEvents()
    splitter = w.pages.widget(6).findChild(QSplitter)
    side, plot = splitter.sizes()
    assert 240 <= side <= 360, f"sidebar should stay narrow, got {side}"
    assert plot > side, "the plot must get the bulk of the width"
    clipped = [l for l in w.summary_inner.findChildren(QLabel)
               if l.isVisible() and l.wordWrap() and l.height() + 1 < l.heightForWidth(l.width())]
    assert not clipped, f"{len(clipped)} summary labels are cut off"
    w.close()


def test_no_label_on_any_step_is_clipped(app):
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import QLabel
    from diffusor.gui.main_window import STEPS
    w = _window(app)
    for i in range(w.lst_examples.count()):
        if w.lst_examples.item(i).data(Qt.UserRole) == "opx_shinmoedake":
            w.lst_examples.setCurrentRow(i)
            break
    w.load_example()
    for page in range(len(STEPS)):
        w._go(page)
        app.processEvents()
        app.processEvents()
        clipped = [l for l in w.pages.widget(page).findChildren(QLabel)
                   if l.isVisible() and l.wordWrap()
                   and l.height() + 1 < l.heightForWidth(l.width())]
        assert not clipped, f"step '{STEPS[page]}' clips {len(clipped)} labels"
    w.close()


def test_superseded_coefficients_are_marked_in_the_chooser(app):
    w = _window(app)
    keys = [w.cmb_mineral.itemData(i) for i in range(w.cmb_mineral.count())]
    w.cmb_mineral.setCurrentIndex(keys.index("plagioclase"))
    w.cmb_species.setCurrentText("Sr")
    app.processEvents()
    texts = [w.lst_coef.item(i).text() for i in range(w.lst_coef.count())]
    assert any("superseded" in t for t in texts), texts
    w.close()


def test_a_fit_runs_end_to_end_through_the_interface(app):
    from PySide6.QtCore import Qt
    from diffusor.constants import SEC_PER_YEAR
    from diffusor.fitting import fit_time
    w = _window(app)
    for i in range(w.lst_examples.count()):
        if w.lst_examples.item(i).data(Qt.UserRole) == "opx_shinmoedake":
            w.lst_examples.setCurrentRow(i)
            break
    w.load_example()
    w._go(6)
    model = w._model(w._checked_keys()[0])
    r = fit_time(model, w.profile.x, w.profile.C, w.profile.sigma, w._free_parameters())
    w._fit_done(r)
    app.processEvents()
    assert 0.5 < r.t_seconds / (1.5 * SEC_PER_YEAR) < 2.0
    assert w.fit_result is not None
    w.close()


def test_kizimen_dataset_is_the_real_ostorero_traverse():
    d = ds.get("opx_kizimen")
    assert d.kind == "measured" and d.citation == "ostorero2022"
    df = read_table(d.path)
    assert len(df) == 109
    assert df["Distance_from_rim_um"].iloc[0] == pytest.approx(2.28)
    prof = build_profile(df, ProfileSpec(**d.spec))
    assert prof.x.min() >= 4.0 and prof.x.max() <= 60.0
    # reverse zone: Fe-poor rim band against a Fe-richer core
    assert prof.C[prof.x < 11].max() < 0.30 < prof.C[prof.x > 20].min()


def test_kizimen_fit_lands_inside_ostorero_uncertainty(app):
    """Ostorero et al. (2022) Supplementary Data 4: K9_L10C4 = 2.32 yr (+7.16 / -1.75)."""
    from PySide6.QtCore import Qt
    from diffusor.constants import SEC_PER_YEAR
    from diffusor.fitting import fit_time
    w = _window(app)
    for i in range(w.lst_examples.count()):
        if w.lst_examples.item(i).data(Qt.UserRole) == "opx_kizimen":
            w.lst_examples.setCurrentRow(i)
            break
    w.load_example()
    app.processEvents()
    assert w._checked_keys() == ["opx_FeMg_ganguly_tazzoli1994_nofo2"]
    assert w.sp_T.value() == pytest.approx(850.0)
    model = w._model(w._checked_keys()[0])
    assert model.comp_key == "XFe"
    r = fit_time(model, w.profile.x, w.profile.C, w.profile.sigma, w._free_parameters())
    years = r.t_seconds / SEC_PER_YEAR
    assert 2.32 - 1.75 < years < 2.32 + 7.16, years
    w.close()
