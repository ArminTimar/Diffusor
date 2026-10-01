"""The bundled dataset catalogue, the superseded-coefficient flags, and the interface."""
import os
from pathlib import Path

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


# --- regressions from user testing ---------------------------------------------------
def _load_example(app, key):
    from PySide6.QtCore import Qt
    w = _window(app)
    for i in range(w.lst_examples.count()):
        if w.lst_examples.item(i).data(Qt.UserRole) == key:
            w.lst_examples.setCurrentRow(i)
    w.load_example()
    app.processEvents()
    return w


def test_magnetite_example_fits_its_plateaus_and_recovers_8_days(app):
    """Two regressions. The initial-profile guess ran after the example settings
    and replaced x_Ti = 0.1 by the mean TiO2 (5.3 wt%, clamped to 1), which made
    D 250 times too large and the fit about 0.02 d. And the traverse ends before
    the profile flattens: with the plateaus held at the outer points the fit gave
    about 5.3 d for a true 8 d."""
    from diffusor.constants import SEC_PER_DAY
    from diffusor.fitting import fit_time
    w = _load_example(app, "magnetite_shinmoedake")
    assert w.sp_xcomp.value() == pytest.approx(0.1)
    assert w.chk_free_plateaus.isChecked()
    assert set(w._free_parameters()) == {"t", "x0", "C_left", "C_right"}
    r = fit_time(w._model(w._checked_keys()[0]), w.profile.x, w.profile.C, w.profile.sigma,
                 w._free_parameters())
    assert 0.75 < r.t_seconds / (8 * SEC_PER_DAY) < 1.33
    w.close()
    w = _load_example(app, "opx_shinmoedake")
    assert not w.chk_free_plateaus.isChecked(), "the next example must not inherit the setting"
    w.close()


def test_the_coefficient_is_chosen_before_the_conditions(app):
    """The law decides what the later steps need, so it comes right after the
    mineral; the resolution sits with the data, the cooling path and the host
    composition with the conditions."""
    from diffusor.gui import main_window as mw
    assert mw.STEPS == ["Data", "Mineral", "Coefficient", "Conditions", "Model",
                        "Uncertainty", "Results"]
    w = _window(app)
    page = w.pages.widget
    assert page(mw.DATA).isAncestorOf(w.cmb_resolution)
    assert page(mw.DATA).isAncestorOf(w.sp_xscale_sig)
    assert page(mw.COEFFICIENT).isAncestorOf(w.lst_coef)
    for wdg in (w.sp_T, w.chk_cooling, w.sp_Tend, w.sp_dbuf, w.sp_xcomp, w.chk_comp_dep,
                w.box_unused):
        assert page(mw.CONDITIONS).isAncestorOf(wdg), wdg
    assert page(mw.MODEL).isAncestorOf(w.cmb_solver)
    assert page(mw.MINERAL).isAncestorOf(w.cmb_axis)
    w.close()


def test_conditions_the_law_ignores_are_switched_off(app):
    from PySide6.QtCore import Qt
    from diffusor.gui.main_window import CONDITIONS

    def tick(w, key):
        for i in range(w.lst_coef.count()):
            it = w.lst_coef.item(i)
            it.setCheckState(Qt.Checked if it.data(Qt.UserRole) == key else Qt.Unchecked)
        w._go(CONDITIONS)

    w = _load_example(app, "cpx_stromboli")         # Mueller et al. (2013): no P, no fO2
    tick(w, "cpx_FeMg_muller2013")
    assert not w.sp_P.isEnabled() and not w.sp_dbuf.isEnabled() and not w.cmb_buffer.isEnabled()
    assert w.box_unused.isVisibleTo(w)
    w.close()

    w = _load_example(app, "magnetite_shinmoedake")  # fO2 term, no pressure term
    tick(w, "mt_Ti_vanorman_crispin2010")
    assert w.sp_dbuf.isEnabled()
    assert w.cmb_fo2_mode.currentIndex() == 1 and not w.sp_P.isEnabled(), \
        "absolute fO2: pressure enters nowhere"
    w.cmb_fo2_mode.setCurrentIndex(0)
    assert w.sp_P.isEnabled() and w.cmb_buffer.isEnabled(), "a buffer moves with pressure"
    w.close()


def _load_own_file(w, filename, spec):
    """Load a bundled file the way a user's own file comes in: no dataset settings."""
    path = ds.EXAMPLES_DIR / filename
    w._use_table(read_table(path), spec, path, None, (None, False))


def _choose(w, mineral, species):
    keys = [w.cmb_mineral.itemData(i) for i in range(w.cmb_mineral.count())]
    w.cmb_mineral.setCurrentIndex(keys.index(mineral))
    w.cmb_species.setCurrentText(species)


TIO2 = ProfileSpec("Distance_um", "TiO2_wt", None, "TiO2_err", mode="A")
FEMG = ProfileSpec("Distance_um", "FeO_wt", "MgO_wt", "FeO_err", "MgO_err",
                   mode="A/(A+B)", oxide_a="FeO", oxide_b="MgO")


def test_a_concentration_profile_never_becomes_the_host_composition(app):
    """A TiO2 wt% mean (about 5) used to land in x_Ti, clamped to 1, with D then
    about 250 times too large. The profile is usually loaded before the mineral
    is chosen, so the value has to follow the mineral and species as well."""
    w = _window(app)
    before = w.sp_xcomp.value()
    _load_own_file(w, "magnetite_ti.csv", TIO2)
    _choose(w, "magnetite", "Ti")
    assert w.sp_xcomp.value() == pytest.approx(before)
    w.sp_xcomp.setValue(0.1)                     # typed in
    w._guess_initial(from_button=True)
    assert w.sp_xcomp.value() == pytest.approx(0.1)
    m = w._model("mt_Ti_vanorman_crispin2010")
    assert m.D_bulk() == pytest.approx(4.35e-16, rel=0.01)   # Tomiya et al. (2013)
    w.close()


def test_the_profile_mean_is_taken_back_when_the_species_changes(app):
    w = _window(app)
    _choose(w, "opx", "Fe-Mg")
    before = w.sp_xcomp.value()
    _load_own_file(w, "opx_femg_step.csv", FEMG)
    xfe = float(np.mean(w.profile.C))
    assert 0.2 < xfe < 0.35 and w.sp_xcomp.value() == pytest.approx(xfe, abs=1e-4)
    _choose(w, "magnetite", "Ti")                # the X_Fe mean means nothing for x_Ti
    assert w.sp_xcomp.value() == pytest.approx(before)
    _choose(w, "opx", "Fe-Mg")
    assert w.sp_xcomp.value() == pytest.approx(xfe, abs=1e-4)
    w.close()


def test_a_typed_composition_survives_a_mineral_change_but_not_guess(app):
    w = _window(app)
    _choose(w, "opx", "Fe-Mg")
    _load_own_file(w, "opx_femg_step.csv", FEMG)
    w.sp_xcomp.setValue(0.4)
    _choose(w, "cpx", "Fe-Mg")
    assert w.sp_xcomp.value() == pytest.approx(0.4)
    w._guess_initial(from_button=True)           # asked for: the profile mean
    assert w.sp_xcomp.value() == pytest.approx(float(np.mean(w.profile.C)), abs=1e-4)
    w.close()


def test_olivine_composition_goes_through_the_chosen_coordinate(app):
    w = _window(app)
    _choose(w, "olivine", "Fe-Mg")
    before = w.sp_xcomp.value()
    _load_own_file(w, "olivine_fo.csv", ProfileSpec("Distance_um", "Fo_mol", None, "Fo_err", mode="A"))
    assert w.sp_xcomp.value() == pytest.approx(before), "Fo ~85 is no X_Fe"
    w.cmb_ol_coordinate.setCurrentIndex(2)       # Forsterite mol%
    xfe = 1.0 - float(np.mean(w.profile.C)) / 100.0
    assert 0.1 < xfe < 0.25 and w.sp_xcomp.value() == pytest.approx(xfe, abs=1e-4)
    w.close()


def test_compare_runs_to_completion_and_ignores_a_second_click(app, monkeypatch):
    """Compare used to crash: its progress lambda touched widgets from the worker thread."""
    import inspect
    import time
    from PySide6.QtCore import Qt
    from diffusor.gui import main_window as mw
    for fn in (mw.MainWindow.run_compare, mw.MainWindow.run_mc, mw.MainWindow.run_fit):
        src = inspect.getsource(fn)
        assert "lambda" not in src, f"{fn.__name__} connects a worker signal to a lambda"
    monkeypatch.setattr(mw.richtext, "show", lambda *a, **k: None)
    w = _load_example(app, "opx_kizimen")
    for i in range(w.lst_coef.count()):
        it = w.lst_coef.item(i)
        if it.data(Qt.UserRole) == "opx_FeMg_dias2025":
            it.setCheckState(Qt.Checked)
    w._go(6)
    w.run_compare()
    w.run_fit()              # a second job while one runs must be refused, not crash
    assert len(w._jobs) == 1
    t0 = time.time()
    while w.compare_results is None and time.time() - t0 < 300:
        app.processEvents()
        time.sleep(0.02)
    assert w.compare_results is not None and len(w.compare_results) == 2
    w.close()


def test_mouse_wheel_does_not_change_numbers(app):
    from PySide6.QtCore import QPoint, QPointF, Qt
    from PySide6.QtGui import QWheelEvent
    from PySide6.QtWidgets import QApplication
    from diffusor.gui.main_window import CONDITIONS
    w = _load_example(app, "opx_kizimen")
    w._go(CONDITIONS)
    app.processEvents()
    before = w.sp_T.value()
    ev = QWheelEvent(QPointF(5, 5), QPointF(5, 5), QPoint(0, 0), QPoint(0, 120),
                     Qt.NoButton, Qt.NoModifier, Qt.NoScrollPhase, False)
    QApplication.sendEvent(w.sp_T, ev)
    QApplication.sendEvent(w.cmb_fo2_mode, ev)
    assert w.sp_T.value() == before
    assert w.cmb_fo2_mode.currentIndex() == 0
    w.close()


def test_example_says_what_it_filled_in(app):
    from diffusor.gui.main_window import COEFFICIENT, CONDITIONS, DATA, MINERAL, MODEL
    w = _load_example(app, "opx_kizimen")
    assert w.lbl_data.isVisibleTo(w.pages.widget(DATA))
    assert "filled in" in w.lbl_data.text()
    for step in (MINERAL, COEFFICIENT, CONDITIONS, MODEL):
        assert not w._prefill[step].isHidden(), step
    assert "850" in w._prefill[CONDITIONS].text()
    assert "Ganguly" in w._prefill[COEFFICIENT].text()
    w.close()


def test_anorthite_initial_profile_is_offered_only_for_plagioclase(app):
    w = _load_example(app, "opx_kizimen")
    assert [w.cmb_ic.itemData(i) for i in range(w.cmb_ic.count())] == ["step"]
    w.close()
    w = _load_example(app, "plag_santorini")
    assert w.cmb_ic.currentData() == "equilibrium_plag"
    w.close()


def test_resolution_presets_set_sigma_from_the_spot_size(app):
    from diffusor.gui.main_window import RESOLUTION_PRESETS
    w = _window(app)
    labels = [p[0] for p in RESOLUTION_PRESETS]
    w.cmb_resolution.setCurrentIndex(labels.index("Microprobe, focused beam"))
    assert w.sp_beam.value() == pytest.approx(0.6)
    w.cmb_resolution.setCurrentIndex(labels.index("LA-ICP-MS spot"))
    w.sp_width.setValue(10.0)
    assert w.sp_beam.value() == pytest.approx(2.5)
    w.cmb_resolution.setCurrentIndex(labels.index("LA-ICP-MS line scan"))
    w.sp_width.setValue(7.5)
    assert w.sp_beam.value() == pytest.approx(7.5 / np.sqrt(12), abs=0.01)
    assert "Grocolas" in w.lbl_resolution.text()
    w.close()


def test_column_dialog_guesses_and_fits_a_laptop_screen(app):
    from diffusor.gui.main_window import ColumnDialog
    from diffusor.dataio.profiles import suggest_spec
    d = ds.get("opx_kizimen")
    df = read_table(d.path)
    spec = suggest_spec(df)
    assert (spec.column_a, spec.column_b, spec.oxide_a, spec.oxide_b) == ("FeO_wt", "MgO_wt", "FeO", "MgO")
    assert spec.distance_column == "Distance_from_rim_um" and spec.distance_unit == "um"
    dlg = ColumnDialog(df, spec)
    dlg.show()
    app.processEvents()
    assert dlg.sizeHint().height() < 640, dlg.sizeHint().height()
    got = dlg.spec()
    assert (got.column_a, got.column_b, got.oxide_a, got.mode) == ("FeO_wt", "MgO_wt", "FeO", "A/(A+B)")
    dlg.close()


def test_interface_text_has_no_semicolons(app):
    from PySide6.QtWidgets import QAbstractButton, QLabel
    w = _load_example(app, "plag_santorini")
    for page in range(w.pages.count()):
        w._go(page)
        app.processEvents()
        for lab in w.pages.widget(page).findChildren(QLabel):
            assert ";" not in lab.text().replace("&middot;", ""), lab.text()
        for b in w.pages.widget(page).findChildren(QAbstractButton):
            assert ";" not in b.text()
    for d in ds.DATASETS:
        assert ";" not in d.provenance + d.expected + d.notes, d.key
    w.close()


# --- Help > All references ---------------------------------------------------------
def test_all_references_lists_every_reference_once_or_more(app):
    from diffusor.gui.reference_list import reference_groups, references_html
    placed = {k for _, keys in reference_groups() for k in keys}
    assert placed == set(REFERENCES)
    text = references_html()
    for key, ref in REFERENCES.items():
        assert f">{key}<" in text, key
        if ref.doi:
            assert f"https://doi.org/{ref.doi}" in text, key
    assert f"{len(REFERENCES)} references" in text


def test_all_references_shows_what_uses_each_paper(app):
    from diffusor.gui.reference_list import reference_usage
    usage = reference_usage()
    assert "ol_FeMg_dohmen_chakraborty2007_tamed" in usage["dohmen_chakraborty2007"]
    assert "opx_FeMg_dohmen2016 (secondary)" in usage["sato2022"]
    assert "example plag_santorini" in usage["druitt2012"]


def test_all_references_search_narrows_the_list(app):
    from diffusor.gui.reference_list import references_html
    text = references_html("Sakurajima zzz-no-such-paper")
    assert "No references match" in text
    text = references_html("Polo-Sanchez Kameni")
    assert ">polo_sanchez2023<" in text and ">crank1975<" not in text
    assert f"1 of {len(REFERENCES)} references match" in text


def test_all_references_saves_the_bibtex_file(app, tmp_path):
    from pathlib import Path
    from diffusor.gui.reference_list import ReferenceListDialog
    from diffusor.references import bibtex_document
    dlg = ReferenceListDialog()
    out = dlg.save_bibtex(str(tmp_path / "refs.bib"))
    saved = Path(out).read_text(encoding="utf-8")
    assert saved == bibtex_document()
    dlg.search.setText("Dohmen")
    assert "Dohmen" in dlg.view.toPlainText()
    dlg.close()


def test_help_menu_opens_all_references(app):
    w = _window(app)
    help_menu = [a.menu() for a in w.menuBar().actions() if a.text() == "&Help"][0]
    assert "All references" in [a.text() for a in help_menu.actions()]
    w.close()


# --- recent profiles on the Data step --------------------------------------------------
def _fake_column_dialog(monkeypatch, on_exec=None):
    """Accept the column dialog without showing it; count the times it would open."""
    from PySide6.QtWidgets import QDialog
    from diffusor.gui import main_window as mw
    seen = []

    def fake_exec(dlg):
        seen.append(dlg)
        if on_exec:
            on_exec(dlg, len(seen))
        return QDialog.Accepted
    monkeypatch.setattr(mw.ColumnDialog, "exec", fake_exec)
    return seen


def test_recent_profiles_remember_the_file_and_its_columns(app, tmp_path, monkeypatch):
    import shutil
    src = tmp_path / "mt.csv"
    shutil.copy(ds.EXAMPLES_DIR / "magnetite_ti.csv", src)

    def first_time(dlg, n):
        if n == 1:
            dlg.xmax.setText("60")          # a choice the guess would not make
    seen = _fake_column_dialog(monkeypatch, first_time)
    w = _window(app)
    assert w.lst_recent.count() == 0 and w.data_stack.currentIndex() == 0
    w._load(src, dataset=None)
    assert len(seen) == 1 and w.profile.spec.x_max == 60
    assert w.lst_recent.count() == 1 and w.data_stack.currentIndex() == 0, \
        "after loading, the profile is in front"
    w.btn_data_view.click()
    assert w.data_stack.currentIndex() == 1 and w.lbl_data_title.text() == "Recent profiles"
    w._load_recent_item(w.lst_recent.item(0))
    assert len(seen) == 1, "an unchanged recent file loads without the column dialog"
    assert w.profile.spec.x_max == 60 and w.data_stack.currentIndex() == 0
    w.close()

    # a new window opens on the recent list, and a changed file asks again,
    # starting from last time's choices
    df = read_table(src)
    df["extra"] = 1.0
    df.to_csv(src, index=False)
    w = _window(app)
    assert w.data_stack.currentIndex() == 1 and w.lst_recent.count() == 1
    w._load_recent_item(w.lst_recent.item(0))
    assert len(seen) == 2 and seen[-1].xmax.text() == "60"
    w.close()


def test_a_missing_recent_file_is_shown_but_not_loaded(app, tmp_path, monkeypatch):
    import shutil
    from PySide6.QtCore import Qt
    src = tmp_path / "gone.csv"
    shutil.copy(ds.EXAMPLES_DIR / "magnetite_ti.csv", src)
    _fake_column_dialog(monkeypatch)
    w = _window(app)
    w._load(src, dataset=None)
    src.unlink()
    w._refresh_recent()
    item = w.lst_recent.item(0)
    assert "not found" in item.text() and not item.flags() & Qt.ItemIsEnabled
    before = w.profile
    w._load_recent_item(item)
    assert w.profile is before
    w.btn_recent_clear.click()
    assert w.lst_recent.count() == 0
    w.close()


def test_the_file_dialog_opens_in_the_last_folder(app, tmp_path, monkeypatch):
    import shutil
    from diffusor.gui import main_window as mw
    src = tmp_path / "mt.csv"
    shutil.copy(ds.EXAMPLES_DIR / "magnetite_ti.csv", src)
    _fake_column_dialog(monkeypatch)
    asked = []
    monkeypatch.setattr(mw.QFileDialog, "getOpenFileName",
                        staticmethod(lambda parent, title, start, filt: asked.append(start)
                                     or (str(src), "")))
    w = _window(app)
    w.load_file()
    w.load_file()
    assert Path(asked[1]) == tmp_path
    w.close()


def test_examples_do_not_enter_the_recent_list(app):
    w = _load_example(app, "magnetite_shinmoedake")
    assert w.lst_recent.count() == 0
    w.close()
