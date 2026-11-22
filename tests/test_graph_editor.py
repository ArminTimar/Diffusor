"""The graph editor of the results plot, and the interface fixes of November 2026."""
import os

import numpy as np
import pytest


@pytest.fixture(scope="module")
def app():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    import matplotlib
    matplotlib.use("QtAgg")
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _fit():
    """A real fit of the sanidine Ba example, small and quick."""
    from diffusor import datasets as ds
    from diffusor.coefficients import Conditions, get
    from diffusor.dataio import ProfileSpec, build_profile, read_table
    from diffusor.fitting import DiffusionModel, fit_time
    from diffusor.solvers import InitialCondition, dirichlet
    d = ds.get("sanidine_ba")
    p = build_profile(read_table(d.path), ProfileSpec(**d.spec))
    ic = InitialCondition("step", {"x0": 10.0, "C_left": 4800.0, "C_right": 1600.0})
    m = DiffusionModel(coefficient=get("kfs_Ba_cherniak2002"), conditions=Conditions(T_K=1063.15),
                       initial=ic, bc_left=dirichlet(4800.0), bc_right=dirichlet(1600.0))
    return fit_time(m, p.x, p.C, p.sigma, ["t", "x0"])


def _curve(plot, ident):
    from diffusor.gui.graph_editor import curves
    for ax in plot.figure.axes:
        for c in curves(ax):
            if c.ident == ident:
                return c
    raise AssertionError(f"{ident} not drawn")


def test_the_toolbars_drop_matplotlibs_subplot_and_options_dialogs(app):
    from diffusor.gui.plot_widget import GraphToolbar, SaveToolbar
    plain = [cb for _t, _tip, _i, cb in SaveToolbar.toolitems if cb]
    graph = [cb for _t, _tip, _i, cb in GraphToolbar.toolitems if cb]
    assert "configure_subplots" not in plain + graph
    assert "edit_parameters" not in plain
    assert "edit_parameters" in graph
    names = [t for t, _tip, _i, cb in GraphToolbar.toolitems if cb == "edit_parameters"]
    assert names == ["Edit graph"]


def test_cut_points_tooltip_has_several_lines(app):
    from diffusor.gui.plot_widget import ProfilePlot
    plot = ProfilePlot()
    tip = plot.act_cut.toolTip()
    assert tip.count("<br>") >= 3, tip
    plot.set_points(np.arange(5.0), np.arange(5.0), np.array([1, 0, 0, 0, 1], bool))
    assert "2 points cut" in plot.act_cut.toolTip()
    plot.deleteLater()


def test_the_editor_lists_every_axis_once_and_every_curve(app):
    from diffusor.gui.graph_editor import GraphEditor, axis_rows, curves
    from diffusor.gui.plot_widget import ProfilePlot
    plot = ProfilePlot()
    plot.show_fit(_fit())
    rows = [name for _ax, _w, name in axis_rows(plot.figure)]
    # the profile and residual panels share the distance axis, listed once
    assert rows == ["Distance", "Composition", "Residual"]
    keys = {c.key for ax in plot.figure.axes for c in curves(ax)}
    assert {"Measured points", "Best fit", "Initial profile", "Residuals", "Zero line"} <= keys
    editor = GraphEditor(plot)
    editor.rebuild()
    editor.close()
    editor.deleteLater()
    plot.deleteLater()


def test_changes_apply_at_once_and_survive_a_new_fit(app):
    from matplotlib.colors import to_hex
    from diffusor.gui.plot_widget import ProfilePlot
    plot = ProfilePlot()
    r = _fit()
    plot.show_fit(r)
    plot.overrides["curves"]["profile/Best fit"] = {"color": "#1f77b4", "linewidth": 3.5,
                                                   "label": "my fit"}
    plot.overrides["curves"]["profile/Initial profile"] = {"visible": False}
    plot.overrides["axes"]["residual"] = {"xlabel": "distance from the rim (µm)",
                                          "xlim": (2.0, 20.0)}
    plot.overrides["figure"]["title"] = "Crystal 4"
    plot.restyle()

    def check():
        fit = _curve(plot, "profile/Best fit").first
        assert to_hex(fit.get_color()) == "#1f77b4" and fit.get_linewidth() == 3.5
        assert not _curve(plot, "profile/Initial profile").visible()
        assert plot.ax_res.get_xlabel() == "distance from the rim (µm)"
        assert plot.ax.get_xlim() == pytest.approx((2.0, 20.0))     # shared with the residuals
        assert plot.ax.get_title() == "Crystal 4"
        legend = [t.get_text() for t in plot.ax.get_legend().get_texts()]
        assert "my fit" in legend and "initial profile" not in legend, legend
    check()
    plot.show_fit(r)            # the next fit draws everything anew
    check()
    plot.forget_limits()        # what loading another profile does
    plot.show_fit(r)
    assert plot.ax.get_xlim() != pytest.approx((2.0, 20.0))
    plot.reset_overrides()
    assert not _curve(plot, "profile/Initial profile").first is None
    assert _curve(plot, "profile/Initial profile").visible()
    assert plot.ax.get_title() == ""
    plot.deleteLater()


def test_the_legend_lists_the_data_first(app):
    from diffusor.gui.plot_widget import ProfilePlot
    plot = ProfilePlot()
    plot.show_fit(_fit())
    texts = [t.get_text() for t in plot.ax.get_legend().get_texts()]
    assert texts[0] == "measured" and texts[1].startswith("fit:"), texts
    # the upper panel shares the distance axis, whose numbers are on the residual panel
    assert not any(t.get_visible() and t.get_text() for t in plot.ax.get_xticklabels())
    plot.deleteLater()


def _window(app):
    from diffusor.gui.main_window import MainWindow
    w = MainWindow()
    w.resize(1250, 620)
    w.show()
    app.processEvents()
    return w


def test_the_mineral_note_follows_the_species(app):
    w = _window(app)
    w.cmb_mineral.setCurrentIndex(w.cmb_mineral.findData("opx"))
    w.cmb_species.setCurrentText("Fe-Mg")
    fe_mg = w.lbl_mineral_note.text()
    w.cmb_species.setCurrentText("Lu")
    lu = w.lbl_mineral_note.text()
    assert fe_mg != lu
    assert "The profile is X_Fe" in fe_mg
    assert "measured Lu content" in lu
    w.close()


def test_the_coefficient_step_makes_no_checked_against_claim(app):
    w = _window(app)
    w.cmb_mineral.setCurrentIndex(w.cmb_mineral.findData("opx"))
    w.cmb_species.setCurrentText("Fe-Mg")
    w.lst_coef.setCurrentRow(0)
    app.processEvents()
    assert "checked against" not in w.lbl_coef.text().lower()
    assert " -- " not in w.lst_coef.item(0).text()
    w.close()


def test_the_footer_counts_all_seven_steps(app):
    w = _window(app)
    assert w.lbl_footer.text() == "Step 1 of 7"
    w.close()


def test_the_multicomponent_fit_buttons_stay_in_view(app):
    """They were at the end of a scrolling column, under the fold on a laptop screen."""
    from PySide6.QtWidgets import QScrollArea
    from diffusor.gui.multicomponent_study import MulticomponentStudyDialog
    w = _window(app)
    dlg = MulticomponentStudyDialog(w)
    for button in (dlg.g_run, dlg.i_run):
        assert not any(a.widget().isAncestorOf(button) for a in dlg.findChildren(QScrollArea)
                       if a.widget() is not None)
    dlg.close()
    w.close()


def test_the_joint_study_enables_fit_with_two_profiles_only(app):
    from diffusor.gui.joint_study import JointStudyDialog
    w = _window(app)
    dlg = JointStudyDialog(w)
    assert not dlg.fit.isEnabled() and not dlg.remove.isEnabled()
    assert "two or more" in dlg.status.text()
    dlg.close()
    w.close()


def test_the_validation_library_names_studies_and_loads_on_double_click(app):
    from diffusor.gui.validation_library import ValidationLibraryDialog, study_name
    from diffusor.validation import records
    assert study_name("iovine2017") == "Iovine et al. (2017)"
    assert study_name("sundermeyer2020_eifel") == "Sundermeyer 2020, Eifel"
    dlg = ValidationLibraryDialog(records())
    item = dlg.items.item(0)
    assert not item.text().startswith(item.data(0x0100)["study"] + " ")    # no raw key
    dlg.items.itemDoubleClicked.emit(item)
    from PySide6.QtWidgets import QDialog
    assert dlg.result() == QDialog.Accepted and dlg.selected is not None


def test_long_drop_down_entries_are_not_cut_off(app):
    """The outlier method next to the k box read "Median +/- k M"; the garnet tracer
    laws ran off the end of their box; the editor's "No marker" lost its last letter."""
    from PySide6.QtWidgets import QComboBox
    from diffusor.gui.image_extractor import ImageExtractorDialog
    from diffusor.gui.multicomponent_study import _family_text
    from diffusor.gui.graph_editor import MARKERS, _combo
    assert _family_text("Garnet Fe, Mg, Mn, Ca tracer laws, Carlson (2006)") == \
        "Fe, Mg, Mn, Ca, Carlson (2006)"
    assert _family_text("Garnet tracer laws, Carlson (2006) with the Mn law of Chen & Chu "
                        "(2024)") == "Carlson (2006) with the Mn law of Chen & Chu (2024)"
    dlg = ImageExtractorDialog()
    assert dlg.cmb_clip.sizeAdjustPolicy() != QComboBox.AdjustToMinimumContentsLengthWithIcon
    dlg.close()
    cmb = _combo(MARKERS, "None")
    assert cmb.minimumWidth() >= cmb.sizeHint().width() - 1
