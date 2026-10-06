"""Export: the user chooses which outputs are written, before choosing the folder."""
import os

import numpy as np
import pytest

from diffusor.constants import SEC_PER_YEAR


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
    mc = run_montecarlo(m, x, m.profile(SEC_PER_YEAR, x), np.full_like(x, 0.003), n_draws=12,
                        seed=1, budget=UncertaintyBudget(sigma_T_K=20.0, buffer="NNO",
                                                         delta_buffer=1.0))
    return r, mc


def test_only_the_chosen_outputs_are_written(run, tmp_path):
    from diffusor.dataio import save_results
    r, mc = run
    written = save_results(tmp_path, r, mc, include=["json", "methods"], basename="t")
    assert set(written) == {"json", "methods"}
    assert sorted(p.name for p in tmp_path.iterdir()) == ["t_methods.txt", "t_results.json"]


def test_the_workbook_and_figures_do_not_need_the_text_files(run, tmp_path):
    from diffusor.dataio import save_results
    r, mc = run
    written = save_results(tmp_path, r, mc, include=["workbook", "figure_svg"], basename="t")
    assert set(written) == {"workbook", "figure_svg"}
    assert sorted(p.name for p in tmp_path.iterdir()) == ["t_figure.svg", "t_results.xlsx"]


def test_no_include_writes_everything_available_and_skips_monte_carlo_without_one(run, tmp_path):
    from diffusor.dataio import save_results
    from diffusor.dataio.export import export_keys
    r, mc = run
    assert set(save_results(tmp_path / "a", r, mc)) == set(export_keys(True))
    plain = save_results(tmp_path / "b", r, None)
    assert set(plain) == set(export_keys(False)) and "montecarlo" not in plain
    assert "montecarlo" not in save_results(tmp_path / "c", r, None, include=["montecarlo", "json"])


def test_an_unknown_output_is_refused(run, tmp_path):
    from diffusor.dataio import save_results
    r, _ = run
    with pytest.raises(ValueError, match="unknown export items"):
        save_results(tmp_path, r, include=["json", "nonsense"])


# --- the dialog ----------------------------------------------------------------------
@pytest.fixture(scope="module")
def app():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def test_the_dialog_lists_every_output_and_ticks_what_is_available(app):
    from diffusor.dataio.export import EXPORT_ITEMS
    from diffusor.gui.export_dialog import ExportDialog
    dlg = ExportDialog(None, has_monte_carlo=False)
    assert list(dlg._boxes) == [k for k, *_ in EXPORT_ITEMS]
    assert not dlg._boxes["montecarlo"].isEnabled() and not dlg._boxes["montecarlo"].isChecked()
    assert dlg.selected() == [k for k, *_ in EXPORT_ITEMS if k != "montecarlo"]
    with_mc = ExportDialog(None, has_monte_carlo=True)
    assert "montecarlo" in with_mc.selected()


def test_the_dialog_remembers_choices_and_cannot_continue_with_nothing_ticked(app):
    from diffusor.gui.export_dialog import ExportDialog
    dlg = ExportDialog(None, has_monte_carlo=True, chosen=["json", "figure"])
    assert dlg.selected() == ["json", "figure"] and dlg.btn_ok.isEnabled()
    dlg._set_all(False)
    assert dlg.selected() == [] and not dlg.btn_ok.isEnabled()
    dlg._set_all(True)
    assert len(dlg.selected()) == 8 and dlg.btn_ok.isEnabled()


def test_about_names_the_author(app):
    import matplotlib
    matplotlib.use("QtAgg")
    from diffusor.gui import main_window as mw
    shown = []
    original = mw.QMessageBox.about
    mw.QMessageBox.about = staticmethod(lambda parent, title, text: shown.append((title, text)))
    try:
        w = mw.MainWindow()
        w.show_about()
        w.close()
    finally:
        mw.QMessageBox.about = original
    title, text = shown[0]
    assert title == "About Diffusor" and "Created by Armin Timar" in text and "MIT" in text
