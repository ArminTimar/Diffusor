"""Editable vector exports (Inkscape, CorelDRAW) and the plot toolbar that saves them."""
import io
import os
import re
import xml.dom.minidom

import numpy as np
import pytest
from matplotlib.figure import Figure

from diffusor.dataio.vector import editable_copy, save_vector


def _figure():
    fig = Figure(figsize=(6, 4), layout="constrained")
    ax = fig.add_subplot(111)
    x = np.linspace(0, 30, 12)
    ax.errorbar(x, np.tanh(x / 10), yerr=0.05, fmt="o", ms=4, capsize=2, label="measured")
    ax.plot(x, np.tanh(x / 10) + 0.02, "-", label="fit")
    ax.plot(x[[2, 3]], np.tanh(x[[2, 3]] / 10), "x", label="cut from the fit")
    ax.set_xlabel("distance (um)")
    ax.legend()
    return fig


def test_svg_keeps_text_and_makes_every_point_its_own_group(tmp_path):
    path = save_vector(_figure(), tmp_path / "f.svg")
    svg = open(path, encoding="utf-8").read()
    xml.dom.minidom.parseString(svg.encode("utf-8"))               # still well-formed
    assert "distance (um)" in svg and "<text" in svg, "text was turned into outlines"
    for i in range(1, 13):
        group = re.search(rf'<g id="measured_point_{i:03d}">(.*?)</g>\s*(?=<g id=|</g>)', svg, re.S)
        assert group, f"point {i} has no group"
    assert len(re.findall(r'<g id="measured_point_\d+">', svg)) == 12
    assert len(re.findall(r'<g id="cut-from-the-fit_point_\d+">', svg)) == 2
    # a point carries its error bar and caps in the same group
    first = svg[svg.index('<g id="measured_point_001">'):svg.index('<g id="measured_point_002">')]
    assert "measured_point_001_marker" in first and "measured_point_001_errorbar" in first
    assert 'id="legend"' in svg and 'id="fit"' in svg


def test_markers_are_plain_paths_not_shared_copies(tmp_path):
    svg = open(save_vector(_figure(), tmp_path / "f.svg"), encoding="utf-8").read()
    marker = re.search(r'<g id="measured_point_001_marker">(.*?)</g>', svg, re.S).group(1)
    assert "<path" in marker and "<use" not in marker


def test_the_figure_on_screen_is_not_changed(tmp_path):
    fig = _figure()
    before = [(len(a.lines), len(a.collections)) for a in fig.axes]
    save_vector(fig, tmp_path / "f.svg")
    assert [(len(a.lines), len(a.collections)) for a in fig.axes] == before


def test_split_points_look_like_the_original():
    from PIL import Image
    fig = _figure()

    def render(f):
        buf = io.BytesIO()
        f.savefig(buf, format="png", dpi=100)
        return np.asarray(Image.open(io.BytesIO(buf.getvalue())).convert("RGB"), float)
    a, b = render(fig), render(editable_copy(fig))
    assert a.shape == b.shape and np.abs(a - b).mean() < 0.5


def test_pdf_and_eps_are_written_and_other_suffixes_refused(tmp_path):
    assert os.path.getsize(save_vector(_figure(), tmp_path / "f.pdf")) > 1000
    assert os.path.getsize(save_vector(_figure(), tmp_path / "f.eps")) > 1000
    with pytest.raises(ValueError):
        save_vector(_figure(), tmp_path / "f.png")


# --- the toolbar --------------------------------------------------------------------
@pytest.fixture(scope="module")
def app():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    import matplotlib
    matplotlib.use("QtAgg")
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _plot(app):
    from diffusor.gui.plot_widget import ProfilePlot
    p = ProfilePlot()
    p.resize(900, 400)
    p.show()
    x = np.linspace(0, 30, 12)
    p.set_points(x, np.tanh(x / 10), None)
    p.show_data(x, np.tanh(x / 10), np.full_like(x, 0.05))
    app.processEvents()
    return p


def test_only_one_tool_is_active_and_the_active_one_is_marked(app):
    p = _plot(app)
    assert "QToolButton:checked" in p.toolbar.styleSheet()
    p.act_cut.setChecked(True)
    assert p.act_cut.isChecked()
    p.toolbar.zoom()
    assert not p.act_cut.isChecked() and p.toolbar._actions["zoom"].isChecked()
    p.act_cut.setChecked(True)
    assert not p.toolbar._actions["zoom"].isChecked() and not p.toolbar.mode
    p.toolbar.pan()
    assert not p.act_cut.isChecked()
    p.close()


def test_cut_and_restore_have_icons_like_the_other_buttons(app):
    p = _plot(app)
    for act in (p.act_cut, p.act_restore):
        assert not act.icon().isNull()
    p.close()


def test_the_legend_has_a_frame_so_it_is_not_taken_for_a_point(app):
    p = _plot(app)
    leg = p.ax.get_legend()
    assert leg is not None and leg.get_frame_on()
    assert leg.get_frame().get_alpha() > 0.8
    p.close()


def test_saving_svg_from_the_toolbar_writes_the_editable_form(app, tmp_path, monkeypatch):
    from matplotlib.backends import backend_qt
    p = _plot(app)
    target = tmp_path / "toolbar.svg"

    class Dialog:          # the real dialog is modal and would wait for a person
        @staticmethod
        def getSaveFileName(*_args, **_kwargs):
            return str(target), "Scalable Vector Graphics (*.svg)"
    monkeypatch.setattr(backend_qt.QtWidgets, "QFileDialog", Dialog)
    p.toolbar.save_figure()
    svg = target.read_text(encoding="utf-8")
    assert 'id="measured_point_001"' in svg and "<text" in svg
    assert "savefig" not in p.figure.__dict__, "the override must not outlive the save"
    p.close()
