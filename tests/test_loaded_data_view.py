"""What the Loaded data card says about the points, and the y label staying inside the plot."""
import os
import time

import numpy as np
import pytest

from diffusor import datasets as ds
from diffusor.dataio import ProfileSpec, build_profile, read_table


def _kizimen():
    d = {x.key: x for x in ds.DATASETS}["opx_kizimen"]
    return d, build_profile(read_table(d.path), ProfileSpec(**d.spec), str(d.path))


def test_a_profile_knows_how_many_points_the_fit_window_left_out():
    d, p = _kizimen()
    assert len(p) == 28 and p.n_outside_window == 81 and len(p.raw) == 109
    assert p.sorted().n_outside_window == 81 and p.kept().n_outside_window == 81
    assert p.spec.window_text() == "4 to 60 um"
    whole = build_profile(read_table(d.path), ProfileSpec(**{**d.spec, "x_min": None, "x_max": None}), "x")
    assert len(whole) == 109 and whole.n_outside_window == 0 and whole.spec.window_text() == ""


def test_the_points_phrase_gives_both_counts_and_the_window():
    from diffusor.gui.main_window import _points_phrase
    _, p = _kizimen()
    text = _points_phrase(p)
    assert text.startswith("28 of 109 points, from 4.2 to 59.5 um.")
    assert "fit window (4 to 60 um) leaves out the other 81" in text
    p.n_outside_window = 0
    assert _points_phrase(p) == "28 points, from 4.2 to 59.5 um."


# --- the plot ----------------------------------------------------------------------
@pytest.fixture(scope="module")
def app():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    import matplotlib
    matplotlib.use("QtAgg")
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _settle(app, n=8):
    for _ in range(n):
        app.processEvents()
        time.sleep(0.01)


def test_a_y_label_that_sticks_out_is_redrawn_once_and_never_in_a_loop(app):
    from diffusor.gui.plot_widget import DataPreview
    p = DataPreview()
    p.resize(500, 220)
    p.show()
    _settle(app)
    draws = []
    p.canvas.mpl_connect("draw_event", lambda e: draws.append(1))
    x = np.linspace(0, 10, 12)
    p.show_data(x, x * 0.01, y_label="X_Fe (Fe-Mg)")
    _settle(app)
    assert 1 <= len(draws) <= 2, "a label inside the figure needs no extra draw"
    # a label pinned outside the figure, which no layout can bring back
    draws.clear()
    p.figure.axes[0].yaxis.set_label_coords(-0.5, 0.5)
    p.canvas.draw_idle()
    _settle(app, 40)
    assert 1 <= len(draws) <= 4, f"{len(draws)} draws: the guard must not loop"
    p.close()


def test_the_kizimen_text_explains_source_part_used_and_result_in_three_paragraphs():
    from diffusor.gui import richtext
    d, _ = _kizimen()
    paras = richtext.paragraphs(d.provenance)
    assert len(paras) == 3
    source, used, result = paras
    assert "K9_L10C4" in source and "Supplementary Data 2" in source and "zenodo.7307563" in source
    assert "109 analyses" in used and "28 of them" in used and "4.2 to 59.5 um" in used
    assert "80 points beyond 60 um" in used
    assert "2.32 years" in result and "about 2.5 years" in result
    assert "Real measured" not in d.provenance and "real data" not in d.provenance.lower()
    html = richtext.example_html(d)
    assert html.count("<p>") >= 4 and "Where the data come from" in html
