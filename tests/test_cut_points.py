"""Cutting points out of the fitted profile: the mask, the plot gestures and what the fit sees."""
from types import SimpleNamespace

import os

import numpy as np
import pytest
from PySide6.QtCore import Qt

from diffusor.dataio.profiles import Profile


def _profile(n=8):
    x = np.arange(n, dtype=float)[::-1]          # unsorted on purpose
    return Profile(x=x, C=x * 0.1, sigma=np.full(n, 0.01), source="a.csv").sorted()


def test_kept_drops_the_cut_points_and_says_so():
    p = _profile()
    assert p.kept() is p and p.n_excluded == 0
    mask = np.zeros(len(p), bool)
    mask[[2, 5]] = True
    p.set_excluded(mask)
    k = p.kept()
    assert len(k) == len(p) - 2 and p.n_excluded == 2
    assert not np.isin(p.x[mask], k.x).any()
    assert k.sigma.size == k.x.size
    assert any("cut by hand" in n for n in k.notes)
    assert len(p) == 8, "the cut points stay in the loaded profile"


def test_an_empty_mask_means_nothing_is_cut_and_a_bad_one_is_refused():
    p = _profile()
    p.set_excluded(np.zeros(len(p), bool))
    assert p.excluded is None
    with pytest.raises(ValueError):
        p.set_excluded(np.zeros(3, bool))


def test_the_cut_follows_the_points_when_the_profile_is_sorted():
    p = Profile(x=np.array([3.0, 1.0, 2.0]), C=np.array([30.0, 10.0, 20.0]))
    p.excluded = np.array([True, False, False])      # the point at x = 3
    q = p.sorted()
    assert q.x[q.excluded].tolist() == [3.0]


# --- the window ------------------------------------------------------------------
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


def _open_example(w, app, key="opx_shinmoedake"):
    for i in range(w.lst_examples.count()):
        if w.lst_examples.item(i).data(Qt.UserRole) == key:
            w.lst_examples.setCurrentRow(i)
            break
    w.load_example()
    w._go(6)
    w.plot.canvas.draw()
    app.processEvents()


def _event(plot, x, y, key=None, inside=True):
    px, py = plot.ax.transData.transform((x, y))
    return SimpleNamespace(x=px, y=py, xdata=x, ydata=y, button=1, key=key,
                           inaxes=plot.ax if inside else None)


def _drag(plot, a, b, key=None):
    plot._cut_press(_event(plot, *a, key=key))
    plot._cut_motion(_event(plot, *b, key=key))
    plot._cut_release(_event(plot, *b, key=key))


def test_a_click_cuts_the_nearest_point_and_a_second_click_restores_it(app):
    w = _window(app)
    _open_example(w, app)
    p, plot = w.profile, w.plot
    n = len(p)
    plot.act_cut.setChecked(True)
    i = 4
    plot._cut_press(_event(plot, p.x[i], p.C[i]))
    plot._cut_release(_event(plot, p.x[i], p.C[i]))
    assert p.n_excluded == 1 and p.excluded[i]
    assert len(w._fit_profile()) == n - 1
    assert w.plot._pts[2][i], "the plot marks the point as cut"
    plot._cut_press(_event(plot, p.x[i], p.C[i]))
    plot._cut_release(_event(plot, p.x[i], p.C[i]))
    assert p.n_excluded == 0 and p.excluded is None
    w.close()


def test_a_box_cuts_every_point_inside_and_shift_restores_them(app):
    w = _window(app)
    _open_example(w, app)
    p, plot = w.profile, w.plot
    plot.act_cut.setChecked(True)
    lo, hi = p.x[3] - 1e-9, p.x[9] + 1e-9
    ymin, ymax = float(p.C.min()), float(p.C.max())
    pad = 0.05 * (ymax - ymin)
    _drag(plot, (lo, ymin - pad), (hi, ymax + pad))
    assert p.n_excluded == 7 and not p.excluded[:3].any() and p.excluded[3:10].all()
    assert len(w._fit_profile()) == len(p) - 7
    _drag(plot, (lo, ymin - pad), (hi, ymax + pad), key="shift")
    assert p.n_excluded == 0
    w.close()


def test_gestures_do_nothing_unless_cut_mode_is_on(app):
    w = _window(app)
    _open_example(w, app)
    p, plot = w.profile, w.plot
    plot._cut_press(_event(plot, p.x[2], p.C[2]))
    plot._cut_release(_event(plot, p.x[2], p.C[2]))
    assert p.n_excluded == 0
    w.close()


def test_a_click_far_from_every_point_cuts_nothing(app):
    w = _window(app)
    _open_example(w, app)
    p, plot = w.profile, w.plot
    plot.act_cut.setChecked(True)
    mid = 0.5 * (p.x[0] + p.x[1])
    far_y = float(p.C.max()) + 3 * (float(p.C.max()) - float(p.C.min()))
    plot.ax.set_ylim(float(p.C.min()), far_y)
    plot.canvas.draw()
    plot._cut_press(_event(plot, mid, far_y * 0.999))
    plot._cut_release(_event(plot, mid, far_y * 0.999))
    assert p.n_excluded == 0
    w.close()


def test_at_least_three_points_must_stay_in_the_fit(app):
    w = _window(app)
    _open_example(w, app)
    p = w.profile
    mask = np.ones(len(p), bool)
    mask[:2] = False
    w._on_points_cut(mask)
    assert p.n_excluded == 0, "leaving two points is refused"
    mask[:3] = False
    w._on_points_cut(mask)
    assert len(w._fit_profile()) == 3
    w.close()


def test_cutting_discards_a_result_made_on_the_old_points_and_refits(app, monkeypatch):
    from diffusor.fitting import fit_time
    w = _window(app)
    _open_example(w, app)
    model = w._model(w._checked_keys()[0])
    w.fit_result = fit_time(model, w.profile.x, w.profile.C, w.profile.sigma,
                            w._free_parameters())
    called = []
    monkeypatch.setattr(w, "run_fit", lambda: called.append(True))
    mask = np.zeros(len(w.profile), bool)
    mask[0] = True
    w._on_points_cut(mask)
    assert w.fit_result is None and called == [True]
    w.close()


def test_the_fit_uses_the_kept_points_only(app):
    from diffusor.fitting import fit_time
    w = _window(app)
    _open_example(w, app)
    p = w.profile
    bad = len(p) // 2
    C = p.C.copy()
    C[bad] += 5 * (C.max() - C.min())               # a wild point on the step
    p.C = C
    model = w._model(w._checked_keys()[0])
    full = fit_time(model, p.x, p.C, p.sigma, w._free_parameters())
    mask = np.zeros(len(p), bool)
    mask[bad] = True
    w._on_points_cut(mask)
    k = w._fit_profile()
    assert not np.isclose(k.C, C[bad]).any()
    cut = fit_time(model, k.x, k.C, k.sigma, w._free_parameters())
    assert cut.x_data.size == len(p) - 1
    assert cut.stats.reduced_chi2 < full.stats.reduced_chi2
    w.close()


def test_a_new_file_starts_with_every_point_in_the_fit(app):
    w = _window(app)
    _open_example(w, app)
    mask = np.zeros(len(w.profile), bool)
    mask[1] = True
    w._on_points_cut(mask)
    assert w.profile.n_excluded == 1
    assert "1 point cut." in w.plot.act_cut.toolTip()
    _open_example(w, app, "opx_shinmoedake")
    assert w.profile.n_excluded == 0 and " cut. " not in w.plot.act_cut.toolTip()
    w.close()
