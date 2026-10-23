"""Window sizing and the typeset equation of the coefficient step."""
import os

import pytest
from PySide6.QtCore import QPoint, Qt


@pytest.fixture(scope="module")
def app():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    import matplotlib
    matplotlib.use("QtAgg")
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _window(app, width=1250, height=620):
    from diffusor.gui.main_window import MainWindow
    w = MainWindow()
    w.resize(width, height)
    w.show()
    app.processEvents()
    return w


def _load_example(w, app, key="opx_shinmoedake"):
    for i in range(w.lst_examples.count()):
        if w.lst_examples.item(i).data(Qt.UserRole) == key:
            w.lst_examples.setCurrentRow(i)
            break
    w.load_example()
    app.processEvents()


def _select(w, text):
    for i in range(w.lst_coef.count()):
        if text in w.lst_coef.item(i).text():
            w.lst_coef.setCurrentRow(i)
            return
    raise AssertionError(text)


def _settled(app, view, seconds=3.0):
    """Wait until the equation is laid out for the width the view has now.

    The view is cut for a stand-in width before it is first shown and again 90 ms
    after its width last changed, so its height is final only once that timer has
    fired and no new resize has followed. Font metrics differ between systems, so
    whether a height measured earlier happens to be right differs too.
    """
    import time
    end = time.monotonic() + seconds
    while time.monotonic() < end:
        app.processEvents()
        if not view._settle.isActive() and view._laid_width == view.width():
            return
        time.sleep(0.01)
    raise AssertionError(f"the equation was not laid out again for width {view.width()}")


def _expected_height(view, latex):
    from diffusor.coefficients.latex import equation_pieces, pixel_ratio
    pieces = equation_pieces(latex, max(view.width() - 2, view.MIN_WIDTH), pixel_ratio())
    return sum(h for _png, _w, h in pieces) + view.GAP * (len(pieces) - 1)


def test_a_long_equation_does_not_raise_the_window_minimum(app):
    """The longest law once made the window ask for more than a laptop screen has."""
    w = _window(app)
    _load_example(w, app)
    w._go(2)
    heights = {}
    for text in ("Ganguly & Tazzoli (1994) / Allan", "Dias, Dohmen & Behrens (2025)"):
        _select(w, text)
        app.processEvents()
        heights[text] = w.minimumSizeHint().height()
    assert len(set(heights.values())) == 1, heights
    # 672 px is what a 1920 x 1080 screen at 150 % scaling leaves, less the title bar
    assert max(heights.values()) <= 640, heights
    w.close()


def test_the_buttons_stay_in_view_with_the_longest_equation(app):
    w = _window(app, height=620)
    _load_example(w, app)
    w._go(2)
    _select(w, "Dias, Dohmen & Behrens (2025)")
    app.processEvents()
    assert w.height() <= 625, "the window grew to fit the page"
    bottom = w.btn_next.mapTo(w, QPoint(0, 0)).y() + w.btn_next.height()
    assert bottom <= w.height()
    w.close()


def test_the_equation_card_scrolls_in_a_short_window_instead_of_growing(app):
    from PySide6.QtWidgets import QScrollArea
    w = _window(app, height=620)
    _load_example(w, app)
    w._go(2)
    _select(w, "Dias, Dohmen & Behrens (2025)")
    app.processEvents()
    area = w.lbl_coef_eq.parentWidget().parentWidget().parentWidget()
    assert isinstance(area, QScrollArea)
    assert area.minimumSizeHint().height() <= 140
    w.close()


def test_a_wrapped_label_gives_back_height_once_it_is_wide(app):
    """A label first laid out narrow once kept the tall minimum it needed then."""
    from PySide6.QtWidgets import QVBoxLayout, QWidget
    from diffusor.gui.widgets import WrapLabel
    host = QWidget()
    lay = QVBoxLayout(host)
    label = WrapLabel("A sentence long enough to need several lines in a narrow column " * 3)
    lay.addWidget(label)
    host.resize(120, 600)
    host.show()
    app.processEvents()
    narrow = label.minimumHeight()
    host.resize(900, 600)
    app.processEvents()
    assert label.minimumHeight() < narrow / 3, (narrow, label.minimumHeight())
    host.close()


def test_the_equation_is_exactly_as_tall_as_its_lines(app):
    from diffusor.coefficients.latex import coefficient_latex
    from diffusor.coefficients import get as get_coefficient
    w = _window(app)
    _load_example(w, app)
    w._go(2)
    heights = {}
    for text in ("Dias, Dohmen & Behrens (2025)", "Ganguly & Tazzoli (1994) / Allan"):
        _select(w, text)
        eq = w.lbl_coef_eq
        _settled(app, eq)
        key = w.lst_coef.currentItem().data(Qt.UserRole)
        expected = _expected_height(eq, coefficient_latex(get_coefficient(key)))
        assert eq.height() == expected, f"{text}: {eq.height()} against {expected}"
        heights[text] = eq.height()
    # the shorter law is not padded out to the height of the longer one
    assert heights["Ganguly & Tazzoli (1994) / Allan"] < heights["Dias, Dohmen & Behrens (2025)"]
    w.close()


def test_the_equation_is_cut_again_when_the_view_narrows(app):
    """A narrower view gets lines cut for its new width, and the height that goes with them.

    The height need not grow: a line with nowhere to cut is drawn smaller instead,
    so a narrower view can be a few pixels shorter. The view sits in a plain host
    here, because the main window's smallest width depends on the system's fonts.
    """
    from PySide6.QtWidgets import QVBoxLayout, QWidget
    from diffusor.coefficients.latex import coefficient_latex
    from diffusor.coefficients import get as get_coefficient
    from diffusor.gui.widgets import EquationView
    latex = coefficient_latex(get_coefficient("opx_FeMg_dias2025"))
    host = QWidget()
    lay = QVBoxLayout(host)
    view = EquationView()
    lay.addWidget(view)
    lay.addStretch(1)
    view.set_latex(latex)
    host.resize(900, 700)
    host.show()
    _settled(app, view)
    wide = (view.width(), view.height())
    host.resize(300, 700)
    _settled(app, view)
    assert view.width() < wide[0]
    assert view.height() == _expected_height(view, latex)
    assert view.height() != wide[1], "the lines were not cut again for the narrow width"
    host.close()
