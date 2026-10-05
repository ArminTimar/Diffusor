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
    from diffusor.coefficients.latex import coefficient_latex, equation_pieces, pixel_ratio
    from diffusor.coefficients import get as get_coefficient
    w = _window(app)
    _load_example(w, app)
    w._go(2)
    view = w.lst_coef
    heights = {}
    for text in ("Dias, Dohmen & Behrens (2025)", "Ganguly & Tazzoli (1994) / Allan"):
        _select(w, text)
        app.processEvents()
        eq = w.lbl_coef_eq
        key = w.lst_coef.currentItem().data(Qt.UserRole)
        pieces = equation_pieces(coefficient_latex(get_coefficient(key)),
                                 max(eq.width() - 2, eq.MIN_WIDTH), pixel_ratio())
        expected = sum(h for _png, _w, h in pieces) + eq.GAP * (len(pieces) - 1)
        assert eq.height() == expected, f"{text}: {eq.height()} against {expected}"
        heights[text] = eq.height()
    # the shorter law is not padded out to the height of the longer one
    assert heights["Ganguly & Tazzoli (1994) / Allan"] < heights["Dias, Dohmen & Behrens (2025)"]
    w.close()


def test_the_equation_is_cut_again_when_the_window_narrows(app):
    w = _window(app, width=1250)
    _load_example(w, app)
    w._go(2)
    _select(w, "Dias, Dohmen & Behrens (2025)")
    app.processEvents()
    wide = w.lbl_coef_eq.height()
    w.resize(1140, 620)
    app.processEvents()
    w.lbl_coef_eq._relayout()
    assert w.lbl_coef_eq.height() >= wide
    w.close()
