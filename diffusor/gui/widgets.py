"""Small layout helpers shared by the step pages."""
from __future__ import annotations

from typing import Optional

from PySide6.QtCore import QEvent, QObject, QSize, Qt, QTimer
from PySide6.QtGui import QGuiApplication, QPainter, QPixmap
from PySide6.QtWidgets import (QAbstractScrollArea, QAbstractSpinBox, QApplication, QComboBox,
                               QFrame, QHBoxLayout, QLabel, QPushButton, QScrollArea,
                               QSizePolicy, QVBoxLayout, QWidget)

from . import theme


class WrapLabel(QLabel):
    """A word-wrapped label that actually reserves the height it needs.

    A plain wrapped QLabel reports a one-line minimum height, so a vertical
    layout hands it one line and the rest of the text paints over whatever sits
    below it. Reserving heightForWidth on every resize is the standard cure.
    """

    def __init__(self, text: str = "", parent=None):
        super().__init__(text, parent)
        self.setWordWrap(True)
        self.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Minimum)
        self.setAlignment(Qt.AlignTop | Qt.AlignLeft)

    def sizeHint(self):
        # A wrapped QLabel guesses its height for a narrow, golden-ratio width and
        # asks for far more than it needs. Once it has a real width, use that.
        s = super().sizeHint()
        if self.width() > 0:
            s.setHeight(self._height_for(self.width()))
        return s

    def _height_for(self, width: int) -> int:
        # Qt answers -1 ("no preference") for a label with no text yet; a minimum
        # height of -1 is refused with a "Negative sizes" warning, so use 0
        return max(self.heightForWidth(width), 0)

    def _sync(self):
        w = self.width()
        if w > 0:
            old = self.minimumHeight()
            if old:
                # QLabel never answers heightForWidth with less than its own minimum
                # height, so measure with it cleared: otherwise a label first laid out
                # narrow keeps the tall height it needed then, even once it is wide
                self.setMinimumHeight(0)
            h = self._height_for(w)
            if h != old:
                self.setMinimumHeight(h)
                self.updateGeometry()
            elif old:
                self.setMinimumHeight(old)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._sync()

    def setText(self, text):
        super().setText(text)
        self._sync()
        self.updateGeometry()


class EquationView(QWidget):
    """A typeset equation, drawn line by line from pictures, exactly as tall as it is.

    The equation used to sit in a rich-text label, whose line boxes add room above
    and below every picture and whose height was measured before the label had its
    final width. Here each line is a picture drawn at the screen's own resolution,
    the lines are cut to the width the widget really has, and the widget's height is
    the sum of the pictures, so nothing is left over below it.
    """
    GAP = 5             # logical pixels between two lines
    MIN_WIDTH = 160     # the narrowest the lines are cut for

    def __init__(self, parent=None):
        super().__init__(parent)
        self._latex = ""
        self._pixmaps: list = []
        self._laid_width = 0
        self.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        self.setFixedHeight(0)
        # redraw once a resize has settled, not at every step of a drag
        self._settle = QTimer(self)
        self._settle.setSingleShot(True)
        self._settle.setInterval(90)
        self._settle.timeout.connect(self._relayout)

    def set_latex(self, latex: str):
        self._latex = latex or ""
        self._relayout()

    def _relayout(self):
        from ..coefficients.latex import equation_pieces, pixel_ratio
        width = max(self.width() - 2, self.MIN_WIDTH) if self.width() > 50 else 560
        ratio = pixel_ratio()
        pixmaps, height = [], 0
        for png, _w, h in (equation_pieces(self._latex, width, ratio) if self._latex else []):
            pm = QPixmap()
            pm.loadFromData(png)
            pm.setDevicePixelRatio(ratio)
            pixmaps.append((pm, height))
            height += h + self.GAP
        self._pixmaps = pixmaps
        self._laid_width = self.width()
        self.setFixedHeight(max(height - self.GAP, 0))
        self.update()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        if self.width() != self._laid_width:
            self._settle.start()

    def paintEvent(self, _event):
        p = QPainter(self)
        for pm, y in self._pixmaps:
            p.drawPixmap(0, y, pm)


def card(title: str = "", subtitle: str = "") -> tuple[QFrame, QVBoxLayout]:
    """A white rounded panel with an optional heading; returns (frame, body layout)."""
    frame = QFrame()
    frame.setObjectName("Card")
    outer = QVBoxLayout(frame)
    outer.setContentsMargins(16, 12, 16, 14)
    outer.setSpacing(8)
    if title:
        lab = QLabel(title)
        lab.setObjectName("H2")
        outer.addWidget(lab)
    if subtitle:
        sub = WrapLabel(subtitle)
        sub.setObjectName("Sub")
        outer.addWidget(sub)
    body = QVBoxLayout()
    body.setSpacing(8)
    outer.addLayout(body, 1)
    return frame, body


def field(label: str, widget: QWidget, hint: str = "") -> QWidget:
    """A labelled control laid out vertically, so long labels never clip."""
    w = QWidget()
    v = QVBoxLayout(w)
    v.setContentsMargins(0, 0, 0, 0)
    v.setSpacing(3)
    lab = QLabel(label)
    lab.setObjectName("FieldLabel")
    v.addWidget(lab)
    v.addWidget(widget)
    if hint:
        h = WrapLabel(hint)
        h.setObjectName("Hint")
        v.addWidget(h)
    return w


def row(*widgets: QWidget, spacing: int = 12, stretch_last: bool = False) -> QWidget:
    w = QWidget()
    h = QHBoxLayout(w)
    h.setContentsMargins(0, 0, 0, 0)
    h.setSpacing(spacing)
    for i, x in enumerate(widgets):
        h.addWidget(x, 1 if (stretch_last and i == len(widgets) - 1) else 0)
    if not stretch_last:
        h.addStretch(1)
    return w


def pair(a: QWidget, b: QWidget) -> QWidget:
    """Two controls side by side, each taking half the width."""
    w = QWidget()
    h = QHBoxLayout(w)
    h.setContentsMargins(0, 0, 0, 0)
    h.setSpacing(12)
    h.addWidget(a, 1)
    h.addWidget(b, 1)
    return w


def divider() -> QFrame:
    f = QFrame()
    f.setObjectName("Divider")
    f.setFixedHeight(1)
    return f


def note(text: str, kind: str = "Hint") -> QLabel:
    lab = WrapLabel(text)
    lab.setObjectName(kind)
    return lab


def callout(kind: str = "Info") -> QFrame:
    """A tinted box holding a wrapped label, available as ``box.label``.

    The padding lives on the frame. Padding set on a wrapped QLabel itself makes
    Qt overestimate the label's height and leaves an empty band under the text.
    """
    box = QFrame()
    box.setObjectName(kind + "Box")
    lay = QVBoxLayout(box)
    lay.setContentsMargins(11, 8, 11, 8)
    lab = WrapLabel("")
    lab.setObjectName("CalloutText")
    lay.addWidget(lab)
    box.label = lab
    return box


def fit_to_screen(window: QWidget, width: int, height: int,
                  width_share: float = 0.95, height_share: float = 0.88) -> None:
    """Resize to ``width`` x ``height``, or less so the whole window fits the screen.

    ``height_share`` leaves room for the title bar and the taskbar, so the
    buttons at the bottom of a dialog are always on screen.
    """
    screen = (window.screen() if window.screen() is not None
              else QGuiApplication.primaryScreen()).availableGeometry()
    window.resize(min(width, int(screen.width() * width_share)),
                  min(height, int(screen.height() * height_share)))


class FitScrollArea(QScrollArea):
    """A scroll area that asks for the whole height of its content, and can give it up.

    A plain QScrollArea asks for at most about 24 lines of text, so a tall panel
    would scroll even in a large window. This one asks for its content's own height,
    which a layout grants when there is room, and shrinks below it when there is not,
    so a long equation scrolls inside its card instead of making the window taller
    than the screen.
    """

    def sizeHint(self):
        inner = self.widget()
        if inner is None:
            return super().sizeHint()
        s = inner.sizeHint()
        f = 2 * self.frameWidth()
        return QSize(s.width() + f, s.height() + f)

    def minimumSizeHint(self):
        return QSize(120, 120)


def scrollable(inner: QWidget) -> QScrollArea:
    area = QScrollArea()
    area.setWidget(inner)
    area.setWidgetResizable(True)
    area.setFrameShape(QFrame.NoFrame)
    area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
    return area


def page_body(*widgets: QWidget, max_width: int = 720) -> QScrollArea:
    """A scrolling, width-capped column of cards, centred in the page."""
    inner = QWidget()
    inner.setObjectName("Page")
    outer = QHBoxLayout(inner)
    outer.setContentsMargins(24, 20, 24, 24)
    outer.addStretch(1)
    col = QWidget()
    col.setMaximumWidth(max_width)
    col.setMinimumWidth(420)
    col.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.MinimumExpanding)
    v = QVBoxLayout(col)
    v.setContentsMargins(0, 0, 0, 0)
    v.setSpacing(14)
    for x in widgets:
        v.addWidget(x)
    v.addStretch(1)
    outer.addWidget(col, 10)
    outer.addStretch(1)
    return scrollable(inner)


def page_columns(left, right, top=(), max_width: int = 1200) -> QWidget:
    """A page that never scrolls: an optional full-width top row and two columns.

    ``left`` and ``right`` hold widgets or ``(widget, stretch)`` pairs. A column
    without a stretching widget gets its spare height at the bottom. Only lists
    and text panes inside the cards scroll.
    """
    page = QWidget()
    page.setObjectName("Page")
    outer = QHBoxLayout(page)
    outer.setContentsMargins(24, 10, 24, 10)
    outer.addStretch(1)
    col = QWidget()
    col.setMaximumWidth(max_width)
    v = QVBoxLayout(col)
    v.setContentsMargins(0, 0, 0, 0)
    v.setSpacing(12)
    for w in top:
        v.addWidget(w)
    h = QHBoxLayout()
    h.setSpacing(14)
    for items in (left, right):
        c = QVBoxLayout()
        c.setSpacing(12)
        stretched = False
        for it in items:
            w, s = it if isinstance(it, tuple) else (it, 0)
            c.addWidget(w, s)
            stretched = stretched or s > 0
        if not stretched:
            c.addStretch(1)
        h.addLayout(c, 1)
    v.addLayout(h, 1)
    outer.addWidget(col, 20)
    outer.addStretch(1)
    return page


def badge(text: str, kind: str = "muted") -> QLabel:
    lab = QLabel(text)
    colours = {
        "measured": (theme.OK, "#E8F2EC"),
        "synthetic": (theme.WARN, theme.WARN_SOFT),
        "muted": (theme.TEXT_MUTED, theme.SURFACE_ALT),
        "accent": (theme.ACCENT, theme.ACCENT_SOFT),
        "danger": (theme.DANGER, "#F8EBEB"),
    }
    fg, bg = colours.get(kind, colours["muted"])
    lab.setStyleSheet(
        f"color:{fg}; background:{bg}; border-radius:9px; padding:2px 9px; font-size:11px;")
    lab.setAlignment(Qt.AlignCenter)
    lab.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Maximum)
    return lab


def ghost_button(text: str) -> QPushButton:
    b = QPushButton(text)
    b.setObjectName("Ghost")
    b.setCursor(Qt.PointingHandCursor)
    return b


def primary_button(text: str) -> QPushButton:
    b = QPushButton(text)
    b.setObjectName("Primary")
    b.setCursor(Qt.PointingHandCursor)
    return b


class NoWheelOnInputs(QObject):
    """Stop the mouse wheel and touchpad from changing numbers and choices.

    Scrolling over a spin box or a drop-down scrolls the page instead, so a
    value can only change when it is typed or clicked.
    """

    def eventFilter(self, obj, event):
        if event.type() == QEvent.Wheel and isinstance(obj, (QAbstractSpinBox, QComboBox)):
            w = obj.parentWidget()
            while w is not None and not isinstance(w, QAbstractScrollArea):
                w = w.parentWidget()
            if w is not None:
                QApplication.sendEvent(w.verticalScrollBar(), event)
            return True
        return False


def install_no_wheel(app) -> NoWheelOnInputs:
    f = NoWheelOnInputs(app)
    app.installEventFilter(f)
    return f


def collapsible(title: str, content: QWidget, expanded: bool = False) -> QWidget:
    """A text button that shows or hides ``content`` underneath it."""
    box = QWidget()
    v = QVBoxLayout(box)
    v.setContentsMargins(0, 0, 0, 0)
    v.setSpacing(8)
    btn = ghost_button("")
    btn.setStyleSheet("text-align:left; padding-left:0;")

    def _set(on):
        content.setVisible(on)
        btn.setText(("▾  " if on else "▸  ") + title)
    btn.clicked.connect(lambda: _set(not content.isVisible()))
    v.addWidget(btn, 0, Qt.AlignLeft)
    v.addWidget(content)
    _set(expanded)
    box.toggle_button = btn
    return box
