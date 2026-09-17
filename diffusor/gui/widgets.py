"""Small layout helpers shared by the step pages."""
from __future__ import annotations

from typing import Optional

from PySide6.QtCore import QEvent, QObject, Qt
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
        self.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.MinimumExpanding)
        self.setAlignment(Qt.AlignTop | Qt.AlignLeft)

    def _sync(self):
        w = self.width()
        if w > 0:
            self.setMinimumHeight(self.heightForWidth(w))

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._sync()

    def setText(self, text):
        super().setText(text)
        self._sync()
        self.updateGeometry()


def card(title: str = "", subtitle: str = "") -> tuple[QFrame, QVBoxLayout]:
    """A white rounded panel with an optional heading; returns (frame, body layout)."""
    frame = QFrame()
    frame.setObjectName("Card")
    outer = QVBoxLayout(frame)
    outer.setContentsMargins(18, 16, 18, 16)
    outer.setSpacing(10)
    if title:
        lab = QLabel(title)
        lab.setObjectName("H2")
        outer.addWidget(lab)
    if subtitle:
        sub = WrapLabel(subtitle)
        sub.setObjectName("Sub")
        outer.addWidget(sub)
    body = QVBoxLayout()
    body.setSpacing(10)
    outer.addLayout(body)
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
