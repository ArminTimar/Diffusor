"""The graph editor: every axis and curve of the plot in one column beside the
graph, applied as typed.

It replaces matplotlib's own figure-options dialog, which asks first for one set
of axes, then shows the axes and the curves on separate tabs, one curve at a time.

What the user changes is kept as *overrides* on the plot, keyed by stable names
(the gid of an axes or a curve, or a curve's legend text), not on the matplotlib
objects themselves. The plot draws itself again after every fit, and the overrides
are laid over each new drawing, so a colour or an axis label chosen once survives
the next fit. Axis limits are dropped when a different profile is loaded.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional

from matplotlib.collections import Collection, LineCollection, PolyCollection
from matplotlib.colors import to_hex, to_rgba
from matplotlib.container import BarContainer, ErrorbarContainer
from matplotlib.lines import Line2D
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor, QDoubleValidator
from PySide6.QtWidgets import (QCheckBox, QColorDialog, QComboBox, QDoubleSpinBox,
                               QFrame, QHBoxLayout, QLabel, QLineEdit,
                               QPushButton, QScrollArea, QVBoxLayout, QWidget)

from . import theme

# names shown for the axes Diffusor draws, by the gid each one is given
PANEL_NAMES = {
    "profile": "Profile", "histogram": "Monte Carlo",
    "mc_time": "Time histogram", "mc_T": "Temperature", "mc_f": "Oxygen fugacity",
    "mc_D": "Diffusion coefficient",
}
# what the user sees for one axis of one panel; anything missing is "<panel>, x" or "y"
AXIS_NAMES = {
    ("profile", "x"): "Distance", ("profile", "y"): "Composition",
    ("histogram", "x"): "Time", ("histogram", "y"): "Draws",
    ("mc_time", "x"): "Time", ("mc_time", "y"): "Draws",
}
LINESTYLES = [("Solid", "-"), ("Dashed", "--"), ("Dotted", ":"), ("Dash-dot", "-."),
              ("No line", "None")]
MARKERS = [("No marker", "None"), ("Circle", "o"), ("Square", "s"), ("Triangle", "^"),
           ("Diamond", "D"), ("Cross", "x"), ("Plus", "+"), ("Dot", ".")]
LEGEND_PLACES = [("Best place", "best"), ("Top right", "upper right"),
                 ("Top left", "upper left"), ("Bottom left", "lower left"),
                 ("Bottom right", "lower right"), ("Right, middle", "center right"),
                 ("Left, middle", "center left")]
BASE_TEXT_SIZE = 10.0


def new_overrides() -> dict:
    return {"figure": {}, "axes": {}, "curves": {}}


# ====================================================================== curves
@dataclass
class Curve:
    """One thing drawn in a panel, possibly made of several matplotlib artists."""
    key: str                         # stable name: a gid, or else the legend text
    axes_gid: str
    members: List[object] = field(default_factory=list)

    @property
    def ident(self) -> str:
        return f"{self.axes_gid}/{self.key}"

    @property
    def first(self):
        return self.members[0]

    def kind(self) -> str:
        """"line" (style, width and marker apply), "lines" (width applies) or "fill"."""
        if any(isinstance(m, (Line2D, ErrorbarContainer)) for m in self.members):
            return "line"
        if any(isinstance(m, LineCollection) for m in self.members):
            return "lines"
        return "fill"

    def legend_text(self) -> str:
        for m in self.members:
            label = m.get_label() if hasattr(m, "get_label") else ""
            if label and not str(label).startswith("_"):
                return str(label)
        return ""

    # ------------------------------------------------------------ current look
    def visible(self) -> bool:
        return _artist_visible(self.first)

    def color(self) -> str:
        m = self.first
        try:
            if isinstance(m, ErrorbarContainer):
                return to_hex(m.lines[0].get_color())
            if isinstance(m, Line2D):
                return to_hex(m.get_color())
            if isinstance(m, BarContainer):
                return to_hex(m.patches[0].get_facecolor()) if m.patches else "#000000"
            if isinstance(m, LineCollection):
                return to_hex(m.get_color()[0])
            if isinstance(m, Collection):
                return to_hex(m.get_facecolor()[0])
        except (IndexError, ValueError, TypeError):
            pass
        return "#000000"

    def alpha(self) -> float:
        m = self.first
        if isinstance(m, ErrorbarContainer):
            m = m.lines[0]
        elif isinstance(m, BarContainer):
            m = m.patches[0] if m.patches else None
        a = m.get_alpha() if m is not None else None
        return 1.0 if a is None else float(a)

    def line(self) -> Optional[Line2D]:
        for m in self.members:
            if isinstance(m, ErrorbarContainer):
                return m.lines[0]
            if isinstance(m, Line2D):
                return m
        return None

    def linewidth(self) -> float:
        ln = self.line()
        if ln is not None:
            return float(ln.get_linewidth())
        if isinstance(self.first, LineCollection):
            return float(self.first.get_linewidth()[0])
        return 1.0


def _artist_visible(m) -> bool:
    if isinstance(m, ErrorbarContainer):
        return m.lines[0].get_visible()
    if isinstance(m, BarContainer):
        return bool(m.patches) and m.patches[0].get_visible()
    return m.get_visible()


def _key(artist) -> Optional[str]:
    gid = artist.get_gid() if hasattr(artist, "get_gid") else None
    if gid:
        return str(gid)
    label = artist.get_label() if hasattr(artist, "get_label") else None
    if label and not str(label).startswith("_"):
        return str(label)
    return None


def curves(ax) -> List[Curve]:
    """Everything in ``ax`` the editor can restyle, in drawing order, same names merged."""
    gid = ax.get_gid()
    found: Dict[str, Curve] = {}
    claimed = set()
    for c in ax.containers:
        if isinstance(c, ErrorbarContainer):
            data, caps, bars = c.lines
            claimed.update(id(a) for a in (data, *caps, *bars))
            k = _key(data) or _key(c)
        elif isinstance(c, BarContainer):
            claimed.update(id(p) for p in c.patches)
            k = _key(c) or "Histogram bars"
        else:
            continue
        if k:
            found.setdefault(k, Curve(k, gid)).members.append(c)
    for artist in list(ax.lines) + list(ax.collections):
        if id(artist) in claimed:
            continue
        k = _key(artist)
        if k:
            found.setdefault(k, Curve(k, gid)).members.append(artist)
    return list(found.values())


def editable_axes(fig) -> list:
    """The panels with a gid and real axes (the 'not varied' notes have no ticks)."""
    return [ax for ax in fig.axes if ax.get_gid() in PANEL_NAMES and ax.axison
            and (len(ax.get_xticks()) or len(ax.get_yticks()))]


def _shares_x(ax, other) -> bool:
    try:
        return ax.get_shared_x_axes().joined(ax, other)
    except AttributeError:
        return False


def axis_rows(fig) -> List[tuple]:
    """(axes, "x" or "y", name) for every axis to show, a shared x axis once.

    A shared x axis is edited on the panel that carries its label, which is the
    bottom one, so the label typed is the one seen. Each group of panels sharing
    an x axis gives that axis first, then the vertical axis of each panel."""
    rows = []
    axes = editable_axes(fig)
    done = set()
    for ax in axes:
        if id(ax) in done:
            continue
        group = [a for a in axes if a is ax or _shares_x(ax, a)]
        done.update(id(a) for a in group)
        owner = next((a for a in reversed(group) if a.xaxis.get_label().get_text()),
                     group[-1])
        rows.append((owner, "x", _axis_name(owner, "x")))
        rows.extend((a, "y", _axis_name(a, "y")) for a in group)
    return rows


def _axis_name(ax, which: str) -> str:
    gid = ax.get_gid()
    return AXIS_NAMES.get((gid, which)) or f"{PANEL_NAMES.get(gid, gid)}, {which} axis"


# ====================================================================== applying
def apply_overrides(fig, ov: dict) -> None:
    """Lay the user's changes over the figure as it is drawn now.

    Every value is set outright, so this can run again on the same drawing after
    each change in the editor."""
    axes = editable_axes(fig)
    for ax in axes:
        # remembered at the first pass over a new drawing, before any legend is removed
        if not hasattr(ax, "_diffusor_legend"):
            ax._diffusor_legend = ax.get_legend() is not None
    if not ov or not any(ov.values()):
        return
    figure = ov.get("figure", {})
    size = figure.get("text_size")
    if axes and "title" in figure:
        t = _title_of(axes[0])
        t.set_text(figure["title"])
        t.set_color(theme.TEXT)
        t.set_fontsize((size or BASE_TEXT_SIZE) + 1)
    for ax in axes:
        gid = ax.get_gid()
        a = ov.get("axes", {}).get(gid, {})
        for which in ("x", "y"):
            axis = ax.xaxis if which == "x" else ax.yaxis
            if f"{which}label" in a:
                axis.set_label_text(a[f"{which}label"])
            if f"{which}scale" in a:
                try:
                    (ax.set_xscale if which == "x" else ax.set_yscale)(a[f"{which}scale"])
                except ValueError:
                    pass
            lim = a.get(f"{which}lim")
            if lim and any(v is not None for v in lim):
                (ax.set_xlim if which == "x" else ax.set_ylim)(*lim)
        if size:
            ax.tick_params(labelsize=size - 1)
            ax.xaxis.label.set_fontsize(size)
            ax.yaxis.label.set_fontsize(size)
            for loc in ("center", "left", "right"):
                t = ax._left_title if loc == "left" else (ax._right_title if loc == "right"
                                                          else ax.title)
                if t.get_text():
                    t.set_fontsize(size + 1)
        if "grid" in figure:
            ax.grid(bool(figure["grid"]), color=theme.PLOT_GRID, lw=0.8, alpha=0.9)
        for c in curves(ax):
            style = ov.get("curves", {}).get(c.ident)
            if style:
                _style_curve(c, style)
    for ax in axes:
        _legend_again(ax, figure)


def _title_of(ax):
    """The title in use: Diffusor puts some on the left."""
    if ax._left_title.get_text() and not ax.title.get_text():
        return ax._left_title
    return ax.title


def _children(m) -> list:
    if isinstance(m, ErrorbarContainer):
        data, caps, bars = m.lines
        return [data, *caps, *bars]
    if isinstance(m, BarContainer):
        return list(m.patches)
    return [m]


def _style_curve(c: Curve, s: dict) -> None:
    for m in c.members:
        if "visible" in s:
            for ch in _children(m):
                ch.set_visible(bool(s["visible"]))
        if "label" in s:
            if isinstance(m, ErrorbarContainer) or isinstance(m, BarContainer):
                m.set_label(s["label"])
            elif m.get_label() and not str(m.get_label()).startswith("_") or m is c.first:
                m.set_label(s["label"] or "_nolegend_")
        color, alpha = s.get("color"), s.get("alpha")
        if isinstance(m, ErrorbarContainer):
            data, caps, bars = m.lines
            _style_line(data, s)
            if color:
                for ch in (*caps, *bars):
                    ch.set_color(to_rgba(color, 0.45))
        elif isinstance(m, Line2D):
            _style_line(m, s)
        elif isinstance(m, BarContainer):
            for p in m.patches:
                if color:
                    p.set_facecolor(color)
                if alpha is not None:
                    p.set_alpha(alpha)
        elif isinstance(m, LineCollection):
            if color:
                m.set_color(color)
            if alpha is not None:
                m.set_alpha(alpha)
            if "linewidth" in s:
                m.set_linewidth(s["linewidth"])
        elif isinstance(m, (PolyCollection, Collection)):
            if color:
                m.set_facecolor(color)
                m.set_edgecolor("none")
            if alpha is not None:
                m.set_alpha(alpha)


def _style_line(ln: Line2D, s: dict) -> None:
    if s.get("color"):
        ln.set_color(s["color"])
        ln.set_markerfacecolor(s["color"])
        ln.set_markeredgecolor(s["color"])
    if s.get("alpha") is not None:
        ln.set_alpha(s["alpha"])
    if "linestyle" in s:
        ln.set_linestyle(s["linestyle"])
    if "linewidth" in s:
        ln.set_linewidth(s["linewidth"])
    if "marker" in s:
        ln.set_marker(s["marker"])
    if "markersize" in s:
        ln.set_markersize(s["markersize"])


# legend order: the data, then the model, then what surrounds it
LEGEND_ORDER = ["Measured points", "Cut points", "Best fit", "Initial profile", "68% band",
                "95% band"]


def legend_entries(ax) -> list:
    """(handle, text) pairs in reading order: matplotlib lists error-bar data last."""
    handles, labels = ax.get_legend_handles_labels()

    def rank(h):
        k = _key(h.lines[0] if isinstance(h, ErrorbarContainer) else h)
        return LEGEND_ORDER.index(k) if k in LEGEND_ORDER else len(LEGEND_ORDER)
    pairs = list(zip(handles, labels))
    return sorted(pairs, key=lambda p: rank(p[0]))


def _legend_again(ax, figure: dict) -> None:
    """Rebuild a legend, so it shows the new looks and leaves out hidden curves.

    Only panels that were drawn with a legend get one back."""
    if not getattr(ax, "_diffusor_legend", False):
        return
    show = figure.get("legend", True)
    old = ax.get_legend()
    if old is not None:
        old.remove()
    if not show:
        return
    keep = [(h, lab) for h, lab in legend_entries(ax) if _artist_visible(h)]
    if not keep:
        return
    size = (figure.get("text_size") or BASE_TEXT_SIZE) - 1
    leg = ax.legend([h for h, _ in keep], [lab for _, lab in keep], fontsize=size,
                    loc=figure.get("legend_loc", "best"), frameon=True, fancybox=False,
                    framealpha=0.95, facecolor=theme.SURFACE,
                    edgecolor=theme.BORDER_STRONG, borderpad=0.6)
    leg.get_frame().set_linewidth(0.9)


# ====================================================================== dialog
def _number(text: str) -> Optional[float]:
    t = text.strip().replace(",", ".")
    if not t:
        return None
    try:
        v = float(t)
    except ValueError:
        return None
    return v if math.isfinite(v) else None


def _fmt(v: float) -> str:
    return f"{v:.6g}"


class ColorButton(QPushButton):
    """A swatch that opens a colour picker."""
    picked = Signal(str)

    def __init__(self, color: str, parent=None):
        super().__init__(parent)
        self.setFixedSize(34, 26)
        self.setCursor(Qt.PointingHandCursor)
        self.setToolTip("Choose a colour")
        self.clicked.connect(self._choose)
        self.set_color(color)

    def set_color(self, color: str):
        self._color = color
        self.setStyleSheet(f"QPushButton {{ background:{color}; border:1px solid "
                           f"{theme.BORDER_STRONG}; border-radius:5px; padding:0; }}"
                           f"QPushButton:hover {{ border:2px solid {theme.ACCENT}; }}")

    def _choose(self):
        c = QColorDialog.getColor(QColor(self._color), self, "Curve colour")
        if c.isValid():
            self.set_color(c.name())
            self.picked.emit(c.name())


def _small_spin(lo, hi, value, step, decimals=1) -> QDoubleSpinBox:
    s = QDoubleSpinBox()
    s.setRange(lo, hi)
    s.setDecimals(decimals)
    s.setSingleStep(step)
    s.setValue(value)
    s.setKeyboardTracking(False)
    s.setFixedWidth(84)
    return s


def _combo(items, current) -> QComboBox:
    cmb = QComboBox()
    for text, data in items:
        cmb.addItem(text, data)
    i = cmb.findData(current)
    cmb.setCurrentIndex(i if i >= 0 else 0)
    # never narrower than the longest entry: the grid squeezed "No marker" to "No marke"
    cmb.setSizeAdjustPolicy(QComboBox.AdjustToContents)
    cmb.setMinimumWidth(cmb.sizeHint().width())
    return cmb


def _head(text: str) -> QLabel:
    lab = QLabel(text)
    lab.setObjectName("FieldLabel")
    return lab


def _section(title: str):
    frame = QFrame()
    frame.setObjectName("Card")
    v = QVBoxLayout(frame)
    v.setContentsMargins(14, 10, 14, 12)
    v.setSpacing(6)
    lab = QLabel(title)
    lab.setObjectName("H2")
    v.addWidget(lab)
    return frame, v


# denser inputs than the steps use, so the whole graph fits on one page
COMPACT = (
    "QLineEdit, QDoubleSpinBox, QSpinBox, QComboBox { padding: 2px 7px; min-height: 17px; }")


class GraphEditor(QWidget):
    """Axes, curves, legend and text of the plot, in a column beside the graph.

    The column sits next to the plot, so every change shows on the graph as it is
    made. Changes apply at once; ``done`` is emitted when the user closes it."""
    done = Signal()
    WIDTH = 400

    def __init__(self, plot, parent=None):
        super().__init__(parent)
        self.plot = plot
        self.setObjectName("Page")
        self.setStyleSheet(COMPACT)
        self.setFixedWidth(self.WIDTH)
        root = QVBoxLayout(self)
        root.setContentsMargins(12, 10, 8, 10)
        root.setSpacing(8)
        top = QHBoxLayout()
        title = QLabel("Edit graph")
        title.setObjectName("H2")
        top.addWidget(title)
        top.addStretch(1)
        reset = QPushButton("Reset all")
        reset.setToolTip("Undo every change made here and draw the graph as Diffusor made it")
        reset.clicked.connect(self._reset)
        top.addWidget(reset)
        close = QPushButton("Done")
        close.setObjectName("Primary")
        close.setToolTip("Close the editor. The changes stay on the graph.")
        close.clicked.connect(self.done.emit)
        top.addWidget(close)
        root.addLayout(top)
        self.body = QWidget()
        self.body.setObjectName("Page")
        self.body_layout = QVBoxLayout(self.body)
        self.body_layout.setContentsMargins(0, 0, 6, 0)
        self.body_layout.setSpacing(10)
        area = QScrollArea()
        area.setWidget(self.body)
        area.setWidgetResizable(True)
        area.setFrameShape(QFrame.NoFrame)
        area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        root.addWidget(area, 1)
        hint = QLabel("Changes show on the graph at once and are kept when it is drawn "
                      "again after a fit. Axis ranges reset when another profile is loaded.")
        hint.setObjectName("Hint")
        hint.setWordWrap(True)
        root.addWidget(hint)
        self.rebuild()
        # a new fit draws new curves; the open editor follows
        plot.redrawn.connect(self._plot_redrawn)

    def _plot_redrawn(self):
        if self.isVisible():
            self.rebuild()

    # ---------------------------------------------------------------- building
    def rebuild(self):
        """Read the plot as it is now and lay out a block for each axis and curve."""
        while self.body_layout.count():
            w = self.body_layout.takeAt(0).widget()
            if w is not None:
                w.setParent(None)
                w.deleteLater()
        fig = self.plot.figure
        axes = editable_axes(fig)
        if not axes:
            lab = QLabel("There is nothing to edit in this view yet.")
            lab.setObjectName("Hint")
            self.body_layout.addWidget(lab)
            self.body_layout.addStretch(1)
            return
        self.body_layout.addWidget(self._text_section(axes))
        self.body_layout.addWidget(self._axes_section(fig))
        self.body_layout.addWidget(self._curves_section(axes))
        self.body_layout.addStretch(1)

    def _text_section(self, axes) -> QWidget:
        frame, v = _section("Graph")
        fig_ov = self.plot.overrides["figure"]
        title = QLineEdit(_title_of(axes[0]).get_text())
        title.setPlaceholderText("No title")
        title.editingFinished.connect(lambda: self._set_figure("title", title.text()))
        v.addWidget(_head("Title"))
        v.addWidget(title)
        legend = QCheckBox("Legend")
        legend.setToolTip("Show or hide the legend")
        legend.setChecked(fig_ov.get("legend", any(getattr(a, "_diffusor_legend", False)
                                                    for a in axes)))
        place = _combo(LEGEND_PLACES, fig_ov.get("legend_loc", "best"))
        place.setToolTip("Where the legend goes")
        place.setEnabled(legend.isChecked())
        legend.toggled.connect(lambda on: (place.setEnabled(on),
                                           self._set_figure("legend", bool(on))))
        place.currentIndexChanged.connect(
            lambda _i: self._set_figure("legend_loc", place.currentData()))
        row = QHBoxLayout()
        row.addWidget(legend)
        row.addWidget(place, 1)
        v.addLayout(row)
        size = _small_spin(6, 24, fig_ov.get("text_size", BASE_TEXT_SIZE), 0.5)
        size.setPrefix("text ")
        size.setSuffix(" pt")
        size.setFixedWidth(110)
        size.setToolTip("Axis labels at this size, tick numbers and legend one point smaller")
        size.valueChanged.connect(lambda val: self._set_figure("text_size", float(val)))
        grid_on = QCheckBox("Grid lines")
        grid_on.setChecked(fig_ov.get("grid", True))
        grid_on.toggled.connect(lambda on: self._set_figure("grid", bool(on)))
        row = QHBoxLayout()
        row.addWidget(size)
        row.addSpacing(12)
        row.addWidget(grid_on)
        row.addStretch(1)
        v.addLayout(row)
        return frame

    def _axes_section(self, fig) -> QWidget:
        frame, v = _section("Axes")
        for n, (ax, which, name) in enumerate(axis_rows(fig)):
            if n:
                v.addSpacing(4)
            gid = ax.get_gid()
            a = self.plot.overrides["axes"].get(gid, {})
            axis = ax.xaxis if which == "x" else ax.yaxis
            lo, hi = ax.get_xlim() if which == "x" else ax.get_ylim()
            fixed = bool(a.get(f"{which}lim"))
            lab = QLabel(name)
            lab.setStyleSheet("font-weight:600;")
            scale = _combo([("Linear", "linear"), ("Logarithmic", "log")],
                           ax.get_xscale() if which == "x" else ax.get_yscale())
            scale.setToolTip("Axis scale")
            head = QHBoxLayout()
            head.addWidget(lab, 1)
            head.addWidget(scale)
            v.addLayout(head)
            text = QLineEdit(axis.get_label().get_text())
            text.setPlaceholderText("No label")
            text.setToolTip("Axis label")
            text.editingFinished.connect(
                lambda g=gid, w=which, t=text: self._set_axes(g, f"{w}label", t.text()))
            v.addWidget(text)
            e_lo, e_hi = QLineEdit(_fmt(lo)), QLineEdit(_fmt(hi))
            for e, tip in ((e_lo, "From"), (e_hi, "To")):
                e.setValidator(QDoubleValidator())
                e.setMinimumWidth(70)
                e.setToolTip(tip)
                e.setEnabled(fixed)
            auto = QCheckBox("Automatic")
            auto.setToolTip("Let the range follow the data. Untick to type your own.")
            auto.setChecked(not fixed)

            def set_range(g=gid, w=which, a_lo=e_lo, a_hi=e_hi, chk=auto):
                if chk.isChecked():
                    return
                self._set_axes(g, f"{w}lim", (_number(a_lo.text()), _number(a_hi.text())))

            def toggle_auto(on, g=gid, w=which, a_lo=e_lo, a_hi=e_hi):
                a_lo.setEnabled(not on)
                a_hi.setEnabled(not on)
                if on:
                    self._set_axes(g, f"{w}lim", None)
                else:
                    self._set_axes(g, f"{w}lim", (_number(a_lo.text()), _number(a_hi.text())))

            e_lo.editingFinished.connect(set_range)
            e_hi.editingFinished.connect(set_range)
            auto.toggled.connect(toggle_auto)
            scale.currentIndexChanged.connect(
                lambda _i, g=gid, w=which, s=scale: self._set_axes(g, f"{w}scale", s.currentData()))
            rng = QHBoxLayout()
            rng.addWidget(e_lo, 1)
            to = QLabel("to")
            to.setObjectName("Hint")
            rng.addWidget(to)
            rng.addWidget(e_hi, 1)
            rng.addWidget(auto)
            v.addLayout(rng)
        return frame

    def _curves_section(self, axes) -> QWidget:
        frame, v = _section("Curves")
        several = sum(1 for ax in axes if curves(ax)) > 1
        count = 0
        for ax in axes:
            found = curves(ax)
            if not found:
                continue
            if several:
                sub = QLabel(PANEL_NAMES.get(ax.get_gid(), ax.get_gid()).upper())
                sub.setObjectName("SummaryHead")
                v.addWidget(sub)
            for c in found:
                if count:
                    line = QFrame()
                    line.setObjectName("Divider")
                    line.setFixedHeight(1)
                    v.addWidget(line)
                self._curve_block(v, c)
                count += 1
        if not count:
            lab = QLabel("No curves in this view.")
            lab.setObjectName("Hint")
            v.addWidget(lab)
        return frame

    def _curve_block(self, v: QVBoxLayout, c: Curve):
        """One curve: shown or not, name and colour; legend text; line; marker."""
        kind = c.kind()
        ident = c.ident
        show = QCheckBox()
        show.setChecked(c.visible())
        show.setToolTip("Show or hide this curve")
        show.toggled.connect(lambda on: self._set_curve(ident, "visible", bool(on)))
        name = QLabel(c.key)
        name.setToolTip(c.key)
        name.setStyleSheet("font-weight:600;")
        name.setMinimumWidth(10)
        color = ColorButton(c.color())
        color.picked.connect(lambda col: self._set_curve(ident, "color", col))
        opacity = _small_spin(0.0, 1.0, c.alpha(), 0.05, 2)
        opacity.setPrefix("opacity ")
        opacity.setFixedWidth(124)
        opacity.setToolTip("0 is invisible, 1 is solid")
        opacity.valueChanged.connect(lambda val: self._set_curve(ident, "alpha", float(val)))
        head = QHBoxLayout()
        head.addWidget(show)
        head.addWidget(name, 1)
        head.addWidget(opacity)
        head.addWidget(color)
        v.addLayout(head)

        legend = QLineEdit(c.legend_text())
        legend.setPlaceholderText("not in the legend")
        legend.setToolTip("Legend text")
        legend.editingFinished.connect(lambda: self._set_curve(ident, "label", legend.text()))
        v.addWidget(legend)

        ln = c.line()
        style = _combo(LINESTYLES, ln.get_linestyle() if ln is not None else "-")
        style.setToolTip("Line style")
        style.currentIndexChanged.connect(
            lambda _i: self._set_curve(ident, "linestyle", style.currentData()))
        width = _small_spin(0.0, 12.0, c.linewidth(), 0.5)
        width.setPrefix("width ")
        width.setFixedWidth(100)
        width.setToolTip("Line width (pt)")
        width.valueChanged.connect(lambda val: self._set_curve(ident, "linewidth", float(val)))
        marker = _combo(MARKERS, str(ln.get_marker()) if ln is not None else "None")
        marker.setToolTip("Marker")
        marker.currentIndexChanged.connect(
            lambda _i: self._set_curve(ident, "marker", marker.currentData()))
        size = _small_spin(0.0, 30.0, float(ln.get_markersize()) if ln is not None else 0.0, 0.5)
        size.setPrefix("size ")
        size.setFixedWidth(100)
        size.setToolTip("Marker size (pt)")
        size.valueChanged.connect(lambda val: self._set_curve(ident, "markersize", float(val)))
        style.setEnabled(kind == "line")
        marker.setEnabled(kind == "line")
        size.setEnabled(kind == "line")
        width.setEnabled(kind in ("line", "lines"))
        if kind == "fill":
            for w in (style, width, marker, size):
                w.setToolTip("A filled area has a colour and an opacity only")
        row = QHBoxLayout()
        row.addWidget(style, 1)
        row.addWidget(width)
        v.addLayout(row)
        row = QHBoxLayout()
        row.addWidget(marker, 1)
        row.addWidget(size)
        v.addLayout(row)

    # ---------------------------------------------------------------- changes
    def _changed(self):
        self.plot.restyle()

    def _set_figure(self, key, value):
        self.plot.overrides["figure"][key] = value
        self._changed()

    def _set_axes(self, gid, key, value):
        a = self.plot.overrides["axes"].setdefault(gid, {})
        if value is None:
            a.pop(key, None)
            if key.endswith("lim"):
                # back to automatic: the range follows the data again
                for ax in self.plot.figure.axes:
                    if ax.get_gid() == gid:
                        ax.autoscale(enable=True, axis=key[0])
                        ax.autoscale_view()
        else:
            a[key] = value
        self._changed()

    def _set_curve(self, ident, key, value):
        self.plot.overrides["curves"].setdefault(ident, {})[key] = value
        self._changed()

    def _reset(self):
        self.plot.reset_overrides()
        self.rebuild()
