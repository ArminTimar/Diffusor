"""Draw a boundary on an image and read profiles perpendicular to it.

A nonmodal window, opened from File > Extract profile from image or the Data
step. The image is on the left with the profile lines, rejected values and
exclusion areas drawn over it, and the averaged profile underneath; the
settings are on the right. Every change re-runs the extraction, which takes a
fraction of a second, so what is saved is always what is shown.

Mouse, with the matching tool selected (the plot toolbar's pan and zoom tools
take the mouse while they are on):

* Boundary: click to add guideline points, right-click to remove the last. The
  lines start on the right of the guideline, looking from its first point.
* Scale bar: click both ends of the scale bar and type its length.
* Legend: click the low-value end of the colour legend, then the high end.
* Exclude area: click around an inclusion or lamella, right-click to close.
"""
from __future__ import annotations

import io
import traceback
from pathlib import Path
from typing import List, Optional

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from matplotlib.patches import Polygon
from PySide6.QtCore import Qt, QTimer
from PySide6.QtWidgets import (QButtonGroup, QCheckBox, QComboBox, QDialog, QDialogButtonBox,
                               QDoubleSpinBox, QFileDialog, QHBoxLayout, QInputDialog,
                               QLineEdit, QMessageBox, QPushButton, QSpinBox, QSplitter,
                               QVBoxLayout, QWidget)

from ..dataio import images as imgio
from ..dataio.image_profiles import (CLIP_METHODS, KEPT, OUTSIDE, CleaningSettings,
                                     ExtractionSettings, extract_profiles, is_extraction_workbook,
                                     suggest_limits, table_from_extraction, write_workbook)
from . import theme
from .widgets import card, field, ghost_button, note, pair, primary_button, row, scrollable

TOOLS = {
    "boundary": "Click along the zone boundary to add guideline points. Right-click removes "
                "the last point. The lines start on the right of the guideline, looking from "
                "its first point.",
    "scale": "Click both ends of the scale bar, then type its length.",
    "legend": "Click the low-value end of the colour legend, then the high-value end. Stay "
              "inside the coloured bar.",
    "exclude": "Click around an inclusion, crack or lamella. Right-click closes the area.",
}
COLORMAPS = ["jet", "turbo", "rainbow", "nipy_spectral", "viridis", "plasma", "inferno",
             "magma", "cividis", "hot", "gray"]
MAX_DRAWN_LINES = 40
MAX_REJECTED_DOTS = 20000


def _spin(lo, hi, val, dec=2, step=1.0, suffix=""):
    s = QDoubleSpinBox()
    s.setRange(lo, hi)
    s.setDecimals(dec)
    s.setSingleStep(step)
    s.setValue(val)
    if suffix:
        s.setSuffix(suffix)
    # a huge range would otherwise size the box for its longest possible number
    s.setMinimumWidth(90)
    return s


class RawDialog(QDialog):
    """Width, height and type of a headerless binary image."""

    def __init__(self, path: Path, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Read as raw binary")
        v = QVBoxLayout(self)
        size = path.stat().st_size
        v.addWidget(note(f"{path.name} is {size:,} bytes and has no header Diffusor knows. "
                         "Give its layout to read it as raw numbers.", "Sub"))
        self.sp_w, self.sp_h = QSpinBox(), QSpinBox()
        for s in (self.sp_w, self.sp_h):
            s.setRange(1, 1_000_000)
        self.sp_w.setValue(512)
        self.sp_h.setValue(512)
        v.addWidget(pair(field("Width (px)", self.sp_w), field("Height (px)", self.sp_h)))
        self.cmb_type = QComboBox()
        self.cmb_type.addItems(["uint8", "uint16", "int16", "uint32", "int32", "float32",
                                "float64"])
        self.cmb_type.setCurrentText("uint16")
        self.cmb_order = QComboBox()
        self.cmb_order.addItems(["little", "big"])
        v.addWidget(pair(field("Data type", self.cmb_type), field("Byte order", self.cmb_order)))
        self.sp_off, self.sp_bands = QSpinBox(), QSpinBox()
        self.sp_off.setRange(0, 1 << 30)
        self.sp_bands.setRange(1, 64)
        v.addWidget(pair(field("Header bytes to skip", self.sp_off),
                         field("Channels", self.sp_bands)))
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        v.addWidget(bb)

    def options(self) -> dict:
        return dict(width=self.sp_w.value(), height=self.sp_h.value(),
                    dtype=self.cmb_type.currentText(), byte_order=self.cmb_order.currentText(),
                    offset=self.sp_off.value(), bands=self.sp_bands.value())


class ImageExtractorDialog(QDialog):
    def __init__(self, owner=None):
        super().__init__(owner)
        self.owner = owner
        self.setWindowTitle("Extract a profile from an image")
        self.setWindowFlag(Qt.WindowMaximizeButtonHint, True)
        self.resize(1360, 880)
        self.image: Optional[imgio.LoadedImage] = None
        self.values = None                     # (values, no_value, label)
        self.color_scale: Optional[imgio.ColorScale] = None
        self.guideline: List[tuple] = []
        self.exclusions: List[List[tuple]] = []
        self.current_poly: List[tuple] = []
        self.clicks: List[tuple] = []          # scale bar or legend ends being drawn
        self.legend_pts: Optional[List[tuple]] = None
        self.extraction = None
        self._overlay = []
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(120)
        self._timer.timeout.connect(self.recompute)

        outer = QHBoxLayout(self)
        outer.setContentsMargins(12, 12, 12, 12)
        outer.setSpacing(12)
        split = QSplitter(Qt.Vertical)
        split.addWidget(self._image_panel())
        split.addWidget(self._preview_panel())
        split.setSizes([620, 240])
        outer.addWidget(split, 1)
        side = scrollable(self._controls())
        side.setFixedWidth(400)
        outer.addWidget(side)
        self._set_tool("boundary")
        self._mode_widgets()

    # ================================================================ layout
    def _image_panel(self) -> QWidget:
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setContentsMargins(0, 0, 0, 0)
        self.fig = Figure(figsize=(8, 6), layout="constrained", facecolor=theme.SURFACE)
        self.canvas = FigureCanvasQTAgg(self.fig)
        self.canvas.setFocusPolicy(Qt.StrongFocus)
        self.toolbar = NavigationToolbar2QT(self.canvas, w)
        lay.addWidget(self.toolbar)
        lay.addWidget(self.canvas, 1)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_axis_off()
        self.ax.text(0.5, 0.5, "Open an image to begin", ha="center", va="center",
                     color=theme.TEXT_MUTED, transform=self.ax.transAxes)
        self.canvas.mpl_connect("button_press_event", self._on_click)
        self.canvas.mpl_connect("key_press_event", self._on_key)
        return w

    def _preview_panel(self) -> QWidget:
        self.pfig = Figure(figsize=(8, 2.4), layout="constrained", facecolor=theme.SURFACE)
        self.pcanvas = FigureCanvasQTAgg(self.pfig)
        self.pax = self.pfig.add_subplot(111)
        theme.apply_plot_style(self.pfig, [self.pax])
        return self.pcanvas

    def _controls(self) -> QWidget:
        box = QWidget()
        box.setObjectName("Page")
        v = QVBoxLayout(box)
        v.setContentsMargins(0, 0, 6, 0)
        v.setSpacing(10)

        c, b = card("Image")
        btn = primary_button("Open image...")
        btn.clicked.connect(self.open_image)
        b.addWidget(btn)
        self.lbl_file = note("TIFF, PNG, JPEG, BMP, ENVI, raw binary or a text grid of counts.",
                             "Hint")
        b.addWidget(self.lbl_file)
        self.cmb_mode = QComboBox()
        self.cmb_channel = QComboBox()
        self.cmb_mode.currentIndexChanged.connect(self._mode_changed)
        self.cmb_channel.currentIndexChanged.connect(self._values_changed)
        b.addWidget(field("Values from", self.cmb_mode))
        self.fld_channel = field("Channel", self.cmb_channel)
        b.addWidget(self.fld_channel)
        v.addWidget(c)

        c, b = card("Draw")
        self.tool_group = QButtonGroup(self)
        self.tool_buttons = {}
        btns = []
        for key, text in (("boundary", "Boundary"), ("scale", "Scale bar"), ("legend", "Legend"),
                          ("exclude", "Exclude area")):
            t = QPushButton(text)
            t.setCheckable(True)
            t.clicked.connect(lambda _=False, k=key: self._set_tool(k))
            self.tool_group.addButton(t)
            self.tool_buttons[key] = t
            btns.append(t)
        b.addWidget(row(*btns[:2]))
        b.addWidget(row(*btns[2:]))
        self.lbl_tool = note("", "Hint")
        b.addWidget(self.lbl_tool)
        undo = ghost_button("Undo point")
        undo.clicked.connect(self.undo_point)
        clr = ghost_button("Clear boundary")
        clr.clicked.connect(self.clear_boundary)
        clx = ghost_button("Clear areas")
        clx.clicked.connect(self.clear_exclusions)
        b.addWidget(row(undo, clr, clx, spacing=4))
        v.addWidget(c)

        c, b = card("Scale")
        self.sp_px = _spin(0.0, 1e4, 0.0, 6, 0.01, " um/px")
        self.sp_px.setSpecialValueText("not set (pixels)")
        self.sp_px.valueChanged.connect(self._schedule)
        self.lbl_px = note("", "Hint")
        b.addWidget(field("Pixel size", self.sp_px))
        b.addWidget(self.lbl_px)
        v.addWidget(c)

        self.card_legend, b = card("Colour legend")
        self.cmb_legend = QComboBox()
        self.cmb_legend.addItem("Drawn on the image", "drawn")
        for name in COLORMAPS:
            self.cmb_legend.addItem(f"Colour map '{name}'", name)
        self.cmb_legend.currentIndexChanged.connect(self._legend_changed)
        b.addWidget(field("Legend", self.cmb_legend,
                          "Use a named map only when you know the image was drawn with it."))
        self.sp_v0 = _spin(-1e12, 1e12, 0.0, 4)
        self.sp_v1 = _spin(-1e12, 1e12, 100.0, 4)
        b.addWidget(pair(field("Value at the low end", self.sp_v0),
                         field("Value at the high end", self.sp_v1)))
        self.chk_log = QCheckBox("Logarithmic legend")
        self.le_units = QLineEdit()
        self.le_units.setPlaceholderText("units, e.g. wt% MgO or counts")
        b.addWidget(self.chk_log)
        b.addWidget(self.le_units)
        self.sp_de = _spin(1.0, 200.0, 20.0, 1, 1.0)
        b.addWidget(field("Largest colour difference (Delta E)", self.sp_de,
                          "Colours farther than this from every legend colour get no value: "
                          "cracks, epoxy, labels."))
        self.lbl_legend = note("", "Hint")
        b.addWidget(self.lbl_legend)
        for w in (self.sp_v0, self.sp_v1, self.sp_de):
            w.valueChanged.connect(self._legend_changed)
        self.chk_log.toggled.connect(self._legend_changed)
        self.le_units.editingFinished.connect(self._legend_changed)
        v.addWidget(self.card_legend)

        c, b = card("Profile lines")
        self.sp_before = _spin(0.0, 1e5, 50.0, 1, 5.0, " px")
        self.sp_after = _spin(0.0, 1e5, 50.0, 1, 5.0, " px")
        b.addWidget(pair(field("Length before", self.sp_before),
                         field("Length after", self.sp_after)))
        self.sp_spacing = _spin(0.1, 1e4, 1.0, 2, 0.5, " px")
        self.sp_step = _spin(0.1, 1e4, 1.0, 2, 0.5, " px")
        b.addWidget(pair(field("Line spacing", self.sp_spacing),
                         field("Sample step", self.sp_step)))
        self.cmb_interp = QComboBox()
        self.cmb_interp.addItem("Bilinear", "bilinear")
        self.cmb_interp.addItem("Nearest pixel", "nearest")
        self.chk_flip = QCheckBox("Start on the left")
        b.addWidget(row(self.cmb_interp, self.chk_flip))
        self.lbl_lines = note("NIDIS uses 50 px either side and one line per pixel.", "Hint")
        b.addWidget(self.lbl_lines)
        for w in (self.sp_before, self.sp_after, self.sp_spacing, self.sp_step):
            w.valueChanged.connect(self._schedule)
        self.cmb_interp.currentIndexChanged.connect(self._schedule)
        self.chk_flip.toggled.connect(self._schedule)
        v.addWidget(c)

        c, b = card("Cleaning")
        self.chk_limits = QCheckBox("Reject values outside limits")
        self.sp_low = _spin(-1e12, 1e12, 0.0, 3)
        self.sp_high = _spin(-1e12, 1e12, 255.0, 3)
        sug = ghost_button("Suggest")
        sug.setToolTip("0.5th and 99.5th percentiles of the image")
        sug.clicked.connect(self.suggest_limits)
        b.addWidget(row(self.chk_limits, sug))
        b.addWidget(pair(field("Below", self.sp_low), field("Above", self.sp_high)))
        b.addWidget(note("Dark cracks and holes fall below, bright inclusions above.", "Hint"))
        self.sp_grow = QSpinBox()
        self.sp_grow.setRange(0, 50)
        self.sp_grow.setValue(1)
        self.sp_grow.setSuffix(" px")
        b.addWidget(field("Grow rejected areas by", self.sp_grow,
                          "The edge of a crack is part crack, part crystal."))
        self.cmb_clip = QComboBox()
        for key, label in CLIP_METHODS.items():
            self.cmb_clip.addItem(label, key)
        self.cmb_clip.setCurrentIndex(1)
        self.sp_k = _spin(0.5, 10.0, 3.0, 1, 0.5)
        self.cmb_clip.currentIndexChanged.connect(self._clip_changed)
        b.addWidget(field("Outliers across the lines", row(self.cmb_clip, self.sp_k, spacing=6)))
        self.chk_drop = QCheckBox("Drop whole lines with more rejected than")
        self.chk_drop.setChecked(True)
        self.sp_drop = _spin(1.0, 100.0, 50.0, 0, 5.0, " %")
        b.addWidget(self.chk_drop)
        b.addWidget(row(self.sp_drop))
        self.chk_show = QCheckBox("Show rejected values on the image")
        self.chk_show.setChecked(True)
        b.addWidget(self.chk_show)
        for w in (self.chk_limits, self.chk_drop, self.chk_show):
            w.toggled.connect(self._schedule)
        for w in (self.sp_low, self.sp_high, self.sp_k, self.sp_drop):
            w.valueChanged.connect(self._schedule)
        self.sp_grow.valueChanged.connect(self._schedule)
        v.addWidget(c)

        c, b = card("Result")
        self.lbl_result = note("Draw a boundary to see the profile.", "Hint")
        b.addWidget(self.lbl_result)
        self.btn_save = primary_button("Save workbook...")
        self.btn_save.clicked.connect(self.save_workbook)
        self.btn_use = QPushButton("Use in Diffusor...")
        self.btn_use.clicked.connect(self.use_in_diffusor)
        self.btn_use.setVisible(self.owner is not None and hasattr(self.owner, "load_image_table"))
        b.addWidget(row(self.btn_save, self.btn_use))
        reuse = ghost_button("Reuse settings from a workbook...")
        reuse.clicked.connect(self.reuse_settings)
        b.addWidget(reuse)
        v.addWidget(c)
        v.addStretch(1)
        return box

    # ================================================================ image
    def open_image(self, path: Optional[str] = None):
        if not path:
            path, _ = QFileDialog.getOpenFileName(self, "Open image", "", imgio.FILE_FILTER)
            if not path:
                return
        p = Path(path)
        try:
            image = imgio.load_image(p)
        except ValueError as exc:
            if QMessageBox.question(self, "Unknown format",
                                    f"{exc}\n\nRead it as raw binary?") != QMessageBox.Yes:
                return
            dlg = RawDialog(p, self)
            if dlg.exec() != QDialog.Accepted:
                return
            try:
                image = imgio.load_image(p, raw=dlg.options())
            except Exception as exc2:
                QMessageBox.critical(self, "Could not read the file", str(exc2))
                return
        except Exception:
            QMessageBox.critical(self, "Could not read the image", traceback.format_exc())
            return
        self.set_image(image)

    def set_image(self, image: imgio.LoadedImage):
        same_size = (self.image is not None and self.image.data.shape[:2] == image.data.shape[:2])
        self.image = image
        if not same_size:
            self.guideline, self.exclusions, self.legend_pts = [], [], None
        self.current_poly, self.clicks = [], []
        self.color_scale = None
        self.lbl_file.setText(f"<b>{Path(image.source).name}</b><br>{image.describe()}"
                              + "".join(f"<br>{n}" for n in image.notes))
        if image.pixel_size_um:
            self.sp_px.setValue(image.pixel_size_um)
            self.lbl_px.setText(f"From {image.pixel_size_source}. Check it against the scale "
                                "bar.")
        else:
            self.lbl_px.setText("Not found in the file. Type it, or measure the scale bar with "
                                "the Scale bar tool.")
        self.cmb_mode.blockSignals(True)
        self.cmb_mode.clear()
        if image.n_channels == 1:
            self.cmb_mode.addItem("Grey value", "luminance")
        else:
            if image.is_rgb:
                self.cmb_mode.addItem(imgio.VALUE_MODES["luminance"], "luminance")
                self.cmb_mode.addItem(imgio.VALUE_MODES["mean"], "mean")
            self.cmb_mode.addItem(imgio.VALUE_MODES["channel"], "channel")
            if image.is_rgb:
                self.cmb_mode.addItem(imgio.VALUE_MODES["colour_scale"], "colour_scale")
        self.cmb_mode.blockSignals(False)
        self.cmb_channel.blockSignals(True)
        self.cmb_channel.clear()
        self.cmb_channel.addItems(image.channel_names)
        self.cmb_channel.blockSignals(False)
        lo, hi = image.value_range
        self.sp_low.setValue(lo)
        self.sp_high.setValue(hi)
        self._mode_changed()

    def _mode(self) -> str:
        return self.cmb_mode.currentData() or "luminance"

    def _mode_widgets(self):
        mode = self._mode()
        self.fld_channel.setVisible(mode == "channel")
        self.card_legend.setVisible(mode == "colour_scale")
        self.tool_buttons["legend"].setEnabled(mode == "colour_scale")

    def _mode_changed(self, *_):
        self._mode_widgets()
        if self._mode() == "colour_scale":
            self._legend_changed()
        else:
            self._values_changed()

    def _legend_changed(self, *_):
        if self.image is None or self._mode() != "colour_scale":
            return
        src = self.cmb_legend.currentData()
        v0, v1 = self.sp_v0.value(), self.sp_v1.value()
        try:
            if src == "drawn":
                if not self.legend_pts:
                    self.color_scale = None
                    self.lbl_legend.setText("Draw the legend with the Legend tool.")
                    self._set_tool("legend")
                    self._values_changed()
                    return
                cs = imgio.ColorScale.from_legend(self.image, *self.legend_pts, v0, v1,
                                                  units=self.le_units.text().strip(),
                                                  log=self.chk_log.isChecked())
            else:
                cs = imgio.ColorScale.from_colormap(src, v0, v1, units=self.le_units.text().strip(),
                                                    log=self.chk_log.isChecked())
            cs.max_distance = self.sp_de.value()
            self.color_scale = cs
            self.lbl_legend.setText(f"{len(cs.values)} legend colours, {cs.source}.")
        except Exception as exc:
            self.color_scale = None
            self.lbl_legend.setText(f"<span style='color:{theme.DANGER}'>{exc}</span>")
        self._values_changed()

    def _values_changed(self, *_):
        if self.image is None:
            return
        mode = self._mode()
        try:
            if mode == "colour_scale" and self.color_scale is None:
                self.values = None
            else:
                self.values = imgio.value_map(self.image, mode, max(self.cmb_channel.currentIndex(), 0),
                                              self.color_scale)
        except Exception as exc:
            self.values = None
            self.lbl_result.setText(str(exc))
        if mode == "colour_scale" and self.values is not None:
            off = float(np.mean(self.values[1]))
            self.lbl_legend.setText(self.lbl_legend.text().split("<br>")[0]
                                    + f"<br>{100 * off:.1f}% of the image is off the legend.")
        self._draw_image()
        self._schedule()

    def _draw_image(self):
        self.fig.clear()
        self.ax = self.fig.add_subplot(111)
        self.ax.set_axis_off()
        self._overlay = []
        img = self.image
        if img is None:
            self.canvas.draw_idle()
            return
        # a colour map is shown as it is, so the legend can be found; anything else is
        # shown as the grey values the profile will actually read
        if img.is_rgb and self._mode() == "colour_scale":
            self.ax.imshow(img.rgb01(), interpolation="nearest")
        else:
            data = self.values[0] if self.values is not None else img.channel(0)
            finite = data[np.isfinite(data)]
            lo, hi = (np.percentile(finite, [0.5, 99.5]) if finite.size else (0, 1))
            self.ax.imshow(data, cmap="gray", vmin=lo, vmax=hi if hi > lo else lo + 1,
                           interpolation="nearest")
        self.ax.set_autoscale_on(False)
        self._draw_overlay()

    # ================================================================ mouse
    def _set_tool(self, key: str):
        self.tool = key
        self.clicks = []
        for k, b in self.tool_buttons.items():
            b.setChecked(k == key)
        self.lbl_tool.setText(TOOLS[key])

    def _on_click(self, ev):
        if self.image is None or ev.inaxes is not self.ax or self.toolbar.mode:
            return
        if ev.xdata is None:
            return
        pt = (float(ev.xdata), float(ev.ydata))
        right = ev.button == 3
        if self.tool == "boundary":
            if right:
                self.undo_point()
                return
            if ev.dblclick:
                return
            self.guideline.append(pt)
            self._schedule()
        elif self.tool == "exclude":
            if right or ev.dblclick:
                if len(self.current_poly) >= 3:
                    self.exclusions.append(list(self.current_poly))
                self.current_poly = []
                self._schedule()
            else:
                self.current_poly.append(pt)
                self._draw_overlay()
        elif self.tool in ("scale", "legend") and not right:
            self.clicks.append(pt)
            if len(self.clicks) == 2:
                a, b = self.clicks
                self.clicks = []
                if self.tool == "scale":
                    self._scale_from(a, b)
                else:
                    self.legend_pts = [a, b]
                    self.cmb_legend.setCurrentIndex(0)
                    self._legend_changed()
                    self._set_tool("boundary")
            self._draw_overlay()

    def _on_key(self, ev):
        if ev.key in ("backspace", "ctrl+z"):
            self.undo_point()
        elif ev.key in ("enter", "escape") and self.current_poly:
            if len(self.current_poly) >= 3 and ev.key == "enter":
                self.exclusions.append(list(self.current_poly))
            self.current_poly = []
            self._schedule()

    def _scale_from(self, a, b):
        px = float(np.hypot(b[0] - a[0], b[1] - a[1]))
        if px < 2:
            return
        um, ok = QInputDialog.getDouble(self, "Scale bar", f"The line is {px:.1f} px long. "
                                        "How long is the scale bar (um)?", 100.0, 1e-6, 1e7, 4)
        if ok:
            self.sp_px.setValue(um / px)
            self.lbl_px.setText(f"Measured: {um:g} um over {px:.1f} px.")
            self._set_tool("boundary")

    def undo_point(self):
        if self.tool == "exclude" and self.current_poly:
            self.current_poly.pop()
            self._draw_overlay()
        elif self.guideline:
            self.guideline.pop()
            self._schedule()

    def clear_boundary(self):
        self.guideline = []
        self.extraction = None
        self._schedule()

    def clear_exclusions(self):
        self.exclusions, self.current_poly = [], []
        self._schedule()

    # ================================================================ extraction
    def _clip_changed(self, *_):
        self.sp_k.setValue(1.0 if self.cmb_clip.currentData() == "nidis" else 3.0)
        self.sp_k.setEnabled(self.cmb_clip.currentData() != "none")
        self._schedule()

    def suggest_limits(self):
        if self.values is None:
            return
        lo, hi = suggest_limits(self.values[0])
        if lo is not None:
            self.sp_low.setValue(lo)
            self.sp_high.setValue(hi)
            self.chk_limits.setChecked(True)

    def settings(self) -> ExtractionSettings:
        cl = CleaningSettings(
            low=self.sp_low.value() if self.chk_limits.isChecked() else None,
            high=self.sp_high.value() if self.chk_limits.isChecked() else None,
            grow_px=self.sp_grow.value(),
            exclusions=[list(p) for p in self.exclusions],
            clip_method=self.cmb_clip.currentData(), clip_k=self.sp_k.value(),
            max_rejected_fraction=self.sp_drop.value() / 100.0 if self.chk_drop.isChecked() else None)
        return ExtractionSettings(
            guideline=list(self.guideline), length_before_px=self.sp_before.value(),
            length_after_px=self.sp_after.value(), line_spacing_px=self.sp_spacing.value(),
            sample_step_px=self.sp_step.value(), interpolation=self.cmb_interp.currentData(),
            flip=self.chk_flip.isChecked(), value_mode=self._mode(),
            channel=max(self.cmb_channel.currentIndex(), 0),
            pixel_size_um=self.sp_px.value() or None, cleaning=cl)

    def apply_settings(self, s: ExtractionSettings):
        """Put saved settings back into the controls (NIDIS's "use existing coordinates")."""
        self.guideline = [tuple(p) for p in s.guideline]
        self.exclusions = [list(map(tuple, p)) for p in s.cleaning.exclusions]
        for sp, val in ((self.sp_before, s.length_before_px), (self.sp_after, s.length_after_px),
                        (self.sp_spacing, s.line_spacing_px), (self.sp_step, s.sample_step_px),
                        (self.sp_k, s.cleaning.clip_k)):
            sp.setValue(val)
        if s.pixel_size_um:
            self.sp_px.setValue(s.pixel_size_um)
        self.cmb_interp.setCurrentIndex(max(self.cmb_interp.findData(s.interpolation), 0))
        self.chk_flip.setChecked(s.flip)
        self.cmb_clip.setCurrentIndex(max(self.cmb_clip.findData(s.cleaning.clip_method), 0))
        self.sp_k.setValue(s.cleaning.clip_k)
        self.chk_limits.setChecked(s.cleaning.low is not None or s.cleaning.high is not None)
        if s.cleaning.low is not None:
            self.sp_low.setValue(s.cleaning.low)
        if s.cleaning.high is not None:
            self.sp_high.setValue(s.cleaning.high)
        self.sp_grow.setValue(s.cleaning.grow_px)
        self.chk_drop.setChecked(s.cleaning.max_rejected_fraction is not None)
        if s.cleaning.max_rejected_fraction is not None:
            self.sp_drop.setValue(100 * s.cleaning.max_rejected_fraction)
        i = self.cmb_mode.findData(s.value_mode)
        if i >= 0:
            self.cmb_mode.setCurrentIndex(i)
        if s.value_mode == "channel" and s.channel < self.cmb_channel.count():
            self.cmb_channel.setCurrentIndex(s.channel)
        self._schedule()

    def _schedule(self, *_):
        self._timer.start()

    def recompute(self):
        self.extraction = None
        if self.image is not None and self.values is not None and len(self.guideline) >= 2:
            try:
                self.extraction = extract_profiles(self.image, self.settings(), self.color_scale,
                                                   values=self.values)
            except Exception as exc:
                self.lbl_result.setText(f"<span style='color:{theme.DANGER}'>{exc}</span>")
        self._update_line_hint()
        self._draw_overlay()
        self._draw_preview()
        self._update_result()

    def _update_line_hint(self):
        px = self.sp_px.value()
        if px > 0:
            self.lbl_lines.setText(
                f"{self.sp_before.value() * px:.3g} um before and {self.sp_after.value() * px:.3g} "
                f"um after the boundary, a value every {self.sp_step.value() * px:.3g} um.")
        else:
            self.lbl_lines.setText("NIDIS uses 50 px either side and one line per pixel.")

    def _update_result(self):
        ex = self.extraction
        self.btn_save.setEnabled(ex is not None)
        self.btn_use.setEnabled(ex is not None)
        if ex is None:
            if self.image is not None and len(self.guideline) < 2 and \
                    not self.lbl_result.text().startswith("<span"):
                self.lbl_result.setText("Draw a boundary (at least two points) to see the "
                                        "profile.")
            return
        counts = ex.reason_counts()
        inside = ex.raw.size - counts["outside the image"]
        kept = counts["kept"]
        parts = [f"<b>{ex.n_lines}</b> lines of <b>{ex.n_samples}</b> values. "
                 f"{kept:,} of {inside:,} values kept ({100 * kept / max(inside, 1):.1f}%)."]
        for k, n in counts.items():
            if n and k not in ("kept",):
                parts.append(f"{n:,} {k}")
        parts += ex.notes
        self.lbl_result.setText("<br>".join(parts))

    # ================================================================ drawing
    def _clear_overlay(self):
        for a in self._overlay:
            try:
                a.remove()
            except (ValueError, NotImplementedError):
                pass
        self._overlay = []

    def _draw_overlay(self):
        self._clear_overlay()
        if self.image is None:
            self.canvas.draw_idle()
            return
        ax, o = self.ax, self._overlay
        kw = dict(scalex=False, scaley=False)
        for poly in self.exclusions:
            o.append(ax.add_patch(Polygon(poly, closed=True, fc=theme.DANGER, alpha=0.25,
                                          ec=theme.DANGER, lw=1)))
        if self.current_poly:
            xs, ys = zip(*self.current_poly)
            o += ax.plot(xs, ys, "--o", color=theme.DANGER, ms=3, lw=1, **kw)
        ex = self.extraction
        if ex is not None:
            n = ex.n_lines
            idx = np.unique(np.linspace(0, n - 1, min(n, MAX_DRAWN_LINES)).astype(int))
            for i in idx:
                o += ax.plot([ex.xs[0, i], ex.xs[-1, i]], [ex.ys[0, i], ex.ys[-1, i]],
                             color="#F2C14E", lw=0.6, alpha=0.7, **kw)
            o += ax.plot(ex.xs[0], ex.ys[0], color="#F2C14E", lw=1.6, **kw)      # start side
            o += ax.plot(ex.xs[-1], ex.ys[-1], color="#F2C14E", lw=0.8, ls=":", **kw)
            if self.chk_show.isChecked():
                bad = (ex.reason != KEPT) & (ex.reason != OUTSIDE)
                ys, xs = ex.ys[bad], ex.xs[bad]
                if xs.size > MAX_REJECTED_DOTS:
                    pick = np.random.default_rng(0).choice(xs.size, MAX_REJECTED_DOTS, replace=False)
                    xs, ys = xs[pick], ys[pick]
                if xs.size:
                    o.append(ax.scatter(xs, ys, s=2, c="#FF3B30", lw=0, alpha=0.8))
        if self.guideline:
            gx, gy = zip(*self.guideline)
            o += ax.plot(gx, gy, "-o", color="#2FD1C5", ms=4, lw=1.6, **kw)
            o += ax.plot(gx[:1], gy[:1], "o", color="#2FD1C5", ms=9, mfc="none", mew=2, **kw)
            o.append(ax.annotate("start", (gx[0], gy[0]), xytext=(6, 6),
                                 textcoords="offset points", color="#2FD1C5", fontsize=9))
        if self.legend_pts and self._mode() == "colour_scale" \
                and self.cmb_legend.currentData() == "drawn":
            (x0, y0), (x1, y1) = self.legend_pts
            o += ax.plot([x0, x1], [y0, y1], "-", color="white", lw=2.5, **kw)
            o += ax.plot([x0, x1], [y0, y1], "-", color="black", lw=1, **kw)
            o.append(ax.annotate("low", (x0, y0), xytext=(5, -10), textcoords="offset points",
                                 fontsize=8, color="black", backgroundcolor="white"))
        for p in self.clicks:
            o += ax.plot([p[0]], [p[1]], "+", color="#FF3B30", ms=12, mew=2, **kw)
        self.canvas.draw_idle()

    def _draw_preview(self):
        self.pfig.clear()
        ax = self.pax = self.pfig.add_subplot(111)
        ex = self.extraction
        if ex is None:
            ax.text(0.5, 0.5, "The averaged profile appears here", ha="center", va="center",
                    color=theme.TEXT_MUTED, transform=ax.transAxes)
            theme.apply_plot_style(self.pfig, [ax])
            self.pcanvas.draw_idle()
            return
        um = ex.distance_um()
        x = um if um is not None else ex.distance_px
        raw, clean = ex.statistics("raw"), ex.statistics("clean")
        ax.plot(x, raw["Mean"], "o", ms=3, color=theme.PLOT_INITIAL, alpha=0.7, label="raw mean")
        ax.fill_between(x, clean["Mean"] - clean["SD"], clean["Mean"] + clean["SD"],
                        color=theme.ACCENT, alpha=0.15, lw=0, label="cleaned, +/- 1 SD")
        ax.plot(x, clean["Mean"], "-", color=theme.ACCENT, lw=1.6, label="cleaned mean")
        b = ex.settings.length_before_px * (ex.settings.pixel_size_um or 1.0)
        ax.axvline(b, color=theme.BORDER_STRONG, ls="--", lw=1)
        ax.set_xlabel("distance (um)" if um is not None else "distance (px)")
        ax.set_ylabel(ex.value_label)
        ax.legend(fontsize=8, frameon=False, loc="best")
        theme.apply_plot_style(self.pfig, [ax])
        self.pcanvas.draw_idle()

    def overlay_png(self) -> bytes:
        buf = io.BytesIO()
        self.fig.savefig(buf, format="png", dpi=110, facecolor=theme.SURFACE)
        return buf.getvalue()

    # ================================================================ output
    def save_workbook(self, path: Optional[str] = None) -> Optional[str]:
        if self.extraction is None:
            return None
        if not path:
            src = Path(self.image.source)
            n = 1
            while (src.parent / f"{src.stem}_profile{n}.xlsx").exists():
                n += 1
            path, _ = QFileDialog.getSaveFileName(self, "Save workbook",
                                                  str(src.parent / f"{src.stem}_profile{n}.xlsx"),
                                                  "Excel workbook (*.xlsx)")
            if not path:
                return None
        try:
            out = write_workbook(path, self.extraction, self.overlay_png())
        except Exception:
            QMessageBox.critical(self, "Could not save", traceback.format_exc())
            return None
        self.lbl_result.setText(self.lbl_result.text() + f"<br>Saved to {Path(out).name}.")
        return out

    def use_in_diffusor(self):
        if self.extraction is None or self.owner is None:
            return
        self.owner.load_image_table(table_from_extraction(self.extraction))

    def reuse_settings(self):
        path, _ = QFileDialog.getOpenFileName(self, "Settings from an earlier extraction", "",
                                              "Image-profile workbook (*.xlsx)")
        if not path:
            return
        try:
            if not is_extraction_workbook(path):
                raise ValueError("this is not a Diffusor image-profile workbook")
            from ..dataio.image_profiles import read_extraction
            s = read_extraction(path).settings
            if s is None:
                raise ValueError("the workbook holds no settings")
            self.apply_settings(s)
        except Exception as exc:
            QMessageBox.warning(self, "Could not reuse the settings", str(exc))
