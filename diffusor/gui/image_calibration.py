"""Turn an image profile into a composition profile: pixel size and value map.

Opened when an image-profile workbook is loaded, or from the extractor's "Use
in Diffusor" button. The dialog asks for the two things an image cannot tell
Diffusor: how many micrometres a pixel is, and which composition a grey (or
channel, or legend) value stands for. It then hands the main window an ordinary
table and column mapping, so everything after loading works as for any file.
"""
from __future__ import annotations

import traceback
from pathlib import Path
from typing import Optional

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox,
                               QFileDialog, QHBoxLayout, QHeaderView, QInputDialog, QLineEdit,
                               QMessageBox, QPushButton, QStackedWidget, QTableWidget,
                               QTableWidgetItem, QVBoxLayout, QWidget)

from ..dataio import ProfileSpec, anchors_from_microprobe, calibrate, read_table
from ..dataio.image_profiles import ExtractionTable, composition_table
from .widgets import callout, field, fit_to_screen, ghost_button, note, pair, row, scrollable

STATISTICS = {"Clean_Mean": "Cleaned mean (recommended)", "Clean_Median": "Cleaned median",
              "Raw_Mean": "Raw mean, nothing removed"}
UNCERTAINTIES = [("SE", "Standard error of the mean (NIDIS)"),
                 ("SD", "Standard deviation across the lines"), (None, "None")]
MAPS = [("two", "Linear, through two reference points"),
        ("anchors", "Fitted to microprobe anchor points"),
        ("none", "None, the values already are compositions")]


def _spin(lo, hi, val, dec=4, suffix=""):
    s = QDoubleSpinBox()
    s.setRange(lo, hi)
    s.setDecimals(dec)
    s.setValue(val)
    if suffix:
        s.setSuffix(suffix)
    return s


class ImageCalibrationDialog(QDialog):
    """Pixel size and value-to-composition map for an :class:`ExtractionTable`."""

    def __init__(self, table: ExtractionTable, parent=None):
        super().__init__(parent)
        self.table = table
        self.frame = None
        self.spec: Optional[ProfileSpec] = None
        self.calibration = None
        self.setWindowTitle("Use an image profile")
        fit_to_screen(self, 640, 760)
        outer = QVBoxLayout(self)
        outer.setContentsMargins(12, 12, 12, 12)
        # the form scrolls, the buttons stay put, so they are reachable on any screen
        form = QWidget()
        form.setObjectName("Page")
        v = QVBoxLayout(form)
        v.setContentsMargins(4, 0, 8, 0)
        v.setSpacing(10)
        src = Path(table.source).name or Path(table.path).name
        v.addWidget(note(f"<b>{len(table.profile)}</b> positions from <b>{src}</b>. "
                         f"The values are {table.value_label}s.", "Sub"))

        self.sp_px = _spin(0.0, 1e4, table.pixel_size_um or 0.0, 6, " um/px")
        self.sp_px.setSpecialValueText("not set")
        hint = ("Taken from the workbook." if table.pixel_size_um else
                "Not in the workbook. Read it from the instrument, or measure the scale bar "
                "in the extractor.")
        v.addWidget(field("Pixel size", self.sp_px, hint))

        self.cmb_stat = QComboBox()
        for key in table.statistics:
            self.cmb_stat.addItem(STATISTICS[key], key)
        self.cmb_unc = QComboBox()
        for key, label in UNCERTAINTIES:
            self.cmb_unc.addItem(label, key)
        v.addWidget(pair(field("Profile value", self.cmb_stat),
                         field("Uncertainty", self.cmb_unc)))
        v.addWidget(note("Neighbouring lines share pixels, so the standard error is a lower "
                         "bound. The standard deviation across the lines is the cautious "
                         "choice.", "Hint"))

        self.cmb_map = QComboBox()
        for key, label in MAPS:
            self.cmb_map.addItem(label, key)
        v.addWidget(field("Value to composition", self.cmb_map))
        self.stack = QStackedWidget()
        self.stack.addWidget(self._two_point_page())
        self.stack.addWidget(self._anchor_page())
        self.stack.addWidget(note("The profile values are modelled as they are. Use this for a "
                                  "quantitative map, or a colour legend drawn in composition "
                                  "units.", "Hint"))
        v.addWidget(self.stack)

        self.le_name = QLineEdit("X_Fe")
        v.addWidget(field("Name of the composition", self.le_name,
                          "Only a label for the column and the plot."))
        self.box = callout("Info")
        v.addWidget(self.box)
        v.addStretch(1)
        outer.addWidget(scrollable(form), 1)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.button(QDialogButtonBox.Ok).setText("Load profile")
        bb.accepted.connect(self._accept)
        bb.rejected.connect(self.reject)
        outer.addWidget(bb)

        self.cmb_map.currentIndexChanged.connect(self._map_changed)
        for w in (self.sp_px, self.sp_g1, self.sp_c1, self.sp_g2, self.sp_c2):
            w.valueChanged.connect(self._update)
        for w in (self.cmb_stat, self.cmb_unc, self.cmb_degree):
            w.currentIndexChanged.connect(self._update)
        self.tbl.itemChanged.connect(self._update)
        if table.value_label.startswith("legend value"):
            # a legend is usually drawn in composition units already
            self.cmb_map.setCurrentIndex(self.cmb_map.findData("none"))
        self._map_changed()

    # ------------------------------------------------------------ pages
    def _values(self) -> np.ndarray:
        return np.asarray(self.table.profile[self.cmb_stat.currentData()], float)

    def _two_point_page(self) -> QWidget:
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setContentsMargins(0, 0, 0, 0)
        g = self.table.profile[self.table.statistics[0]]
        lo, hi = float(np.nanmin(g)), float(np.nanmax(g))
        self.sp_g1, self.sp_g2 = _spin(-1e9, 1e9, lo, 3), _spin(-1e9, 1e9, hi, 3)
        self.sp_c1, self.sp_c2 = _spin(-1e9, 1e9, 0.0, 5), _spin(-1e9, 1e9, 1.0, 5)
        lay.addWidget(pair(field("Value 1", self.sp_g1), field("is composition", self.sp_c1)))
        lay.addWidget(pair(field("Value 2", self.sp_g2), field("is composition", self.sp_c2)))
        lay.addWidget(note("Two points fix a straight line but say nothing about its "
                           "uncertainty. Use anchor points when you have more.", "Hint"))
        return w

    def _anchor_page(self) -> QWidget:
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setContentsMargins(0, 0, 0, 0)
        self.tbl = QTableWidget(6, 2)
        self.tbl.setHorizontalHeaderLabels(["Profile value", "Composition"])
        self.tbl.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.tbl.setMinimumHeight(170)
        lay.addWidget(self.tbl)
        self.cmb_degree = QComboBox()
        self.cmb_degree.addItem("Linear", 1)
        self.cmb_degree.addItem("Quadratic", 2)
        add = ghost_button("Add a row")
        add.clicked.connect(lambda: self.tbl.insertRow(self.tbl.rowCount()))
        probe = QPushButton("Fill from a microprobe traverse...")
        probe.clicked.connect(self._fill_from_probe)
        lay.addWidget(row(self.cmb_degree, add, probe))
        lay.addWidget(note("A microprobe traverse along the same line gives anchors: the "
                           "profile is averaged around each spot. Its distances must use the "
                           "same origin as the image profile.", "Hint"))
        return w

    def _map_changed(self):
        self.stack.setCurrentIndex(self.cmb_map.currentIndex())
        self._update()

    def _anchor_pairs(self):
        g, c = [], []
        for r in range(self.tbl.rowCount()):
            items = [self.tbl.item(r, k) for k in (0, 1)]
            try:
                gv, cv = (float(it.text()) for it in items)
            except (AttributeError, ValueError):
                continue
            g.append(gv)
            c.append(cv)
        return np.array(g), np.array(c)

    def _fill_from_probe(self):
        px = self.sp_px.value()
        if px <= 0:
            QMessageBox.information(self, "Pixel size needed",
                                    "Set the pixel size first, so the probe distances can be "
                                    "matched to the image profile.")
            return
        path, _ = QFileDialog.getOpenFileName(self, "Microprobe traverse", "",
                                              "Tables (*.csv *.txt *.tsv *.xlsx *.xls)")
        if not path:
            return
        try:
            df = read_table(path)
            numeric = [str(c) for c in df.columns if np.issubdtype(df[c].dtype, np.number)]
            if len(numeric) < 2:
                raise ValueError("the table needs a distance and a composition column")
            dist, ok = QInputDialog.getItem(self, "Distance column",
                                            "Distance along the profile (um):", numeric, 0, False)
            if not ok:
                return
            comp, ok = QInputDialog.getItem(self, "Composition column", "Composition:",
                                            [c for c in numeric if c != dist], 0, False)
            if not ok:
                return
            win, ok = QInputDialog.getDouble(self, "Averaging window",
                                             "Average the image profile within +/- (um):",
                                             1.0, 0.0, 100.0, 2)
            if not ok:
                return
            x_img = np.asarray(self.table.profile["Distance_px"], float) * px
            g, c = anchors_from_microprobe(x_img, self._values(), df[dist], df[comp], win)
            keep = np.isfinite(g) & np.isfinite(c)
            self.tbl.blockSignals(True)
            self.tbl.setRowCount(max(int(keep.sum()), 2))
            for r, (gv, cv) in enumerate(zip(g[keep], c[keep])):
                self.tbl.setItem(r, 0, QTableWidgetItem(f"{gv:.4f}"))
                self.tbl.setItem(r, 1, QTableWidgetItem(f"{cv:.6g}"))
            self.tbl.blockSignals(False)
            self.le_name.setText(comp)
            self._update()
        except Exception as exc:
            QMessageBox.warning(self, "Could not read the traverse", str(exc))

    # ------------------------------------------------------------ result
    def _calibration(self):
        kind = self.cmb_map.currentData()
        if kind == "none":
            return None
        if kind == "two":
            g1, g2 = self.sp_g1.value(), self.sp_g2.value()
            if g1 == g2:
                raise ValueError("the two reference values must differ")
            return calibrate([g1, g2], [self.sp_c1.value(), self.sp_c2.value()], 1)
        g, c = self._anchor_pairs()
        return calibrate(g, c, int(self.cmb_degree.currentData()))

    def _uncertainty_column(self) -> Optional[str]:
        kind = self.cmb_unc.currentData()
        if kind is None:
            return None
        prefix = "Raw_" if self.cmb_stat.currentData().startswith("Raw") else "Clean_"
        return prefix + kind

    def _update(self, *_):
        try:
            cal = self._calibration()
            if cal is None:
                text = "No calibration: the profile values are used as compositions."
            else:
                text = cal.describe() + "".join(f"<br>Note: {n}" for n in cal.notes)
            if self.sp_px.value() <= 0:
                text += "<br><b>Set the pixel size to continue.</b>"
            self.box.label.setText(text)
            self.box.setObjectName("InfoBox")
        except Exception as exc:
            self.box.label.setText(str(exc))
            self.box.setObjectName("WarnBox")
        self.box.style().unpolish(self.box)
        self.box.style().polish(self.box)

    def _accept(self):
        try:
            cal = self._calibration()
            name = self.le_name.text().strip() or "Composition"
            name = "".join(ch if ch.isalnum() or ch in "_#" else "_" for ch in name)
            frame = composition_table(self.table, self.sp_px.value() or None,
                                      self.cmb_stat.currentData(), self._uncertainty_column(),
                                      cal, name)
            if len(frame) < 3:
                raise ValueError("fewer than three positions have a value")
        except Exception as exc:
            QMessageBox.warning(self, "Cannot use the profile", str(exc))
            return
        self.frame, self.calibration = frame, cal
        err = name + "_err"
        self.spec = ProfileSpec(distance_column="Distance_um", column_a=name,
                                sigma_a_column=err if err in frame else None,
                                distance_unit="um", mode="A",
                                label=f"image profile, {self.cmb_stat.currentText().lower()}")
        self.accept()

    def summary(self) -> str:
        """What was chosen, for the log."""
        lines = [f"image profile {self.table.path or self.table.source}",
                 f"pixel size {self.sp_px.value():.6g} um/px, "
                 f"{self.cmb_stat.currentText().lower()}, uncertainty "
                 f"{self.cmb_unc.currentText().lower()}"]
        lines.append("calibration: " + (self.calibration.describe() if self.calibration
                                        else "none, values used as compositions"))
        return "\n".join(lines)
