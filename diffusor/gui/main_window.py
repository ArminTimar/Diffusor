"""Diffusor's main window: a stepped flow ending in a results view.

The settings are collected one group at a time, then summarised in a narrow
sidebar next to the plot. Any group in that summary can be clicked to go back
and change it.
"""
from __future__ import annotations

import traceback
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
from PySide6.QtCore import Qt, Slot
from PySide6.QtGui import QAction, QFont, QGuiApplication
from PySide6.QtWidgets import (QApplication, QCheckBox, QComboBox, QDialog, QDialogButtonBox,
                               QDoubleSpinBox, QFileDialog, QFrame, QHBoxLayout, QHeaderView, QLabel,
                               QLineEdit, QListWidget, QListWidgetItem, QMainWindow,
                               QMessageBox, QProgressBar, QPushButton, QSizePolicy, QSpinBox,
                               QSplitter, QStackedWidget, QTableWidget, QTableWidgetItem,
                               QTextEdit, QVBoxLayout, QWidget)

from .. import datasets as ds
from ..coefficients import Conditions, get as get_coefficient, list_coefficients
from ..coefficients.plagioclase import ACTIVITY_A, activity_theta
from ..dataio import (ProfileSpec, build_profile, methods_paragraph, read_table,
                      save_results, suggest_spec)
from ..fitting import DiffusionModel, UncertaintyBudget
from ..minerals import MINERALS, get_mineral
from ..references import cite, format_reference
from ..solvers import Geometry, InitialCondition, dirichlet, zero_flux
from ..solvers.history import ThermalHistory
from ..solvers.initial import guess_step_from_data
from ..thermo import available_buffers, log_fo2_from_delta
from ..thermo.units import human_time
from . import format_help, theme
from .plot_widget import ProfilePlot
from .widgets import (WrapLabel, card, collapsible, divider, field, ghost_button,
                      install_no_wheel, note, pair, page_body, primary_button, row,
                      scrollable)
from .workers import CompareWorker, FitWorker, MonteCarloWorker, start

STEPS = ["Data", "Mineral", "Conditions", "Model", "Coefficient", "Uncertainty", "Results"]
RESULTS = len(STEPS) - 1

# Analytical resolution presets: (label, width kind, default width in um, sigma for a
# fixed preset, hint). Width kind is "spot" (sigma = d / 4 for an evenly lit round
# spot), "slit" (sigma = w / sqrt 12 for an evenly lit slit), or None.
RESOLUTION_PRESETS = [
    ("No correction", None, 0.0, 0.0,
     "The model is compared with the data as it is."),
    ("BSE image profile", None, 0.0, 0.0,
     "Grey-value profiles resolve better than 0.5 um (Petrone et al. 2016), so a "
     "correction is rarely needed."),
    ("Microprobe, focused beam", None, 0.0, 0.6,
     "Ganguly et al. (1988) found sigma rarely exceeds 0.6 um on a modern microprobe. "
     "Profiles longer than about 15 um are barely affected."),
    ("Microprobe, defocused beam", "spot", 5.0, None,
     "Feldspar is often measured with a 5 um beam (Chamberlain et al. 2014, "
     "Grocolas et al. 2025)."),
    ("LA-ICP-MS spot", "spot", 10.0, None,
     "Druitt et al. (2012) used a 10 um laser spot."),
    ("LA-ICP-MS line scan", "slit", 7.5, None,
     "Grocolas et al. (2025) scanned with a 7.5 um wide slit."),
    ("SIMS spot", "spot", 12.0, None,
     "Druitt et al. (2012) used 10 to 15 um ion beams."),
    ("Custom sigma", None, 0.0, None,
     "Standard deviation of the Gaussian beam profile (Ganguly et al. 1988)."),
]


def _spin(lo, hi, val, dec=3, step=1.0, suffix=""):
    s = QDoubleSpinBox()
    s.setRange(lo, hi)
    s.setDecimals(dec)
    s.setSingleStep(step)
    s.setValue(val)
    if suffix:
        s.setSuffix(suffix)
    return s


def _text_dialog(parent, title, text, width=900, height=600):
    dlg = QDialog(parent)
    dlg.setWindowTitle(title)
    screen = QGuiApplication.primaryScreen().availableGeometry()
    dlg.resize(min(width, int(screen.width() * 0.9)), min(height, int(screen.height() * 0.85)))
    lay = QVBoxLayout(dlg)
    lay.setContentsMargins(18, 18, 18, 18)
    te = QTextEdit()
    te.setReadOnly(True)
    te.setPlainText(text)
    te.setFont(QFont(theme.MONO_STACK.split(",")[0].strip('"'), 9))
    lay.addWidget(te)
    bb = QDialogButtonBox(QDialogButtonBox.Close)
    bb.rejected.connect(dlg.reject)
    lay.addWidget(bb)
    dlg.exec()


def _repolish(w, name):
    w.setObjectName(name)
    w.style().unpolish(w)
    w.style().polish(w)


# ==================================================================== columns
class ColumnDialog(QDialog):
    """Confirm the guessed column mapping. Everything else sits under Advanced."""

    def __init__(self, df, spec: ProfileSpec, parent=None):
        super().__init__(parent)
        self.df = df
        self.setWindowTitle("Check the columns")
        cols = [str(c) for c in df.columns]
        numeric = [str(c) for c in df.columns if np.issubdtype(df[c].dtype, np.number)]
        root = QVBoxLayout(self)
        root.setContentsMargins(20, 16, 20, 16)
        root.setSpacing(10)
        head = QLabel("Is this right?")
        head.setObjectName("H2")
        root.addWidget(head)

        self.dist = QComboBox(); self.dist.addItems(numeric or cols)
        self.unit = QComboBox(); self.unit.addItems(["um", "mm", "nm", "m"])
        self.a = QComboBox(); self.a.addItems(numeric or cols)
        self.b = QComboBox(); self.b.addItems(["(none)"] + (numeric or cols))
        self.dist.setCurrentText(str(spec.distance_column))
        self.unit.setCurrentText(spec.distance_unit)
        self.a.setCurrentText(str(spec.column_a))
        self.b.setCurrentText(str(spec.column_b) if spec.column_b else "(none)")
        root.addWidget(pair(field("Distance", self.dist), field("Unit", self.unit)))
        root.addWidget(pair(field("Element A", self.a), field("Element B", self.b)))
        self.lbl_what = note("", "Good")
        root.addWidget(self.lbl_what)

        self.preview = QTableWidget(3, 3)
        self.preview.verticalHeader().setVisible(False)
        self.preview.setEditTriggers(QTableWidget.NoEditTriggers)
        self.preview.setFocusPolicy(Qt.NoFocus)
        self.preview.setFixedHeight(118)
        self.preview.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        root.addWidget(self.preview)

        # --- advanced
        adv = QWidget()
        av = QVBoxLayout(adv)
        av.setContentsMargins(0, 0, 8, 0)
        av.setSpacing(10)
        self.sa = QComboBox(); self.sa.addItems(["(none)"] + numeric)
        self.sb = QComboBox(); self.sb.addItems(["(none)"] + numeric)
        self.sa.setCurrentText(spec.sigma_a_column or "(none)")
        self.sb.setCurrentText(spec.sigma_b_column or "(none)")
        av.addWidget(pair(field("Uncertainty of A", self.sa), field("Uncertainty of B", self.sb)))
        self.slevel = QComboBox(); self.slevel.addItems(["1s", "2s"])
        self.slevel.setCurrentText(spec.sigma_level)
        self.mode = QComboBox(); self.mode.addItems(["A/(A+B)", "B/(A+B)", "A", "A-B"])
        self.mode.setCurrentText(spec.mode)
        av.addWidget(pair(field("Uncertainty level", self.slevel),
                          field("Modelled variable", self.mode)))
        self.oxa = QLineEdit(spec.oxide_a or ""); self.oxb = QLineEdit(spec.oxide_b or "")
        self.oxa.setPlaceholderText("e.g. FeO")
        self.oxb.setPlaceholderText("e.g. MgO")
        av.addWidget(pair(field("Oxide of A", self.oxa), field("Oxide of B", self.oxb)))
        av.addWidget(note("Name oxides only for wt% oxide columns.", "Hint"))
        self.xmin = QLineEdit("" if spec.x_min is None else f"{spec.x_min:g}")
        self.xmax = QLineEdit("" if spec.x_max is None else f"{spec.x_max:g}")
        self.xmin.setPlaceholderText("first point")
        self.xmax.setPlaceholderText("last point")
        av.addWidget(pair(field("Fit from", self.xmin), field("Fit to", self.xmax)))
        av.addWidget(note("Points outside this range stay in the file but are left out of "
                          "the fit.", "Hint"))
        self.an = QComboBox(); self.an.addItems(["(none)"] + numeric)
        self.an_percent = QCheckBox("in mol%")
        self.an_percent.setChecked(True)
        av.addWidget(pair(field("Anorthite column (plagioclase)", self.an), self.an_percent))
        av.addStretch(1)
        area = scrollable(adv)
        area.setMaximumHeight(260)
        area.setMinimumHeight(200)
        root.addWidget(collapsible("Advanced settings", area))

        bb = QDialogButtonBox()
        ok = bb.addButton("Looks right", QDialogButtonBox.AcceptRole)
        ok.setObjectName("Primary")
        bb.addButton(QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        root.addWidget(bb)

        for w in (self.dist, self.unit, self.a, self.b, self.mode):
            w.currentIndexChanged.connect(self._update)
        self.oxa.textChanged.connect(self._update)
        self.oxb.textChanged.connect(self._update)
        self._update()
        screen = QGuiApplication.primaryScreen().availableGeometry()
        self.setMinimumWidth(min(560, int(screen.width() * 0.9)))
        self.resize(min(620, int(screen.width() * 0.9)), 10)
        self.setMaximumHeight(int(screen.height() * 0.92))

    @staticmethod
    def _none(text):
        text = text.strip()
        return None if text in ("", "(none)") else text

    def _update(self):
        a, b = self.a.currentText(), self._none(self.b.currentText())
        if b and self.mode.currentText() == "A":
            self.mode.setCurrentText("A/(A+B)")
        if not b and self.mode.currentText() != "A":
            self.mode.setCurrentText("A")
        mode = self.mode.currentText()
        if mode == "A":
            what = f"Models {a} as it stands."
        else:
            ox = (" in cation moles" if self.oxa.text().strip() and self.oxb.text().strip()
                  else "")
            forms = {"A/(A+B)": f"{a}/({a}+{b})", "B/(A+B)": f"{b}/({a}+{b})",
                     "A-B": f"{a} - {b}"}
            what = f"Models {forms.get(mode, mode)}{ox}."
        self.lbl_what.setText(f"{what} {len(self.df)} rows.")
        shown = [self.dist.currentText(), a] + ([b] if b else [])
        self.preview.setColumnCount(len(shown))
        self.preview.setHorizontalHeaderLabels(shown)
        for r in range(3):
            for c, name in enumerate(shown):
                v = self.df[name].iloc[r] if r < len(self.df) and name in self.df else ""
                self.preview.setItem(r, c, QTableWidgetItem(f"{v:g}" if isinstance(v, (int, float, np.number)) else str(v)))

    def spec(self) -> ProfileSpec:
        def num(edit):
            t = edit.text().strip()
            return float(t) if t else None
        return ProfileSpec(
            distance_column=self.dist.currentText(), column_a=self.a.currentText(),
            column_b=self._none(self.b.currentText()),
            sigma_a_column=self._none(self.sa.currentText()),
            sigma_b_column=self._none(self.sb.currentText()),
            distance_unit=self.unit.currentText(), mode=self.mode.currentText(),
            oxide_a=self.oxa.text().strip() or None, oxide_b=self.oxb.text().strip() or None,
            sigma_level=self.slevel.currentText(),
            x_min=num(self.xmin), x_max=num(self.xmax))

    def an_column(self):
        return self._none(self.an.currentText()), self.an_percent.isChecked()


# ==================================================================== window
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Diffusor")
        screen = QGuiApplication.primaryScreen().availableGeometry()
        self.resize(min(1440, int(screen.width() * 0.95)), min(920, int(screen.height() * 0.92)))
        self.setStyleSheet(theme.stylesheet())
        self._no_wheel = install_no_wheel(QApplication.instance())

        self.profile = None
        self.dataset: Optional[ds.ExampleDataset] = None
        self.an_values = None
        self.fit_result = None
        self.mc_result = None
        self.compare_results = None
        self._jobs: List = []           # (thread, worker) pairs kept alive until finished
        self._log_lines: List[str] = []
        self._prefill: Dict[int, QLabel] = {}
        self.step = 0

        self.plot = ProfilePlot()

        root = QWidget(); root.setObjectName("Page")
        outer = QVBoxLayout(root)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)
        outer.addWidget(self._build_header())
        outer.addWidget(divider())
        self.pages = QStackedWidget()
        for builder in (self._page_data, self._page_mineral, self._page_conditions,
                        self._page_model, self._page_coefficient, self._page_uncertainty,
                        self._page_results):
            self.pages.addWidget(builder())
        outer.addWidget(self.pages, 1)
        outer.addWidget(divider())
        outer.addWidget(self._build_footer())
        self.setCentralWidget(root)

        self.progress = QProgressBar()
        self.progress.setVisible(False)
        self.progress.setFixedWidth(200)
        self.statusBar().addPermanentWidget(self.progress)
        self._build_menu()
        self._on_mineral_changed()
        self._on_resolution_changed()
        self._go(0)

    # ================================================================ chrome
    def _build_header(self) -> QWidget:
        w = QWidget(); w.setObjectName("Page")
        h = QHBoxLayout(w)
        h.setContentsMargins(24, 12, 24, 12)
        h.setSpacing(16)
        title = QLabel("Diffusor")
        title.setObjectName("H1")
        h.addWidget(title)
        h.addStretch(1)
        self.step_labels = []
        rail = QWidget()
        rh = QHBoxLayout(rail)
        rh.setContentsMargins(0, 0, 0, 0)
        rh.setSpacing(4)
        for i, name in enumerate(STEPS):
            lab = QLabel(f"{i + 1}  {name}")
            lab.setObjectName("StepDot")
            lab.setCursor(Qt.PointingHandCursor)
            lab.mousePressEvent = (lambda e, idx=i: self._go(idx))
            self.step_labels.append(lab)
            rh.addWidget(lab)
        h.addWidget(rail)
        return w

    def _build_footer(self) -> QWidget:
        w = QWidget(); w.setObjectName("Page")
        h = QHBoxLayout(w)
        h.setContentsMargins(24, 10, 24, 10)
        h.setSpacing(10)
        self.btn_back = QPushButton("Back")
        self.btn_back.clicked.connect(lambda: self._go(self.step - 1))
        h.addWidget(self.btn_back)
        self.lbl_footer = QLabel("")
        self.lbl_footer.setObjectName("Sub")
        h.addWidget(self.lbl_footer)
        h.addStretch(1)
        self.btn_next = primary_button("Continue")
        self.btn_next.clicked.connect(self._next)
        h.addWidget(self.btn_next)
        return w

    def _build_menu(self):
        m = self.menuBar().addMenu("&File")
        for text, fn in (("Load profile...", self.load_file),
                         ("Export results...", self.export_results),
                         (None, None),
                         ("Quit", self.close)):
            if text is None:
                m.addSeparator(); continue
            a = QAction(text, self); a.triggered.connect(fn); m.addAction(a)
        v = self.menuBar().addMenu("&View")
        for text, fn in (("Methods and references", self.show_methods),
                         ("Diffusion coefficient details", lambda: self.show_coefficient_info()),
                         ("Monte Carlo histogram", self.show_histogram),
                         ("Log", self.show_log)):
            a = QAction(text, self); a.triggered.connect(fn); v.addAction(a)
        h = self.menuBar().addMenu("&Help")
        a = QAction("Input file format", self); a.triggered.connect(self.show_format); h.addAction(a)
        a = QAction("About", self); a.triggered.connect(self.show_about); h.addAction(a)

    def _prefill_label(self, step: int) -> QLabel:
        lab = note("", "Good")
        lab.setVisible(False)
        self._prefill[step] = lab
        return lab

    # ================================================================ step 1
    def _page_data(self) -> QWidget:
        load_card, body = card("Load a profile")
        btn = primary_button("Choose a file...")
        btn.clicked.connect(self.load_file)
        fmt = ghost_button("File format")
        fmt.clicked.connect(self.show_format)
        body.addWidget(row(btn, fmt))
        body.addWidget(note(format_help.SHORT, "Hint"))
        self.lbl_data = note("", "Good")
        self.lbl_data.setVisible(False)
        body.addWidget(self.lbl_data)

        ex_card, ebody = card("Or try an example")
        self.lst_examples = QListWidget()
        self.lst_examples.setMinimumHeight(200)
        for d in ds.DATASETS:
            kind = "measured" if d.kind == "measured" else "synthetic"
            it = QListWidgetItem(f"{d.name.split(' (')[0]}   ·   {kind}")
            it.setData(Qt.UserRole, d.key)
            if not d.exists:
                it.setFlags(it.flags() & ~Qt.ItemIsEnabled)
            self.lst_examples.addItem(it)
        self.lst_examples.currentItemChanged.connect(self._example_selected)
        self.lst_examples.itemDoubleClicked.connect(lambda _it: self.load_example())
        ebody.addWidget(self.lst_examples)
        self.lbl_example = note("", "Hint")
        ebody.addWidget(self.lbl_example)
        b = QPushButton("Load example")
        b.clicked.connect(self.load_example)
        info = ghost_button("Where the data come from")
        info.clicked.connect(self.show_example_details)
        ebody.addWidget(row(b, info))
        return page_body(load_card, ex_card, max_width=640)

    def _example_selected(self, item, _prev=None):
        if item is None:
            return
        d = ds.get(item.data(Qt.UserRole))
        src = (f"Measured. {cite(d.citation)}." if d.kind == "measured"
               else "Synthetic, made by Diffusor with a known answer.")
        self.lbl_example.setText(f"{src} {d.mineral}, {d.species}.")

    def show_example_details(self):
        item = self.lst_examples.currentItem()
        if item is None:
            return
        d = ds.get(item.data(Qt.UserRole))
        parts = [d.provenance_banner(), "", "Expected answer", d.expected]
        if d.notes:
            parts += ["", "Notes", d.notes]
        _text_dialog(self, d.name, "\n".join(parts), 760, 520)

    def show_format(self):
        _text_dialog(self, "Input file format", format_help.as_plain_text(), 760, 560)

    # ================================================================ step 2
    def _page_mineral(self) -> QWidget:
        c, body = card("Mineral")
        self.cmb_mineral = QComboBox()
        for k, mn in MINERALS.items():
            self.cmb_mineral.addItem(mn.name, k)
        self.cmb_mineral.currentIndexChanged.connect(self._on_mineral_changed)
        self.cmb_species = QComboBox()
        self.cmb_species.currentIndexChanged.connect(self._refresh_coefficients)
        body.addWidget(pair(field("Mineral", self.cmb_mineral),
                            field("Diffusing species", self.cmb_species)))
        self.lbl_mineral_note = note("", "Hint")
        body.addWidget(self.lbl_mineral_note)

        o, obody = card("Traverse direction")
        self.cmb_axis = QComboBox()
        self.cmb_axis.addItems(["Same axis as the experiments",
                                "a-axis", "b-axis", "c-axis", "Angles to a, b and c"])
        self.cmb_axis.currentIndexChanged.connect(self._on_axis_changed)
        obody.addWidget(self.cmb_axis)
        self.sp_alpha = _spin(0, 180, 90, 1, 1, " deg")
        self.sp_beta = _spin(0, 180, 90, 1, 1, " deg")
        self.sp_gamma = _spin(0, 180, 0, 1, 1, " deg")
        self._angle_row = QWidget()
        ar = QVBoxLayout(self._angle_row)
        ar.setContentsMargins(0, 0, 0, 0)
        ar.addWidget(row(field("to a", self.sp_alpha), field("to b", self.sp_beta),
                         field("to c", self.sp_gamma), stretch_last=True))
        ar.addWidget(note("D = Da cos²α + Db cos²β + Dc cos²γ (Costa & Chakraborty 2004)",
                          "Hint"))
        self._angle_row.setVisible(False)
        obody.addWidget(self._angle_row)

        x, xbody = card("Composition")
        self.sp_xcomp = _spin(0, 1, 0.15, 4, 0.01)
        xbody.addWidget(field("Representative value", self.sp_xcomp))
        self.chk_comp_dep = QCheckBox("D follows the composition along the profile")
        self.chk_comp_dep.setChecked(True)
        xbody.addWidget(self.chk_comp_dep)
        return page_body(self._prefill_label(1), c, o, x)

    # ================================================================ step 3
    def _page_conditions(self) -> QWidget:
        c, body = card("Temperature and pressure")
        self.sp_T = _spin(300, 2000, 950, 1, 5, " °C")
        self.sp_T.valueChanged.connect(self._update_fo2_label)
        self.sp_T_sig = _spin(0, 300, 20, 1, 1, " K")
        body.addWidget(pair(field("Temperature", self.sp_T), field("± 1σ", self.sp_T_sig)))
        self.sp_P = _spin(0, 5000, 200, 1, 10, " MPa")
        self.sp_P.valueChanged.connect(self._update_fo2_label)
        self.sp_P_sig = _spin(0, 2000, 100, 1, 10, " MPa")
        body.addWidget(pair(field("Pressure", self.sp_P), field("± 1σ", self.sp_P_sig)))

        f, fbody = card("Oxygen fugacity")
        self.cmb_fo2_mode = QComboBox()
        self.cmb_fo2_mode.addItems(["Relative to a buffer", "Absolute log fO2 (bar)"])
        self.cmb_fo2_mode.currentIndexChanged.connect(self._update_fo2_label)
        self.cmb_buffer = QComboBox()
        self.cmb_buffer.addItems(available_buffers())
        self.cmb_buffer.setCurrentText("NNO")
        self.cmb_buffer.currentIndexChanged.connect(self._update_fo2_label)
        fbody.addWidget(pair(field("Given as", self.cmb_fo2_mode),
                             field("Buffer", self.cmb_buffer)))
        self.sp_dbuf = _spin(-14, 8, 1.0, 2, 0.1)
        self.sp_dbuf.valueChanged.connect(self._update_fo2_label)
        self.sp_dbuf_sig = _spin(0, 5, 0.3, 2, 0.1)
        fbody.addWidget(pair(field("Value", self.sp_dbuf), field("± 1σ", self.sp_dbuf_sig)))
        self.lbl_fo2 = note("", "Hint")
        fbody.addWidget(self.lbl_fo2)

        r, rbody = card("Analytical resolution")
        self.cmb_resolution = QComboBox()
        for label, *_ in RESOLUTION_PRESETS:
            self.cmb_resolution.addItem(label)
        self.cmb_resolution.currentIndexChanged.connect(self._on_resolution_changed)
        self.sp_width = _spin(0, 200, 5.0, 1, 0.5, " um")
        self.sp_width.valueChanged.connect(self._update_beam_sigma)
        self.sp_beam = _spin(0, 50, 0.0, 2, 0.1, " um")
        rbody.addWidget(field("How the profile was measured", self.cmb_resolution))
        self._width_field = field("Spot or slit width", self.sp_width)
        rbody.addWidget(pair(self._width_field, field("Beam σ", self.sp_beam)))
        self.lbl_resolution = note("", "Hint")
        rbody.addWidget(self.lbl_resolution)
        self.sp_xscale_sig = _spin(0, 0.5, 0.0, 3, 0.005)
        rbody.addWidget(field("Distance scale error (relative, 1σ)", self.sp_xscale_sig))
        return page_body(self._prefill_label(2), c, f, r)

    # ================================================================ step 4
    def _page_model(self) -> QWidget:
        g, gbody = card("Geometry")
        self.cmb_geom = QComboBox()
        self.cmb_geom.addItems(["plane", "cylinder", "sphere"])
        self.cmb_geom.currentIndexChanged.connect(self._update_geometry_note)
        gbody.addWidget(self.cmb_geom)
        self.lbl_geom = note("", "Hint")
        gbody.addWidget(self.lbl_geom)

        b, bbody = card("Boundaries")
        self.cmb_bcl = QComboBox(); self.cmb_bcl.addItems(["Fixed concentration", "Zero flux"])
        self.cmb_bcr = QComboBox(); self.cmb_bcr.addItems(["Fixed concentration", "Zero flux"])
        bbody.addWidget(pair(field("Left end", self.cmb_bcl), field("Right end", self.cmb_bcr)))
        bbody.addWidget(note("Fixed: held by the melt or a large reservoir. Zero flux: a "
                             "closed system or the crystal centre.", "Hint"))

        i, ibody = card("Initial profile")
        self.cmb_ic = QComboBox()
        self.cmb_ic.currentIndexChanged.connect(self._on_ic_changed)
        ibody.addWidget(self.cmb_ic)
        self.lbl_ic = note("", "Hint")
        ibody.addWidget(self.lbl_ic)
        self.sp_x0 = _spin(-1e5, 1e5, 0.0, 3, 1, " um")
        self.sp_cl = _spin(-1e6, 1e6, 0.3, 5, 0.01)
        self.sp_cr = _spin(-1e6, 1e6, 0.18, 5, 0.01)
        self.sp_smooth = _spin(0, 50, 0.0, 2, 0.1, " um")
        self._step_controls = QWidget()
        sc = QVBoxLayout(self._step_controls)
        sc.setContentsMargins(0, 0, 0, 0)
        sc.setSpacing(10)
        sc.addWidget(pair(field("Step position", self.sp_x0),
                          field("Initial smoothing", self.sp_smooth)))
        sc.addWidget(pair(field("Left plateau", self.sp_cl), field("Right plateau", self.sp_cr)))
        bg = QPushButton("Guess from the data")
        bg.clicked.connect(self._guess_initial)
        sc.addWidget(row(bg))
        ibody.addWidget(self._step_controls)
        self.chk_free_x0 = QCheckBox("Fit the step position")
        self.chk_free_x0.setChecked(True)
        self.chk_free_plateaus = QCheckBox("Fit the plateaus")
        ibody.addWidget(pair(self.chk_free_x0, self.chk_free_plateaus))

        n, nbody = card("Solver and thermal history")
        self.sp_nodes = QSpinBox(); self.sp_nodes.setRange(51, 4001); self.sp_nodes.setValue(301)
        self.chk_force_num = QCheckBox("Always use the numerical solver")
        nbody.addWidget(pair(field("Grid nodes", self.sp_nodes), self.chk_force_num))
        self.chk_cooling = QCheckBox("Linear cooling")
        self.sp_Tend = _spin(300, 2000, 900, 1, 5, " °C")
        self.sp_Tend.setEnabled(False)
        self.chk_cooling.toggled.connect(self.sp_Tend.setEnabled)
        nbody.addWidget(pair(self.chk_cooling, field("Final temperature", self.sp_Tend)))
        self._update_geometry_note()
        return page_body(self._prefill_label(3), g, b, i, n)

    # ================================================================ step 5
    def _page_coefficient(self) -> QWidget:
        c, body = card("Diffusion coefficient")
        body.addWidget(note("Tick one to fit. Tick several to compare them.", "Hint"))
        self.lst_coef = QListWidget()
        self.lst_coef.setMinimumHeight(230)
        self.lst_coef.currentItemChanged.connect(self._coefficient_selected)
        body.addWidget(self.lst_coef)
        self.lbl_coef = note("", "Hint")
        body.addWidget(self.lbl_coef)
        b = ghost_button("Equation, ranges and citation")
        b.clicked.connect(lambda: self.show_coefficient_info())
        body.addWidget(row(b))
        return page_body(self._prefill_label(4), c)

    # ================================================================ step 6
    def _page_uncertainty(self) -> QWidget:
        c, body = card("Monte Carlo")
        self.sp_draws = QSpinBox(); self.sp_draws.setRange(20, 100000); self.sp_draws.setValue(500)
        self.sp_seed = QSpinBox(); self.sp_seed.setRange(0, 2 ** 31 - 1); self.sp_seed.setValue(12345)
        body.addWidget(pair(field("Draws", self.sp_draws, "500 for a quick look, 1000 for a paper"),
                            field("Random seed", self.sp_seed)))

        s, sbody = card("Sample")
        self.chk_mc_T = QCheckBox("Temperature"); self.chk_mc_T.setChecked(True)
        self.chk_mc_f = QCheckBox("Oxygen fugacity"); self.chk_mc_f.setChecked(True)
        self.chk_mc_P = QCheckBox("Pressure"); self.chk_mc_P.setChecked(True)
        self.chk_mc_D = QCheckBox("Diffusion coefficient"); self.chk_mc_D.setChecked(True)
        self.chk_mc_n = QCheckBox("Measurement noise"); self.chk_mc_n.setChecked(True)
        self.chk_mc_x = QCheckBox("Distance scale")
        sbody.addWidget(pair(self.chk_mc_T, self.chk_mc_f))
        sbody.addWidget(pair(self.chk_mc_P, self.chk_mc_D))
        sbody.addWidget(pair(self.chk_mc_n, self.chk_mc_x))
        self.chk_contrib = QCheckBox("Rank what each source contributes (slower)")
        sbody.addWidget(self.chk_contrib)

        d, dbody = card("How the coefficient is sampled")
        self.cmb_dmode = QComboBox()
        for label, key in (("Best available", "auto"), ("Published covariance", "covariance"),
                           ("Scatter of log D at T", "logD_at_T"),
                           ("Each parameter independently", "independent")):
            self.cmb_dmode.addItem(label, key)
        self.cmb_dmode.currentIndexChanged.connect(self._on_dmode_changed)
        dbody.addWidget(self.cmb_dmode)
        self.lbl_dmode = note("", "Hint")
        dbody.addWidget(self.lbl_dmode)
        self._on_dmode_changed()
        return page_body(c, s, d)

    def _on_dmode_changed(self):
        texts = {
            "auto": "Uses the covariance when the paper gives one, otherwise the scatter of log D.",
            "covariance": "Draws the Arrhenius parameters together from their covariance.",
            "logD_at_T": "Draws log D at the working temperature from the reported scatter.",
            "independent": "Ignores the correlation between D0 and Q and overstates the "
                           "uncertainty. Use it only to reproduce older estimates.",
        }
        self.lbl_dmode.setText(texts[self.cmb_dmode.currentData()])

    # ================================================================ results
    def _page_results(self) -> QWidget:
        w = QWidget(); w.setObjectName("Page")
        v = QVBoxLayout(w)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(0)

        split = QSplitter(Qt.Horizontal)
        self.summary_inner = QWidget()
        self.summary_inner.setObjectName("Summary")
        self.summary_layout = QVBoxLayout(self.summary_inner)
        self.summary_layout.setContentsMargins(18, 16, 18, 16)
        self.summary_layout.setSpacing(4)
        side = scrollable(self.summary_inner)
        side.setMinimumWidth(250)
        side.setMaximumWidth(360)
        side.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)
        side.setStyleSheet(f"QScrollArea {{ background:{theme.SURFACE}; "
                           f"border-right:1px solid {theme.BORDER}; }}")
        split.addWidget(side)

        right = QWidget(); right.setObjectName("Page")
        rv = QVBoxLayout(right)
        rv.setContentsMargins(0, 0, 0, 0)
        rv.setSpacing(0)
        rv.addWidget(self.plot, 1)
        rv.addWidget(divider())

        bar = QWidget(); bar.setObjectName("Page")
        bh = QHBoxLayout(bar)
        bh.setContentsMargins(18, 10, 18, 10)
        bh.setSpacing(9)
        self.btn_fit = primary_button("Fit")
        self.btn_fit.clicked.connect(self.run_fit)
        self.btn_cmp = QPushButton("Compare coefficients")
        self.btn_cmp.clicked.connect(self.run_compare)
        self.btn_mc = QPushButton("Monte Carlo")
        self.btn_mc.clicked.connect(self.run_mc)
        self.btn_stop = QPushButton("Stop")
        self.btn_stop.clicked.connect(self.stop_work)
        self.btn_stop.setEnabled(False)
        for b in (self.btn_fit, self.btn_cmp, self.btn_mc, self.btn_stop):
            bh.addWidget(b)
        bh.addStretch(1)
        for text, fn in (("Methods", self.show_methods),
                         ("Histogram", self.show_histogram),
                         ("Log", self.show_log)):
            gb = ghost_button(text); gb.clicked.connect(fn); bh.addWidget(gb)
        self.btn_export = QPushButton("Export...")
        self.btn_export.clicked.connect(self.export_results)
        bh.addWidget(self.btn_export)
        rv.addWidget(bar)
        split.addWidget(right)
        split.setStretchFactor(0, 0)
        split.setStretchFactor(1, 1)
        split.setSizes([300, 1140])
        v.addWidget(split)
        return w

    # ---------------------------------------------------------------- summary
    def _summary_line(self, key: str, value: str):
        k = QLabel(key); k.setObjectName("SummaryKey")
        val = WrapLabel(value); val.setObjectName("SummaryVal")
        self.summary_layout.addWidget(k)
        self.summary_layout.addWidget(val)
        self.summary_layout.addSpacing(6)

    def _summary_head(self, text: str, step: int):
        w = QWidget()
        h = QHBoxLayout(w)
        h.setContentsMargins(0, 10, 0, 2)
        h.setSpacing(6)
        lab = QLabel(text.upper())
        lab.setObjectName("SummaryHead")
        h.addWidget(lab)
        h.addStretch(1)
        e = ghost_button("edit")
        e.setStyleSheet("font-size:11px; padding:1px 6px;")
        e.clicked.connect(lambda: self._go(step))
        h.addWidget(e)
        self.summary_layout.addWidget(w)

    def _rebuild_summary(self):
        # setParent(None) detaches at once. deleteLater alone leaves the old widgets
        # painted until the event loop catches up.
        while self.summary_layout.count():
            item = self.summary_layout.takeAt(0)
            w = item.widget()
            if w is not None:
                w.setParent(None)
                w.deleteLater()

        if self.fit_result is not None:
            head = QLabel("RESULT")
            head.setObjectName("SummaryHead")
            self.summary_layout.addWidget(head)
            res = WrapLabel("")
            txt = (f"<b style='font-size:19px;color:{theme.ACCENT}'>"
                   f"{human_time(self.fit_result.t_seconds)}</b>")
            if self.mc_result is not None:
                txt += (f"<br><span style='color:{theme.TEXT}'>68%: "
                        f"{human_time(self.mc_result.p16)} to "
                        f"{human_time(self.mc_result.p84)}</span>"
                        f"<br><span style='color:{theme.TEXT_MUTED};font-size:11.5px'>95%: "
                        f"{human_time(self.mc_result.p2_5)} to "
                        f"{human_time(self.mc_result.p97_5)}</span>")
            txt += (f"<br><span style='color:{theme.TEXT_FAINT};font-size:11px'>"
                    f"reduced chi2 {self.fit_result.stats.reduced_chi2:.2f} &middot; "
                    f"R2 {self.fit_result.stats.r_squared:.4f}</span>")
            res.setText(txt)
            self.summary_layout.addWidget(res)
            self.summary_layout.addSpacing(6)
            self.summary_layout.addWidget(divider())

        self._summary_head("Data", 0)
        if self.profile is not None:
            p = self.profile
            name = Path(p.source).name if p.source else "loaded profile"
            self._summary_line(name, f"{len(p)} points, {p.x.min():.1f} to {p.x.max():.1f} um")
            if self.dataset is not None:
                self._summary_line("provenance",
                                   "measured, published" if self.dataset.kind == "measured"
                                   else "synthetic")
        else:
            self._summary_line("", "no data loaded")

        self._summary_head("Mineral", 1)
        self._summary_line(get_mineral(self.cmb_mineral.currentData()).name,
                           f"{self.cmb_species.currentText()}, "
                           f"{self.cmb_axis.currentText().lower()}")

        self._summary_head("Conditions", 2)
        self._summary_line("temperature", f"{self.sp_T.value():.0f} °C ± {self.sp_T_sig.value():.0f}")
        self._summary_line("pressure", f"{self.sp_P.value():.0f} MPa ± {self.sp_P_sig.value():.0f}")
        if self.cmb_fo2_mode.currentIndex() == 0:
            self._summary_line("oxygen fugacity",
                               f"{self.cmb_buffer.currentText()} {self.sp_dbuf.value():+.2f} "
                               f"± {self.sp_dbuf_sig.value():.2f}")
        else:
            self._summary_line("oxygen fugacity", f"log fO2 {self.sp_dbuf.value():.2f} bar")
        if self.sp_beam.value() > 0:
            self._summary_line("beam σ", f"{self.sp_beam.value():.2f} um "
                               f"({self.cmb_resolution.currentText().lower()})")

        self._summary_head("Model", 3)
        self._summary_line("geometry", self.cmb_geom.currentText())
        if self._equilibrium_ic():
            self._summary_line("initial profile", "equilibrium with the anorthite zoning")
        else:
            self._summary_line("initial profile",
                               f"step at {self.sp_x0.value():.1f} um, "
                               f"{self.sp_cl.value():.4g} to {self.sp_cr.value():.4g}")
        self._summary_line("boundaries",
                           f"{self.cmb_bcl.currentText().lower()} / "
                           f"{self.cmb_bcr.currentText().lower()}")

        self._summary_head("Coefficient", 4)
        for k in self._checked_keys():
            c = get_coefficient(k)
            flag = ""
            if not c.verified:
                flag = "  [unverified]"
            if c.superseded_by or c.superseded_note:
                flag += "  [superseded]"
            self._summary_line("", c.label + flag)

        self._summary_head("Uncertainty", 5)
        b = self._budget()
        self._summary_line(f"{self.sp_draws.value()} draws, seed {self.sp_seed.value()}",
                           ", ".join(b.active_sources()) or "nothing sampled")

        self.summary_layout.addStretch(1)

    # ================================================================ navigation
    def _go(self, index: int):
        index = max(0, min(index, RESULTS))
        if index > 0 and self.profile is None:
            self.statusBar().showMessage("Load a profile first")
            index = 0
        self.step = index
        self.pages.setCurrentIndex(index)
        for i, lab in enumerate(self.step_labels):
            _repolish(lab, "StepDotActive" if i == index
                      else ("StepDotDone" if i < index else "StepDot"))
        self.btn_back.setVisible(index > 0)
        if index == RESULTS:
            self.btn_next.setVisible(False)
            self.lbl_footer.setText("Click edit in the summary to change a setting.")
            self._rebuild_summary()
        else:
            self.btn_next.setVisible(True)
            self.btn_next.setText("Review and fit" if index == RESULTS - 1 else "Continue")
            self.lbl_footer.setText(f"Step {index + 1} of {RESULTS}")
        self.btn_next.setEnabled(self.profile is not None or index == 0)

    def _next(self):
        if self.step == 0 and self.profile is None:
            QMessageBox.information(self, "No data", "Load a file or an example first.")
            return
        self._go(self.step + 1)

    # ================================================================ data
    def load_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Load profile", "",
            "Tables (*.csv *.txt *.tsv *.xlsx *.xls);;All files (*)")
        if path:
            self._load(Path(path), dataset=None)

    def load_example(self):
        item = self.lst_examples.currentItem()
        if item is None:
            QMessageBox.information(self, "Nothing selected", "Pick an example from the list.")
            return
        d = ds.get(item.data(Qt.UserRole))
        if not d.exists:
            QMessageBox.warning(self, "File missing", f"{d.filename} is not in the examples folder.")
            return
        self._load(d.path, dataset=d)

    def _load(self, path: Path, dataset: Optional[ds.ExampleDataset]):
        try:
            df = read_table(path)
            an = (None, True)
            if dataset is not None:
                spec = ProfileSpec(**dataset.spec)
            else:
                dlg = ColumnDialog(df, suggest_spec(df), self)
                if dlg.exec() != QDialog.Accepted:
                    return
                spec = dlg.spec()
                an = dlg.an_column()
            self.profile = build_profile(df, spec, source=str(path))
            self.dataset = dataset
            self.an_values = None
            for lab in self._prefill.values():
                lab.setVisible(False)
            p = self.profile
            if dataset is not None:
                self._apply_dataset_settings(dataset, df)
                msg = (f"Loaded {dataset.name.split(' (')[0]}, {len(p)} points. The example "
                       "also filled in the mineral, conditions, model and coefficient. "
                       "Each step shows what it set.")
            else:
                msg = f"Loaded {path.name}, {len(p)} points from {p.x.min():.1f} to {p.x.max():.1f} um."
                if an[0]:
                    vals = self.profile.column(an[0])
                    self.an_values = vals / 100.0 if an[1] else vals
                self._refresh_ic_options()
            self.lbl_data.setText(msg)
            self.lbl_data.setVisible(True)
            self._log(f"loaded {path}")
            self._log("  " + p.spec.describe())
            for n in p.notes:
                self._log(f"  note: {n}")
            self._guess_initial()
            self.plot.show_data(p.x, p.C, p.sigma, y_label=self._y_label())
            self.fit_result = self.mc_result = self.compare_results = None
            self.btn_next.setEnabled(True)
            self.statusBar().showMessage(f"Loaded {len(p)} points")
        except Exception:
            QMessageBox.critical(self, "Could not load the file", traceback.format_exc())

    def _apply_dataset_settings(self, d: ds.ExampleDataset, df):
        """Fill in the later steps from a bundled dataset and say so on each step."""
        s = d.settings
        keys = [self.cmb_mineral.itemData(i) for i in range(self.cmb_mineral.count())]
        if d.mineral in keys:
            self.cmb_mineral.setCurrentIndex(keys.index(d.mineral))
        self.cmb_species.setCurrentText(d.species)
        self.sp_T.setValue(s.get("T_C", 950.0))
        self.sp_T_sig.setValue(s.get("sigma_T_K", 20.0))
        self.sp_P.setValue(s.get("P_MPa", 200.0))
        self.sp_P_sig.setValue(s.get("sigma_P_MPa", 100.0))
        if "fo2_absolute" in s:
            self.cmb_fo2_mode.setCurrentIndex(1)
            self.sp_dbuf.setValue(s["fo2_absolute"])
            fo2 = f"log fO2 {s['fo2_absolute']:g} bar"
        else:
            self.cmb_fo2_mode.setCurrentIndex(0)
            self.cmb_buffer.setCurrentText(s.get("buffer", "NNO"))
            self.sp_dbuf.setValue(s.get("delta_buffer", 0.0))
            fo2 = f"{s.get('buffer', 'NNO')} {s.get('delta_buffer', 0.0):+g}"
        self.sp_dbuf_sig.setValue(s.get("sigma_delta", 0.3))
        axis = s.get("axis")
        self.cmb_axis.setCurrentIndex({"a": 1, "b": 2, "c": 3}.get(axis, 0))
        self.cmb_geom.setCurrentText(s.get("geometry", "plane"))
        self.chk_comp_dep.setChecked(bool(s.get("composition_dependent", False)))
        if "x_composition" in s:
            self.sp_xcomp.setValue(s["x_composition"])
        want = s.get("coefficient")
        self._refresh_coefficients()
        for i in range(self.lst_coef.count()):
            it = self.lst_coef.item(i)
            it.setCheckState(Qt.Checked if it.data(Qt.UserRole) == want else Qt.Unchecked)
            if it.data(Qt.UserRole) == want:
                self.lst_coef.setCurrentRow(i)
        an_col = s.get("an_column")
        if an_col and an_col in df.columns:
            vals = self.profile.column(an_col)
            if s.get("an_is_percent"):
                vals = vals / 100.0
            self.an_values = vals
            self.sp_xcomp.setValue(float(np.nanmean(vals)))
            self._log(f"  anorthite read from '{an_col}' "
                      f"(X_An {np.nanmin(vals):.2f} to {np.nanmax(vals):.2f})")
        self._refresh_ic_options()
        self._select_ic("equilibrium_plag" if s.get("initial_condition") == "equilibrium_plag"
                        else "step")

        name = d.name.split(" (")[0]
        head = f"Set by the {name} example: "
        self._show_prefill(1, head + f"{self.cmb_mineral.currentText()}, {d.species}, "
                           f"{self.cmb_axis.currentText().lower()}.")
        self._show_prefill(2, head + f"{s.get('T_C', 950):g} ± {s.get('sigma_T_K', 20):g} °C, "
                           f"{s.get('P_MPa', 200):g} ± {s.get('sigma_P_MPa', 100):g} MPa, {fo2}.")
        self._show_prefill(3, head + f"{s.get('geometry', 'plane')} geometry, "
                           f"{self.cmb_ic.currentText().lower()}.")
        if want:
            self._show_prefill(4, head + get_coefficient(want).label + ".")

    def _show_prefill(self, step: int, text: str):
        lab = self._prefill.get(step)
        if lab is not None:
            lab.setText(text)
            lab.setVisible(True)

    # ================================================================ reactions
    def _log(self, msg: str):
        self._log_lines.append(str(msg))
        self.statusBar().showMessage(str(msg).splitlines()[0][:150])

    def _on_mineral_changed(self):
        mineral = get_mineral(self.cmb_mineral.currentData())
        self.cmb_species.blockSignals(True)
        self.cmb_species.clear()
        self.cmb_species.addItems(mineral.species_keys())
        self.cmb_species.blockSignals(False)
        self.cmb_axis.setEnabled(not mineral.isotropic)
        cv = mineral.composition_variable
        self.lbl_mineral_note.setText(f"Modelled as {cv.label} = {cv.definition}." if cv else "")
        self._refresh_coefficients()
        self._refresh_ic_options()
        self._update_fo2_label()

    def _on_axis_changed(self):
        self._angle_row.setVisible(self.cmb_axis.currentIndex() == 4)

    def _refresh_coefficients(self):
        mineral = self.cmb_mineral.currentData()
        species = self.cmb_species.currentText()
        self.lst_coef.clear()
        for c in list_coefficients(mineral, species):
            tags = []
            if c.recommended:
                tags.append("recommended")
            if not c.verified:
                tags.append("unverified")
            if c.superseded_by or c.superseded_note:
                tags.append("superseded")
            label = c.label + (("\n" + "  ·  ".join(tags)) if tags else "")
            it = QListWidgetItem(label)
            it.setData(Qt.UserRole, c.key)
            it.setFlags(it.flags() | Qt.ItemIsUserCheckable)
            it.setCheckState(Qt.Checked if c.recommended else Qt.Unchecked)
            if c.superseded_by or c.superseded_note:
                it.setForeground(Qt.darkRed)
            elif not c.verified:
                it.setForeground(Qt.darkYellow)
            self.lst_coef.addItem(it)
        if self.lst_coef.count():
            self.lst_coef.setCurrentRow(0)

    def _coefficient_selected(self, item, _prev=None):
        if item is None:
            return
        c = get_coefficient(item.data(Qt.UserRole))
        parts = [c.equation_text]
        if c.superseded_by:
            parts.append(f"Superseded by {cite(c.superseded_by)}.")
        elif c.superseded_note:
            parts.append("Superseded. See the details.")
        if not c.verified:
            parts.append("Not checked against the original paper. See the details.")
        self.lbl_coef.setText("\n".join(parts))
        _repolish(self.lbl_coef, "Warn" if (c.superseded_by or c.superseded_note
                                            or not c.verified) else "Hint")

    # --- initial profile -----------------------------------------------------
    def _refresh_ic_options(self):
        current = self.cmb_ic.currentData()
        self.cmb_ic.blockSignals(True)
        self.cmb_ic.clear()
        self.cmb_ic.addItem("Sharp step between two plateaus", "step")
        if self.cmb_mineral.currentData() == "plagioclase" and self.an_values is not None:
            self.cmb_ic.addItem("Equilibrium with the anorthite zoning", "equilibrium_plag")
        self.cmb_ic.blockSignals(False)
        self._select_ic(current or "step")

    def _select_ic(self, key: str):
        i = self.cmb_ic.findData(key)
        self.cmb_ic.setCurrentIndex(i if i >= 0 else 0)
        self._on_ic_changed()

    def _equilibrium_ic(self) -> bool:
        return self.cmb_ic.currentData() == "equilibrium_plag" and self.an_values is not None

    def _on_ic_changed(self):
        equil = self.cmb_ic.currentData() == "equilibrium_plag"
        self._step_controls.setVisible(not equil)
        self.chk_free_x0.setEnabled(not equil)
        self.chk_free_plateaus.setEnabled(not equil)
        if equil:
            self.lbl_ic.setText("The trace element starts in equilibrium with the measured "
                                "anorthite profile, C = C0 exp(A X_An / RT) (Dohmen et al. "
                                "2017, eq. A13).")
        else:
            self.lbl_ic.setText("Crank (1975) eq. 2.14")

    def _update_geometry_note(self):
        self.lbl_geom.setText(Geometry(self.cmb_geom.currentText()).describe())

    # --- resolution --------------------------------------------------------------
    def _on_resolution_changed(self):
        label, kind, width, sigma, hint = RESOLUTION_PRESETS[self.cmb_resolution.currentIndex()]
        self._width_field.setVisible(kind is not None)
        custom = label == "Custom sigma"
        self.sp_beam.setReadOnly(not custom)
        if kind is not None:
            self.sp_width.blockSignals(True)
            self.sp_width.setValue(width)
            self.sp_width.blockSignals(False)
            rule = ("σ = diameter / 4 for an evenly lit round spot." if kind == "spot"
                    else "σ = width / √12 for an evenly lit slit.")
            self.lbl_resolution.setText(f"{hint} {rule}")
        else:
            self.lbl_resolution.setText(hint)
        if sigma is not None:
            self.sp_beam.setValue(sigma)
        self._update_beam_sigma()

    def _update_beam_sigma(self):
        kind = RESOLUTION_PRESETS[self.cmb_resolution.currentIndex()][1]
        if kind == "spot":
            self.sp_beam.setValue(self.sp_width.value() / 4.0)
        elif kind == "slit":
            self.sp_beam.setValue(self.sp_width.value() / np.sqrt(12.0))

    def _update_fo2_label(self):
        try:
            T = self.sp_T.value() + 273.15
            P = self.sp_P.value() * 1e6
            if self.cmb_fo2_mode.currentIndex() == 0:
                lf = log_fo2_from_delta(self.cmb_buffer.currentText(), self.sp_dbuf.value(), T, P)
                txt = f"log fO2 = {lf:.2f} bar ({lf + 5:.2f} Pa), buffer from Frost (1991)"
            else:
                lf = self.sp_dbuf.value()
                txt = f"log fO2 = {lf:.2f} bar ({lf + 5:.2f} Pa)"
            self.lbl_fo2.setText(txt)
            self.cmb_buffer.setEnabled(self.cmb_fo2_mode.currentIndex() == 0)
        except Exception as exc:
            self.lbl_fo2.setText(str(exc))

    def _guess_initial(self):
        if self.profile is None:
            return
        ic = guess_step_from_data(self.profile.x, self.profile.C)
        self.sp_x0.setValue(ic.params["x0"])
        self.sp_cl.setValue(ic.params["C_left"])
        self.sp_cr.setValue(ic.params["C_right"])
        if self.an_values is None:
            self.sp_xcomp.setValue(float(np.mean(self.profile.C)))
        self._log("initial profile: " + ic.describe())

    # ================================================================ model
    def _checked_keys(self) -> List[str]:
        keys = [self.lst_coef.item(i).data(Qt.UserRole)
                for i in range(self.lst_coef.count())
                if self.lst_coef.item(i).checkState() == Qt.Checked]
        if not keys and self.lst_coef.currentItem():
            keys = [self.lst_coef.currentItem().data(Qt.UserRole)]
        return keys

    def _conditions(self, coef) -> Conditions:
        T = self.sp_T.value() + 273.15
        P = self.sp_P.value() * 1e6
        lf = (log_fo2_from_delta(self.cmb_buffer.currentText(), self.sp_dbuf.value(), T, P)
              if self.cmb_fo2_mode.currentIndex() == 0 else self.sp_dbuf.value())
        mineral = get_mineral(self.cmb_mineral.currentData())
        X = {mineral.composition_variable.key: self.sp_xcomp.value()}
        for k in coef.requires:
            X.setdefault(k, self.sp_xcomp.value())
        axis = {1: "a", 2: "b", 3: "c"}.get(self.cmb_axis.currentIndex())
        angles = ((self.sp_alpha.value(), self.sp_beta.value(), self.sp_gamma.value())
                  if self.cmb_axis.currentIndex() == 4 else None)
        return Conditions(T_K=T, P_Pa=P, log_fo2_bar=lf, X=X, axis=axis, angles_deg=angles)

    def _model(self, coef_key: str) -> DiffusionModel:
        coef = get_coefficient(coef_key)
        cond = self._conditions(coef)
        mineral = get_mineral(self.cmb_mineral.currentData())
        species = self.cmb_species.currentText()
        if self._equilibrium_ic():
            ic = InitialCondition(
                "equilibrium_plag",
                {"x_an_x": np.asarray(self.profile.x, dtype=float),
                 "x_an_values": np.asarray(self.an_values, dtype=float),
                 "T_K": cond.T_K, "species": species,
                 "C_ref": float(self.sp_cl.value()),
                 "X_An_ref": float(self.an_values[0])},
                description=("equilibrium profile C0 exp(A X_An / RT) anchored at "
                             + format(self.sp_cl.value(), ".4g")
                             + " where X_An = " + format(float(self.an_values[0]), ".3f")
                             + " (Dohmen et al. 2017 App. eq. A13)"))
            c_left = float(self.sp_cl.value())
            c_right = float(self.sp_cr.value())
        else:
            params = {"x0": self.sp_x0.value(), "C_left": self.sp_cl.value(),
                      "C_right": self.sp_cr.value()}
            if self.sp_smooth.value() > 0:
                params["smooth"] = self.sp_smooth.value()
            ic = InitialCondition("step", params)
            c_left, c_right = params["C_left"], params["C_right"]
        bcl = dirichlet(c_left) if self.cmb_bcl.currentIndex() == 0 else zero_flux()
        bcr = dirichlet(c_right) if self.cmb_bcr.currentIndex() == 0 else zero_flux()
        hist = (ThermalHistory.linear(cond.T_K, self.sp_Tend.value() + 273.15, 1.0)
                if self.chk_cooling.isChecked() else None)
        comp_key = mineral.composition_variable.key if self.chk_comp_dep.isChecked() else None
        if comp_key and comp_key not in coef.requires:
            comp_key = None

        theta = 0.0
        an_grid = None
        if mineral.key == "plagioclase" and species in ACTIVITY_A and self.an_values is not None:
            theta = activity_theta(species, cond.T_K)
        n_nodes = self.sp_nodes.value()
        x_grid = None
        if self.profile is not None:
            x_grid = np.linspace(float(self.profile.x.min()), float(self.profile.x.max()), n_nodes)
            if theta and self.an_values is not None:
                an_grid = np.interp(x_grid, self.profile.x, self.an_values)
        return DiffusionModel(
            coefficient=coef, conditions=cond, initial=ic,
            geometry=Geometry(self.cmb_geom.currentText()),
            bc_left=bcl, bc_right=bcr, history=hist,
            beam_sigma_um=self.sp_beam.value(), n_nodes=n_nodes, x_grid=x_grid,
            comp_key=comp_key, composition_dependent=bool(comp_key),
            an_profile=an_grid, activity_theta=theta if an_grid is not None else 0.0,
            force_numerical=self.chk_force_num.isChecked())

    def _free_parameters(self):
        free = ["t"]
        if self._equilibrium_ic():
            return free
        if self.chk_free_x0.isChecked():
            free.append("x0")
        if self.chk_free_plateaus.isChecked():
            free += ["C_left", "C_right"]
        return free

    def _budget(self) -> UncertaintyBudget:
        return UncertaintyBudget(
            sigma_T_K=self.sp_T_sig.value() if self.chk_mc_T.isChecked() else 0.0,
            fo2_mode="buffer" if self.cmb_fo2_mode.currentIndex() == 0 else "absolute",
            buffer=self.cmb_buffer.currentText(),
            delta_buffer=self.sp_dbuf.value(),
            sigma_delta_buffer=self.sp_dbuf_sig.value() if self.chk_mc_f.isChecked() else 0.0,
            sigma_log_fo2=self.sp_dbuf_sig.value() if self.chk_mc_f.isChecked() else 0.0,
            sigma_P_Pa=self.sp_P_sig.value() * 1e6 if self.chk_mc_P.isChecked() else 0.0,
            sample_coefficient=self.chk_mc_D.isChecked(),
            coefficient_mode=self.cmb_dmode.currentData(),
            sample_measurement_noise=self.chk_mc_n.isChecked(),
            sigma_distance_scale=self.sp_xscale_sig.value() if self.chk_mc_x.isChecked() else 0.0)

    def _y_label(self) -> str:
        mineral = get_mineral(self.cmb_mineral.currentData())
        if self.profile is not None and self.profile.spec and self.profile.spec.mode == "A":
            return str(self.profile.spec.column_a)
        return f"{mineral.composition_variable.label} ({self.cmb_species.currentText()})"

    # ================================================================ running
    # Worker signals are connected only to methods of this window. Qt then delivers
    # them on the interface thread. A lambda has no thread of its own, so it would
    # run on the worker thread and touch widgets from there, which crashes.
    def _running(self) -> bool:
        return any(t.isRunning() for t, _ in self._jobs)

    def _launch(self, worker, busy_max: int = 0) -> bool:
        if self._running():
            self.statusBar().showMessage("Still working. Wait or press Stop.")
            return False
        self._busy(True, busy_max)
        thread = start(worker)
        self._jobs.append((thread, worker))
        thread.finished.connect(self._reap_jobs)
        return True

    @Slot()
    def _reap_jobs(self):
        self._jobs = [(t, w) for t, w in self._jobs if not t.isFinished()]

    def _busy(self, on: bool, maximum: int = 0):
        for b in (self.btn_fit, self.btn_cmp, self.btn_mc):
            b.setEnabled(not on)
        self.btn_stop.setEnabled(on)
        self.progress.setVisible(on)
        self.progress.setMaximum(maximum)
        self.progress.setValue(0)

    def stop_work(self):
        for _, w in self._jobs:
            if hasattr(w, "abort"):
                w.abort()
        self._log("stop requested")

    def closeEvent(self, event):
        self.stop_work()
        for t, _ in self._jobs:
            t.quit()
            t.wait(3000)
        super().closeEvent(event)

    def _ready(self) -> bool:
        if self.profile is None:
            QMessageBox.information(self, "No data", "Load a profile first.")
            return False
        if not self._checked_keys():
            QMessageBox.information(self, "No coefficient", "Choose a coefficient on step 5.")
            return False
        return True

    def run_fit(self):
        if not self._ready():
            return
        try:
            model = self._model(self._checked_keys()[0])
        except Exception:
            QMessageBox.critical(self, "Model error", traceback.format_exc())
            return
        p = self.profile
        w = FitWorker(model, p.x, p.C, p.sigma, self._free_parameters(), 1e2, 3.2e12)
        w.finished.connect(self._fit_done)
        w.failed.connect(self._work_failed)
        self._launch(w)

    @Slot(object)
    def _fit_done(self, res):
        self._busy(False)
        self.fit_result = res
        self.mc_result = None
        self.compare_results = None
        self.plot.show_fit(res, None, y_label=self._y_label())
        self._log("fit complete: " + human_time(res.t_seconds))
        for wmsg in res.warnings:
            self._log("  ! " + wmsg)
        self._rebuild_summary()
        self._show_warnings(res.warnings)

    def run_compare(self):
        if not self._ready():
            return
        keys = self._checked_keys()
        if len(keys) < 2:
            QMessageBox.information(self, "Tick at least two",
                                    "Tick two or more coefficients on step 5.")
            return
        try:
            models = {get_coefficient(k).label: self._model(k) for k in keys}
        except Exception:
            QMessageBox.critical(self, "Model error", traceback.format_exc())
            return
        p = self.profile
        w = CompareWorker(models, p.x, p.C, p.sigma, self._free_parameters(), 1e2, 3.2e12)
        w.progress.connect(self._compare_progress)
        w.finished.connect(self._compare_done)
        w.failed.connect(self._work_failed)
        self._launch(w, len(models))

    @Slot(int, int, str)
    def _compare_progress(self, i, n, key):
        self.progress.setValue(i)
        self.statusBar().showMessage(f"{i}/{n}: {key}")

    @Slot(object)
    def _compare_done(self, results):
        self._busy(False)
        self.compare_results = results
        ok = {k: r for k, r in results.items() if not isinstance(r, Exception)}
        lines = []
        for k, r in results.items():
            if isinstance(r, Exception):
                lines.append(f"{k}\n   failed: {r}")
                continue
            lines.append(f"{k}\n   {human_time(r.t_seconds):<14s} "
                         f"reduced chi2 {r.stats.reduced_chi2:.3g}   R2 {r.stats.r_squared:.4f}")
            c = r.model.coefficient
            if c.superseded_by:
                lines.append(f"   superseded by {cite(c.superseded_by)}")
            elif not c.verified:
                lines.append("   not checked against the original paper")
        if len(ok) > 1:
            ts = [r.t_seconds for r in ok.values()]
            lines += ["", f"Range {human_time(min(ts))} to {human_time(max(ts))}, "
                          f"a factor of {max(ts) / min(ts):.1f}."]
        self.plot.show_comparison(ok, y_label=self._y_label())
        if ok:
            self.fit_result = next(iter(ok.values()))
        self._rebuild_summary()
        self._log("comparison complete")
        _text_dialog(self, "Coefficient comparison", "\n".join(lines), 820, 460)

    def run_mc(self):
        if not self._ready():
            return
        try:
            model = self._model(self._checked_keys()[0])
        except Exception:
            QMessageBox.critical(self, "Model error", traceback.format_exc())
            return
        p = self.profile
        n = self.sp_draws.value()
        w = MonteCarloWorker(model, p.x, p.C, p.sigma, self._budget(), n,
                             self.sp_seed.value(), self._free_parameters(), 1e2, 3.2e12,
                             do_contributions=self.chk_contrib.isChecked(),
                             contribution_draws=max(40, n // 5))
        w.progress.connect(self._mc_progress)
        w.stage.connect(self._mc_stage)
        w.finished.connect(self._mc_done)
        w.failed.connect(self._work_failed)
        if self._launch(w, n):
            self._log(f"Monte Carlo: {n} draws, seed {self.sp_seed.value()}, "
                      f"sampling {', '.join(self._budget().active_sources())}")

    @Slot(int, int)
    def _mc_progress(self, i, total):
        self.progress.setValue(i)

    @Slot(str)
    def _mc_stage(self, stage):
        self.statusBar().showMessage(stage)

    @Slot(object)
    def _mc_done(self, res):
        self._busy(False)
        self.mc_result = res
        if self.fit_result is None:
            from ..fitting import fit_time
            self.fit_result = fit_time(self._model(self._checked_keys()[0]), self.profile.x,
                                       self.profile.C, self.profile.sigma,
                                       self._free_parameters())
        self.plot.show_fit(self.fit_result, res, y_label=self._y_label())
        self._rebuild_summary()
        self._log("Monte Carlo complete: median " + human_time(res.median))
        _text_dialog(self, "Monte Carlo result", res.summary(), 780, 460)

    @Slot(str)
    def _work_failed(self, msg: str):
        self._busy(False)
        self._log(msg)
        QMessageBox.critical(self, "Calculation failed", msg)

    def _show_warnings(self, warnings):
        serious = [w for w in warnings
                   if "SUPERSEDED" in w or "NOT verified" in w or "semi-infinite" in w
                   or "resolution limit" in w]
        if serious:
            QMessageBox.warning(self, "Check before using the result", "\n\n".join(serious))

    # ================================================================ dialogs
    def show_methods(self):
        if self.fit_result is None:
            QMessageBox.information(self, "Nothing yet", "Run a fit first.")
            return
        _text_dialog(self, "Methods and references",
                     methods_paragraph(self.fit_result, self.mc_result, self.profile))

    def show_log(self):
        _text_dialog(self, "Log", "\n".join(self._log_lines) or "(nothing yet)", 820, 500)

    def show_histogram(self):
        if self.mc_result is None:
            QMessageBox.information(self, "No Monte Carlo yet", "Run a Monte Carlo first.")
            return
        self.plot.show_histogram(self.mc_result)

    def show_coefficient_info(self):
        it = self.lst_coef.currentItem()
        if it is None:
            return
        c = get_coefficient(it.data(Qt.UserRole))
        text = c.describe() + "\n\nFull reference:\n  " + format_reference(c.citation)
        for k in c.secondary_citations:
            text += "\n  see also: " + format_reference(k)
        _text_dialog(self, c.label, text)

    def export_results(self):
        if self.fit_result is None:
            QMessageBox.information(self, "Nothing to export", "Run a fit first.")
            return
        d = QFileDialog.getExistingDirectory(self, "Export into folder")
        if not d:
            return
        try:
            written = save_results(d, self.fit_result, self.mc_result, self.profile,
                                   figure=self.plot.figure)
            self._log("exported:\n  " + "\n  ".join(f"{k}: {v}" for k, v in written.items()))
            QMessageBox.information(self, "Exported", "Written:\n" + "\n".join(written.values()))
        except Exception:
            QMessageBox.critical(self, "Export failed", traceback.format_exc())

    def show_about(self):
        from .. import __version__
        QMessageBox.about(self, "About Diffusor", (
            f"<b>Diffusor {__version__}</b><br><br>"
            "Diffusion chronometry with closed-form and Crank-Nicolson solvers, a registry "
            "of published diffusion coefficients, and Monte Carlo error propagation.<br><br>"
            "Every equation and coefficient carries its citation. View > Methods lists the "
            "references used by the current run."))
