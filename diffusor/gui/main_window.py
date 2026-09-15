"""Diffusor's main window: a stepped flow ending in a results view.

The settings are collected one group at a time, then summarised in a narrow
sidebar next to the plot. Any group in that summary can be clicked to go back
and change it.
"""
from __future__ import annotations

import copy
import traceback
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtGui import QAction, QFont
from PySide6.QtWidgets import (QButtonGroup, QCheckBox, QComboBox, QDialog,
                               QDialogButtonBox, QDoubleSpinBox, QFileDialog, QFrame,
                               QHBoxLayout, QLabel, QLineEdit, QListWidget,
                               QListWidgetItem, QMainWindow, QMessageBox, QProgressBar,
                               QPushButton, QRadioButton, QSizePolicy, QSpinBox,
                               QSplitter, QStackedWidget, QTextEdit, QVBoxLayout,
                               QWidget)

from .. import datasets as ds
from ..coefficients import Conditions, get as get_coefficient, list_coefficients
from ..coefficients.plagioclase import ACTIVITY_A, activity_theta
from ..dataio import (ProfileSpec, build_profile, methods_paragraph, read_table,
                      save_results, suggest_spec)
from ..fitting import DiffusionModel, UncertaintyBudget
from ..minerals import MINERALS, get_mineral
from ..references import format_reference
from ..solvers import Geometry, InitialCondition, dirichlet, zero_flux
from ..solvers.history import ThermalHistory
from ..solvers.initial import guess_step_from_data
from ..thermo import available_buffers, log_fo2_from_delta
from ..thermo.units import human_time
from . import format_help, theme
from .plot_widget import ProfilePlot
from .widgets import (WrapLabel, badge, card, divider, field, ghost_button, note,
                      pair, page_body, primary_button, row, scrollable)
from .workers import CompareWorker, FitWorker, MonteCarloWorker, start

STEPS = ["Data", "Mineral", "Conditions", "Model", "Coefficient", "Uncertainty", "Results"]
RESULTS = len(STEPS) - 1


def _spin(lo, hi, val, dec=3, step=1.0, suffix=""):
    s = QDoubleSpinBox()
    s.setRange(lo, hi)
    s.setDecimals(dec)
    s.setSingleStep(step)
    s.setValue(val)
    if suffix:
        s.setSuffix(suffix)
    return s


def _text_dialog(parent, title, text, width=940, height=640):
    dlg = QDialog(parent)
    dlg.setWindowTitle(title)
    dlg.resize(width, height)
    lay = QVBoxLayout(dlg)
    lay.setContentsMargins(18, 18, 18, 18)
    te = QTextEdit()
    te.setReadOnly(True)
    te.setPlainText(text)
    te.setFont(QFont(theme.MONO_STACK.split(",")[0].strip('"'), 9))
    lay.addWidget(te)
    bb = QDialogButtonBox(QDialogButtonBox.Close)
    bb.rejected.connect(dlg.reject)
    bb.accepted.connect(dlg.accept)
    lay.addWidget(bb)
    dlg.exec()


class ColumnDialog(QDialog):
    """Map the columns of a loaded table onto the model."""

    def __init__(self, df, spec: ProfileSpec, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Map your columns")
        self.resize(560, 640)
        cols = [""] + [str(c) for c in df.columns]
        root = QVBoxLayout(self)
        root.setContentsMargins(20, 18, 20, 18)
        root.setSpacing(12)
        head = QLabel("Map your columns")
        head.setObjectName("H2")
        root.addWidget(head)
        root.addWidget(note("Diffusor has guessed from the column names. Check each one.", "Sub"))

        self.dist = QComboBox(); self.dist.addItems(cols)
        self.unit = QComboBox(); self.unit.addItems(["um", "mm", "nm", "m"])
        self.a = QComboBox(); self.a.addItems(cols)
        self.b = QComboBox(); self.b.addItems(cols)
        self.sa = QComboBox(); self.sa.addItems(cols)
        self.sb = QComboBox(); self.sb.addItems(cols)
        self.mode = QComboBox(); self.mode.addItems(["A/(A+B)", "B/(A+B)", "A", "A-B"])
        self.oxa = QLineEdit(); self.oxb = QLineEdit()
        self.oxa.setPlaceholderText("e.g. FeO; blank if not wt% oxide")
        self.oxb.setPlaceholderText("e.g. MgO")
        self.slevel = QComboBox(); self.slevel.addItems(["1s", "2s"])

        def setc(combo, val):
            combo.setCurrentText(str(val) if val is not None else "")

        setc(self.dist, spec.distance_column); setc(self.a, spec.column_a)
        setc(self.b, spec.column_b); setc(self.sa, spec.sigma_a_column)
        setc(self.sb, spec.sigma_b_column)
        self.mode.setCurrentText(spec.mode)

        root.addWidget(pair(field("Distance column", self.dist),
                            field("Distance unit", self.unit)))
        root.addWidget(pair(field("Element A", self.a), field("Element B (optional)", self.b)))
        root.addWidget(pair(field("Uncertainty of A", self.sa),
                            field("Uncertainty of B", self.sb)))
        root.addWidget(pair(field("Uncertainty level", self.slevel),
                            field("Modelled variable", self.mode)))
        root.addWidget(pair(field("Oxide of A", self.oxa), field("Oxide of B", self.oxb)))
        root.addWidget(note(
            "Name the oxides only if columns A and B are weight per cent oxide. Diffusor then "
            "converts to cation moles before taking the ratio.", "Hint"))
        root.addStretch(1)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.button(QDialogButtonBox.Ok).setObjectName("Primary")
        bb.accepted.connect(self.accept); bb.rejected.connect(self.reject)
        root.addWidget(bb)

    def spec(self) -> ProfileSpec:
        def g(c):
            return c.currentText().strip() or None
        return ProfileSpec(
            distance_column=g(self.dist), column_a=g(self.a), column_b=g(self.b),
            sigma_a_column=g(self.sa), sigma_b_column=g(self.sb),
            distance_unit=self.unit.currentText(), mode=self.mode.currentText(),
            oxide_a=self.oxa.text().strip() or None, oxide_b=self.oxb.text().strip() or None,
            sigma_level=self.slevel.currentText())


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Diffusor")
        self.resize(1440, 920)
        self.setStyleSheet(theme.stylesheet())

        self.profile = None
        self.dataset: Optional[ds.ExampleDataset] = None
        self.an_values = None
        self.fit_result = None
        self.mc_result = None
        self.compare_results = None
        self._threads: List = []
        self._workers: List = []
        self._log_lines: List[str] = []
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
        self._on_ic_changed()
        self._go(0)

    # ================================================================ chrome
    def _build_header(self) -> QWidget:
        w = QWidget(); w.setObjectName("Page")
        h = QHBoxLayout(w)
        h.setContentsMargins(24, 14, 24, 14)
        h.setSpacing(16)
        title = QLabel("Diffusor")
        title.setObjectName("H1")
        h.addWidget(title)
        sub = QLabel("diffusion chronometry")
        sub.setObjectName("Sub")
        h.addWidget(sub)
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
        h.setContentsMargins(24, 12, 24, 12)
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
        a = QAction("Input file format", self,
                    triggered=lambda: _text_dialog(self, "Input file format",
                                                   format_help.as_plain_text()))
        h.addAction(a)
        a = QAction("About", self); a.triggered.connect(self.show_about); h.addAction(a)

    # ================================================================ step 1
    def _page_data(self) -> QWidget:
        load_card, body = card(
            "Load a profile",
            "A measured traverse across the zone boundary you want to date.")
        btn = primary_button("Choose a file...")
        btn.clicked.connect(self.load_file)
        body.addWidget(row(btn, note(f"Accepts {format_help.FILE_TYPES}", "Hint"),
                           stretch_last=True))
        self.lbl_data = note("No file loaded yet.", "Sub")
        body.addWidget(self.lbl_data)

        fmt_card, fbody = card("How the file has to look", format_help.SHORT)
        sample = QLabel(format_help.EXAMPLE_TABLE)
        sample.setObjectName("Mono")
        sample.setStyleSheet(
            f"background:{theme.SURFACE_ALT}; border:1px solid {theme.BORDER};"
            f"border-radius:7px; padding:11px; font-family:{theme.MONO_STACK};")
        fbody.addWidget(sample)
        self._rules_box = QWidget()
        rb = QVBoxLayout(self._rules_box)
        rb.setContentsMargins(0, 4, 0, 0)
        rb.setSpacing(9)
        for title, text in format_help.RULES:
            t = QLabel(title); t.setObjectName("H2")
            t.setStyleSheet("font-size:12.5px; font-weight:600;")
            rb.addWidget(t)
            rb.addWidget(note(text, "Hint"))
        rb.addWidget(note("Greyscale profiles: " + format_help.GREYSCALE, "Hint"))
        self._rules_box.setVisible(False)
        toggle = ghost_button("Show the detailed rules")

        def _toggle():
            vis = not self._rules_box.isVisible()
            self._rules_box.setVisible(vis)
            toggle.setText("Hide the detailed rules" if vis else "Show the detailed rules")
        toggle.clicked.connect(_toggle)
        fbody.addWidget(row(toggle, stretch_last=True))
        fbody.addWidget(self._rules_box)

        ex_card, ebody = card("Or start with an example", format_help.NO_DATA)
        self.lst_examples = QListWidget()
        self.lst_examples.setMinimumHeight(190)
        for d in ds.DATASETS:
            it = QListWidgetItem(f"{d.name}\n{d.mineral} / {d.species}")
            it.setData(Qt.UserRole, d.key)
            if not d.exists:
                it.setFlags(it.flags() & ~Qt.ItemIsEnabled)
                it.setText(it.text() + "   (file missing)")
            self.lst_examples.addItem(it)
        self.lst_examples.currentItemChanged.connect(self._example_selected)
        ebody.addWidget(self.lst_examples)
        self.lbl_example = note("Select one to see where its numbers come from.", "Hint")
        ebody.addWidget(self.lbl_example)
        b = QPushButton("Load the selected example")
        b.clicked.connect(self.load_example)
        ebody.addWidget(row(b, stretch_last=True))
        return page_body(load_card, fmt_card, ex_card)

    def _example_selected(self, item, _prev=None):
        if item is None:
            return
        d = ds.get(item.data(Qt.UserRole))
        self.lbl_example.setText(d.provenance_banner() +
                                 (f"\n\nExpected answer: {d.expected}" if d.expected else ""))
        self.lbl_example.setObjectName("Warn" if d.kind == "synthetic" else "Hint")
        self.lbl_example.setStyleSheet("")
        self.lbl_example.style().unpolish(self.lbl_example)
        self.lbl_example.style().polish(self.lbl_example)

    # ================================================================ step 2
    def _page_mineral(self) -> QWidget:
        c, body = card("Mineral and diffusing species",
                       "This decides which published diffusion coefficients apply.")
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

        o, obody = card("Crystal orientation",
                        "Diffusion is anisotropic in most of these minerals, so the direction "
                        "of the traverse matters.")
        self.cmb_axis = QComboBox()
        self.cmb_axis.addItems(["Use the reference axis of the publication",
                                "a-axis", "b-axis", "c-axis", "Angles to a, b and c"])
        self.cmb_axis.currentIndexChanged.connect(self._on_axis_changed)
        obody.addWidget(field("Traverse direction", self.cmb_axis))
        self.sp_alpha = _spin(0, 180, 90, 1, 1, " deg")
        self.sp_beta = _spin(0, 180, 90, 1, 1, " deg")
        self.sp_gamma = _spin(0, 180, 0, 1, 1, " deg")
        self._angle_row = row(field("to a", self.sp_alpha), field("to b", self.sp_beta),
                              field("to c", self.sp_gamma), stretch_last=True)
        self._angle_row.setVisible(False)
        obody.addWidget(self._angle_row)
        obody.addWidget(note(
            "With angles, Diffusor applies D = D_a cos2(alpha) + D_b cos2(beta) + D_c cos2(gamma) "
            "after Costa & Chakraborty (2004). The three cosines squared must sum to one.",
            "Hint"))

        x, xbody = card("Composition",
                        "The diffusion coefficient depends on composition in most of these "
                        "minerals.")
        self.sp_xcomp = _spin(0, 1, 0.15, 4, 0.01)
        xbody.addWidget(field("Representative composition", self.sp_xcomp))
        self.chk_comp_dep = QCheckBox("Let the fitted profile set the local composition")
        self.chk_comp_dep.setChecked(True)
        xbody.addWidget(self.chk_comp_dep)
        xbody.addWidget(note(
            "With this on, D is re-evaluated at every grid node from the profile itself, which "
            "is the physically correct treatment and forces the numerical solver. With it off, "
            "D is constant at the value above and a closed-form solution can be used.", "Hint"))
        return page_body(c, o, x)

    # ================================================================ step 3
    def _page_conditions(self) -> QWidget:
        c, body = card("Temperature and pressure",
                       "Give the uncertainty on each: it is propagated by the Monte Carlo, and "
                       "temperature is almost always the dominant term.")
        self.sp_T = _spin(300, 2000, 950, 1, 5, " C")
        self.sp_T.valueChanged.connect(self._update_fo2_label)
        self.sp_T_sig = _spin(0, 300, 20, 1, 1, " K")
        body.addWidget(pair(field("Temperature", self.sp_T),
                            field("1 sigma", self.sp_T_sig)))
        self.sp_P = _spin(0, 5000, 200, 1, 10, " MPa")
        self.sp_P.valueChanged.connect(self._update_fo2_label)
        self.sp_P_sig = _spin(0, 2000, 100, 1, 10, " MPa")
        body.addWidget(pair(field("Pressure", self.sp_P), field("1 sigma", self.sp_P_sig)))

        f, fbody = card("Oxygen fugacity",
                        "Given as an offset from a mineral buffer, the buffer is re-evaluated "
                        "at every sampled temperature, so fO2 and T stay correlated.")
        self.cmb_fo2_mode = QComboBox()
        self.cmb_fo2_mode.addItems(["Offset from a buffer", "Absolute log10 fO2 (bar)"])
        self.cmb_fo2_mode.currentIndexChanged.connect(self._update_fo2_label)
        self.cmb_buffer = QComboBox()
        self.cmb_buffer.addItems(available_buffers())
        self.cmb_buffer.setCurrentText("NNO")
        self.cmb_buffer.currentIndexChanged.connect(self._update_fo2_label)
        fbody.addWidget(pair(field("Specified as", self.cmb_fo2_mode),
                             field("Buffer", self.cmb_buffer)))
        self.sp_dbuf = _spin(-14, 8, 1.0, 2, 0.1)
        self.sp_dbuf.valueChanged.connect(self._update_fo2_label)
        self.sp_dbuf_sig = _spin(0, 5, 0.3, 2, 0.1)
        fbody.addWidget(pair(field("Offset / value", self.sp_dbuf),
                             field("1 sigma", self.sp_dbuf_sig)))
        self.lbl_fo2 = note("", "Good")
        fbody.addWidget(self.lbl_fo2)

        r, rbody = card("Analytical resolution",
                        "A microbeam averages over its interaction volume, which makes a short "
                        "profile look more diffused than it is.")
        self.sp_beam = _spin(0, 20, 0.0, 2, 0.1, " um")
        self.sp_xscale_sig = _spin(0, 0.5, 0.0, 3, 0.005)
        rbody.addWidget(pair(
            field("Beam sigma", self.sp_beam, "0 disables the convolution correction"),
            field("Distance scale 1 sigma", self.sp_xscale_sig,
                  "relative, e.g. 0.02 for a 2% image calibration")))
        return page_body(c, f, r)

    # ================================================================ step 4
    def _page_model(self) -> QWidget:
        g, gbody = card("Geometry",
                        "A petrological choice, not a numerical one. Modelling a 3-D crystal "
                        "in 1-D always returns a maximum estimate of the time.")
        self.cmb_geom = QComboBox()
        self.cmb_geom.addItems(["plane", "cylinder", "sphere"])
        self.cmb_geom.currentIndexChanged.connect(self._update_geometry_note)
        gbody.addWidget(field("Shape", self.cmb_geom))
        self.lbl_geom = note("", "Hint")
        gbody.addWidget(self.lbl_geom)

        b, bbody = card("Boundaries",
                        "Fixed concentration means an open system held by an infinite melt "
                        "reservoir; zero flux means a closed system or a symmetry plane.")
        self.cmb_bcl = QComboBox(); self.cmb_bcl.addItems(["Fixed concentration", "Zero flux"])
        self.cmb_bcr = QComboBox(); self.cmb_bcr.addItems(["Fixed concentration", "Zero flux"])
        bbody.addWidget(pair(field("Left end of the traverse", self.cmb_bcl),
                             field("Right end", self.cmb_bcr)))

        i, ibody = card("Initial condition",
                        "The largest single source of systematic error. A profile that grew "
                        "partly by crystal growth will give a spuriously long time if all of "
                        "it is attributed to diffusion.")
        self.cmb_ic = QComboBox()
        self.cmb_ic.addItems(["Step between two plateaus",
                              "Equilibrium profile from the anorthite gradient"])
        self.cmb_ic.currentIndexChanged.connect(self._on_ic_changed)
        ibody.addWidget(field("Shape of the initial profile", self.cmb_ic))
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
        sc.addWidget(pair(field("Interface position x0", self.sp_x0),
                          field("Initial smoothing", self.sp_smooth)))
        sc.addWidget(pair(field("Plateau on the left", self.sp_cl),
                          field("Plateau on the right", self.sp_cr)))
        bg = QPushButton("Guess from the data")
        bg.clicked.connect(self._guess_initial)
        sc.addWidget(row(bg, stretch_last=True))
        ibody.addWidget(self._step_controls)
        self.chk_free_x0 = QCheckBox("Let the fit adjust x0")
        self.chk_free_x0.setChecked(True)
        self.chk_free_plateaus = QCheckBox("Let the fit adjust the plateaus")
        ibody.addWidget(self.chk_free_x0)
        ibody.addWidget(self.chk_free_plateaus)

        n, nbody = card("Solver and thermal history")
        self.sp_nodes = QSpinBox(); self.sp_nodes.setRange(51, 4001); self.sp_nodes.setValue(301)
        self.chk_force_num = QCheckBox("Always use the numerical solver")
        nbody.addWidget(pair(field("Grid nodes", self.sp_nodes), self.chk_force_num))
        self.chk_cooling = QCheckBox("Linear cooling rather than isothermal")
        self.sp_Tend = _spin(300, 2000, 900, 1, 5, " C")
        self.sp_Tend.setEnabled(False)
        self.chk_cooling.toggled.connect(self.sp_Tend.setEnabled)
        nbody.addWidget(self.chk_cooling)
        nbody.addWidget(field("Final temperature", self.sp_Tend))
        nbody.addWidget(note(
            "Diffusor uses a closed-form solution whenever one is exactly valid and the "
            "Crank-Nicolson solver otherwise, and tells you which it used.", "Hint"))
        self._update_geometry_note()
        return page_body(g, b, i, n)

    # ================================================================ step 5
    def _page_coefficient(self) -> QWidget:
        c, body = card("Diffusion coefficient",
                       "Tick one to fit with. Tick several and Diffusor will fit them all and "
                       "show the spread, which is usually larger than the analytical "
                       "uncertainty.")
        self.lst_coef = QListWidget()
        self.lst_coef.setMinimumHeight(230)
        self.lst_coef.currentItemChanged.connect(self._coefficient_selected)
        body.addWidget(self.lst_coef)
        b = ghost_button("Show the full equation, ranges and citation")
        b.clicked.connect(lambda: self.show_coefficient_info())
        body.addWidget(row(b, stretch_last=True))
        self.lbl_coef = note("", "Hint")
        body.addWidget(self.lbl_coef)
        return page_body(c)

    # ================================================================ step 6
    def _page_uncertainty(self) -> QWidget:
        c, body = card("Monte Carlo",
                       "Every draw re-runs the whole fit, so correlations between the sampled "
                       "quantities are honoured exactly rather than assumed away.")
        self.sp_draws = QSpinBox(); self.sp_draws.setRange(20, 100000); self.sp_draws.setValue(500)
        self.sp_seed = QSpinBox(); self.sp_seed.setRange(0, 2 ** 31 - 1); self.sp_seed.setValue(12345)
        body.addWidget(pair(field("Draws", self.sp_draws, "500 is usually enough; 1000 for a paper"),
                            field("Random seed", self.sp_seed, "recorded, so runs are reproducible")))

        s, sbody = card("What to sample")
        self.chk_mc_T = QCheckBox("Temperature"); self.chk_mc_T.setChecked(True)
        self.chk_mc_f = QCheckBox("Oxygen fugacity"); self.chk_mc_f.setChecked(True)
        self.chk_mc_P = QCheckBox("Pressure"); self.chk_mc_P.setChecked(True)
        self.chk_mc_D = QCheckBox("Diffusion coefficient"); self.chk_mc_D.setChecked(True)
        self.chk_mc_n = QCheckBox("Measurement noise"); self.chk_mc_n.setChecked(True)
        self.chk_mc_x = QCheckBox("Distance scale")
        for w in (self.chk_mc_T, self.chk_mc_f, self.chk_mc_P, self.chk_mc_D,
                  self.chk_mc_n, self.chk_mc_x):
            sbody.addWidget(w)
        self.chk_contrib = QCheckBox("Also rank how much each source contributes (slower)")
        sbody.addWidget(self.chk_contrib)

        d, dbody = card("How the diffusion coefficient is sampled")
        self.cmb_dmode = QComboBox()
        self.cmb_dmode.addItems(["auto", "covariance", "logD_at_T", "independent"])
        dbody.addWidget(field("Sampling mode", self.cmb_dmode))
        dbody.addWidget(note(
            "auto picks the best available for the chosen coefficient. covariance samples the "
            "published parameter covariance. logD_at_T samples ln D directly at the working "
            "temperature using the scatter the paper reports. independent samples each "
            "Arrhenius parameter on its own, which ignores the strong correlation between "
            "ln D0 and Q and inflates the uncertainty; it is offered only so older published "
            "estimates can be reproduced.", "Hint"))
        return page_body(c, s, d)

    # ================================================================ results
    def _page_results(self) -> QWidget:
        w = QWidget(); w.setObjectName("Page")
        v = QVBoxLayout(w)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(0)

        split = QSplitter(Qt.Horizontal)
        # --- summary sidebar
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
        bh.setContentsMargins(18, 11, 18, 11)
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
        # setParent(None) detaches immediately; deleteLater alone leaves the old
        # widgets on screen until the event loop catches up, which paints stale
        # text over the new summary.
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
                                   else "synthetic, not a measurement")
        else:
            self._summary_line("", "no data loaded")

        self._summary_head("Mineral", 1)
        self._summary_line(get_mineral(self.cmb_mineral.currentData()).name,
                           f"{self.cmb_species.currentText()}, "
                           f"{self.cmb_axis.currentText().lower()}")

        self._summary_head("Conditions", 2)
        self._summary_line("temperature", f"{self.sp_T.value():.0f} C  +/- {self.sp_T_sig.value():.0f}")
        self._summary_line("pressure", f"{self.sp_P.value():.0f} MPa  +/- {self.sp_P_sig.value():.0f}")
        if self.cmb_fo2_mode.currentIndex() == 0:
            self._summary_line("oxygen fugacity",
                               f"{self.cmb_buffer.currentText()} {self.sp_dbuf.value():+.2f} "
                               f"+/- {self.sp_dbuf_sig.value():.2f}")
        else:
            self._summary_line("oxygen fugacity", f"log10 fO2 = {self.sp_dbuf.value():.2f} bar")
        if self.sp_beam.value() > 0:
            self._summary_line("beam sigma", f"{self.sp_beam.value():.2f} um")

        self._summary_head("Model", 3)
        self._summary_line("geometry", self.cmb_geom.currentText())
        if self.cmb_ic.currentIndex() == 1 and self.an_values is not None:
            self._summary_line("initial condition",
                               "equilibrium profile from the anorthite gradient")
        else:
            self._summary_line("initial condition",
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
            lab.setObjectName("StepDotActive" if i == index
                              else ("StepDotDone" if i < index else "StepDot"))
            lab.style().unpolish(lab); lab.style().polish(lab)
        self.btn_back.setVisible(index > 0)
        if index == RESULTS:
            self.btn_next.setVisible(False)
            self.lbl_footer.setText("Change anything from the summary on the left.")
            self._rebuild_summary()
        else:
            self.btn_next.setVisible(True)
            self.btn_next.setText("Review and fit" if index == RESULTS - 1 else "Continue")
            self.lbl_footer.setText(f"Step {index + 1} of {RESULTS}")
        self.btn_next.setEnabled(self.profile is not None or index == 0)

    def _next(self):
        if self.step == 0 and self.profile is None:
            QMessageBox.information(self, "No data",
                                    "Load a file or one of the example datasets first.")
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
            QMessageBox.warning(self, "File missing",
                                f"{d.filename} is not in the examples folder. "
                                "Run scripts/make_examples.py to regenerate it.")
            return
        self._load(d.path, dataset=d)

    def _load(self, path: Path, dataset: Optional[ds.ExampleDataset]):
        try:
            df = read_table(path)
            if dataset is not None:
                spec = ProfileSpec(**dataset.spec)
            else:
                dlg = ColumnDialog(df, suggest_spec(df), self)
                if dlg.exec() != QDialog.Accepted:
                    return
                spec = dlg.spec()
            self.profile = build_profile(df, spec, source=str(path))
            self.dataset = dataset
            self.an_values = None
            p = self.profile
            head = f"{path.name}: {len(p)} points, {p.x.min():.2f} to {p.x.max():.2f} um"
            if dataset is not None:
                head += "\n\n" + dataset.provenance_banner()
                self._apply_dataset_settings(dataset, df)
            self.lbl_data.setText(head + "\n" + p.spec.describe())
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
        """Pre-fill the steps with the settings that go with a bundled dataset."""
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
        else:
            self.cmb_fo2_mode.setCurrentIndex(0)
            self.cmb_buffer.setCurrentText(s.get("buffer", "NNO"))
            self.sp_dbuf.setValue(s.get("delta_buffer", 0.0))
        self.sp_dbuf_sig.setValue(s.get("sigma_delta", 0.3))
        axis = s.get("axis")
        self.cmb_axis.setCurrentIndex({"a": 1, "b": 2, "c": 3}.get(axis, 0))
        self.cmb_geom.setCurrentText(s.get("geometry", "plane"))
        self.chk_comp_dep.setChecked(bool(s.get("composition_dependent", False)))
        if "x_composition" in s:
            self.sp_xcomp.setValue(s["x_composition"])
        self.cmb_ic.setCurrentIndex(
            1 if s.get("initial_condition") == "equilibrium_plag" else 0)
        want = s.get("coefficient")
        self._refresh_coefficients()
        for i in range(self.lst_coef.count()):
            it = self.lst_coef.item(i)
            it.setCheckState(Qt.Checked if it.data(Qt.UserRole) == want else Qt.Unchecked)
            if it.data(Qt.UserRole) == want:
                self.lst_coef.setCurrentRow(i)
        an_col = s.get("an_column")
        if an_col and an_col in df.columns:
            vals = np.asarray(df[an_col], dtype=float)
            if s.get("an_is_percent"):
                vals = vals / 100.0
            self.an_values = vals
            self.sp_xcomp.setValue(float(np.nanmean(vals)))
            self._log(f"  anorthite profile read from '{an_col}' "
                      f"(X_An {np.nanmin(vals):.2f} to {np.nanmax(vals):.2f})")

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
        self.lbl_mineral_note.setText(
            (f"Modelled variable: {mineral.composition_variable}. " if mineral.composition_variable
             else "") + mineral.notes)
        self._refresh_coefficients()
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
        if c.superseded_by or c.superseded_note:
            parts.append("SUPERSEDED. " + c.superseded_note)
        if not c.verified:
            parts.append("Coefficients not verified against the primary publication: "
                         + c.verified_from)
        self.lbl_coef.setText("\n\n".join(parts))
        self.lbl_coef.setObjectName(
            "Warn" if (c.superseded_by or c.superseded_note or not c.verified) else "Hint")
        self.lbl_coef.style().unpolish(self.lbl_coef)
        self.lbl_coef.style().polish(self.lbl_coef)

    def _on_ic_changed(self):
        equil = self.cmb_ic.currentIndex() == 1
        self._step_controls.setVisible(not equil)
        self.chk_free_x0.setEnabled(not equil)
        self.chk_free_plateaus.setEnabled(not equil)
        if equil:
            ok = self.an_values is not None
            msg = ("C(x) = C0 exp(A_i X_An(x) / R T), the quasi-steady state a fast trace "
                   "element relaxes to while the anorthite profile stays frozen (Dohmen, Faak "
                   "& Blundy 2017, Appendix eq. A13; Costa et al. 2003; Zellmer et al. 1999). "
                   "Diffusion then runs from this profile towards the boundary values.")
            if not ok:
                msg += ("\n\nNo anorthite column is loaded, so this cannot be used. Load a "
                        "dataset that carries one.")
            self.lbl_ic.setText(msg)
            self.lbl_ic.setObjectName("Hint" if ok else "Warn")
        else:
            self.lbl_ic.setText("Two plateaus meeting at a sharp interface: the classic "
                                "diffusion couple (Crank 1975 eq. 2.14).")
            self.lbl_ic.setObjectName("Hint")
        self.lbl_ic.style().unpolish(self.lbl_ic)
        self.lbl_ic.style().polish(self.lbl_ic)

    def _update_geometry_note(self):
        self.lbl_geom.setText(Geometry(self.cmb_geom.currentText()).describe())

    def _update_fo2_label(self):
        try:
            T = self.sp_T.value() + 273.15
            P = self.sp_P.value() * 1e6
            if self.cmb_fo2_mode.currentIndex() == 0:
                lf = log_fo2_from_delta(self.cmb_buffer.currentText(), self.sp_dbuf.value(), T, P)
                txt = (f"log10 fO2 = {lf:.3f} bar = {lf + 5:.3f} Pa, from Frost (1991) Table 1")
            else:
                lf = self.sp_dbuf.value()
                txt = f"log10 fO2 = {lf:.3f} bar = {lf + 5:.3f} Pa (absolute)"
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
        self._log("initial condition: " + ic.describe())

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
        if self.cmb_ic.currentIndex() == 1 and self.an_values is not None:
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
        if self.cmb_ic.currentIndex() == 1 and self.an_values is not None:
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
            coefficient_mode=self.cmb_dmode.currentText(),
            sample_measurement_noise=self.chk_mc_n.isChecked(),
            sigma_distance_scale=self.sp_xscale_sig.value() if self.chk_mc_x.isChecked() else 0.0)

    def _y_label(self) -> str:
        mineral = get_mineral(self.cmb_mineral.currentData())
        if self.profile is not None and self.profile.spec and self.profile.spec.mode == "A":
            return str(self.profile.spec.column_a)
        return f"{mineral.composition_variable.label} ({self.cmb_species.currentText()})"

    # ================================================================ running
    def _busy(self, on: bool, maximum: int = 0):
        for b in (self.btn_fit, self.btn_cmp, self.btn_mc):
            b.setEnabled(not on)
        self.btn_stop.setEnabled(on)
        self.progress.setVisible(on)
        self.progress.setMaximum(maximum)
        self.progress.setValue(0)

    def stop_work(self):
        for w in self._workers:
            if hasattr(w, "abort"):
                w.abort()
        self._log("stop requested")

    def _ready(self) -> bool:
        if self.profile is None:
            QMessageBox.information(self, "No data", "Load a profile first.")
            return False
        if not self._checked_keys():
            QMessageBox.information(self, "No coefficient",
                                    "Choose a diffusion coefficient on step 5.")
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
        self._busy(True)
        w = FitWorker(model, p.x, p.C, p.sigma, self._free_parameters(), 1e2, 3.2e12)
        self._workers = [w]
        w.finished.connect(self._fit_done)
        w.failed.connect(self._work_failed)
        self._threads.append(start(w))

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
                                    "Go back to step 5 and tick two or more coefficients.")
            return
        try:
            models = {get_coefficient(k).label: self._model(k) for k in keys}
        except Exception:
            QMessageBox.critical(self, "Model error", traceback.format_exc())
            return
        p = self.profile
        self._busy(True, len(models))
        w = CompareWorker(models, p.x, p.C, p.sigma, self._free_parameters(), 1e2, 3.2e12)
        self._workers = [w]
        w.progress.connect(lambda i, n, k: (self.progress.setValue(i),
                                            self.statusBar().showMessage(f"{i}/{n}: {k}")))
        w.finished.connect(self._compare_done)
        w.failed.connect(self._work_failed)
        self._threads.append(start(w))

    def _compare_done(self, results):
        self._busy(False)
        self.compare_results = results
        ok = {k: r for k, r in results.items() if not isinstance(r, Exception)}
        lines = ["Comparison of diffusion coefficients", "=" * 62]
        for k, r in results.items():
            if isinstance(r, Exception):
                lines.append(f"{k}\n   FAILED: {r}")
                continue
            lines.append(f"{k}\n   t = {human_time(r.t_seconds):<14s} "
                         f"reduced chi2 = {r.stats.reduced_chi2:.3g}   "
                         f"R2 = {r.stats.r_squared:.4f}")
            c = r.model.coefficient
            if c.superseded_by or c.superseded_note:
                lines.append("   [superseded] " + c.superseded_note[:160])
            elif not c.verified:
                lines.append("   [coefficients not verified against the primary publication]")
        if len(ok) > 1:
            ts = [r.t_seconds for r in ok.values()]
            lines += ["", f"Spread: {human_time(min(ts))} to {human_time(max(ts))}, "
                          f"a factor of {max(ts)/min(ts):.1f} between the extremes.",
                      "That spread is the real uncertainty on the choice of coefficient, and it "
                      "is usually larger than the analytical uncertainty on any one of them."]
        self.plot.show_comparison(ok, y_label=self._y_label())
        if ok:
            self.fit_result = next(iter(ok.values()))
        self._rebuild_summary()
        self._log("comparison complete")
        _text_dialog(self, "Coefficient comparison", "\n".join(lines), 820, 520)

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
        self._busy(True, n)
        w = MonteCarloWorker(model, p.x, p.C, p.sigma, self._budget(), n,
                             self.sp_seed.value(), self._free_parameters(), 1e2, 3.2e12,
                             do_contributions=self.chk_contrib.isChecked(),
                             contribution_draws=max(40, n // 5))
        self._workers = [w]
        w.progress.connect(lambda i, tot: self.progress.setValue(i))
        w.stage.connect(lambda s: self.statusBar().showMessage(s))
        w.finished.connect(self._mc_done)
        w.failed.connect(self._work_failed)
        self._threads.append(start(w))
        self._log(f"Monte Carlo: {n} draws, seed {self.sp_seed.value()}, "
                  f"sampling {', '.join(self._budget().active_sources())}")

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

    def _work_failed(self, msg: str):
        self._busy(False)
        self._log(msg)
        QMessageBox.critical(self, "Calculation failed", msg)

    def _show_warnings(self, warnings):
        serious = [w for w in warnings
                   if "SUPERSEDED" in w or "NOT verified" in w or "semi-infinite" in w
                   or "resolution limit" in w]
        if serious:
            QMessageBox.warning(self, "Check these before you use the result",
                                "\n\n".join(serious))

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
            "Diffusion chronometry with analytical (Crank 1975) and numerical "
            "(Crank-Nicolson) solvers, a registry of literature diffusion coefficients, "
            "and Monte Carlo error propagation that keeps temperature, oxygen fugacity "
            "and the Arrhenius parameters correlated.<br><br>"
            "Every equation and coefficient carries its citation, and the Methods view "
            "lists the references used by the current run."))
