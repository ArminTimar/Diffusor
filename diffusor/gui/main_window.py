"""Diffusor's main window: a stepped flow ending in a results view.

The settings are collected one group at a time, then summarised in a narrow
sidebar next to the plot. Any group in that summary can be clicked to go back
and change it. Pages never scroll. Only lists and reading panes do.
"""
from __future__ import annotations

import json
import os
import traceback
from dataclasses import asdict, fields
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from PySide6.QtCore import QLocale, QSettings, Qt, QTimer, QUrl, Slot
from PySide6.QtGui import QAction, QDesktopServices, QFont, QGuiApplication
from PySide6.QtWidgets import (QApplication, QCheckBox, QComboBox, QDialog, QDialogButtonBox,
                               QDoubleSpinBox, QFileDialog, QFrame, QHBoxLayout, QHeaderView,
                               QLabel, QLineEdit, QListWidget, QListWidgetItem, QMainWindow,
                               QMessageBox, QProgressBar, QPushButton, QSizePolicy, QSpinBox,
                               QSplitter, QStackedWidget, QTableWidget, QTableWidgetItem,
                               QVBoxLayout, QWidget)

from .. import datasets as ds
from .. import updates
from ..coefficients import Conditions, get as get_coefficient, list_coefficients
from ..coefficients.plagioclase import (ACTIVITY_SETS, ACTIVITY_SPECIES, DEFAULT_ACTIVITY_SET,
                                         activity_A, activity_note)
from ..dataio import ProfileSpec, build_profile, read_table, save_results, suggest_spec
from ..dataio.image_profiles import is_extraction_workbook, read_extraction
from ..dataio.images import PILLOW_SUFFIXES, RAW_SUFFIXES
from ..fitting import DiffusionModel, UncertaintyBudget
from ..fitting.fit import T_MAX_DEFAULT as T_MAX, T_MIN_DEFAULT as T_MIN
from ..fitting.montecarlo import default_workers
from ..minerals import MINERALS, get_mineral
from ..references import cite
from ..solvers import Geometry, InitialCondition, dirichlet, zero_flux
from ..solvers.history import ThermalHistory
from ..solvers.initial import guess_step_from_data
from ..thermo import available_buffers, log_fo2_from_delta
from ..thermo.units import human_time
from . import format_help, richtext, theme
from .plot_widget import DataPreview, ProfilePlot
from .widgets import (EquationView, FitScrollArea, WrapLabel, callout, card, collapsible,
                      divider, field, ghost_button, install_no_wheel, note, page_columns, pair,
                      primary_button, row, scrollable)
from .workers import CompareWorker, FitWorker, MonteCarloWorker, UpdateWorker, start

# The coefficient comes straight after the mineral because it decides what the later
# steps need: which host composition, whether pressure and fO2 enter at all, and
# whether the solver has to follow the composition along the profile.
STEPS = ["Data", "Mineral", "Coefficient", "Conditions", "Model", "Uncertainty", "Results"]
DATA, MINERAL, COEFFICIENT, CONDITIONS, MODEL, UNCERTAINTY, RESULTS = range(len(STEPS))

# For these species the profile is the mineral's own composition variable (X_Fe, X_An,
# ...). For every other species it is a concentration (TiO2 wt%, Sr ppm) that says
# nothing about the host composition the law needs.
EXCHANGE_SPECIES = ("Fe-Mg", "Fe-Ti", "NaSi-CaAl", "Na-K")

# how many of the user's own profile files the Data step remembers
RECENT_MAX = 10

# picked in "Load profile", these open the image extractor instead
IMAGE_ONLY_SUFFIXES = set(PILLOW_SUFFIXES + RAW_SUFFIXES + (".hdr", ".npy"))

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
    ("Microprobe, stated beam diameter", "spot", 5.0, None,
     "Enter the beam diameter given in the analytical methods of the study. The "
     "starting value is only a placeholder."),
    ("LA-ICP-MS spot", "spot", 10.0, None,
     "Druitt et al. (2012) used a 10 um laser spot."),
    ("LA-ICP-MS line scan", "slit", 7.5, None,
     "Grocolas et al. (2025) scanned with a 7.5 um wide slit."),
    ("SIMS spot", "spot", 12.0, None,
     "Druitt et al. (2012) used 10 to 15 um ion beams."),
    ("Custom sigma", None, 0.0, None,
     "Standard deviation of the Gaussian beam profile (Ganguly et al. 1988)."),
]

# What each end of the profile is. The key is stored in dataset settings.
BOUNDARIES = [
    ("far", "Plateau continues"),
    ("rim_melt", "Crystal rim, held by the melt"),
    ("rim_closed", "Crystal rim, closed"),
    ("centre", "Crystal centre"),
]

SOURCE_NAMES = {
    "temperature": "temperature", "fo2": "oxygen fugacity", "pressure": "pressure",
    "diffusion_coefficient": "diffusion coefficient", "measurement_noise": "measurement noise",
    "distance_scale": "distance scale", "boundary_compositions": "plateau compositions",
}

SAMPLING_MODES = [
    ("auto", "Best available for this coefficient"),
    ("covariance", "Correlated D0 and Q"),
    ("logD_at_T", "Scatter of log D at the temperature"),
    ("independent", "Each parameter independently"),
    ("none", "Hold coefficient parameters fixed"),
]


def number_locale() -> QLocale:
    """Decimal point and no digit grouping, whatever the system locale."""
    loc = QLocale(QLocale.English, QLocale.UnitedStates)
    loc.setNumberOptions(QLocale.OmitGroupSeparator)
    return loc


def _spin(lo, hi, val, dec=3, step=1.0, suffix=""):
    s = QDoubleSpinBox()
    s.setLocale(number_locale())
    s.setRange(lo, hi)
    s.setDecimals(dec)
    s.setSingleStep(step)
    s.setValue(val)
    if suffix:
        s.setSuffix(suffix)
    return s


_APP_SETUP = {}


def _style_application(app):
    """Theme the application and stop the wheel changing values, once per process.

    The stylesheet goes on the application so dialogs, message boxes and tooltips
    match. Both are done once: re-applying a stylesheet re-polishes every widget
    that already exists, which gets slow when several windows are opened.
    """
    if app not in _APP_SETUP:
        app.setStyleSheet(theme.stylesheet())
        _APP_SETUP[app] = install_no_wheel(app)
    return _APP_SETUP[app]


def _repolish(w, name):
    w.setObjectName(name)
    w.style().unpolish(w)
    w.style().polish(w)


# ==================================================================== columns
class ColumnDialog(QDialog):
    """Confirm the guessed column mapping. Everything else sits under Advanced."""

    def __init__(self, df, spec: ProfileSpec, parent=None, an=(None, True)):
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
        self.an.setCurrentText(an[0] or "(none)")
        self.an_percent = QCheckBox("in mol%")
        self.an_percent.setChecked(bool(an[1]))
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
                self.preview.setItem(r, c, QTableWidgetItem(
                    f"{v:g}" if isinstance(v, (int, float, np.number)) else str(v)))

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
        # Numbers use a decimal point everywhere, whatever the system locale, so the
        # boxes match the input files and the plots.
        QLocale.setDefault(number_locale())
        self.setWindowTitle("Diffusor")
        screen = QGuiApplication.primaryScreen().availableGeometry()
        self.resize(min(1440, int(screen.width() * 0.95)), min(900, int(screen.height() * 0.92)))
        self._no_wheel = _style_application(QApplication.instance())

        self.profile = None
        self.dataset: Optional[ds.ExampleDataset] = None
        self.an_values = None
        # Where the representative composition came from: "default", "profile" (the
        # profile mean, written by the program and taken back when it stops applying)
        # or "explicit" (typed in, or set by an example or an anorthite column).
        self._xcomp_source = "default"
        self._xcomp_fallback = None     # the value before the profile mean replaced it
        self._xcomp_writing = False
        self.fit_result = None
        self.mc_result = None
        self.compare_results = None
        self._stale = False             # a setting changed after the result was computed
        self._jobs: List = []           # (thread, worker) pairs kept alive until finished
        self._update_job = None         # the (thread, worker) of a running update check
        self._update_manual = False     # the user asked for it, so say the outcome
        self._update_release: Optional[updates.ReleaseInfo] = None
        self._log_entries: List[Tuple[datetime, str, str]] = []
        self._prefill: Dict[int, QLabel] = {}
        self._mc_draws: List[dict] = []
        self._mc_meta: Optional[dict] = None
        self._loaded_key: Optional[str] = None
        self.step = 0

        self.plot = ProfilePlot()
        self.plot.points_cut.connect(self._on_points_cut)

        root = QWidget(); root.setObjectName("Page")
        outer = QVBoxLayout(root)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)
        outer.addWidget(self._build_header())
        outer.addWidget(divider())
        outer.addWidget(self._build_update_banner())
        self.pages = QStackedWidget()
        # Built in this order because later pages call into widgets of earlier ones
        # while they are made; added to the stack in the order of STEPS.
        built = {}
        for step, builder in ((DATA, self._page_data), (MINERAL, self._page_mineral),
                              (CONDITIONS, self._page_conditions), (MODEL, self._page_model),
                              (COEFFICIENT, self._page_coefficient),
                              (UNCERTAINTY, self._page_uncertainty),
                              (RESULTS, self._page_results)):
            built[step] = builder()
        for step in range(len(STEPS)):
            self.pages.addWidget(built[step])
        outer.addWidget(self.pages, 1)
        outer.addWidget(divider())
        outer.addWidget(self._build_footer())
        self.setCentralWidget(root)

        self._build_menu()
        self._on_mineral_changed()
        self._on_resolution_changed()
        self._watch_inputs()
        self._go(DATA)

    # ================================================================ chrome
    def _build_header(self) -> QWidget:
        w = QWidget(); w.setObjectName("Page")
        h = QHBoxLayout(w)
        h.setContentsMargins(24, 8, 24, 8)
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
            lab.setAttribute(Qt.WA_Hover, True)
            lab.setCursor(Qt.PointingHandCursor)
            lab.mousePressEvent = (lambda e, idx=i: self._go(idx))
            self.step_labels.append(lab)
            rh.addWidget(lab)
        h.addWidget(rail)
        return w

    def _build_update_banner(self) -> QWidget:
        """A strip under the header, hidden until a newer release is found."""
        wrap = QWidget(); wrap.setObjectName("Page")
        lay = QHBoxLayout(wrap)
        lay.setContentsMargins(24, 8, 24, 0)
        box = QFrame(); box.setObjectName("InfoBox")
        h = QHBoxLayout(box)
        h.setContentsMargins(11, 6, 11, 6)
        h.setSpacing(10)
        self.lbl_update = WrapLabel("")
        self.lbl_update.setObjectName("CalloutText")
        self.lbl_update.setAlignment(Qt.AlignVCenter | Qt.AlignLeft)
        h.addWidget(self.lbl_update, 1)
        self.btn_update_get = primary_button("Download")
        self.btn_update_get.clicked.connect(self._open_update_page)
        h.addWidget(self.btn_update_get, 0, Qt.AlignVCenter)
        self.btn_update_later = ghost_button("Later")
        self.btn_update_later.clicked.connect(lambda: self.update_banner.setVisible(False))
        h.addWidget(self.btn_update_later, 0, Qt.AlignVCenter)
        lay.addWidget(box)
        wrap.setVisible(False)
        self.update_banner = wrap
        return wrap

    def _build_footer(self) -> QWidget:
        w = QWidget(); w.setObjectName("Page")
        h = QHBoxLayout(w)
        h.setContentsMargins(24, 8, 24, 8)
        h.setSpacing(12)
        self.btn_back = QPushButton("Back")
        self.btn_back.clicked.connect(lambda: self._go(self.step - 1))
        h.addWidget(self.btn_back)
        self.lbl_footer = QLabel("")
        self.lbl_footer.setObjectName("Sub")
        h.addWidget(self.lbl_footer)
        h.addStretch(1)
        self.lbl_status = QLabel("")
        self.lbl_status.setObjectName("Hint")
        h.addWidget(self.lbl_status)
        self.progress = QProgressBar()
        self.progress.setVisible(False)
        self.progress.setFixedWidth(180)
        h.addWidget(self.progress)
        self.btn_next = primary_button("Continue")
        self.btn_next.clicked.connect(self._next)
        h.addWidget(self.btn_next)
        return w

    def _status(self, msg: str):
        self.lbl_status.setText(str(msg).splitlines()[0][:120] if msg else "")

    def _build_menu(self):
        m = self.menuBar().addMenu("&File")
        for text, fn in (("Load profile...", self.load_file),
                         ("Extract profile from image...", self.show_image_extractor),
                         ("Export results...", self.export_results),
                         ("Export figure for Inkscape or CorelDRAW...", self.export_figure),
                         ("Shared-duration study...", self.show_joint_study),
                         ("Multicomponent and isotope study...", self.show_multicomponent_study),
                         (None, None),
                         ("Quit", self.close)):
            if text is None:
                m.addSeparator(); continue
            a = QAction(text, self); a.triggered.connect(fn); m.addAction(a)
        v = self.menuBar().addMenu("&View")
        for text, fn in (("Methods and references", self.show_methods),
                         ("Diffusion coefficient details", lambda: self.show_coefficient_info()),
                         ("Monte Carlo", self.show_histogram),
                         ("Log", self.show_log)):
            a = QAction(text, self); a.triggered.connect(fn); v.addAction(a)
        h = self.menuBar().addMenu("&Help")
        for text, fn in (("Input file format", self.show_format),
                         ("Choosing the boundaries", self.show_boundary_help),
                         ("All references", self.show_all_references),
                         ("About", self.show_about)):
            a = QAction(text, self); a.triggered.connect(fn); h.addAction(a)
        h.addSeparator()
        a = QAction("Check for updates...", self)
        a.triggered.connect(lambda: self.check_for_updates(manual=True))
        h.addAction(a)
        self.act_update_startup = QAction("Check for updates at startup", self)
        self.act_update_startup.setCheckable(True)
        self.act_update_startup.setChecked(self._settings().value(
            "updates/check_at_startup", True, type=bool))
        self.act_update_startup.toggled.connect(
            lambda on: self._settings().setValue("updates/check_at_startup", bool(on)))
        h.addAction(self.act_update_startup)

    def _prefill_bar(self, step: int, with_sources: bool = False) -> QWidget:
        """A one-line note saying what a loaded example set on this step."""
        bar = QWidget()
        h = QHBoxLayout(bar)
        h.setContentsMargins(0, 0, 0, 0)
        h.setSpacing(8)
        lab = WrapLabel("")
        lab.setObjectName("Good")
        lab.setVisible(False)
        self._prefill[step] = lab
        h.addWidget(lab, 1)
        if with_sources:
            self.btn_sources = ghost_button("Where these values come from")
            self.btn_sources.clicked.connect(self.show_example_details)
            self.btn_sources.setVisible(False)
            h.addWidget(self.btn_sources, 0, Qt.AlignTop)
        return bar

    # ================================================================ step 1
    def _page_data(self) -> QWidget:
        load_card, body = card("Load a profile")
        btn = primary_button("Choose a file...")
        btn.clicked.connect(self.load_file)
        img = QPushButton("From an image...")
        img.setToolTip("Draw a boundary on a BSE image or element map and read profiles "
                       "perpendicular to it")
        img.clicked.connect(self.show_image_extractor)
        fmt = ghost_button("File format")
        fmt.clicked.connect(self.show_format)
        body.addWidget(row(btn, img, fmt))
        body.addWidget(note(format_help.SHORT, "Hint"))

        # One card, two views: the loaded profile, or the recently loaded files.
        # Loading anything brings the profile to the front; the button switches.
        data_card, dbody = card()
        self.lbl_data_title = QLabel("Loaded data")
        self.lbl_data_title.setObjectName("H2")
        self.btn_data_view = ghost_button("Recent profiles")
        self.btn_data_view.clicked.connect(
            lambda: self._show_data_view(recent=self.data_stack.currentIndex() == 0))
        head = QWidget()
        hh = QHBoxLayout(head)
        hh.setContentsMargins(0, 0, 0, 0)
        hh.addWidget(self.lbl_data_title)
        hh.addStretch(1)
        hh.addWidget(self.btn_data_view)
        dbody.addWidget(head)
        self.data_stack = QStackedWidget()
        shown = QWidget()
        sv = QVBoxLayout(shown)
        sv.setContentsMargins(0, 0, 0, 0)
        sv.setSpacing(8)
        self.lbl_data = note("Nothing loaded yet.", "Hint")
        sv.addWidget(self.lbl_data)
        self.preview = DataPreview()
        self.preview.setMinimumHeight(180)
        sv.addWidget(self.preview, 1)
        self.data_stack.addWidget(shown)
        recent = QWidget()
        rv = QVBoxLayout(recent)
        rv.setContentsMargins(0, 0, 0, 0)
        rv.setSpacing(8)
        self.lbl_recent = note("", "Hint")
        rv.addWidget(self.lbl_recent)
        self.lst_recent = QListWidget()
        self.lst_recent.setMinimumHeight(120)
        # long folder paths are shortened in the middle rather than scrolled sideways
        self.lst_recent.setTextElideMode(Qt.ElideMiddle)
        self.lst_recent.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.lst_recent.itemDoubleClicked.connect(self._load_recent_item)
        rv.addWidget(self.lst_recent, 1)
        self.btn_recent_open = QPushButton("Load")
        self.btn_recent_open.clicked.connect(
            lambda: self._load_recent_item(self.lst_recent.currentItem()))
        self.btn_recent_clear = ghost_button("Clear the list")
        self.btn_recent_clear.clicked.connect(self._clear_recent)
        rv.addWidget(row(self.btn_recent_open, self.btn_recent_clear))
        self.data_stack.addWidget(recent)
        dbody.addWidget(self.data_stack, 1)
        self._refresh_recent()
        self._show_data_view(recent=bool(self._recent()))

        ex_card, ebody = card("Examples")
        self.lst_examples = QListWidget()
        self.lst_examples.setMouseTracking(True)
        for d in ds.DATASETS:
            it = QListWidgetItem(self._example_label(d, loaded=False))
            it.setData(Qt.UserRole, d.key)
            if not d.exists:
                it.setFlags(it.flags() & ~Qt.ItemIsEnabled)
            self.lst_examples.addItem(it)
        self.lst_examples.currentItemChanged.connect(self._example_selected)
        self.lst_examples.itemDoubleClicked.connect(lambda _it: self.load_example())
        ebody.addWidget(self.lst_examples, 1)
        self.lbl_example = note("Pick an example to see where it comes from.", "Hint")
        ebody.addWidget(self.lbl_example)
        b = QPushButton("Load example")
        b.clicked.connect(self.load_example)
        info = ghost_button("Where the data come from")
        info.clicked.connect(self.show_example_details)
        ebody.addWidget(row(b, info))

        # how the profile was measured belongs with the profile
        r, rbody = card("Analytical resolution")
        self.cmb_resolution = QComboBox()
        for label, *_ in RESOLUTION_PRESETS:
            self.cmb_resolution.addItem(label)
        self.cmb_resolution.currentIndexChanged.connect(self._on_resolution_changed)
        self.sp_width = _spin(0, 200, 5.0, 1, 0.5, " um")
        self.sp_width.valueChanged.connect(self._update_beam_sigma)
        self.sp_beam = _spin(0, 50, 0.0, 2, 0.1, " um")
        self.sp_xscale_sig = _spin(0, 0.5, 0.0, 3, 0.005)
        # two fields a row: this card shares the Data step with the preview
        rbody.addWidget(pair(field("How the profile was measured", self.cmb_resolution),
                             field("Distance scale error (relative, 1σ)", self.sp_xscale_sig)))
        self._width_field = field("Spot or slit width", self.sp_width)
        rbody.addWidget(pair(self._width_field, field("Beam σ", self.sp_beam)))
        self.lbl_resolution = note("", "Hint")
        rbody.addWidget(self.lbl_resolution)
        return page_columns([load_card, (data_card, 1)], [(ex_card, 1), r])

    @staticmethod
    def _example_label(d, loaded: bool) -> str:
        kind = "measured" if d.kind == "measured" else "synthetic"
        text = f"{d.name.split(' (')[0]}   ·   {kind}"
        return text + ("   ·   loaded" if loaded else "")

    def _mark_loaded_example(self, key: Optional[str]):
        self._loaded_key = key
        for i in range(self.lst_examples.count()):
            it = self.lst_examples.item(i)
            d = ds.get(it.data(Qt.UserRole))
            on = d.key == key
            it.setText(self._example_label(d, on))
            f = it.font()
            f.setBold(on)
            it.setFont(f)

    def _example_selected(self, item, _prev=None):
        if item is None:
            return
        d = ds.get(item.data(Qt.UserRole))
        src = (f"Measured. {cite(d.citation)}." if d.kind == "measured"
               else "Synthetic, made by Diffusor with a known answer.")
        mineral = get_mineral(d.mineral).name
        self.lbl_example.setText(f"{src} {mineral}, {d.species}.")

    def show_example_details(self):
        d = None
        item = self.lst_examples.currentItem()
        if self.step != DATA and self.dataset is not None:
            d = self.dataset
        elif item is not None:
            d = ds.get(item.data(Qt.UserRole))
        elif self.dataset is not None:
            d = self.dataset
        if d is None:
            return
        coef = None
        try:
            coef = get_coefficient(d.settings.get("coefficient"))
        except Exception:
            pass
        richtext.show(self, d.name, richtext.example_html(d, coef), 780, 640)

    def show_format(self):
        richtext.show(self, "Input file format", richtext.format_html(), 720, 560)

    def show_boundary_help(self):
        richtext.show(self, "Choosing the boundaries", richtext.boundaries_html(), 860, 680)

    def show_all_references(self):
        from .reference_list import ReferenceListDialog
        ReferenceListDialog(self).exec()

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
        self.lbl_axis = note("", "Hint")
        obody.addWidget(self.lbl_axis)
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

        return page_columns([c], [o], top=[self._prefill_bar(MINERAL)])

    # ================================================================ step 4
    def _page_conditions(self) -> QWidget:
        c, body = card("Temperature and pressure")
        self.sp_T = _spin(300, 2000, 950, 1, 5, " °C")
        self.sp_T.valueChanged.connect(self._update_fo2_label)
        self.sp_T_sig = _spin(0, 300, 20, 1, 1, " K")
        body.addWidget(pair(field("Temperature", self.sp_T), field("± 1σ", self.sp_T_sig)))
        self.chk_cooling = QCheckBox("Linear cooling to")
        self.sp_Tend = _spin(300, 2000, 900, 1, 5, " °C")
        self.sp_Tend.setEnabled(False)
        self.chk_cooling.toggled.connect(self.sp_Tend.setEnabled)
        self.chk_cooling.toggled.connect(self._update_solver_note)
        body.addWidget(row(self.chk_cooling, self.sp_Tend, spacing=10))
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
        self.box_unused = callout("Info")
        self.lbl_unused = self.box_unused.label
        self.box_unused.setVisible(False)

        x, xbody = card("Composition")
        self.sp_xcomp = _spin(0, 1, 0.15, 4, 0.01)
        self.sp_xcomp.valueChanged.connect(self._on_xcomp_changed)
        xbody.addWidget(field("Representative value", self.sp_xcomp))
        self.chk_comp_dep = QCheckBox("D follows the composition along the profile")
        self.chk_comp_dep.setChecked(True)
        self.chk_comp_dep.toggled.connect(self._update_solver_note)
        xbody.addWidget(self.chk_comp_dep)
        self.cmb_ol_coordinate = QComboBox()
        for text, data in (("Fe fraction: XFe", (1., 0.)),
                           ("Forsterite fraction: XFo", (-1., 1.)),
                           ("Forsterite mol%: Fo", (-.01, 1.))):
            self.cmb_ol_coordinate.addItem(text, data)
        self.cmb_ol_coordinate.setToolTip("Defines how measured olivine Fe-Mg values map to XFe in the diffusion law.")
        self.cmb_ol_coordinate.currentIndexChanged.connect(self._update_solver_note)
        self.cmb_ol_coordinate.currentIndexChanged.connect(lambda: self._sync_composition())
        xbody.addWidget(self.cmb_ol_coordinate)
        self.lbl_comp = note("", "Hint")
        xbody.addWidget(self.lbl_comp)
        return page_columns([c, f], [x, self.box_unused],
                            top=[self._prefill_bar(CONDITIONS, with_sources=True)])

    # ================================================================ step 5
    def _page_model(self) -> QWidget:
        g, gbody = card("Geometry")
        self.cmb_geom = QComboBox()
        self.cmb_geom.addItems(["plane", "cylinder", "sphere"])
        self.cmb_geom.currentIndexChanged.connect(self._update_geometry_note)
        gbody.addWidget(self.cmb_geom)
        self.lbl_geom = note("", "Hint")
        gbody.addWidget(self.lbl_geom)

        b, bbody = card("Ends of the profile")
        self.cmb_bcl = QComboBox()
        self.cmb_bcr = QComboBox()
        for cmb in (self.cmb_bcl, self.cmb_bcr):
            for key, label in BOUNDARIES:
                cmb.addItem(label, key)
            cmb.currentIndexChanged.connect(self._on_boundaries_changed)
        bbody.addWidget(pair(field("Left end", self.cmb_bcl), field("Right end", self.cmb_bcr)))
        self.lbl_bc = note("", "Hint")
        bbody.addWidget(self.lbl_bc)
        how = ghost_button("How to choose")
        how.clicked.connect(self.show_boundary_help)
        bbody.addWidget(row(how))

        i, ibody = card("Initial profile")
        self.cmb_ic = QComboBox()
        self.cmb_ic.currentIndexChanged.connect(self._on_ic_changed)
        ibody.addWidget(self.cmb_ic)
        self.lbl_ic = note("", "Hint")
        ibody.addWidget(self.lbl_ic)
        # which published A_i set couples the trace element to the anorthite gradient
        self.cmb_activity = QComboBox()
        for key, s in ACTIVITY_SETS.items():
            self.cmb_activity.addItem(s["label"], key)
        self.cmb_activity.setCurrentIndex(self.cmb_activity.findData(DEFAULT_ACTIVITY_SET))
        self.cmb_activity.currentIndexChanged.connect(self._on_ic_changed)
        self._activity_field = field("Anorthite activity factors A_i", self.cmb_activity)
        ibody.addWidget(self._activity_field)
        self._activity_field.setVisible(False)
        self.sp_x0 = _spin(-1e5, 1e5, 0.0, 3, 1, " um")
        self.sp_cl = _spin(-1e6, 1e6, 0.3, 5, 0.01)
        self.sp_cr = _spin(-1e6, 1e6, 0.18, 5, 0.01)
        self.sp_smooth = _spin(0, 50, 0.0, 2, 0.1, " um")
        self.sp_smooth.valueChanged.connect(self._update_solver_note)
        self._step_controls = QWidget()
        sc = QVBoxLayout(self._step_controls)
        sc.setContentsMargins(0, 0, 0, 0)
        sc.setSpacing(8)
        sc.addWidget(pair(field("Step position", self.sp_x0),
                          field("Initial smoothing", self.sp_smooth)))
        sc.addWidget(pair(field("Left plateau", self.sp_cl), field("Right plateau", self.sp_cr)))
        ibody.addWidget(self._step_controls)
        bg = QPushButton("Guess from the data")
        bg.clicked.connect(lambda: self._guess_initial(from_button=True))
        self.chk_free_x0 = QCheckBox("Fit the step position")
        self.chk_free_x0.setChecked(True)
        self.chk_free_plateaus = QCheckBox("Fit the plateaus")
        ibody.addWidget(row(bg, self.chk_free_x0, self.chk_free_plateaus, spacing=14))

        n, nbody = card("Solver")
        self.cmb_solver = QComboBox()
        self.cmb_solver.addItems(["Automatic", "Always numerical"])
        self.cmb_solver.currentIndexChanged.connect(self._update_solver_note)
        self.sp_nodes = QSpinBox(); self.sp_nodes.setRange(51, 4001); self.sp_nodes.setValue(301)
        self.sp_nodes.setLocale(number_locale())
        nbody.addWidget(pair(field("Solver", self.cmb_solver), field("Grid nodes", self.sp_nodes)))
        box = callout("Info")
        self.lbl_solver = box.label
        nbody.addWidget(box)
        self._update_geometry_note()
        self._on_boundaries_changed()
        return page_columns([g, n], [i, b], top=[self._prefill_bar(MODEL)])

    # ================================================================ step 3
    def _page_coefficient(self) -> QWidget:
        c, body = card("Diffusion coefficient")
        body.addWidget(note("Tick one to fit. Tick several to compare them.", "Hint"))
        self.lst_coef = QListWidget()
        self.lst_coef.setMouseTracking(True)
        self.lst_coef.currentItemChanged.connect(self._coefficient_selected)
        self.lst_coef.itemChanged.connect(lambda _it: self._on_coefficients_ticked())
        body.addWidget(self.lst_coef, 1)

        d, dbody = card("Selected")
        # The name, equation and notes scroll inside the card when the window is short,
        # so a long equation never makes the window taller than the screen. The button
        # stays below them, in view.
        inner = QWidget()
        ilay = QVBoxLayout(inner)
        ilay.setContentsMargins(0, 0, 0, 0)
        ilay.setSpacing(8)
        self.lbl_coef_name = note("", "H2")
        ilay.addWidget(self.lbl_coef_name)
        ilay.addSpacing(6)
        self.lbl_coef_eq = EquationView()
        ilay.addWidget(self.lbl_coef_eq)
        self.lbl_coef = note("", "Hint")
        ilay.addWidget(self.lbl_coef)
        self.lbl_coef_range = note("", "Hint")
        ilay.addWidget(self.lbl_coef_range)
        ilay.addStretch(1)
        area = FitScrollArea()
        area.setWidget(inner)
        area.setWidgetResizable(True)
        area.setFrameShape(QFrame.NoFrame)
        area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        dbody.addWidget(area, 1)
        b = ghost_button("Full details and references")
        b.clicked.connect(lambda: self.show_coefficient_info())
        dbody.addWidget(row(b))
        return page_columns([(c, 1)], [d], top=[self._prefill_bar(COEFFICIENT)])

    # ================================================================ step 6
    def _page_uncertainty(self) -> QWidget:
        c, body = card("Monte Carlo")
        self.sp_draws = QSpinBox(); self.sp_draws.setRange(20, 100000); self.sp_draws.setValue(500)
        self.sp_seed = QSpinBox(); self.sp_seed.setRange(0, 2 ** 31 - 1); self.sp_seed.setValue(12345)
        self.sp_draws.setLocale(number_locale())
        self.sp_seed.setLocale(number_locale())
        body.addWidget(pair(field("Draws", self.sp_draws), field("Random seed", self.sp_seed)))
        body.addWidget(note("500 for a quick look, 1000 for a paper.", "Hint"))
        cores = os.cpu_count() or 1
        self.sp_cores = QSpinBox(); self.sp_cores.setRange(1, max(1, min(cores, 61)))
        self.sp_cores.setValue(default_workers())
        body.addWidget(field("Processor cores", self.sp_cores,
                             f"Draws are fitted in parallel on this many of the {cores} cores. "
                             "The seed gives the same answer on any number. Short runs use one."))

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
        for key, label in SAMPLING_MODES:
            self.cmb_dmode.addItem(label, key)
        self.cmb_dmode.currentIndexChanged.connect(self._on_dmode_changed)
        dbody.addWidget(self.cmb_dmode)
        self.lbl_dmode = note("", "Hint")
        dbody.addWidget(self.lbl_dmode)
        self._on_dmode_changed()
        return page_columns([c, s], [d])

    def _on_dmode_changed(self):
        """Describe the sampling mode for the coefficient that will be run."""
        mode = self.cmb_dmode.currentData()
        coef = None
        if hasattr(self, "lst_coef"):
            keys = self._checked_keys()
            coef = get_coefficient(keys[0]) if keys else None
        name = coef.label.split(",", 1)[1].strip() if coef and "," in coef.label else "this law"
        texts = {
            "covariance": ("Draws D0 and Q together so a high D0 comes with a high Q. Only "
                           "Grocolas et al. (2025) Sr and Ba in plagioclase provide this. The "
                           "authors state that log D0 and Q are strongly correlated and treat "
                           "them as perfectly correlated. Diffusor copies that assumption. No "
                           "other paper in Diffusor publishes a covariance."),
            "logD_at_T": ("Draws log D at the working temperature from the scatter the paper "
                          "reports about its fit. D0 and Q stay on the published line."),
            "independent": ("Ignores the correlation between D0 and Q and overstates the "
                            "uncertainty. Use it only to reproduce older estimates."),
            "none": "Coefficient parameters are held fixed. Their uncertainty is not propagated.",
        }
        if mode == "auto":
            if coef is None:
                text = "Uses the most faithful option the chosen coefficient supports."
            else:
                best = coef.default_sampling_mode()
                detail = {"covariance": "correlated D0 and Q",
                          "logD_at_T": f"the scatter of log D, {coef.sigma_logD} log units",
                          "independent": "each parameter independently",
                          "none": ("fixed coefficient parameters. The paper gives neither a "
                                   "covariance nor a scatter of log D, so the coefficient "
                                   "uncertainty is not propagated")}[best]
                text = f"For {name} this is {detail}."
        else:
            text = texts[mode]
        self.lbl_dmode.setText(text)
        # grey out what the coefficient cannot do
        model = self.cmb_dmode.model()
        for i, (key, _label) in enumerate(SAMPLING_MODES):
            enabled = True
            if coef is not None and key == "covariance":
                enabled = coef.covariance is not None
            if coef is not None and key == "logD_at_T":
                enabled = coef.sigma_logD is not None
            model.item(i).setEnabled(enabled)
        if not model.item(self.cmb_dmode.currentIndex()).isEnabled():
            self.cmb_dmode.setCurrentIndex(0)

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
        self.summary_layout.setContentsMargins(18, 14, 18, 14)
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
        bh.setContentsMargins(18, 8, 18, 8)
        bh.setSpacing(8)
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
        self.btn_view_fit = ghost_button("Profile")
        self.btn_view_fit.clicked.connect(self.show_profile_view)
        self.btn_view_mc = ghost_button("Monte Carlo view")
        self.btn_view_mc.clicked.connect(self.show_histogram)
        bh.addWidget(self.btn_view_fit)
        bh.addWidget(self.btn_view_mc)
        bh.addWidget(self._vline())
        for text, fn in (("Details", self.show_run_details),
                         ("Methods", self.show_methods),
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

    @staticmethod
    def _vline() -> QFrame:
        f = QFrame()
        f.setFixedWidth(1)
        f.setFixedHeight(20)
        f.setStyleSheet(f"background:{theme.BORDER};")
        return f

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
                   f"{human_time(self.fit_result.t_seconds)}</b>"
                   f" <span style='color:{theme.TEXT_MUTED};font-size:11.5px'>best fit</span>")
            if self.mc_result is not None:
                txt += (f"<br><span style='color:{theme.TEXT}'>Monte Carlo median "
                        f"{human_time(self.mc_result.median)}</span>"
                        f"<br><span style='color:{theme.TEXT}'>68%: "
                        f"{human_time(self.mc_result.p16)} to "
                        f"{human_time(self.mc_result.p84)}</span>"
                        f"<br><span style='color:{theme.TEXT_MUTED};font-size:11.5px'>95%: "
                        f"{human_time(self.mc_result.p2_5)} to "
                        f"{human_time(self.mc_result.p97_5)}</span>")
            txt += (f"<br><span style='color:{theme.TEXT_FAINT};font-size:11px'>"
                    f"reduced chi2 {self.fit_result.stats.reduced_chi2:.2f} &middot; "
                    f"R2 {self.fit_result.stats.r_squared:.4f}</span>")
            if self._stale:
                txt += (f"<br><span style='color:{theme.WARN};font-size:11.5px'>"
                        "A setting changed after this was computed. Run Fit again.</span>")
            res.setText(txt)
            self.summary_layout.addWidget(res)
            self.summary_layout.addSpacing(6)
            self.summary_layout.addWidget(divider())

        self._summary_head("Data", DATA)
        if self.profile is not None:
            p = self.profile
            name = Path(p.source).name if p.source else "loaded profile"
            self._summary_line(name, f"{len(p) - p.n_excluded} points, {p.x.min():.1f} to "
                                     f"{p.x.max():.1f} um")
            if p.n_excluded:
                self._summary_line("cut by hand", f"{p.n_excluded} of {len(p)} points")
            if self.dataset is not None:
                self._summary_line("provenance",
                                   "measured, published" if self.dataset.kind == "measured"
                                   else "synthetic")
        else:
            self._summary_line("", "no data loaded")
        if self.sp_beam.value() > 0:
            self._summary_line("beam σ", f"{self.sp_beam.value():.2f} um "
                               f"({self.cmb_resolution.currentText().lower()})")

        self._summary_head("Mineral", MINERAL)
        self._summary_line(get_mineral(self.cmb_mineral.currentData()).name,
                           f"{self.cmb_species.currentText()}, "
                           f"{self.cmb_axis.currentText().lower()}")

        self._summary_head("Coefficient", COEFFICIENT)
        for k in self._checked_keys():
            c = get_coefficient(k)
            flag = ""
            if not c.verified:
                flag = "  [unverified]"
            if c.superseded_by or c.superseded_note:
                flag += "  [superseded]"
            self._summary_line("", c.label + flag)

        self._summary_head("Conditions", CONDITIONS)
        self._summary_line("temperature", f"{self.sp_T.value():.0f} °C ± {self.sp_T_sig.value():.0f}")
        if self.chk_cooling.isChecked():
            self._summary_line("cooling", f"linear to {self.sp_Tend.value():.0f} °C")
        self._summary_line("pressure", f"{self.sp_P.value():.0f} MPa ± {self.sp_P_sig.value():.0f}")
        if self.cmb_fo2_mode.currentIndex() == 0:
            self._summary_line("oxygen fugacity",
                               f"{self.cmb_buffer.currentText()} {self.sp_dbuf.value():+.2f} "
                               f"± {self.sp_dbuf_sig.value():.2f}")
        else:
            self._summary_line("oxygen fugacity", f"log fO2 {self.sp_dbuf.value():.2f} bar")
        keys = self._checked_keys()
        needs = get_coefficient(keys[0]).requires if keys else ()
        if needs:
            self._summary_line("host composition",
                               f"{', '.join(needs)} {self.sp_xcomp.value():.4g}"
                               + (", following the profile" if self.chk_comp_dep.isChecked()
                                  and self._model(keys[0]).comp_key else ""))

        self._summary_head("Model", MODEL)
        self._summary_line("geometry", self.cmb_geom.currentText())
        if self._equilibrium_ic():
            self._summary_line("initial profile", "equilibrium with the anorthite zoning")
        else:
            self._summary_line("initial profile",
                               f"step at {self.sp_x0.value():.1f} um, "
                               f"{self.sp_cl.value():.4g} to {self.sp_cr.value():.4g}")
        self._summary_line("ends",
                           f"{self.cmb_bcl.currentText().lower()} / "
                           f"{self.cmb_bcr.currentText().lower()}")
        self._summary_line("solver", self._solver_short())

        self._summary_head("Uncertainty", UNCERTAINTY)
        b = self._budget()
        self._summary_line(f"{self.sp_draws.value()} draws, seed {self.sp_seed.value()}",
                           ", ".join(SOURCE_NAMES.get(s, s) for s in b.active_sources())
                           or "nothing sampled")

        self.summary_layout.addStretch(1)

    # ================================================================ navigation
    def _go(self, index: int):
        height_before = self.height()
        index = max(DATA, min(index, RESULTS))
        if index > DATA and self.profile is None:
            self._status("Load a profile first")
            index = DATA
        self.step = index
        self.pages.setCurrentIndex(index)
        for i, lab in enumerate(self.step_labels):
            _repolish(lab, "StepDotActive" if i == index
                      else ("StepDotDone" if i < index else "StepDot"))
        self.btn_back.setVisible(index > DATA)
        if index == CONDITIONS:
            self._update_unused_note()
        if index == MODEL:
            self._update_solver_note()
        if index == UNCERTAINTY:
            self._on_dmode_changed()
        if index == RESULTS:
            self.btn_next.setVisible(False)
            self.lbl_footer.setText("Click edit in the summary to change a setting.")
            self._rebuild_summary()
        else:
            self.btn_next.setVisible(True)
            self.btn_next.setText("Review and fit" if index == RESULTS - 1 else "Continue")
            self.lbl_footer.setText(f"Step {index + 1} of {RESULTS}")
        self.btn_next.setEnabled(self.profile is not None or index == DATA)
        if self.isVisible():
            # the window as context drops the call if the window is deleted first
            QTimer.singleShot(60, self, lambda h=height_before: self._restore_height(h))

    def _restore_height(self, height: int):
        """Undo a growth caused by wrapped text measuring itself before it had a width.

        A newly shown page briefly asks for more height than it needs. Once its
        labels know their width the page fits again, so the window goes back to
        the height the user gave it.
        """
        if self.isMaximized() or self.isFullScreen():
            return
        if self.height() > height and self.minimumSizeHint().height() <= height:
            self.resize(self.width(), height)

    def _next(self):
        if self.step == DATA and self.profile is None:
            QMessageBox.information(self, "No data", "Load a file or an example first.")
            return
        self._go(self.step + 1)

    # ================================================================ data
    def load_file(self):
        start = str(self._settings().value("recent/last_dir", "") or "")
        path, _ = QFileDialog.getOpenFileName(
            self, "Load profile", start if start and Path(start).is_dir() else "",
            "Tables and image profiles (*.csv *.txt *.tsv *.xlsx *.xls);;All files (*)")
        if path:
            self._settings().setValue("recent/last_dir", str(Path(path).resolve().parent))
        if path and Path(path).suffix.lower() in IMAGE_ONLY_SUFFIXES:
            self.show_image_extractor()
            self.image_extractor.open_image(path)
            return
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

    def _load(self, path: Path, dataset: Optional[ds.ExampleDataset],
              from_recent: bool = False):
        if dataset is None and is_extraction_workbook(path):
            try:
                table = read_extraction(path)
            except Exception:
                QMessageBox.critical(self, "Could not load the file", traceback.format_exc())
                return
            if self.load_image_table(table):
                self._remember_file(path)
            return
        try:
            df = read_table(path)
            an = (None, True)
            if dataset is not None:
                spec = ProfileSpec(**dataset.spec)
            else:
                spec, an = self._column_mapping(df, path, from_recent)
                if spec is None:
                    return
        except Exception:
            self._log(traceback.format_exc(), "error")
            QMessageBox.critical(self, "Could not load the file", traceback.format_exc())
            return
        if self._use_table(df, spec, path, dataset, an) and dataset is None:
            self._remember_file(path, df, spec, an)

    # ---------------------------------------------------------------- recent profiles
    def _recent(self) -> List[dict]:
        """The remembered files, newest first: path, time, columns and mapping."""
        try:
            items = json.loads(str(self._settings().value("recent/profiles", "[]") or "[]"))
        except (TypeError, ValueError):
            return []
        return [r for r in items if isinstance(r, dict) and r.get("path")][:RECENT_MAX]

    @staticmethod
    def _same_path(a: str, b: str) -> bool:
        return os.path.normcase(os.path.abspath(a)) == os.path.normcase(os.path.abspath(b))

    def _remembered(self, path) -> Optional[dict]:
        return next((r for r in self._recent() if self._same_path(r["path"], str(path))), None)

    def _remember_file(self, path, df=None, spec: Optional[ProfileSpec] = None, an=(None, True)):
        entry = {"path": str(Path(path).resolve()),
                 "loaded": datetime.now().strftime("%Y-%m-%d %H:%M")}
        if df is not None and spec is not None:
            entry["columns"] = [str(c) for c in df.columns]
            entry["spec"] = asdict(spec)
            entry["an"] = [an[0], bool(an[1])]
        items = [r for r in self._recent() if not self._same_path(r["path"], entry["path"])]
        s = self._settings()
        s.setValue("recent/profiles", json.dumps([entry] + items[:RECENT_MAX - 1]))
        s.setValue("recent/last_dir", str(Path(entry["path"]).parent))
        self._refresh_recent()

    def _column_mapping(self, df, path, from_recent: bool):
        """The column mapping for a user's file: remembered, or confirmed in the dialog.

        A file opened from the recent list whose columns are unchanged loads with
        last time's mapping and no dialog. Otherwise the dialog opens, starting
        from last time's mapping when every column it names is still there.
        """
        cols = [str(c) for c in df.columns]
        known = self._remembered(path) or {}
        spec, an = None, (None, True)
        if known.get("spec"):
            try:
                names = {f.name for f in fields(ProfileSpec)}
                spec = ProfileSpec(**{k: v for k, v in known["spec"].items() if k in names})
                an = (known.get("an") or [None, True])[0], bool((known.get("an") or [None, True])[1])
            except Exception:
                spec, an = None, (None, True)
        if spec is not None:
            used = [spec.distance_column, spec.column_a, spec.column_b, spec.sigma_a_column,
                    spec.sigma_b_column, an[0]]
            if not all(c in cols for c in used if c):
                spec, an = None, (None, True)
        if spec is not None and from_recent and known.get("columns") == cols:
            self._log(f"columns as last time: {spec.describe()}")
            return spec, an
        dlg = ColumnDialog(df, spec or suggest_spec(df), self, an=an)
        if dlg.exec() != QDialog.Accepted:
            return None, (None, True)
        return dlg.spec(), dlg.an_column()

    def _refresh_recent(self):
        self.lst_recent.clear()
        items = self._recent()
        for r in items:
            path = Path(r["path"])
            here = path.exists()
            it = QListWidgetItem(f"{path.name}   ·   {r.get('loaded', '')}"
                                 + ("" if here else "   ·   not found") + f"\n{path.parent}")
            it.setData(Qt.UserRole, str(path))
            it.setToolTip(str(path))
            if not here:
                it.setFlags(it.flags() & ~Qt.ItemIsEnabled)
            self.lst_recent.addItem(it)
        self.lbl_recent.setText(
            "Double-click to load. Each file opens with the columns chosen last time."
            if items else "Files you load appear here, with the columns you chose for them.")
        self.btn_recent_open.setEnabled(bool(items))
        self.btn_recent_clear.setEnabled(bool(items))

    def _load_recent_item(self, item):
        if item is None or not (item.flags() & Qt.ItemIsEnabled):
            return
        path = Path(item.data(Qt.UserRole))
        if not path.exists():
            self._refresh_recent()
            return
        self._load(path, dataset=None, from_recent=True)

    def _clear_recent(self):
        self._settings().setValue("recent/profiles", "[]")
        self._refresh_recent()

    def _show_data_view(self, recent: bool):
        """The profile, or the recent files, in the Data step's left card."""
        self.data_stack.setCurrentIndex(1 if recent else 0)
        self.lbl_data_title.setText("Recent profiles" if recent else "Loaded data")
        self.btn_data_view.setText("Show the profile" if recent else "Recent profiles")

    def load_image_table(self, table) -> bool:
        """Ask for pixel size and value-to-composition map, then load the profile."""
        from .image_calibration import ImageCalibrationDialog
        dlg = ImageCalibrationDialog(table, self)
        if dlg.exec() != QDialog.Accepted:
            return False
        self._log(dlg.summary())
        path = Path(table.path or table.source or "image profile")
        self._use_table(dlg.frame, dlg.spec, path, None, (None, True))
        self.raise_()
        self.activateWindow()
        return self.profile is not None

    def show_image_extractor(self):
        from .image_extractor import ImageExtractorDialog
        if not hasattr(self, "image_extractor"):
            self.image_extractor = ImageExtractorDialog(self)
        self.image_extractor.show()
        self.image_extractor.raise_()

    def _use_table(self, df, spec, path: Path, dataset: Optional[ds.ExampleDataset],
                   an) -> bool:
        try:
            self.profile = build_profile(df, spec, source=str(path))
            self.dataset = dataset
            self.an_values = None
            self._release_composition()     # the last profile's mean is not this one's
            for lab in self._prefill.values():
                lab.setVisible(False)
            self.btn_sources.setVisible(False)
            p = self.profile
            if dataset is not None:
                self._apply_dataset_settings(dataset, df)
                msg = (f"<b>{dataset.name.split(' (')[0]}</b>, {len(p)} points from "
                       f"{p.x.min():.1f} to {p.x.max():.1f} um. The example also filled in the "
                       "mineral, conditions, model, coefficient and analytical resolution. "
                       "Each step shows what it set.")
                self._mark_loaded_example(dataset.key)
            else:
                msg = (f"<b>{path.name}</b>, {len(p)} points from {p.x.min():.1f} to "
                       f"{p.x.max():.1f} um.")
                self._mark_loaded_example(None)
                if an[0]:
                    vals = self.profile.column(an[0])
                    self.an_values = vals / 100.0 if an[1] else vals
                self._refresh_ic_options()
            self.lbl_data.setText(msg)
            _repolish(self.lbl_data, "Good")
            # graphs saved from the plot toolbar go beside the data, like exports
            self.plot.toolbar.start_dir = Path(self.output_dir()) if self.output_dir() else None
            self.plot.toolbar.default_name = f"{Path(p.source).stem or 'profile'}_diffusor.png"
            self._log(f"loaded {path}\n{p.spec.describe()}"
                      + "".join(f"\nnote: {n}" for n in p.notes))
            self._guess_initial()
            self.preview.show_data(p.x, p.C, p.sigma, y_label=self._y_label())
            self.plot.set_points(p.x, p.C, None)
            self.plot.show_data(p.x, p.C, p.sigma, y_label=self._y_label())
            self.fit_result = self.mc_result = self.compare_results = None
            self._mc_draws = []
            self.btn_next.setEnabled(True)
            self._status(f"Loaded {len(p)} points")
            self._show_data_view(recent=False)
            return True
        except Exception:
            self._log(traceback.format_exc(), "error")
            QMessageBox.critical(self, "Could not load the file", traceback.format_exc())
            return False

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
        fo2 += f" ± {s.get('sigma_delta', 0.3):g}"
        axis = s.get("axis")
        self.cmb_axis.setCurrentIndex({"a": 1, "b": 2, "c": 3}.get(axis, 0))
        self.cmb_geom.setCurrentText(s.get("geometry", "plane"))
        for cmb, key in ((self.cmb_bcl, "bc_left"), (self.cmb_bcr, "bc_right")):
            i = cmb.findData(s.get(key, "far"))
            cmb.setCurrentIndex(i if i >= 0 else 0)
        self.chk_comp_dep.setChecked(bool(s.get("composition_dependent", False)))
        self.chk_free_plateaus.setChecked(bool(s.get("fit_plateaus", False)))
        self.cmb_ol_coordinate.setCurrentIndex(s.get("olivine_coordinate", 0))
        if "x_composition" in s:
            self._set_xcomp(s["x_composition"], "explicit")
        want = s.get("coefficient")
        self._refresh_coefficients()
        self.lst_coef.blockSignals(True)
        for i in range(self.lst_coef.count()):
            it = self.lst_coef.item(i)
            it.setCheckState(Qt.Checked if it.data(Qt.UserRole) == want else Qt.Unchecked)
        self.lst_coef.blockSignals(False)
        for i in range(self.lst_coef.count()):
            if self.lst_coef.item(i).data(Qt.UserRole) == want:
                self.lst_coef.setCurrentRow(i)
        an_col = s.get("an_column")
        if an_col and an_col in df.columns:
            vals = self.profile.column(an_col)
            if s.get("an_is_percent"):
                vals = vals / 100.0
            self.an_values = vals
            self._set_xcomp(float(np.nanmean(vals)), "explicit")
            self._log(f"anorthite read from '{an_col}' "
                      f"(X_An {np.nanmin(vals):.2f} to {np.nanmax(vals):.2f})")
        self._refresh_ic_options()
        self._select_ic("equilibrium_plag" if s.get("initial_condition") == "equilibrium_plag"
                        else "step")
        labels = [p[0] for p in RESOLUTION_PRESETS]
        self.cmb_resolution.setCurrentIndex(labels.index(s.get("resolution", "No correction")))
        if "beam_width_um" in s:
            self.sp_width.setValue(s["beam_width_um"])

        name = d.name.split(" (")[0]
        head = f"Set by the {name} example: "
        self._show_prefill(MINERAL, head + f"{self.cmb_mineral.currentText()}, {d.species}, "
                           f"{self.cmb_axis.currentText().lower()}.")
        self._show_prefill(CONDITIONS, head + f"{s.get('T_C', 950):g} ± {s.get('sigma_T_K', 20):g} °C, "
                           f"{s.get('P_MPa', 200):g} ± {s.get('sigma_P_MPa', 100):g} MPa, {fo2}.")
        self.btn_sources.setVisible(True)
        if self._boundaries_far():
            ends = "plateaus continuing past both ends"
        else:
            ends = (f"left end {self.cmb_bcl.currentText().lower()}, right end "
                    f"{self.cmb_bcr.currentText().lower()}")
        fitted = ", plateaus fitted" if s.get("fit_plateaus") else ""
        self._show_prefill(MODEL, head + f"{s.get('geometry', 'plane')} geometry, "
                           f"{self.cmb_ic.currentText().lower()}, {ends}{fitted}.")
        if want:
            self._show_prefill(COEFFICIENT, head + get_coefficient(want).label + ".")

    def _show_prefill(self, step: int, text: str):
        lab = self._prefill.get(step)
        if lab is not None:
            lab.setText(text)
            lab.setVisible(True)

    # ================================================================ reactions
    def _log(self, msg: str, level: str = "info"):
        self._log_entries.append((datetime.now(), level, str(msg)))
        self._status(msg)

    def _on_mineral_changed(self):
        mineral = get_mineral(self.cmb_mineral.currentData())
        self.cmb_species.blockSignals(True)
        self.cmb_species.clear()
        self.cmb_species.addItems([sp for sp in mineral.species_keys()
                                   if list_coefficients(mineral.key, sp)])
        self.cmb_species.blockSignals(False)
        self.cmb_axis.setEnabled(not mineral.isotropic)
        self.lbl_axis.setText("Treated as isotropic, so the direction does not matter."
                              if mineral.isotropic else "")
        self.lbl_axis.setVisible(mineral.isotropic)
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
        self.lst_coef.blockSignals(True)
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
        self.lst_coef.blockSignals(False)
        if self.lst_coef.count():
            self.lst_coef.setCurrentRow(0)
        self._on_coefficients_ticked()
        self._sync_composition()

    def _on_coefficients_ticked(self):
        if not hasattr(self, "cmb_dmode"):
            return
        self._update_composition_note()
        self._update_unused_note()
        self._on_dmode_changed()
        self._update_solver_note()

    def _coefficient_selected(self, item, _prev=None):
        if item is None:
            return
        c = get_coefficient(item.data(Qt.UserRole))
        self.lbl_coef_name.setText(c.label)
        from ..coefficients.latex import coefficient_latex
        self.lbl_coef_eq.set_latex(coefficient_latex(c))
        self.lbl_coef_eq.setToolTip(c.equation_text)
        parts = []
        if c.superseded_by:
            parts.append(f"Superseded by {cite(c.superseded_by)}.")
        elif c.superseded_note:
            parts.append("Superseded. See the details.")
        if not c.verified:
            parts.append("Not checked against the original paper. See the details.")
        if parts:
            self.lbl_coef.setText(f"<span style='color:{theme.WARN}'>{' '.join(parts)}</span>")
        else:
            self.lbl_coef.setText("Checked against the original paper.")
        self.lbl_coef_range.setText(f"Calibrated for {c.T_range}."
                                    if c.T_range.lo is not None else "")

    def _update_composition_note(self):
        """Grey out the composition box for a law that does not use it."""
        keys = self._checked_keys() if self.lst_coef.count() else []
        coef = get_coefficient(keys[0]) if keys else None
        needed = coef is None or bool(coef.requires)
        for wdg in (self.sp_xcomp, self.chk_comp_dep):
            wdg.setEnabled(needed)
        self.lbl_comp.setText("" if needed else
                              f"{cite(coef.citation)} has no composition term, so this value "
                              "does not enter the model.")
        ol_exchange = coef is not None and coef.mineral == "olivine" and coef.species == "Fe-Mg"
        self.cmb_ol_coordinate.setVisible(ol_exchange)
        if needed and coef:
            self.lbl_comp.setText("Host composition required by this law: " + ", ".join(coef.requires) +
                                 ". Use mole fractions. Trace-element concentration is a separate variable.")

    def _update_unused_note(self):
        keys = self._checked_keys() if self.lst_coef.count() else []
        coef = get_coefficient(keys[0]) if keys else None
        text = richtext.unused_conditions(coef) if coef else ""
        self.lbl_unused.setText(text)
        self.box_unused.setVisible(bool(text))
        self._unused_conditions = richtext.unused_condition_names(coef) if coef else []
        self._update_condition_inputs()

    def _update_condition_inputs(self):
        """Switch off the pressure and fO2 inputs that the chosen law ignores.

        Pressure stays on when fO2 is given relative to a buffer and the law uses
        fO2, because the buffer itself moves with pressure.
        """
        unused = getattr(self, "_unused_conditions", [])
        fo2_used = "oxygen fugacity" not in unused
        buffer = self.cmb_fo2_mode.currentIndex() == 0
        for wdg in (self.cmb_fo2_mode, self.sp_dbuf, self.sp_dbuf_sig):
            wdg.setEnabled(fo2_used)
        self.cmb_buffer.setEnabled(fo2_used and buffer)
        p_used = "pressure" not in unused or (fo2_used and buffer)
        for wdg in (self.sp_P, self.sp_P_sig):
            wdg.setEnabled(p_used)

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
        plag_an = self.cmb_mineral.currentData() == "plagioclase" and self.an_values is not None
        self._activity_field.setVisible(plag_an)
        texts = []
        if equil:
            texts.append("The trace element starts in equilibrium with the measured anorthite "
                         "profile, C = C0 exp(A X_An / RT) (Dohmen et al. 2017, eq. A13). That is the "
                         "state diffusion ends in, so it is an initial state only if the crystal had "
                         "already equilibrated with one melt before the event being timed.")
        if plag_an:
            set_key = self.cmb_activity.currentData() or DEFAULT_ACTIVITY_SET
            texts.append(activity_note(set_key, self.sp_T.value() + 273.15) + " The two published "
                         "sets differ in sign for Mg (Dohmen et al. 2017, Table 1).")
        self.lbl_ic.setText(" ".join(texts))
        self.lbl_ic.setVisible(bool(texts))
        self._update_solver_note()

    def _update_geometry_note(self):
        self.lbl_geom.setText(Geometry(self.cmb_geom.currentText()).describe())
        self._update_solver_note()

    # --- boundaries and solver ------------------------------------------------------
    def _boundaries_far(self) -> bool:
        return self.cmb_bcl.currentData() == "far" and self.cmb_bcr.currentData() == "far"

    def _on_boundaries_changed(self):
        texts = {
            "far": "a plateau that continues past the traverse",
            "rim_melt": "a rim whose composition the melt holds fixed",
            "rim_closed": "a rim closed to exchange (zero flux)",
            "centre": "the crystal centre (zero flux by symmetry)",
        }
        l, r = self.cmb_bcl.currentData(), self.cmb_bcr.currentData()
        if l == r == "far":
            self.lbl_bc.setText("Both ends are plateaus that continue. This is the usual "
                                "error-function set-up, and the far ends do not affect the "
                                "time while both plateaus survive.")
        else:
            self.lbl_bc.setText(f"Left: {texts[l]}. Right: {texts[r]}.")
        self._update_solver_note()

    def _solver_short(self) -> str:
        try:
            ok, _ = self._model(self._checked_keys()[0]).can_use_analytical()
            return "closed form, Crank (1975) eq. 2.14" if ok else "numerical, Crank-Nicolson"
        except Exception:
            return "numerical, Crank-Nicolson" if self.cmb_solver.currentIndex() else "automatic"

    def _update_solver_note(self, *_):
        if not hasattr(self, "lbl_solver"):
            return
        ok, why = None, ""
        if self.profile is not None and hasattr(self, "lst_coef") and self.lst_coef.count():
            try:
                ok, why = self._model(self._checked_keys()[0]).can_use_analytical()
            except Exception:
                ok = None
        if ok is None:
            self.lbl_solver.setText("Diffusor picks the closed-form solution when it is exact "
                                    "and the numerical solver otherwise. Load data to see which.")
            self.sp_nodes.setEnabled(True)
            return
        if ok:
            self.lbl_solver.setText("<b>Closed form</b> for this run: " + why + ".")
        else:
            self.lbl_solver.setText("<b>Numerical</b> (Crank-Nicolson) for this run, because "
                                    + why + ".")
        self.sp_nodes.setEnabled(not ok)

    # --- results that no longer match the settings -------------------------------
    def _watch_inputs(self):
        """Note when any setting changes after a result was computed.

        A fit stays on screen until the next one, so without this a changed setting
        leaves a number that belongs to the old settings next to the new ones."""
        # these only configure a Monte Carlo run, not the fitted time
        skip = {self.sp_draws, self.sp_seed, self.sp_cores, self.chk_contrib}
        for kind, signal in ((QDoubleSpinBox, "valueChanged"), (QSpinBox, "valueChanged"),
                             (QComboBox, "currentIndexChanged"), (QCheckBox, "toggled"),
                             (QListWidget, "itemChanged")):
            for w in self.pages.findChildren(kind):
                if w not in skip:
                    getattr(w, signal).connect(self._inputs_changed)

    def _inputs_changed(self, *_):
        if self.fit_result is not None and not self._stale:
            self._stale = True
            self._rebuild_summary()

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
            self._update_condition_inputs()
        except Exception as exc:
            self.lbl_fo2.setText(str(exc))

    def _fit_profile(self):
        """The loaded profile without the points the user cut: what every fit uses."""
        return self.profile.kept()

    @Slot(object)
    def _on_points_cut(self, mask):
        """Points were cut from, or restored to, the fit on the plot."""
        p = self.profile
        if p is None:
            return
        if self._running():
            self._status("Wait for the running calculation or press Stop before cutting points.")
            return
        mask = np.asarray(mask, dtype=bool)
        if int(np.count_nonzero(~mask)) < 3:
            self._status("Keep at least three points in the fit.")
            return
        refit = self.fit_result is not None and self.mc_result is None
        p.set_excluded(mask)
        kept = p.kept()
        cut = p.excluded
        self.fit_result = self.mc_result = self.compare_results = None
        self._mc_draws, self._mc_meta = [], None
        self._stale = False
        self.plot.set_points(p.x, p.C, cut)
        self.plot.show_data(kept.x, kept.C, kept.sigma, y_label=self._y_label())
        self.preview.show_data(kept.x, kept.C, kept.sigma, y_label=self._y_label(),
                               cut=None if cut is None else (p.x[cut], p.C[cut]))
        if cut is None:
            self._log("no points cut: the whole profile is fitted")
        else:
            self._log(f"{p.n_excluded} of {len(p)} points cut from the fit, at x = "
                      + ", ".join(f"{v:.2f}" for v in p.x[cut][:12])
                      + (", ..." if p.n_excluded > 12 else "") + " um")
        self._rebuild_summary()
        if refit and self._ready():
            self.run_fit()

    def _guess_initial(self, from_button: bool = False):
        """Step and plateaus from the data, and the composition if the profile is one.

        Pressing **Guess** also replaces a typed-in composition; the guess made on
        loading leaves a typed-in or example value alone.
        """
        if self.profile is None:
            return
        kept = self._fit_profile()
        ic = guess_step_from_data(kept.x, kept.C)
        self.sp_x0.setValue(ic.params["x0"])
        self.sp_cl.setValue(ic.params["C_left"])
        self.sp_cr.setValue(ic.params["C_right"])
        self._log("initial profile: " + ic.describe())
        self._sync_composition(replace_explicit=from_button, report=True)

    # -- representative composition --------------------------------------------
    def _profile_as_composition(self) -> Optional[np.ndarray]:
        """The loaded profile in the law's composition coordinate, or None.

        Only for an exchange species is the profile the mineral's composition
        variable; olivine Fe-Mg goes through the chosen coordinate (X_Fe, X_Fo or
        Fo mol%). A TiO2 wt% or ppm profile is a concentration, and its mean put
        into x_Ti or X_An would be nonsense (5 wt% TiO2 became x_Ti = 1, and D 250
        times too large). Values outside 0-1 are refused for the same reason.
        """
        if self.profile is None or self.an_values is not None:
            return None
        mineral = get_mineral(self.cmb_mineral.currentData())
        species = self.cmb_species.currentText()
        if species not in EXCHANGE_SPECIES or mineral.composition_variable.key == "none":
            return None
        scale, offset = 1.0, 0.0
        if mineral.key == "olivine" and species == "Fe-Mg":
            scale, offset = self.cmb_ol_coordinate.currentData()
        X = offset + scale * np.asarray(self._fit_profile().C, dtype=float)
        if not (np.all(np.isfinite(X)) and 0.0 <= X.min() and X.max() <= 1.0):
            return None
        return X

    def _set_xcomp(self, value: float, source: str):
        self._xcomp_writing = True
        try:
            self.sp_xcomp.setValue(float(value))
        finally:
            self._xcomp_writing = False
        self._xcomp_source = source

    def _on_xcomp_changed(self, _value):
        if not self._xcomp_writing:
            self._xcomp_source = "explicit"
            self._xcomp_fallback = None

    def _release_composition(self):
        """Take back a value that came from the profile mean."""
        if self._xcomp_source == "profile" and self._xcomp_fallback is not None:
            self._set_xcomp(self._xcomp_fallback, "default")
        self._xcomp_source = "default"
        self._xcomp_fallback = None

    def _sync_composition(self, replace_explicit: bool = False, report: bool = False):
        """Keep the representative composition consistent with the profile.

        Runs on loading and whenever the mineral, species or olivine coordinate
        changes, because a profile is usually loaded before the mineral is chosen.
        """
        if self._xcomp_source == "explicit" and not replace_explicit:
            return
        X = self._profile_as_composition()
        if X is not None:
            if self._xcomp_source != "profile":
                self._xcomp_fallback = self.sp_xcomp.value()
            self._set_xcomp(float(np.mean(X)), "profile")
            if report:
                self._log(f"representative composition: {np.mean(X):.4g}, the profile mean")
            return
        if self._xcomp_source == "profile":
            self._release_composition()
        if report and self.profile is not None and self.an_values is None:
            self._log("representative composition left at "
                      f"{self.sp_xcomp.value():.4g}: the profile is not the host composition. "
                      "Set it on the Conditions step if the law needs one.")

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
        X = {}
        for k in coef.requires:
            X.setdefault(k, self.sp_xcomp.value())
        axis = {1: "a", 2: "b", 3: "c"}.get(self.cmb_axis.currentIndex())
        angles = ((self.sp_alpha.value(), self.sp_beta.value(), self.sp_gamma.value())
                  if self.cmb_axis.currentIndex() == 4 else None)
        return Conditions(T_K=T, P_Pa=P, log_fo2_bar=lf, X=X, axis=axis, angles_deg=angles)

    def _boundary(self, cmb, plateau: float):
        return zero_flux() if cmb.currentData() in ("rim_closed", "centre") else dirichlet(plateau)

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
                 "activity_set": self.cmb_activity.currentData() or DEFAULT_ACTIVITY_SET,
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
        bcl = self._boundary(self.cmb_bcl, c_left)
        bcr = self._boundary(self.cmb_bcr, c_right)
        hist = (ThermalHistory.linear(cond.T_K, self.sp_Tend.value() + 273.15, 1.0)
                if self.chk_cooling.isChecked() else None)
        exchange = species in ("Fe-Mg", "Fe-Ti", "NaSi-CaAl", "Na-K")
        comp_key = mineral.composition_variable.key if self.chk_comp_dep.isChecked() and exchange else None
        comp_scale, comp_offset = 1.0, 0.0
        if mineral.key == "olivine" and species == "Fe-Mg":
            comp_scale, comp_offset = self.cmb_ol_coordinate.currentData()
            comp_key = "XFe" if self.chk_comp_dep.isChecked() else None
        if comp_key and comp_key not in coef.requires:
            comp_key = None

        A_kJ = 0.0
        an_grid = None
        activity_set = self.cmb_activity.currentData() or DEFAULT_ACTIVITY_SET
        if mineral.key == "plagioclase" and species in ACTIVITY_SPECIES and self.an_values is not None:
            A_kJ = activity_A(species, cond.T_K, activity_set)
        n_nodes = self.sp_nodes.value()
        x_grid = None
        if self.profile is not None:
            x_grid = np.linspace(float(self.profile.x.min()), float(self.profile.x.max()), n_nodes)
            if A_kJ and self.an_values is not None:
                an_grid = np.interp(x_grid, self.profile.x, self.an_values)
        return DiffusionModel(
            coefficient=coef, conditions=cond, initial=ic,
            geometry=Geometry(self.cmb_geom.currentText()),
            bc_left=bcl, bc_right=bcr, history=hist,
            beam_sigma_um=self.sp_beam.value(), n_nodes=n_nodes, x_grid=x_grid,
            comp_key=comp_key, composition_dependent=bool(comp_key),
            comp_scale=comp_scale, comp_offset=comp_offset,
            an_profile=an_grid, activity_A_kJ=A_kJ if an_grid is not None else 0.0,
            activity_set=activity_set if an_grid is not None else "",
            force_numerical=self.cmb_solver.currentIndex() == 1,
            boundaries_far=self._boundaries_far(),
            fo2_buffer=((self.cmb_buffer.currentText(), self.sp_dbuf.value())
                        if self.cmb_fo2_mode.currentIndex() == 0 else None))

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
            self._status("Still working. Wait or press Stop.")
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
        for name in ("joint_study", "multicomponent_study"):
            study = getattr(self, name, None)
            if study is not None and study.thread is not None and study.thread.isRunning():
                study.abort()
                self._status("Cancelling the study. Close once its current calculation finishes.")
                event.ignore()
                return
        self.stop_work()
        for t, _ in self._jobs:
            t.quit()
            t.wait(3000)
        if self._update_job is not None:
            # a web request cannot be interrupted; it ends within its own timeout
            self._update_job[0].wait(int(updates.TIMEOUT_SECONDS * 1000) + 1000)
        super().closeEvent(event)

    def _ready(self) -> bool:
        if self.profile is None:
            QMessageBox.information(self, "No data", "Load a profile first.")
            return False
        if not self._checked_keys():
            QMessageBox.information(self, "No coefficient", "Choose a coefficient on step 3.")
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
        p = self._fit_profile()
        w = FitWorker(model, p.x, p.C, p.sigma, self._free_parameters(), T_MIN, T_MAX)
        w.finished.connect(self._fit_done)
        w.failed.connect(self._work_failed)
        if self._launch(w):
            self._status("Fitting...")

    @Slot(object)
    def _fit_done(self, res):
        self._busy(False)
        self.fit_result = res
        self.mc_result = None
        self.compare_results = None
        self._stale = False
        self.plot.show_fit(res, None, y_label=self._y_label())
        self._log(f"fit complete: {human_time(res.t_seconds)}\nsolver: {res.route}"
                  + "".join(f"\n{w}" for w in res.warnings))
        self._rebuild_summary()
        self._show_warnings(res.warnings)

    def run_compare(self):
        if not self._ready():
            return
        keys = self._checked_keys()
        if len(keys) < 2:
            QMessageBox.information(self, "Tick at least two",
                                    "Tick two or more coefficients on step 3.")
            return
        try:
            models = {get_coefficient(k).label: self._model(k) for k in keys}
        except Exception:
            QMessageBox.critical(self, "Model error", traceback.format_exc())
            return
        p = self._fit_profile()
        w = CompareWorker(models, p.x, p.C, p.sigma, self._free_parameters(), T_MIN, T_MAX)
        w.progress.connect(self._compare_progress)
        w.finished.connect(self._compare_done)
        w.failed.connect(self._work_failed)
        self._launch(w, len(models))

    @Slot(int, int, str)
    def _compare_progress(self, i, n, key):
        self.progress.setValue(i)
        self._status(f"{i}/{n}: {key}")

    @Slot(object)
    def _compare_done(self, results):
        self._busy(False)
        self.compare_results = results
        ok = {k: r for k, r in results.items() if not isinstance(r, Exception)}
        self.plot.show_comparison(ok, y_label=self._y_label())
        if ok:
            self.fit_result = next(iter(ok.values()))
            self._stale = False
        self._rebuild_summary()
        self._log("comparison complete\n" + "\n".join(
            f"{k}: {'failed' if isinstance(r, Exception) else human_time(r.t_seconds)}"
            for k, r in results.items()))
        self._compare_dialog = richtext.show(self, "Coefficient comparison",
                                             richtext.compare_html(results), 860, 480,
                                             modal=False)

    def run_mc(self):
        if not self._ready():
            return
        try:
            model = self._model(self._checked_keys()[0])
        except Exception:
            QMessageBox.critical(self, "Model error", traceback.format_exc())
            return
        p = self._fit_profile()
        n = self.sp_draws.value()
        budget = self._budget()
        w = MonteCarloWorker(model, p.x, p.C, p.sigma, budget, n,
                             self.sp_seed.value(), self._free_parameters(), T_MIN, T_MAX,
                             do_contributions=self.chk_contrib.isChecked(),
                             contribution_draws=max(40, n // 5),
                             workers=self.sp_cores.value())
        w.progress.connect(self._mc_progress)
        w.stage.connect(self._mc_stage)
        w.draws.connect(self._mc_draws_arrived)
        w.finished.connect(self._mc_done)
        w.failed.connect(self._work_failed)
        if self._launch(w, n):
            active = budget.active_sources()
            fixed = {}
            if "temperature" not in active:
                fixed["T"] = f"{self.sp_T.value():g} °C"
            if "fo2" not in active or not model.coefficient.needs_fo2:
                fixed["f"] = ("no fO2 term" if not model.coefficient.needs_fo2
                              else self._fo2_text())
            self._mc_draws = []
            self._mc_meta = dict(x=p.x, C=p.C, sigma=p.sigma, y_label=self._y_label(),
                                 n_total=n, sampled=active, fixed=fixed)
            self.plot.start_monte_carlo(**self._mc_meta)
            self._log(f"Monte Carlo: {n} draws, seed {self.sp_seed.value()}, "
                      f"up to {self.sp_cores.value()} cores, varying "
                      + ", ".join(SOURCE_NAMES.get(a, a) for a in active))

    def _fo2_text(self) -> str:
        if self.cmb_fo2_mode.currentIndex() == 0:
            return f"{self.cmb_buffer.currentText()} {self.sp_dbuf.value():+.2f}"
        return f"log fO2 {self.sp_dbuf.value():.2f}"

    @Slot(object)
    def _mc_draws_arrived(self, batch):
        self._mc_draws.extend(batch)
        if self._mc_meta is not None and getattr(self.plot, "_mc", None) is not None:
            self.plot.add_monte_carlo_draws(batch)

    @Slot(int, int)
    def _mc_progress(self, i, total):
        self.progress.setValue(i)
        self._status(f"Monte Carlo {i} of {total}")

    @Slot(str)
    def _mc_stage(self, stage):
        self._status(stage)

    @Slot(object)
    def _mc_done(self, res):
        self._busy(False)
        self.mc_result = res
        # The run refits the data it was given, so its best fit is the one for these
        # settings. An earlier Fit may have been made with other settings, and the
        # headline, the plot legend and the histograms must all show the same time.
        if res.base_fit is not None:
            self.fit_result = res.base_fit
        elif self.fit_result is None:
            from ..fitting import fit_time
            p = self._fit_profile()
            self.fit_result = fit_time(self._model(self._checked_keys()[0]), p.x, p.C, p.sigma,
                                       self._free_parameters())
        self._stale = False
        self.plot.finish_monte_carlo(res, self._best_curve(res))
        self._rebuild_summary()
        self._log(f"Monte Carlo complete: median {human_time(res.median)}, 68% "
                  f"{human_time(res.p16)} to {human_time(res.p84)}, "
                  f"{res.workers} {'core' if res.workers == 1 else 'cores'}")
        self._status("Monte Carlo done. Details shows the numbers.")

    def _best_curve(self, res):
        try:
            x = self.profile.x
            xf = np.linspace(float(x.min()), float(x.max()), 300)
            return xf, self.fit_result.model.profile(res.t_best, xf)
        except Exception:
            return None

    @Slot(str)
    def _work_failed(self, msg: str):
        self._busy(False)
        self._log(msg, "error")
        QMessageBox.critical(self, "Calculation failed", msg)

    def _show_warnings(self, warnings):
        serious = [w for w in warnings
                   if "SUPERSEDED" in w or "NOT verified" in w or "semi-infinite" in w
                   or "resolution limit" in w]
        if serious:
            QMessageBox.warning(self, "Check before using the result", "\n\n".join(serious))

    # ================================================================ views and dialogs
    def show_profile_view(self):
        if self.compare_results and self.fit_result is None:
            return
        if self.fit_result is not None:
            self.plot.show_fit(self.fit_result, self.mc_result, y_label=self._y_label())
        elif self.profile is not None:
            p = self._fit_profile()
            self.plot.show_data(p.x, p.C, p.sigma, y_label=self._y_label())

    def show_histogram(self):
        if self._mc_meta is None or not self._mc_draws:
            QMessageBox.information(self, "No Monte Carlo yet", "Run a Monte Carlo first.")
            return
        if self.step != RESULTS:
            self._go(RESULTS)
        self.plot.start_monte_carlo(**self._mc_meta)
        self.plot.add_monte_carlo_draws(self._mc_draws)
        if self.mc_result is not None:
            self.plot.finish_monte_carlo(self.mc_result, self._best_curve(self.mc_result))

    def show_run_details(self):
        if self.mc_result is not None:
            richtext.show(self, "Monte Carlo result", richtext.mc_html(self.mc_result), 760, 620,
                          modal=False)
        elif self.compare_results:
            richtext.show(self, "Coefficient comparison",
                          richtext.compare_html(self.compare_results), 860, 480, modal=False)
        elif self.fit_result is not None:
            self.show_methods()
        else:
            QMessageBox.information(self, "Nothing yet", "Run a fit first.")

    def show_methods(self):
        if self.fit_result is None:
            QMessageBox.information(self, "Nothing yet", "Run a fit first.")
            return
        richtext.show(self, "Methods and references",
                      richtext.methods_html(self.fit_result, self.mc_result, self._fit_profile()),
                      900, 720)

    def show_log(self):
        richtext.show(self, "Log", richtext.log_html(self._log_entries), 860, 560)

    def show_coefficient_info(self):
        it = self.lst_coef.currentItem()
        if it is None:
            return
        c = get_coefficient(it.data(Qt.UserRole))
        richtext.show(self, c.label, richtext.coefficient_html(c), 860, 700)

    def show_multicomponent_study(self):
        from .multicomponent_study import MulticomponentStudyDialog
        if not hasattr(self, "multicomponent_study"):
            self.multicomponent_study = MulticomponentStudyDialog(self)
        self.multicomponent_study.show()
        self.multicomponent_study.raise_()

    def show_joint_study(self):
        from .joint_study import JointStudyDialog
        if not hasattr(self, "joint_study"):
            self.joint_study = JointStudyDialog(self)
        self.joint_study.show()
        self.joint_study.raise_()

    def export_results(self):
        if self.fit_result is None:
            QMessageBox.information(self, "Nothing to export", "Run a fit first.")
            return
        d = QFileDialog.getExistingDirectory(self, "Export into folder", self.output_dir())
        if not d:
            return
        try:
            written = save_results(d, self.fit_result, self.mc_result, self._fit_profile(),
                                   figure=self.plot.figure)
            self._log("exported\n" + "\n".join(f"{k}: {v}" for k, v in written.items()))
            QMessageBox.information(self, "Exported", "Written:\n" + "\n".join(written.values()))
        except Exception:
            QMessageBox.critical(self, "Export failed", traceback.format_exc())

    def export_figure(self):
        """The plot as an editable vector file: text as text, one object per data point."""
        from ..dataio.vector import save_vector
        source = self.profile.source if self.profile is not None else None
        name = f"{Path(source).stem if source else 'profile'}_diffusor"
        start = str(Path(self.output_dir()) / name) if self.output_dir() else name
        path, _ = QFileDialog.getSaveFileName(
            self, "Export figure for editing", start + ".svg",
            "SVG, for Inkscape and CorelDRAW (*.svg);;PDF (*.pdf);;EPS (*.eps)")
        if not path:
            return
        try:
            save_vector(self.plot.figure, path, bbox_inches="tight")
            self._log("figure exported\n" + path)
            QMessageBox.information(
                self, "Figure exported",
                path + "\n\nText is editable text and every data point is its own object, "
                "named for its series and number (for example measured_point_007), "
                "with its error bar in the same group.")
        except Exception:
            QMessageBox.critical(self, "Export failed", traceback.format_exc())

    def output_dir(self) -> str:
        """The loaded profile's folder: where save and export dialogs open."""
        if self.profile is not None and self.profile.source:
            folder = Path(self.profile.source).parent
            if folder.is_dir():
                return str(folder)
        return ""

    # ================================================================ updates
    # The check asks GitHub for the newest release and, if it is newer than this
    # one, shows the banner. It never downloads or changes anything itself.
    @staticmethod
    def _settings() -> QSettings:
        # DIFFUSOR_SETTINGS names an .ini file to use instead, which keeps the test
        # suite away from the user's own settings and recent files
        path = os.environ.get("DIFFUSOR_SETTINGS")
        if path:
            return QSettings(path, QSettings.IniFormat)
        return QSettings("Diffusor", "Diffusor")

    def startup_update_check(self):
        """Check quietly, if the user allows it and none succeeded in the last day."""
        s = self._settings()
        if not s.value("updates/check_at_startup", True, type=bool):
            return
        if not updates.due(str(s.value("updates/last_check", "") or "")):
            return
        self.check_for_updates(manual=False)

    def check_for_updates(self, manual: bool = False):
        if self._update_job is not None and self._update_job[0].isRunning():
            if manual:
                self._status("Already checking for updates.")
            return
        from .. import __version__
        self._update_manual = manual
        worker = UpdateWorker(__version__)
        worker.finished.connect(self._update_found)
        worker.failed.connect(self._update_failed)
        self._update_job = (start(worker), worker)
        if manual:
            self._status("Checking for updates...")

    @Slot(object)
    def _update_found(self, info):
        from .. import __version__
        self._settings().setValue("updates/last_check",
                                  datetime.now(timezone.utc).isoformat())
        manual, self._update_manual = self._update_manual, False
        if updates.is_newer(info.version, __version__):
            self._update_release = info
            self.lbl_update.setText(f"Diffusor {info.version} is available. "
                                    f"This is {__version__}.")
            self.lbl_update.setToolTip(info.notes[:800])
            self.update_banner.setVisible(True)
            self._log(f"update available: {info.version}")
            if manual:
                self._status("")
        elif manual:
            self._status("")
            QMessageBox.information(self, "Check for updates",
                                    f"Diffusor {__version__} is the newest version.")

    @Slot(str)
    def _update_failed(self, message: str):
        manual, self._update_manual = self._update_manual, False
        self._log(f"update check: {message}")
        if manual:
            self._status("")
            QMessageBox.information(self, "Check for updates", message)

    def _open_update_page(self):
        if self._update_release is not None:
            QDesktopServices.openUrl(QUrl(self._update_release.url))
        self.update_banner.setVisible(False)

    def show_about(self):
        from .. import __version__
        QMessageBox.about(self, "About Diffusor", (
            f"<b>Diffusor {__version__}</b><br><br>"
            "Diffusion chronometry with closed-form and Crank-Nicolson solvers, a registry "
            "of published diffusion coefficients, and Monte Carlo error propagation.<br><br>"
            "Every equation and coefficient carries its citation. View > Methods lists the "
            "references used by the current run, and Help > All references lists every "
            "source in the app."))
