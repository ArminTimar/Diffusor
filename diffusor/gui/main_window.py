"""Diffusor's main window: a stepped flow ending in a results view.

The settings are collected one group at a time, then summarised in a narrow
sidebar next to the plot. Any group in that summary can be clicked to go back
and change it. Pages never scroll. Only lists and reading panes do.
"""
from __future__ import annotations

import traceback
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from PySide6.QtCore import QLocale, Qt, QTimer, Slot
from PySide6.QtGui import QAction, QFont, QGuiApplication
from PySide6.QtWidgets import (QApplication, QCheckBox, QComboBox, QDialog, QDialogButtonBox,
                               QDoubleSpinBox, QFileDialog, QFrame, QHBoxLayout, QHeaderView,
                               QLabel, QLineEdit, QListWidget, QListWidgetItem, QMainWindow,
                               QMessageBox, QProgressBar, QPushButton, QSizePolicy, QSpinBox,
                               QSplitter, QStackedWidget, QTableWidget, QTableWidgetItem,
                               QVBoxLayout, QWidget)

from .. import datasets as ds
from ..coefficients import Conditions, get as get_coefficient, list_coefficients
from ..coefficients.plagioclase import ACTIVITY_A, activity_theta
from ..dataio import ProfileSpec, build_profile, read_table, save_results, suggest_spec
from ..dataio.image_profiles import is_extraction_workbook, read_extraction
from ..dataio.images import PILLOW_SUFFIXES, RAW_SUFFIXES
from ..fitting import DiffusionModel, UncertaintyBudget
from ..fitting.fit import T_MAX_DEFAULT as T_MAX, T_MIN_DEFAULT as T_MIN
from ..minerals import MINERALS, get_mineral
from ..references import cite
from ..solvers import Geometry, InitialCondition, dirichlet, zero_flux
from ..solvers.history import ThermalHistory
from ..solvers.initial import guess_step_from_data
from ..thermo import available_buffers, log_fo2_from_delta
from ..thermo.units import human_time
from . import format_help, richtext, theme
from .plot_widget import DataPreview, ProfilePlot
from .widgets import (WrapLabel, callout, card, collapsible, divider, field, ghost_button,
                      install_no_wheel, note, page_columns, pair, primary_button, row,
                      scrollable)
from .workers import CompareWorker, FitWorker, MonteCarloWorker, start

STEPS = ["Data", "Mineral", "Conditions", "Model", "Coefficient", "Uncertainty", "Results"]
RESULTS = len(STEPS) - 1

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
        self.fit_result = None
        self.mc_result = None
        self.compare_results = None
        self._jobs: List = []           # (thread, worker) pairs kept alive until finished
        self._log_entries: List[Tuple[datetime, str, str]] = []
        self._prefill: Dict[int, QLabel] = {}
        self._mc_draws: List[dict] = []
        self._mc_meta: Optional[dict] = None
        self._loaded_key: Optional[str] = None
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

        self._build_menu()
        self._on_mineral_changed()
        self._on_resolution_changed()
        self._go(0)

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
                         ("Shared-duration study...", self.show_joint_study),
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
                         ("About", self.show_about)):
            a = QAction(text, self); a.triggered.connect(fn); h.addAction(a)

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

        data_card, dbody = card("Loaded data")
        self.lbl_data = note("Nothing loaded yet.", "Hint")
        dbody.addWidget(self.lbl_data)
        self.preview = DataPreview()
        self.preview.setMinimumHeight(180)
        dbody.addWidget(self.preview, 1)

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
        return page_columns([load_card, (data_card, 1)], [(ex_card, 1)])

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
        if self.step != 0 and self.dataset is not None:
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

        x, xbody = card("Composition")
        self.sp_xcomp = _spin(0, 1, 0.15, 4, 0.01)
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
        xbody.addWidget(self.cmb_ol_coordinate)
        self.lbl_comp = note("", "Hint")
        xbody.addWidget(self.lbl_comp)
        return page_columns([c, x], [o], top=[self._prefill_bar(1)])

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
        self.box_unused = callout("Info")
        self.lbl_unused = self.box_unused.label
        self.box_unused.setVisible(False)

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
        return page_columns([c, f], [r, self.box_unused],
                            top=[self._prefill_bar(2, with_sources=True)])

    # ================================================================ step 4
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
        bg.clicked.connect(self._guess_initial)
        self.chk_free_x0 = QCheckBox("Fit the step position")
        self.chk_free_x0.setChecked(True)
        self.chk_free_plateaus = QCheckBox("Fit the plateaus")
        ibody.addWidget(row(bg, self.chk_free_x0, self.chk_free_plateaus, spacing=14))

        n, nbody = card("Solver and thermal history")
        self.cmb_solver = QComboBox()
        self.cmb_solver.addItems(["Automatic", "Always numerical"])
        self.cmb_solver.currentIndexChanged.connect(self._update_solver_note)
        self.sp_nodes = QSpinBox(); self.sp_nodes.setRange(51, 4001); self.sp_nodes.setValue(301)
        self.sp_nodes.setLocale(number_locale())
        nbody.addWidget(pair(field("Solver", self.cmb_solver), field("Grid nodes", self.sp_nodes)))
        box = callout("Info")
        self.lbl_solver = box.label
        nbody.addWidget(box)
        self.chk_cooling = QCheckBox("Linear cooling to")
        self.sp_Tend = _spin(300, 2000, 900, 1, 5, " °C")
        self.sp_Tend.setEnabled(False)
        self.chk_cooling.toggled.connect(self.sp_Tend.setEnabled)
        self.chk_cooling.toggled.connect(self._update_solver_note)
        nbody.addWidget(row(self.chk_cooling, self.sp_Tend, spacing=10))
        self._update_geometry_note()
        self._on_boundaries_changed()
        return page_columns([g, n], [i, b], top=[self._prefill_bar(3)])

    # ================================================================ step 5
    def _page_coefficient(self) -> QWidget:
        c, body = card("Diffusion coefficient")
        body.addWidget(note("Tick one to fit. Tick several to compare them.", "Hint"))
        self.lst_coef = QListWidget()
        self.lst_coef.setMouseTracking(True)
        self.lst_coef.currentItemChanged.connect(self._coefficient_selected)
        self.lst_coef.itemChanged.connect(lambda _it: self._on_coefficients_ticked())
        body.addWidget(self.lst_coef, 1)

        d, dbody = card("Selected")
        self.lbl_coef_name = note("", "H2")
        dbody.addWidget(self.lbl_coef_name)
        self.lbl_coef_eq = note("", "Equation")
        # set in code, so the wrapped height is measured with the font actually drawn
        self.lbl_coef_eq.setFont(QFont("Consolas", 9))
        dbody.addWidget(self.lbl_coef_eq)
        self.lbl_coef = note("", "Hint")
        dbody.addWidget(self.lbl_coef)
        self.lbl_coef_range = note("", "Hint")
        dbody.addWidget(self.lbl_coef_range)
        b = ghost_button("Full details and references")
        b.clicked.connect(lambda: self.show_coefficient_info())
        dbody.addWidget(row(b))
        return page_columns([(c, 1)], [d], top=[self._prefill_bar(4)])

    # ================================================================ step 6
    def _page_uncertainty(self) -> QWidget:
        c, body = card("Monte Carlo")
        self.sp_draws = QSpinBox(); self.sp_draws.setRange(20, 100000); self.sp_draws.setValue(500)
        self.sp_seed = QSpinBox(); self.sp_seed.setRange(0, 2 ** 31 - 1); self.sp_seed.setValue(12345)
        self.sp_draws.setLocale(number_locale())
        self.sp_seed.setLocale(number_locale())
        body.addWidget(pair(field("Draws", self.sp_draws), field("Random seed", self.sp_seed)))
        body.addWidget(note("500 for a quick look, 1000 for a paper.", "Hint"))

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
                          "independent": "each parameter independently, because the paper "
                                         "gives neither a covariance nor a scatter",
                          "none": "fixed coefficient parameters; coefficient uncertainty is not quantified"}[best]
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
        self._summary_line("ends",
                           f"{self.cmb_bcl.currentText().lower()} / "
                           f"{self.cmb_bcr.currentText().lower()}")
        self._summary_line("solver", self._solver_short())

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
                           ", ".join(SOURCE_NAMES.get(s, s) for s in b.active_sources())
                           or "nothing sampled")

        self.summary_layout.addStretch(1)

    # ================================================================ navigation
    def _go(self, index: int):
        height_before = self.height()
        index = max(0, min(index, RESULTS))
        if index > 0 and self.profile is None:
            self._status("Load a profile first")
            index = 0
        self.step = index
        self.pages.setCurrentIndex(index)
        for i, lab in enumerate(self.step_labels):
            _repolish(lab, "StepDotActive" if i == index
                      else ("StepDotDone" if i < index else "StepDot"))
        self.btn_back.setVisible(index > 0)
        if index == 2:
            self._update_unused_note()
        if index == 3:
            self._update_solver_note()
        if index == 5:
            self._on_dmode_changed()
        if index == RESULTS:
            self.btn_next.setVisible(False)
            self.lbl_footer.setText("Click edit in the summary to change a setting.")
            self._rebuild_summary()
        else:
            self.btn_next.setVisible(True)
            self.btn_next.setText("Review and fit" if index == RESULTS - 1 else "Continue")
            self.lbl_footer.setText(f"Step {index + 1} of {RESULTS}")
        self.btn_next.setEnabled(self.profile is not None or index == 0)
        if self.isVisible():
            QTimer.singleShot(60, lambda h=height_before: self._restore_height(h))

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
        if self.step == 0 and self.profile is None:
            QMessageBox.information(self, "No data", "Load a file or an example first.")
            return
        self._go(self.step + 1)

    # ================================================================ data
    def load_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Load profile", "",
            "Tables and image profiles (*.csv *.txt *.tsv *.xlsx *.xls);;All files (*)")
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

    def _load(self, path: Path, dataset: Optional[ds.ExampleDataset]):
        if dataset is None and is_extraction_workbook(path):
            try:
                table = read_extraction(path)
            except Exception:
                QMessageBox.critical(self, "Could not load the file", traceback.format_exc())
                return
            self.load_image_table(table)
            return
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
        except Exception:
            self._log(traceback.format_exc(), "error")
            QMessageBox.critical(self, "Could not load the file", traceback.format_exc())
            return
        self._use_table(df, spec, path, dataset, an)

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

    def _use_table(self, df, spec, path: Path, dataset: Optional[ds.ExampleDataset], an):
        try:
            self.profile = build_profile(df, spec, source=str(path))
            self.dataset = dataset
            self.an_values = None
            for lab in self._prefill.values():
                lab.setVisible(False)
            self.btn_sources.setVisible(False)
            p = self.profile
            if dataset is not None:
                self._apply_dataset_settings(dataset, df)
                msg = (f"<b>{dataset.name.split(' (')[0]}</b>, {len(p)} points from "
                       f"{p.x.min():.1f} to {p.x.max():.1f} um. The example also filled in the "
                       "mineral, conditions, model and coefficient. Each step shows what it set.")
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
            self._log(f"loaded {path}\n{p.spec.describe()}"
                      + "".join(f"\nnote: {n}" for n in p.notes))
            self._guess_initial()
            self.preview.show_data(p.x, p.C, p.sigma, y_label=self._y_label())
            self.plot.show_data(p.x, p.C, p.sigma, y_label=self._y_label())
            self.fit_result = self.mc_result = self.compare_results = None
            self._mc_draws = []
            self.btn_next.setEnabled(True)
            self._status(f"Loaded {len(p)} points")
        except Exception:
            self._log(traceback.format_exc(), "error")
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
        fo2 += f" ± {s.get('sigma_delta', 0.3):g}"
        axis = s.get("axis")
        self.cmb_axis.setCurrentIndex({"a": 1, "b": 2, "c": 3}.get(axis, 0))
        self.cmb_geom.setCurrentText(s.get("geometry", "plane"))
        for cmb, key in ((self.cmb_bcl, "bc_left"), (self.cmb_bcr, "bc_right")):
            i = cmb.findData(s.get(key, "far"))
            cmb.setCurrentIndex(i if i >= 0 else 0)
        self.chk_comp_dep.setChecked(bool(s.get("composition_dependent", False)))
        self.cmb_ol_coordinate.setCurrentIndex(s.get("olivine_coordinate", 0))
        if "x_composition" in s:
            self.sp_xcomp.setValue(s["x_composition"])
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
            self.sp_xcomp.setValue(float(np.nanmean(vals)))
            self._log(f"anorthite read from '{an_col}' "
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
        self.btn_sources.setVisible(True)
        if self._boundaries_far():
            ends = "plateaus continuing past both ends"
        else:
            ends = (f"left end {self.cmb_bcl.currentText().lower()}, right end "
                    f"{self.cmb_bcr.currentText().lower()}")
        self._show_prefill(3, head + f"{s.get('geometry', 'plane')} geometry, "
                           f"{self.cmb_ic.currentText().lower()}, {ends}.")
        if want:
            self._show_prefill(4, head + get_coefficient(want).label + ".")

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
        self.lbl_coef_eq.setText(c.equation_text)
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
        text = richtext.unused_conditions(get_coefficient(keys[0])) if keys else ""
        self.lbl_unused.setText(text)
        self.box_unused.setVisible(bool(text))

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
            self.lbl_ic.setText("")
        self.lbl_ic.setVisible(equil)
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
            comp_scale=comp_scale, comp_offset=comp_offset,
            an_profile=an_grid, activity_theta=theta if an_grid is not None else 0.0,
            force_numerical=self.cmb_solver.currentIndex() == 1,
            boundaries_far=self._boundaries_far())

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
        if hasattr(self, "joint_study") and self.joint_study.thread is not None and self.joint_study.thread.isRunning():
            self.joint_study.abort()
            self._status("Cancelling the joint study. Close once its current calculation finishes.")
            event.ignore()
            return
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
                                    "Tick two or more coefficients on step 5.")
            return
        try:
            models = {get_coefficient(k).label: self._model(k) for k in keys}
        except Exception:
            QMessageBox.critical(self, "Model error", traceback.format_exc())
            return
        p = self.profile
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
        p = self.profile
        n = self.sp_draws.value()
        budget = self._budget()
        w = MonteCarloWorker(model, p.x, p.C, p.sigma, budget, n,
                             self.sp_seed.value(), self._free_parameters(), T_MIN, T_MAX,
                             do_contributions=self.chk_contrib.isChecked(),
                             contribution_draws=max(40, n // 5))
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
            self._log(f"Monte Carlo: {n} draws, seed {self.sp_seed.value()}, varying "
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
        if self.fit_result is None:
            from ..fitting import fit_time
            self.fit_result = fit_time(self._model(self._checked_keys()[0]), self.profile.x,
                                       self.profile.C, self.profile.sigma,
                                       self._free_parameters())
        self.plot.finish_monte_carlo(res, self._best_curve(res))
        self._rebuild_summary()
        self._log(f"Monte Carlo complete: median {human_time(res.median)}, 68% "
                  f"{human_time(res.p16)} to {human_time(res.p84)}")
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
            self.plot.show_data(self.profile.x, self.profile.C, self.profile.sigma,
                                y_label=self._y_label())

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
                      richtext.methods_html(self.fit_result, self.mc_result, self.profile),
                      900, 720)

    def show_log(self):
        richtext.show(self, "Log", richtext.log_html(self._log_entries), 860, 560)

    def show_coefficient_info(self):
        it = self.lst_coef.currentItem()
        if it is None:
            return
        c = get_coefficient(it.data(Qt.UserRole))
        richtext.show(self, c.label, richtext.coefficient_html(c), 860, 700)

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
        d = QFileDialog.getExistingDirectory(self, "Export into folder")
        if not d:
            return
        try:
            written = save_results(d, self.fit_result, self.mc_result, self.profile,
                                   figure=self.plot.figure)
            self._log("exported\n" + "\n".join(f"{k}: {v}" for k, v in written.items()))
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
