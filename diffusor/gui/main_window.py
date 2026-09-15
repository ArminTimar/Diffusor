"""The Diffusor main window."""
from __future__ import annotations

import copy
import traceback
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtGui import QAction, QFont
from PySide6.QtWidgets import (QCheckBox, QComboBox, QDialog, QDialogButtonBox,
                               QDoubleSpinBox, QFileDialog, QFormLayout, QGroupBox,
                               QHBoxLayout, QLabel, QLineEdit, QListWidget,
                               QListWidgetItem, QMainWindow, QMessageBox, QProgressBar,
                               QPushButton, QScrollArea, QSpinBox, QSplitter,
                               QTabWidget, QTextEdit, QVBoxLayout, QWidget)

from ..coefficients import Conditions, list_coefficients, get as get_coefficient
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
from .plot_widget import ProfilePlot
from .workers import CompareWorker, FitWorker, MonteCarloWorker, start

MONO = "Consolas, Menlo, monospace"


def _spin(lo, hi, val, dec=3, step=1.0, suffix=""):
    s = QDoubleSpinBox()
    s.setRange(lo, hi)
    s.setDecimals(dec)
    s.setSingleStep(step)
    s.setValue(val)
    if suffix:
        s.setSuffix(suffix)
    return s


class ColumnDialog(QDialog):
    """Map the columns of a loaded table onto the model."""

    def __init__(self, df, spec: ProfileSpec, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Map columns")
        self.df = df
        cols = [""] + [str(c) for c in df.columns]
        form = QFormLayout(self)

        self.dist = QComboBox(); self.dist.addItems(cols)
        self.unit = QComboBox(); self.unit.addItems(["um", "mm", "nm", "m"])
        self.a = QComboBox(); self.a.addItems(cols)
        self.b = QComboBox(); self.b.addItems(cols)
        self.sa = QComboBox(); self.sa.addItems(cols)
        self.sb = QComboBox(); self.sb.addItems(cols)
        self.mode = QComboBox(); self.mode.addItems(["A/(A+B)", "B/(A+B)", "A", "A-B"])
        self.oxa = QLineEdit(); self.oxb = QLineEdit()
        self.oxa.setPlaceholderText("e.g. FeO (leave blank if not wt% oxide)")
        self.oxb.setPlaceholderText("e.g. MgO")
        self.slevel = QComboBox(); self.slevel.addItems(["1s", "2s"])

        def setc(combo, val):
            combo.setCurrentText(str(val) if val is not None else "")

        setc(self.dist, spec.distance_column); setc(self.a, spec.column_a)
        setc(self.b, spec.column_b); setc(self.sa, spec.sigma_a_column)
        setc(self.sb, spec.sigma_b_column)
        self.mode.setCurrentText(spec.mode)

        form.addRow("Distance column", self.dist)
        form.addRow("Distance unit", self.unit)
        form.addRow("Element A", self.a)
        form.addRow("Element B", self.b)
        form.addRow("Uncertainty of A", self.sa)
        form.addRow("Uncertainty of B", self.sb)
        form.addRow("Uncertainty level", self.slevel)
        form.addRow("Modelled variable", self.mode)
        form.addRow("Oxide of A", self.oxa)
        form.addRow("Oxide of B", self.oxb)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept); bb.rejected.connect(self.reject)
        form.addRow(bb)

    def spec(self) -> ProfileSpec:
        def g(c):
            t = c.currentText().strip()
            return t or None
        return ProfileSpec(
            distance_column=g(self.dist), column_a=g(self.a), column_b=g(self.b),
            sigma_a_column=g(self.sa), sigma_b_column=g(self.sb),
            distance_unit=self.unit.currentText(), mode=self.mode.currentText(),
            oxide_a=self.oxa.text().strip() or None, oxide_b=self.oxb.text().strip() or None,
            sigma_level=self.slevel.currentText())


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Diffusor - diffusion chronometry")
        self.resize(1500, 950)
        self.profile = None
        self.fit_result = None
        self.mc_result = None
        self.compare_results = None
        self._threads = []

        self.plot = ProfilePlot()
        self.log = QTextEdit(); self.log.setReadOnly(True)
        self.log.setFont(QFont(MONO.split(",")[0], 9))
        self.results = QTextEdit(); self.results.setReadOnly(True)
        self.results.setFont(QFont(MONO.split(",")[0], 9))
        self.methods = QTextEdit(); self.methods.setReadOnly(True)
        self.methods.setFont(QFont(MONO.split(",")[0], 9))

        right = QTabWidget()
        right.addTab(self.results, "Results")
        right.addTab(self.methods, "Methods and references")
        right.addTab(self.log, "Log")

        centre = QSplitter(Qt.Horizontal)
        centre.addWidget(self.plot)
        centre.addWidget(right)
        centre.setSizes([950, 550])

        outer = QSplitter(Qt.Horizontal)
        outer.addWidget(self._build_controls())
        outer.addWidget(centre)
        outer.setSizes([430, 1070])
        self.setCentralWidget(outer)

        self.progress = QProgressBar(); self.progress.setVisible(False)
        self.statusBar().addPermanentWidget(self.progress)
        self.statusBar().showMessage("Load a profile to begin")
        self._build_menu()
        self._on_mineral_changed()

    # ------------------------------------------------------------------
    def _build_menu(self):
        m = self.menuBar().addMenu("&File")
        a = QAction("Load profile...", self); a.triggered.connect(self.load_profile); m.addAction(a)
        a = QAction("Export results...", self); a.triggered.connect(self.export_results); m.addAction(a)
        m.addSeparator()
        a = QAction("Quit", self); a.triggered.connect(self.close); m.addAction(a)
        h = self.menuBar().addMenu("&Help")
        a = QAction("Coefficient details", self); a.triggered.connect(self.show_coefficient_info); h.addAction(a)
        a = QAction("About", self); a.triggered.connect(self.show_about); h.addAction(a)

    def _build_controls(self) -> QWidget:
        panel = QWidget()
        v = QVBoxLayout(panel)
        v.setContentsMargins(6, 6, 6, 6)

        # --- data ---------------------------------------------------------
        g = QGroupBox("1. Data"); f = QFormLayout(g)
        self.btn_load = QPushButton("Load profile (CSV / Excel)...")
        self.btn_load.clicked.connect(self.load_profile)
        f.addRow(self.btn_load)
        self.lbl_data = QLabel("no data loaded"); self.lbl_data.setWordWrap(True)
        f.addRow(self.lbl_data)
        self.chk_greyscale = QCheckBox("column A is a BSE grey value")
        self.chk_greyscale.setToolTip(
            "Grey values must be calibrated against microprobe anchor points before fitting "
            "(Petrone et al. 2016). Use File > Load profile for the anchors.")
        f.addRow(self.chk_greyscale)
        v.addWidget(g)

        # --- mineral ------------------------------------------------------
        g = QGroupBox("2. Mineral and species"); f = QFormLayout(g)
        self.cmb_mineral = QComboBox()
        for k, mn in MINERALS.items():
            self.cmb_mineral.addItem(mn.name, k)
        self.cmb_mineral.currentIndexChanged.connect(self._on_mineral_changed)
        f.addRow("Mineral", self.cmb_mineral)
        self.cmb_species = QComboBox()
        self.cmb_species.currentIndexChanged.connect(self._refresh_coefficients)
        f.addRow("Species", self.cmb_species)
        self.cmb_axis = QComboBox()
        self.cmb_axis.addItems(["reference axis of the paper", "a", "b", "c", "angles to a,b,c"])
        self.cmb_axis.currentIndexChanged.connect(self._on_axis_changed)
        f.addRow("Traverse", self.cmb_axis)
        row = QWidget(); rl = QHBoxLayout(row); rl.setContentsMargins(0, 0, 0, 0)
        self.sp_alpha = _spin(0, 180, 90, 1, 1, " deg")
        self.sp_beta = _spin(0, 180, 90, 1, 1, " deg")
        self.sp_gamma = _spin(0, 180, 0, 1, 1, " deg")
        for s in (self.sp_alpha, self.sp_beta, self.sp_gamma):
            s.setEnabled(False); rl.addWidget(s)
        f.addRow("alpha/beta/gamma", row)
        self.sp_xcomp = _spin(0, 1, 0.15, 4, 0.01)
        f.addRow("Composition for D (X)", self.sp_xcomp)
        self.chk_comp_dep = QCheckBox("profile itself sets X (composition-dependent D)")
        self.chk_comp_dep.setChecked(True)
        f.addRow(self.chk_comp_dep)
        v.addWidget(g)

        # --- conditions ----------------------------------------------------
        g = QGroupBox("3. Conditions"); f = QFormLayout(g)
        self.sp_T = _spin(300, 2000, 950, 1, 5, " C")
        self.sp_T_sig = _spin(0, 300, 20, 1, 1, " K")
        f.addRow("Temperature", self.sp_T); f.addRow("  1 sigma", self.sp_T_sig)
        self.sp_P = _spin(0, 5000, 200, 1, 10, " MPa")
        self.sp_P_sig = _spin(0, 2000, 100, 1, 10, " MPa")
        f.addRow("Pressure", self.sp_P); f.addRow("  1 sigma", self.sp_P_sig)
        self.cmb_fo2_mode = QComboBox(); self.cmb_fo2_mode.addItems(["buffer offset", "absolute log fO2"])
        self.cmb_fo2_mode.currentIndexChanged.connect(self._update_fo2_label)
        f.addRow("fO2 given as", self.cmb_fo2_mode)
        self.cmb_buffer = QComboBox(); self.cmb_buffer.addItems(available_buffers())
        self.cmb_buffer.setCurrentText("NNO")
        self.cmb_buffer.currentIndexChanged.connect(self._update_fo2_label)
        f.addRow("Buffer", self.cmb_buffer)
        self.sp_dbuf = _spin(-8, 8, 1.0, 2, 0.1)
        self.sp_dbuf.valueChanged.connect(self._update_fo2_label)
        self.sp_dbuf_sig = _spin(0, 5, 0.3, 2, 0.1)
        f.addRow("Offset / log fO2 (bar)", self.sp_dbuf); f.addRow("  1 sigma", self.sp_dbuf_sig)
        self.lbl_fo2 = QLabel(""); self.lbl_fo2.setWordWrap(True)
        f.addRow(self.lbl_fo2)
        self.sp_beam = _spin(0, 20, 0.0, 2, 0.1, " um")
        self.sp_beam.setToolTip("Gaussian sigma of the analytical spatial resolution; 0 disables "
                                "the convolution correction (Ganguly et al. 1988).")
        f.addRow("Beam sigma", self.sp_beam)
        self.sp_xscale_sig = _spin(0, 0.5, 0.0, 3, 0.005)
        self.sp_xscale_sig.setToolTip("Relative 1 sigma of the distance calibration, e.g. 0.02 "
                                      "for a 2% image scale uncertainty.")
        f.addRow("Distance scale 1 sigma", self.sp_xscale_sig)
        v.addWidget(g)

        # --- model ----------------------------------------------------------
        g = QGroupBox("4. Model"); f = QFormLayout(g)
        self.cmb_geom = QComboBox(); self.cmb_geom.addItems(["plane", "cylinder", "sphere"])
        f.addRow("Geometry", self.cmb_geom)
        self.cmb_bcl = QComboBox(); self.cmb_bcl.addItems(["fixed (Dirichlet)", "zero flux"])
        self.cmb_bcr = QComboBox(); self.cmb_bcr.addItems(["fixed (Dirichlet)", "zero flux"])
        f.addRow("Left boundary", self.cmb_bcl); f.addRow("Right boundary", self.cmb_bcr)
        self.cmb_ic = QComboBox(); self.cmb_ic.addItems(["step", "plateau + rim"])
        f.addRow("Initial condition", self.cmb_ic)
        self.sp_x0 = _spin(-1e5, 1e5, 0.0, 3, 1, " um")
        self.sp_cl = _spin(-1e6, 1e6, 0.3, 5, 0.01)
        self.sp_cr = _spin(-1e6, 1e6, 0.18, 5, 0.01)
        self.sp_smooth = _spin(0, 50, 0.0, 2, 0.1, " um")
        f.addRow("Interface x0", self.sp_x0)
        f.addRow("Plateau left", self.sp_cl); f.addRow("Plateau right", self.sp_cr)
        f.addRow("Initial smoothing", self.sp_smooth)
        self.btn_guess = QPushButton("Guess from data")
        self.btn_guess.clicked.connect(self._guess_initial)
        f.addRow(self.btn_guess)
        self.sp_nodes = QSpinBox(); self.sp_nodes.setRange(51, 4001); self.sp_nodes.setValue(301)
        f.addRow("Grid nodes", self.sp_nodes)
        self.chk_force_num = QCheckBox("force the numerical solver")
        f.addRow(self.chk_force_num)
        self.chk_cooling = QCheckBox("linear cooling instead of isothermal")
        self.sp_Tend = _spin(300, 2000, 900, 1, 5, " C"); self.sp_Tend.setEnabled(False)
        self.chk_cooling.toggled.connect(self.sp_Tend.setEnabled)
        f.addRow(self.chk_cooling); f.addRow("  final temperature", self.sp_Tend)
        self.chk_free_x0 = QCheckBox("also fit x0")
        self.chk_free_plateaus = QCheckBox("also fit the plateaus")
        f.addRow(self.chk_free_x0); f.addRow(self.chk_free_plateaus)
        v.addWidget(g)

        # --- coefficients -----------------------------------------------------
        g = QGroupBox("5. Diffusion coefficients"); vb = QVBoxLayout(g)
        self.lst_coef = QListWidget()
        self.lst_coef.setSelectionMode(QListWidget.SingleSelection)
        self.lst_coef.currentItemChanged.connect(lambda *_: self.show_coefficient_info(brief=True))
        vb.addWidget(self.lst_coef)
        b = QPushButton("Show full details and citation")
        b.clicked.connect(self.show_coefficient_info)
        vb.addWidget(b)
        vb.addWidget(QLabel("Tick several and use Compare to run them all."))
        v.addWidget(g)

        # --- Monte Carlo --------------------------------------------------------
        g = QGroupBox("6. Monte Carlo"); f = QFormLayout(g)
        self.sp_draws = QSpinBox(); self.sp_draws.setRange(20, 100000); self.sp_draws.setValue(500)
        f.addRow("Draws", self.sp_draws)
        self.sp_seed = QSpinBox(); self.sp_seed.setRange(0, 2 ** 31 - 1); self.sp_seed.setValue(12345)
        f.addRow("Seed", self.sp_seed)
        self.chk_mc_T = QCheckBox("temperature"); self.chk_mc_T.setChecked(True)
        self.chk_mc_f = QCheckBox("oxygen fugacity"); self.chk_mc_f.setChecked(True)
        self.chk_mc_P = QCheckBox("pressure"); self.chk_mc_P.setChecked(True)
        self.chk_mc_D = QCheckBox("diffusion coefficient"); self.chk_mc_D.setChecked(True)
        self.chk_mc_n = QCheckBox("measurement noise"); self.chk_mc_n.setChecked(True)
        self.chk_mc_x = QCheckBox("distance scale")
        self.chk_contrib = QCheckBox("variance contributions (slower)")
        for c in (self.chk_mc_T, self.chk_mc_f, self.chk_mc_P, self.chk_mc_D,
                  self.chk_mc_n, self.chk_mc_x, self.chk_contrib):
            f.addRow(c)
        self.cmb_dmode = QComboBox()
        self.cmb_dmode.addItems(["auto", "covariance", "logD_at_T", "independent"])
        self.cmb_dmode.setToolTip(
            "How the diffusion-law parameters are sampled. 'independent' reproduces the "
            "assumption made by NIDIS and overestimates the uncertainty.")
        f.addRow("D sampling", self.cmb_dmode)
        v.addWidget(g)

        # --- actions ---------------------------------------------------------
        row = QWidget(); rl = QHBoxLayout(row); rl.setContentsMargins(0, 0, 0, 0)
        self.btn_fit = QPushButton("Fit"); self.btn_fit.clicked.connect(self.run_fit)
        self.btn_cmp = QPushButton("Compare"); self.btn_cmp.clicked.connect(self.run_compare)
        self.btn_mc = QPushButton("Monte Carlo"); self.btn_mc.clicked.connect(self.run_mc)
        self.btn_stop = QPushButton("Stop"); self.btn_stop.clicked.connect(self.stop_work)
        self.btn_stop.setEnabled(False)
        for b in (self.btn_fit, self.btn_cmp, self.btn_mc, self.btn_stop):
            rl.addWidget(b)
        v.addWidget(row)
        self.btn_export = QPushButton("Export results...")
        self.btn_export.clicked.connect(self.export_results)
        v.addWidget(self.btn_export)
        v.addStretch(1)

        scroll = QScrollArea(); scroll.setWidget(panel); scroll.setWidgetResizable(True)
        scroll.setMinimumWidth(410)
        return scroll

    # ------------------------------------------------------------------
    def _log(self, msg: str):
        self.log.append(str(msg))
        self.statusBar().showMessage(str(msg).splitlines()[0][:140])

    def _on_mineral_changed(self):
        key = self.cmb_mineral.currentData()
        mineral = get_mineral(key)
        self.cmb_species.blockSignals(True)
        self.cmb_species.clear()
        self.cmb_species.addItems(mineral.species_keys())
        self.cmb_species.blockSignals(False)
        self.cmb_axis.setEnabled(not mineral.isotropic)
        self._refresh_coefficients()
        self._update_fo2_label()

    def _on_axis_changed(self):
        ang = self.cmb_axis.currentText().startswith("angles")
        for s in (self.sp_alpha, self.sp_beta, self.sp_gamma):
            s.setEnabled(ang)

    def _refresh_coefficients(self):
        mineral = self.cmb_mineral.currentData()
        species = self.cmb_species.currentText()
        self.lst_coef.clear()
        for c in list_coefficients(mineral, species):
            item = QListWidgetItem(("* " if c.recommended else "  ") + c.label +
                                   ("" if c.verified else "   [UNVERIFIED]"))
            item.setData(Qt.UserRole, c.key)
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
            item.setCheckState(Qt.Checked if c.recommended else Qt.Unchecked)
            if not c.verified:
                item.setForeground(Qt.darkRed)
            self.lst_coef.addItem(item)
        if self.lst_coef.count():
            self.lst_coef.setCurrentRow(0)

    def _update_fo2_label(self):
        try:
            T = self.sp_T.value() + 273.15
            P = self.sp_P.value() * 1e6
            if self.cmb_fo2_mode.currentIndex() == 0:
                lf = log_fo2_from_delta(self.cmb_buffer.currentText(), self.sp_dbuf.value(), T, P)
                self.lbl_fo2.setText(
                    f"log10 fO2 = {lf:.3f} bar = {lf + 5:.3f} Pa  "
                    f"(buffer from Frost 1991 Table 1, re-evaluated at each sampled T)")
            else:
                lf = self.sp_dbuf.value()
                self.lbl_fo2.setText(f"log10 fO2 = {lf:.3f} bar = {lf + 5:.3f} Pa (absolute)")
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
        self.sp_xcomp.setValue(float(np.mean(self.profile.C)))
        self._log("initial condition: " + ic.describe())

    # ------------------------------------------------------------------
    def load_profile(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Load profile", "", "Tables (*.csv *.txt *.tsv *.xlsx *.xls);;All files (*)")
        if not path:
            return
        try:
            df = read_table(path)
            dlg = ColumnDialog(df, suggest_spec(df), self)
            if dlg.exec() != QDialog.Accepted:
                return
            self.profile = build_profile(df, dlg.spec(), source=path)
            p = self.profile
            self.lbl_data.setText(f"{Path(path).name}: {len(p)} points, "
                                  f"{p.x.min():.2f} to {p.x.max():.2f} um\n{p.spec.describe()}")
            self._log(f"loaded {path}\n  {p.spec.describe()}")
            for n in p.notes:
                self._log(f"  note: {n}")
            self._guess_initial()
            self.plot.show_data(p.x, p.C, p.sigma)
            self.fit_result = self.mc_result = self.compare_results = None
        except Exception:
            QMessageBox.critical(self, "Could not load the file", traceback.format_exc())

    # ------------------------------------------------------------------
    def _checked_keys(self) -> List[str]:
        keys = []
        for i in range(self.lst_coef.count()):
            it = self.lst_coef.item(i)
            if it.checkState() == Qt.Checked:
                keys.append(it.data(Qt.UserRole))
        if not keys and self.lst_coef.currentItem():
            keys = [self.lst_coef.currentItem().data(Qt.UserRole)]
        return keys

    def _conditions(self, coef) -> Conditions:
        T = self.sp_T.value() + 273.15
        P = self.sp_P.value() * 1e6
        if self.cmb_fo2_mode.currentIndex() == 0:
            lf = log_fo2_from_delta(self.cmb_buffer.currentText(), self.sp_dbuf.value(), T, P)
        else:
            lf = self.sp_dbuf.value()
        mineral = get_mineral(self.cmb_mineral.currentData())
        X = {mineral.composition_variable.key: self.sp_xcomp.value()}
        for k in coef.requires:
            X.setdefault(k, self.sp_xcomp.value())
        axis = None
        angles = None
        t = self.cmb_axis.currentText()
        if t in ("a", "b", "c"):
            axis = t
        elif t.startswith("angles"):
            angles = (self.sp_alpha.value(), self.sp_beta.value(), self.sp_gamma.value())
        return Conditions(T_K=T, P_Pa=P, log_fo2_bar=lf, X=X, axis=axis, angles_deg=angles)

    def _model(self, coef_key: str) -> DiffusionModel:
        coef = get_coefficient(coef_key)
        cond = self._conditions(coef)
        mineral = get_mineral(self.cmb_mineral.currentData())
        params = {"x0": self.sp_x0.value(), "C_left": self.sp_cl.value(),
                  "C_right": self.sp_cr.value()}
        if self.sp_smooth.value() > 0:
            params["smooth"] = self.sp_smooth.value()
        ic = InitialCondition("step", params)
        bcl = dirichlet(params["C_left"]) if self.cmb_bcl.currentIndex() == 0 else zero_flux()
        bcr = dirichlet(params["C_right"]) if self.cmb_bcr.currentIndex() == 0 else zero_flux()
        hist = None
        if self.chk_cooling.isChecked():
            hist = ThermalHistory.linear(cond.T_K, self.sp_Tend.value() + 273.15, 1.0)
        comp_key = mineral.composition_variable.key if self.chk_comp_dep.isChecked() else None
        if comp_key and comp_key not in coef.requires:
            comp_key = None
        theta = 0.0
        species = self.cmb_species.currentText()
        if mineral.key == "plagioclase" and species in ACTIVITY_A:
            theta = activity_theta(species, cond.T_K)
        return DiffusionModel(
            coefficient=coef, conditions=cond, initial=ic,
            geometry=Geometry(self.cmb_geom.currentText()),
            bc_left=bcl, bc_right=bcr, history=hist,
            beam_sigma_um=self.sp_beam.value(), n_nodes=self.sp_nodes.value(),
            comp_key=comp_key, composition_dependent=bool(comp_key),
            activity_theta=theta, force_numerical=self.chk_force_num.isChecked())

    def _free_parameters(self):
        free = ["t"]
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

    def _busy(self, on: bool, maximum: int = 0):
        for b in (self.btn_fit, self.btn_cmp, self.btn_mc):
            b.setEnabled(not on)
        self.btn_stop.setEnabled(on)
        self.progress.setVisible(on)
        self.progress.setMaximum(maximum)
        self.progress.setValue(0)

    def stop_work(self):
        for w in getattr(self, "_workers", []):
            if hasattr(w, "abort"):
                w.abort()
        self._log("stop requested")

    # ------------------------------------------------------------------
    def run_fit(self):
        if self.profile is None:
            QMessageBox.information(self, "No data", "Load a profile first.")
            return
        keys = self._checked_keys()
        if not keys:
            QMessageBox.information(self, "No coefficient", "Select a diffusion coefficient.")
            return
        try:
            model = self._model(keys[0])
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
        self.results.setPlainText(res.summary())
        self.plot.show_fit(res, None, y_label=self._y_label())
        self.methods.setPlainText(methods_paragraph(res, None, self.profile))
        self._log("fit complete: " + human_time(res.t_seconds))

    def run_compare(self):
        if self.profile is None:
            QMessageBox.information(self, "No data", "Load a profile first.")
            return
        keys = self._checked_keys()
        if len(keys) < 2:
            QMessageBox.information(self, "Select at least two",
                                    "Tick two or more coefficients to compare them.")
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
        lines = ["Comparison of diffusion coefficients", "=" * 60]
        ok = {}
        for k, r in results.items():
            if isinstance(r, Exception):
                lines.append(f"{k}\n   FAILED: {r}")
                continue
            ok[k] = r
            lines.append(f"{k}\n   t = {human_time(r.t_seconds):<14s} "
                         f"reduced chi2 = {r.stats.reduced_chi2:.3g}  R2 = {r.stats.r_squared:.4f}")
            if not r.model.coefficient.verified:
                lines.append("   [coefficients NOT verified against the primary publication]")
        if len(ok) > 1:
            ts = [r.t_seconds for r in ok.values()]
            lines += ["", f"spread: {human_time(min(ts))} to {human_time(max(ts))} "
                          f"(factor {max(ts)/min(ts):.1f} between the extremes)"]
        self.results.setPlainText("\n".join(lines))
        self.plot.show_comparison(ok, y_label=self._y_label())
        if ok:
            first = next(iter(ok.values()))
            self.fit_result = first
            self.methods.setPlainText(methods_paragraph(first, None, self.profile))
        self._log("comparison complete")

    def run_mc(self):
        if self.profile is None:
            QMessageBox.information(self, "No data", "Load a profile first.")
            return
        keys = self._checked_keys()
        if not keys:
            return
        try:
            model = self._model(keys[0])
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
        self._log(f"Monte Carlo started: {n} draws, seed {self.sp_seed.value()}, "
                  f"sources {self._budget().active_sources()}")

    def _mc_done(self, res):
        self._busy(False)
        self.mc_result = res
        from ..fitting import fit_time
        if self.fit_result is None:
            keys = self._checked_keys()
            self.fit_result = fit_time(self._model(keys[0]), self.profile.x, self.profile.C,
                                       self.profile.sigma, self._free_parameters())
        self.results.setPlainText(res.summary())
        self.plot.show_fit(self.fit_result, res, y_label=self._y_label())
        self.methods.setPlainText(methods_paragraph(self.fit_result, res, self.profile))
        self._log("Monte Carlo complete: median " + human_time(res.median))

    def _work_failed(self, msg: str):
        self._busy(False)
        self._log(msg)
        QMessageBox.critical(self, "Calculation failed", msg)

    def _y_label(self) -> str:
        mineral = get_mineral(self.cmb_mineral.currentData())
        return f"{mineral.composition_variable.label} ({self.cmb_species.currentText()})"

    # ------------------------------------------------------------------
    def show_coefficient_info(self, brief: bool = False):
        it = self.lst_coef.currentItem()
        if it is None:
            return
        c = get_coefficient(it.data(Qt.UserRole))
        text = c.describe() + "\n\nFull reference:\n  " + format_reference(c.citation)
        for k in c.secondary_citations:
            text += "\n  see also: " + format_reference(k)
        if brief:
            self.results.setPlainText(text)
        else:
            dlg = QDialog(self)
            dlg.setWindowTitle(c.label)
            dlg.resize(900, 600)
            lay = QVBoxLayout(dlg)
            te = QTextEdit(); te.setReadOnly(True); te.setPlainText(text)
            te.setFont(QFont(MONO.split(",")[0], 9))
            lay.addWidget(te)
            bb = QDialogButtonBox(QDialogButtonBox.Close)
            bb.rejected.connect(dlg.reject); bb.accepted.connect(dlg.accept)
            lay.addWidget(bb)
            dlg.exec()

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
            QMessageBox.information(self, "Exported",
                                    "Written:\n" + "\n".join(written.values()))
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
            "Every equation and coefficient carries its citation; the Methods tab lists "
            "the references actually used by the current run."))
