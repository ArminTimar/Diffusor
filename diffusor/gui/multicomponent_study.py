"""File > Multicomponent and isotope study.

Two tabs:

* **Garnet (Fe-Mg-Mn-Ca)** reads its own table with one column per cation,
  models the four components together with a tracer family
  (:mod:`diffusor.coefficients.garnet`) and fits one duration to all of them.
* **Isotopes** uses the profile and model set up in the main window and adds
  the delta-value columns of the same table: a dilute element (Li) with its
  isotopes diffusing separately, or Fe-Mg in olivine with the seven isotopes
  coupled (Oeser et al. 2026).
"""
from __future__ import annotations

import copy
import traceback
from pathlib import Path

import numpy as np
from PySide6.QtCore import QObject, Qt, Signal, Slot
from PySide6.QtWidgets import (QCheckBox, QComboBox, QDialog, QDoubleSpinBox, QFileDialog, QFormLayout,
                               QGridLayout, QGroupBox, QHBoxLayout, QLabel, QMessageBox, QPushButton,
                               QScrollArea, QSpinBox, QTabWidget, QTextEdit, QVBoxLayout, QWidget)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg

from ..coefficients.families import families_for
from ..coefficients.garnet import COMPONENTS
from ..coefficients.isotopes import ABUNDANCES, DELTA_LABEL, betas_for
from ..coefficients.registry import get as get_coefficient
from ..dataio.multicomponent import (isotope_figure, multicomponent_figure, save_isotope_results,
                                     save_multicomponent_results)
from ..dataio.profiles import read_table
from ..fitting.isotopes import (CoupledFeMgIsotopeModel, DiluteIsotopeModel, fit_coupled_isotopes,
                                fit_dilute_isotopes)
from ..fitting.multicomponent import (MulticomponentModel, fit_multicomponent_time,
                                      monte_carlo_multicomponent)
from ..solvers.geometry import Geometry
from ..solvers.history import ThermalHistory
from ..solvers.initial import InitialCondition
from ..thermo.units import human_time, length_to_m
from .widgets import fit_to_screen
from .workers import start

NONE = "(none)"


class _Worker(QObject):
    finished = Signal(object)
    failed = Signal(str)

    def __init__(self, job):
        super().__init__()
        self.job = job
        self.cancelled = False

    def abort(self):
        self.cancelled = True

    def run(self):
        try:
            self.finished.emit(self.job(lambda _f: self.cancelled))
        except InterruptedError:
            self.failed.emit("Cancelled.")
        except Exception:
            self.failed.emit(traceback.format_exc())


def _spin(lo, hi, value, decimals=3, step=None, suffix=""):
    s = QDoubleSpinBox()
    s.setRange(lo, hi)
    s.setDecimals(decimals)
    s.setValue(value)
    if step:
        s.setSingleStep(step)
    if suffix:
        s.setSuffix(suffix)
    return s


def _ends(x, y, frac=0.15):
    order = np.argsort(x)
    y = np.asarray(y, dtype=float)[order]
    n = max(2, int(frac * y.size))
    return float(np.median(y[:n])), float(np.median(y[-n:]))


class MulticomponentStudyDialog(QDialog):
    def __init__(self, owner):
        super().__init__(owner)
        self.owner = owner
        self.thread = None
        self.worker = None
        self.garnet_result = None
        self.iso_result = None
        self.iso_description = ""
        self.table = None
        self.table_path = None
        self.setWindowTitle("Multicomponent and isotope study")
        fit_to_screen(self, 1100, 860)
        layout = QVBoxLayout(self)
        self.tabs = QTabWidget()
        layout.addWidget(self.tabs, 1)
        self.tabs.addTab(self._garnet_tab(), "Garnet (Fe-Mg-Mn-Ca)")
        self.tabs.addTab(self._isotope_tab(), "Isotopes")
        bar = QHBoxLayout()
        self.status = QLabel("")
        self.status.setWordWrap(True)
        self.stop = QPushButton("Stop")
        self.stop.setEnabled(False)
        self.stop.clicked.connect(self.abort)
        bar.addWidget(self.status, 1)
        bar.addWidget(self.stop)
        layout.addLayout(bar)
        # long entries (tracer families, beta descriptions) must not widen the settings column
        for cb in self.findChildren(QComboBox):
            cb.setSizeAdjustPolicy(QComboBox.AdjustToMinimumContentsLengthWithIcon)
            cb.setMinimumContentsLength(8)
            cb.updateGeometry()

    # ================================================================== garnet
    def _garnet_tab(self):
        w = QWidget()
        outer = QHBoxLayout(w)
        left = self._settings_column(outer)
        info = QLabel("Load a table with the distance and one column per cation (Fe, Mg, Mn, Ca as "
                      "apfu, mole fractions or end-member fractions). Each row is normalised to "
                      "Fe + Mg + Mn + Ca = 1. The four components are modelled together with the "
                      "ideal ionic diffusion matrix of the chosen tracer family.")
        info.setWordWrap(True)
        left.addWidget(info)
        load = QPushButton("Load table…")
        load.clicked.connect(lambda: self.load_garnet_table())
        left.addWidget(load)
        self.g_file = QLabel("No table loaded.")
        self.g_file.setWordWrap(True)
        left.addWidget(self.g_file)

        cols = QGroupBox("Columns")
        cf = QFormLayout(cols)
        self.g_dist = QComboBox()
        self.g_unit = QComboBox()
        self.g_unit.addItems(["um", "mm", "nm", "m"])
        cf.addRow("Distance", self.g_dist)
        cf.addRow("Distance unit", self.g_unit)
        self.g_cols, self.g_sig = {}, {}
        for c in COMPONENTS:
            cb, sb = QComboBox(), QComboBox()
            self.g_cols[c], self.g_sig[c] = cb, sb
            row = QHBoxLayout()
            row.addWidget(cb, 2)
            row.addWidget(QLabel("σ"))
            row.addWidget(sb, 1)
            cf.addRow(c, row)
        self.g_sig_const = _spin(0.0, 1.0, 0.005, 4, 0.001)
        cf.addRow("σ where no column (mole fraction)", self.g_sig_const)
        left.addWidget(cols)

        cond = QGroupBox("Model")
        mf = QFormLayout(cond)
        self.g_family = QComboBox()
        self.g_family.setSizeAdjustPolicy(QComboBox.AdjustToMinimumContentsLengthWithIcon)
        self.g_family.setMinimumContentsLength(8)
        for f in families_for("garnet"):
            self.g_family.addItem(f.label, f.key)
        self.g_family.currentIndexChanged.connect(self._garnet_family_changed)
        mf.addRow("Tracer laws", self.g_family)
        self.g_T = _spin(300, 1600, 750, 1, 5, " °C")
        self.g_Tend = _spin(300, 1600, 650, 1, 5, " °C")
        self.g_cool = QCheckBox("cool linearly to")
        mf.addRow("Temperature", self.g_T)
        cr = QHBoxLayout()
        cr.addWidget(self.g_cool)
        cr.addWidget(self.g_Tend)
        mf.addRow("", cr)
        self.g_P = _spin(0.0, 10.0, 0.8, 3, 0.1, " GPa")
        mf.addRow("Pressure", self.g_P)
        self.g_fo2 = _spin(-10, 10, 0.0, 2, 0.5)
        mf.addRow("log fO2 − log fO2(graphite)", self.g_fo2)
        self.g_dep = QComboBox()
        self.g_dep.addItems(list(COMPONENTS))
        self.g_dep.setCurrentText("Ca")
        mf.addRow("Dependent component", self.g_dep)
        self.g_x0 = _spin(-1e6, 1e6, 0.0, 2, 1.0, " µm")
        self.g_fitx0 = QCheckBox("fit the step position")
        mf.addRow("Initial step at", self.g_x0)
        mf.addRow("", self.g_fitx0)
        self.g_left = QComboBox()
        self.g_right = QComboBox()
        for cmb in (self.g_left, self.g_right):
            cmb.addItem("closed (zero flux)", False)
            cmb.addItem("fixed at the plateau", True)
        mf.addRow("Left end", self.g_left)
        mf.addRow("Right end", self.g_right)
        self.g_geom = QComboBox()
        self.g_geom.addItems(["plane", "sphere", "cylinder"])
        mf.addRow("Geometry", self.g_geom)
        self.g_beam = _spin(0, 50, 0.0, 2, 0.5, " µm")
        mf.addRow("Beam σ", self.g_beam)
        self.g_fit = {c: QCheckBox(c) for c in COMPONENTS}
        fr = QHBoxLayout()
        for cb in self.g_fit.values():
            cb.setChecked(True)
            fr.addWidget(cb)
        mf.addRow("Fit to", fr)
        self.g_mc = QSpinBox()
        self.g_mc.setRange(0, 1000)
        self.g_mc.setValue(0)
        self.g_Tsig = _spin(0, 200, 15, 1, 5, " °C")
        mf.addRow("Monte Carlo draws", self.g_mc)
        mf.addRow("σ of temperature", self.g_Tsig)
        left.addWidget(cond)
        btns = QHBoxLayout()
        self.g_run = QPushButton("Fit duration")
        self.g_run.clicked.connect(self.run_garnet)
        self.g_export = QPushButton("Export…")
        self.g_export.setEnabled(False)
        self.g_export.clicked.connect(self.export_garnet)
        btns.addWidget(self.g_run)
        btns.addWidget(self.g_export)
        left.addLayout(btns)
        left.addStretch(1)

        right = QVBoxLayout()
        outer.addLayout(right, 1)
        self.g_desc = QTextEdit()
        self.g_desc.setReadOnly(True)
        self.g_desc.setMaximumHeight(170)
        right.addWidget(self.g_desc)
        self.g_scroll = QScrollArea()
        self.g_scroll.setWidgetResizable(True)
        right.addWidget(self.g_scroll, 1)
        self._garnet_family_changed()
        return w

    def _garnet_family_changed(self):
        from ..coefficients.families import get_family
        fam = get_family(self.g_family.currentData())
        self.g_fo2.setEnabled("dlogfo2_graphite" in fam.options)
        from ..coefficients.latex import equation_html
        from html import escape
        self.g_desc.setHtml(equation_html(fam.latex(), max_width=560) + "<pre>" + escape(fam.describe()) + "</pre>")

    def load_garnet_table(self, path=None):
        if path is None:
            start_dir = self.owner.output_dir() if hasattr(self.owner, "output_dir") else ""
            path, _ = QFileDialog.getOpenFileName(self, "Garnet profile", start_dir,
                                                  "Tables (*.csv *.tsv *.txt *.xlsx *.xls)")
        if not path:
            return
        try:
            self.table = read_table(path)
        except Exception as exc:
            QMessageBox.warning(self, "Cannot read the table", str(exc))
            return
        self.table_path = Path(path)
        names = [str(c) for c in self.table.columns]
        self.g_dist.clear()
        self.g_dist.addItems(names)
        guess = {"Fe": ("fe", "alm"), "Mg": ("mg", "prp", "py"), "Mn": ("mn", "sps"), "Ca": ("ca", "grs")}
        for c in COMPONENTS:
            cb, sb = self.g_cols[c], self.g_sig[c]
            cb.clear()
            sb.clear()
            cb.addItems(names)
            sb.addItem(NONE)
            sb.addItems(names)
            low = [n.lower() for n in names]
            hit = next((n for n, l in zip(names, low) if any(l.startswith(g) for g in guess[c])
                        and not any(e in l for e in ("err", "sig", "sd", "unc"))), None)
            if hit:
                cb.setCurrentText(hit)
            err = next((n for n, l in zip(names, low) if any(l.startswith(g) for g in guess[c])
                        and any(e in l for e in ("err", "sig", "sd", "unc"))), None)
            if err:
                sb.setCurrentText(err)
        dist = next((n for n in names if n.lower().startswith(("dist", "x", "position"))), names[0])
        self.g_dist.setCurrentText(dist)
        self.g_file.setText(f"{self.table_path.name}: {len(self.table)} rows.")
        try:
            x, data, _ = self._garnet_data()
            contrast = max(COMPONENTS, key=lambda c: abs(np.subtract(*_ends(x, data[c]))))
            lo, hi = _ends(x, data[contrast])
            order = np.argsort(x)
            xs, ys = x[order], data[contrast][order]
            mid = 0.5 * (lo + hi)
            i = int(np.argmin(np.abs(ys - mid)))
            self.g_x0.setValue(float(xs[i]))
        except Exception:
            pass

    def _garnet_data(self):
        df = self.table
        x = length_to_m(np.asarray(df[self.g_dist.currentText()], dtype=float), self.g_unit.currentText()) * 1e6
        raw = {c: np.asarray(df[self.g_cols[c].currentText()], dtype=float) for c in COMPONENTS}
        total = sum(raw.values())
        ok = np.isfinite(x) & np.isfinite(total) & (total > 0)
        data = {c: raw[c][ok] / total[ok] for c in COMPONENTS}
        sig = {}
        for c in COMPONENTS:
            col = self.g_sig[c].currentText()
            if col and col != NONE:
                sig[c] = np.asarray(df[col], dtype=float)[ok] / total[ok]
            elif self.g_sig_const.value() > 0:
                sig[c] = np.full(int(ok.sum()), self.g_sig_const.value())
            else:
                sig[c] = None
        x = x[ok]
        order = np.argsort(x)
        return x[order], {c: v[order] for c, v in data.items()}, {c: None if v is None else v[order] for c, v in sig.items()}

    def _garnet_model(self, x, data):
        from ..coefficients.families import get_family
        fam = get_family(self.g_family.currentData())
        x0 = self.g_x0.value()
        ends = {c: _ends(x, data[c]) for c in COMPONENTS}
        initial = {c: InitialCondition("step", {"x0": x0, "C_left": ends[c][0], "C_right": ends[c][1]})
                   for c in COMPONENTS}
        T = self.g_T.value() + 273.15
        hist = ThermalHistory.linear(T, self.g_Tend.value() + 273.15, 1.0) if self.g_cool.isChecked() else None
        opts = {"dlogfo2_graphite": self.g_fo2.value()} if "dlogfo2_graphite" in fam.options else {}
        return MulticomponentModel(
            fam, T, initial, P_Pa=self.g_P.value() * 1e9, options=opts,
            geometry=Geometry(self.g_geom.currentText()),
            left_fixed={c: ends[c][0] for c in COMPONENTS} if self.g_left.currentData() else None,
            right_fixed={c: ends[c][1] for c in COMPONENTS} if self.g_right.currentData() else None,
            history=hist, beam_sigma_um=self.g_beam.value(), dependent=self.g_dep.currentText(),
            x_grid=np.linspace(x.min(), x.max(), 201))

    def run_garnet(self):
        if self.table is None:
            self.status.setText("Load a garnet table first.")
            return
        try:
            x, data, sig = self._garnet_data()
            fitted = [c for c, cb in self.g_fit.items() if cb.isChecked()]
            model = self._garnet_model(x, data)
        except Exception as exc:
            QMessageBox.warning(self, "Cannot set up the model", str(exc))
            return
        n_mc, T_sig, fit_x0 = self.g_mc.value(), self.g_Tsig.value(), self.g_fitx0.isChecked()

        def job(cancel):
            res = fit_multicomponent_time(model, x, data, {c: sig[c] for c in fitted}, fitted,
                                          fit_x0=fit_x0, progress=cancel)
            if n_mc:
                monte_carlo_multicomponent(res, n_mc, T_sigma_K=T_sig, progress=cancel)
            return res
        self.status.setText("Fitting the duration to all components (each forward model solves the coupled system)…")
        self._launch(job, self._garnet_done)

    @Slot(object)
    def _garnet_done(self, result):
        self.garnet_result = result
        canvas = FigureCanvasQTAgg(multicomponent_figure(result))
        canvas.setMinimumHeight(260 * len(result.fitted))
        self.g_scroll.setWidget(canvas)
        self.g_desc.setPlainText(result.summary())
        self.g_export.setEnabled(True)
        self.status.setText(f"Duration {human_time(result.t_seconds)} · χ²/dof = {result.stats.reduced_chi2:.3g}.")

    def export_garnet(self):
        self._export(self.garnet_result, save_multicomponent_results)

    # ================================================================ isotopes
    def _isotope_tab(self):
        w = QWidget()
        outer = QHBoxLayout(w)
        left = self._settings_column(outer)
        info = QLabel("Set up the concentration profile, coefficient, conditions and initial state in "
                      "the main window, then read them here and choose the delta-value columns of the "
                      "same table. Li and other dilute elements: each isotope diffuses with "
                      "D_m = D (m_ref/m)^β. Olivine Fe-Mg: the three Mg and four Fe isotopes are "
                      "coupled in one exchange (Oeser et al. 2026).")
        info.setWordWrap(True)
        left.addWidget(info)
        take = QPushButton("Use the main-window profile and model")
        take.clicked.connect(self.take_main)
        left.addWidget(take)
        self.i_what = QLabel("Nothing read yet.")
        self.i_what.setWordWrap(True)
        left.addWidget(self.i_what)

        box = QGroupBox("Isotope data and β")
        self.i_form = QGridLayout(box)
        self.i_rows = {}
        for r, (key, label) in enumerate((("a", "δ column"), ("b", "second δ column"))):
            cb, sb = QComboBox(), QComboBox()
            const = _spin(0.0, 100.0, 0.1, 3, 0.05, " ‰")
            const.setMaximumWidth(110)
            beta = QComboBox()
            custom = _spin(0.0, 1.0, 0.25, 3, 0.01)
            custom.setMaximumWidth(110)
            lab = QLabel(label)
            self.i_rows[key] = (lab, cb, sb, const, beta, custom)
            self.i_form.addWidget(lab, 3 * r, 0)
            self.i_form.addWidget(cb, 3 * r, 1)
            self.i_form.addWidget(QLabel("σ"), 3 * r, 2)
            self.i_form.addWidget(sb, 3 * r, 3)
            self.i_form.addWidget(const, 3 * r, 4)
            self.i_form.addWidget(QLabel("β"), 3 * r + 1, 0)
            self.i_form.addWidget(beta, 3 * r + 1, 1, 1, 3)
            self.i_form.addWidget(custom, 3 * r + 1, 4)
            beta.currentIndexChanged.connect(lambda _i, b=beta, c=custom: c.setEnabled(b.currentData() is None))
        left.addWidget(box)

        opts = QGroupBox("Model")
        of = QFormLayout(opts)
        self.i_csig = _spin(0.0, 1e6, 0.0, 4, 0.01)
        of.addRow("Concentration σ (if none)", self.i_csig)
        self.i_source = QComboBox()
        self.i_source.setSizeAdjustPolicy(QComboBox.AdjustToMinimumContentsLengthWithIcon)
        self.i_source.setMinimumContentsLength(8)
        self.i_source.addItem("Oeser et al. (2026) Fe and Mg tracer laws", "oeser")
        self.i_source.addItem("main-window Fe-Mg law with D*Fe/D*Mg", "ratio")
        of.addRow("Olivine tracer coefficients", self.i_source)
        self.i_ratio = _spin(0.05, 50.0, 2.0, 2, 0.1)
        of.addRow("D*Fe / D*Mg", self.i_ratio)
        self.i_fitbeta = QCheckBox("also fit β (dilute element)")
        of.addRow("", self.i_fitbeta)
        left.addWidget(opts)
        btns = QHBoxLayout()
        self.i_run = QPushButton("Fit duration")
        self.i_run.clicked.connect(self.run_isotopes)
        self.i_export = QPushButton("Export…")
        self.i_export.setEnabled(False)
        self.i_export.clicked.connect(self.export_isotopes)
        btns.addWidget(self.i_run)
        btns.addWidget(self.i_export)
        left.addLayout(btns)
        left.addStretch(1)
        right = QVBoxLayout()
        outer.addLayout(right, 1)
        self.i_desc = QTextEdit()
        self.i_desc.setReadOnly(True)
        self.i_desc.setMaximumHeight(170)
        right.addWidget(self.i_desc)
        self.i_scroll = QScrollArea()
        self.i_scroll.setWidgetResizable(True)
        right.addWidget(self.i_scroll, 1)
        self.i_model = None
        self.i_mode = None
        return w

    def take_main(self):
        o = self.owner
        if not o._ready():
            return
        try:
            self.i_model = o._model(o._checked_keys()[0])
        except Exception as exc:
            QMessageBox.warning(self, "Cannot read the main-window model", str(exc))
            return
        self.i_profile = o._fit_profile()
        c = self.i_model.coefficient
        if c.mineral == "olivine" and c.species == "Fe-Mg":
            self.i_mode, elements = "coupled", ("Fe", "Mg")
        elif c.species in ABUNDANCES:
            self.i_mode, elements = "dilute", (c.species,)
        else:
            self.i_mode = None
            self.i_what.setText(f"{c.species} in {c.mineral}: no isotope model. Isotope data are included "
                                f"for {', '.join(sorted(ABUNDANCES))}, and for Fe-Mg in olivine.")
            return
        names = [NONE] + ([str(n) for n in self.i_profile.raw.columns] if self.i_profile.raw is not None else [])
        axis = self.i_model.conditions.axis
        for key, el in zip(("a", "b"), elements + (None,) * (2 - len(elements))):
            lab, cb, sb, const, beta, custom = self.i_rows[key]
            for wdg in (lab, cb, sb, const, beta, custom):
                wdg.setVisible(el is not None)
            if el is None:
                continue
            lab.setText(f"{DELTA_LABEL[el]} column")
            cb.clear(); sb.clear(); beta.clear()
            cb.addItems(names); sb.addItems(names)
            hit = next((n for n in names if el.lower() in n.lower() and ("d" in n.lower() or "δ" in n)
                        and not any(e in n.lower() for e in ("err", "sig", "sd"))), NONE)
            cb.setCurrentText(hit)
            for b in betas_for(el, c.mineral):
                beta.addItem(b.describe()[:110], b.key)
                if axis and b.direction == axis:
                    beta.setCurrentIndex(beta.count() - 1)
            beta.addItem("custom β", None)
            custom.setEnabled(beta.currentData() is None)
        self.i_source.setEnabled(self.i_mode == "coupled")
        self.i_ratio.setEnabled(self.i_mode == "coupled")
        self.i_fitbeta.setEnabled(self.i_mode == "dilute")
        if self.i_mode == "coupled":
            try:
                fe, mg = get_coefficient("ol_Fe_oeser2026_tracer"), get_coefficient("ol_Mg_oeser2026_tracer")
                cond = self.i_model.conditions.replace(X={"XFe": 0.1})
                if cond.axis is not None or cond.angles_deg is not None:
                    self.i_ratio.setValue(float(fe.D(cond) / mg.D(cond)))
            except Exception:
                pass
        self.i_what.setText(f"{c.label}: {len(self.i_profile)} points, "
                            f"{'coupled Fe-Mg isotope model' if self.i_mode == 'coupled' else 'dilute-element isotope model'}.")
        self.i_desc.setPlainText("\n".join(self.i_model.warnings()))

    def _beta(self, key):
        from ..coefficients.isotopes import BETAS
        beta, custom = self.i_rows[key][4], self.i_rows[key][5]
        k = beta.currentData()
        if k is None:
            return custom.value(), "custom"
        b = BETAS[k]
        return b.beta, k + (", a published maximum" if b.upper_bound else "") + (
            ", assumed in the source" if b.assumed else "")

    def _delta(self, key):
        _, cb, sb, const, _, _ = self.i_rows[key]
        name = cb.currentText()
        if not name or name == NONE:
            return None, None
        d = self.i_profile.column(name)
        s_name = sb.currentText()
        s = self.i_profile.column(s_name) if s_name and s_name != NONE else np.full(d.shape, const.value())
        return d, s

    def run_isotopes(self):
        if self.i_mode is None or self.i_model is None:
            self.status.setText("Read the main-window profile and model first.")
            return
        p, dm = self.i_profile, self.i_model
        C_sig = p.sigma if p.sigma is not None else (np.full(p.x.shape, self.i_csig.value()) if self.i_csig.value() > 0 else None)
        try:
            if self.i_mode == "dilute":
                el = dm.coefficient.species
                d, ds = self._delta("a")
                beta, bkey = self._beta("a")
                dl, dr = _ends(p.x, d) if d is not None else (0.0, 0.0)
                iso = DiluteIsotopeModel(dm, el, beta, dl, dr)
                fit_beta = self.i_fitbeta.isChecked()
                self.iso_description = (f"dilute {el} isotopes, beta = {beta:g} ({bkey}), D of {iso.reference_mass}{el} "
                                        f"from {dm.coefficient.key}; delta plateaus {dl:.3g} and {dr:.3g} per mil")

                def job(cancel):
                    return fit_dilute_isotopes(iso, p.x, p.C, C_sig, d, ds, fit_beta=fit_beta, progress=cancel)
            else:
                (dfe, sfe), (dmg, smg) = self._delta("a"), self._delta("b")
                (bfe, kfe), (bmg, kmg) = self._beta("a"), self._beta("b")
                to_xfe = lambda C: dm.comp_offset + dm.comp_scale * np.asarray(C, dtype=float)
                ic = copy.deepcopy(dm.initial)
                if ic.kind != "step":
                    raise ValueError("the coupled model needs a step initial state")
                ic.params["C_left"], ic.params["C_right"] = (float(to_xfe(ic.params["C_left"])),
                                                             float(to_xfe(ic.params["C_right"])))
                dfe_lr = _ends(p.x, dfe) if dfe is not None else (0.0, 0.0)
                dmg_lr = _ends(p.x, dmg) if dmg is not None else (0.0, 0.0)
                src = self.i_source.currentData()
                kw = (dict(tracer_fe=get_coefficient("ol_Fe_oeser2026_tracer"), tracer_mg=get_coefficient("ol_Mg_oeser2026_tracer"))
                      if src == "oeser" else dict(interdiffusion=dm.coefficient, ratio_fe_mg=self.i_ratio.value()))
                cm = CoupledFeMgIsotopeModel(
                    dm.conditions, ic, bfe, bmg, delta_fe=dfe_lr, delta_mg=dmg_lr, geometry=dm.geometry,
                    left_fixed=dm.bc_left.kind == "dirichlet", right_fixed=dm.bc_right.kind == "dirichlet",
                    history=dm.history, fo2_buffer=dm.fo2_buffer, beam_sigma_um=dm.beam_sigma_um, **kw)
                obs, sig = {}, {}
                if C_sig is not None:
                    obs["XFe"], sig["XFe"] = to_xfe(p.C), np.abs(dm.comp_scale) * np.asarray(C_sig, dtype=float)
                if dfe is not None:
                    obs["d56Fe"], sig["d56Fe"] = dfe, sfe
                if dmg is not None:
                    obs["d26Mg"], sig["d26Mg"] = dmg, smg
                self.iso_description = (f"coupled Fe-Mg isotopes in olivine (Oeser et al. 2026), beta_Fe = {bfe:g} ({kfe}), "
                                        f"beta_Mg = {bmg:g} ({kmg}), tracer coefficients: "
                                        + ("Oeser et al. (2026) Table 4" if src == "oeser" else
                                           f"{dm.coefficient.key} with D*Fe/D*Mg = {self.i_ratio.value():g}"))

                def job(cancel):
                    return fit_coupled_isotopes(cm, p.x, obs, sig, progress=cancel)
        except Exception as exc:
            QMessageBox.warning(self, "Cannot set up the isotope model", str(exc))
            return
        self.status.setText("Fitting the duration to the concentration and isotope profiles…")
        self._launch(job, self._isotopes_done)

    @Slot(object)
    def _isotopes_done(self, result):
        self.iso_result = result
        canvas = FigureCanvasQTAgg(isotope_figure(result))
        canvas.setMinimumHeight(270 * len(result.observed))
        self.i_scroll.setWidget(canvas)
        self.i_desc.setPlainText(self.iso_description + "\n\n" + result.summary())
        self.i_export.setEnabled(True)
        self.status.setText(f"Duration {human_time(result.t_seconds)} · χ²/dof = {result.stats.reduced_chi2:.3g}.")

    def export_isotopes(self):
        self._export(self.iso_result, lambda d, r: save_isotope_results(d, r, self.iso_description))

    # ================================================================ shared
    @staticmethod
    def _settings_column(outer):
        """A fixed-width, scrollable settings column, so the plots get the rest of the width."""
        holder = QWidget()
        column = QVBoxLayout(holder)
        column.setContentsMargins(0, 0, 6, 0)
        area = QScrollArea()
        area.setWidgetResizable(True)
        area.setWidget(holder)
        area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        area.setFixedWidth(500)
        outer.addWidget(area, 0)
        return column

    def _launch(self, job, done):
        if self.thread is not None and self.thread.isRunning():
            self.status.setText("Still working. Wait or press Stop.")
            return
        for b in (self.g_run, self.i_run):
            b.setEnabled(False)
        self.stop.setEnabled(True)
        self.worker = _Worker(job)
        self.worker.finished.connect(done)
        self.worker.failed.connect(self._failed)
        self.thread = start(self.worker)
        self.thread.finished.connect(self._idle)

    @Slot()
    def _idle(self):
        for b in (self.g_run, self.i_run):
            b.setEnabled(True)
        self.stop.setEnabled(False)

    @Slot(str)
    def _failed(self, message):
        self.status.setText(message.strip().splitlines()[-1])

    def abort(self):
        if self.worker is not None:
            self.worker.abort()

    def closeEvent(self, event):
        if self.thread is not None and self.thread.isRunning():
            self.abort()
            self.status.setText("Cancelling. Close this window once the current calculation finishes.")
            event.ignore()
        else:
            super().closeEvent(event)

    def _export(self, result, save):
        if result is None:
            return
        start_dir = self.owner.output_dir() if hasattr(self.owner, "output_dir") else ""
        directory = QFileDialog.getExistingDirectory(self, "Export study", start_dir)
        if directory:
            try:
                save(directory, result)
                self.status.setText(f"Exported to {directory}")
            except Exception as exc:
                QMessageBox.critical(self, "Export failed", str(exc))
