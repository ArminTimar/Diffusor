"""A nonmodal spreadsheet of crystals that share settings, each keeping its own values.

The first row holds the settings of all crystals; every row under it is one crystal.
A value typed in the first row applies to every crystal. A value typed in a crystal's
row belongs to that crystal only and is shown in bold. The model behind the table is
:class:`crystal_setup.Batch`; fitting goes through :func:`crystal_setup.run_crystal`,
so each crystal is fitted exactly as the main window would fit it.
"""
import traceback
from pathlib import Path

import numpy as np
from PySide6.QtCore import QObject, Qt, QTimer, Signal, Slot
from PySide6.QtGui import QAction, QBrush, QColor, QFont
from PySide6.QtWidgets import (QAbstractItemView, QComboBox, QDialog, QFileDialog, QHBoxLayout,
                               QLabel, QLineEdit, QMenu, QMessageBox, QPushButton,
                               QStyledItemDelegate, QTableWidget, QTableWidgetItem, QVBoxLayout)

from .. import crystal_setup as cs
from ..coefficients import get as get_coefficient, list_coefficients
from ..dataio import read_table, suggest_spec
from ..minerals import MINERALS, get_mineral
from ..thermo import available_buffers
from ..thermo.units import human_time
from . import theme
from .widgets import fit_to_screen
from .workers import start

SHARED_ROW = 0
AXIS_NAMES = {"experiments": "as in the experiments", "a": "a-axis", "b": "b-axis",
              "c": "c-axis", "angles": "angles to a, b, c"}
BOUNDARY_NAMES = {"far": "plateau continues", "rim_melt": "rim, held by melt",
                  "rim_closed": "rim, closed", "centre": "crystal centre"}
ANGLE_KEYS = {"angle_a": 0, "angle_b": 1, "angle_c": 2}
FO2_NAMES = {"buffer": "offset from a buffer", "absolute": "log fO2 (bar)"}

# (key, header, kind). Keys starting with "_" are not settings; angle_* edit angles_deg.
COLUMNS = [
    ("_label", "Crystal", "text"), ("_group", "Group", "text"),
    ("_status", "Status", "info"), ("_time", "Best-fit time", "info"),
    ("_interval", "Monte Carlo 68 %", "info"), ("_chi2", "χ² / dof", "info"),
    ("_file", "File", "info"), ("mineral", "Mineral", "choice"), ("species", "Species", "choice"),
    ("coefficients", "Diffusion law", "choice"),
    ("T_C", "T (°C)", "float"), ("sigma_T_K", "± T (K)", "float"),
    ("P_MPa", "P (MPa)", "float"), ("sigma_P_MPa", "± P (MPa)", "float"),
    ("fo2_mode", "fO2 given as", "choice"), ("buffer", "fO2 buffer", "choice"), ("fo2_value", "fO2 value", "float"),
    ("sigma_fo2", "± fO2", "float"), ("x_composition", "Host composition", "float"),
    ("axis", "Traverse direction", "choice"),
    ("angle_a", "∠ to a (°)", "float"), ("angle_b", "∠ to b (°)", "float"),
    ("angle_c", "∠ to c (°)", "float"),
    ("geometry", "Geometry", "choice"), ("bc_left", "Left end", "choice"),
    ("bc_right", "Right end", "choice"),
    ("x0_um", "Step position (µm)", "float"), ("C_left", "Left plateau", "float"),
    ("C_right", "Right plateau", "float"), ("fit_x0", "Fit step", "bool"),
    ("fit_plateaus", "Fit plateaus", "bool"), ("beam_sigma_um", "Beam σ (µm)", "float"),
    ("mc_draws", "MC draws", "int"),
]
KEYS = [c[0] for c in COLUMNS]
HIDDEN_BY_DEFAULT = {"sigma_T_K", "sigma_P_MPa", "sigma_fo2", "geometry", "bc_left", "bc_right",
                     "beam_sigma_um", "mc_draws"}
OWN = QColor(theme.ACCENT_SOFT)
SHARED = QColor("#EFEBE4")


def _setting_key(key: str) -> str:
    return "angles_deg" if key in ANGLE_KEYS else key


def _get(settings, key):
    if key in ANGLE_KEYS:
        return settings["angles_deg"][ANGLE_KEYS[key]]
    if key == "coefficients":
        return settings["coefficients"][0] if settings["coefficients"] else ""
    return settings[key]


def _new_value(settings, key, value):
    """The value to store under ``_setting_key(key)`` when ``key`` is set to ``value``."""
    if key in ANGLE_KEYS:
        angles = list(settings["angles_deg"])
        angles[ANGLE_KEYS[key]] = float(value)
        return angles
    if key == "coefficients":
        return [value]
    return value


def choices(key: str, settings) -> list:
    """(stored value, shown text) pairs for a choice column, for these settings."""
    if key == "mineral":
        return [(k, m.name) for k, m in MINERALS.items()]
    if key == "species":
        m = get_mineral(settings["mineral"])
        return [(sp, sp) for sp in m.species_keys() if list_coefficients(m.key, sp)]
    if key == "coefficients":
        return [(c.key, c.label) for c in list_coefficients(settings["mineral"], settings["species"])]
    if key == "buffer":
        return [(b, b) for b in available_buffers()]
    if key == "fo2_mode":
        return list(FO2_NAMES.items())
    if key == "axis":
        return list(AXIS_NAMES.items())
    if key == "geometry":
        return [(g, g) for g in ("plane", "cylinder", "sphere")]
    if key in ("bc_left", "bc_right"):
        return list(BOUNDARY_NAMES.items())
    return []


def shown(key: str, value) -> str:
    if isinstance(value, bool):
        return ""
    if key == "mineral":
        return MINERALS[value].name if value in MINERALS else str(value)
    if key == "coefficients":
        try:
            return get_coefficient(value).label if value else "none chosen"
        except Exception:
            return f"{value} (not in this version)"
    if key == "axis":
        return AXIS_NAMES.get(value, str(value))
    if key == "fo2_mode":
        return FO2_NAMES.get(value, str(value))
    if key in ("bc_left", "bc_right"):
        return BOUNDARY_NAMES.get(value, str(value))
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def consistent(settings) -> dict:
    """Species and law that suit the mineral: what to change so a mineral or species change is valid."""
    fix = {}
    species = [sp for sp, _ in choices("species", settings)]
    sp = settings["species"] if settings["species"] in species else (species[0] if species else "")
    if sp != settings["species"]:
        fix["species"] = sp
    laws = list_coefficients(settings["mineral"], sp)
    keys = [c.key for c in laws]
    if not settings["coefficients"] or settings["coefficients"][0] not in keys:
        rec = [c.key for c in laws if c.recommended]
        fix["coefficients"] = [(rec or keys or [""])[0]] if keys else []
    return fix


class ChoiceDelegate(QStyledItemDelegate):
    """A drop-down list in the choice columns, filled for the row's own mineral and species."""

    def __init__(self, dialog):
        super().__init__(dialog)
        self.dialog = dialog

    def createEditor(self, parent, option, index):
        key = KEYS[index.column()]
        combo = QComboBox(parent)
        for value, text in choices(key, self.dialog.row_settings(index.row())):
            combo.addItem(text, value)
        combo.activated.connect(lambda _i, c=combo: self.commitData.emit(c))
        return combo

    def setEditorData(self, editor, index):
        key = KEYS[index.column()]
        i = editor.findData(_get(self.dialog.row_settings(index.row()), key))
        editor.setCurrentIndex(max(i, 0))

    def setModelData(self, editor, model, index):
        row, key, value = index.row(), KEYS[index.column()], editor.currentData()
        # the table is rebuilt after an edit, so it happens once the editor has closed
        QTimer.singleShot(0, lambda: self.dialog.set_value(row, key, value))


class BatchWorker(QObject):
    progress = Signal(int, int, str)
    crystal_done = Signal(int, object, object)
    crystal_failed = Signal(int, str)
    finished = Signal()

    def __init__(self, jobs, monte_carlo: bool, workers: int):
        super().__init__()
        self.jobs = jobs                  # (row index, label, settings, data)
        self.monte_carlo = monte_carlo
        self.workers = workers
        self.cancelled = False

    def abort(self):
        self.cancelled = True

    def run(self):
        n = len(self.jobs)
        for k, (i, label, settings, data) in enumerate(self.jobs):
            if self.cancelled:
                break
            self.progress.emit(k, n, label)
            try:
                fit, mc = cs.run_crystal(settings, data, self.monte_carlo, self.workers,
                                         progress=lambda *_: self.cancelled)
                if self.cancelled and mc is not None:
                    break                 # a stopped Monte Carlo is not a result
                self.crystal_done.emit(i, fit, mc)
            except Exception as exc:
                self.crystal_failed.emit(i, f"{type(exc).__name__}: {exc}\n\n{traceback.format_exc()}")
        self.finished.emit()


class BatchStudyDialog(QDialog):
    def __init__(self, owner):
        super().__init__(owner)
        self.owner = owner
        self.batch = cs.Batch()
        self.thread = None
        self.worker = None
        self._filling = False
        self._in_main = None              # the crystal last opened in the main window
        self.setWindowTitle("Batch of crystals")
        fit_to_screen(self, 1260, 760)
        layout = QVBoxLayout(self)
        info = QLabel(
            "Set up one crystal in the main window, then add profiles here: each takes the main "
            "window's settings and finds its own step and plateaus. The first row, <b>All "
            "crystals</b>, sets a value for every crystal. A value typed in a crystal's row is "
            "that crystal's own and is shown in <b>bold</b>; right-click a cell to give it back "
            "the shared value. Open a crystal in the main window to look at its fit, then take "
            "its settings back.")
        info.setWordWrap(True)
        layout.addWidget(info)

        row = QHBoxLayout()
        row.addWidget(QLabel("Group for new crystals:"))
        self.group_edit = QLineEdit()
        self.group_edit.setPlaceholderText("for example 2011 sub-Plinian")
        self.group_edit.setMaximumWidth(260)
        row.addWidget(self.group_edit)
        self.add_files = QPushButton("Add profiles…")
        self.add_current = QPushButton("Add the main window's crystal")
        self.add_saved = QPushButton("Add saved setups…")
        self.remove = QPushButton("Remove selected")
        for b in (self.add_files, self.add_current, self.add_saved, self.remove):
            row.addWidget(b)
        row.addStretch(1)
        layout.addLayout(row)

        row = QHBoxLayout()
        self.open_main = QPushButton("Open in main window")
        self.take_main = QPushButton("Take back from main window")
        self.fit = QPushButton("Fit")
        self.fit_mc = QPushButton("Fit and Monte Carlo")
        self.stop = QPushButton("Stop")
        self.stop.setEnabled(False)
        for b in (self.open_main, self.take_main, self.fit, self.fit_mc, self.stop):
            row.addWidget(b)
        row.addStretch(1)
        self.columns_btn = QPushButton("Columns")
        self.to_population = QPushButton("Send to population")
        self.save = QPushButton("Save batch…")
        self.load = QPushButton("Open batch…")
        for b in (self.columns_btn, self.to_population, self.save, self.load):
            row.addWidget(b)
        layout.addLayout(row)

        self.table = QTableWidget(0, len(COLUMNS))
        self.table.setHorizontalHeaderLabels([c[1] for c in COLUMNS])
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.DoubleClicked | QAbstractItemView.EditKeyPressed
                                   | QAbstractItemView.AnyKeyPressed)
        self.table.verticalHeader().setVisible(False)
        self.table.setContextMenuPolicy(Qt.CustomContextMenu)
        self.table.customContextMenuRequested.connect(self._context_menu)
        delegate = ChoiceDelegate(self)
        for j, (_, _, kind) in enumerate(COLUMNS):
            if kind == "choice":
                self.table.setItemDelegateForColumn(j, delegate)
        self.table.itemChanged.connect(self._edited)
        layout.addWidget(self.table, 1)
        self.status = QLabel("No crystals yet.")
        self.status.setWordWrap(True)
        layout.addWidget(self.status)

        menu = QMenu(self)
        hidden = self._hidden_columns()
        for j, (key, header, _) in enumerate(COLUMNS):
            if key in ("_label",):
                continue
            a = QAction(header, menu, checkable=True)
            a.setChecked(key not in hidden)
            a.toggled.connect(lambda on, j=j: self._show_column(j, on))
            menu.addAction(a)
            self.table.setColumnHidden(j, key in hidden)
        self.columns_btn.setMenu(menu)

        self.add_files.clicked.connect(self.add_profiles)
        self.add_current.clicked.connect(self.add_from_main)
        self.add_saved.clicked.connect(self.add_setups)
        self.remove.clicked.connect(self.remove_selected)
        self.open_main.clicked.connect(self.open_in_main)
        self.take_main.clicked.connect(self.take_from_main)
        self.fit.clicked.connect(lambda: self.run(False))
        self.fit_mc.clicked.connect(lambda: self.run(True))
        self.stop.clicked.connect(self.abort)
        self.to_population.clicked.connect(self.send_to_population)
        self.save.clicked.connect(self.save_batch)
        self.load.clicked.connect(self.open_batch)
        self.refresh()

    # --- columns ---------------------------------------------------------------------
    def _hidden_columns(self):
        stored = self.owner._settings().value("batch/hidden_columns", None)
        if stored is None:
            return set(HIDDEN_BY_DEFAULT)
        return {k for k in str(stored).split(",") if k}

    def _show_column(self, j, on):
        self.table.setColumnHidden(j, not on)
        hidden = [KEYS[i] for i in range(len(KEYS)) if self.table.isColumnHidden(i)]
        self.owner._settings().setValue("batch/hidden_columns", ",".join(hidden))

    # --- the table -------------------------------------------------------------------
    def crystal_index(self, row: int):
        """The crystal of a table row, or None for the shared row."""
        return None if row == SHARED_ROW else row - 1

    def row_settings(self, row: int):
        i = self.crystal_index(row)
        return self.batch.shared if i is None else self.batch.settings(i)

    def refresh(self):
        self._filling = True
        try:
            self.table.setRowCount(len(self.batch) + 1)
            for row in range(len(self.batch) + 1):
                for j, (key, _, kind) in enumerate(COLUMNS):
                    self.table.setItem(row, j, self._item(row, key, kind))
            self.table.resizeColumnsToContents()
            for j, (key, _, _) in enumerate(COLUMNS):
                if key == "coefficients":
                    self.table.setColumnWidth(j, min(self.table.columnWidth(j), 260))
        finally:
            self._filling = False
        busy = self.thread is not None and self.thread.isRunning()
        n = len(self.batch)
        fitted = sum(c.result is not None for c in self.batch.crystals)
        if not busy:
            self.status.setText(f"{n} crystals, {fitted} with a result." if n else
                                "No crystals yet. Add profiles, saved setups, or the main window's crystal.")
        for b in (self.fit, self.fit_mc, self.open_main, self.take_main, self.remove, self.save,
                  self.to_population):
            b.setEnabled(bool(n) and not busy)

    def _item(self, row, key, kind):
        i = self.crystal_index(row)
        item = QTableWidgetItem()
        flags = Qt.ItemIsSelectable | Qt.ItemIsEnabled
        if i is None:
            item.setBackground(QBrush(SHARED))
            if key == "_label":
                item.setText("All crystals")
                f = item.font(); f.setBold(True); item.setFont(f)
                item.setFlags(flags)
                return item
            if key.startswith("_"):
                item.setFlags(flags)
                return item
            settings = self.batch.shared
        else:
            c = self.batch.crystals[i]
            settings = self.batch.settings(i)
            info = {"_label": c.label, "_group": c.group,
                    "_file": Path(c.data.source).name if c.data.source else "",
                    "_status": c.status}
            if key in info:
                item.setText(info[key])
                if key == "_file":
                    item.setToolTip(c.data.source)
                if key == "_status" and c.result and c.result.get("warnings"):
                    item.setToolTip("\n\n".join(c.result["warnings"]))
                    item.setForeground(QBrush(QColor(theme.WARN)))
            elif key in ("_time", "_interval", "_chi2"):
                item.setText(self._result_text(c.result, key))
            if key.startswith("_"):
                item.setFlags(flags | (Qt.ItemIsEditable if kind == "text" else Qt.NoItemFlags))
                return item
            if self.batch.is_own(i, _setting_key(key)):
                item.setBackground(QBrush(OWN))
                f = item.font(); f.setBold(True); item.setFont(f)
        value = _get(settings, key)
        if kind == "bool":
            item.setFlags(flags | Qt.ItemIsUserCheckable)
            item.setCheckState(Qt.Checked if value else Qt.Unchecked)
        else:
            item.setText(shown(key, value))
            item.setFlags(flags | Qt.ItemIsEditable)
        if key in ANGLE_KEYS and settings["axis"] != "angles":
            item.setForeground(QBrush(QColor(theme.TEXT_FAINT)))
            item.setToolTip("Used only when the traverse direction is given as angles.")
        elif key == "buffer" and settings["fo2_mode"] != "buffer":
            item.setForeground(QBrush(QColor(theme.TEXT_FAINT)))
            item.setToolTip("Not used: fO2 is given as an absolute value.")
        elif key == "coefficients":
            item.setToolTip(shown(key, value))
        return item

    @staticmethod
    def _result_text(r, key):
        if not r:
            return ""
        if key == "_time":
            return human_time(r["t_seconds"])
        if key == "_chi2":
            v = r.get("reduced_chi2")
            return "" if v is None else f"{v:.3g}"
        t = (r.get("montecarlo") or {}).get("times_s")
        if not t:
            return ""
        lo, hi = np.percentile(np.asarray(t, dtype=float), [16, 84])
        return f"{human_time(lo)} – {human_time(hi)}"

    # --- editing ---------------------------------------------------------------------
    @Slot(object)
    def _edited(self, item):
        if self._filling:
            return
        row, key = item.row(), KEYS[item.column()]
        kind = COLUMNS[item.column()][2]
        i = self.crystal_index(row)
        if key in ("_label", "_group") and i is not None:
            text = item.text().strip()
            c = self.batch.crystals[i]
            if key == "_group":
                c.group = text
            elif text and text not in {x.label for x in self.batch.crystals if x is not c}:
                c.label = text
            QTimer.singleShot(0, self.refresh)
            return
        if kind == "bool":
            value = item.checkState() == Qt.Checked
        elif kind in ("float", "int"):
            try:
                value = float(item.text().replace(",", "."))
                if kind == "int":
                    value = int(round(value))
            except ValueError:
                self.status.setText(f"'{item.text()}' is not a number.")
                QTimer.singleShot(0, self.refresh)
                return
        else:
            return
        QTimer.singleShot(0, lambda: self.set_value(row, key, value))

    def set_value(self, row: int, key: str, value):
        """Set a value in the shared row (for every crystal) or in one crystal's row."""
        i = self.crystal_index(row)
        skey = _setting_key(key)
        if i is None:
            new = _new_value(self.batch.shared, key, value)
            keep_own = False
            n_own = self.batch.own_count(skey)
            if n_own:
                answer = QMessageBox.question(
                    self, "Crystals with their own value",
                    f"{n_own} crystal(s) have their own value here. Give them the new value too?",
                    QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel, QMessageBox.Yes)
                if answer == QMessageBox.Cancel:
                    self.refresh()
                    return
                keep_own = answer == QMessageBox.No
            self.batch.set_shared(skey, new, keep_own=keep_own)
            if key in ("mineral", "species"):
                for k, v in consistent(self.batch.shared).items():
                    self.batch.set_shared(k, v, keep_own=keep_own)
        else:
            self.batch.set_own(i, skey, _new_value(self.batch.settings(i), key, value))
            if key in ("mineral", "species"):
                for k, v in consistent(self.batch.settings(i)).items():
                    self.batch.set_own(i, k, v)
        self.refresh()

    def _context_menu(self, pos):
        index = self.table.indexAt(pos)
        if not index.isValid():
            return
        row, key = index.row(), KEYS[index.column()]
        i = self.crystal_index(row)
        menu = QMenu(self)
        skey = _setting_key(key)
        if not key.startswith("_"):
            if i is None:
                a = menu.addAction("Give every crystal this value")
                a.setEnabled(self.batch.own_count(skey) > 0)
                a.triggered.connect(lambda: self._clear_column(skey))
            else:
                a = menu.addAction("Use the value of All crystals")
                a.setEnabled(self.batch.is_own(i, skey))
                a.triggered.connect(lambda: (self.batch.clear_own(i, skey), self.refresh()))
                b = menu.addAction("Use this value for all crystals")
                b.triggered.connect(lambda: self._copy_to_all(i, skey))
        if i is not None:
            menu.addSeparator()
            g = menu.addAction("Guess step and plateaus from the data")
            g.triggered.connect(lambda: self.guess_selected())
            o = menu.addAction("Open in main window")
            o.triggered.connect(self.open_in_main)
        if not menu.isEmpty():
            menu.exec(self.table.viewport().mapToGlobal(pos))

    def _clear_column(self, skey):
        self.batch.set_shared(skey, self.batch.shared[skey], keep_own=False)
        self.refresh()

    def _copy_to_all(self, i, skey):
        self.batch.set_shared(skey, self.batch.settings(i)[skey], keep_own=False)
        self.refresh()

    def _selected(self):
        rows = sorted({ix.row() for ix in self.table.selectedIndexes()})
        return [self.crystal_index(r) for r in rows if r != SHARED_ROW]

    def guess_selected(self):
        for i in self._selected():
            for k, v in cs.guess_initial(self.batch.settings(i), self.batch.crystals[i].data).items():
                self.batch.set_own(i, k, v)
        self.refresh()

    # --- adding and removing ---------------------------------------------------------
    def _template(self):
        """The main window's settings become the shared settings of an empty batch."""
        if not len(self.batch):
            self.batch.shared = cs.complete(self.owner.current_settings())

    def add_profiles(self):
        from .main_window import ColumnDialog
        start = str(self.owner._settings().value("recent/last_dir", "") or "")
        paths, _ = QFileDialog.getOpenFileNames(
            self, "Add profiles", start if start and Path(start).is_dir() else "",
            "Tables (*.csv *.txt *.tsv *.xlsx *.xls);;All files (*)")
        if not paths:
            return
        self._template()
        spec, an, added, failed = None, (None, True), 0, []
        for path in paths:
            try:
                df = read_table(path)
                cols = {str(c) for c in df.columns}
                used = [] if spec is None else [spec.distance_column, spec.column_a, spec.column_b,
                                                spec.sigma_a_column, spec.sigma_b_column, an[0]]
                if spec is None or not all(c in cols for c in used if c):
                    dlg = ColumnDialog(df, spec or suggest_spec(df), self, an=an)
                    dlg.setWindowTitle(f"Columns of {Path(path).name}")
                    if dlg.exec() != QDialog.Accepted:
                        break
                    spec, an = dlg.spec(), dlg.an_column()
                data = cs.CrystalData(df, spec, str(path), an_column=an[0], an_is_percent=an[1])
                if len(data.profile()) < 3:
                    raise ValueError("fewer than three points")
                self.batch.add(Path(path).stem, data, group=self.group_edit.text().strip())
                added += 1
            except Exception as exc:
                failed.append(f"{Path(path).name}: {exc}")
        self.refresh()
        if failed:
            QMessageBox.warning(self, "Some profiles were not added", "\n".join(failed))
        if added:
            self.status.setText(f"{added} profile(s) added with the shared settings; step and "
                                "plateaus were found from each profile.")

    def add_from_main(self):
        data = self.owner.current_data()
        if data is None:
            QMessageBox.information(self, "Nothing to add", "Load a profile in the main window first.")
            return
        self._template()
        result = None if self.owner._stale else cs.result_to_json(self.owner.fit_result,
                                                                  self.owner.mc_result)
        self.batch.add(data.label, data, self.owner.current_settings(),
                       self.group_edit.text().strip(), result)
        self.refresh()

    def add_setups(self):
        paths, _ = QFileDialog.getOpenFileNames(self, "Add saved setups", self.owner.output_dir(),
                                                "Diffusor setups (*.diffusor)")
        failed = []
        for path in paths:
            try:
                crystals = cs.crystals_in(cs.read_file(path))
                if crystals and not len(self.batch):
                    self.batch.shared = cs.complete(crystals[0].settings)
                for c in crystals:
                    self.batch.add(c.label, c.data, c.settings,
                                   c.group or self.group_edit.text().strip(), c.result)
            except Exception as exc:
                failed.append(f"{Path(path).name}: {exc}")
        self.refresh()
        if failed:
            QMessageBox.warning(self, "Some files were not added", "\n".join(failed))

    def remove_selected(self):
        for i in sorted(self._selected(), reverse=True):
            self.batch.crystals.pop(i)
            if self._in_main == i:
                self._in_main = None
        self.refresh()

    # --- the main window ---------------------------------------------------------------
    def open_in_main(self):
        rows = self._selected()
        if len(rows) != 1:
            self.status.setText("Select one crystal to open in the main window.")
            return
        i = rows[0]
        c = self.batch.crystals[i]
        crystal = cs.Crystal(c.label, self.batch.settings(i), c.data, c.group, c.result)
        if self.owner.open_crystal(crystal, c.label, refit=c.result is not None):
            self._in_main = i
            self.status.setText(f"{c.label} is open in the main window. Change it there, then "
                                "press Take back from main window.")
            self.owner.raise_()
            self.owner.activateWindow()

    def take_from_main(self):
        rows = self._selected()
        i = rows[0] if len(rows) == 1 else self._in_main
        if i is None:
            self.status.setText("Select the crystal that should take the main window's settings.")
            return
        data = self.owner.current_data()
        if data is None:
            return
        c = self.batch.crystals[i]
        if Path(data.source).name != Path(c.data.source).name:
            answer = QMessageBox.question(
                self, "A different profile",
                f"The main window shows {Path(data.source).name}, not {Path(c.data.source).name}. "
                f"Give {c.label} that profile and its settings?")
            if answer != QMessageBox.Yes:
                return
        c.data = data
        c.overrides = cs.differences(self.owner.current_settings(), self.batch.shared)
        c.result = None if self.owner._stale else cs.result_to_json(self.owner.fit_result,
                                                                   self.owner.mc_result)
        c.status = "taken from the main window" + (" with its result" if c.result else "")
        self.refresh()

    # --- fitting ---------------------------------------------------------------------
    def run(self, monte_carlo: bool):
        if self.thread is not None and self.thread.isRunning():
            return
        rows = self._selected() or list(range(len(self.batch)))
        jobs = [(i, self.batch.crystals[i].label, self.batch.settings(i), self.batch.crystals[i].data)
                for i in rows]
        if not jobs:
            return
        for i in rows:
            self.batch.crystals[i].status = "waiting"
        self.refresh()
        for b in (self.fit, self.fit_mc, self.add_files, self.add_current, self.add_saved,
                  self.remove, self.open_main, self.take_main, self.load, self.save):
            b.setEnabled(False)
        self.stop.setEnabled(True)
        cores = self.owner.sp_cores.value() if hasattr(self.owner, "sp_cores") else 1
        self.worker = BatchWorker(jobs, monte_carlo, cores)
        self.worker.progress.connect(self._progress)
        self.worker.crystal_done.connect(self._crystal_done)
        self.worker.crystal_failed.connect(self._crystal_failed)
        self.thread = start(self.worker)
        self.thread.finished.connect(self._idle)

    def abort(self):
        if self.worker is not None:
            self.worker.abort()
        self.status.setText("Stopping after the current crystal…")

    @Slot(int, int, str)
    def _progress(self, k, n, label):
        self.status.setText(f"Crystal {k + 1} of {n}: {label}")

    @Slot(int, object, object)
    def _crystal_done(self, i, fit, mc):
        if i >= len(self.batch):
            return
        c = self.batch.crystals[i]
        c.result = cs.result_to_json(fit, mc)
        c.status = (f"check: {len(fit.warnings)} warning(s)" if fit.warnings else "fitted")
        self.refresh()

    @Slot(int, str)
    def _crystal_failed(self, i, msg):
        if i >= len(self.batch):
            return
        c = self.batch.crystals[i]
        c.result = None
        c.status = "failed: " + msg.splitlines()[0][:120]
        self.refresh()

    @Slot()
    def _idle(self):
        for c in self.batch.crystals:
            if c.status == "waiting":
                c.status = "stopped"
        self.stop.setEnabled(False)
        for b in (self.add_files, self.add_current, self.add_saved, self.load):
            b.setEnabled(True)
        self.refresh()

    # --- out -------------------------------------------------------------------------
    def send_to_population(self):
        from ..population import Timescale
        records = []
        for i, c in enumerate(self.batch.crystals):
            r = c.result
            if not r:
                continue
            mc = r.get("montecarlo") or {}
            records.append(Timescale(
                label=c.label, group=c.group, t_best=r["t_seconds"],
                samples=np.asarray(mc["times_s"], dtype=float) if mc.get("times_s") else None,
                coefficient=r.get("coefficient", ""), T_C=self.batch.settings(i)["T_C"],
                source=c.data.source))
        if not records:
            self.status.setText("No crystal has a result yet. Fit them first.")
            return
        self.owner.show_population_study()
        self.owner.population_study._add(records)
        self.status.setText(f"{len(records)} timescales sent to the population window.")

    def save_batch(self):
        start = str(Path(self.owner.output_dir()) / "batch") if self.owner.output_dir() else "batch"
        path, _ = QFileDialog.getSaveFileName(self, "Save batch", start + cs.SUFFIX,
                                              "Diffusor setups (*.diffusor)")
        if not path:
            return
        try:
            written = self.batch.save(path)
            self.status.setText(f"Saved {len(self.batch)} crystals to {written.name}.")
        except Exception:
            QMessageBox.critical(self, "Could not save the batch", traceback.format_exc())

    def open_batch(self, path=None):
        if path is None:
            path, _ = QFileDialog.getOpenFileName(self, "Open batch", self.owner.output_dir(),
                                                  "Diffusor setups (*.diffusor)")
            if not path:
                return
        try:
            d = cs.read_file(path)
            if d["format"] != cs.FORMAT_BATCH:
                raise ValueError(f"{Path(path).name} is a single crystal; use Add saved setups")
            batch = cs.Batch.from_json(d)
        except Exception as exc:
            QMessageBox.warning(self, "Could not open the batch", str(exc))
            return
        if len(self.batch):
            answer = QMessageBox.question(self, "Replace the batch",
                                          f"Replace the {len(self.batch)} crystals in the table?")
            if answer != QMessageBox.Yes:
                return
        self.batch = batch
        self._in_main = None
        self.refresh()

    def closeEvent(self, event):
        if self.thread is not None and self.thread.isRunning():
            self.abort()
        super().closeEvent(event)
