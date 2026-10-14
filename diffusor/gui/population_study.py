"""A nonmodal workspace that plots many crystals' timescales together."""
import traceback
from pathlib import Path

from PySide6.QtCore import Qt, Slot
from PySide6.QtWidgets import (QAbstractItemView, QCheckBox, QComboBox, QDialog, QFileDialog,
                               QHBoxLayout, QHeaderView, QLabel, QLineEdit, QMessageBox,
                               QPushButton, QScrollArea, QSplitter, QTableWidget,
                               QTableWidgetItem, QVBoxLayout, QWidget)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg

from ..dataio.figures import GROUP_COLOURS, population_figure
from ..population import (UNIT_NAMES, extend_unique, from_run, group_name, groups_of,
                          load_result_json, load_setup_file, load_table,
                          save_population)
from ..thermo.units import human_time
from .widgets import fit_to_screen

COLUMNS = ("Label", "Group", "Best-fit time", "68% interval", "Uncertainty")
UNITS = ("Automatic", "h", "d", "yr", "kyr")


class PopulationStudyDialog(QDialog):
    def __init__(self, owner):
        super().__init__(owner)
        self.owner = owner
        self.records = []
        self.figure = None
        self.setWindowTitle("Timescale population")
        fit_to_screen(self, 1220, 780)
        layout = QVBoxLayout(self)
        info = QLabel("Collect the timescales of many crystals and compare them in one figure: the "
                      "times ranked, a summed probability density per group, and the cumulative "
                      "distribution. Add the run in the main window, saved crystal setups or batches "
                      "(*.diffusor) and exported results (*_results.json), or a table of times. Groups, for "
                      "example eruptions or samples, are compared in colour; edit a label or group "
                      "in the table.")
        info.setWordWrap(True)
        layout.addWidget(info)

        split = QSplitter(Qt.Horizontal)
        layout.addWidget(split, 1)

        left = QWidget()
        lv = QVBoxLayout(left)
        lv.setContentsMargins(0, 0, 0, 0)
        row = QHBoxLayout()
        row.addWidget(QLabel("Group for new entries:"))
        self.group_edit = QLineEdit()
        self.group_edit.setPlaceholderText("for example 2011 sub-Plinian")
        row.addWidget(self.group_edit, 1)
        lv.addLayout(row)
        row = QHBoxLayout()
        self.add_run = QPushButton("Add current run")
        self.add_files = QPushButton("Add saved results…")
        self.add_table = QPushButton("Add table…")
        for b in (self.add_run, self.add_files, self.add_table):
            row.addWidget(b)
        lv.addLayout(row)
        self.table = QTableWidget(0, len(COLUMNS))
        self.table.setHorizontalHeaderLabels(COLUMNS)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.verticalHeader().setVisible(False)
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.Interactive)
        header.setStretchLastSection(True)
        lv.addWidget(self.table, 1)
        row = QHBoxLayout()
        self.set_group = QPushButton("Set group of selected")
        self.remove = QPushButton("Remove selected")
        self.clear = QPushButton("Clear")
        for b in (self.set_group, self.remove, self.clear):
            row.addWidget(b)
        lv.addLayout(row)
        split.addWidget(left)

        right = QWidget()
        rv = QVBoxLayout(right)
        rv.setContentsMargins(0, 0, 0, 0)
        row = QHBoxLayout()
        self.panel_boxes = {}
        for key, text in (("ranked", "Ranked"), ("density", "Density"), ("cumulative", "Cumulative")):
            box = QCheckBox(text)
            box.setChecked(True)
            box.toggled.connect(self.redraw)
            self.panel_boxes[key] = box
            row.addWidget(box)
        row.addSpacing(12)
        row.addWidget(QLabel("Time axis:"))
        self.unit = QComboBox()
        for u in UNITS:
            self.unit.addItem(u if u == "Automatic" else UNIT_NAMES[u], u)
        self.unit.currentIndexChanged.connect(self.redraw)
        row.addWidget(self.unit)
        row.addStretch(1)
        self.export = QPushButton("Export…")
        self.export.setEnabled(False)
        row.addWidget(self.export)
        rv.addLayout(row)
        self.status = QLabel("No timescales added.")
        self.status.setWordWrap(True)
        rv.addWidget(self.status)
        self.scroll = QScrollArea()
        self.scroll.setWidgetResizable(True)
        rv.addWidget(self.scroll, 1)
        split.addWidget(right)
        split.setSizes([480, 740])

        self.add_run.clicked.connect(self.add_current)
        self.add_files.clicked.connect(self.add_result_files)
        self.add_table.clicked.connect(self.add_table_file)
        self.set_group.clicked.connect(self.set_selected_group)
        self.remove.clicked.connect(self.remove_selected)
        self.clear.clicked.connect(self.clear_all)
        self.export.clicked.connect(self.export_population)
        self.table.itemChanged.connect(self._edited)

    # --- adding --------------------------------------------------------------------
    def _group(self) -> str:
        return self.group_edit.text().strip()

    def _add(self, new):
        added = extend_unique(self.records, new)
        self._fill()
        self.redraw()
        return added

    def add_current(self):
        o = self.owner
        if o.fit_result is None or getattr(o, "_stale", False):
            QMessageBox.information(self, "No result", "Fit a profile in the main window first. "
                                    "A Monte Carlo run adds its uncertainty.")
            return
        source = o.profile.source if o.profile is not None and o.profile.source else ""
        label = Path(source).stem if source else "Run"
        try:
            self._add([from_run(o.fit_result, o.mc_result, label, self._group(), source)])
        except Exception as exc:
            QMessageBox.warning(self, "Cannot add the run", str(exc))

    def add_result_files(self):
        paths, _ = QFileDialog.getOpenFileNames(
            self, "Add Diffusor results", self.owner.output_dir(),
            "Saved setups and exported results (*.diffusor *_results.json);;All files (*)")
        new, failed = [], []
        for p in paths:
            try:
                if p.lower().endswith(".diffusor"):
                    new.extend(load_setup_file(p, self._group()))
                else:
                    r = load_result_json(p)
                    r.group = self._group()
                    new.append(r)
            except Exception as exc:
                failed.append(f"{Path(p).name}: {exc}")
        if new:
            self._add(new)
        if failed:
            QMessageBox.warning(self, "Some files were not added", "\n".join(failed))

    def add_table_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Add a table of timescales", self.owner.output_dir(),
            "Tables (*.csv *.txt *.tsv *.xlsx *.xls)")
        if not path:
            return
        unit = self._ask_unit()
        if unit is None:
            return
        try:
            self._add(load_table(path, default_unit=unit, group=self._group()))
        except Exception as exc:
            QMessageBox.warning(self, "Cannot read the table", str(exc))

    def _ask_unit(self):
        """The unit of time columns whose header does not name one."""
        from PySide6.QtWidgets import QInputDialog
        choices = [UNIT_NAMES[u] for u in ("h", "d", "yr", "kyr")]
        text, ok = QInputDialog.getItem(
            self, "Time unit", "Unit of time columns whose name gives none\n"
            "(a header such as 'time (days)' or 'p16_yr' overrides this):", choices, 1, False)
        if not ok:
            return None
        return {UNIT_NAMES[u]: u for u in ("h", "d", "yr", "kyr")}[text]

    # --- editing -------------------------------------------------------------------
    def _selected_rows(self):
        return sorted({i.row() for i in self.table.selectedIndexes()})

    def set_selected_group(self):
        rows = self._selected_rows()
        if not rows:
            self.status.setText("Select rows in the table first.")
            return
        for i in rows:
            self.records[i].group = self._group()
        self._fill()
        self.redraw()

    def remove_selected(self):
        for i in reversed(self._selected_rows()):
            self.records.pop(i)
        self._fill()
        self.redraw()

    def clear_all(self):
        self.records.clear()
        self._fill()
        self.redraw()

    @Slot(object)
    def _edited(self, item):
        r = self.records[item.row()]
        text = item.text().strip()
        if item.column() == 0:
            if text and text not in {x.label for x in self.records if x is not r}:
                r.label = text
            else:
                self._fill()          # an empty or repeated label is put back
                return
        elif item.column() == 1:
            r.group = "" if text == group_name("") else text
        self.redraw()

    def _fill(self):
        self.table.blockSignals(True)
        self.table.setRowCount(len(self.records))
        for i, r in enumerate(self.records):
            interval = (f"{human_time(r.p16)} – {human_time(r.p84)}" if r.has_interval else "")
            basis = ("Monte Carlo" if r.samples is not None and r.samples.size >= 20
                     else "68% interval" if r.has_interval else "none")
            for j, text in enumerate((r.label, group_name(r.group), human_time(r.t_best), interval, basis)):
                item = QTableWidgetItem(text)
                if j > 1:
                    item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                if j == 4 and r.notes:
                    item.setToolTip("\n".join(r.notes))
                if j == 0 and r.source:
                    item.setToolTip(r.source)
                self.table.setItem(i, j, item)
        self.table.resizeColumnsToContents()
        self.table.blockSignals(False)

    # --- figure --------------------------------------------------------------------
    @Slot()
    def redraw(self):
        old = self.scroll.takeWidget()
        if old is not None:
            old.deleteLater()
        self.figure = None
        self.export.setEnabled(False)
        if not self.records:
            self.status.setText("No timescales added.")
            return
        panels = [k for k, b in self.panel_boxes.items() if b.isChecked()]
        unit = self.unit.currentData()
        try:
            self.figure = population_figure(self.records, panels=panels,
                                            unit=None if unit == "Automatic" else unit)
        except Exception as exc:
            self.status.setText(str(exc))
            return
        canvas = FigureCanvasQTAgg(self.figure)
        w, h = self.figure.get_size_inches()
        canvas.setMinimumHeight(int(h * 72))
        self.scroll.setWidget(canvas)
        self.export.setEnabled(True)
        groups = groups_of(self.records)
        self.status.setText(f"{len(self.records)} timescales in {len(groups)} group(s)"
                            + (f"; at most {len(GROUP_COLOURS)} groups can be drawn." if len(groups) > len(GROUP_COLOURS) else "."))

    def export_population(self):
        if self.figure is None:
            return
        d = QFileDialog.getExistingDirectory(self, "Export the population into folder", self.owner.output_dir())
        if not d:
            return
        try:
            written = save_population(d, self.records, self.figure)
            QMessageBox.information(self, "Exported", "Written:\n" + "\n".join(written.values()))
        except Exception:
            QMessageBox.critical(self, "Export failed", traceback.format_exc())
