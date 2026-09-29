"""A nonmodal workspace for independent profiles sharing one duration."""
import copy
import traceback
from PySide6.QtCore import QObject, Signal, Slot
from PySide6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QPushButton,
    QLabel, QListWidget, QMessageBox, QFileDialog, QScrollArea)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from ..fitting.joint import ProfileConstraint, fit_joint_time
from ..dataio.joint import save_joint_results
from ..dataio.figures import joint_figure
from ..thermo.units import human_time
from .widgets import fit_to_screen
from .workers import start


class JointWorker(QObject):
    finished = Signal(object)
    failed = Signal(str)
    def __init__(self, profiles):
        super().__init__()
        self.profiles = profiles
        self.cancelled = False
    def abort(self):
        self.cancelled = True
    def run(self):
        try:
            self.finished.emit(fit_joint_time(self.profiles, progress=lambda _: self.cancelled))
        except Exception:
            self.failed.emit(traceback.format_exc())


class JointStudyDialog(QDialog):
    def __init__(self, owner):
        super().__init__(owner)
        self.owner = owner
        self.profiles = []
        self.result = None
        self.thread = None
        self.worker = None
        self.setWindowTitle("Shared-duration study")
        fit_to_screen(self, 1040, 780)
        layout = QVBoxLayout(self)
        info = QLabel("Configure a profile in the main window, then add it here. Repeat for each element or crystal. "
                      "Each profile needs measurement uncertainties. Only duration is fitted; other settings stay fixed. "
                      "This study combines independent scalar models and assumes independent measurement errors.")
        info.setWordWrap(True)
        layout.addWidget(info)
        self.list = QListWidget()
        self.list.setMaximumHeight(125)
        layout.addWidget(self.list)
        buttons = QHBoxLayout()
        self.add = QPushButton("Add current profile")
        self.remove = QPushButton("Remove selected")
        self.fit = QPushButton("Fit shared duration")
        self.stop = QPushButton("Stop")
        self.export = QPushButton("Export study…")
        self.stop.setEnabled(False)
        self.export.setEnabled(False)
        for button in (self.add, self.remove, self.fit, self.stop, self.export):
            buttons.addWidget(button)
        layout.addLayout(buttons)
        self.status = QLabel("No profiles added.")
        self.status.setWordWrap(True)
        layout.addWidget(self.status)
        self.scroll = QScrollArea()
        self.scroll.setWidgetResizable(True)
        layout.addWidget(self.scroll, 1)
        self.add.clicked.connect(self.add_current)
        self.remove.clicked.connect(self.remove_selected)
        self.fit.clicked.connect(self.run_fit)
        self.stop.clicked.connect(self.abort)
        self.export.clicked.connect(self.export_study)

    def _invalidate(self):
        self.result = None
        self.export.setEnabled(False)
        self.status.setText(f"{len(self.profiles)} profiles in this study. Fit to update results.")
        old = self.scroll.takeWidget()
        if old is not None:
            old.deleteLater()

    def add_current(self):
        if not self.owner._ready():
            return
        try:
            p = self.owner.profile
            model = copy.deepcopy(self.owner._model(self.owner._checked_keys()[0]))
            index = 1
            names = {c.name for c in self.profiles}
            while f"Profile {index}: {model.coefficient.mineral} {model.coefficient.species}" in names:
                index += 1
            name = f"Profile {index}: {model.coefficient.mineral} {model.coefficient.species}"
            constraint = ProfileConstraint(name, model, p.x, p.C, p.sigma)
            self.profiles.append(constraint)
            self.list.addItem(f"{name} · {len(p.x)} points · {model.coefficient.key}")
            self._invalidate()
        except Exception as exc:
            QMessageBox.warning(self, "Cannot add profile", str(exc))

    def remove_selected(self):
        i = self.list.currentRow()
        if i >= 0:
            self.profiles.pop(i)
            self.list.takeItem(i)
            self._invalidate()

    def run_fit(self):
        if len(self.profiles) < 2:
            self.status.setText("Add at least two independent profiles first.")
            return
        if self.thread is not None and self.thread.isRunning():
            return
        self._invalidate()
        for b in (self.add, self.remove, self.fit):
            b.setEnabled(False)
        self.stop.setEnabled(True)
        self.status.setText("Fitting the shared duration…")
        self.worker = JointWorker(copy.deepcopy(self.profiles))
        self.worker.finished.connect(self._done)
        self.worker.failed.connect(self._failed)
        self.thread = start(self.worker)
        self.thread.finished.connect(self._idle)

    @Slot()
    def _idle(self):
        for b in (self.add, self.remove, self.fit):
            b.setEnabled(True)
        self.stop.setEnabled(False)

    @Slot(object)
    def _done(self, result):
        self.result = result
        canvas = FigureCanvasQTAgg(joint_figure(result))
        canvas.setMinimumHeight(300*len(result.profiles))
        self.scroll.setWidget(canvas)
        self.export.setEnabled(True)
        self.status.setText(f"Shared duration: {human_time(result.t_seconds)} · χ²/dof = {result.reduced_chi2:.3g}. "
                            "Inspect each residual panel. This result holds conditions and coefficients fixed.")

    @Slot(str)
    def _failed(self, message):
        self.status.setText(message.splitlines()[-1])

    def abort(self):
        if self.worker is not None:
            self.worker.abort()

    def closeEvent(self, event):
        if self.thread is not None and self.thread.isRunning():
            self.abort()
            self.status.setText("Cancelling. Close this window once the current forward calculation finishes.")
            event.ignore()
        else:
            super().closeEvent(event)

    def export_study(self):
        if self.result is None:
            return
        directory = QFileDialog.getExistingDirectory(self, "Export shared-duration study")
        if directory:
            try:
                save_joint_results(directory, self.result)
                self.status.setText(f"Study exported to {directory}")
            except Exception as exc:
                QMessageBox.critical(self, "Export failed", str(exc))
