"""File > Export results: choose which outputs to write, before choosing the folder."""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QCheckBox, QDialog, QGridLayout, QHBoxLayout, QLabel, QPushButton,
                               QVBoxLayout)

from ..dataio.export import EXPORT_ITEMS
from .widgets import divider, ghost_button, note, primary_button


class ExportDialog(QDialog):
    """One tick box per output, with a few words beside it on what it holds.

    ``chosen`` is what the user ticked last time; None ticks everything that is
    available. The Monte Carlo times are listed but cannot be ticked until a Monte
    Carlo has been run, so the dialog always shows everything an export can hold.
    """

    def __init__(self, parent=None, has_monte_carlo: bool = False,
                 chosen: Optional[Sequence[str]] = None):
        super().__init__(parent)
        self.setWindowTitle("Export results")
        self.setModal(True)
        self.setMinimumWidth(600)
        self._boxes: Dict[str, QCheckBox] = {}

        lay = QVBoxLayout(self)
        lay.setContentsMargins(20, 18, 20, 16)
        lay.setSpacing(10)
        head = QLabel("What do you want to export?")
        head.setObjectName("H2")
        lay.addWidget(head)
        lay.addWidget(note("Tick the outputs to write. You choose the folder next.", "Hint"))

        quick = QHBoxLayout()
        all_ = ghost_button("Select all")
        none = ghost_button("Select none")
        all_.clicked.connect(lambda: self._set_all(True))
        none.clicked.connect(lambda: self._set_all(False))
        quick.addWidget(all_)
        quick.addWidget(none)
        quick.addStretch(1)
        lay.addLayout(quick)
        lay.addWidget(divider())

        grid = QGridLayout()
        grid.setHorizontalSpacing(14)
        grid.setVerticalSpacing(8)
        for row, (key, label, text, needs_mc) in enumerate(EXPORT_ITEMS):
            available = has_monte_carlo or not needs_mc
            box = QCheckBox(label)
            box.setEnabled(available)
            box.setChecked(available and (chosen is None or key in chosen))
            box.toggled.connect(self._update_buttons)
            self._boxes[key] = box
            detail = QLabel(text if available else text + ". Run a Monte Carlo first")
            detail.setObjectName("Hint")
            grid.addWidget(box, row, 0)
            grid.addWidget(detail, row, 1)
        grid.setColumnStretch(1, 1)
        lay.addLayout(grid)

        lay.addWidget(divider())
        buttons = QHBoxLayout()
        buttons.addStretch(1)
        cancel = QPushButton("Cancel")
        cancel.clicked.connect(self.reject)
        self.btn_ok = primary_button("Choose folder...")
        self.btn_ok.setDefault(True)
        self.btn_ok.clicked.connect(self.accept)
        buttons.addWidget(cancel)
        buttons.addWidget(self.btn_ok)
        lay.addLayout(buttons)
        self._update_buttons()

    def _set_all(self, on: bool) -> None:
        for box in self._boxes.values():
            if box.isEnabled():
                box.setChecked(on)

    def _update_buttons(self, *_):
        self.btn_ok.setEnabled(bool(self.selected()))

    def selected(self) -> List[str]:
        """The keys of the ticked outputs, in the order they are listed."""
        return [k for k, b in self._boxes.items() if b.isEnabled() and b.isChecked()]
