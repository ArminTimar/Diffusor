"""Measured validation data browser, with explicit model-readiness labels."""
import re

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QDialog, QVBoxLayout, QLineEdit, QListWidget, QListWidgetItem,
                               QLabel, QDialogButtonBox)

HINT = ("Pick a profile, then press Load profile or double-click it. "
        "\"Settings included\" also sets the study's mineral, conditions and model; "
        "\"measurements only\" loads the data and leaves the settings to you.")


def study_name(key: str) -> str:
    """A study's key as a reader would write it: "Iovine et al. (2017)", or, for a study
    without an entry in the reference list, "Sundermeyer 2020, Eifel". The second is a
    label, not a citation: the key does not say whether the paper had co-authors."""
    from ..references import REFERENCES, cite
    if key in REFERENCES:
        return cite(key)
    m = re.fullmatch(r"([a-z]+)(?:_([a-z]+))?(\d{4})(?:_([a-z]+))?", key)
    if not m:
        return key
    author, tag, year, place = m.groups()
    name = f"{author.capitalize()} {year}"
    for extra in (tag, place):
        if extra:
            name += f", {extra.upper() if len(extra) <= 3 else extra.capitalize()}"
    return name


class ValidationLibraryDialog(QDialog):
    def __init__(self, records, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Published validation library')
        self.resize(740, 570)
        self.selected = None
        layout = QVBoxLayout(self)
        self.search = QLineEdit()
        self.search.setPlaceholderText('Filter by study, crystal or mineral')
        layout.addWidget(self.search)
        self.items = QListWidget()
        from ..datasets import BY_KEY
        for r in records:
            status = 'settings included' if r['key'] in BY_KEY else 'measurements only'
            item = QListWidgetItem(f"{study_name(r['study'])}  ·  {r['sample']}  ·  "
                                   f"{r.get('mineral', 'olivine')}  ·  {r['n_points']} points  ·  "
                                   f"{status}")
            item.setData(Qt.UserRole, r)
            self.items.addItem(item)
        layout.addWidget(self.items, 1)
        self.details = QLabel(HINT)
        self.details.setWordWrap(True)
        self.details.setTextFormat(Qt.PlainText)
        layout.addWidget(self.details)
        buttons = QDialogButtonBox(QDialogButtonBox.Open | QDialogButtonBox.Cancel)
        buttons.button(QDialogButtonBox.Open).setText('Load profile')
        buttons.button(QDialogButtonBox.Open).setEnabled(False)
        buttons.button(QDialogButtonBox.Open).setObjectName("Primary")
        layout.addWidget(buttons)

        def changed(item):
            buttons.button(QDialogButtonBox.Open).setEnabled(item is not None)
            self.selected = item.data(Qt.UserRole) if item else None
            if self.selected:
                r = self.selected
                self.details.setText(f"Source: {r['source_file']} / {r['source_sheet']}\n"
                                     f"{r['notes']}")
            else:
                self.details.setText(HINT)

        def filtered(text):
            for i in range(self.items.count()):
                item = self.items.item(i)
                item.setHidden(text.casefold() not in item.text().casefold())
            if self.items.currentItem() and self.items.currentItem().isHidden():
                self.items.setCurrentRow(-1)

        def double_clicked(item):
            self.items.setCurrentItem(item)
            self.accept()

        self.items.currentItemChanged.connect(lambda item, previous: changed(item))
        self.items.itemDoubleClicked.connect(double_clicked)
        self.search.textChanged.connect(filtered)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
