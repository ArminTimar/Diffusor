"""Measured validation data browser, with explicit model-readiness labels."""
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QDialog,QVBoxLayout,QLineEdit,QListWidget,QListWidgetItem,
                              QLabel,QDialogButtonBox)


class ValidationLibraryDialog(QDialog):
    def __init__(self, records, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Published validation library')
        self.resize(740,570)
        self.selected=None
        layout=QVBoxLayout(self)
        self.search=QLineEdit()
        self.search.setPlaceholderText('Filter by study or crystal name')
        layout.addWidget(self.search)
        self.items=QListWidget()
        from ..datasets import BY_KEY
        for r in records:
            status='model preset' if r['key'] in BY_KEY else 'data only, setup required'
            item=QListWidgetItem(f"{r['study']} / {r['sample']} — {r.get('mineral','olivine')} — "
                                 f"{r['n_points']} points — {status}")
            item.setData(Qt.UserRole,r)
            self.items.addItem(item)
        layout.addWidget(self.items,1)
        self.details=QLabel('Choose a profile. Measured data and published model curves are separate columns.')
        self.details.setWordWrap(True)
        self.details.setTextFormat(Qt.PlainText)
        layout.addWidget(self.details)
        buttons=QDialogButtonBox(QDialogButtonBox.Open|QDialogButtonBox.Cancel)
        buttons.button(QDialogButtonBox.Open).setText('Load profile')
        buttons.button(QDialogButtonBox.Open).setEnabled(False)
        layout.addWidget(buttons)
        def changed(item):
            buttons.button(QDialogButtonBox.Open).setEnabled(item is not None)
            self.selected=item.data(Qt.UserRole) if item else None
            if self.selected:
                r=self.selected
                self.details.setText(f"Source: {r['source_file']} / {r['source_sheet']}\n{r['notes']}")
        def filtered(text):
            for i in range(self.items.count()):
                item=self.items.item(i)
                item.setHidden(text.casefold() not in item.text().casefold())
            if self.items.currentItem() and self.items.currentItem().isHidden():
                self.items.setCurrentRow(-1)
        self.items.currentItemChanged.connect(lambda item,previous:changed(item))
        self.search.textChanged.connect(filtered)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
