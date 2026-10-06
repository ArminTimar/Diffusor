"""The examples list as a tree: mineral, then study, then example.

Twenty-odd examples from one study would bury the others in a flat list, so the
list is grouped, and every mineral and study starts closed. Only the last level is
an example; the mineral and study rows open and close and cannot be loaded.

The tree keeps the small surface a flat list had, so the rest of the window and
its tests can still ask for the n-th example: ``count()``, ``item(i)``,
``setCurrentRow(i)``, and items that take no column in ``data``, ``setText`` and
``setFont``. Those all mean the examples, in the order they appear.
"""
from __future__ import annotations

from typing import List, Optional, Sequence, Tuple

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QAbstractItemView, QTreeWidget, QTreeWidgetItem

KEY_ROLE = Qt.UserRole

# (key, label, enabled)
Example = Tuple[str, str, bool]
# (study label, examples)
Study = Tuple[str, Sequence[Example]]
# (mineral label, studies)
Mineral = Tuple[str, Sequence[Study]]


class ExampleItem(QTreeWidgetItem):
    """One example. Reads like a list item: no column to give."""

    def data(self, a, b=None):
        # Qt itself calls data(column, role); callers of a list item give only the role
        return super().data(0, a) if b is None else super().data(a, b)

    def setText(self, a, b=None):
        if b is None:
            super().setText(0, a)
        else:
            super().setText(a, b)

    def font(self, column=0):
        return super().font(column)

    def setFont(self, a, b=None):
        if b is None:
            super().setFont(0, a)
        else:
            super().setFont(a, b)


class ExampleTree(QTreeWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setHeaderHidden(True)
        self.setIndentation(16)
        self.setMouseTracking(True)
        self.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.setSelectionMode(QAbstractItemView.SingleSelection)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.setTextElideMode(Qt.ElideRight)
        self._leaves: List[ExampleItem] = []

    # -- building -------------------------------------------------------------
    def populate(self, minerals: Sequence[Mineral]) -> None:
        self.clear()
        self._leaves = []
        for mineral, studies in minerals:
            n = sum(len(examples) for _, examples in studies)
            top = self._group(self.invisibleRootItem(), f"{mineral}   ({n})")
            top.setExpanded(False)
            for study, examples in studies:
                node = self._group(top, f"{study}   ({len(examples)})")
                for key, label, enabled in examples:
                    leaf = ExampleItem(node)
                    leaf.setText(label)
                    leaf.setData(0, KEY_ROLE, key)
                    if not enabled:
                        leaf.setFlags(leaf.flags() & ~Qt.ItemIsEnabled)
                    self._leaves.append(leaf)
                node.setExpanded(False)

    def _group(self, parent, text: str) -> QTreeWidgetItem:
        item = QTreeWidgetItem(parent)
        item.setText(0, text)
        f = item.font(0)
        f.setBold(True)
        item.setFont(0, f)
        item.setFlags(Qt.ItemIsEnabled)             # opens and closes, cannot be selected
        return item

    # -- the flat view of the examples ---------------------------------------
    def count(self) -> int:
        return len(self._leaves)

    def item(self, i: int) -> Optional[ExampleItem]:
        return self._leaves[i] if 0 <= i < len(self._leaves) else None

    def setCurrentRow(self, i: int) -> None:
        leaf = self.item(i)
        if leaf is not None:
            self.reveal(leaf)
            self.setCurrentItem(leaf)

    def reveal(self, leaf: QTreeWidgetItem) -> None:
        """Open the study and mineral above an example and scroll it into view."""
        parent = leaf.parent()
        while parent is not None:
            parent.setExpanded(True)
            parent = parent.parent()
        self.scrollToItem(leaf)

    @staticmethod
    def key_of(item) -> Optional[str]:
        """The dataset key of an item, or None for a mineral or study row."""
        return None if item is None else item.data(0, KEY_ROLE)
