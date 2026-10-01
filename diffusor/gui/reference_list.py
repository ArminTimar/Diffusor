"""Help > All references: every source Diffusor cites, searchable, with BibTeX export.

The groups and the "Used by" lines are built from the code itself (the coefficient
registry, the buffer citations and the example datasets), so the window cannot fall
out of step with what the app actually uses. A paper cited for several minerals
appears under each of them.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Tuple

from PySide6.QtWidgets import (QDialog, QDialogButtonBox, QFileDialog, QLineEdit, QMessageBox,
                               QPushButton, QTextBrowser, QVBoxLayout)

from ..coefficients import list_coefficients
from ..datasets import DATASETS
from ..minerals.definitions import MINERALS
from ..references import REFERENCES, bibtex_document
from ..thermo.buffers import BUFFER_CITATIONS
from .richtext import Doc, esc, reference_items
from .widgets import fit_to_screen

OTHER_GROUP = "Methods, constants, reviews and applications"
BUFFER_GROUP = "Oxygen buffers"


def _sort_key(key: str):
    r = REFERENCES[key]
    return (r.authors.lower(), r.year, key)


def reference_usage() -> Dict[str, List[str]]:
    """For each citation key, what in Diffusor uses it."""
    used: Dict[str, List[str]] = {}
    for c in list_coefficients():
        used.setdefault(c.citation, []).append(c.key)
        for k in c.secondary_citations:
            used.setdefault(k, []).append(f"{c.key} (secondary)")
    for name, key in BUFFER_CITATIONS.items():
        used.setdefault(key, []).append(f"oxygen buffers ({name})")
    for d in DATASETS:
        if d.citation:
            used.setdefault(d.citation, []).append(f"example {d.key}")
    return {k: sorted(set(v)) for k, v in used.items()}


def reference_groups() -> List[Tuple[str, List[str]]]:
    """(group title, citation keys), minerals first, in the app's mineral order."""
    by_mineral: Dict[str, set] = {}
    for c in list_coefficients():
        by_mineral.setdefault(c.mineral, set()).update((c.citation, *c.secondary_citations))
    groups = [(MINERALS[m].name if m in MINERALS else m, sorted(by_mineral[m], key=_sort_key))
              for m in MINERALS if m in by_mineral]
    groups += [(m, sorted(keys, key=_sort_key)) for m, keys in by_mineral.items()
               if m not in MINERALS]
    groups.append((BUFFER_GROUP, sorted(set(BUFFER_CITATIONS.values()), key=_sort_key)))
    placed = {k for _, keys in groups for k in keys}
    groups.append((OTHER_GROUP, sorted(set(REFERENCES) - placed, key=_sort_key)))
    return [(title, keys) for title, keys in groups if keys]


def _matches(query: str, key: str, group: str, usage: List[str]) -> bool:
    """Every word of the query must start a word somewhere in the entry.

    Matching at word starts keeps "Ni" from hitting "units". Keys are searched
    both whole and split at underscores, so "ol_Ni_petry2004" and "Ni" both work.
    """
    if not query.strip():
        return True
    text = " ".join((key, REFERENCES[key].full(), group, *usage))
    hay = (text + " " + text.replace("_", " ")).lower()
    return all(re.search(r"(?<![a-z0-9])" + re.escape(w), hay) for w in query.lower().split())


def references_html(query: str = "") -> str:
    """The whole list, or only the entries matching every word of ``query``."""
    usage = reference_usage()
    groups = reference_groups()
    shown = {title: [k for k in keys if _matches(query, k, title, usage.get(k, []))]
             for title, keys in groups}
    n_total = len(REFERENCES)
    n_shown = len({k for keys in shown.values() for k in keys})
    if query.strip():
        subtitle = esc(f"{n_shown} of {n_total} references match “{query.strip()}”.")
    else:
        subtitle = esc(f"{n_total} references. Every coefficient, equation, constant and "
                       f"buffer in Diffusor cites one of these.")
    d = Doc("All references", subtitle)
    if not n_shown:
        d.p("No references match. Search looks in authors, titles, journals, DOIs, "
            "citation keys, minerals and the coefficients that use each paper.", "muted")
        return d.html()
    if not query.strip():
        d.p("Grouped by the mineral whose diffusion laws cite them. A paper cited for "
            "several minerals appears under each. Under each entry is what in Diffusor "
            "uses it.", "muted")
    for title, keys in groups:
        if not shown[title]:
            continue
        d.h(title)
        items = []
        for key, text in zip(shown[title], reference_items(shown[title])):
            line = f"{text}<br><span class='faint'>Key <span class='mono'>{esc(key)}</span>"
            if usage.get(key):
                line += " · used by " + esc(", ".join(usage[key]))
            items.append(line + "</span>")
        d.ul(items)
    return d.html()


class ReferenceListDialog(QDialog):
    """Searchable list of every reference, with a button to save them as BibTeX."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("All references")
        fit_to_screen(self, 920, 760)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(16, 16, 16, 12)
        lay.setSpacing(10)
        self.search = QLineEdit()
        self.search.setPlaceholderText("Search authors, titles, minerals, coefficient keys or DOIs")
        self.search.setClearButtonEnabled(True)
        self.search.textChanged.connect(self._refresh)
        lay.addWidget(self.search)
        self.view = QTextBrowser()
        self.view.setOpenExternalLinks(True)
        self.view.document().setDocumentMargin(22)
        lay.addWidget(self.view)
        bb = QDialogButtonBox(QDialogButtonBox.Close)
        self.btn_bib = QPushButton("Save as BibTeX...")
        self.btn_bib.clicked.connect(self.save_bibtex)
        bb.addButton(self.btn_bib, QDialogButtonBox.ActionRole)
        bb.rejected.connect(self.reject)
        lay.addWidget(bb)
        self._refresh()

    def _refresh(self, *_):
        self.view.setHtml(references_html(self.search.text()))

    def save_bibtex(self, path: str = "") -> str:
        """Write every reference to a .bib file. Returns the path, or "" if cancelled."""
        if not path:
            path, _ = QFileDialog.getSaveFileName(
                self, "Save references as BibTeX", str(Path.home() / "diffusor_references.bib"),
                "BibTeX (*.bib)")
        if not path:
            return ""
        try:
            Path(path).write_text(bibtex_document(), encoding="utf-8")
        except OSError as exc:
            QMessageBox.warning(self, "Could not save", f"{path}\n\n{exc}")
            return ""
        return path
