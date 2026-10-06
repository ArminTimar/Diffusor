"""The examples list: mineral, then study, then example, and no 'measured' label."""
import os

import pytest
from PySide6.QtCore import Qt

from diffusor import datasets as ds


@pytest.fixture(scope="module")
def app():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    import matplotlib
    matplotlib.use("QtAgg")
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def test_every_example_is_listed_once_under_its_own_mineral_and_study():
    from diffusor.gui.main_window import MainWindow
    groups = MainWindow._example_groups()
    seen = []
    for mineral, studies in groups:
        for study, examples in studies:
            for key, label, _enabled in examples:
                d = ds.get(key)
                assert d.study_name == study and label.startswith(d.short_name)
                seen.append(key)
    assert sorted(seen) == sorted(d.key for d in ds.DATASETS)
    names = [m for m, _ in groups]
    assert names.index("Olivine") < names.index("Orthopyroxene"), "the order follows the mineral registry"
    for _, studies in groups:
        order = [s for s, _ in studies]
        assert order == sorted(order, key=lambda s: (s == ds.SYNTHETIC_STUDY, s.lower()))


def test_everything_starts_closed_and_choosing_an_example_opens_its_path(app):
    from diffusor.gui.example_tree import ExampleTree
    from diffusor.gui.main_window import MainWindow
    tree = ExampleTree()
    tree.populate(MainWindow._example_groups())
    assert tree.count() == len(ds.DATASETS)
    for m in range(tree.topLevelItemCount()):
        mineral = tree.topLevelItem(m)
        assert not mineral.isExpanded(), mineral.text(0)
        for s in range(mineral.childCount()):
            assert not mineral.child(s).isExpanded(), mineral.child(s).text(0)
    lynn = next(i for i in range(tree.count()) if tree.item(i).data(Qt.UserRole).startswith("lynn2024"))
    leaf = tree.item(lynn)
    tree.setCurrentRow(lynn)
    assert leaf.parent().isExpanded() and leaf.parent().parent().isExpanded()
    assert tree.currentItem() is leaf, "choosing one opens its study and its mineral"


def test_mineral_and_study_rows_cannot_be_chosen_or_loaded(app):
    from diffusor.gui.example_tree import ExampleTree
    from diffusor.gui.main_window import MainWindow
    tree = ExampleTree()
    tree.populate(MainWindow._example_groups())
    top = tree.topLevelItem(0)
    study = top.child(0)
    for row in (top, study):
        assert not (row.flags() & Qt.ItemIsSelectable)
        assert ExampleTree.key_of(row) is None
    assert ExampleTree.key_of(tree.item(0)) == tree.item(0).data(Qt.UserRole)


def test_the_window_shows_no_measured_label(app):
    from diffusor.gui import richtext
    from diffusor.gui.main_window import MainWindow
    w = MainWindow()
    texts = [w.lst_examples.item(i).text(0) for i in range(w.lst_examples.count())]
    root = w.lst_examples.invisibleRootItem()
    stack = [root.child(i) for i in range(root.childCount())]
    while stack:
        it = stack.pop()
        texts.append(it.text(0))
        stack += [it.child(i) for i in range(it.childCount())]
    assert texts and not any("measured" in t.lower() for t in texts)
    for d in ds.DATASETS:
        assert "(measured)" not in d.name
        head = richtext.example_html(d)[:900].lower()
        if d.kind != "synthetic":
            assert "measured" not in head, d.key
    w.close()
