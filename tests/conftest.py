"""Shared test setup.

Interface tests build many windows in one process. A window that is only
closed stays alive until Python's garbage collector destroys it at some later,
arbitrary moment, possibly while Qt still has events queued for its children,
which crashed the suite on Windows now and then (heap corruption, 0xc0000374).
After each test every leftover top-level window is closed and deleted through
Qt's own event loop instead, so destruction happens at a known point.
"""
import gc

import pytest


@pytest.fixture(autouse=True)
def _delete_leftover_windows():
    yield
    import sys
    if "PySide6.QtWidgets" not in sys.modules:
        return
    from PySide6.QtCore import QCoreApplication, QEvent
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance()
    if app is None:
        return
    for w in QApplication.topLevelWidgets():
        w.close()
        w.deleteLater()
    QCoreApplication.sendPostedEvents(None, QEvent.DeferredDelete)
    app.processEvents()
    gc.collect()
