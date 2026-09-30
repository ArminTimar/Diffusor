"""Entry point for the desktop application: ``python -m diffusor``."""
from __future__ import annotations

import sys
from pathlib import Path

ICON = Path(__file__).resolve().parent / "icons" / "diffusor.png"


def _own_taskbar_entry() -> None:
    """On Windows, group the window under Diffusor's icon, not pythonw's."""
    if sys.platform == "win32":
        try:
            import ctypes
            ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("Diffusor.Diffusor")
        except (AttributeError, OSError):
            pass


def main(argv=None) -> int:
    from PySide6.QtCore import QTimer
    from PySide6.QtGui import QIcon
    from PySide6.QtWidgets import QApplication
    import matplotlib
    matplotlib.use("QtAgg")
    from .main_window import MainWindow

    _own_taskbar_entry()
    app = QApplication(argv if argv is not None else sys.argv)
    app.setApplicationName("Diffusor")
    app.setWindowIcon(QIcon(str(ICON)))
    win = MainWindow()
    win.show()
    # after the window is up, so a slow network never delays the start
    QTimer.singleShot(1500, win.startup_update_check)
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
