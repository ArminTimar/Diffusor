"""Entry point for the desktop application: ``python -m diffusor``."""
from __future__ import annotations

import sys
from pathlib import Path

ICONS = Path(__file__).resolve().parent / "icons"
ICON = ICONS / "diffusor.png"
ICON_WIN = ICONS / "diffusor.ico"
APP_ID = "Diffusor.Diffusor"


def _own_taskbar_entry() -> None:
    """On Windows, group the window under Diffusor's icon, not pythonw's."""
    if sys.platform == "win32":
        try:
            import ctypes
            ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(APP_ID)
        except (AttributeError, OSError):
            pass


def _taskbar_icon(window) -> None:
    """Give the window's taskbar button Diffusor's icon on Windows.

    A window with its own application ID gets a plain default icon unless the
    shell is told which icon belongs to it. Two things do that: the icon sent to
    the window itself, and the relaunch properties, which the taskbar also uses
    when the button is pinned. Anything that goes wrong here leaves the window
    with the icon Qt already set.
    """
    if sys.platform != "win32":
        return
    try:
        import ctypes
        from ctypes import wintypes
        hwnd = int(window.winId())
        user32 = ctypes.windll.user32
        user32.LoadImageW.restype = wintypes.HANDLE
        user32.LoadImageW.argtypes = [wintypes.HINSTANCE, wintypes.LPCWSTR, wintypes.UINT,
                                      ctypes.c_int, ctypes.c_int, wintypes.UINT]
        user32.SendMessageW.argtypes = [wintypes.HWND, wintypes.UINT, wintypes.WPARAM,
                                        wintypes.HANDLE]
        user32.GetSystemMetrics.restype = ctypes.c_int
        IMAGE_ICON, LR_LOADFROMFILE, WM_SETICON = 1, 0x10, 0x80
        for which, metric in ((0, (49, 50)), (1, (11, 12))):    # SM_CXSMICON, SM_CXICON
            w, h = user32.GetSystemMetrics(metric[0]), user32.GetSystemMetrics(metric[1])
            icon = user32.LoadImageW(None, str(ICON_WIN), IMAGE_ICON, w, h, LR_LOADFROMFILE)
            if icon:
                user32.SendMessageW(hwnd, WM_SETICON, which, icon)
        _set_window_properties(hwnd, {
            5: APP_ID,                                                  # AppUserModel.ID
            2: f'"{sys.executable}" -m diffusor',                       # RelaunchCommand
            3: f"{ICON_WIN},0",                                         # RelaunchIconResource
            4: "Diffusor"})                                             # RelaunchDisplayName
    except Exception:
        pass


def _set_window_properties(hwnd: int, values: dict) -> None:
    """Set string properties of the System.AppUserModel group on a window."""
    import ctypes
    from ctypes import wintypes

    class GUID(ctypes.Structure):
        _fields_ = [("a", wintypes.DWORD), ("b", wintypes.WORD), ("c", wintypes.WORD),
                    ("d", ctypes.c_ubyte * 8)]

    def guid(text: str) -> GUID:
        g = GUID()
        ctypes.oledll.ole32.CLSIDFromString(text, ctypes.byref(g))
        return g

    class PropKey(ctypes.Structure):
        _fields_ = [("fmtid", GUID), ("pid", wintypes.DWORD)]

    class PropVariant(ctypes.Structure):
        _fields_ = [("vt", wintypes.WORD), ("r1", wintypes.WORD), ("r2", wintypes.WORD),
                    ("r3", wintypes.WORD), ("ptr", ctypes.c_void_p), ("pad", ctypes.c_void_p)]

    VT_LPWSTR = 31
    shell32 = ctypes.windll.shell32
    store = ctypes.c_void_p()
    shell32.SHGetPropertyStoreForWindow.argtypes = [wintypes.HWND, ctypes.POINTER(GUID),
                                                    ctypes.POINTER(ctypes.c_void_p)]
    shell32.SHGetPropertyStoreForWindow(
        hwnd, ctypes.byref(guid("{886D8EEB-8CF2-4446-8D02-CDBA1DBDCF99}")),   # IPropertyStore
        ctypes.byref(store))
    if not store:
        return
    vtable = ctypes.cast(ctypes.cast(store, ctypes.POINTER(ctypes.c_void_p))[0],
                         ctypes.POINTER(ctypes.c_void_p))
    proto = ctypes.WINFUNCTYPE
    set_value = proto(ctypes.HRESULT, ctypes.c_void_p, ctypes.POINTER(PropKey),
                      ctypes.POINTER(PropVariant))(vtable[6])
    commit = proto(ctypes.HRESULT, ctypes.c_void_p)(vtable[7])
    release = proto(ctypes.c_ulong, ctypes.c_void_p)(vtable[2])
    fmtid = guid("{9F4C2855-9F79-4B39-A8D0-E1D42DE1D5F3}")
    buffers = []
    try:
        for pid, text in values.items():
            buf = ctypes.create_unicode_buffer(text)
            buffers.append(buf)         # alive until the call has copied the string
            var = PropVariant(vt=VT_LPWSTR, ptr=ctypes.cast(buf, ctypes.c_void_p).value)
            set_value(store, ctypes.byref(PropKey(fmtid, pid)), ctypes.byref(var))
        commit(store)
    finally:
        release(store)


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
    _taskbar_icon(win)
    # after the window is up, so a slow network never delays the start
    QTimer.singleShot(1500, win.startup_update_check)
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
