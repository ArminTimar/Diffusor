"""The start-up window: the icon, the name and a bar, shown while Diffusor loads.

Starting takes several seconds, almost all of it spent importing the numerical and
plotting libraries. This window opens as soon as Qt is up, before any of that, so a
click on the shortcut is answered at once. The heavy imports then run on a second
thread while the bar moves, and the main window is built on the interface thread
when they are done.

The splash paints itself and uses nothing from the application style sheet, which
is not set yet when it appears.
"""
from __future__ import annotations

import threading
import traceback
from typing import Callable, List, Optional, Sequence, Tuple

from PySide6.QtCore import QRectF, Qt, QTimer
from PySide6.QtGui import QColor, QFont, QGuiApplication, QPainter, QPainterPath, QPen, QPixmap
from PySide6.QtWidgets import QWidget

from . import theme

WIDTH, HEIGHT = 440, 292
ICON_SIZE = 92
BAR_WIDTH, BAR_HEIGHT = 264, 6

# One entry per stage: what the user is told, how far the bar goes when it ends, and
# a guess at how long it takes, which only sets how fast the bar creeps meanwhile.
Stage = Tuple[str, float, float, Callable[[], None]]


def _import_numerics() -> None:
    import numpy  # noqa: F401
    import scipy.linalg  # noqa: F401
    import scipy.optimize  # noqa: F401


def _import_plotting() -> None:
    import matplotlib
    matplotlib.use("QtAgg")
    import matplotlib.backends.backend_qtagg  # noqa: F401
    import matplotlib.figure  # noqa: F401


def _import_science() -> None:
    from .. import coefficients, fitting, minerals  # noqa: F401


def _import_interface() -> None:
    from . import main_window  # noqa: F401


# (message, bar value when the stage is done, seconds the stage usually takes)
STAGES: List[Tuple[str, float, float, Callable[[], None]]] = [
    ("Loading numerical libraries", 0.50, 2.2, _import_numerics),
    ("Loading plotting", 0.70, 0.8, _import_plotting),
    ("Loading diffusion laws", 0.78, 0.3, _import_science),
    ("Loading the interface", 0.90, 0.7, _import_interface),
]


class Preloader:
    """Runs the stages on a thread and reports where it is, without Qt signals.

    The splash reads ``stage``, ``done`` and ``error`` from its timer. They are plain
    attributes written by one thread and read by another, which is safe here because
    each is only ever assigned a whole new value.
    """

    def __init__(self, stages: Sequence[Stage] = ()):
        self.stages = list(stages or STAGES)
        self.stage = 0                       # index of the stage running now
        self.done = False
        self.error: Optional[str] = None
        self._thread = threading.Thread(target=self._run, name="diffusor-preload", daemon=True)

    def start(self) -> None:
        self._thread.start()

    def _run(self) -> None:
        try:
            for i, (_, _, _, work) in enumerate(self.stages):
                self.stage = i
                work()
            self.stage = len(self.stages)
        except Exception:
            self.error = traceback.format_exc()
        self.done = True


class StartupSplash(QWidget):
    """A small frameless window with Diffusor's icon, its name, a status line and a bar."""

    def __init__(self, preloader: Optional[Preloader] = None, version: str = ""):
        super().__init__(None, Qt.SplashScreen | Qt.FramelessWindowHint | Qt.WindowStaysOnTopHint)
        self.setAttribute(Qt.WA_TranslucentBackground, True)
        self.setFixedSize(WIDTH, HEIGHT)
        self.preloader = preloader
        self.version = version
        self.message = "Starting"
        self.value = 0.0                     # what the bar shows
        self.target = 0.04                   # where it is heading
        self._icon: Optional[QPixmap] = None
        self._completing = False
        self._timer = QTimer(self)
        self._timer.setInterval(16)
        self._timer.timeout.connect(self._tick)
        screen = QGuiApplication.primaryScreen()
        if screen is not None:
            area = screen.availableGeometry()
            self.move(area.center().x() - WIDTH // 2, area.center().y() - HEIGHT // 2)

    # -- driving the bar ------------------------------------------------------
    def begin(self) -> None:
        self.show()
        self._timer.start()

    def set_progress(self, target: float, message: Optional[str] = None) -> None:
        """Send the bar towards ``target`` (0 to 1) and optionally change the line under it."""
        self.target = max(self.target, min(1.0, float(target)))
        if message is not None:
            self.message = message
        self.update()

    def _stage_goal(self) -> Tuple[float, str]:
        """Where the bar may go while the loading thread is still working."""
        pre = self.preloader
        if pre is None:
            return self.target, self.message
        stages = pre.stages
        i = min(pre.stage, len(stages) - 1)
        start = stages[i - 1][1] if i > 0 else 0.0
        message, end, _, _ = stages[i]
        if pre.done:
            return stages[-1][1], self.message
        # stay a little short of the end of the stage until the stage reports it is over
        return start + 0.92 * (end - start), message

    def _tick(self) -> None:
        goal, message = self._stage_goal()
        if self.preloader is not None and not self.preloader.done:
            self.target = max(self.target, goal)
            self.message = message
        speed = 0.3 if self._completing else 0.045      # share of the remaining distance per frame
        gap = self.target - self.value
        if gap > 0:
            self.value += max(gap * speed, 0.0006)
            self.value = min(self.value, self.target)
        self.update()
        if self._completing and self.value >= 0.995:
            self._completing = False
            QTimer.singleShot(140, self.finish)         # a moment with the bar full

    def complete(self, message: str = "Ready") -> None:
        """Run the bar to the end quickly, then close."""
        self.set_progress(1.0, message)
        self._completing = True
        QTimer.singleShot(1500, self.finish)            # in case the timer never gets there

    def finish(self) -> None:
        self._timer.stop()
        self.close()

    # -- painting -------------------------------------------------------------
    def _icon_pixmap(self) -> QPixmap:
        if self._icon is None:
            ratio = self.devicePixelRatioF() or 1.0
            px = QPixmap(str(theme.ICONS / "diffusor.png"))
            size = int(ICON_SIZE * ratio)
            self._icon = px.scaled(size, size, Qt.KeepAspectRatio, Qt.SmoothTransformation)
            self._icon.setDevicePixelRatio(ratio)
        return self._icon

    def paintEvent(self, _event) -> None:
        p = QPainter(self)
        p.setRenderHints(QPainter.Antialiasing | QPainter.TextAntialiasing
                         | QPainter.SmoothPixmapTransform)
        card = QRectF(0.5, 0.5, WIDTH - 1, HEIGHT - 1)
        path = QPainterPath()
        path.addRoundedRect(card, 16, 16)
        p.fillPath(path, QColor(theme.BG))
        p.setPen(QPen(QColor(theme.BORDER_STRONG), 1))
        p.drawPath(path)

        icon = self._icon_pixmap()
        p.drawPixmap(int((WIDTH - ICON_SIZE) / 2), 34, icon)

        title = QFont("Segoe UI", 22)
        title.setWeight(QFont.DemiBold)
        p.setFont(title)
        p.setPen(QColor(theme.TEXT))
        p.drawText(QRectF(0, 138, WIDTH, 40), Qt.AlignHCenter | Qt.AlignVCenter, "Diffusor")

        small = QFont("Segoe UI", 10)
        p.setFont(small)
        p.setPen(QColor(theme.TEXT_MUTED))
        p.drawText(QRectF(0, 182, WIDTH, 22), Qt.AlignHCenter | Qt.AlignVCenter,
                   self.message + "...")

        left = (WIDTH - BAR_WIDTH) / 2
        track = QRectF(left, 220, BAR_WIDTH, BAR_HEIGHT)
        p.setPen(Qt.NoPen)
        p.setBrush(QColor(theme.BORDER))
        p.drawRoundedRect(track, BAR_HEIGHT / 2, BAR_HEIGHT / 2)
        fill = QRectF(left, 220, max(BAR_HEIGHT, BAR_WIDTH * self.value), BAR_HEIGHT)
        p.setBrush(QColor(theme.ACCENT))
        p.drawRoundedRect(fill, BAR_HEIGHT / 2, BAR_HEIGHT / 2)

        if self.version:
            tiny = QFont("Segoe UI", 8)
            p.setFont(tiny)
            p.setPen(QColor(theme.TEXT_FAINT))
            p.drawText(QRectF(0, 250, WIDTH, 22), Qt.AlignHCenter | Qt.AlignVCenter,
                       f"Version {self.version}")
        p.end()
