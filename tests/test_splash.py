"""The start-up splash: it opens first, its bar only moves forward, it reports a failed import."""
import os
import time

import pytest


@pytest.fixture(scope="module")
def app():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _wait(app, condition, seconds=5.0):
    end = time.time() + seconds
    while time.time() < end and not condition():
        app.processEvents()
        time.sleep(0.005)
    return condition()


def test_the_preloader_runs_every_stage_in_order_on_another_thread():
    import threading
    from diffusor.gui.splash import Preloader
    seen = []
    stages = [("a", 0.4, 0.1, lambda: seen.append(("a", threading.current_thread().name))),
              ("b", 0.9, 0.1, lambda: seen.append(("b", threading.current_thread().name)))]
    loader = Preloader(stages)
    loader.start()
    loader._thread.join(5)
    assert loader.done and loader.error is None and loader.stage == 2
    assert [s for s, _ in seen] == ["a", "b"]
    assert all(name == "diffusor-preload" for _, name in seen)


def test_a_failed_import_is_reported_and_does_not_hang_the_start():
    from diffusor.gui.splash import Preloader

    def broken():
        raise ImportError("no module named nothing")
    loader = Preloader([("a", 0.5, 0.1, broken), ("b", 1.0, 0.1, lambda: None)])
    loader.start()
    loader._thread.join(5)
    assert loader.done and "no module named nothing" in loader.error
    assert loader.stage == 0, "it stopped at the stage that failed"


def test_the_bar_only_moves_forward_and_stays_short_until_the_work_is_done(app):
    from diffusor.gui.splash import Preloader, StartupSplash
    release = []
    loader = Preloader([("slow", 0.5, 1.0, lambda: _wait(app, lambda: bool(release), 10))])
    s = StartupSplash(loader, "9.9.9")
    s.begin()
    loader.start()
    last = 0.0
    for _ in range(60):
        app.processEvents()
        time.sleep(0.01)
        assert s.value >= last - 1e-12
        last = s.value
    assert 0.0 < s.value < 0.5, "it creeps but does not claim the stage is over"
    assert s.message == "slow"
    release.append(1)
    loader._thread.join(5)
    s.finish()


def test_complete_runs_the_bar_to_the_end_and_closes_the_window(app):
    from diffusor.gui.splash import StartupSplash
    s = StartupSplash(None, "9.9.9")
    s.begin()
    assert s.isVisible()
    s.complete()
    assert _wait(app, lambda: not s.isVisible(), 3.0)
    assert s.value >= 0.995


def test_the_splash_paints_without_the_application_style_sheet(app):
    from diffusor.gui.splash import HEIGHT, WIDTH, StartupSplash
    s = StartupSplash(None, "9.9.9")
    image = s.grab().toImage()
    assert image.width() >= WIDTH and image.height() >= HEIGHT
    s.close()
