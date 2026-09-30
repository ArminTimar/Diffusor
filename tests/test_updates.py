"""The update notice: version comparison, the GitHub reader and the banner.

No test touches the network. ``updates.urlopen`` is replaced by a fake, and the
window's settings are redirected to a temporary folder.
"""
import io
import json
import os
from datetime import datetime, timedelta, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError

import pytest

import diffusor
from diffusor import updates

ROOT = Path(__file__).resolve().parents[1]


# --- versions ------------------------------------------------------------------
def test_parse_version_reads_plain_versions_only():
    assert updates.parse_version("v0.2.0") == (0, 2, 0)
    assert updates.parse_version("1.10") == (1, 10)
    assert updates.parse_version(" V2.0.1 ") == (2, 0, 1)
    for bad in ("", "latest", "0.2.0rc1", "v1.x", "1..2", None):
        assert updates.parse_version(bad) is None, bad


def test_is_newer_compares_numbers_not_text():
    assert updates.is_newer("v0.2.0", "0.1.0")
    assert updates.is_newer("0.10.0", "0.9.0")          # text order would say the opposite
    assert updates.is_newer("1.0.1", "1.0")
    assert not updates.is_newer("1.0", "1.0.0")          # the same release
    assert not updates.is_newer("0.1.0", "0.1.0")
    assert not updates.is_newer("0.0.9", "0.1.0")
    assert not updates.is_newer("nightly", "0.1.0")      # unreadable: no notice
    assert not updates.is_newer("0.2.0", "dev")


def test_the_version_is_written_in_one_place():
    assert updates.parse_version(diffusor.__version__) is not None
    toml = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert 'dynamic = ["version"]' in toml
    assert 'attr = "diffusor.__version__"' in toml
    assert '\nversion = "' not in toml


# --- the release record ----------------------------------------------------------
def test_release_from_json_reads_tag_link_and_notes():
    info = updates.release_from_json({
        "tag_name": "v0.2.0", "body": "  Faster fits.\n",
        "html_url": f"https://github.com/{updates.REPO}/releases/tag/v0.2.0"})
    assert (info.version, info.tag) == ("0.2.0", "v0.2.0")
    assert info.url.endswith("/releases/tag/v0.2.0")
    assert info.notes == "Faster fits."


def test_a_link_to_another_site_is_replaced_by_the_projects_own_page():
    for evil in ("https://example.com/download", "http://github.com/x/y/releases/1",
                 f"https://github.com/{updates.REPO}-evil/releases/1", "", None):
        info = updates.release_from_json({"tag_name": "v0.2.0", "html_url": evil})
        assert info.url == f"https://github.com/{updates.REPO}/releases/tag/v0.2.0"


def test_a_tag_that_is_not_a_version_is_refused():
    with pytest.raises(updates.UpdateError):
        updates.release_from_json({"tag_name": "nightly"})
    with pytest.raises(updates.UpdateError):
        updates.release_from_json({})


# --- reading GitHub --------------------------------------------------------------
class _Response(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def _serve(monkeypatch, body=None, error=None):
    seen = {}

    def fake(req, timeout=None):
        seen["url"], seen["timeout"] = req.full_url, timeout
        seen["agent"] = req.get_header("User-agent")
        if error is not None:
            raise error
        return _Response(body if isinstance(body, bytes) else json.dumps(body).encode())

    monkeypatch.setattr(updates, "urlopen", fake)
    return seen


def test_fetch_latest_asks_github_for_this_project(monkeypatch):
    seen = _serve(monkeypatch, {"tag_name": "v9.0.0", "body": "x"})
    info = updates.fetch_latest("0.1.0")
    assert info.version == "9.0.0"
    assert seen["url"] == f"https://api.github.com/repos/{updates.REPO}/releases/latest"
    assert seen["agent"] == "Diffusor/0.1.0"            # GitHub refuses requests without one
    assert seen["timeout"] == updates.TIMEOUT_SECONDS


@pytest.mark.parametrize("error, fragment", [
    (HTTPError("u", 404, "Not Found", {}, None), "No release"),
    (HTTPError("u", 403, "Forbidden", {}, None), "limiting"),
    (HTTPError("u", 500, "Oops", {}, None), "error 500"),
    (URLError("no route"), "Could not reach"),
    (TimeoutError(), "Could not reach"),
])
def test_every_failure_becomes_a_sentence(monkeypatch, error, fragment):
    _serve(monkeypatch, error=error)
    with pytest.raises(updates.UpdateError, match=fragment):
        updates.fetch_latest("0.1.0")


@pytest.mark.parametrize("body", [
    pytest.param(b"<html>not json</html>", id="html"),
    pytest.param(b"[1, 2]", id="list"),
    pytest.param(b"\xff\xfe", id="not-utf8"),
    pytest.param(b"x" * (updates.MAX_BYTES + 5), id="oversize")])
def test_a_reply_that_is_not_a_release_is_refused(monkeypatch, body):
    _serve(monkeypatch, body=body)
    with pytest.raises(updates.UpdateError):
        updates.fetch_latest("0.1.0")


# --- how often ------------------------------------------------------------------
def test_the_startup_check_runs_at_most_about_once_a_day():
    now = datetime(2026, 10, 1, 12, 0, tzinfo=timezone.utc)
    assert updates.due("", now)
    assert updates.due("garbage", now)
    assert not updates.due((now - timedelta(hours=2)).isoformat(), now)
    assert updates.due((now - timedelta(hours=21)).isoformat(), now)
    assert updates.due((now + timedelta(days=3)).isoformat(), now)   # clock was set back
    assert not updates.due((now - timedelta(hours=2)).replace(tzinfo=None).isoformat(), now)


# --- the worker -----------------------------------------------------------------
def test_the_worker_reports_a_release_or_a_reason(monkeypatch):
    from diffusor.gui.workers import UpdateWorker
    got, failed = [], []
    w = UpdateWorker("0.1.0")
    w.finished.connect(got.append)
    w.failed.connect(failed.append)

    info = updates.ReleaseInfo("0.2.0", "v0.2.0", "https://github.com/x", "")
    monkeypatch.setattr(updates, "fetch_latest", lambda current: info)
    w.run()
    assert got == [info] and not failed

    def refuse(current):
        raise updates.UpdateError("Could not reach GitHub.")
    monkeypatch.setattr(updates, "fetch_latest", refuse)
    w.run()
    assert failed == ["Could not reach GitHub."]

    def crash(current):
        raise RuntimeError("bug")
    monkeypatch.setattr(updates, "fetch_latest", crash)
    w.run()
    assert len(failed) == 2 and "unexpectedly" in failed[1]


# --- the window -----------------------------------------------------------------
@pytest.fixture()
def window(tmp_path, monkeypatch):
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    import matplotlib
    matplotlib.use("QtAgg")
    from PySide6.QtCore import QSettings
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])
    from diffusor.gui import main_window as mw
    ini = str(tmp_path / "settings.ini")        # a file of its own: the real settings stay untouched
    monkeypatch.setattr(mw.MainWindow, "_settings",
                        staticmethod(lambda: QSettings(ini, QSettings.IniFormat)))
    w = mw.MainWindow()
    w.resize(1440, 900)
    w.show()
    app.processEvents()
    yield w
    w.close()


def _release(version):
    return updates.ReleaseInfo(version, "v" + version,
                               f"https://github.com/{updates.REPO}/releases/tag/v{version}",
                               "Notes for the release.")


def test_the_banner_is_hidden_until_a_newer_release_is_found(window):
    assert not window.update_banner.isVisible()
    window._update_found(_release(diffusor.__version__))
    assert not window.update_banner.isVisible()
    window._update_found(_release("99.0.0"))
    assert window.update_banner.isVisible()
    assert "99.0.0" in window.lbl_update.text()
    assert diffusor.__version__ in window.lbl_update.text()
    assert "Notes for the release." in window.lbl_update.toolTip()


def test_later_hides_the_banner_and_download_opens_the_release_page(window, monkeypatch):
    from diffusor.gui import main_window as mw
    opened = []
    monkeypatch.setattr(mw.QDesktopServices, "openUrl", lambda url: opened.append(url.toString()))
    window._update_found(_release("99.0.0"))
    window.btn_update_later.click()
    assert not window.update_banner.isVisible() and not opened
    window._update_found(_release("99.0.0"))
    window.btn_update_get.click()
    assert opened == [f"https://github.com/{updates.REPO}/releases/tag/v99.0.0"]
    assert not window.update_banner.isVisible()


def test_a_manual_check_always_answers_and_a_startup_check_stays_quiet(window, monkeypatch):
    from diffusor.gui import main_window as mw
    shown = []
    monkeypatch.setattr(mw.QMessageBox, "information",
                        lambda parent, title, text: shown.append(text))
    window._update_manual = False
    window._update_failed("Could not reach GitHub.")
    window._update_found(_release(diffusor.__version__))
    assert shown == []                                   # a silent check never interrupts

    window._update_manual = True
    window._update_failed("Could not reach GitHub.")
    window._update_manual = True
    window._update_found(_release(diffusor.__version__))
    assert shown[0] == "Could not reach GitHub."
    assert "newest version" in shown[1]
    window._update_manual = True
    window._update_found(_release("99.0.0"))
    assert len(shown) == 2 and window.update_banner.isVisible()   # the banner is the answer


def test_a_successful_check_is_remembered_and_a_failed_one_is_not(window):
    assert window._settings().value("updates/last_check", "") == ""
    window._update_failed("Could not reach GitHub.")
    assert window._settings().value("updates/last_check", "") == ""
    window._update_found(_release(diffusor.__version__))
    assert not updates.due(str(window._settings().value("updates/last_check", "")))


def test_startup_check_respects_the_switch_and_the_daily_limit(window, monkeypatch):
    calls = []
    monkeypatch.setattr(window, "check_for_updates", lambda manual=False: calls.append(manual))
    window.startup_update_check()
    assert calls == [False]                              # on by default, never checked

    window.act_update_startup.setChecked(False)
    window.startup_update_check()
    assert calls == [False]                              # switched off
    window.act_update_startup.setChecked(True)

    window._settings().setValue("updates/last_check", datetime.now(timezone.utc).isoformat())
    window.startup_update_check()
    assert calls == [False]                              # checked a moment ago


def test_the_startup_switch_is_remembered_between_launches(window):
    window.act_update_startup.setChecked(False)
    from diffusor.gui import main_window as mw
    second = mw.MainWindow()
    try:
        assert not second.act_update_startup.isChecked()
    finally:
        second.close()


def test_the_help_menu_offers_the_check(window):
    help_menu = [a.menu() for a in window.menuBar().actions() if a.text() == "&Help"][0]
    texts = [a.text() for a in help_menu.actions()]
    assert "Check for updates..." in texts
    assert "Check for updates at startup" in texts


def test_a_check_started_from_the_window_goes_through_the_worker(window, monkeypatch):
    from PySide6.QtWidgets import QApplication
    import time
    monkeypatch.setattr(updates, "fetch_latest", lambda current: _release("99.0.0"))
    window.check_for_updates(manual=True)
    for _ in range(100):
        QApplication.processEvents()
        if window.update_banner.isVisible():
            break
        time.sleep(0.05)
    assert window.update_banner.isVisible()
    window._update_job[0].wait(2000)
