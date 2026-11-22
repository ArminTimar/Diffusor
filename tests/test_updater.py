"""The one-click update: release assets, how this copy installs, the download and hand-over.

No test touches the network or runs an installer: ``updater.urlopen`` and the
process launcher are replaced by fakes, and every file lives in a temporary folder.
"""
import hashlib
import io
import os
import shutil
import subprocess
from pathlib import Path

import pytest

import diffusor
from diffusor import updater, updates

ROOT = Path(__file__).resolve().parents[1]
PREFIX = f"https://github.com/{updates.REPO}/releases/download/v9.0.0/"
NAMES = ("Diffusor-9.0.0-Linux-x86_64.AppImage", "Diffusor-9.0.0-macOS-AppleSilicon.dmg",
         "Diffusor-9.0.0-macOS-Intel.dmg", "Diffusor-9.0.0-Windows-Setup.exe",
         "diffusor_9.0.0_amd64.deb")


def _asset_json(name, data=b"x", digest=True):
    item = {"name": name, "browser_download_url": PREFIX + name, "size": len(data)}
    if digest:
        item["digest"] = "sha256:" + hashlib.sha256(data).hexdigest()
    return item


def _release(assets=NAMES, digest=True):
    return updates.release_from_json({
        "tag_name": "v9.0.0", "html_url": f"https://github.com/{updates.REPO}/releases/tag/v9.0.0",
        "assets": [_asset_json(n, digest=digest) for n in assets]})


# --- the release's files ------------------------------------------------------------
def test_the_installers_are_read_from_the_release_record():
    info = _release()
    assert [a.name for a in info.assets] == list(NAMES)
    a = info.assets[3]
    assert a.url == PREFIX + a.name and a.size == 1
    assert a.sha256 == hashlib.sha256(b"x").hexdigest()


def test_doubtful_files_are_left_out():
    good = _asset_json("Diffusor-9.0.0-Windows-Setup.exe")
    bad = [
        dict(good, browser_download_url="https://example.com/Diffusor-9.0.0-Windows-Setup.exe"),
        dict(good, browser_download_url=f"https://github.com/{updates.REPO}-evil/releases/"
                                        "download/v9.0.0/Diffusor-9.0.0-Windows-Setup.exe"),
        dict(good, browser_download_url=PREFIX + "other.exe"),     # the name and link differ
        dict(good, name="../Diffusor-9.0.0-Windows-Setup.exe"),
        dict(good, name=".hidden"),
        dict(good, size=0), dict(good, size="12"), dict(good, size=True),
        "not a record", None,
    ]
    assert updates.assets_from_json(bad + [good]) == updates.assets_from_json([good])
    assert len(updates.assets_from_json([good])) == 1
    assert updates.assets_from_json(None) == () and updates.assets_from_json("x") == ()
    no_digest = updates.assets_from_json([dict(good, digest="md5:abc")])
    assert no_digest[0].sha256 == ""


def test_a_release_without_assets_still_reads():
    info = updates.release_from_json({"tag_name": "v9.0.0"})
    assert info.assets == ()
    assert updates.ReleaseInfo("1.0", "v1.0", "u", "").assets == ()


# --- how this copy is installed -------------------------------------------------------
def test_run_from_source_has_no_plan():
    assert updater.detect() is None             # the tests run from source
    assert updater.detect(frozen=False, system="Windows") is None


def test_windows_needs_the_inno_setup_uninstaller(tmp_path):
    exe = tmp_path / "Diffusor.exe"
    exe.write_bytes(b"")
    assert updater.detect(True, str(exe), "Windows") is None       # a copied folder
    (tmp_path / "unins000.exe").write_bytes(b"")
    plan = updater.detect(True, str(exe), "Windows")
    assert plan == updater.Plan("windows", tmp_path, True)


def test_macos_needs_a_writable_app_that_is_not_translocated(tmp_path):
    exe = tmp_path / "Applications" / "Diffusor.app" / "Contents" / "MacOS" / "Diffusor"
    exe.parent.mkdir(parents=True)
    plan = updater.detect(True, str(exe), "Darwin")
    assert plan == updater.Plan("macos", tmp_path / "Applications" / "Diffusor.app", True)
    moved = tmp_path / "AppTranslocation" / "X" / "d" / "Diffusor.app" / "Contents" / "MacOS"
    moved.mkdir(parents=True)
    assert updater.detect(True, str(moved / "Diffusor"), "Darwin") is None
    assert updater.detect(True, str(tmp_path / "Diffusor"), "Darwin") is None


def test_linux_appimage_and_deb(tmp_path):
    image = tmp_path / "Diffusor-0.8.0-Linux-x86_64.AppImage"
    image.write_bytes(b"")
    plan = updater.detect(True, "/tmp/.mount_x/usr/lib/diffusor/Diffusor", "Linux",
                          {"APPIMAGE": str(image)})
    assert plan == updater.Plan("appimage", image, True)
    assert updater.detect(True, "/x/Diffusor", "Linux",
                          {"APPIMAGE": str(tmp_path / "gone.AppImage")}) is None
    deb = updater.detect(True, "/opt/diffusor/Diffusor", "Linux", {})
    assert deb is not None and deb.kind == "deb" and not deb.restarts
    assert updater.detect(True, "/home/u/Diffusor/Diffusor", "Linux", {}) is None


def test_each_installation_picks_its_own_installer():
    info = _release()
    pick = lambda kind, machine=None: updater.pick_asset(info, updater.Plan(kind, Path("."), True),
                                                         machine)
    assert pick("windows").name == "Diffusor-9.0.0-Windows-Setup.exe"
    assert pick("macos", "arm64").name == "Diffusor-9.0.0-macOS-AppleSilicon.dmg"
    assert pick("macos", "x86_64").name == "Diffusor-9.0.0-macOS-Intel.dmg"
    assert pick("appimage").name == "Diffusor-9.0.0-Linux-x86_64.AppImage"
    assert pick("deb").name == "diffusor_9.0.0_amd64.deb"


def test_no_installer_without_a_digest_or_before_the_files_are_attached():
    plan = updater.Plan("windows", Path("."), True)
    assert updater.pick_asset(_release(digest=False), plan) is None
    assert updater.pick_asset(_release(assets=()), plan) is None
    # a file of another version attached by mistake is not used
    other = _release(assets=("Diffusor-8.0.0-Windows-Setup.exe",))
    assert updater.pick_asset(other, plan) is None


def test_the_file_names_match_what_the_build_produces():
    build = (ROOT / "packaging" / "build.py").read_text(encoding="utf-8")
    for fragment in ("-Windows-Setup", "-macOS-{kind}.dmg", "-Linux-x86_64.AppImage",
                     "_amd64.deb", "AppleSilicon", "Intel"):
        assert fragment in build, fragment


# --- the download --------------------------------------------------------------------
class _Body(io.BytesIO):
    def __init__(self, data, url="https://objects.githubusercontent.com/x"):
        super().__init__(data)
        self._url = url

    def geturl(self):
        return self._url

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def _asset(data, name="Diffusor-9.0.0-Windows-Setup.exe", size=None, sha=None):
    return updates.Asset(name, PREFIX + name, len(data) if size is None else size,
                         hashlib.sha256(data).hexdigest() if sha is None else sha)


def _serve(monkeypatch, data, url="https://objects.githubusercontent.com/x"):
    seen = {}

    def fake(req, timeout=None):
        seen["url"], seen["agent"] = req.full_url, req.get_header("User-agent")
        return _Body(data, url)

    monkeypatch.setattr(updater, "urlopen", fake)
    monkeypatch.setattr(updater, "CHUNK", 1000)
    return seen


def test_a_matching_download_is_kept(monkeypatch, tmp_path):
    data = os.urandom(5500)
    seen = _serve(monkeypatch, data)
    steps = []
    path = updater.download(_asset(data), tmp_path, "Diffusor/1.0",
                            progress=lambda got, total: steps.append((got, total)))
    assert path == tmp_path / "Diffusor-9.0.0-Windows-Setup.exe"
    assert path.read_bytes() == data
    assert seen["url"].startswith(PREFIX) and seen["agent"] == "Diffusor/1.0"
    assert steps[-1] == (5500, 5500) and len(steps) == 6
    assert sorted(p.name for p in tmp_path.iterdir()) == [path.name]


@pytest.mark.parametrize("asset_for, fragment", [
    (lambda d: _asset(d, sha="0" * 64), "does not match"),
    (lambda d: _asset(d, size=len(d) - 1), "larger than"),
    (lambda d: _asset(d, size=len(d) + 1), "before the end"),
])
def test_a_wrong_download_is_deleted(monkeypatch, tmp_path, asset_for, fragment):
    data = os.urandom(3000)
    _serve(monkeypatch, data)
    with pytest.raises(updates.UpdateError, match=fragment):
        updater.download(asset_for(data), tmp_path)
    assert list(tmp_path.iterdir()) == []


def test_a_cancelled_download_leaves_nothing(monkeypatch, tmp_path):
    data = os.urandom(3000)
    _serve(monkeypatch, data)
    with pytest.raises(updater.Cancelled):
        updater.download(_asset(data), tmp_path, cancelled=lambda: True)
    assert list(tmp_path.iterdir()) == []


def test_only_this_projects_files_are_downloaded(monkeypatch, tmp_path):
    data = b"abc"
    _serve(monkeypatch, data)
    elsewhere = updates.Asset("x.exe", "https://example.com/x.exe", 3,
                              hashlib.sha256(data).hexdigest())
    with pytest.raises(updates.UpdateError):
        updater.download(elsewhere, tmp_path)
    with pytest.raises(updates.UpdateError):                    # no digest to check against
        updater.download(_asset(data, sha=""), tmp_path)
    _serve(monkeypatch, data, url="http://objects.githubusercontent.com/x")
    with pytest.raises(updates.UpdateError, match="secure"):
        updater.download(_asset(data), tmp_path)
    assert list(tmp_path.iterdir()) == []


def test_a_network_failure_becomes_a_sentence(monkeypatch, tmp_path):
    from urllib.error import HTTPError, URLError

    for error in (URLError("down"), TimeoutError(), HTTPError("u", 503, "x", None, None)):
        def fake(req, timeout=None, error=error):
            raise error
        monkeypatch.setattr(updater, "urlopen", fake)
        with pytest.raises(updates.UpdateError) as exc:
            updater.download(_asset(b"abc"), tmp_path)
        assert str(exc.value).endswith(".")


def test_the_download_folder(tmp_path):
    image = tmp_path / "Diffusor.AppImage"
    assert updater.download_folder(updater.Plan("appimage", image, True)) == tmp_path
    win = updater.download_folder(updater.Plan("windows", tmp_path, True))
    assert win.name == updater.TEMP_FOLDER


# --- handing over to the installer ----------------------------------------------------
class _Launcher:
    def __init__(self):
        self.calls = []

    def __call__(self, args, **kwargs):
        self.calls.append((args, kwargs))


def test_windows_runs_the_setup_silently_and_asks_it_to_restart(tmp_path):
    setup = tmp_path / "Diffusor-9.0.0-Windows-Setup.exe"
    run = _Launcher()
    updater.start_install(updater.Plan("windows", tmp_path, True), setup, popen=run)
    (args, kwargs), = run.calls
    assert args[0] == str(setup)
    for flag in ("/VERYSILENT", "/SUPPRESSMSGBOXES", "/NORESTART", "/relaunch=1"):
        assert flag in args
    assert kwargs["stdin"] == subprocess.DEVNULL


def test_the_windows_installer_restarts_diffusor_only_for_an_in_app_update():
    iss = (ROOT / "packaging" / "windows" / "diffusor.iss").read_text(encoding="utf-8")
    assert "Flags: nowait postinstall skipifsilent" in iss       # a normal install asks
    assert "Check: RelaunchAfterUpdate" in iss
    assert "WizardSilent and (ExpandConstant('{param:relaunch|0}') = '1')" in iss
    assert "runasoriginaluser" in iss
    assert "CloseApplications=yes" in iss and "PrivilegesRequired=lowest" in iss


def test_macos_writes_the_swap_script_and_waits_for_this_process(tmp_path):
    dmg = tmp_path / "Diffusor-9.0.0-macOS-Intel.dmg"
    app = tmp_path / "Applications" / "Diffusor.app"
    run = _Launcher()
    updater.start_install(updater.Plan("macos", app, True), dmg, pid=4321, popen=run)
    (args, kwargs), = run.calls
    script = tmp_path / "install-update.sh"
    assert args == ["/bin/sh", str(script), "4321", str(dmg), str(app)]
    assert kwargs["start_new_session"]
    text = script.read_bytes()
    assert b"\r" not in text and text.decode() == updater.MACOS_SCRIPT


@pytest.mark.skipif(shutil.which("sh") is None, reason="needs a POSIX shell")
def test_the_macos_script_is_valid_shell(tmp_path):
    script = tmp_path / "s.sh"
    script.write_bytes(updater.MACOS_SCRIPT.encode())
    assert subprocess.run([shutil.which("sh"), "-n", str(script)]).returncode == 0


def test_the_macos_script_removes_the_old_app_only_after_a_full_copy():
    s = updater.MACOS_SCRIPT
    assert s.index('ditto "$src" "$app.new"') < s.index('rm -rf "$app" && mv')
    assert 'kill -0 "$pid"' in s and s.rstrip().endswith('open "$app"')


def test_appimage_starts_the_new_file_then_removes_the_old_one(tmp_path):
    old = tmp_path / "Diffusor-0.8.0-Linux-x86_64.AppImage"
    new = tmp_path / "Diffusor-9.0.0-Linux-x86_64.AppImage"
    old.write_bytes(b"old"); new.write_bytes(b"new")
    run = _Launcher()
    updater.start_install(updater.Plan("appimage", old, True), new, popen=run)
    assert run.calls[0][0] == [str(new)]
    assert new.exists() and not old.exists()


def test_appimage_keeps_the_old_file_if_the_new_one_does_not_start(tmp_path):
    old = tmp_path / "old.AppImage"
    new = tmp_path / "new.AppImage"
    old.write_bytes(b"old"); new.write_bytes(b"new")

    def refuse(args, **kwargs):
        raise OSError("not executable")
    with pytest.raises(OSError):
        updater.start_install(updater.Plan("appimage", old, True), new, popen=refuse)
    assert old.exists()


def test_deb_opens_the_package_in_the_software_installer(tmp_path):
    deb = tmp_path / "diffusor_9.0.0_amd64.deb"
    run = _Launcher()
    updater.start_install(updater.Plan("deb", Path("/opt/diffusor"), False), deb, popen=run)
    assert run.calls[0][0] == ["xdg-open", str(deb)]


# --- the banner --------------------------------------------------------------------
@pytest.fixture()
def window(tmp_path, monkeypatch):
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    import matplotlib
    matplotlib.use("QtAgg")
    from PySide6.QtCore import QSettings
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])
    from diffusor.gui import main_window as mw
    ini = str(tmp_path / "settings.ini")
    monkeypatch.setattr(mw.MainWindow, "_settings",
                        staticmethod(lambda: QSettings(ini, QSettings.IniFormat)))
    w = mw.MainWindow()
    w.resize(1440, 900)
    w.show()
    app.processEvents()
    yield w
    w.close()


def _wait(predicate, seconds=10.0):
    import time
    from PySide6.QtWidgets import QApplication
    end = time.monotonic() + seconds
    while time.monotonic() < end:
        QApplication.processEvents()
        if predicate():
            return True
        time.sleep(0.01)
    return False


def _installable(monkeypatch, tmp_path, kind="windows"):
    plan = updater.Plan(kind, tmp_path, kind != "deb")
    monkeypatch.setattr(updater, "detect", lambda: plan)
    monkeypatch.setattr(updater, "download_folder", lambda p: tmp_path / "dl")
    return plan


def test_a_copy_that_can_update_offers_install(window, monkeypatch, tmp_path):
    _installable(monkeypatch, tmp_path)
    window._update_found(_release())
    assert window.update_banner.isVisible()
    assert window.btn_update_get.text() == "Install update"
    assert not window.btn_update_notes.isHidden()


def test_without_an_installer_the_banner_opens_the_page(window, monkeypatch, tmp_path):
    from diffusor.gui import main_window as mw
    opened = []
    monkeypatch.setattr(mw.QDesktopServices, "openUrl", lambda url: opened.append(url.toString()))
    _installable(monkeypatch, tmp_path)
    window._update_found(_release(assets=()))                  # not attached yet
    assert window.btn_update_get.text() == "Download"
    assert window.btn_update_notes.isHidden()
    window.btn_update_get.click()
    assert opened and not window.update_banner.isVisible()


def test_install_downloads_then_restarts_on_request(window, monkeypatch, tmp_path):
    from diffusor.gui import main_window as mw
    plan = _installable(monkeypatch, tmp_path)
    saved = tmp_path / "dl" / NAMES[3]

    def fake_download(asset, folder, agent, progress=None, cancelled=None):
        assert asset.name == NAMES[3] and agent == f"Diffusor/{diffusor.__version__}"
        progress(1, 1)
        folder.mkdir(parents=True, exist_ok=True)
        saved.write_bytes(b"x")
        return saved
    monkeypatch.setattr(updater, "download", fake_download)
    started, closed = [], []
    monkeypatch.setattr(updater, "start_install", lambda p, path: started.append((p, path)))
    monkeypatch.setattr(mw.MainWindow, "close", lambda self: closed.append(True) or True)

    window._update_found(_release())
    window.btn_update_get.click()
    assert _wait(lambda: window._update_state == "ready")
    assert window.btn_update_get.text() == "Restart and install"
    assert not started                      # nothing happens until the user says so
    window.btn_update_get.click()
    assert started == [(plan, saved)] and closed == [True]


def test_cancel_stops_the_download(window, monkeypatch, tmp_path):
    import threading
    _installable(monkeypatch, tmp_path)
    go = threading.Event()

    def slow_download(asset, folder, agent, progress=None, cancelled=None):
        while not cancelled():
            go.wait(0.01)
        raise updater.Cancelled("cancelled")
    monkeypatch.setattr(updater, "download", slow_download)
    window._update_found(_release())
    window.btn_update_get.click()
    assert window._update_state == "downloading"
    assert window.btn_update_later.text() == "Cancel" and window.btn_update_get.isHidden()
    window.btn_update_later.click()
    assert _wait(lambda: window._update_state == "offer")
    assert window.update_banner.isVisible()
    assert window.btn_update_get.text() == "Install update"


def test_a_failed_download_says_why(window, monkeypatch, tmp_path):
    from diffusor.gui import main_window as mw
    _installable(monkeypatch, tmp_path)

    def broken(asset, folder, agent, progress=None, cancelled=None):
        raise updates.UpdateError("The download stopped.")
    monkeypatch.setattr(updater, "download", broken)
    said = []
    monkeypatch.setattr(mw.QMessageBox, "information", lambda *a: said.append(a[2]))
    window._update_found(_release())
    window.btn_update_get.click()
    assert _wait(lambda: bool(said))
    assert said == ["The download stopped."] and window._update_state == "offer"


def test_an_installer_that_does_not_start_keeps_diffusor_open(window, monkeypatch, tmp_path):
    from diffusor.gui import main_window as mw
    _installable(monkeypatch, tmp_path)
    window._update_found(_release())
    window._update_file = tmp_path / "setup.exe"
    window._update_file.write_bytes(b"x")

    def refuse(plan, path):
        raise OSError("blocked")
    monkeypatch.setattr(updater, "start_install", refuse)
    warned, closed = [], []
    monkeypatch.setattr(mw.QMessageBox, "warning", lambda *a: warned.append(a[2]))
    monkeypatch.setattr(mw.MainWindow, "close", lambda self: closed.append(True) or True)
    window._install_update()
    assert warned and "blocked" in warned[0] and "setup.exe" in warned[0]
    assert not closed


def test_deb_hands_over_and_stays_open(window, monkeypatch, tmp_path):
    from diffusor.gui import main_window as mw
    _installable(monkeypatch, tmp_path, kind="deb")
    saved = tmp_path / "dl" / NAMES[4]

    def fake_download(asset, folder, agent, progress=None, cancelled=None):
        folder.mkdir(parents=True, exist_ok=True)
        saved.write_bytes(b"x")
        return saved
    monkeypatch.setattr(updater, "download", fake_download)
    started, closed = [], []
    monkeypatch.setattr(updater, "start_install", lambda p, path: started.append(path))
    monkeypatch.setattr(mw.MainWindow, "close", lambda self: closed.append(True) or True)
    window._update_found(_release())
    window.btn_update_get.click()
    assert _wait(lambda: window._update_state == "handed")
    assert started == [saved] and not closed
    assert window.btn_update_later.text() == "Close"
