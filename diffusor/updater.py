"""Download the installer of a newer release and hand over to it.

``updates`` finds out whether a newer release exists; this module installs it, so
that an update is one click and a restart instead of a download and a setup run.
How it is installed depends on how this copy of Diffusor was installed, which
``detect`` works out:

windows   the Inno Setup installation (``unins000.exe`` next to ``Diffusor.exe``).
          The new setup runs silently over it, keeps the folder, the shortcuts and
          the settings, and starts Diffusor again: ``/relaunch=1`` is read by
          packaging/windows/diffusor.iss.
macos     a ``Diffusor.app`` in a folder the user can write to. A short shell script
          waits for this process to end, copies the new ``Diffusor.app`` out of the
          disk image in place of the old one and opens it.
appimage  the AppImage named by ``$APPIMAGE``. The new AppImage is saved next to it
          and started, and the old file is removed.
deb       the package installed in /opt/diffusor. Installing a package needs an
          administrator password, so the downloaded package is opened in the
          system's software installer and Diffusor stays open.

Run from source, from pip, or from a copy that is none of these, ``detect`` returns
None and the update notice opens the release page instead, as before.

A file is downloaded only from this project's own release assets, only when GitHub
gives its SHA-256 digest, and it is kept only if its size and digest both match.
"""
from __future__ import annotations

import hashlib
import os
import platform
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Mapping, Optional
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from .updates import DOWNLOAD_PREFIX, Asset, ReleaseInfo, UpdateError

TIMEOUT_SECONDS = 30.0          # for each read; a slow connection still finishes
CHUNK = 256 * 1024
TEMP_FOLDER = "Diffusor-update"

# What each kind of installation downloads, by the file names packaging/build.py gives.
SUFFIXES = {
    "windows": ("-Windows-Setup.exe",),
    "macos-arm": ("-macOS-AppleSilicon.dmg",),
    "macos-intel": ("-macOS-Intel.dmg",),
    "appimage": ("-Linux-x86_64.AppImage",),
    "deb": ("_amd64.deb",),
}

# Waits for Diffusor to quit, puts the new app in place of the old one and opens it.
# The old app is removed only after the new one has been copied in full, and the
# app is opened at the end whether or not the copy worked.
MACOS_SCRIPT = """#!/bin/sh
pid="$1"; dmg="$2"; app="$3"
while kill -0 "$pid" 2>/dev/null; do sleep 0.5; done
mnt="$(mktemp -d)"
if hdiutil attach -nobrowse -readonly -noautoopen -mountpoint "$mnt" "$dmg" >/dev/null; then
    src="$(ls -d "$mnt"/*.app 2>/dev/null | head -n 1)"
    if [ -n "$src" ]; then
        rm -rf "$app.new"
        if ditto "$src" "$app.new"; then
            rm -rf "$app" && mv "$app.new" "$app"
        fi
        rm -rf "$app.new"
    fi
    hdiutil detach "$mnt" -quiet
fi
rmdir "$mnt" 2>/dev/null
rm -f "$dmg"
open "$app"
"""


class Cancelled(UpdateError):
    """The user stopped the download."""


@dataclass(frozen=True)
class Plan:
    """How this copy of Diffusor can be updated."""
    kind: str           # "windows", "macos", "appimage" or "deb"
    target: Path        # the install folder, the .app bundle, the AppImage or /opt/diffusor
    restarts: bool      # True: Diffusor closes and the new version starts by itself


def detect(frozen: Optional[bool] = None, executable: Optional[str] = None,
           system: Optional[str] = None,
           environ: Optional[Mapping[str, str]] = None) -> Optional[Plan]:
    """The way this copy can be updated, or None when only the release page can help."""
    frozen = bool(getattr(sys, "frozen", False)) if frozen is None else frozen
    if not frozen:
        return None
    exe = Path(executable or sys.executable)
    system = system or platform.system()
    env = os.environ if environ is None else environ
    if system == "Windows":
        if (exe.parent / "unins000.exe").is_file():
            return Plan("windows", exe.parent, True)
        return None
    if system == "Darwin":
        app = next((p for p in exe.parents if p.suffix == ".app"), None)
        # a quarantined app runs from a read-only copy (App Translocation)
        if app is None or "AppTranslocation" in str(app) or not os.access(app.parent, os.W_OK):
            return None
        return Plan("macos", app, True)
    if system == "Linux":
        image = env.get("APPIMAGE")
        if image:
            p = Path(image)
            if p.is_file() and os.access(p.parent, os.W_OK):
                return Plan("appimage", p, True)
            return None
        if exe.parent == Path("/opt/diffusor"):
            return Plan("deb", exe.parent, False)
    return None


def pick_asset(release: ReleaseInfo, plan: Plan, machine: Optional[str] = None) -> Optional[Asset]:
    """The installer in ``release`` for this installation, or None if it has none yet.

    The installers are attached some minutes after a release is published, so a
    release found early can have none. An asset without a digest is not used.
    """
    key = plan.kind
    if key == "macos":
        m = (machine or platform.machine()).lower()
        # under Rosetta an Intel build reports x86_64 and keeps the Intel build
        key = "macos-arm" if m in ("arm64", "aarch64") else "macos-intel"
    for asset in release.assets:
        if asset.sha256 and asset.name.endswith(SUFFIXES[key]) and release.version in asset.name:
            return asset
    return None


def download_folder(plan: Plan) -> Path:
    """Where the installer is saved: next to the AppImage, so it can replace it in place."""
    if plan.kind == "appimage":
        return plan.target.parent
    if plan.kind == "deb":
        downloads = Path.home() / "Downloads"
        if downloads.is_dir():
            return downloads
    return Path(tempfile.gettempdir()) / TEMP_FOLDER


def remove_leftovers() -> None:
    """Delete installers left from an earlier update (a running setup cannot delete itself)."""
    shutil.rmtree(Path(tempfile.gettempdir()) / TEMP_FOLDER, ignore_errors=True)


def download(asset: Asset, folder: Path, user_agent: str = "Diffusor",
             progress: Optional[Callable[[int, int], None]] = None,
             cancelled: Optional[Callable[[], bool]] = None,
             timeout: float = TIMEOUT_SECONDS) -> Path:
    """Save ``asset`` in ``folder`` and return its path; raise ``UpdateError`` otherwise.

    The file is written under a temporary name and renamed only once its size and
    SHA-256 digest match the release record, so a broken download never sits where
    an installer is expected.
    """
    if not asset.url.startswith(DOWNLOAD_PREFIX) or not asset.sha256:
        raise UpdateError("The installer is not one of this project's release files.")
    folder = Path(folder)
    try:
        folder.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise UpdateError(f"Could not create the folder {folder} for the download.") from exc
    final, part = folder / asset.name, folder / (asset.name + ".part")
    req = Request(asset.url, headers={"Accept": "application/octet-stream",
                                      "User-Agent": user_agent})
    digest, got = hashlib.sha256(), 0
    try:
        with urlopen(req, timeout=timeout) as resp, open(part, "wb") as out:
            if not str(resp.geturl()).startswith("https://"):
                raise UpdateError("The download was redirected away from a secure address.")
            while True:
                if cancelled is not None and cancelled():
                    raise Cancelled("The download was cancelled.")
                block = resp.read(CHUNK)
                if not block:
                    break
                got += len(block)
                if got > asset.size:
                    raise UpdateError("The download was larger than the release says. "
                                      "It was deleted.")
                digest.update(block)
                out.write(block)
                if progress is not None:
                    progress(got, asset.size)
    except UpdateError:
        part.unlink(missing_ok=True)
        raise
    except HTTPError as exc:
        part.unlink(missing_ok=True)
        raise UpdateError(f"GitHub answered with error {exc.code}. Try again later.") from exc
    except (URLError, TimeoutError, OSError) as exc:
        part.unlink(missing_ok=True)
        raise UpdateError("The download stopped. Check the internet connection and "
                          "try again.") from exc
    if got != asset.size:
        part.unlink(missing_ok=True)
        raise UpdateError("The download stopped before the end. Try again.")
    if digest.hexdigest() != asset.sha256:
        part.unlink(missing_ok=True)
        raise UpdateError("The downloaded file does not match the release record, so it "
                          "was deleted. Try again, or download it from the release page.")
    os.replace(part, final)
    return final


def windows_command(installer: Path) -> list:
    """The silent setup run over the existing installation, starting Diffusor at the end."""
    return [str(installer), "/VERYSILENT", "/SUPPRESSMSGBOXES", "/NORESTART",
            "/FORCECLOSEAPPLICATIONS", "/NORESTARTAPPLICATIONS",
            f"/LOG={installer.with_suffix('.log')}", "/relaunch=1"]


def start_install(plan: Plan, installer: Path, pid: Optional[int] = None,
                  popen=subprocess.Popen) -> None:
    """Start the installation. For a plan that restarts, Diffusor must quit next.

    Raises ``OSError`` if the installer cannot be started; nothing has been
    changed then, and the running copy can carry on.
    """
    installer = Path(installer)
    quiet = dict(stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                 stderr=subprocess.DEVNULL, close_fds=True)
    if plan.kind == "windows":
        flags = (getattr(subprocess, "DETACHED_PROCESS", 0)
                 | getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0))
        popen(windows_command(installer), creationflags=flags, **quiet)
    elif plan.kind == "macos":
        script = installer.parent / "install-update.sh"
        script.write_text(MACOS_SCRIPT, encoding="utf-8", newline="\n")
        popen(["/bin/sh", str(script), str(pid or os.getpid()), str(installer),
               str(plan.target)], start_new_session=True, **quiet)
    elif plan.kind == "appimage":
        installer.chmod(0o755)
        popen([str(installer)], start_new_session=True, **quiet)
        # the new copy is running, so the old one can go; the running process keeps
        # its open file until it quits
        if installer.resolve() != plan.target.resolve():
            plan.target.unlink(missing_ok=True)
    elif plan.kind == "deb":
        popen(["xdg-open", str(installer)], start_new_session=True, **quiet)
    else:
        raise OSError(f"no way to install an update for '{plan.kind}'")
