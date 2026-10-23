"""Build the Diffusor installer for the system this runs on.

    python packaging/build.py                    # bundle, check and installer
    python packaging/build.py installer          # installer from an existing bundle

Steps, in this order:

bundle     PyInstaller makes dist/Diffusor (dist/Diffusor.app on macOS): Python,
           the libraries and Diffusor in one folder that needs nothing installed.
check      starts the bundled program with --self-check: the window is built off
           screen and a Monte Carlo fit runs on two worker processes.
installer  Windows: Inno Setup makes Diffusor-<version>-Windows-Setup.exe.
           macOS:   a disk image with the app and a link to Applications.
           Linux:   an AppImage and a .deb package.

The installers land in dist/. Needs: the build extra (``pip install -e .[build,images]``);
on Windows Inno Setup 6 (ISCC.exe on the PATH, in its default folder or named by the
ISCC variable); on Linux appimagetool (on the PATH or named by APPIMAGETOOL) and
dpkg-deb.

macOS signing and notarisation run when these environment variables are set:
DIFFUSOR_MACOS_SIGN_IDENTITY (the "Developer ID Application: ..." name of a
certificate in the keychain), and for notarisation APPLE_ID, APPLE_TEAM_ID and
APPLE_APP_PASSWORD (an app-specific password of that Apple ID).
"""
from __future__ import annotations

import argparse
import os
import platform
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DIST = ROOT / "dist"
WORK = ROOT / "build" / "pyinstaller"
STEPS = ("bundle", "check", "installer")


def version() -> str:
    text = (ROOT / "diffusor" / "__init__.py").read_text(encoding="utf8")
    return re.search(r'^__version__ = "([^"]+)"', text, re.M).group(1)


def run(cmd, shown=None, **kwargs) -> subprocess.CompletedProcess:
    """Run a command, stop the build if it fails. ``shown`` replaces it in the log."""
    print("+", " ".join(str(c) for c in (shown or cmd)), flush=True)
    return subprocess.run([str(c) for c in cmd], check=True, **kwargs)


def app_folder() -> Path:
    return DIST / ("Diffusor.app" if sys.platform == "darwin" else "Diffusor")


def executable() -> Path:
    if sys.platform == "win32":
        return DIST / "Diffusor" / "Diffusor.exe"
    if sys.platform == "darwin":
        return DIST / "Diffusor.app" / "Contents" / "MacOS" / "Diffusor"
    return DIST / "Diffusor" / "Diffusor"


# --- steps --------------------------------------------------------------------------
def bundle() -> None:
    run([sys.executable, "-m", "PyInstaller", "--noconfirm", "--clean",
         "--distpath", DIST, "--workpath", WORK, ROOT / "packaging" / "diffusor.spec"])


def check() -> None:
    WORK.mkdir(parents=True, exist_ok=True)
    report = WORK / "self-check.txt"
    report.unlink(missing_ok=True)
    env = dict(os.environ, QT_QPA_PLATFORM="offscreen")
    print(f"+ {executable()} --self-check", flush=True)
    proc = subprocess.run([str(executable()), "--self-check", str(report)], env=env,
                          timeout=900)
    print(report.read_text(encoding="utf8") if report.is_file() else "(no report written)")
    if proc.returncode != 0:
        raise SystemExit(f"the bundled program failed its self-check (exit code "
                         f"{proc.returncode})")


def installer() -> None:
    if not app_folder().exists():
        raise SystemExit(f"{app_folder()} does not exist; run the bundle step first")
    {"win32": installer_windows, "darwin": installer_macos}.get(
        sys.platform, installer_linux)()


# --- Windows ------------------------------------------------------------------------
def find_iscc() -> Path:
    candidates = [os.environ.get("ISCC"), shutil.which("iscc"),
                  Path(os.environ.get("ProgramFiles(x86)", r"C:\Program Files (x86)"))
                  / "Inno Setup 6" / "ISCC.exe",
                  Path(os.environ.get("LOCALAPPDATA", "")) / "Programs" / "Inno Setup 6"
                  / "ISCC.exe"]
    for c in candidates:
        if c and Path(c).is_file():
            return Path(c)
    raise SystemExit("Inno Setup 6 was not found. Install it from https://jrsoftware.org "
                     "or name ISCC.exe in the ISCC environment variable.")


def installer_windows() -> None:
    v = version()
    run([find_iscc(), "/Qp", f"/DAppVersion={v}", f"/DSourceDir={DIST / 'Diffusor'}",
         f"/DOutputDir={DIST}", f"/DOutputName=Diffusor-{v}-Windows-Setup",
         ROOT / "packaging" / "windows" / "diffusor.iss"])


# --- macOS --------------------------------------------------------------------------
def installer_macos() -> None:
    kind = "AppleSilicon" if platform.machine() == "arm64" else "Intel"
    dmg = DIST / f"Diffusor-{version()}-macOS-{kind}.dmg"
    stage = WORK / "dmg"
    shutil.rmtree(stage, ignore_errors=True)
    stage.mkdir(parents=True)
    shutil.copytree(app_folder(), stage / "Diffusor.app", symlinks=True)
    (stage / "Applications").symlink_to("/Applications")
    dmg.unlink(missing_ok=True)
    run(["hdiutil", "create", "-volname", "Diffusor", "-srcfolder", stage, "-ov",
         "-format", "UDZO", dmg])

    identity = os.environ.get("DIFFUSOR_MACOS_SIGN_IDENTITY")
    if not identity:
        print("DIFFUSOR_MACOS_SIGN_IDENTITY is not set: the app is signed ad hoc and "
              "macOS will ask the user to allow it under Privacy & Security.")
        return
    run(["codesign", "--sign", identity, "--timestamp", dmg])
    apple_id, team, password = (os.environ.get(k) for k in
                                ("APPLE_ID", "APPLE_TEAM_ID", "APPLE_APP_PASSWORD"))
    if not (apple_id and team and password):
        print("APPLE_ID, APPLE_TEAM_ID or APPLE_APP_PASSWORD is not set: the disk image "
              "is signed but not notarised.")
        return
    result = run(["xcrun", "notarytool", "submit", dmg, "--apple-id", apple_id,
                  "--team-id", team, "--password", password, "--wait"],
                 shown=["xcrun", "notarytool", "submit", dmg.name, "--wait"],
                 capture_output=True, text=True)
    print(result.stdout)
    if "status: Accepted" not in result.stdout:
        raise SystemExit("Apple did not accept the disk image for notarisation; "
                         "`xcrun notarytool log <id>` gives the reasons.")
    run(["xcrun", "stapler", "staple", dmg])


# --- Linux --------------------------------------------------------------------------
LINUX = ROOT / "packaging" / "linux"
ICON = ROOT / "diffusor" / "gui" / "icons" / "diffusor.png"
APPRUN = """#!/bin/sh
HERE="$(dirname "$(readlink -f "$0")")"
exec "$HERE/usr/lib/diffusor/Diffusor" "$@"
"""


def installer_linux() -> None:
    appimage()
    deb()


def appimage() -> None:
    tool = os.environ.get("APPIMAGETOOL") or shutil.which("appimagetool")
    if not tool:
        raise SystemExit("appimagetool was not found: put it on the PATH or name it in "
                         "the APPIMAGETOOL environment variable.")
    appdir = WORK / "Diffusor.AppDir"
    shutil.rmtree(appdir, ignore_errors=True)
    shutil.copytree(app_folder(), appdir / "usr" / "lib" / "diffusor", symlinks=True)
    shutil.copy(LINUX / "diffusor.desktop", appdir / "diffusor.desktop")
    shutil.copy(ICON, appdir / "diffusor.png")
    (appdir / "AppRun").write_text(APPRUN, encoding="utf8", newline="\n")
    (appdir / "AppRun").chmod(0o755)
    out = DIST / f"Diffusor-{version()}-Linux-x86_64.AppImage"
    out.unlink(missing_ok=True)
    # appimagetool is itself an AppImage; extracting it avoids needing FUSE here
    run([tool, appdir, out], env=dict(os.environ, ARCH="x86_64",
                                      APPIMAGE_EXTRACT_AND_RUN="1"))


def deb() -> None:
    v = version()
    root = WORK / "deb"
    shutil.rmtree(root, ignore_errors=True)
    shutil.copytree(app_folder(), root / "opt" / "diffusor", symlinks=True)
    (root / "usr" / "bin").mkdir(parents=True)
    (root / "usr" / "bin" / "diffusor").symlink_to("/opt/diffusor/Diffusor")
    apps = root / "usr" / "share" / "applications"
    icons = root / "usr" / "share" / "icons" / "hicolor" / "256x256" / "apps"
    apps.mkdir(parents=True)
    icons.mkdir(parents=True)
    shutil.copy(LINUX / "diffusor.desktop", apps / "diffusor.desktop")
    shutil.copy(ICON, icons / "diffusor.png")
    size_kb = sum(p.stat().st_size for p in root.rglob("*") if p.is_file()) // 1024
    (root / "DEBIAN").mkdir()
    (root / "DEBIAN" / "control").write_text(f"""\
Package: diffusor
Version: {v}
Architecture: amd64
Maintainer: Armin Timar <timar.armin@gmail.com>
Installed-Size: {size_kb}
Depends: libc6 (>= 2.35), libegl1, libgl1, libfontconfig1, libxkbcommon0, libxkbcommon-x11-0, libdbus-1-3, libxcb-cursor0
Section: science
Priority: optional
Homepage: https://github.com/ArminTimar/Diffusor
Description: diffusion chronometry of mineral zoning profiles
 Diffusor fits analytical and numerical diffusion models to measured
 concentration profiles in minerals, using diffusion coefficients from the
 published literature, and propagates their uncertainties by Monte Carlo
 sampling to give a timescale with its interval.
""", encoding="utf8", newline="\n")
    for p in [root, *root.rglob("*")]:
        if p.is_dir() and not p.is_symlink():
            p.chmod(0o755)
    run(["dpkg-deb", "--root-owner-group", "--build", root,
         DIST / f"diffusor_{v}_amd64.deb"])


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("steps", nargs="*", choices=STEPS,
                        help="steps to run (default: all, in order)")
    args = parser.parse_args(argv)
    for step in STEPS:
        if step in (args.steps or STEPS):
            globals()[step]()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
