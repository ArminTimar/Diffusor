# PyInstaller recipe for the Diffusor application, used on Windows, macOS and Linux.
# Run it through packaging/build.py, which also makes the installer for the system.
# The result is one folder (dist/Diffusor, or dist/Diffusor.app on macOS) holding
# Python, the libraries and Diffusor, so the user needs no Python of their own.
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(SPECPATH).resolve().parent
sys.path.insert(0, str(ROOT))
from diffusor import __version__  # noqa: E402

ICONS = ROOT / "diffusor" / "gui" / "icons"


def public_files(folder):
    """The files of a folder that git tracks; every file when git is not available.

    A local checkout holds development notes beside the examples that are not
    distributed, and the bundle should carry what the public repository carries.
    """
    try:
        listed = subprocess.run(["git", "ls-files", "-z", folder], cwd=ROOT, check=True,
                                capture_output=True).stdout.decode("utf8").split("\0")
        files = [ROOT / p for p in listed if p]
    except (OSError, subprocess.CalledProcessError):
        files = [p for p in (ROOT / folder).rglob("*")
                 if p.is_file() and "__pycache__" not in p.parts]
    return [(str(p), str(p.parent.relative_to(ROOT))) for p in files]


# Paths inside the bundle mirror the source tree, so diffusor/datasets.py and the
# icon lookups find their files next to the modules exactly as in a checkout.
datas = [(str(p), "diffusor/gui/icons") for p in ICONS.iterdir() if p.is_file()]
datas += public_files("examples")

a = Analysis(
    [str(ROOT / "packaging" / "diffusor_app.py")],
    pathex=[str(ROOT)],
    datas=datas,
    hiddenimports=["matplotlib.backends.backend_qtagg"],
    excludes=["tkinter", "_tkinter", "PyQt5", "PyQt6", "PySide2", "IPython", "pytest",
              "matplotlib.backends.backend_tkagg"],
    noarchive=False,
)
pyz = PYZ(a.pure)

version_file = None
if sys.platform == "win32":
    from PyInstaller.utils.win32.versioninfo import (FixedFileInfo, StringFileInfo,
                                                     StringStruct, StringTable,
                                                     VarFileInfo, VarStruct, VSVersionInfo)
    parts = [int(n) for n in __version__.split(".")][:4]
    numbers = tuple(parts + [0] * (4 - len(parts)))
    info = VSVersionInfo(
        ffi=FixedFileInfo(filevers=numbers, prodvers=numbers),
        kids=[StringFileInfo([StringTable("040904B0", [
                  StringStruct("CompanyName", "Armin Timar"),
                  StringStruct("FileDescription", "Diffusor"),
                  StringStruct("FileVersion", __version__),
                  StringStruct("InternalName", "Diffusor"),
                  StringStruct("LegalCopyright", "MIT License, Armin Timar"),
                  StringStruct("OriginalFilename", "Diffusor.exe"),
                  StringStruct("ProductName", "Diffusor"),
                  StringStruct("ProductVersion", __version__)])]),
              VarFileInfo([VarStruct("Translation", [0x0409, 1200])])])
    version_file = os.path.join(workpath, "version_info.txt")
    with open(version_file, "w", encoding="utf8") as f:
        f.write(str(info))

# macOS: PyInstaller signs every binary with this identity (ad hoc when it is empty)
# and turns on the hardened runtime that notarisation requires.
identity = os.environ.get("DIFFUSOR_MACOS_SIGN_IDENTITY") or None
entitlements = str(ROOT / "packaging" / "macos" / "entitlements.plist") if identity else None

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="Diffusor",
    console=False,
    icon=str(ICONS / ("diffusor.ico" if sys.platform == "win32" else "diffusor.png")),
    version=version_file,
    codesign_identity=identity,
    entitlements_file=entitlements,
)
coll = COLLECT(exe, a.binaries, a.datas, name="Diffusor")

if sys.platform == "darwin":
    app = BUNDLE(
        coll,
        name="Diffusor.app",
        icon=str(ICONS / "diffusor.png"),
        bundle_identifier="io.github.armintimar.diffusor",
        version=__version__,
        info_plist={
            "CFBundleName": "Diffusor",
            "CFBundleDisplayName": "Diffusor",
            "CFBundleShortVersionString": __version__,
            "CFBundleVersion": __version__,
            "NSHighResolutionCapable": True,
            "NSHumanReadableCopyright": "MIT License, Armin Timar",
            "LSApplicationCategoryType": "public.app-category.education",
        },
        codesign_identity=identity,
        entitlements_file=entitlements,
    )
