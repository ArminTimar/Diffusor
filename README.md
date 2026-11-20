# Diffusor

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23074989.svg)](https://doi.org/10.5281/zenodo.23074989)
[![tests](https://github.com/ArminTimar/Diffusor/actions/workflows/tests.yml/badge.svg)](https://github.com/ArminTimar/Diffusor/actions/workflows/tests.yml)

Diffusor is a desktop program for diffusion chronometry of zoned crystals. It
fits a diffusion time to a measured compositional profile and propagates the
uncertainties of temperature, oxygen fugacity, pressure, the diffusion
coefficient and the measurements to that time by Monte Carlo sampling. It
contains 161 published diffusion laws for 14 minerals and for water in silicate
melt. Each law is stored with its source, the equation as printed and its
calibration range.

## Install

Download the file for your system from the **Assets** list of the
[latest release](https://github.com/ArminTimar/Diffusor/releases/latest). The
installers include Python and every library Diffusor needs, so nothing else has
to be installed.

| System | File |
| ------ | ---- |
| Windows 10 and 11 (64-bit) | `Diffusor-<version>-Windows-Setup.exe` |
| macOS, Apple silicon (M1 and later) | `Diffusor-<version>-macOS-AppleSilicon.dmg` |
| macOS, Intel processor | `Diffusor-<version>-macOS-Intel.dmg` |
| Linux, any distribution (x86-64) | `Diffusor-<version>-Linux-x86_64.AppImage` |
| Ubuntu, Debian, Linux Mint | `diffusor_<version>_amd64.deb` |

**Windows.** Run the setup file and follow the wizard. It installs for the
current user and needs no administrator rights. The installer is not yet
code-signed, so SmartScreen shows *Windows protected your PC* the first time.
Click **More info**, then **Run anyway**.

**macOS.** Open the `.dmg` and drag Diffusor onto Applications. Apple has not
notarised the app yet, so macOS blocks the first start. On macOS 15 and later,
open **System Settings > Privacy & Security** and click **Open Anyway**. On
earlier versions, right-click Diffusor in Applications and choose **Open**.

**Linux.** Allow the AppImage to run as a program under **Properties >
Permissions** and double-click it, or install the `.deb` from the software
centre.

Diffusor checks GitHub for a newer release once a day and shows a bar when one
exists. It never downloads anything by itself.

## First run

The window leads through six steps: Data, Mineral, Coefficient, Conditions,
Model and Uncertainty. **Continue** and **Back** at the bottom move between
them, and **Review and fit** on the last step opens the results.

To see a complete run before using your own data, pick an example on the Data
step and click **Load example**. It fills in every step with the conditions of
a published study or a synthetic test case. Step through to see what was set,
then click **Fit**. **Monte Carlo** repeats the fit with sampled inputs and
reports the time as a median with 68 % and 95 % intervals.

Your own profile goes in with **Choose a file...** on the Data step. CSV, TSV, text
and Excel files are read, with one row per analysed point:

```
Distance_um,FeO_wt,MgO_wt,FeO_err,MgO_err
0.0,20.9,22.4,0.15,0.20
1.7,20.8,22.5,0.15,0.20
```

Column names are free. Diffusor suggests which column is the distance and which
are compositions, and you confirm the choice. Two oxide columns are converted to
the molar ratio the Fe-Mg laws are calibrated against. A single column is
modelled as given.

**File > Export results** writes the fitted time, the profiles, the figure and a
methods paragraph that names the diffusion law, solver, boundary conditions and
uncertainty model of the run, followed by its references.

## Running from source

This needs Python 3.10 or newer. Download the source with **Code > Download
ZIP**, unpack it and double-click **Start Diffusor.bat** on Windows or **Start
Diffusor.command** on macOS and Linux. The first start sets up a private Python
environment in a `.venv` folder, which takes a few minutes and needs an internet
connection. As a Python package:

```bash
python -m pip install -e .
python -m diffusor
```

## Documentation

The [user guide](manual/USER_GUIDE.md) covers each step in detail, profiles
from BSE images and element maps, batches of crystals, timescale populations,
multicomponent garnet and isotope profiles, exports and the
[code signing policy](manual/USER_GUIDE.md#code-signing-policy).
[METHODS.md](manual/METHODS.md) describes the solvers, the uncertainty model,
the registry of diffusion laws and the known limitations.
[VERIFICATION.md](manual/VERIFICATION.md) lists the test suite, the comparison
of each law with its source and the refits of published studies. Changes
between versions are in [CHANGELOG.md](CHANGELOG.md).

## Citation

If Diffusor contributes to published work, cite the version you used. It is
stated in every exported methods paragraph.
[doi:10.5281/zenodo.23074989](https://doi.org/10.5281/zenodo.23074989) points to
the newest release, and each version's Zenodo page gives its own DOI. GitHub's
**Cite this repository** button gives the citation in APA and BibTeX format.
Please cite the diffusion laws as well. **Help > All references** exports them
as BibTeX.

## Licence

MIT, see [LICENSE](LICENSE).
