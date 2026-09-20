"""WRITEUP.md must keep up with the code.

The writeup is the framework document: a module-by-module account of what runs
when, and a paper-by-paper account of what every equation, coefficient and
convention is taken from.  It is written by hand, so the only thing a test can
enforce is coverage: every module, coefficient, citation, dataset and buffer
that exists in the code must be named somewhere in the document.

A failure here means the document has fallen behind, not that the code is
wrong.  Fix the document.
"""
from pathlib import Path

import pytest

from diffusor import datasets as ds
from diffusor.coefficients import list_coefficients
from diffusor.references import REFERENCES
from diffusor.thermo import available_buffers

ROOT = Path(__file__).resolve().parents[1]
WRITEUP = ROOT / "WRITEUP.md"

# Modules deliberately not given their own entry: package __init__ files that
# only re-export, and the icon folder.
SKIP_MODULES = {"diffusor/coefficients/__init__.py", "diffusor/dataio/__init__.py",
                "diffusor/fitting/__init__.py", "diffusor/gui/__init__.py",
                "diffusor/minerals/__init__.py", "diffusor/solvers/__init__.py",
                "diffusor/thermo/__init__.py"}


@pytest.fixture(scope="module")
def text() -> str:
    assert WRITEUP.exists(), "WRITEUP.md is missing"
    return WRITEUP.read_text(encoding="utf-8")


def _modules():
    out = []
    for pattern in ("diffusor/**/*.py", "scripts/*.py", "tests/*.py"):
        for p in sorted(ROOT.glob(pattern)):
            rel = p.relative_to(ROOT).as_posix()
            if rel not in SKIP_MODULES:
                out.append(rel)
    return out


def test_every_module_is_described(text):
    missing = [m for m in _modules() if m not in text and Path(m).name not in text]
    assert not missing, ("WRITEUP.md does not mention these modules: "
                         + ", ".join(missing))


def test_every_coefficient_key_appears(text):
    missing = [c.key for c in list_coefficients() if c.key not in text]
    assert not missing, ("WRITEUP.md does not mention these coefficients: "
                         + ", ".join(missing))


def test_every_citation_key_appears(text):
    missing = [k for k in sorted(REFERENCES) if k not in text]
    assert not missing, ("WRITEUP.md section 3 does not cover these references: "
                         + ", ".join(missing))


def test_every_example_dataset_appears(text):
    missing = [d.key for d in ds.DATASETS if d.key not in text]
    assert not missing, ("WRITEUP.md does not mention these datasets: "
                         + ", ".join(missing))


def test_every_buffer_appears(text):
    missing = [b for b in available_buffers() if b not in text]
    assert not missing, ("WRITEUP.md does not mention these buffers: "
                         + ", ".join(missing))


def test_the_headline_counts_are_right(text):
    n_coef = len(list_coefficients())
    n_verified = sum(c.verified for c in list_coefficients())
    n_refs = len(REFERENCES)
    assert f"{n_coef} entries" in text or f"**{n_coef}**" in text, \
        f"the coefficient count ({n_coef}) is not stated in WRITEUP.md"
    assert str(n_verified) in text, f"the verified count ({n_verified}) is not stated"
    assert f"{n_refs} keys" in text, f"the reference count ({n_refs}) is not stated"
    n_mod = len(_modules())
    assert f"{n_mod} Python modules" in text, \
        f"the module count ({n_mod}) is not stated in WRITEUP.md"
