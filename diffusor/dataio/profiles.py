"""Loading measured profiles.

The expected input is a table with a distance column and one or two
composition columns, optionally with uncertainty columns:

    distance_um, FeO_wt, MgO_wt, FeO_err, MgO_err

CSV, TSV and Excel are supported.  Nothing is guessed silently: the loader
reports which columns it matched, and :class:`ProfileSpec` records the choice
so it can be written into the exported methods block.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from ..thermo.units import composition_variable, length_to_m

DISTANCE_HINTS = ("distance", "dist", "x", "micron", "um", "µm", "position", "real_distance")
ERROR_HINTS = ("err", "error", "sigma", "sd", "std", "2s", "1s", "stdev", "standard_error")


@dataclass
class ProfileSpec:
    """How the columns of the file map onto the model."""
    distance_column: str
    column_a: str
    column_b: Optional[str] = None
    sigma_a_column: Optional[str] = None
    sigma_b_column: Optional[str] = None
    distance_unit: str = "um"
    mode: str = "A/(A+B)"            # see thermo.units.composition_variable
    oxide_a: Optional[str] = None    # e.g. 'FeO' to convert wt% to cation moles
    oxide_b: Optional[str] = None
    sigma_level: str = "1s"
    label: str = ""

    def describe(self) -> str:
        s = f"x = '{self.distance_column}' [{self.distance_unit}], "
        s += f"C = {self.mode} of '{self.column_a}'"
        if self.column_b:
            s += f" and '{self.column_b}'"
        if self.oxide_a and self.oxide_b:
            s += f" (as {self.oxide_a}/{self.oxide_b} wt%, converted to cation moles)"
        if self.sigma_a_column:
            s += f", uncertainty from '{self.sigma_a_column}' ({self.sigma_level})"
        return s


@dataclass
class Profile:
    x: np.ndarray                     # micrometres
    C: np.ndarray                     # modelled composition variable
    sigma: Optional[np.ndarray] = None
    raw: Optional[pd.DataFrame] = None
    spec: Optional[ProfileSpec] = None
    source: str = ""
    notes: List[str] = field(default_factory=list)

    def __len__(self) -> int:
        return int(self.x.size)

    def sorted(self) -> "Profile":
        idx = np.argsort(self.x)
        return Profile(self.x[idx], self.C[idx],
                       None if self.sigma is None else self.sigma[idx],
                       self.raw, self.spec, self.source, list(self.notes))

    def zeroed(self) -> "Profile":
        """Shift x so the traverse starts at zero."""
        p = self.sorted()
        p.x = p.x - p.x[0]
        return p


def read_table(path) -> pd.DataFrame:
    p = Path(path)
    suf = p.suffix.lower()
    if suf in (".xlsx", ".xlsm", ".xls"):
        return pd.read_excel(p)
    if suf in (".tsv", ".tab"):
        return pd.read_csv(p, sep="\t")
    return pd.read_csv(p, sep=None, engine="python")


def _match(columns: Sequence[str], hints: Sequence[str]) -> Optional[str]:
    low = {c: str(c).strip().lower() for c in columns}
    for c, lc in low.items():
        if lc in hints:
            return c
    for c, lc in low.items():
        if any(h in lc for h in hints):
            return c
    return None


def suggest_spec(df: pd.DataFrame) -> ProfileSpec:
    """A first guess at the column mapping; always shown to the user for confirmation."""
    cols = list(df.columns)
    numeric = [c for c in cols if pd.api.types.is_numeric_dtype(df[c])]
    dist = _match(numeric, DISTANCE_HINTS) or (numeric[0] if numeric else cols[0])
    rest = [c for c in numeric if c != dist]
    errs = [c for c in rest if any(h in str(c).lower() for h in ERROR_HINTS)]
    vals = [c for c in rest if c not in errs]
    a = vals[0] if vals else None
    b = vals[1] if len(vals) > 1 else None
    return ProfileSpec(distance_column=dist, column_a=a, column_b=b,
                       sigma_a_column=errs[0] if errs else None,
                       sigma_b_column=errs[1] if len(errs) > 1 else None,
                       mode="A/(A+B)" if b else "A")


def build_profile(df: pd.DataFrame, spec: ProfileSpec, source: str = "") -> Profile:
    """Apply a :class:`ProfileSpec` to a table and return the modelled profile."""
    x_raw = np.asarray(df[spec.distance_column], dtype=float)
    x = length_to_m(x_raw, spec.distance_unit) * 1.0e6          # work in micrometres
    a = np.asarray(df[spec.column_a], dtype=float)
    b = np.asarray(df[spec.column_b], dtype=float) if spec.column_b else None
    C = composition_variable(a, b, spec.mode, spec.oxide_a, spec.oxide_b)

    sigma = None
    notes: List[str] = []
    if spec.sigma_a_column:
        sa = np.asarray(df[spec.sigma_a_column], dtype=float)
        if spec.sigma_level == "2s":
            sa = sa / 2.0
            notes.append("uncertainty column interpreted as 2 sigma and halved")
        if spec.mode in ("A/(A+B)", "B/(A+B)") and b is not None:
            # propagate to the ratio; if only one uncertainty column is given the
            # other is assumed to have the same relative uncertainty
            sb = (np.asarray(df[spec.sigma_b_column], dtype=float)
                  if spec.sigma_b_column else sa / np.maximum(np.abs(a), 1e-30) * np.abs(b))
            if spec.sigma_b_column and spec.sigma_level == "2s":
                sb = sb / 2.0
            tot = a + b
            with np.errstate(divide="ignore", invalid="ignore"):
                # d(a/(a+b)) = (b da - a db)/(a+b)^2, errors added in quadrature
                sigma = np.sqrt((b * sa) ** 2 + (a * sb) ** 2) / tot ** 2
            if spec.mode == "B/(A+B)":
                pass                       # the ratio uncertainty is symmetric
            if not spec.sigma_b_column:
                notes.append("only one uncertainty column was given; the second element was "
                             "assumed to carry the same relative uncertainty")
        else:
            sigma = sa
    if sigma is not None:
        bad = ~np.isfinite(sigma) | (sigma <= 0)
        if np.any(bad):
            good = sigma[~bad]
            fill = float(np.median(good)) if good.size else float(np.std(C)) or 1e-6
            sigma = np.where(bad, fill, sigma)
            notes.append(f"{int(bad.sum())} non-positive or missing uncertainties replaced by "
                         f"the median ({fill:.3g})")

    ok = np.isfinite(x) & np.isfinite(C)
    if not np.all(ok):
        notes.append(f"{int((~ok).sum())} rows with missing values were dropped")
    prof = Profile(x=x[ok], C=C[ok], sigma=None if sigma is None else sigma[ok],
                   raw=df, spec=spec, source=source, notes=notes)
    return prof.sorted()


def load_profile(path, spec: Optional[ProfileSpec] = None) -> Profile:
    df = read_table(path)
    spec = spec or suggest_spec(df)
    return build_profile(df, spec, source=str(path))
