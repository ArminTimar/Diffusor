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
import re
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from ..constants import OXIDE_CATIONS, OXIDE_MOLAR_MASS
from ..thermo.units import composition_variable, length_to_m


def composition_sigma(a, b, sa, sb, mode: str,
                      oxide_a: Optional[str] = None, oxide_b: Optional[str] = None):
    """1-sigma uncertainty of the modelled variable built by ``composition_variable``.

    First-order propagation for independent inputs (JCGM 100:2008, the Guide to
    the expression of uncertainty in measurement, eq. 10 of clause 5.1.2).  With
    oxides named, the modelled variable is formed from cation moles,
    ``A = n_a a / M_a`` and ``B = n_b b / M_b`` (a, b in wt%, M the oxide molar
    mass, n the number of cations per formula unit), so the sensitivity
    coefficients with respect to the measured wt% are

        d[A/(A+B)]/da =  k_a B / (A+B)^2        k = n / M
        d[A/(A+B)]/db = -k_b A / (A+B)^2

    and the combined standard uncertainty is ``sqrt((k_a B sa)^2 + (k_b A sb)^2)
    / (A+B)^2``.  Without oxides k_a = k_b = 1 (the raw columns are the moles).
    ``B/(A+B)`` has the same magnitude (it is 1 - A/(A+B)); ``A-B`` gives
    ``sqrt((k_a sa)^2 + (k_b sb)^2)``; mode ``A`` returns ``sa``.
    """
    sa = np.asarray(sa, dtype=float)
    if mode == "A" or b is None:
        return sa
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    sb = np.asarray(sb, dtype=float)
    ka = kb = 1.0
    if oxide_a and oxide_b:
        ka = OXIDE_CATIONS[oxide_a] / OXIDE_MOLAR_MASS[oxide_a]
        kb = OXIDE_CATIONS[oxide_b] / OXIDE_MOLAR_MASS[oxide_b]
    if mode == "A-B":
        return np.sqrt((ka * sa) ** 2 + (kb * sb) ** 2)
    A, B = ka * a, kb * b
    tot = A + B
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.sqrt((ka * B * sa) ** 2 + (kb * A * sb) ** 2) / tot ** 2



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
    x_min: Optional[float] = None    # keep only points with x_min <= x <= x_max,
    x_max: Optional[float] = None    # in distance_unit; None means no limit

    def describe(self) -> str:
        s = f"x = '{self.distance_column}' [{self.distance_unit}], "
        s += f"C = {self.mode} of '{self.column_a}'"
        if self.column_b:
            s += f" and '{self.column_b}'"
        if self.oxide_a and self.oxide_b:
            s += f" (as {self.oxide_a}/{self.oxide_b} wt%, converted to cation moles)"
        if self.sigma_a_column:
            s += f", uncertainty from '{self.sigma_a_column}' ({self.sigma_level})"
        if self.x_min is not None or self.x_max is not None:
            lo = "start" if self.x_min is None else f"{self.x_min:g}"
            hi = "end" if self.x_max is None else f"{self.x_max:g}"
            s += f", fitted window {lo} to {hi} {self.distance_unit}"
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
    row_index: Optional[np.ndarray] = None   # rows of ``raw`` behind each point

    def __len__(self) -> int:
        return int(self.x.size)

    def sorted(self) -> "Profile":
        idx = np.argsort(self.x, kind="stable")
        return Profile(self.x[idx], self.C[idx],
                       None if self.sigma is None else self.sigma[idx],
                       self.raw, self.spec, self.source, list(self.notes),
                       None if self.row_index is None else self.row_index[idx])

    def column(self, name: str) -> np.ndarray:
        """Another column of the source table, aligned point for point with x.

        Use this for anything that must follow the same rows as the profile, such
        as an anorthite column, so that windowing and sorting can never put it out
        of step with the modelled composition.
        """
        if self.raw is None or self.row_index is None:
            raise ValueError("this profile does not keep its source table")
        return np.asarray(self.raw[name], dtype=float)[self.row_index]

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


DISTANCE_WORDS = ("distance", "dist", "position", "pos", "x", "depth", "micron", "microns")
ERROR_WORDS = ("err", "error", "sigma", "sd", "std", "stdev", "1s", "2s", "unc", "uncertainty",
               "standard_error", "se")
UNIT_WORDS = {"um": "um", "µm": "um", "micron": "um", "microns": "um", "mum": "um",
              "mm": "mm", "nm": "nm", "m": "m"}
# pairs that are normally modelled as a molar ratio, most common first
OXIDE_PAIRS = (("FeO", "MgO"), ("CaO", "Na2O"), ("FeO", "MnO"))
KNOWN_OXIDES = ("SiO2", "TiO2", "Al2O3", "Cr2O3", "FeO", "Fe2O3", "MnO", "MgO", "NiO", "CaO",
                "Na2O", "K2O", "Li2O", "SrO", "BaO")


CATION_WORDS = ("fe", "mg", "ca", "na", "mn")
CATION_PAIRS = (("fe", "mg"), ("ca", "na"), ("fe", "mn"))
FRACTION_WORDS = ("xan", "an", "xfe", "xmg", "fo", "fa", "en", "fs", "mg#", "mgnumber", "grey",
                  "gray", "greyvalue", "grayvalue")


def _tokens(name) -> List[str]:
    return [t for t in re.split(r"[^0-9a-zµ]+", str(name).strip().lower()) if t]


def _oxide_of(name) -> Optional[str]:
    """The oxide a column holds, from names like 'FeO', 'FeO_wt' or 'MgO (wt%)'."""
    low = str(name).strip().lower()
    for ox in sorted(KNOWN_OXIDES, key=len, reverse=True):
        o = ox.lower()
        if low == o or (low.startswith(o) and not low[len(o)].isalnum()):
            return ox
    return None


def _is_error(name) -> bool:
    return any(t in ERROR_WORDS for t in _tokens(name))


def _unit_of(name) -> Optional[str]:
    for t in reversed(_tokens(name)):
        if t in UNIT_WORDS:
            return UNIT_WORDS[t]
    if "µm" in str(name):
        return "um"
    return None


def guess_distance_unit(name) -> str:
    return _unit_of(name) or "um"


def suggest_spec(df: pd.DataFrame) -> ProfileSpec:
    """A first guess at the column mapping, shown to the user for confirmation.

    Distance is the first numeric column named like a distance. If the table
    holds a common oxide pair (FeO and MgO, CaO and Na2O) that pair is used as a
    molar ratio with the oxides named. Otherwise a table with exactly two
    composition columns is modelled as their ratio, and anything else models
    the first composition column on its own. Uncertainty columns are paired to
    their value column by name.
    """
    cols = list(df.columns)
    numeric = [c for c in cols if pd.api.types.is_numeric_dtype(df[c])]
    if not numeric:
        raise ValueError("the table has no numeric columns")
    dist = next((c for c in numeric if any(t in DISTANCE_WORDS for t in _tokens(c))), None)
    dist = dist or numeric[0]
    rest = [c for c in numeric if c != dist]
    errs = [c for c in rest if _is_error(c)]
    vals = [c for c in rest if c not in errs]

    by_oxide = {}
    for c in vals:
        ox = _oxide_of(c)
        if ox and ox not in by_oxide:
            by_oxide[ox] = c
    a = b = ox_a = ox_b = None
    for pa, pb in OXIDE_PAIRS:
        if pa in by_oxide and pb in by_oxide:
            a, b, ox_a, ox_b = by_oxide[pa], by_oxide[pb], pa, pb
            break
    if a is None:
        by_cation = {}
        for c in vals:
            el = _tokens(c)[0] if _tokens(c) else ""
            if el in CATION_WORDS and el not in by_cation:
                by_cation[el] = c
        for pa, pb in CATION_PAIRS:
            if pa in by_cation and pb in by_cation:
                a, b = by_cation[pa], by_cation[pb]
                break
    if a is None and vals:
        a = vals[0]
        # two bare columns with no unit or name clue are taken as elements A and B
        if len(vals) == 2 and all(len(_tokens(c)) <= 1 and _tokens(c)[0] not in FRACTION_WORDS
                                  for c in vals):
            b = vals[1]

    def error_for(col):
        if col is None or not _tokens(col):
            return None
        stem = _tokens(col)[0]
        for e in errs:
            t = _tokens(e)
            if t and (t[0] == stem or stem.startswith(t[0]) or t[0].startswith(stem)):
                return e
        return None

    sa, sb = error_for(a), error_for(b)
    level = "2s" if any("2s" in _tokens(e) for e in (sa, sb) if e) else "1s"
    return ProfileSpec(distance_column=dist, column_a=a, column_b=b,
                       sigma_a_column=sa, sigma_b_column=sb,
                       distance_unit=guess_distance_unit(dist),
                       mode="A/(A+B)" if b else "A",
                       oxide_a=ox_a, oxide_b=ox_b, sigma_level=level)


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
        if spec.mode in ("A/(A+B)", "B/(A+B)", "A-B") and b is not None:
            # propagate to the modelled variable; if only one uncertainty column is
            # given the other is assumed to have the same relative uncertainty
            sb = (np.asarray(df[spec.sigma_b_column], dtype=float)
                  if spec.sigma_b_column else sa / np.maximum(np.abs(a), 1e-30) * np.abs(b))
            if spec.sigma_b_column and spec.sigma_level == "2s":
                sb = sb / 2.0
            # first-order propagation through the same wt% -> cation-mole conversion
            # that builds C (JCGM 100:2008 eq. 10); see composition_sigma
            sigma = composition_sigma(a, b, sa, sb, spec.mode, spec.oxide_a, spec.oxide_b)
            if not spec.sigma_b_column:
                notes.append("only one uncertainty column was given. The second element was "
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
    if spec.x_min is not None or spec.x_max is not None:
        lo = -np.inf if spec.x_min is None else float(length_to_m(spec.x_min, spec.distance_unit)) * 1.0e6
        hi = np.inf if spec.x_max is None else float(length_to_m(spec.x_max, spec.distance_unit)) * 1.0e6
        inside = (x >= lo) & (x <= hi)
        n_out = int((ok & ~inside).sum())
        ok = ok & inside
        if n_out:
            notes.append(f"{n_out} points outside the fitted window were excluded")
        if not np.any(ok):
            raise ValueError("no points fall inside the fitted window")
    rows = np.arange(len(df))
    prof = Profile(x=x[ok], C=C[ok], sigma=None if sigma is None else sigma[ok],
                   raw=df, spec=spec, source=source, notes=notes, row_index=rows[ok])
    return prof.sorted()


def load_profile(path, spec: Optional[ProfileSpec] = None) -> Profile:
    df = read_table(path)
    spec = spec or suggest_spec(df)
    return build_profile(df, spec, source=str(path))
