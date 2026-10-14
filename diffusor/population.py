"""Timescale populations: many crystals' times compared in one figure.

A population is a list of :class:`Timescale` records, one per crystal or
profile, each with an optional group (an eruption, a sample, a mineral or a
zoning type). Records come from the current run, from saved results JSON
files (with their Monte Carlo times when those were exported next to them), or
from a table of times typed or collected elsewhere.

Each record's uncertainty becomes a probability density in log10 time:

* with Monte Carlo times, a Gaussian kernel density estimate of log10 t, with
  the rule-of-thumb bandwidth of Silverman (1986, section 3.4.2),
  h = 0.9 min(sd, IQR / 1.34) n^(-1/5);
* with only a 68% interval, two half-normals in log10 t joined at the median,
  each holding half the probability, whose widths put the 16th and 84th
  percentiles where they were given;
* with neither, no density; the record still counts in the ranked and
  cumulative plots.

A group's density is the mean of its members' densities, each integrating to
one, so every crystal weighs the same however narrow its interval (a summed
probability density). Densities are per unit log10 time, so they do not depend
on the unit the axis is drawn in.
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from .constants import SEC_PER_DAY, SEC_PER_HOUR, SEC_PER_MIN, SEC_PER_YEAR
from .references import format_reference
from .thermo.units import human_time

# a record with fewer Monte Carlo times than this uses its percentiles instead
MIN_SAMPLES_FOR_KDE = 20

TIME_UNITS: Dict[str, float] = {
    "s": 1.0, "min": SEC_PER_MIN, "h": SEC_PER_HOUR, "d": SEC_PER_DAY,
    "yr": SEC_PER_YEAR, "kyr": 1.0e3 * SEC_PER_YEAR, "Myr": 1.0e6 * SEC_PER_YEAR,
}
UNIT_NAMES = {"s": "seconds", "min": "minutes", "h": "hours", "d": "days",
              "yr": "years", "kyr": "thousand years", "Myr": "million years"}
_UNIT_ALIASES = {
    "s": "s", "sec": "s", "secs": "s", "second": "s", "seconds": "s",
    "min": "min", "mins": "min", "minute": "min", "minutes": "min",
    "h": "h", "hr": "h", "hrs": "h", "hour": "h", "hours": "h",
    "d": "d", "day": "d", "days": "d",
    "y": "yr", "yr": "yr", "yrs": "yr", "a": "yr", "year": "yr", "years": "yr",
    "kyr": "kyr", "ka": "kyr", "myr": "Myr", "ma": "Myr",
}


def unit_seconds(unit: str) -> float:
    """Seconds in one ``unit``; accepts the usual spellings (days, yr, a, ka...)."""
    key = _UNIT_ALIASES.get(str(unit).strip().lower())
    if key is None:
        raise ValueError(f"unknown time unit {unit!r}; use s, min, h, d, yr, kyr or Myr")
    return TIME_UNITS[key]


@dataclass
class Timescale:
    """One crystal's or profile's time. All times are in seconds."""
    label: str
    t_best: float
    group: str = ""
    p16: Optional[float] = None
    p84: Optional[float] = None
    p2_5: Optional[float] = None
    p97_5: Optional[float] = None
    median: Optional[float] = None
    samples: Optional[np.ndarray] = None
    coefficient: str = ""
    T_C: Optional[float] = None
    source: str = ""
    notes: List[str] = field(default_factory=list)

    def __post_init__(self):
        self.t_best = float(self.t_best)
        if not np.isfinite(self.t_best) or self.t_best <= 0:
            raise ValueError(f"{self.label}: the time must be a positive number, not {self.t_best}")
        if self.samples is not None:
            s = np.asarray(self.samples, dtype=float)
            s = s[np.isfinite(s) & (s > 0)]
            self.samples = s if s.size else None
        if self.samples is not None and self.samples.size >= MIN_SAMPLES_FOR_KDE:
            # percentiles follow the samples, so the table and the density agree
            q = np.percentile(self.samples, [2.5, 16, 50, 84, 97.5])
            self.p2_5, self.p16, self.median, self.p84, self.p97_5 = (float(v) for v in q)

    @property
    def has_interval(self) -> bool:
        return self.p16 is not None and self.p84 is not None and self.p84 > self.p16 > 0

    @property
    def density_basis(self) -> str:
        if self.samples is not None and self.samples.size >= MIN_SAMPLES_FOR_KDE:
            return f"kernel density of {self.samples.size} Monte Carlo times"
        if self.has_interval:
            return "two half-normals in log10 t through the 16th and 84th percentiles"
        return "none (no uncertainty)"

    def log_density(self, grid_log10: np.ndarray) -> Optional[np.ndarray]:
        """Probability density per unit log10(t / s) on ``grid_log10``, or None."""
        g = np.asarray(grid_log10, dtype=float)
        if self.samples is not None and self.samples.size >= MIN_SAMPLES_FOR_KDE:
            z = np.log10(self.samples)
            iqr = np.subtract(*np.percentile(z, [75, 25]))
            spread = min(np.std(z, ddof=1), iqr / 1.34) if iqr > 0 else np.std(z, ddof=1)
            h = max(0.9 * spread * z.size ** -0.2, 0.01)
            u = (g[:, None] - z[None, :]) / h
            return np.exp(-0.5 * u * u).sum(axis=1) / (z.size * h * np.sqrt(2 * np.pi))
        if self.has_interval:
            mid = self.median if self.median is not None else self.t_best
            mu = np.log10(mid)
            lo = max(mu - np.log10(self.p16), 1e-3)
            hi = max(np.log10(self.p84) - mu, 1e-3)
            # half the probability on each side of the median, so the 16th and 84th
            # percentiles fall exactly where they were given; with an asymmetric
            # interval the density steps at the median
            s = np.where(g < mu, lo, hi)
            return np.exp(-0.5 * ((g - mu) / s) ** 2) / (np.sqrt(2 * np.pi) * s)
        return None


# --- building records ------------------------------------------------------------
def from_run(fit_result, mc_result=None, label: str = "", group: str = "",
             source: str = "") -> Timescale:
    """A record from the fit (and Monte Carlo) currently in the main window."""
    model = fit_result.model
    return Timescale(
        label=label or "Run", group=group, t_best=fit_result.t_seconds,
        samples=None if mc_result is None else np.asarray(mc_result.times, dtype=float),
        coefficient=model.coefficient.key, T_C=model.conditions.T_K - 273.15, source=source)


def load_result_json(path) -> Timescale:
    """A record from a ``*_results.json`` written by File > Export results.

    The Monte Carlo times are read from the ``*_montecarlo_times.csv`` beside it
    when that was exported too; otherwise the percentiles in the JSON are used.
    """
    path = Path(path)
    d = json.loads(path.read_text(encoding="utf-8"))
    if "profiles" in d and "fit" not in d:
        raise ValueError(f"{path.name} is a shared-duration study, which has one time for all "
                         "its profiles; add the profiles' own results files instead")
    if "fit" not in d or "t_seconds" not in d["fit"]:
        raise ValueError(f"{path.name} is not a Diffusor results file (no fitted time)")
    stem = path.name[:-len("_results.json")] if path.name.endswith("_results.json") else path.stem
    data = d.get("data") or {}
    label = Path(data["source"]).stem if data.get("source") else stem
    mc = d.get("montecarlo") or {}
    samples, notes = None, []
    csv = path.with_name(f"{stem}_montecarlo_times.csv")
    if csv.exists():
        import pandas as pd
        samples = pd.read_csv(csv)["t_seconds"].to_numpy(dtype=float)
        notes.append(f"Monte Carlo times from {csv.name}")
    cond = d.get("conditions") or {}
    return Timescale(
        label=label, t_best=d["fit"]["t_seconds"],
        p16=mc.get("p16_s"), p84=mc.get("p84_s"), p2_5=mc.get("p2.5_s"),
        p97_5=mc.get("p97.5_s"), median=mc.get("median_s"), samples=samples,
        coefficient=(d.get("coefficient") or {}).get("key", ""),
        T_C=cond.get("T_C"), source=str(path), notes=notes)


def load_setup_file(path, group: str = "") -> List[Timescale]:
    """Records from a crystal setup or batch file (``.diffusor``): every crystal with a saved result.

    A crystal's own group is kept; crystals without one get ``group``.
    """
    from .crystal_setup import crystals_in, read_file
    path = Path(path)
    out = []
    for c in crystals_in(read_file(path)):
        r = c.result
        if not r:
            continue
        mc = r.get("montecarlo") or {}
        out.append(Timescale(
            label=c.label, group=c.group or group, t_best=r["t_seconds"],
            samples=np.asarray(mc["times_s"], dtype=float) if mc.get("times_s") else None,
            coefficient=r.get("coefficient", ""), T_C=c.settings.get("T_C"),
            source=str(path), notes=[f"from {path.name}"]))
    if not out:
        raise ValueError(f"{path.name} holds no fitted result; fit the crystals and save again")
    return out


def _norm(name) -> str:
    return re.sub(r"[^a-z0-9.]+", "_", str(name).strip().lower()).strip("_")


_LABEL = ("label", "crystal", "name", "id", "profile", "sample")
_GROUP = ("group", "eruption", "population", "event", "unit")
_TIME = ("t_best", "best", "time", "t", "timescale", "duration", "median")
_LOWER = ("p16", "lower", "lo", "low", "min_68")
_UPPER = ("p84", "upper", "hi", "high", "max_68")
_LOWER95 = ("p2.5", "p2_5", "lower_95", "lo95")
_UPPER95 = ("p97.5", "p97_5", "upper_95", "hi95")
_SIGMA_LOG = ("sigma_log10", "sigma_log", "log_sigma", "1s_log10")


def _split_unit(header: str) -> Tuple[str, Optional[str]]:
    """'Time (days)' -> ('time', 'd'); 'p16_yr' -> ('p16', 'yr'); 'time' -> ('time', None)."""
    m = re.match(r"^(.*?)\s*[\(\[]\s*([A-Za-z]+)\s*[\)\]]\s*$", str(header))
    if m and m.group(2).lower() in _UNIT_ALIASES:
        return _norm(m.group(1)), _UNIT_ALIASES[m.group(2).lower()]
    n = _norm(header)
    if "_" in n:
        base, last = n.rsplit("_", 1)
        if last in _UNIT_ALIASES and base:
            return base, _UNIT_ALIASES[last]
    return n, None


def load_table(path, default_unit: str = "d", group: str = "") -> List[Timescale]:
    """Records from a CSV or Excel table, one row per crystal.

    Columns are recognised by name, ignoring case: a label (label, crystal,
    name, id), an optional group (group, eruption, population, event), the time
    (time, t, best, median...), and optionally a 68% interval (p16 and p84, or
    lower and upper), a 95% interval (p2.5 and p97.5) or a 1 sigma in log10
    units (sigma_log10). A unit may follow a column name, as in ``time (days)``
    or ``p16_yr``, or come from a ``unit`` column; otherwise ``default_unit``.
    """
    import pandas as pd
    path = Path(path)
    if path.suffix.lower() in (".xlsx", ".xlsm", ".xls"):
        df = pd.read_excel(path)
    else:
        df = pd.read_csv(path, sep=None, engine="python")
    cols: Dict[str, Tuple[str, Optional[str]]] = {}
    for c in df.columns:
        base, unit = _split_unit(c)
        cols.setdefault(base, (c, unit))
    unit_col = cols.get("unit", (None, None))[0] if "unit" in cols else None

    def pick(names):
        for n in names:
            if n in cols and cols[n][0] != unit_col:
                return cols[n]
        return None, None

    time_col, time_unit = pick(_TIME)
    if time_col is None:
        raise ValueError(f"{path.name}: no time column. Name one 'time', 't' or 'median', "
                         f"for example 'time (days)'. Found: {', '.join(map(str, df.columns))}")
    label_col, _ = pick(_LABEL)
    group_col, _ = pick([g for g in _GROUP if not (g == "unit" and unit_col)])
    lo_col, lo_unit = pick(_LOWER)
    hi_col, hi_unit = pick(_UPPER)
    lo95_col, lo95_unit = pick(_LOWER95)
    hi95_col, hi95_unit = pick(_UPPER95)
    sig_col, _ = pick(_SIGMA_LOG)

    out: List[Timescale] = []
    for i, row in df.iterrows():
        t = pd.to_numeric(row[time_col], errors="coerce")
        if not np.isfinite(t):
            continue
        row_unit = str(row[unit_col]) if unit_col is not None and pd.notna(row[unit_col]) else None

        def sec(col, unit):
            if col is None or pd.isna(row[col]):
                return None
            v = pd.to_numeric(row[col], errors="coerce")
            if not np.isfinite(v):
                return None
            return float(v) * unit_seconds(unit or row_unit or default_unit)

        t_s = sec(time_col, time_unit)
        p16, p84 = sec(lo_col, lo_unit), sec(hi_col, hi_unit)
        if (p16 is None or p84 is None) and sig_col is not None and pd.notna(row[sig_col]):
            s = float(row[sig_col])
            p16, p84 = t_s * 10 ** -s, t_s * 10 ** s
        label = str(row[label_col]) if label_col is not None and pd.notna(row[label_col]) else f"Row {i + 2}"
        grp = str(row[group_col]) if group_col is not None and pd.notna(row[group_col]) else group
        out.append(Timescale(label=label, group=grp, t_best=t_s, p16=p16, p84=p84,
                             p2_5=sec(lo95_col, lo95_unit), p97_5=sec(hi95_col, hi95_unit),
                             source=f"{path}, row {i + 2}"))
    if not out:
        raise ValueError(f"{path.name}: no row has a numeric time in column {time_col!r}")
    return out


# --- statistics ------------------------------------------------------------------
def groups_of(records: Sequence[Timescale]) -> List[str]:
    """Group names in order of first appearance; ungrouped records form ''."""
    return list(dict.fromkeys(r.group for r in records))


def group_name(group: str) -> str:
    return group or "Ungrouped"


def log_grid(records: Sequence[Timescale], n: int = 600) -> np.ndarray:
    """A log10(t / s) grid covering every record with room for its tails."""
    lows, highs = [], []
    for r in records:
        lows.append(np.log10(min(v for v in (r.t_best, r.p2_5, r.p16) if v)))
        highs.append(np.log10(max(v for v in (r.t_best, r.p97_5, r.p84) if v)))
        if r.has_interval:
            # four sigma of each half, so the tails of the density are not cut
            mu = np.log10(r.median if r.median is not None else r.t_best)
            lows.append(mu - 4 * (mu - np.log10(r.p16)))
            highs.append(mu + 4 * (np.log10(r.p84) - mu))
        if r.samples is not None:
            lows.append(np.log10(r.samples.min()))
            highs.append(np.log10(r.samples.max()))
    lo, hi = min(lows), max(highs)
    pad = max(0.3, 0.05 * (hi - lo))
    return np.linspace(lo - pad, hi + pad, n)


def group_density(records: Sequence[Timescale], grid_log10: np.ndarray) -> Tuple[Optional[np.ndarray], int]:
    """Mean density of the records that have one, and how many did."""
    curves = [d for d in (r.log_density(grid_log10) for r in records) if d is not None]
    if not curves:
        return None, 0
    return np.mean(curves, axis=0), len(curves)


def summary(records: Sequence[Timescale]) -> List[Dict]:
    """One row per group: count, median and spread of the best-fit times."""
    rows = []
    for g in groups_of(records):
        t = np.array([r.t_best for r in records if r.group == g])
        q16, q50, q84 = np.percentile(t, [16, 50, 84])
        rows.append({"group": group_name(g), "n": int(t.size),
                     "n_with_uncertainty": sum(r.has_interval for r in records if r.group == g),
                     "median_s": float(q50), "p16_of_population_s": float(q16),
                     "p84_of_population_s": float(q84), "min_s": float(t.min()),
                     "max_s": float(t.max())})
    return rows


def auto_unit(records: Sequence[Timescale]) -> str:
    """The unit in which the population's median reads most naturally."""
    m = float(np.median([r.t_best for r in records]))
    for unit in ("Myr", "kyr", "yr", "d", "h", "min"):
        if m >= TIME_UNITS[unit]:
            return unit
    return "s"


def describe(records: Sequence[Timescale]) -> str:
    """A methods paragraph for the population figure."""
    n_kde = sum(r.samples is not None and r.samples.size >= MIN_SAMPLES_FOR_KDE for r in records)
    n_split = sum(r.has_interval for r in records) - n_kde
    n_none = len(records) - n_kde - n_split
    lines = [f"{len(records)} timescales in {len(groups_of(records))} group(s). "
             "Markers are best-fit times; bars are 68% (thick) and 95% (thin) intervals.",
             "Each group's density is the mean of its members' densities in log10 time, each "
             "normalised to one, so every crystal has equal weight (summed-probability plot).",
             f"{n_kde} densities are Gaussian kernel density estimates of the Monte Carlo times "
             "(rule-of-thumb bandwidth of Silverman 1986); "
             f"{n_split} are two half-normals in log10 time joined at the median, through the "
             "16th and 84th percentiles; "
             f"{n_none} records without uncertainty appear only in the ranked and cumulative plots.",
             "The cumulative curves are the empirical distributions of the best-fit times."]
    for row in summary(records):
        lines.append(f"{row['group']}: n = {row['n']}, median {human_time(row['median_s'])}, "
                     f"16th-84th percentile of the population {human_time(row['p16_of_population_s'])} "
                     f"to {human_time(row['p84_of_population_s'])}.")
    if n_kde:
        lines += ["", "References", f"- {format_reference('silverman1986')}"]
    return "\n".join(lines)


# --- export ----------------------------------------------------------------------
def save_population(directory, records: Sequence[Timescale], figure=None,
                    basename: str = "population") -> Dict[str, str]:
    """Write the records, group summary, densities, methods text and the figure."""
    import pandas as pd
    out = Path(directory)
    out.mkdir(parents=True, exist_ok=True)
    written: Dict[str, str] = {}
    yr = SEC_PER_YEAR
    rows = [{"label": r.label, "group": group_name(r.group),
             "t_best_s": r.t_best, "t_best_days": r.t_best / SEC_PER_DAY, "t_best_years": r.t_best / yr,
             "median_s": r.median, "p16_s": r.p16, "p84_s": r.p84, "p2.5_s": r.p2_5, "p97.5_s": r.p97_5,
             "density_basis": r.density_basis, "coefficient": r.coefficient, "T_C": r.T_C,
             "source": r.source} for r in records]
    p = out / f"{basename}_timescales.csv"
    pd.DataFrame(rows).to_csv(p, index=False)
    written["timescales"] = str(p)
    p = out / f"{basename}_groups.csv"
    pd.DataFrame(summary(records)).to_csv(p, index=False)
    written["groups"] = str(p)
    grid = log_grid(records)
    dens = {"log10_t_s": grid, "t_s": 10 ** grid}
    for g in groups_of(records):
        d, _ = group_density([r for r in records if r.group == g], grid)
        if d is not None:
            dens[f"density_{group_name(g)}"] = d
    p = out / f"{basename}_density.csv"
    pd.DataFrame(dens).to_csv(p, index=False)
    written["density"] = str(p)
    p = out / f"{basename}_methods.txt"
    p.write_text(describe(records) + "\n", encoding="utf-8")
    written["methods"] = str(p)
    if figure is not None:
        from .dataio.vector import save_vector
        f = out / f"{basename}_figure.png"
        figure.savefig(f, dpi=300, bbox_inches="tight")
        written["figure"] = str(f)
        for suffix in ("svg", "pdf"):
            f = out / f"{basename}_figure.{suffix}"
            save_vector(figure, f, bbox_inches="tight")
            written[f"figure_{suffix}"] = str(f)
    return written


def extend_unique(records: List[Timescale], new: Iterable[Timescale]) -> List[str]:
    """Append ``new`` to ``records``, renaming clashing labels; return the labels added."""
    taken = {r.label for r in records}
    added = []
    for r in new:
        base, k = r.label, 2
        while r.label in taken:
            r.label = f"{base} ({k})"
            k += 1
        taken.add(r.label)
        records.append(r)
        added.append(r.label)
    return added
