"""Crystal setups: the settings of one run, its data, and how the model is built from them.

A setup is everything the main window needs to repeat a run: the profile table with
its column mapping and the points cut by hand, and the settings of every step as a
flat dictionary (:data:`DEFAULTS` lists the keys). The window turns its widgets into
such a dictionary and back; :func:`build_model`, :func:`free_parameters` and
:func:`budget` turn the dictionary and the data into the model, the fitted
parameters and the Monte Carlo budget. The batch study and the setup files use the
same functions, so a crystal fitted from a saved setup or a batch row is fitted
exactly as it would be in the main window.

Setups are saved as JSON with the extension ``.diffusor``. One file holds one crystal
(``"format": "diffusor-setup"``) or a batch of crystals sharing settings
(``"format": "diffusor-batch"``). A saved result (fitted time, Monte Carlo times) is
kept with the setup, so a timescale population can be built from the files.
"""
from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass, field, fields
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from . import __version__
from . import datasets as ds
from .coefficients import Conditions, get as get_coefficient
from .coefficients.plagioclase import ACTIVITY_SPECIES, DEFAULT_ACTIVITY_SET, activity_A
from .dataio.profiles import Profile, ProfileSpec, build_profile
from .fitting import DiffusionModel, UncertaintyBudget
from .minerals import MINERALS, get_mineral
from .solvers import Geometry, InitialCondition, dirichlet, zero_flux
from .solvers.history import ThermalHistory
from .solvers.initial import guess_step_from_data
from .thermo import log_fo2_from_delta

FORMAT_SETUP = "diffusor-setup"
FORMAT_BATCH = "diffusor-batch"
FILE_VERSION = 1
SUFFIX = ".diffusor"

# For these species the profile is the mineral's own composition variable (X_Fe, X_An,
# ...). For every other species it is a concentration (TiO2 wt%, Sr ppm) that says
# nothing about the host composition the law needs.
# species whose profile is the composition variable the law depends on (melt H2O: D grows with C_H2Ot)
EXCHANGE_SPECIES = ("Fe-Mg", "Fe-Ti", "NaSi-CaAl", "Na-K", "H2O")
# composition variables given in wt% rather than as a mole fraction
WT_PERCENT_VARIABLES = ("H2Ot",)
# measured olivine Fe-Mg values to X_Fe: (scale, offset) for XFe, XFo and Fo mol%
OLIVINE_COORDINATES = ((1., 0.), (-1., 1.), (-.01, 1.))

AXES = ("experiments", "a", "b", "c", "angles")
BOUNDARY_KEYS = ("far", "rim_melt", "rim_closed", "centre")
FO2_MODES = ("buffer", "absolute")
SOLVERS = ("auto", "numerical")

# Every setting with the value the window starts with. Monte Carlo cores are left out:
# they belong to the computer, not to the run.
DEFAULTS: Dict[str, Any] = {
    "mineral": next(iter(MINERALS)), "species": "", "coefficients": [],
    "axis": "experiments", "angles_deg": [90.0, 90.0, 0.0],
    "T_C": 950.0, "sigma_T_K": 20.0, "cooling": False, "T_end_C": 900.0,
    "P_MPa": 200.0, "sigma_P_MPa": 100.0,
    "fo2_mode": "buffer", "buffer": "NNO", "fo2_value": 1.0, "sigma_fo2": 0.3,
    "x_composition": 0.15, "x_composition_source": "default",
    "composition_dependent": True, "olivine_coordinate": 0,
    "geometry": "plane", "bc_left": "far", "bc_right": "far",
    "initial_condition": "step", "activity_set": DEFAULT_ACTIVITY_SET,
    "x0_um": 0.0, "C_left": 0.3, "C_right": 0.18, "smooth_um": 0.0,
    "fit_x0": True, "fit_plateaus": False,
    "solver": "auto", "n_nodes": 301,
    "resolution": "No correction", "beam_width_um": 5.0, "beam_sigma_um": 0.0,
    "distance_scale_sigma": 0.0,
    "mc_draws": 500, "mc_seed": 12345, "mc_temperature": True, "mc_fo2": True,
    "mc_pressure": True, "mc_coefficient": True, "mc_noise": True, "mc_distance": False,
    "mc_contributions": False, "coefficient_sampling": "auto",
}


def complete(settings: Dict[str, Any]) -> Dict[str, Any]:
    """``settings`` with every missing key filled from :data:`DEFAULTS`; unknown keys dropped."""
    out = {k: (list(v) if isinstance(v, list) else v) for k, v in DEFAULTS.items()}
    out.update({k: v for k, v in settings.items() if k in DEFAULTS})
    return out


# --- data --------------------------------------------------------------------------
@dataclass
class CrystalData:
    """A profile table and how it was read: what a setup needs besides the settings."""
    table: pd.DataFrame
    spec: ProfileSpec
    source: str = ""
    an_column: Optional[str] = None
    an_is_percent: bool = True
    dataset_key: Optional[str] = None
    excluded: Optional[List[bool]] = None      # one entry per point of the sorted profile

    @property
    def dataset(self) -> Optional[ds.ExampleDataset]:
        return ds.BY_KEY.get(self.dataset_key) if self.dataset_key else None

    def profile(self) -> Profile:
        p = build_profile(self.table, self.spec, source=self.source)
        if self.excluded is not None:
            mask = np.asarray(self.excluded, dtype=bool)
            if mask.shape == p.x.shape:
                p.set_excluded(mask)
        return p

    def an_values(self, profile: Optional[Profile] = None) -> Optional[np.ndarray]:
        if not self.an_column or self.an_column not in self.table.columns:
            return None
        vals = (profile or self.profile()).column(self.an_column)
        return vals / 100.0 if self.an_is_percent else vals

    @property
    def label(self) -> str:
        return Path(self.source).stem if self.source else "crystal"


def _plain(v):
    """A table value JSON can hold: numpy scalars as Python numbers, missing as None."""
    if isinstance(v, np.generic):
        v = v.item()
    if isinstance(v, float) and not math.isfinite(v):
        return None
    if v is None or isinstance(v, (bool, int, float, str)):
        return v
    try:
        if pd.isna(v):
            return None
    except (TypeError, ValueError):
        pass
    return str(v)


def table_to_json(df: pd.DataFrame) -> Dict:
    return {"columns": [str(c) for c in df.columns],
            "rows": [[_plain(v) for v in row] for row in df.itertuples(index=False, name=None)]}


def table_from_json(d: Dict) -> pd.DataFrame:
    df = pd.DataFrame(d["rows"], columns=d["columns"])
    for c in df.columns:
        col = df[c]
        num = pd.to_numeric(col, errors="coerce")
        if int(num.notna().sum()) == int(col.notna().sum()):
            df[c] = num.astype(float)
    return df


def data_to_json(data: CrystalData) -> Dict:
    return {"source": data.source, "spec": asdict(data.spec), "an_column": data.an_column,
            "an_is_percent": bool(data.an_is_percent), "dataset": data.dataset_key,
            "excluded": None if data.excluded is None else [bool(b) for b in data.excluded],
            "table": table_to_json(data.table)}


def data_from_json(d: Dict) -> CrystalData:
    names = {f.name for f in fields(ProfileSpec)}
    return CrystalData(table=table_from_json(d["table"]),
                       spec=ProfileSpec(**{k: v for k, v in d["spec"].items() if k in names}),
                       source=d.get("source", ""), an_column=d.get("an_column"),
                       an_is_percent=bool(d.get("an_is_percent", True)),
                       dataset_key=d.get("dataset"), excluded=d.get("excluded"))


# --- from settings to a model ------------------------------------------------------
def fo2_log_bar(s: Dict[str, Any]) -> float:
    T, P = s["T_C"] + 273.15, s["P_MPa"] * 1e6
    if s["fo2_mode"] == "buffer":
        return log_fo2_from_delta(s["buffer"], s["fo2_value"], T, P)
    return float(s["fo2_value"])


def conditions(s: Dict[str, Any], coef) -> Conditions:
    X = {k: s["x_composition"] for k in coef.requires}
    axis = s["axis"] if s["axis"] in ("a", "b", "c") else None
    angles = tuple(float(v) for v in s["angles_deg"]) if s["axis"] == "angles" else None
    return Conditions(T_K=s["T_C"] + 273.15, P_Pa=s["P_MPa"] * 1e6, log_fo2_bar=fo2_log_bar(s),
                      X=X, axis=axis, angles_deg=angles)


def table_initial_available(data: Optional[CrystalData]) -> bool:
    d = None if data is None else data.dataset
    return d is not None and bool(d.settings.get("initial_column"))


def initial_kind(s: Dict[str, Any], data: Optional[CrystalData], an_values=None) -> str:
    """The initial condition that applies: the chosen one if this crystal allows it, else a step."""
    kind = s["initial_condition"]
    if kind == "table" and table_initial_available(data):
        return "table"
    if kind == "equilibrium_plag" and an_values is not None:
        return "equilibrium_plag"
    return "step"


def coefficient_key(s: Dict[str, Any]) -> str:
    if not s["coefficients"]:
        raise ValueError("no diffusion coefficient is chosen")
    return s["coefficients"][0]


def build_model(s: Dict[str, Any], data: Optional[CrystalData], coef_key: Optional[str] = None,
                profile: Optional[Profile] = None, an_values=None) -> DiffusionModel:
    """The diffusion model for these settings and this crystal.

    ``profile`` and ``an_values`` may be passed when the caller has them already; the
    window does, so that a profile with points cut on the plot is used as it stands.
    """
    s = complete(s)
    coef = get_coefficient(coef_key or coefficient_key(s))
    cond = conditions(s, coef)
    mineral = get_mineral(s["mineral"])
    species = s["species"]
    if profile is None and data is not None:
        profile = data.profile()
    if an_values is None and data is not None and profile is not None:
        an_values = data.an_values(profile)
    kind = initial_kind(s, data, an_values)
    if kind == "table":
        values = profile.column(data.dataset.settings["initial_column"])
        if not np.all(np.isfinite(values)):
            raise ValueError("Published initial profile contains missing values in the fit window")
        ic = InitialCondition("table", {"x_table": profile.x.copy(), "C_table": values.copy()},
                              description="Published tabulated initial profile (linear interpolation)")
        c_left, c_right = float(values[0]), float(values[-1])
    elif kind == "equilibrium_plag":
        activity_set = s["activity_set"] or DEFAULT_ACTIVITY_SET
        ic = InitialCondition(
            "equilibrium_plag",
            {"x_an_x": np.asarray(profile.x, dtype=float),
             "x_an_values": np.asarray(an_values, dtype=float),
             "T_K": cond.T_K, "species": species, "activity_set": activity_set,
             "C_ref": float(s["C_left"]), "X_An_ref": float(an_values[0])},
            description=("equilibrium profile C0 exp(A X_An / RT) anchored at "
                         + format(s["C_left"], ".4g")
                         + " where X_An = " + format(float(an_values[0]), ".3f")
                         + " (Dohmen et al. 2017 App. eq. A13)"))
        c_left, c_right = float(s["C_left"]), float(s["C_right"])
    else:
        params = {"x0": float(s["x0_um"]), "C_left": float(s["C_left"]),
                  "C_right": float(s["C_right"])}
        if s["smooth_um"] > 0:
            params["smooth"] = float(s["smooth_um"])
        ic = InitialCondition("step", params)
        c_left, c_right = params["C_left"], params["C_right"]

    def boundary(kind_, plateau):
        return zero_flux() if kind_ in ("rim_closed", "centre") else dirichlet(plateau)

    hist = (ThermalHistory.linear(cond.T_K, s["T_end_C"] + 273.15, 1.0) if s["cooling"] else None)
    exchange = species in EXCHANGE_SPECIES
    comp_key = (mineral.composition_variable.key
                if s["composition_dependent"] and exchange else None)
    comp_scale, comp_offset = 1.0, 0.0
    if mineral.key == "olivine" and species == "Fe-Mg":
        comp_scale, comp_offset = OLIVINE_COORDINATES[int(s["olivine_coordinate"])]
        comp_key = "XFe" if s["composition_dependent"] else None
    if comp_key and comp_key not in coef.requires:
        comp_key = None

    A_kJ, an_grid = 0.0, None
    activity_set = s["activity_set"] or DEFAULT_ACTIVITY_SET
    if mineral.key == "plagioclase" and species in ACTIVITY_SPECIES and an_values is not None:
        A_kJ = activity_A(species, cond.T_K, activity_set)
    n_nodes = int(s["n_nodes"])
    x_grid = None
    if profile is not None:
        lo = float(profile.x.min())
        dataset = None if data is None else data.dataset
        if dataset is not None and "domain_start_um" in dataset.settings:
            # the crystal rim lies before the first analysis that was kept
            lo = min(lo, float(dataset.settings["domain_start_um"]))
        x_grid = np.linspace(lo, float(profile.x.max()), n_nodes)
        if A_kJ and an_values is not None:
            an_grid = np.interp(x_grid, profile.x, an_values)
    return DiffusionModel(
        coefficient=coef, conditions=cond, initial=ic,
        geometry=Geometry(s["geometry"]),
        bc_left=boundary(s["bc_left"], c_left), bc_right=boundary(s["bc_right"], c_right),
        history=hist, beam_sigma_um=float(s["beam_sigma_um"]), n_nodes=n_nodes, x_grid=x_grid,
        comp_key=comp_key, composition_dependent=bool(comp_key),
        comp_scale=comp_scale, comp_offset=comp_offset,
        an_profile=an_grid, activity_A_kJ=A_kJ if an_grid is not None else 0.0,
        activity_set=activity_set if an_grid is not None else "",
        force_numerical=s["solver"] == "numerical",
        boundaries_far=s["bc_left"] == "far" and s["bc_right"] == "far",
        fo2_buffer=((s["buffer"], float(s["fo2_value"])) if s["fo2_mode"] == "buffer" else None))


def free_parameters(s: Dict[str, Any], data: Optional[CrystalData] = None, an_values=None) -> List[str]:
    s = complete(s)
    free = ["t"]
    if initial_kind(s, data, an_values) != "step":
        return free
    if s["fit_x0"]:
        free.append("x0")
    if s["fit_plateaus"]:
        free += ["C_left", "C_right"]
    return free


def budget(s: Dict[str, Any]) -> UncertaintyBudget:
    s = complete(s)
    return UncertaintyBudget(
        sigma_T_K=s["sigma_T_K"] if s["mc_temperature"] else 0.0,
        fo2_mode=s["fo2_mode"], buffer=s["buffer"], delta_buffer=s["fo2_value"],
        sigma_delta_buffer=s["sigma_fo2"] if s["mc_fo2"] else 0.0,
        sigma_log_fo2=s["sigma_fo2"] if s["mc_fo2"] else 0.0,
        sigma_P_Pa=s["sigma_P_MPa"] * 1e6 if s["mc_pressure"] else 0.0,
        sample_coefficient=bool(s["mc_coefficient"]),
        coefficient_mode=s["coefficient_sampling"],
        sample_measurement_noise=bool(s["mc_noise"]),
        sigma_distance_scale=s["distance_scale_sigma"] if s["mc_distance"] else 0.0)


def profile_composition(s: Dict[str, Any], profile: Optional[Profile], an_values=None) -> Optional[np.ndarray]:
    """The profile in the law's composition coordinate, or None when it is not the host composition.

    Only for an exchange species is the profile the mineral's composition variable;
    olivine Fe-Mg goes through the chosen coordinate (X_Fe, X_Fo or Fo mol%). A TiO2
    wt% or ppm profile is a concentration, and its mean put into x_Ti or X_An would be
    nonsense. Values outside 0-1 (0-100 for melt water) are refused for the same reason.
    """
    s = complete(s)
    if profile is None or an_values is not None:
        return None
    mineral = get_mineral(s["mineral"])
    if s["species"] not in EXCHANGE_SPECIES or mineral.composition_variable.key == "none":
        return None
    scale, offset = 1.0, 0.0
    if mineral.key == "olivine" and s["species"] == "Fe-Mg":
        scale, offset = OLIVINE_COORDINATES[int(s["olivine_coordinate"])]
    X = offset + scale * np.asarray(profile.kept().C, dtype=float)
    top = 100.0 if mineral.composition_variable.key in WT_PERCENT_VARIABLES else 1.0
    if not (np.all(np.isfinite(X)) and 0.0 <= X.min() and X.max() <= top):
        return None
    return X


def guess_initial(s: Dict[str, Any], data: CrystalData, use_dataset: bool = True) -> Dict[str, float]:
    """Step position and plateaus from this crystal's own data (``x0_um``, ``C_left``, ``C_right``)."""
    s = complete(s)
    profile = data.profile()
    if initial_kind(s, data, data.an_values(profile)) == "table":
        values = profile.column(data.dataset.settings["initial_column"])
        return {"C_left": float(values[0]), "C_right": float(values[-1])}
    kept = profile.kept()
    start = dict(guess_step_from_data(kept.x, kept.C).params)
    if use_dataset and data.dataset is not None:
        start.update(data.dataset.settings.get("initial_step", {}))
    return {"x0_um": float(start["x0"]), "C_left": float(start["C_left"]),
            "C_right": float(start["C_right"])}


def crystal_settings(template: Dict[str, Any], data: CrystalData) -> Dict[str, Any]:
    """The template's settings with this crystal's own step, plateaus and composition."""
    s = complete(template)
    s.update(guess_initial(s, data))
    if s["x_composition_source"] == "profile":
        X = profile_composition(s, data.profile(), data.an_values())
        if X is not None:
            s["x_composition"] = float(np.mean(X))
    return s


# --- results -----------------------------------------------------------------------
def result_to_json(fit_result=None, mc_result=None) -> Optional[Dict]:
    if fit_result is None:
        return None
    r = {"coefficient": fit_result.model.coefficient.key,
         "t_seconds": float(fit_result.t_seconds),
         "free": {k: float(v) for k, v in fit_result.free.items()},
         "reduced_chi2": _plain(float(fit_result.stats.reduced_chi2)),
         "r_squared": _plain(float(fit_result.stats.r_squared)),
         "warnings": list(fit_result.warnings), "route": fit_result.route}
    if mc_result is not None:
        r["montecarlo"] = {"n_draws": int(mc_result.n_draws), "n_failed": int(mc_result.n_failed),
                           "seed": int(mc_result.seed),
                           "times_s": [float(t) for t in mc_result.times]}
    return r


# --- files -------------------------------------------------------------------------
def _stamp() -> Dict[str, str]:
    return {"diffusor_version": __version__,
            "saved_utc": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")}


def setup_to_json(settings: Dict[str, Any], data: CrystalData, result: Optional[Dict] = None,
                  label: str = "", group: str = "") -> Dict:
    return {"format": FORMAT_SETUP, "version": FILE_VERSION, **_stamp(),
            "label": label or data.label, "group": group,
            "settings": complete(settings), "data": data_to_json(data), "result": result}


def save_setup(path, settings: Dict[str, Any], data: CrystalData, result: Optional[Dict] = None,
               label: str = "", group: str = "") -> Path:
    path = Path(path)
    if path.suffix.lower() != SUFFIX:
        path = path.with_suffix(SUFFIX)
    path.write_text(json.dumps(setup_to_json(settings, data, result, label, group), indent=1),
                    encoding="utf-8")
    return path


def read_file(path) -> Dict:
    """A setup or batch file, checked for its format and version."""
    path = Path(path)
    try:
        d = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ValueError(f"{path.name} is not a Diffusor setup file: {exc}") from exc
    if not isinstance(d, dict) or d.get("format") not in (FORMAT_SETUP, FORMAT_BATCH):
        raise ValueError(f"{path.name} is not a Diffusor setup or batch file")
    if int(d.get("version", 0)) > FILE_VERSION:
        raise ValueError(f"{path.name} was written by a newer Diffusor "
                         f"({d.get('diffusor_version', '?')}); update Diffusor to open it")
    return d


@dataclass
class Crystal:
    """One crystal of a setup or batch file, with its effective settings."""
    label: str
    settings: Dict[str, Any]
    data: CrystalData
    group: str = ""
    result: Optional[Dict] = None
    overrides: Dict[str, Any] = field(default_factory=dict)


def crystals_in(d: Dict) -> List[Crystal]:
    """Every crystal in a parsed setup or batch file."""
    if d["format"] == FORMAT_SETUP:
        data = data_from_json(d["data"])
        return [Crystal(d.get("label") or data.label, complete(d["settings"]), data,
                        d.get("group", ""), d.get("result"))]
    shared = complete(d.get("shared", {}))
    out = []
    for c in d.get("crystals", []):
        data = data_from_json(c["data"])
        over = {k: v for k, v in (c.get("overrides") or {}).items() if k in DEFAULTS}
        out.append(Crystal(c.get("label") or data.label, complete({**shared, **over}), data,
                           c.get("group", ""), c.get("result"), over))
    return out


def load_setup(path) -> Crystal:
    """The one crystal of a setup file."""
    d = read_file(path)
    if d["format"] != FORMAT_SETUP:
        raise ValueError(f"{Path(path).name} is a batch of crystals; open it in the batch study")
    return crystals_in(d)[0]


# --- running one crystal -----------------------------------------------------------
def run_crystal(settings: Dict[str, Any], data: CrystalData, monte_carlo: bool = False,
                workers: int = 1, progress=None):
    """Fit one crystal, and run its Monte Carlo if asked, as the main window would.

    Returns ``(fit_result, mc_result)``; ``mc_result`` is None without a Monte Carlo.
    ``progress(i, n)`` may return True to stop the Monte Carlo.
    """
    from .fitting import fit_time, run_montecarlo
    from .fitting.fit import T_MAX_DEFAULT, T_MIN_DEFAULT
    s = complete(settings)
    profile = data.profile()
    an = data.an_values(profile)
    model = build_model(s, data, profile=profile, an_values=an)
    free = free_parameters(s, data, an)
    kept = profile.kept()
    if not monte_carlo:
        return fit_time(model, kept.x, kept.C, kept.sigma, free, T_MIN_DEFAULT, T_MAX_DEFAULT), None
    mc = run_montecarlo(model, kept.x, kept.C, kept.sigma, budget=budget(s),
                        n_draws=int(s["mc_draws"]), seed=int(s["mc_seed"]), free_parameters=free,
                        t_min=T_MIN_DEFAULT, t_max=T_MAX_DEFAULT, progress=progress,
                        workers=workers)
    fit = mc.base_fit if mc.base_fit is not None else fit_time(
        model, kept.x, kept.C, kept.sigma, free, T_MIN_DEFAULT, T_MAX_DEFAULT)
    return fit, mc


# --- a batch of crystals sharing settings ---------------------------------------------
def _same(a, b) -> bool:
    if isinstance(a, float) or isinstance(b, float):
        try:
            return math.isclose(float(a), float(b), rel_tol=1e-12, abs_tol=1e-12)
        except (TypeError, ValueError):
            return False
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        return len(a) == len(b) and all(_same(x, y) for x, y in zip(a, b))
    return a == b


def differences(settings: Dict[str, Any], shared: Dict[str, Any]) -> Dict[str, Any]:
    """The settings that differ from ``shared``: what a crystal must keep as its own."""
    s, base = complete(settings), complete(shared)
    return {k: v for k, v in s.items() if not _same(v, base[k])}


@dataclass
class BatchCrystal:
    label: str
    data: CrystalData
    group: str = ""
    overrides: Dict[str, Any] = field(default_factory=dict)
    result: Optional[Dict] = None
    status: str = ""


class Batch:
    """Crystals that share settings, each keeping the values it sets for itself.

    A crystal's settings are the shared settings with its own values on top. Setting a
    shared value sets it for every crystal: by default the crystals' own values for that
    setting are cleared, so the new value really applies to all.
    """

    def __init__(self, shared: Optional[Dict[str, Any]] = None):
        self.shared = complete(shared or {})
        self.crystals: List[BatchCrystal] = []

    def __len__(self) -> int:
        return len(self.crystals)

    def settings(self, i: int) -> Dict[str, Any]:
        return complete({**self.shared, **self.crystals[i].overrides})

    def value(self, i: Optional[int], key: str):
        return self.shared[key] if i is None else self.settings(i)[key]

    def is_own(self, i: int, key: str) -> bool:
        return key in self.crystals[i].overrides

    def own_count(self, key: str) -> int:
        return sum(key in c.overrides for c in self.crystals)

    def set_shared(self, key: str, value, keep_own: bool = False) -> None:
        if key not in DEFAULTS:
            raise KeyError(key)
        self.shared[key] = value
        for c in self.crystals:
            if not keep_own:
                c.overrides.pop(key, None)
            elif key in c.overrides and _same(c.overrides[key], value):
                del c.overrides[key]
        self._invalidate(range(len(self.crystals)))

    def set_own(self, i: int, key: str, value) -> None:
        if key not in DEFAULTS:
            raise KeyError(key)
        c = self.crystals[i]
        if _same(value, self.shared[key]):
            c.overrides.pop(key, None)
        else:
            c.overrides[key] = value
        self._invalidate([i])

    def clear_own(self, i: int, key: str) -> None:
        self.crystals[i].overrides.pop(key, None)
        self._invalidate([i])

    def _invalidate(self, rows) -> None:
        for i in rows:
            c = self.crystals[i]
            if c.result is not None:
                c.result = None
                c.status = "settings changed: fit again"

    def add(self, label: str, data: CrystalData, settings: Optional[Dict[str, Any]] = None,
            group: str = "", result: Optional[Dict] = None) -> BatchCrystal:
        """Add a crystal. Without ``settings``, its step and plateaus come from its own data."""
        s = crystal_settings(self.shared, data) if settings is None else complete(settings)
        taken = {c.label for c in self.crystals}
        base, k = label, 2
        while label in taken:
            label = f"{base} ({k})"
            k += 1
        c = BatchCrystal(label, data, group, differences(s, self.shared), result,
                         "saved result" if result else "")
        self.crystals.append(c)
        return c

    def to_json(self) -> Dict:
        return {"format": FORMAT_BATCH, "version": FILE_VERSION, **_stamp(),
                "shared": complete(self.shared),
                "crystals": [{"label": c.label, "group": c.group, "overrides": c.overrides,
                              "data": data_to_json(c.data), "result": c.result}
                             for c in self.crystals]}

    def save(self, path) -> Path:
        path = Path(path)
        if path.suffix.lower() != SUFFIX:
            path = path.with_suffix(SUFFIX)
        path.write_text(json.dumps(self.to_json(), indent=1), encoding="utf-8")
        return path

    @classmethod
    def from_json(cls, d: Dict) -> "Batch":
        if d["format"] != FORMAT_BATCH:
            raise ValueError("not a batch file")
        b = cls(d.get("shared", {}))
        for c in crystals_in(d):
            b.crystals.append(BatchCrystal(c.label, c.data, c.group, dict(c.overrides), c.result,
                                           "saved result" if c.result else ""))
        return b
