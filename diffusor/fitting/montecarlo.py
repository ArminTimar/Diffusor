"""Monte Carlo error propagation.

Why Monte Carlo rather than analytical propagation
--------------------------------------------------
The time retrieved from a diffusion profile depends on the intensive
parameters through a strongly non-linear function, and those parameters are
*not* independent:

* fO2 is usually known as an offset from a mineral buffer (NNO+1, FMQ-0.5).
  The buffer itself is a function of temperature, so sampling T and fO2
  independently is wrong -- their correlation is fixed by the buffer equation.
  Diffusor re-evaluates the buffer at each sampled temperature, which builds
  the correlation in exactly.
* The Arrhenius parameters ln D0 and Q are typically positively correlated by the
  regression that produced them.  Sampling them independently inflates the
  spread of D at the temperature of interest enormously.  Diffusor samples
  from a covariance matrix when the registry holds one for the coefficient
  (either published, or built by Diffusor from a correlation the source states
  in words; the entry says which), and otherwise samples ln D directly at the
  working temperature using the paper's stated scatter (see
  :mod:`diffusor.coefficients.base`).
* The composition enters both the diffusion coefficient and the profile being
  fitted, so measurement noise propagates through two routes at once.

The linear-propagation formula used by NIDIS (Petrone et al. 2016, Methods)
treats sqrt(4Dt) and T as independent, which neither captures the
buffer-temperature correlation nor the D0-Q correlation.  Each draw here
re-runs the *whole* fit. Only correlations explicitly encoded by the input
model (such as buffer-temperature coupling or a supplied covariance) are
propagated; unmodelled systematic errors remain unquantified.

Reported statistics
-------------------
The empirical median and percentiles describe the sampled times without
assuming a distribution shape. Skewness and multiple modes can matter; these
are propagated-input ensembles, not Bayesian posterior samples.
"""
from __future__ import annotations

import copy
import multiprocessing as mp
import os
import pickle
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures.process import BrokenProcessPool
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence

import numpy as np

from ..solvers.history import ThermalHistory
from ..thermo.buffers import log_fo2_from_delta
from ..thermo.units import human_time
from .fit import T_MAX_DEFAULT, T_MIN_DEFAULT, fit_time
from .model import DiffusionModel

# expected run time below which the draws stay in this process (starting workers
# costs about a second each)
PARALLEL_THRESHOLD_S = 6.0
BLAS_THREAD_VARIABLES = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")

SOURCES = ("temperature", "fo2", "pressure", "diffusion_coefficient",
           "measurement_noise", "distance_scale", "boundary_compositions")


@dataclass
class UncertaintyBudget:
    """1-sigma uncertainties on everything that feeds the timescale."""
    sigma_T_K: float = 0.0
    fo2_mode: str = "buffer"              # 'buffer' (delta from a buffer) or 'absolute'
    buffer: str = "NNO"
    delta_buffer: float = 0.0
    sigma_delta_buffer: float = 0.0
    sigma_log_fo2: float = 0.0            # used when fo2_mode == 'absolute'
    sigma_P_Pa: float = 0.0
    sample_coefficient: bool = True
    coefficient_mode: str = "auto"        # 'auto' | 'covariance' | 'logD_at_T' | 'independent'
    sample_measurement_noise: bool = True
    sigma_distance_scale: float = 0.0     # relative, e.g. 0.02 for a 2% image calibration
    sigma_boundary: float = 0.0           # on the plateau compositions
    refit_each_draw: bool = True

    def active_sources(self) -> List[str]:
        s = []
        if self.sigma_T_K > 0:
            s.append("temperature")
        if (self.fo2_mode == "buffer" and self.sigma_delta_buffer > 0) or \
           (self.fo2_mode == "absolute" and self.sigma_log_fo2 > 0):
            s.append("fo2")
        if self.sigma_P_Pa > 0:
            s.append("pressure")
        if self.sample_coefficient:
            s.append("diffusion_coefficient")
        if self.sample_measurement_noise:
            s.append("measurement_noise")
        if self.sigma_distance_scale > 0:
            s.append("distance_scale")
        if self.sigma_boundary > 0:
            s.append("boundary_compositions")
        return s


@dataclass
class MonteCarloResult:
    times: np.ndarray
    t_best: float
    n_draws: int
    n_failed: int
    seed: int
    budget: UncertaintyBudget
    profiles: Optional[np.ndarray] = None       # (n_kept, n_x) model profiles
    x_profiles: Optional[np.ndarray] = None
    contributions: Dict[str, float] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)
    workers: int = 1                            # processes that fitted the draws

    # -- statistics ---------------------------------------------------------
    @property
    def median(self) -> float:
        return float(np.median(self.times))

    def percentile(self, q) -> float:
        return float(np.percentile(self.times, q))

    @property
    def p16(self) -> float:
        return self.percentile(16)

    @property
    def p84(self) -> float:
        return self.percentile(84)

    @property
    def p2_5(self) -> float:
        return self.percentile(2.5)

    @property
    def p97_5(self) -> float:
        return self.percentile(97.5)

    @property
    def sigma_log10(self) -> float:
        return float(np.std(np.log10(self.times)))

    def envelope(self, q_low: float = 16, q_high: float = 84):
        """Percentile envelope of the modelled profiles."""
        if self.profiles is None:
            return None
        return (np.percentile(self.profiles, q_low, axis=0),
                np.percentile(self.profiles, q_high, axis=0))

    def summary(self) -> str:
        lines = [
            f"best fit          : {human_time(self.t_best)}",
            f"Monte Carlo median: {human_time(self.median)}  ({self.n_draws} draws, seed {self.seed})",
            f"  68% interval    : {human_time(self.p16)} to {human_time(self.p84)}",
            f"  95% interval    : {human_time(self.p2_5)} to {human_time(self.p97_5)}",
            f"  scatter         : {self.sigma_log10:.3f} log10 units (1 sigma)",
        ]
        if self.n_failed:
            lines.append(f"  {self.n_failed} draws failed and were discarded")
        if self.contributions:
            lines.append("  variance contributions (sigma of log10 t with only that source active):")
            for k, v in sorted(self.contributions.items(), key=lambda kv: -kv[1]):
                lines.append(f"     {k:<24s} {v:.3f}")
        for w in self.warnings:
            lines.append(f"  ! {w}")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
def _draw_conditions(model: DiffusionModel, budget: UncertaintyBudget,
                     rng: np.random.Generator, only: Optional[str] = None):
    """One sampled Conditions object, with T-fO2 correlation built in.

    Returns the conditions and the sampled buffer offset (None when fO2 is given
    as an absolute value), which a cooling path needs to follow the buffer down.
    """
    cond = model.conditions
    active = budget.active_sources() if only is None else ([only] if only in budget.active_sources() else [])

    T = cond.T_K
    if "temperature" in active:
        T = float(rng.normal(cond.T_K, budget.sigma_T_K))
    P = cond.P_Pa
    if "pressure" in active:
        P = max(0.0, float(rng.normal(cond.P_Pa, budget.sigma_P_Pa)))

    log_fo2 = cond.log_fo2_bar
    delta = None
    if budget.fo2_mode == "buffer":
        delta = budget.delta_buffer
        if "fo2" in active:
            delta = float(rng.normal(budget.delta_buffer, budget.sigma_delta_buffer))
        # the buffer is re-evaluated at the *sampled* temperature: this is the
        # T-fO2 correlation that independent sampling misses
        log_fo2 = log_fo2_from_delta(budget.buffer, delta, T, P)
    elif "fo2" in active and log_fo2 is not None:
        log_fo2 = float(rng.normal(cond.log_fo2_bar, budget.sigma_log_fo2))

    return cond.replace(T_K=T, P_Pa=P, log_fo2_bar=log_fo2), delta


def _draw_data(x_data, C_data, sigma, budget: UncertaintyBudget,
               rng: np.random.Generator, only: Optional[str] = None):
    active = budget.active_sources() if only is None else ([only] if only in budget.active_sources() else [])
    x = np.asarray(x_data, dtype=float)
    C = np.asarray(C_data, dtype=float)
    if "measurement_noise" in active and sigma is not None:
        C = C + rng.normal(0.0, np.asarray(sigma, dtype=float))
    if "distance_scale" in active:
        x = x * float(rng.normal(1.0, budget.sigma_distance_scale))
    return x, C


@dataclass
class _Draw:
    """Everything random about one draw.

    Draws are sampled in the calling process, in draw order, so the seed alone
    fixes every result however many processes then fit them.
    """
    index: int
    conditions: object
    fo2_buffer: Optional[tuple]
    history: Optional[ThermalHistory]
    initial: object                 # None keeps the model's initial condition
    overrides: dict
    x: np.ndarray
    C: np.ndarray


def _sample_draw(i: int, model: DiffusionModel, budget: UncertaintyBudget,
                 rng: np.random.Generator, x_data, C_data, sigma,
                 only_source: Optional[str], active: Sequence[str]) -> _Draw:
    cond, delta = _draw_conditions(model, budget, rng, only_source)
    xd, Cd = _draw_data(x_data, C_data, sigma, budget, rng, only_source)
    overrides = {}
    if ("diffusion_coefficient" in active
            and (only_source is None or only_source == "diffusion_coefficient")):
        try:
            overrides = model.coefficient.sample(rng, budget.coefficient_mode)
        except ValueError:
            overrides = {}
    initial = None
    if "boundary_compositions" in active and budget.sigma_boundary > 0 \
            and (only_source is None or only_source == "boundary_compositions"):
        initial = copy.deepcopy(model.initial)
        for k in ("C_left", "C_right", "C_core", "C_rim"):
            if k in initial.params:
                initial.params[k] = float(rng.normal(initial.params[k], budget.sigma_boundary))
    fo2_buffer = model.fo2_buffer
    if delta is not None and budget.fo2_mode == "buffer":
        fo2_buffer = (budget.buffer, delta)
    history = model.history
    if history is not None and cond.T_K != model.conditions.T_K:
        # a cooling path moves as a whole with the sampled temperature: the
        # uncertainty is on where the path sits, its shape is kept
        dT = cond.T_K - model.conditions.T_K
        history = ThermalHistory(history.times.copy(), history.temps_K + dT, history.label)
    return _Draw(i, cond, fo2_buffer, history, initial, overrides, xd, Cd)


def _evaluate(model: DiffusionModel, d: _Draw, sigma, free_parameters, t_min, t_max,
              t_guess: float, refit: bool, x_data) -> Optional[dict]:
    """Fit one draw. Returns what ``on_draw`` receives, or None when the draw fails."""
    try:
        m = copy.copy(model)
        m.conditions, m.fo2_buffer, m.history = d.conditions, d.fo2_buffer, d.history
        if d.initial is not None:
            m.initial = d.initial
        if refit:
            r = fit_time(m, d.x, d.C, sigma, free_parameters, t_min, t_max,
                         overrides=d.overrides, t_guess=t_guess)
            t, prof = r.t_seconds, r.C_model
        else:
            t = t_guess
            prof = m.profile(t, x_data, d.overrides)
        if not np.isfinite(t) or t <= 0:
            return None
        try:
            logD = float(np.log10(m.D_bulk(d.overrides, C_ref=float(np.mean(d.C)))))
        except Exception:
            logD = float("nan")
        return {"index": d.index, "t": float(t), "T_K": float(d.conditions.T_K),
                "log_fo2_bar": (None if d.conditions.log_fo2_bar is None
                                else float(d.conditions.log_fo2_bar)),
                "log10_D": logD, "x": np.asarray(d.x, dtype=float),
                "C": np.asarray(d.C, dtype=float), "C_model": np.asarray(prof, dtype=float)}
    except Exception:
        return None


# --- parallel evaluation --------------------------------------------------------
# Worker processes receive the base model once, when they start, and then batches
# of draws. The coefficient travels as its registry key: many laws are closures
# (the magnetite tables, for example), which cannot be pickled, and every worker
# rebuilds the same registry on import anyway.
_WORKER: dict = {}


def _worker_init(payload: bytes, coefficient_key: Optional[str]) -> None:
    context = pickle.loads(payload)
    if coefficient_key is not None:
        from ..coefficients import get
        context["model"].coefficient = get(coefficient_key)
    _WORKER.clear()
    _WORKER.update(context)


def _worker_batch(draws: List[_Draw]) -> List[Optional[dict]]:
    w = _WORKER
    return [_evaluate(w["model"], d, w["sigma"], w["free"], w["t_min"], w["t_max"],
                      w["t_guess"], w["refit"], w["x_data"]) for d in draws]


def _portable(model: DiffusionModel):
    """The model with its coefficient replaced by a registry key, if it is one."""
    from ..coefficients import get
    try:
        if get(model.coefficient.key) is model.coefficient:
            light = copy.copy(model)
            light.coefficient = None
            return light, model.coefficient.key
    except KeyError:
        pass
    pickle.dumps(model)            # raises when the model cannot go to another process
    return model, None


def default_workers() -> int:
    """All cores but one, so the computer stays usable while a run goes on."""
    return max(1, (os.cpu_count() or 2) - 1)


def run(model: DiffusionModel, x_data, C_data, sigma=None, *,
        budget: Optional[UncertaintyBudget] = None,
        n_draws: int = 1000, seed: int = 12345,
        free_parameters: Sequence[str] = ("t",),
        t_min: float = T_MIN_DEFAULT, t_max: float = T_MAX_DEFAULT,
        keep_profiles: int = 200,
        progress: Optional[Callable[[int, int], bool]] = None,
        only_source: Optional[str] = None,
        on_draw: Optional[Callable[[dict], None]] = None,
        workers: Optional[int] = 1,
        parallel_threshold_s: float = PARALLEL_THRESHOLD_S) -> MonteCarloResult:
    """Run the Monte Carlo, re-fitting the time for every draw.

    ``progress(i, n)`` may return True to abort.  ``only_source`` restricts the
    sampling to a single source, which is how the variance contributions are
    computed.  ``on_draw(info)`` receives every successful draw: the fitted time,
    the sampled temperature, fO2 and log10 D, the perturbed data and the fitted
    profile. The interface uses it to draw the Monte Carlo while it runs.

    ``workers`` above 1 fits the draws in that many processes (None or 0: all
    cores but one). A run expected to take less than ``parallel_threshold_s``
    stays in this process, because starting workers costs about a second each.
    The draws are sampled here either way, so the result for a given seed does
    not depend on the number of workers.
    """
    budget = budget or UncertaintyBudget()
    if (model.coefficient.fixed_temperature_K is not None and budget.sigma_T_K > 0
            and only_source in (None, "temperature")):
        raise ValueError("This coefficient has measurements at one temperature only. "
                         "Set temperature uncertainty to zero or select a calibrated temperature-dependent law.")
    rng = np.random.default_rng(seed)
    x_data = np.asarray(x_data, dtype=float)
    C_data = np.asarray(C_data, dtype=float)

    started = time.perf_counter()
    base = fit_time(model, x_data, C_data, sigma, free_parameters, t_min, t_max)
    # a draw is seeded from the best fit and scans a narrower range: about 2.5x faster
    per_draw = (time.perf_counter() - started) / 2.5
    active = budget.active_sources()
    draws = [_sample_draw(i, model, budget, rng, x_data, C_data, sigma, only_source, active)
             for i in range(n_draws)]

    n_workers = default_workers() if not workers else int(workers)
    n_workers = min(n_workers, n_draws, 61)          # 61: the Windows process limit
    if per_draw * n_draws < parallel_threshold_s:
        n_workers = 1

    results: Dict[int, dict] = {}
    notes: List[str] = []
    state = {"done": 0, "aborted": False}

    def collect(r: Optional[dict]) -> None:
        state["done"] += 1
        if r is None:
            return
        results[r["index"]] = r
        if on_draw is not None:
            on_draw({k: v for k, v in r.items() if k != "index"})

    def run_here(pending: List[_Draw]) -> None:
        for d in pending:
            collect(_evaluate(model, d, sigma, free_parameters, t_min, t_max,
                              base.t_seconds, budget.refit_each_draw, x_data))
            if progress is not None and (state["done"] % 5 == 0 or state["done"] == n_draws):
                if progress(state["done"], n_draws):
                    state["aborted"] = True
                    return

    if n_workers > 1:
        try:
            light, key = _portable(model)
            payload = pickle.dumps(dict(model=light, sigma=sigma, free=tuple(free_parameters),
                                        t_min=t_min, t_max=t_max, t_guess=base.t_seconds,
                                        refit=budget.refit_each_draw, x_data=x_data))
        except Exception as exc:
            notes.append(f"the draws ran on one core: the model cannot be sent to other "
                         f"processes ({type(exc).__name__})")
            n_workers = 1
    if n_workers > 1:
        # batches of about half a second each, and several per worker for balance
        size = max(1, min(int(0.5 / max(per_draw, 1e-3)),
                          int(np.ceil(n_draws / (4 * n_workers)))))
        batches = [draws[i:i + size] for i in range(0, n_draws, size)]
        pool = ProcessPoolExecutor(max_workers=n_workers, mp_context=mp.get_context("spawn"),
                                   initializer=_worker_init, initargs=(payload, key))
        finished = set()
        # one BLAS thread per worker: the workers are the parallelism, and each
        # starting its own thread pool would oversubscribe the cores. Workers
        # start when the batches are submitted and copy the environment then.
        saved = {v: os.environ.get(v) for v in BLAS_THREAD_VARIABLES}
        os.environ.update({v: "1" for v in BLAS_THREAD_VARIABLES})
        try:
            try:
                futures = {pool.submit(_worker_batch, b): n for n, b in enumerate(batches)}
            finally:
                for v, old in saved.items():
                    if old is None:
                        os.environ.pop(v, None)
                    else:
                        os.environ[v] = old
            for f in as_completed(futures):
                for r in f.result():
                    collect(r)
                finished.add(futures[f])
                if progress is not None and progress(state["done"], n_draws):
                    state["aborted"] = True
                    break
        except BrokenProcessPool:
            left = [d for n, b in enumerate(batches) if n not in finished for d in b]
            notes.append(f"the worker processes stopped; {len(left)} draws were finished "
                         "on one core")
            run_here(left)
        finally:
            pool.shutdown(wait=not state["aborted"], cancel_futures=True)
    else:
        run_here(draws)

    kept = [results[i] for i in sorted(results)]
    if not kept:
        raise RuntimeError("every Monte Carlo draw failed. Check the model set-up.")
    n_run = state["done"] if state["aborted"] else n_draws
    n_failed = max(0, n_run - len(kept))
    times = np.array([r["t"] for r in kept])
    profiles = [r["C_model"] for r in kept[:keep_profiles]]
    res = MonteCarloResult(
        times=times, t_best=base.t_seconds, n_draws=len(kept),
        n_failed=n_failed, seed=seed, budget=budget,
        profiles=np.array(profiles) if profiles else None,
        x_profiles=x_data, warnings=list(base.warnings) + notes, workers=n_workers)
    if n_failed > 0.1 * n_draws:
        res.warnings.append(f"{n_failed} of {n_draws} draws failed. The result may be biased.")
    at_bound = int(np.sum(res.times >= 0.95 * t_max) + np.sum(res.times <= 1.05 * t_min))
    if at_bound:
        res.warnings.append(
            f"{at_bound} of {len(res.times)} draws stopped at the edge of the search range "
            f"({human_time(t_min)} to {human_time(t_max)}), so the interval is cut off there. "
            "Those draws are the coldest ones. Widen the range or tighten the temperature.")
    return res


def contributions(model: DiffusionModel, x_data, C_data, sigma=None, *,
                  budget: Optional[UncertaintyBudget] = None,
                  n_draws: int = 200, seed: int = 12345,
                  free_parameters: Sequence[str] = ("t",),
                  progress: Optional[Callable[[int, int], bool]] = None,
                  workers: Optional[int] = 1) -> Dict[str, float]:
    """Variance decomposition: sigma(log10 t) with one source active at a time.

    This is a one-at-a-time sensitivity analysis, so the individual values do
    not add up to the total in quadrature when the sources interact; it is
    meant for ranking which measurement to improve first.
    """
    budget = budget or UncertaintyBudget()
    out: Dict[str, float] = {}
    for src in budget.active_sources():
        r = run(model, x_data, C_data, sigma, budget=budget, n_draws=n_draws, seed=seed,
                free_parameters=free_parameters, keep_profiles=0, only_source=src,
                progress=progress, workers=workers)
        out[src] = r.sigma_log10
    return out
