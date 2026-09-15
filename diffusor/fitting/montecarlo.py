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
* The Arrhenius parameters ln D0 and Q are strongly anti-correlated by the
  regression that produced them.  Sampling them independently inflates the
  spread of D at the temperature of interest enormously.  Diffusor samples
  from the published covariance when it exists, and otherwise samples ln D
  directly at the working temperature using the paper's stated scatter (see
  :mod:`diffusor.coefficients.base`).
* The composition enters both the diffusion coefficient and the profile being
  fitted, so measurement noise propagates through two routes at once.

The linear-propagation formula used by NIDIS (Petrone et al. 2016, Methods)
treats sqrt(4Dt) and T as independent, which neither captures the
buffer-temperature correlation nor the D0-Q correlation.  Each draw here
re-runs the *whole* fit, so every correlation is honoured automatically.

Reported statistics
-------------------
Diffusion times are approximately log-normal (as Mutch et al. 2021 also find),
so the median and the 16th/84th percentiles are reported rather than a mean
and a symmetric sigma.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence

import numpy as np

from ..thermo.buffers import log_fo2_from_delta
from ..thermo.units import human_time
from .fit import fit_time
from .model import DiffusionModel

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
    """One sampled Conditions object, with T-fO2 correlation built in."""
    cond = model.conditions
    active = budget.active_sources() if only is None else ([only] if only in budget.active_sources() else [])

    T = cond.T_K
    if "temperature" in active:
        T = float(rng.normal(cond.T_K, budget.sigma_T_K))
    P = cond.P_Pa
    if "pressure" in active:
        P = max(0.0, float(rng.normal(cond.P_Pa, budget.sigma_P_Pa)))

    log_fo2 = cond.log_fo2_bar
    if budget.fo2_mode == "buffer":
        delta = budget.delta_buffer
        if "fo2" in active:
            delta = float(rng.normal(budget.delta_buffer, budget.sigma_delta_buffer))
        # the buffer is re-evaluated at the *sampled* temperature: this is the
        # T-fO2 correlation that independent sampling misses
        log_fo2 = log_fo2_from_delta(budget.buffer, delta, T, P)
    elif "fo2" in active and log_fo2 is not None:
        log_fo2 = float(rng.normal(cond.log_fo2_bar, budget.sigma_log_fo2))

    return cond.replace(T_K=T, P_Pa=P, log_fo2_bar=log_fo2)


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


def run(model: DiffusionModel, x_data, C_data, sigma=None, *,
        budget: Optional[UncertaintyBudget] = None,
        n_draws: int = 1000, seed: int = 12345,
        free_parameters: Sequence[str] = ("t",),
        t_min: float = 1.0e2, t_max: float = 3.2e12,
        keep_profiles: int = 200,
        progress: Optional[Callable[[int, int], bool]] = None,
        only_source: Optional[str] = None) -> MonteCarloResult:
    """Run the Monte Carlo, re-fitting the time for every draw.

    ``progress(i, n)`` may return True to abort.  ``only_source`` restricts the
    sampling to a single source, which is how the variance contributions are
    computed.
    """
    budget = budget or UncertaintyBudget()
    rng = np.random.default_rng(seed)
    x_data = np.asarray(x_data, dtype=float)
    C_data = np.asarray(C_data, dtype=float)

    base = fit_time(model, x_data, C_data, sigma, free_parameters, t_min, t_max)
    times: List[float] = []
    kept: List[np.ndarray] = []
    n_failed = 0
    active = budget.active_sources()

    for i in range(n_draws):
        cond = _draw_conditions(model, budget, rng, only_source)
        xd, Cd = _draw_data(x_data, C_data, sigma, budget, rng, only_source)
        overrides = {}
        if ("diffusion_coefficient" in active
                and (only_source is None or only_source == "diffusion_coefficient")):
            try:
                overrides = model.coefficient.sample(rng, budget.coefficient_mode)
            except ValueError:
                overrides = {}
        import copy
        m = copy.copy(model)
        m.conditions = cond
        if "boundary_compositions" in active and budget.sigma_boundary > 0 \
                and (only_source is None or only_source == "boundary_compositions"):
            m.initial = copy.deepcopy(model.initial)
            for k in ("C_left", "C_right", "C_core", "C_rim"):
                if k in m.initial.params:
                    m.initial.params[k] = float(rng.normal(m.initial.params[k], budget.sigma_boundary))
        try:
            if budget.refit_each_draw:
                r = fit_time(m, xd, Cd, sigma, free_parameters, t_min, t_max,
                             overrides=overrides, t_guess=base.t_seconds)
                t = r.t_seconds
                prof = r.C_model
            else:
                t = base.t_seconds
                prof = m.profile(t, x_data, overrides)
            if not np.isfinite(t) or t <= 0:
                raise ValueError("non-finite time")
            times.append(float(t))
            if len(kept) < keep_profiles:
                kept.append(np.asarray(prof, dtype=float))
        except Exception:
            n_failed += 1
        if progress is not None and (i % 5 == 0):
            if progress(i + 1, n_draws):
                break

    if not times:
        raise RuntimeError("every Monte Carlo draw failed; check the model set-up")

    res = MonteCarloResult(
        times=np.array(times), t_best=base.t_seconds, n_draws=len(times),
        n_failed=n_failed, seed=seed, budget=budget,
        profiles=np.array(kept) if kept else None,
        x_profiles=x_data, warnings=list(base.warnings))
    if n_failed > 0.1 * n_draws:
        res.warnings.append(f"{n_failed} of {n_draws} draws failed; the result may be biased")
    return res


def contributions(model: DiffusionModel, x_data, C_data, sigma=None, *,
                  budget: Optional[UncertaintyBudget] = None,
                  n_draws: int = 200, seed: int = 12345,
                  free_parameters: Sequence[str] = ("t",),
                  progress: Optional[Callable[[int, int], bool]] = None) -> Dict[str, float]:
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
                progress=progress)
        out[src] = r.sigma_log10
    return out
