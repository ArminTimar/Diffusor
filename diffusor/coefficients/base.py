"""Diffusion coefficient framework.

A :class:`DiffusionCoefficient` is a *traceable* implementation of one
published diffusion law.  It carries

* the callable that evaluates D (m^2/s) from :class:`Conditions`,
* the named fit parameters with their published uncertainties and, where the
  source gives it, their covariance,
* the calibration ranges (T, P, fO2, composition) so the GUI can warn when a
  run extrapolates,
* the citation key, the equation number and a verbatim transcription of the
  equation as printed in the source,
* a ``verified`` flag stating whether the numbers were read from the primary
  publication or taken from a secondary compilation.

Uncertainty and Monte Carlo
---------------------------
Diffusion laws are fitted as ``ln D = ln D0 - Q/RT (+ n ln fO2 + m X + ...)``.
The fit parameters are *strongly anti-correlated*: a higher ln D0 always comes
with a higher Q.  Sampling them independently (as NIDIS does, Petrone et al.
2016) inflates the uncertainty at the temperature of interest by a large
factor.  Diffusor therefore offers three sampling modes, in order of
preference:

``covariance``
    Sample the parameter vector from the published covariance matrix
    (this is the approach of Mutch et al. 2021, DFENS).
``logD_at_T``
    Sample ln D directly at the working temperature using the paper's stated
    scatter about the Arrhenius line (e.g. "reproduces the data within 1 log
    unit", Mueller et al. 2013).  This is the honest fallback when no
    covariance is published, and is what most papers actually constrain.
``independent``
    Sample each parameter independently.  Provided only so the user can
    reproduce the (over-)estimates of other software; never the default.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..constants import R_GAS
from ..references import cite, get as get_reference

LN10 = np.log(10.0)


# ---------------------------------------------------------------------------
@dataclass
class Conditions:
    """Intensive conditions at which D is evaluated.

    ``X`` holds the composition variables the law needs, e.g.
    ``{'XFe': 0.12}`` or ``{'XFo': array_of_nodes}``.  Array values are
    allowed so the numerical solver can evaluate D node by node.
    """
    T_K: float
    P_Pa: float = 1.0e5
    log_fo2_bar: Optional[float] = None
    X: Dict[str, object] = field(default_factory=dict)
    axis: Optional[str] = None                      # 'a' | 'b' | 'c' | None
    angles_deg: Optional[Tuple[float, float, float]] = None   # to a, b, c

    def x(self, key: str, default=None):
        v = self.X.get(key, default)
        if v is None:
            raise KeyError(f"composition variable '{key}' is required by this diffusion law "
                           f"but was not supplied (have: {sorted(self.X)})")
        return v

    @property
    def log_fo2_Pa(self) -> float:
        """log10 fO2 in Pascal (1 bar = 1e5 Pa)."""
        if self.log_fo2_bar is None:
            raise ValueError("this diffusion law needs fO2, but none was given")
        return self.log_fo2_bar + 5.0

    def replace(self, **kw) -> "Conditions":
        d = dict(T_K=self.T_K, P_Pa=self.P_Pa, log_fo2_bar=self.log_fo2_bar,
                 X=dict(self.X), axis=self.axis, angles_deg=self.angles_deg)
        d.update(kw)
        return Conditions(**d)


@dataclass
class Parameter:
    """One fitted parameter of a diffusion law."""
    name: str
    value: float
    sigma: float = 0.0
    unit: str = ""
    sigma_level: str = "1s"       # '1s' or '2s' as published
    description: str = ""

    @property
    def sigma_1s(self) -> float:
        return self.sigma / 2.0 if self.sigma_level == "2s" else self.sigma


@dataclass
class Range:
    lo: Optional[float] = None
    hi: Optional[float] = None
    unit: str = ""

    def contains(self, v) -> bool:
        v = float(np.min(v)), float(np.max(v))
        if self.lo is not None and v[0] < self.lo:
            return False
        if self.hi is not None and v[1] > self.hi:
            return False
        return True

    def __str__(self) -> str:
        lo = "-inf" if self.lo is None else f"{self.lo:g}"
        hi = "inf" if self.hi is None else f"{self.hi:g}"
        return f"{lo} to {hi} {self.unit}".strip()


# ---------------------------------------------------------------------------
@dataclass
class DiffusionCoefficient:
    """One published diffusion law, fully attributed."""

    key: str
    mineral: str
    species: str
    label: str
    citation: str                                   # key in diffusor.references
    equation_text: str                              # verbatim, as printed
    equation_number: str = ""
    func: Callable[["DiffusionCoefficient", Conditions, Dict[str, float]], object] = None
    params: Dict[str, Parameter] = field(default_factory=dict)
    covariance: Optional[np.ndarray] = None          # over ``cov_order``
    cov_order: Sequence[str] = ()
    sigma_logD: Optional[float] = None               # 1 sigma scatter of log10 D about the fit
    axis_factors: Dict[str, float] = field(default_factory=dict)   # D_axis / D_reference
    reference_axis: str = ""
    requires: Sequence[str] = ()                     # composition keys needed
    needs_fo2: bool = False
    fo2_unit: str = "Pa"                             # unit the law's fO2 term expects
    T_range: Range = field(default_factory=Range)
    P_range: Range = field(default_factory=Range)
    fo2_range: Range = field(default_factory=Range)
    X_range: Range = field(default_factory=Range)
    verified: bool = False
    verified_from: str = ""
    superseded_by: str = ""          # citation key of a newer calibration
    superseded_note: str = ""        # what changed and by how much
    secondary_citations: Sequence[str] = ()
    notes: str = ""
    recommended: bool = False

    # -- evaluation ---------------------------------------------------------
    def _values(self, overrides: Optional[Dict[str, float]] = None) -> Dict[str, float]:
        v = {k: p.value for k, p in self.params.items()}
        if overrides:
            v.update(overrides)
        return v

    def D_reference_axis(self, cond: Conditions, overrides=None):
        """D (m^2/s) along the law's reference axis, ignoring anisotropy."""
        return self.func(self, cond, self._values(overrides))

    def D(self, cond: Conditions, overrides=None):
        """D (m^2/s) along the traverse described by ``cond``.

        If ``cond.axis`` is given, the axis factor for that axis is applied.
        If ``cond.angles_deg`` is given, the direction-cosine relation of
        Costa & Chakraborty (2004) is used.  With neither, the reference axis
        of the publication is used unchanged.
        """
        D0 = self.D_reference_axis(cond, overrides)
        if cond.angles_deg is not None and self.axis_factors:
            from ..minerals.base import direction_factor
            fa = self.axis_factors.get("a", 1.0)
            fb = self.axis_factors.get("b", 1.0)
            fc = self.axis_factors.get("c", 1.0)
            return direction_factor(D0 * fa, D0 * fb, D0 * fc, *cond.angles_deg)
        if cond.axis and self.axis_factors:
            if cond.axis not in self.axis_factors:
                raise KeyError(f"{self.key} has no anisotropy factor for axis '{cond.axis}' "
                               f"(has {sorted(self.axis_factors)})")
            return D0 * self.axis_factors[cond.axis]
        return D0

    def log10_D(self, cond: Conditions, overrides=None):
        return np.log10(self.D(cond, overrides))

    # -- Monte Carlo --------------------------------------------------------
    def sample(self, rng: np.random.Generator, mode: str = "auto") -> Dict[str, float]:
        """Draw one set of parameter overrides for a Monte Carlo iteration.

        Returns a dict of parameter overrides; for ``logD_at_T`` the special
        key ``__log10_D_offset__`` is returned instead and applied by
        :meth:`D_sampled`.
        """
        mode = self.default_sampling_mode() if mode == "auto" else mode
        if mode == "covariance":
            if self.covariance is None:
                raise ValueError(f"{self.key} has no published covariance matrix")
            mean = np.array([self.params[k].value for k in self.cov_order], dtype=float)
            draw = rng.multivariate_normal(mean, np.asarray(self.covariance, dtype=float))
            return dict(zip(self.cov_order, draw))
        if mode == "logD_at_T":
            s = self.sigma_logD
            if s is None:
                raise ValueError(f"{self.key} has no published log10 D scatter")
            return {"__log10_D_offset__": float(rng.normal(0.0, s))}
        if mode == "independent":
            return {k: float(rng.normal(p.value, p.sigma_1s)) if p.sigma_1s > 0 else p.value
                    for k, p in self.params.items()}
        if mode == "none":
            return {}
        raise ValueError(f"unknown sampling mode '{mode}'")

    def default_sampling_mode(self) -> str:
        if self.covariance is not None:
            return "covariance"
        if self.sigma_logD is not None:
            return "logD_at_T"
        return "independent"

    def D_sampled(self, cond: Conditions, overrides: Dict[str, float]):
        """Evaluate D with a set of overrides that may contain the log-offset key."""
        ov = dict(overrides)
        off = ov.pop("__log10_D_offset__", 0.0)
        D = self.D(cond, ov if ov else None)
        return D * (10.0 ** off) if off else D

    # -- reporting ----------------------------------------------------------
    def check_conditions(self, cond: Conditions) -> List[str]:
        """Warnings for conditions outside the published calibration range."""
        w: List[str] = []
        if self.T_range.lo is not None and not self.T_range.contains(cond.T_K):
            w.append(f"T = {cond.T_K:.0f} K is outside the calibration range "
                     f"({self.T_range}) of {cite(self.citation)}")
        if cond.log_fo2_bar is not None and self.fo2_range.lo is not None:
            lf = cond.log_fo2_Pa if self.fo2_unit == "Pa" else cond.log_fo2_bar
            if not self.fo2_range.contains(lf):
                w.append(f"log fO2 = {lf:.2f} ({self.fo2_unit}) is outside the calibration "
                         f"range ({self.fo2_range}) of {cite(self.citation)}")
        if self.P_range.hi is not None and not self.P_range.contains(cond.P_Pa):
            w.append(f"P = {cond.P_Pa/1e9:.2f} GPa is outside the calibration range "
                     f"({self.P_range}) of {cite(self.citation)}")
        for key in self.requires:
            if key in cond.X and self.X_range.lo is not None:
                if not self.X_range.contains(cond.X[key]):
                    w.append(f"{key} is outside the calibration range ({self.X_range}) "
                             f"of {cite(self.citation)}")
        if not self.verified:
            w.append(f"{self.key}: coefficients NOT verified against the primary publication "
                     f"({self.verified_from or 'no secondary source recorded'}). Check before publishing.")
        if self.superseded_by:
            w.append(f"SUPERSEDED: a newer calibration exists -- {cite(self.superseded_by)}. "
                     f"{self.superseded_note}")
        return w

    def describe(self) -> str:
        ref = get_reference(self.citation)
        lines = [f"{self.label}",
                 f"  mineral/species : {self.mineral} / {self.species}",
                 f"  source          : {ref.full()}",
                 f"  equation        : {self.equation_number or '(unnumbered)'}",
                 f"     {self.equation_text}",
                 f"  D units         : m^2/s"]
        if self.needs_fo2:
            lines.append(f"  fO2 unit in law : {self.fo2_unit}")
        if self.params:
            lines.append("  parameters:")
            for p in self.params.values():
                s = f" +/- {p.sigma:g} ({p.sigma_level})" if p.sigma else ""
                lines.append(f"     {p.name} = {p.value:g}{s} {p.unit}  {p.description}")
        if self.sigma_logD is not None:
            lines.append(f"  scatter about the fit: {self.sigma_logD:g} log10 units (1 sigma)")
        if self.covariance is not None:
            lines.append(f"  covariance published for: {', '.join(self.cov_order)}")
        if self.axis_factors:
            fac = ", ".join(f"D_{k}/D_{self.reference_axis} = {v:g}" for k, v in self.axis_factors.items())
            lines.append(f"  anisotropy      : {fac}")
        lines.append(f"  calibration     : T {self.T_range}; P {self.P_range}; "
                     f"fO2 {self.fo2_range}; X {self.X_range}")
        lines.append(f"  verified        : {'YES -- ' + self.verified_from if self.verified else 'NO -- ' + (self.verified_from or 'unchecked')}")
        if self.superseded_by:
            lines.append(f"  SUPERSEDED BY   : {get_reference(self.superseded_by).full()}")
            lines.append(f"     {self.superseded_note}")
        if self.secondary_citations:
            lines.append("  see also        : " + "; ".join(cite(c) for c in self.secondary_citations))
        if self.notes:
            lines.append(f"  notes           : {self.notes}")
        return "\n".join(lines)

    def methods_sentence(self) -> str:
        """One sentence for the exported methods paragraph."""
        ref = get_reference(self.citation)
        eq = f", eq. {self.equation_number}" if self.equation_number else ""
        return (f"{self.species} diffusion in {self.mineral} was calculated with the "
                f"parameterisation of {ref.short()}{eq}.")


# --- helpers used by the individual laws --------------------------------------
def arrhenius(D0: float, Q_J: float, T_K: float, dV: float = 0.0,
              P_Pa: float = 1.0e5, P0_Pa: float = 1.0e5):
    """``D = D0 exp(-(Q + (P - P0) dV) / (R T))``.

    The standard Arrhenius form with an activation volume (Crank 1975 section 11;
    Costa et al. 2008 eq. 20, where the pressure term is written
    ``Q + P dV``).  ``dV`` in m^3/mol, ``Q_J`` in J/mol.
    """
    return D0 * np.exp(-(Q_J + (P_Pa - P0_Pa) * dV) / (R_GAS * T_K))
