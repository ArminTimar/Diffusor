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
The fit parameters are typically *positively correlated*: a higher ln D0 comes
with a higher Q.  Sampling them independently (as NIDIS does, Petrone et al.
2016) inflates the uncertainty at the temperature of interest by a large
factor.  Diffusor therefore offers three sampling modes, in order of
preference:

``covariance``
    Sample the parameter vector from a covariance matrix (this is the
    approach of Mutch et al. 2021, DFENS). The matrix is either published or
    built from a correlation the source states in words; each entry's
    uncertainty note says which.
``logD_at_T``
    Sample ln D directly at the working temperature using the paper's stated
    scatter about the Arrhenius line (e.g. "reproduces the data within 1 log
    unit", Mueller et al. 2013).  This is the fallback when no covariance
    is published; most papers constrain D at a given temperature in this form.
``independent``
    Sample each parameter independently.  Provided so that results of
    software that samples this way can be reproduced; never the default.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..constants import R_GAS
from ..references import cite, get as get_reference

LN10 = np.log(10.0)

# What a measured coefficient describes. A tracer coefficient and an exchange
# (interdiffusion) coefficient of the same elements can differ by orders of
# magnitude (Schaffer et al. 2014; Oeser et al. 2026), so every law states which
# one it is and the interface shows the definition.
TRANSPORT_KINDS = {
    "tracer": ("Tracer (self-) diffusion: the mobility of one species or isotope in a host of "
               "fixed chemical composition, measured by isotope exchange or tracer in-diffusion "
               "without a gradient in the host composition."),
    "interdiffusion": ("Interdiffusion (exchange): two or more major components exchange on one "
                       "site, so the flux of one is balanced by the others (Fe-Mg, Ca-Mg, Na-K, "
                       "NaSi-CaAl, Fe-Ti). This is the coefficient that relaxes major-element zoning."),
    "chemical": ("Chemical diffusion of a dilute (trace) element down its own concentration "
                 "gradient in a host of fixed major-element composition."),
    "effective": ("An effective coefficient that lumps several species, sites or mechanisms into "
                  "one D, or that the entry applies beyond the conditions it was measured at."),
}

TRACER_ADVICE = (
    "This is a tracer coefficient. It equals the chemical diffusion coefficient only for a "
    "dilute species that mixes ideally. Major-element zoning relaxes by interdiffusion; for an "
    "ideal binary exchange of equally charged ions D_AB = D*_A D*_B / (X_A D*_A + X_B D*_B) "
    "(Chakraborty & Ganguly 1992 eq. 2 in the binary limit; Oeser et al. 2026 eq. 6), which "
    "can differ from either tracer coefficient.")


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
    # '1s' or '2s' as published, or 'unstated' when the source prints a +/- without
    # saying what it is (standard deviation, standard error, 95 % bound ...). An
    # unstated level is sampled as if it were 1 sigma and is reported as unstated.
    sigma_level: str = "1s"
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
        if not np.all(np.isfinite(v)):
            return False
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
    sigma_logD_basis: str = ""                        # published / derived / assumed, see uncertainty_basis.py
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
    kind: str = ""  # a key of TRANSPORT_KINDS; required
    model_family: str = "scalar_fickian"
    transported_variable: str = ""
    reference_state: str = ""
    calibration_notes: Sequence[str] = ()
    uncertainty_note: str = ""
    # Independent principal laws, evaluated BEFORE projection. Ratios may vary with T.
    principal_funcs: Dict[str, Callable] = field(default_factory=dict)
    allowed_axes: Sequence[str] = ()
    orientation_required: bool = False
    fixed_temperature_K: Optional[float] = None
    # Keys of the tracer laws an interdiffusion law was computed from (diffusor.coefficients.transport)
    derived_from: Sequence[str] = ()

    def __post_init__(self):
        if self.kind not in TRANSPORT_KINDS:
            raise ValueError(f"{self.key}: transport kind {self.kind!r} must be one of {sorted(TRANSPORT_KINDS)}")

    @property
    def kind_description(self) -> str:
        return TRANSPORT_KINDS[self.kind]

    @property
    def validation_level(self) -> str:
        return "source_transcription_checked" if self.verified else "transcribed_unverified"

    # -- evaluation ---------------------------------------------------------
    def _values(self, overrides: Optional[Dict[str, float]] = None) -> Dict[str, float]:
        v = {k: p.value for k, p in self.params.items()}
        if overrides:
            v.update(overrides)
        return v

    def D_reference_axis(self, cond: Conditions, overrides=None):
        """D (m^2/s) along the law's reference axis, ignoring anisotropy."""
        self._validate_conditions(cond)
        return self.func(self, cond, self._values(overrides))

    def _validate_conditions(self, cond):
        if not np.isfinite(cond.T_K) or cond.T_K <= 0:
            raise ValueError("temperature must be finite and positive in kelvin")
        if not np.isfinite(cond.P_Pa) or cond.P_Pa < 0:
            raise ValueError("pressure must be finite and nonnegative in Pa")
        if self.needs_fo2 and (cond.log_fo2_bar is None or not np.isfinite(cond.log_fo2_bar)):
            raise ValueError(f"{self.key} requires a finite log10 fO2 (bar)")
        if self.fixed_temperature_K is not None and not np.isclose(
                cond.T_K, self.fixed_temperature_K, rtol=0, atol=1e-6):
            raise ValueError(f"{self.key} is calibrated only at {self.fixed_temperature_K - 273.15:g} C; "
                             "no temperature dependence was measured")

    def D(self, cond: Conditions, overrides=None):
        """D (m^2/s) along the traverse described by ``cond``.

        If ``cond.axis`` is given, the axis factor for that axis is applied.
        If ``cond.angles_deg`` is given, the direction-cosine relation
        D_trav = D_a cos^2(alpha) + D_b cos^2(beta) + D_c cos^2(gamma) is used
        (Costa & Chakraborty 2004, EPSL 227, eq. 5, there taken from Philibert
        1991; also stated in Costa, Dohmen & Chakraborty 2008, RiMG 69, p. 574, who
        note that the principal axes need not coincide with the crystallographic
        axes).  With neither, the reference axis
        of the publication is used unchanged.
        """
        self._validate_conditions(cond)
        if cond.axis is not None and cond.angles_deg is not None:
            raise ValueError("specify either an axis or direction angles, not both")
        if self.orientation_required and cond.axis is None and cond.angles_deg is None:
            raise ValueError(f"{self.key} requires an explicit traverse orientation")
        if self.allowed_axes:
            if cond.angles_deg is not None or (cond.axis and cond.axis not in self.allowed_axes):
                raise ValueError(f"{self.key} is calibrated only for axes {tuple(self.allowed_axes)}; "
                                 "a full diffusion tensor is unavailable")
        if self.principal_funcs:
            p = self._values(overrides)
            if cond.angles_deg is not None:
                from ..minerals.base import direction_factor
                return direction_factor(*(self.principal_funcs[a](self, cond, p) for a in ("a", "b", "c")),
                                        *cond.angles_deg)
            axis = cond.axis or self.reference_axis
            if axis not in self.principal_funcs:
                raise ValueError(f"no principal diffusivity for axis {axis!r}")
            return self.principal_funcs[axis](self, cond, p)
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
                raise ValueError(f"{self.key} has no covariance matrix")
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
        """The sampling used when the user leaves the choice to Diffusor.

        Independent sampling of log D0 and Q is never the default: the two are
        correlated regression parameters, and drawing them separately spreads
        log D by orders of magnitude away from the centroid of the data. A law
        with only marginal errors is held fixed unless the user asks otherwise.
        """
        if self.covariance is not None:
            return "covariance"
        if self.sigma_logD is not None:
            return "logD_at_T"
        return "none"

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
        if not self.T_range.contains(cond.T_K):
            w.append(f"T = {cond.T_K:.0f} K is outside the calibration range "
                     f"({self.T_range}) of {cite(self.citation)}")
        if cond.log_fo2_bar is not None:
            lf = (cond.log_fo2_Pa if self.fo2_unit == "Pa" else
                  cond.log_fo2_bar - np.log10(1.01325) if self.fo2_unit == "atm" else cond.log_fo2_bar)
            if not self.fo2_range.contains(lf):
                w.append(f"log fO2 = {lf:.2f} ({self.fo2_unit}) is outside the calibration "
                         f"range ({self.fo2_range}) of {cite(self.citation)}")
        if not self.P_range.contains(cond.P_Pa):
            w.append(f"P = {cond.P_Pa/1e9:.2f} GPa is outside the calibration range "
                     f"({self.P_range}) of {cite(self.citation)}")
        for key in self.requires:
            if key not in cond.X:
                w.append(f"{self.key}: host composition {key} was not supplied; inspect the reference state before interpreting this fit.")
            if key in cond.X and (self.X_range.lo is not None or self.X_range.hi is not None):
                if not self.X_range.contains(cond.X[key]):
                    w.append(f"{key} is outside the calibration range ({self.X_range}) "
                             f"of {cite(self.citation)}")
        if not self.verified:
            w.append(f"{self.key}: coefficients NOT verified against the primary publication "
                     f"({self.verified_from or 'no secondary source recorded'}). Check before publishing.")
        if self.superseded_by:
            w.append(f"SUPERSEDED: a newer calibration exists -- {cite(self.superseded_by)}. "
                     f"{self.superseded_note}")
        w.extend(self.calibration_notes)
        if self.uncertainty_note:
            w.append(self.uncertainty_note)
        if self.reference_axis and not cond.axis and cond.angles_deg is None:
            w.append(f"Orientation unspecified: using the published {self.reference_axis} direction; "
                     "this is not an orientation average.")
        return w

    def describe(self) -> str:
        ref = get_reference(self.citation)
        lines = [f"{self.label}",
                 f"  mineral/species : {self.mineral} / {self.species}",
                 f"  source          : {ref.full()}",
                 f"  equation        : {self.equation_number or '(unnumbered)'}",
                 f"     {self.equation_text}",
                 f"  D units         : m^2/s"]
        lines.extend([f"  coefficient kind: {self.kind} -- {self.kind_description}",
                      f"  model family    : {self.model_family}",
                      f"  validation      : {self.validation_level}"])
        if self.derived_from:
            lines.append(f"  computed from   : {', '.join(self.derived_from)}")
        if self.transported_variable:
            lines.append(f"  state variable  : {self.transported_variable}")
        if self.reference_state:
            lines.append(f"  reference state : {self.reference_state}")
        if self.principal_funcs:
            lines.append("  anisotropy      : independent temperature-dependent principal laws")
        if self.allowed_axes:
            lines.append(f"  measured axes   : {', '.join(self.allowed_axes)} (no tensor inferred)")
        lines.extend(f"  limitation      : {s}" for s in self.calibration_notes)
        if self.uncertainty_note:
            lines.append(f"  uncertainty     : {self.uncertainty_note}")
        if self.needs_fo2:
            lines.append(f"  fO2 unit in law : {self.fo2_unit}")
        if self.params:
            lines.append("  parameters:")
            for p in self.params.values():
                s = f" +/- {p.sigma:g} ({p.sigma_level})" if p.sigma else ""
                lines.append(f"     {p.name} = {p.value:g}{s} {p.unit}  {p.description}")
        if self.sigma_logD is not None:
            lines.append(f"  scatter about the fit: {self.sigma_logD:g} log10 units (1 sigma); "
                         f"{self.sigma_logD_basis or 'basis not recorded'}")
        if self.covariance is not None:
            lines.append(f"  covariance (see uncertainty note for its basis): {', '.join(self.cov_order)}")
        if self.axis_factors:
            fac = ", ".join(f"D_{k}/D_{self.reference_axis} = {v:g}" for k, v in self.axis_factors.items())
            lines.append(f"  anisotropy      : {fac}")
        lines.append(f"  calibration     : T {self.T_range}, P {self.P_range}, "
                     f"fO2 {self.fo2_range}, X {self.X_range}")
        lines.append(f"  verified        : {'YES -- ' + self.verified_from if self.verified else 'NO -- ' + (self.verified_from or 'unchecked')}")
        if self.superseded_by:
            lines.append(f"  SUPERSEDED BY   : {get_reference(self.superseded_by).full()}")
            lines.append(f"     {self.superseded_note}")
        if self.secondary_citations:
            lines.append("  see also        : " + ", ".join(cite(c) for c in self.secondary_citations))
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

    The standard Arrhenius form with an activation volume (Costa et al. 2008
    eq. 20, where the pressure term is written
    ``Q + P dV``).  ``dV`` in m^3/mol, ``Q_J`` in J/mol.
    """
    return D0 * np.exp(-(Q_J + (P_Pa - P0_Pa) * dV) / (R_GAS * T_K))
