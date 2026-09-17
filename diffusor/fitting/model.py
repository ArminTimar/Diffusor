"""The forward model: everything needed to turn a time into a predicted profile.

A :class:`DiffusionModel` bundles the mineral, coefficient, conditions,
geometry, initial and boundary conditions, thermal history and beam
convolution.  It exposes :meth:`profile`, which maps a time (seconds) to a
concentration profile on the measurement coordinates, and is used identically
by the least-squares fit and by every Monte Carlo draw.

Analytical vs numerical
-----------------------
The analytical route (Crank 1975 closed forms) is used when it is exactly
valid: constant D (no composition dependence), a step or plateau initial
condition, and a geometry with a published series solution.  Otherwise the
Crank-Nicolson solver runs.  ``DiffusionModel.can_use_analytical()`` reports
which applies and why, so the user is never silently given the wrong one.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

import numpy as np

from ..coefficients.base import Conditions, DiffusionCoefficient
from ..solvers import analytical
from ..solvers.boundary import BoundaryCondition, zero_flux
from ..solvers.convolution import gaussian_convolve, resolution_warning
from ..solvers.geometry import Geometry, suggest_grid
from ..solvers.history import ThermalHistory, effective_Dt
from ..solvers.initial import InitialCondition
from ..solvers.numerical import solve_1d

# Diffusor works internally in micrometres and seconds: D is converted from
# m^2/s to um^2/s once, which keeps the grid numbers O(1-1000) and avoids the
# ill-conditioning of a metre-based grid with D ~ 1e-20 m^2/s.
M2_PER_S_TO_UM2_PER_S = 1.0e12


@dataclass
class DiffusionModel:
    coefficient: DiffusionCoefficient
    conditions: Conditions
    initial: InitialCondition
    geometry: Geometry = field(default_factory=Geometry)
    bc_left: BoundaryCondition = field(default_factory=zero_flux)
    bc_right: BoundaryCondition = field(default_factory=zero_flux)
    history: Optional[ThermalHistory] = None
    beam_sigma_um: float = 0.0
    n_nodes: int = 401
    x_grid: Optional[np.ndarray] = None
    composition_dependent: bool = True
    comp_key: Optional[str] = None       # which X key the profile itself is
    an_profile: Optional[np.ndarray] = None       # plagioclase X_An on the grid
    activity_theta: float = 0.0                   # A_i/(RT) for the plag term
    force_numerical: bool = False

    # -- grid ---------------------------------------------------------------
    def grid(self, x_data) -> np.ndarray:
        if self.x_grid is not None:
            return self.x_grid
        return suggest_grid(x_data, self.n_nodes)

    # -- diffusivity --------------------------------------------------------
    def _D_um2s(self, C_nodes, T_K: float, overrides=None):
        """D in um^2/s at the given nodal compositions."""
        cond = self.conditions.replace(T_K=T_K)
        if self.composition_dependent and self.comp_key:
            X = dict(cond.X)
            X[self.comp_key] = C_nodes
            cond = cond.replace(X=X)
        D = self.coefficient.D_sampled(cond, overrides or {})
        return np.asarray(D, dtype=float) * M2_PER_S_TO_UM2_PER_S

    def D_bulk(self, overrides=None, C_ref: Optional[float] = None) -> float:
        """A single representative D (m^2/s) for reporting and for the analytical route."""
        cond = self.conditions
        if self.composition_dependent and self.comp_key and C_ref is not None:
            X = dict(cond.X)
            X[self.comp_key] = C_ref
            cond = cond.replace(X=X)
        return float(np.mean(np.atleast_1d(self.coefficient.D_sampled(cond, overrides or {}))))

    # -- analytical applicability -------------------------------------------
    def can_use_analytical(self) -> Tuple[bool, str]:
        if self.force_numerical:
            return False, "numerical solver requested explicitly"
        if self.composition_dependent and self.comp_key:
            return False, ("D depends on the composition being modelled, so no closed-form "
                           "solution exists (Crank 1975 section 7.2). The numerical solver is used.")
        if self.activity_theta:
            return False, ("the plagioclase activity term couples D to the anorthite gradient "
                           "(Costa et al. 2003 eq. 7), which has no closed form")
        if self.initial.kind not in ("step", "plateau_rim"):
            return False, f"initial condition '{self.initial.kind}' has no closed-form solution"
        if self.geometry.kind != "plane":
            return False, ("closed forms for cylinder and sphere assume a uniform initial "
                           "profile and a fixed surface concentration. Call "
                           "diffusor.solvers.analytical directly for that case.")
        if self.initial.params.get("smooth"):
            return False, "a smoothed initial step has no closed-form solution"
        return True, "step initial condition, constant D, plane geometry: Crank (1975) eq. 2.14"

    # -- forward model ------------------------------------------------------
    def profile(self, t_seconds: float, x_out, overrides=None) -> np.ndarray:
        """Predicted concentration at ``x_out`` after ``t_seconds``."""
        x_out = np.asarray(x_out, dtype=float)
        ok, _ = self.can_use_analytical()
        if ok:
            C = self._profile_analytical(t_seconds, x_out, overrides)
            if self.beam_sigma_um > 0:
                xf = self.grid(x_out)
                Cf = self._profile_analytical(t_seconds, xf, overrides)
                Cf = gaussian_convolve(xf, Cf, self.beam_sigma_um)
                C = np.interp(x_out, xf, Cf)
            return C
        return self._profile_numerical(t_seconds, x_out, overrides)

    def _profile_analytical(self, t_seconds: float, x_out, overrides=None) -> np.ndarray:
        p = self.initial.params
        C_left, C_right = p["C_left"], p["C_right"]
        Dt = self._effective_Dt(t_seconds, overrides, C_ref=0.5 * (C_left + C_right))
        return analytical.step_infinite(x_out, p["x0"], C_left, C_right, Dt)

    def _effective_Dt(self, t_seconds: float, overrides=None, C_ref=None) -> float:
        """Dt in um^2, integrating over the thermal history if there is one."""
        if self.history is None or self.history.is_isothermal:
            return self.D_bulk(overrides, C_ref) * M2_PER_S_TO_UM2_PER_S * t_seconds
        hist = self.history.shifted_to_end(t_seconds)

        def D_of_T(T_array):
            out = []
            for T in np.atleast_1d(T_array):
                cond = self.conditions.replace(T_K=float(T))
                if self.comp_key and C_ref is not None:
                    X = dict(cond.X)
                    X[self.comp_key] = C_ref
                    cond = cond.replace(X=X)
                out.append(float(np.mean(np.atleast_1d(
                    self.coefficient.D_sampled(cond, overrides or {})))))
            return np.array(out) * M2_PER_S_TO_UM2_PER_S

        return effective_Dt(D_of_T, hist, t_seconds)

    def _profile_numerical(self, t_seconds: float, x_out, overrides=None) -> np.ndarray:
        x = self.grid(x_out)
        C0 = self.initial.evaluate(x)
        hist = (self.history.shifted_to_end(t_seconds) if self.history is not None
                else ThermalHistory.isothermal(self.conditions.T_K, t_seconds))

        def D_func(C, T_K):
            return self._D_um2s(C, T_K, overrides)

        res = solve_1d(x, C0, D_func, t_seconds, m=self.geometry.m,
                       bc_left=self.bc_left, bc_right=self.bc_right,
                       history=hist, theta_time=0.5,
                       an_profile=self.an_profile, theta_activity=self.activity_theta)
        C = res.C_final
        if self.beam_sigma_um > 0:
            C = gaussian_convolve(x, C, self.beam_sigma_um)
        return np.interp(np.asarray(x_out, dtype=float), x, C)

    # -- diagnostics ---------------------------------------------------------
    def warnings(self, t_seconds: Optional[float] = None) -> List[str]:
        w = list(self.coefficient.check_conditions(self.conditions))
        if t_seconds is not None and self.beam_sigma_um > 0:
            Dt = self._effective_Dt(t_seconds)
            msg = resolution_warning(Dt, self.beam_sigma_um)
            if msg:
                w.append(msg)
        if t_seconds is not None and self.initial.kind in ("step", "plateau_rim"):
            x0 = self.initial.params.get("x0", self.initial.params.get("rim_start"))
            if x0 is not None and self.x_grid is not None:
                half = min(abs(float(self.x_grid[-1]) - x0), abs(x0 - float(self.x_grid[0])))
                L = 2.0 * np.sqrt(max(self._effective_Dt(t_seconds), 0.0))
                if L > 0.6 * half:
                    w.append(
                        f"the diffusion length 2*sqrt(Dt) = {L:.1f} um is a large fraction of "
                        f"the distance from the interface to the end of the profile "
                        f"({half:.1f} um). The far field has been reached, so the semi-infinite "
                        "assumption is breaking down. Measure a longer traverse or model the "
                        "whole crystal with an explicit geometry.")
        if self.geometry.kind == "plane":
            w.append("1-D modelling of a 3-D crystal gives a maximum estimate of the time. "
                     "Sectioning can bias it further (Shea et al. 2015, "
                     "Krimer & Costa 2017).")
        return w
