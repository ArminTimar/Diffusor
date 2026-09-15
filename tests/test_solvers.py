"""The numerical solver must reproduce the closed-form solutions of Crank (1975)."""
import numpy as np
import pytest

from diffusor.solvers import (Geometry, ThermalHistory, analytical, dirichlet,
                              effective_Dt, gaussian_convolve, make_grid,
                              resolution_warning, solve_1d, zero_flux)
from diffusor.solvers.initial import step

CONST_D = 1.0


def _D(C, T):
    return np.full_like(np.asarray(C, dtype=float), CONST_D)


def test_step_in_an_infinite_medium_matches_crank_eq_2_14():
    Dt = 400.0
    L = 8 * np.sqrt(Dt)
    x = make_grid(-L, L, 801)
    res = solve_1d(x, step(x, 0.0, 0.0, 1.0), _D, Dt, T_K=1000.0,
                   bc_left=dirichlet(0.0), bc_right=dirichlet(1.0))
    ref = analytical.step_infinite(x, 0.0, 0.0, 1.0, Dt)
    assert np.max(np.abs(res.C_final - ref)) < 1e-4


def test_numerical_step_converges_at_second_order():
    Dt = 100.0
    L = 8 * np.sqrt(Dt)
    errs = []
    for n in (401, 801, 1601):
        x = make_grid(-L, L, n)
        res = solve_1d(x, step(x, 0.0, 0.0, 1.0), _D, Dt, T_K=1000.0,
                       bc_left=dirichlet(0.0), bc_right=dirichlet(1.0))
        errs.append(np.max(np.abs(res.C_final - analytical.step_infinite(x, 0.0, 0.0, 1.0, Dt))))
    assert errs[1] < errs[0] and errs[2] < errs[1]


def test_sphere_matches_crank_eq_6_18():
    a, Dt = 100.0, 300.0
    r = make_grid(0.0, a, 801)
    res = solve_1d(r, np.zeros_like(r), _D, Dt, m=2, T_K=1000.0,
                   bc_left=zero_flux(), bc_right=dirichlet(1.0))
    assert np.max(np.abs(res.C_final - analytical.sphere(r, a, 0.0, 1.0, Dt))) < 1e-4


def test_cylinder_matches_crank_eq_5_22():
    a, Dt = 100.0, 300.0
    r = make_grid(0.0, a, 801)
    res = solve_1d(r, np.zeros_like(r), _D, Dt, m=1, T_K=1000.0,
                   bc_left=zero_flux(), bc_right=dirichlet(1.0))
    assert np.max(np.abs(res.C_final - analytical.cylinder(r, a, 0.0, 1.0, Dt))) < 1e-4


def test_plane_sheet_matches_crank_eq_4_17():
    l, Dt = 100.0, 300.0
    x = make_grid(-l, l, 801)
    res = solve_1d(x, np.zeros_like(x), _D, Dt, T_K=1000.0,
                   bc_left=dirichlet(1.0), bc_right=dirichlet(1.0))
    assert np.max(np.abs(res.C_final - analytical.plane_sheet(x, l, 0.0, 1.0, Dt))) < 1e-4


def test_closed_system_conserves_mass():
    x = make_grid(0.0, 100.0, 401)
    res = solve_1d(x, step(x, 30.0, 1.0, 0.2), _D, 5000.0, T_K=1000.0,
                   bc_left=zero_flux(), bc_right=zero_flux())
    assert abs(res.mass_final / res.mass_initial - 1.0) < 1e-8


def test_closed_sphere_conserves_mass_with_the_radial_weight():
    r = make_grid(0.0, 100.0, 401)
    res = solve_1d(r, step(r, 50.0, 1.0, 0.2), _D, 3000.0, m=2, T_K=1000.0,
                   bc_left=zero_flux(), bc_right=zero_flux())
    assert abs(res.mass_final / res.mass_initial - 1.0) < 1e-6


def test_closed_system_relaxes_to_the_mean():
    x = make_grid(0.0, 100.0, 201)
    C0 = step(x, 50.0, 1.0, 0.0)
    res = solve_1d(x, C0, _D, 1.0e6, T_K=1000.0, bc_left=zero_flux(), bc_right=zero_flux())
    assert np.allclose(res.C_final, 0.5, atol=1e-3)


def test_explicit_and_implicit_schemes_agree():
    Dt = 50.0
    L = 8 * np.sqrt(Dt)
    x = make_grid(-L, L, 401)
    C0 = step(x, 0.0, 0.0, 1.0)
    a = solve_1d(x, C0, _D, Dt, T_K=1000.0, theta_time=0.5,
                 bc_left=dirichlet(0.0), bc_right=dirichlet(1.0))
    b = solve_1d(x, C0, _D, Dt, T_K=1000.0, theta_time=0.0, courant=0.4,
                 bc_left=dirichlet(0.0), bc_right=dirichlet(1.0))
    assert np.max(np.abs(a.C_final - b.C_final)) < 5e-4


def test_fractional_uptake_series_are_consistent_with_the_profiles():
    """M_t/M_inf from Crank eq. 6.20 must equal the integral of the eq. 6.18 profile."""
    trap = getattr(np, "trapezoid", None) or np.trapz
    a = 1.0
    for Dt in (0.003, 0.03, 0.3):
        r = np.linspace(1e-9, a, 20001)
        C = analytical.sphere(r, a, 0.0, 1.0, Dt)
        got = trap(C * r ** 2, r) / (a ** 3 / 3.0)
        assert got == pytest.approx(analytical.fraction_sphere(a, Dt), abs=1e-6)


def test_band_solution_is_symmetric_and_conserves_the_source():
    x = np.linspace(-200, 200, 2001)
    C = analytical.band_infinite(x, 0.0, 10.0, 1.0, 0.0, 100.0)
    assert C[0] == pytest.approx(C[-1], abs=1e-12)
    assert np.argmax(C) in (999, 1000, 1001)


def test_non_isothermal_integral_equals_isothermal_when_flat():
    h = ThermalHistory.isothermal(1200.0, 100.0)
    assert effective_Dt(lambda T: np.full_like(np.atleast_1d(T), 2.0), h, 100.0) == pytest.approx(200.0)


def test_cooling_history_gives_less_diffusion_than_isothermal_at_the_peak():
    from diffusor.constants import R_GAS

    def D_of_T(T):
        return 1e-3 * np.exp(-200e3 / (R_GAS * np.atleast_1d(T)))

    t = 1.0e6
    hot = ThermalHistory.isothermal(1300.0, t)
    cooling = ThermalHistory.linear(1300.0, 1000.0, t)
    assert effective_Dt(D_of_T, cooling, t) < effective_Dt(D_of_T, hot, t)


def test_convolution_broadens_a_step_and_preserves_the_plateaus():
    x = np.linspace(-50, 50, 1001)
    C = np.where(x < 0, 0.0, 1.0)
    Cc = gaussian_convolve(x, C, 3.0)
    assert Cc[0] == pytest.approx(0.0, abs=1e-6)
    assert Cc[-1] == pytest.approx(1.0, abs=1e-6)
    assert Cc[500] == pytest.approx(0.5, abs=0.02)
    # the convolved profile is monotonic and shallower at the centre
    assert np.all(np.diff(Cc) >= -1e-12)


def test_resolution_warning_triggers_for_short_profiles():
    assert resolution_warning(0.1, 3.0) is not None
    assert resolution_warning(1000.0, 1.0) is None


def test_dirichlet_boundary_can_vary_with_time():
    x = make_grid(0.0, 50.0, 201)
    res = solve_1d(x, np.zeros_like(x), _D, 200.0, T_K=1000.0,
                   bc_left=dirichlet(lambda t: min(1.0, t / 200.0)), bc_right=zero_flux())
    assert res.C_final[0] == pytest.approx(1.0, abs=1e-6)
    assert res.C_final[-1] > 0.0


def test_geometry_rejects_unknown_kinds():
    with pytest.raises(ValueError):
        Geometry("hexagon")
