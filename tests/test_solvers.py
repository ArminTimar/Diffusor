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


def test_step_between_nodes_gets_the_cell_average():
    """A sharp interface inside a node's cell gives that node the cell-weighted average,
    not the whole plateau, so the discretised step is not offset by up to dx/2."""
    x = np.linspace(0.0, 1.0, 11)
    C = step(x, 0.33, 0.0, 1.0)
    assert C[3] == pytest.approx(0.2)                  # cell [0.25, 0.35] is 20 % right of 0.33
    assert C[2] == 0.0 and C[4] == 1.0
    assert step(x, 0.3, 0.0, 1.0)[3] == pytest.approx(0.5)   # 0.3 is not exactly a float node
    assert np.allclose(step(np.linspace(0, 10, 11), 5.0, 0.0, 1.0)[4:7], [0.0, 0.5, 1.0])


def test_off_node_step_is_as_accurate_as_an_on_node_step():
    Dt = 4.0
    x = make_grid(-40.0, 40.0, 161)
    for x0 in (0.13, 0.37):
        res = solve_1d(x, step(x, x0, 0.0, 1.0), _D, Dt, T_K=1000.0,
                       bc_left=dirichlet(0.0), bc_right=dirichlet(1.0))
        assert np.max(np.abs(res.C_final - analytical.step_infinite(x, x0, 0.0, 1.0, Dt))) < 2e-3


def test_multi_step_edges_between_nodes_are_cell_averaged():
    from diffusor.solvers.initial import multi_step
    x = np.linspace(0.0, 1.0, 11)
    C = multi_step(x, [0.33, 0.77], [0.0, 1.0, 3.0])
    assert C[3] == pytest.approx(0.2) and C[7] == 1.0 and C[8] == pytest.approx(1.0 + 2.0 * 0.8)
    assert C[0] == 0.0 and C[10] == 3.0


def test_source_comments_carry_no_stale_crank_locators():
    """Crank (1975) locators checked against the book: 5.4, 6.3, 7.7, 8.31, 8.33 and 8.35
    are not the equations these comments used to cite."""
    import pathlib
    import re
    import diffusor
    root = pathlib.Path(diffusor.__file__).parent
    stale = re.compile(r"Crank[^\n]{0,60}eqs?\.? (5\.4|6\.3|7\.7|8\.31|8\.33|8\.35)\b")
    files = (list(root.glob("solvers/*.py")) + list(root.glob("dataio/*.py"))
             + list(root.glob("fitting/*.py")))
    hits = [str(p) for p in files if stale.search(p.read_text(encoding="utf-8"))]
    assert not hits, hits


def test_resolution_warning_says_the_threshold_is_diffusors_own():
    msg = resolution_warning(0.1, 3.0)
    assert "rule of thumb" in msg and "gradient width/spot size" in msg
