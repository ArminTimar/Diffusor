"""Fitting and Monte Carlo behaviour on synthetic profiles with a known answer."""
import numpy as np
import pytest

from diffusor.coefficients import Conditions, get
from diffusor.constants import SEC_PER_YEAR
from diffusor.fitting import (DiffusionModel, UncertaintyBudget, contributions,
                              fit_time, run_montecarlo)
from diffusor.solvers import Geometry, InitialCondition, dirichlet
from diffusor.thermo import log_fo2_from_delta


def make_model(comp_dependent=True, n_nodes=201, T_C=950.0):
    T = T_C + 273.15
    coef = get("opx_FeMg_dohmen2016")
    cond = Conditions(T_K=T, P_Pa=2.0e8, log_fo2_bar=log_fo2_from_delta("NNO", 1.0, T),
                      X={"XFe": 0.25}, axis="c")
    ic = InitialCondition("step", {"x0": 0.0, "C_left": 0.30, "C_right": 0.18})
    return DiffusionModel(coefficient=coef, conditions=cond, initial=ic,
                          geometry=Geometry("plane"),
                          bc_left=dirichlet(0.30), bc_right=dirichlet(0.18),
                          comp_key="XFe" if comp_dependent else None,
                          composition_dependent=comp_dependent, n_nodes=n_nodes)


def test_analytical_route_is_used_only_when_it_is_valid():
    assert make_model(comp_dependent=True).can_use_analytical()[0] is False
    ok, why = make_model(comp_dependent=False).can_use_analytical()
    assert ok and "Crank" in why


def test_fit_recovers_a_known_time_without_noise():
    m = make_model(comp_dependent=False)
    t_true = 3.0 * SEC_PER_YEAR
    x = np.linspace(-40, 40, 61)
    C = m.profile(t_true, x)
    r = fit_time(m, x, C)
    assert r.t_seconds == pytest.approx(t_true, rel=0.02)


def test_fit_recovers_a_known_time_with_composition_dependent_D():
    m = make_model(comp_dependent=True)
    t_true = 2.0 * SEC_PER_YEAR
    x = np.linspace(-40, 40, 41)
    C = m.profile(t_true, x)
    r = fit_time(m, x, C)
    assert r.t_seconds == pytest.approx(t_true, rel=0.05)


def test_fit_recovers_a_known_time_through_noise():
    rng = np.random.default_rng(3)
    m = make_model(comp_dependent=False)
    t_true = 1.0 * SEC_PER_YEAR
    x = np.linspace(-40, 40, 61)
    sigma = np.full_like(x, 0.003)
    C = m.profile(t_true, x) + rng.normal(0.0, sigma)
    r = fit_time(m, x, C, sigma)
    assert 0.6 < r.t_seconds / t_true < 1.6
    assert r.stats.reduced_chi2 < 3.0


def test_longer_times_give_wider_profiles():
    m = make_model(comp_dependent=False)
    x = np.linspace(-60, 60, 201)
    widths = []
    for t in (0.5 * SEC_PER_YEAR, 5.0 * SEC_PER_YEAR):
        C = m.profile(t, x)
        widths.append(np.sum((C > 0.19) & (C < 0.29)))
    assert widths[1] > widths[0]


def test_free_x0_is_recovered():
    m = make_model(comp_dependent=False)
    t_true = 2.0 * SEC_PER_YEAR
    x = np.linspace(-40, 40, 81)
    m_shift = make_model(comp_dependent=False)
    m_shift.initial.params["x0"] = 6.0
    C = m_shift.profile(t_true, x)
    r = fit_time(m, x, C, free_parameters=("t", "x0"))
    assert r.free["x0"] == pytest.approx(6.0, abs=1.5)


def test_beam_convolution_lengthens_the_apparent_time_if_ignored():
    """A profile blurred by the beam looks more diffused; ignoring the blur
    inflates the retrieved time (Bradshaw & Kent 2017)."""
    sharp = make_model(comp_dependent=False)
    blurred = make_model(comp_dependent=False)
    blurred.beam_sigma_um = 3.0
    t_true = 0.5 * SEC_PER_YEAR
    x = np.linspace(-40, 40, 81)
    C_blurred = blurred.profile(t_true, x)
    naive = fit_time(sharp, x, C_blurred)
    corrected = fit_time(blurred, x, C_blurred)
    assert naive.t_seconds > corrected.t_seconds
    assert corrected.t_seconds == pytest.approx(t_true, rel=0.1)


def test_montecarlo_spread_grows_with_temperature_uncertainty():
    m = make_model(comp_dependent=False)
    t_true = 2.0 * SEC_PER_YEAR
    x = np.linspace(-40, 40, 41)
    sigma = np.full_like(x, 0.003)
    C = m.profile(t_true, x)
    narrow = run_montecarlo(m, x, C, sigma, n_draws=60, seed=1,
                            budget=UncertaintyBudget(sigma_T_K=5.0, sample_coefficient=False,
                                                     sample_measurement_noise=False))
    wide = run_montecarlo(m, x, C, sigma, n_draws=60, seed=1,
                          budget=UncertaintyBudget(sigma_T_K=40.0, sample_coefficient=False,
                                                   sample_measurement_noise=False))
    assert wide.sigma_log10 > 2.0 * narrow.sigma_log10


def test_montecarlo_median_is_close_to_the_best_fit():
    m = make_model(comp_dependent=False)
    t_true = 2.0 * SEC_PER_YEAR
    x = np.linspace(-40, 40, 41)
    sigma = np.full_like(x, 0.003)
    C = m.profile(t_true, x)
    mc = run_montecarlo(m, x, C, sigma, n_draws=80, seed=5,
                        budget=UncertaintyBudget(sigma_T_K=15.0, sample_coefficient=False,
                                                 sample_measurement_noise=True))
    assert 0.5 < mc.median / mc.t_best < 2.0
    assert mc.p16 < mc.median < mc.p84


def test_zero_uncertainty_budget_gives_no_spread():
    m = make_model(comp_dependent=False)
    x = np.linspace(-40, 40, 41)
    C = m.profile(SEC_PER_YEAR, x)
    mc = run_montecarlo(m, x, C, np.full_like(x, 0.003), n_draws=25, seed=2,
                        budget=UncertaintyBudget(sample_coefficient=False,
                                                 sample_measurement_noise=False))
    assert mc.sigma_log10 < 1e-6


def test_temperature_and_fo2_stay_correlated_through_the_buffer():
    """Sampling T must move fO2 with it when fO2 is given as a buffer offset."""
    from diffusor.fitting.montecarlo import _draw_conditions
    m = make_model(comp_dependent=False)
    budget = UncertaintyBudget(sigma_T_K=30.0, buffer="NNO", delta_buffer=1.0,
                               sigma_delta_buffer=0.0)
    rng = np.random.default_rng(0)
    Ts, fs = [], []
    for _ in range(300):
        c = _draw_conditions(m, budget, rng)
        Ts.append(c.T_K)
        fs.append(c.log_fo2_bar)
    r = np.corrcoef(Ts, fs)[0, 1]
    assert r > 0.99, "fO2 must track temperature through the buffer equation"


def test_independent_sampling_overestimates_relative_to_logD_at_T():
    """Sampling ln D0 and Q independently inflates the spread, which is exactly
    the flaw this project set out to avoid."""
    coef = get("opx_FeMg_dohmen2016")
    rng = np.random.default_rng(4)
    cond = Conditions(T_K=1223.15, log_fo2_bar=-12.0, X={"XFe": 0.15}, axis="c")
    indep = [np.log10(float(coef.D_sampled(cond, coef.sample(rng, "independent"))))
             for _ in range(400)]
    at_T = [np.log10(float(coef.D_sampled(cond, coef.sample(rng, "logD_at_T"))))
            for _ in range(400)]
    assert np.std(indep) > 3.0 * np.std(at_T)


def test_contributions_rank_temperature_first_for_opx():
    m = make_model(comp_dependent=False)
    x = np.linspace(-40, 40, 41)
    C = m.profile(2.0 * SEC_PER_YEAR, x)
    sigma = np.full_like(x, 0.003)
    budget = UncertaintyBudget(sigma_T_K=25.0, buffer="NNO", delta_buffer=1.0,
                               sigma_delta_buffer=0.3, sample_coefficient=False,
                               sample_measurement_noise=False)
    con = contributions(m, x, C, sigma, budget=budget, n_draws=40, seed=9)
    assert con["temperature"] > con["fo2"], (
        "for orthopyroxene the fO2 exponent is only 0.053, so temperature must dominate")


def test_fit_reports_the_solver_route_and_caveats():
    m = make_model(comp_dependent=True)
    x = np.linspace(-40, 40, 41)
    C = m.profile(SEC_PER_YEAR, x)
    r = fit_time(m, x, C)
    assert "Crank-Nicolson" in r.route
    assert any("1-D modelling" in w for w in r.warnings)
