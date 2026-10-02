"""The anorthite activity factors A_i of plagioclase: two published sets, chosen per run."""
import numpy as np
import pytest

from diffusor.coefficients import Conditions, get
from diffusor.coefficients.plagioclase import (ACTIVITY_SETS, activity_A, activity_note,
                                               equilibrium_profile)
from diffusor.constants import R_GAS
from diffusor.fitting import DiffusionModel
from diffusor.solvers.boundary import dirichlet, zero_flux
from diffusor.solvers.history import ThermalHistory
from diffusor.solvers.initial import InitialCondition


def test_values_are_dohmen_et_al_2017_table_1():
    """Dohmen, Faak & Blundy (2017) Table 1, p. 556 (rendered table)."""
    db = ACTIVITY_SETS["dohmen_blundy2014"]["columns"]
    assert db[1473.15] == {"Mg": 15.8, "Sr": -17.4, "Ba": -35.1, "Li": -1.7, "K": -8.0, "Rb": -15.9}
    assert db[1173.15] == {"Mg": 13.7, "Sr": -15.1, "Ba": -30.5, "Li": -2.5, "K": -8.8, "Rb": -16.7}
    bi = ACTIVITY_SETS["bindeman1998"]["columns"][None]
    assert bi == {"Mg": -26.1, "Sr": -30.4, "Ba": -55.0, "Li": -6.9, "K": -25.5, "Rb": -40.0}
    # Costa et al. (2003) quote the same A_Mg = -26.1 kJ/mol from Bindeman et al. (1998)
    assert activity_A("Mg", 1173.15, "bindeman1998") == -26.1


def test_dohmen_blundy_uses_the_nearer_column_without_interpolating():
    assert activity_A("Mg", 900 + 273.15) == 13.7
    assert activity_A("Mg", 1000 + 273.15) == 13.7
    assert activity_A("Mg", 1100 + 273.15) == 15.8
    assert "900 C column" in activity_note("dohmen_blundy2014", 950 + 273.15)
    with pytest.raises(KeyError):
        activity_A("Mg", 1173.15, "no_such_set")


def test_the_sets_give_opposite_equilibrium_slopes_for_mg():
    an = np.array([0.4, 0.8])
    up = equilibrium_profile(an, 1173.15, "Mg", 100.0, 0.4, "dohmen_blundy2014")
    down = equilibrium_profile(an, 1173.15, "Mg", 100.0, 0.4, "bindeman1998")
    assert up[1] > up[0] and down[1] < down[0]


def _model(history=None, A=-26.1):
    x = np.linspace(0, 200, 201)
    an = np.linspace(0.4, 0.8, x.size)
    T = 1173.15
    return DiffusionModel(get("plag_Mg_costa2003"), Conditions(T, X={"XAn": an}),
                          InitialCondition("table", {"x_table": x, "C_table": np.full(x.size, 150.0)}),
                          bc_left=dirichlet(150.0), bc_right=zero_flux(), x_grid=x, boundaries_far=False,
                          an_profile=an, activity_A_kJ=A, activity_set="bindeman1998", history=history)


def test_theta_follows_the_temperature_of_each_step():
    m = _model()
    assert m.activity_theta == pytest.approx(-26.1e3 / (R_GAS * 1173.15))
    x = np.linspace(0, 200, 41)
    yr = 3.15576e7
    iso = m.profile(30 * yr, x)
    # a path that ends at the nominal temperature but starts 100 K hotter differs from the
    # isothermal run in D *and* in A/RT; with A/RT frozen at the nominal T it would differ
    # only through D
    hot = _model(ThermalHistory.linear(1273.15, 1173.15, 1.0)).profile(30 * yr, x)
    assert not np.allclose(iso, hot)
    th = _model()._theta()
    assert th(1273.15) == pytest.approx(-26.1e3 / (R_GAS * 1273.15))


def test_the_set_is_recorded_in_the_methods_and_citations():
    from diffusor.dataio.export import collect_citations, methods_paragraph
    from diffusor.fitting.fit import FitResult
    from diffusor.fitting.objective import statistics
    m = _model()
    x = np.linspace(0, 200, 11)
    C = m.profile(3e8, x)
    fit = FitResult(3e8, m, x, C, None, C, statistics(C, C, None))
    assert "bindeman1998" in collect_citations(fit)
    assert "A = -26.1 kJ/mol" in methods_paragraph(fit)
    assert "Bindeman et al. (1998)" in methods_paragraph(fit)
