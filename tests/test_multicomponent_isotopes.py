"""Tracer versus exchange coefficients, multicomponent garnet diffusion and isotope fractionation."""
import json

import numpy as np
import pytest

from diffusor.coefficients import Conditions, get
from diffusor.coefficients.base import TRANSPORT_KINDS, DiffusionCoefficient, TRACER_ADVICE
from diffusor.coefficients.families import FAMILIES, get_family
from diffusor.coefficients.garnet import A0_NM, CARLSON_TABLE4, unit_cell_edge
from diffusor.coefficients.registry import REGISTRY
from diffusor.coefficients.transport import binary_interdiffusion, ideal_ionic_matrix
from diffusor.constants import R_GAS
from diffusor.fitting import DiffusionModel
from diffusor.fitting.isotopes import (CoupledFeMgIsotopeModel, DiluteIsotopeModel, fit_coupled_isotopes,
                                       fit_dilute_isotopes)
from diffusor.fitting.multicomponent import MulticomponentModel, fit_multicomponent_time
from diffusor.solvers.initial import InitialCondition, step
from diffusor.solvers.multicomponent import solve_multicomponent, step_constant_matrix
from diffusor.solvers.numerical import solve_1d
import diffusor.coefficients.isotopes as iso


# ---------------------------------------------------------------- transport kinds
def test_every_law_states_a_known_transport_kind():
    assert all(c.kind in TRANSPORT_KINDS for c in REGISTRY.values())
    with pytest.raises(ValueError, match="transport kind"):
        DiffusionCoefficient(key="x", mineral="m", species="s", label="l", citation="crank1975",
                             equation_text="", func=lambda *a: 1.0)
    # the sources' own descriptions
    assert REGISTRY["opx_Mg_schwandt1998_c"].kind == "tracer"
    assert REGISTRY["plag_Sr_giletti_casserly1994"].kind == "tracer"
    assert REGISTRY["opx_FeMg_dohmen2016"].kind == "interdiffusion"


def test_tracer_laws_carry_the_exchange_advice():
    m = DiffusionModel(get("opx_Mg_schwandt1998_c"), Conditions(1100.0, axis="c"),
                       InitialCondition("step", {"x0": 0.0, "C_left": 1.0, "C_right": 0.0}))
    assert TRACER_ADVICE in m.warnings()
    m = DiffusionModel(get("opx_FeMg_dohmen2016"), Conditions(1100.0, log_fo2_bar=-12, X={"XFe": .2}),
                       InitialCondition("step", {"x0": 0.0, "C_left": .2, "C_right": .3}))
    assert TRACER_ADVICE not in m.warnings()


def test_binary_interdiffusion_limits():
    assert binary_interdiffusion(2e-18, 2e-18, .3) == pytest.approx(2e-18)
    # the dilute species controls the exchange
    assert binary_interdiffusion(1e-18, 1e-16, 1e-9) == pytest.approx(1e-18, rel=1e-6)
    assert binary_interdiffusion(1e-18, 1e-16, 1 - 1e-9) == pytest.approx(1e-16, rel=1e-6)
    # symmetric in the labels, and between the two tracer values
    D = binary_interdiffusion(1e-18, 5e-17, .4)
    assert D == pytest.approx(binary_interdiffusion(5e-17, 1e-18, .6))
    assert 1e-18 < D < 5e-17
    with pytest.raises(ValueError):
        binary_interdiffusion(-1, 1, .5)


def test_ideal_ionic_matrix_limits():
    # binary: the 1x1 matrix is the binary interdiffusion coefficient
    M = ideal_ionic_matrix(np.array([3.0, 0.5]), np.array([.3, .7]))
    assert M[0, 0] == pytest.approx(binary_interdiffusion(3.0, .5, .3))
    # equal tracers: no coupling
    M = ideal_ionic_matrix(np.full(4, 2.0), np.array([.1, .2, .3, .4]))
    assert np.allclose(M, 2.0 * np.eye(3))


def test_chakraborty_ganguly_1992_eq5_matrix_is_reproduced():
    """Eq. 5 prints nine matrix elements (cm2/s) at X_Fe .61, X_Mn .20, X_Mg .18, X_Ca .01 with
    Ca dependent. This verifies the structure of eq. 2 element by element, not the tracer
    coefficients. The paper does not print the D* behind eq. 5 (its text points to Table 1, which
    holds microprobe analyses, and Table 2 covers only 1200 C and 20 to 35 kb). Only D*Ca = 0.5 D*Fe
    (p. 81) is stated. The values D*Mn 2.495e-12, D*Mg = D*Fe = 6.0e-13 cm2/s were chosen to
    reproduce the printed matrix: they reproduce all nine elements within 0.45 %. Eq. 3 at the
    experiment's 41 kb and 1430 C gives D*Mn 1.52e-12, D*Mg 4.46e-13 and D*Fe 4.52e-13 cm2/s, and
    a matrix 17 to 36 % below the printed one. D[Fe,Mn] (-8.231e-13 against -8.24e-13) and
    D[Fe,Mg] (-1.125e-13 against -1.12e-13) differ in the last printed digit, which is rounding.
    The comparison has no absolute tolerance: with numpy's default atol of 1e-8 every element,
    all of order 1e-12 or smaller, would pass whatever its value."""
    printed = np.array([[1.37e-12, -1.53e-13, -1.53e-13],
                        [-2.43e-13, 5.67e-13, -3.32e-14],
                        [-8.24e-13, -1.12e-13, 4.88e-13]])
    D = np.array([2.495e-12, 6.0e-13, 6.0e-13, 3.0e-13])       # Mn, Mg, Fe, Ca
    M = ideal_ionic_matrix(D, np.array([.20, .18, .61, .01]), dependent=-1)
    assert np.allclose(M, printed, rtol=5e-3, atol=0.0)
    assert not np.allclose(M, 0.7 * printed, rtol=5e-3, atol=0.0)      # the comparison is not vacuous


# ---------------------------------------------------------------- new scalar laws
def test_oeser2026_tracer_laws_and_derived_interdiffusion():
    fe, mg, femg = get("ol_Fe_oeser2026_tracer"), get("ol_Mg_oeser2026_tracer"), get("ol_FeMg_oeser2026")
    R = 8.314462618
    table = {"a": ((287, -6.29), (231, -8.69), (314, -5.37)), "b": ((307, -5.22), (292, -6.11), (363, -3.22)),
             "c": ((328, -4.18), (329, -4.39), (371, -2.55))}
    for T in (1100.0, 1150.0, 1250.0):
        for ax, (tf, tm, tfm) in table.items():
            cond = Conditions(T + 273.15, X={"XFe": .085}, axis=ax)
            log = lambda q, l: l - q * 1e3 / (np.log(10) * R * cond.T_K)
            assert np.log10(fe.D(cond)) == pytest.approx(log(*tf), abs=1e-6)
            assert np.log10(mg.D(cond)) == pytest.approx(log(*tm), abs=1e-6)
            # eq. 6 applied to the tracer fits reproduces the paper's own D_Fe-Mg fit (Table 4)
            assert np.log10(femg.D(cond)) == pytest.approx(log(*tfm), abs=0.2)
    assert femg.kind == "interdiffusion" and femg.derived_from == (fe.key, mg.key)
    with pytest.raises(ValueError):
        femg.D(Conditions(1450.0, X={"XFe": .085}))          # orientation is required
    # composition term 10^(3 (X_Fe - 0.085))
    c1, c2 = Conditions(1400.0, X={"XFe": .085}, axis="c"), Conditions(1400.0, X={"XFe": .185}, axis="c")
    assert fe.D(c2) / fe.D(c1) == pytest.approx(10 ** 0.3)


def test_schaffer2014_na_k_exchange_laws():
    c = get("kfs_NaK_schaffer2014_001_or92")
    assert c.D(Conditions(1273.15)) == pytest.approx(5.18e-8 * np.exp(-179e3 / (R_GAS * 1273.15)))
    with pytest.raises(ValueError):
        c.D(Conditions(1273.15, axis="c"))        # the normal to (001) is not the c axis
    b = get("kfs_NaK_schaffer2014_010_or92")
    assert b.D(Conditions(1273.15, axis="b")) == pytest.approx(5.81e-5 * np.exp(-272e3 / (R_GAS * 1273.15)))
    # normal to (001) about one order of magnitude faster than normal to (010) (abstract)
    assert 3 < c.D(Conditions(1173.15)) / b.D(Conditions(1173.15)) < 30


def test_borinski_garnet_binary_from_tracers():
    fe, mg, femg = get("grt_Fe_borinski2012"), get("grt_Mg_borinski2012"), get("grt_FeMg_borinski2012")
    cond = Conditions(1573.15, P_Pa=2.5e9, X={"XFe": .6})
    assert fe.D(cond) == pytest.approx(1.64e-10 * np.exp(-(226.9e3 + .56 * (25000 - 1)) / (R_GAS * 1573.15)))
    assert femg.D(cond) == pytest.approx(binary_interdiffusion(fe.D(cond), mg.D(cond), .6))


# ---------------------------------------------------------------- garnet families
def test_carlson2006_transcription():
    fam = get_family("grt_carlson2006")
    alm = {"Fe": 1.0, "Mg": 0.0, "Mn": 0.0, "Ca": 0.0}
    D = fam.tracer(1273.15, 1e9, alm)
    for c, (lnD0, k, Q, dV) in CARLSON_TABLE4.items():
        assert np.log(D[c]) == pytest.approx(lnD0 - (Q + dV) * 1e3 / (R_GAS * 1273.15))
    up = fam.tracer(1273.15, 1e9, alm, {"dlogfo2_graphite": 6.0})
    assert np.log10(up["Fe"] / D["Fe"]) == pytest.approx(1.0)
    X = {"Fe": .05, "Mg": .0, "Mn": .0, "Ca": .95}         # a0 = 1.1836 nm
    assert unit_cell_edge(X) == pytest.approx(.05 * A0_NM["Fe"] + .95 * A0_NM["Ca"])
    assert any("a0 < 1.1821" in w for w in fam.check(1000, 1e9, X))
    assert not any("a0 < 1.1821" in w for w in fam.check(1000, 1e9, alm))


def test_chen_chu_mn_law_and_chakraborty_ganguly_ca():
    # Chen & Chu (2024): D_Mn about 1e-24 m2/s at about 510 C and 2 GPa for their eclogite garnet
    D = get_family("grt_carlson2006_chen_chu2024_mn").tracer(783.15, 2e9, {"Fe": .518, "Mg": .068, "Mn": .05, "Ca": .364})
    assert np.log10(D["Mn"]) == pytest.approx(-24, abs=0.5)
    D = get_family("grt_chakraborty_ganguly1992").tracer(1473.15, 2.06e9, {"Fe": .5, "Mg": .1, "Mn": .3, "Ca": .1})
    assert D["Ca"] == pytest.approx(0.5 * D["Fe"])
    assert len(FAMILIES) == 4


# ---------------------------------------------------------------- vector solver
def test_constant_matrix_matches_eigen_solution_and_conserves_mass():
    D = ideal_ionic_matrix(np.array([2.5, .6, .6, .3]), np.array([.4, .14, .43, .03]))
    CL, CR = np.array([.2, .18, .61]), np.array([.6, .1, .25])
    errs = []
    for n in (201, 401):
        x = np.linspace(-50, 50, n)
        C0 = np.array([step(x, 0.0, a, b) for a, b in zip(CL, CR)])
        r = solve_multicomponent(x, C0, lambda C, T: np.repeat(D[:, :, None], n, 2), 100.0, T_K=1000.0)
        assert np.allclose(r.mass_initial, r.mass_final, atol=1e-9)
        errs.append(np.abs(r.C_final - step_constant_matrix(x, 0.0, CL, CR, D, 100.0)).max())
    assert errs[1] < 2e-3 and errs[1] < errs[0]


def test_two_components_reduce_to_the_scalar_solver():
    x = np.linspace(-30, 30, 121)
    C0 = step(x, 0.0, .8, .2)
    DA, DB = 2.0, .1

    def D_bin(C):
        return binary_interdiffusion(DA, DB, np.clip(C, 0, 1))
    scalar = solve_1d(x, C0, lambda C, T: D_bin(C), 50.0, T_K=1000.0).C_final
    vector = solve_multicomponent(x, C0[None], lambda C, T: ideal_ionic_matrix(
        np.array([np.full(x.size, DA), np.full(x.size, DB)]), np.vstack([C, 1 - C]))[..., :], 50.0,
        T_K=1000.0, min_steps=400).C_final[0]
    assert np.abs(scalar - vector).max() < 2e-3


def _garnet_model(dependent="Ca", **kw):
    L = {"Fe": .70, "Mg": .20, "Mn": .05, "Ca": .05}
    R = {"Fe": .60, "Mg": .10, "Mn": .25, "Ca": .05}
    ini = {c: InitialCondition("step", {"x0": 0.0, "C_left": L[c], "C_right": R[c]}) for c in L}
    return MulticomponentModel(get_family("grt_carlson2006"), 1073.15, ini, P_Pa=.8e9, dependent=dependent,
                               n_nodes=161, **kw)


def test_dependent_component_choice_does_not_change_profiles_and_ca_moves_uphill():
    x = np.linspace(-40, 40, 41)
    a = _garnet_model("Ca").profiles(2e11, x)
    b = _garnet_model("Mn").profiles(2e11, x)
    for c in a:
        assert np.allclose(a[c], b[c], atol=2e-4)
    assert np.allclose(sum(a.values()), 1.0)
    # Ca starts flat; the off-diagonal terms move it (Carlson 2006, Fig. 4)
    assert np.ptp(a["Ca"]) > 5e-4


def test_garnet_fit_recovers_a_synthetic_duration():
    x = np.linspace(-40, 40, 41)
    m = _garnet_model()
    data = m.profiles(2e11, x)
    r = fit_multicomponent_time(m, x, data, {c: .003 for c in data}, scan_points=24)
    assert r.t_seconds == pytest.approx(2e11, rel=0.02)
    assert r.stats.chi2 < 1e-3
    with pytest.raises(InterruptedError):
        fit_multicomponent_time(m, x, data, progress=lambda f: True)


# ---------------------------------------------------------------- isotopes
def test_delta_round_trip_and_reference_insensitivity(monkeypatch):
    for el in iso.ABUNDANCES:
        parts = iso.split_isotopes(el, 2.0, -7.5)
        assert sum(parts.values()) == pytest.approx(2.0)
        assert iso.delta_value(el, parts) == pytest.approx(-7.5)
    c = get("plag_Li_pohl2024_interstitial")
    dm = DiffusionModel(c, Conditions(1173.15, X={"XAn": .6}),
                        InitialCondition("step", {"x0": 150.0, "C_left": 2.0, "C_right": 20.0}),
                        composition_dependent=False)
    x = np.linspace(0, 200, 41)
    _, d1 = DiluteIsotopeModel(dm, "Li", .27).profiles(600, x)
    monkeypatch.setitem(iso.ABUNDANCES, "Li", {6: .0759 * 1.01, 7: 1 - .0759 * 1.01})
    _, d2 = DiluteIsotopeModel(dm, "Li", .27).profiles(600, x)
    assert np.abs(d1 - d2).max() < 0.01 * np.abs(d1).max()


def test_dilute_isotopes_beta_zero_and_recovery():
    c = get("plag_Li_pohl2024_interstitial")
    dm = DiffusionModel(c, Conditions(1173.15, X={"XAn": .6}),
                        InitialCondition("step", {"x0": 150.0, "C_left": 2.0, "C_right": 20.0}),
                        composition_dependent=False)
    x = np.linspace(0, 200, 61)
    C0, d0 = DiluteIsotopeModel(dm, "Li", 0.0).profiles(600, x)
    assert np.allclose(d0, 0.0, atol=1e-9)
    assert np.allclose(C0, dm.profile(600, x))
    model = DiluteIsotopeModel(dm, "Li", .27)
    C, d = model.profiles(600, x)
    assert d.min() < -10                   # light 6Li runs ahead of 7Li into the crystal
    r = fit_dilute_isotopes(model, x, C, .1, d, .5, scan_points=30)
    assert r.t_seconds == pytest.approx(600, rel=1e-3)
    r = fit_dilute_isotopes(DiluteIsotopeModel(dm, "Li", .2), x, C, .1, d, .5, fit_beta=True, scan_points=30)
    assert r.beta == pytest.approx(.27, abs=1e-3)
    with pytest.raises(ValueError, match="uncertainties"):
        fit_dilute_isotopes(model, x, C, None, d, .5)


def _coupled(beta_fe=.22, beta_mg=.17, **kw):
    cond = Conditions(1473.15, axis="b", X={})
    return CoupledFeMgIsotopeModel(cond, InitialCondition("step", {"x0": 0.0, "C_left": .085, "C_right": .25}),
                                   beta_fe, beta_mg, right_fixed=True, **kw)


def test_coupled_isotopes_reduce_to_the_scalar_fe_mg_law():
    x = np.linspace(-40, 20, 61)
    t = 2 * 86400.0
    m = _coupled(0.0, 0.0, tracer_fe=get("ol_Fe_oeser2026_tracer"), tracer_mg=get("ol_Mg_oeser2026_tracer"))
    p = m.profiles(t, x)
    assert np.allclose(p["d56Fe"], 0, atol=1e-9) and np.allclose(p["d26Mg"], 0, atol=1e-9)
    from diffusor.solvers.boundary import dirichlet
    scalar = DiffusionModel(get("ol_FeMg_oeser2026"), Conditions(1473.15, axis="b", X={"XFe": .1}),
                            InitialCondition("step", {"x0": 0.0, "C_left": .085, "C_right": .25}),
                            bc_right=dirichlet(.25), comp_key="XFe", boundaries_far=False, n_nodes=201,
                            x_grid=np.linspace(-40, 20, 201))
    assert np.abs(p["XFe"] - scalar.profile(t, x)).max() < 3e-3


def test_coupled_isotope_signs_ratio_mode_and_fit():
    x = np.linspace(-40, 20, 61)
    t = 2 * 86400.0
    m = _coupled(tracer_fe=get("ol_Fe_oeser2026_tracer"), tracer_mg=get("ol_Mg_oeser2026_tracer"))
    p = m.profiles(t, x)
    inside = x < 0
    assert p["d56Fe"][inside].min() < -0.2       # light Fe arrives first where Fe diffuses in
    assert p["d26Mg"][inside].max() > 0.05       # light Mg leaves first, heavy Mg stays
    r = fit_coupled_isotopes(m, x, {k: p[k] for k in ("XFe", "d56Fe", "d26Mg")},
                             {"XFe": .003, "d56Fe": .05, "d26Mg": .05}, scan_points=20)
    assert r.t_seconds == pytest.approx(t, rel=0.02)
    # an interdiffusion law with a ratio satisfies eq. 6 at every composition
    rm = _coupled(interdiffusion=get("ol_FeMg_oeser2026"), ratio_fe_mg=2.5)
    XFe = np.array([.1, .2, .3])
    dfe, dmg = rm.element_tracers(XFe, 1473.15)
    target = get("ol_FeMg_oeser2026").D(Conditions(1473.15, axis="b", X={"XFe": XFe}))
    assert np.allclose(binary_interdiffusion(dfe, dmg, XFe), target)
    assert np.allclose(dfe / dmg, 2.5)
    assert rm.reference_masses == (56, 24) and m.reference_masses == (57, 25)
    with pytest.raises(ValueError):
        _coupled(interdiffusion=get("ol_FeMg_oeser2026"))     # ratio missing


def test_beta_registry_names_its_model():
    from diffusor.coefficients.isotopes import BETAS, betas_for
    assert {b.direction for b in betas_for("Fe", "olivine") if b.citation == "oeser2026"} == {"a", "b", "c"}
    assert BETAS["ol_Fe_oeser2026_c"].upper_bound and BETAS["ol_Fe_teng2011"].assumed
    assert all(b.diffusion_model for b in BETAS.values())


# ---------------------------------------------------------------- export and interface
def test_exports_write_results(tmp_path):
    from diffusor.dataio.multicomponent import save_isotope_results, save_multicomponent_results
    x = np.linspace(-40, 40, 21)
    m = _garnet_model()
    data = m.profiles(2e11, x)
    r = fit_multicomponent_time(m, x, data, {c: .003 for c in data}, scan_points=12)
    save_multicomponent_results(tmp_path / "grt", r)
    j = json.loads((tmp_path / "grt" / "multicomponent_results.json").read_text(encoding="utf-8"))
    assert j["family"] == "grt_carlson2006" and j["t_seconds"] == pytest.approx(r.t_seconds)
    assert (tmp_path / "grt" / "multicomponent_results.xlsx").exists()
    assert "Lasaga" in (tmp_path / "grt" / "methods.txt").read_text(encoding="utf-8")
    c = get("plag_Li_pohl2024_interstitial")
    dm = DiffusionModel(c, Conditions(1173.15, X={"XAn": .6}),
                        InitialCondition("step", {"x0": 150.0, "C_left": 2.0, "C_right": 20.0}),
                        composition_dependent=False)
    xi = np.linspace(0, 200, 31)
    model = DiluteIsotopeModel(dm, "Li", .27)
    C, d = model.profiles(600, xi)
    ri = fit_dilute_isotopes(model, xi, C, .1, d, .5, scan_points=12)
    save_isotope_results(tmp_path / "li", ri, "test")
    assert json.loads((tmp_path / "li" / "isotope_results.json").read_text(encoding="utf-8"))["model"] == "test"


def test_study_dialog_builds(qtbot=None):
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])
    from diffusor.gui.main_window import MainWindow
    w = MainWindow()
    w.show_multicomponent_study()
    dlg = w.multicomponent_study
    assert dlg.tabs.count() == 2
    assert dlg.g_family.count() == 4
    dlg.close()
    w.close()
