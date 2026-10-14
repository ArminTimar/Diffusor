"""Laws in the form published application studies printed and used (coefficients/study_laws.py)."""
import numpy as np
import pytest

from diffusor.coefficients import Conditions, get
from diffusor.coefficients import study_laws
from diffusor.coefficients.registry import REGISTRY
from diffusor.references import REFERENCES


def test_every_study_law_is_traceable_and_never_recommended():
    assert len(study_laws.COEFFICIENTS) == 43
    for c in study_laws.COEFFICIENTS:
        assert c.key in REGISTRY and c.citation in REFERENCES
        assert not c.recommended and c.verified, c.key
        assert not c.variant_of or c.variant_of in REGISTRY, c.key


def test_tomiya_al_and_mg_assignments_reproduce_the_printed_d():
    """Tomiya et al. (2013) p. 14-15: 950 C, log fO2 = -11, X_Usp = 0.3: Ti (and Al) 4.3e-16, Fe (and Mg) 4.4e-15."""
    cond = Conditions(T_K=1223.15, log_fo2_bar=-11.0, X={"xTi": 0.1})
    assert get("mt_Al_tomiya2013").D(cond) == pytest.approx(4.3e-16, rel=0.05)
    assert get("mt_Mg_tomiya2013").D(cond) == pytest.approx(4.4e-15, rel=0.06)
    assert get("mt_Al_tomiya2013").D(cond) == get("mt_Ti_vanorman_crispin2010").D(cond)


@pytest.mark.parametrize("T, dnno, x_usp, printed", [
    (1073, 0.21, 0.55, 2.29e-16), (1029, 0.40, 0.47, 1.09e-16), (1015, 0.43, 0.46, 8.51e-17),
    (1013, 0.61, 0.42, 8.54e-17), (966, 0.45, 0.44, 3.6e-17), (903, 0.21, 0.40, 8.6e-18)])
def test_morgado_eq5_reproduces_the_printed_d(T, dnno, x_usp, printed):
    """Morgado et al. (2019) Tables 2-3, with the NNO expression of their SM2 workbook (atm)."""
    log_atm = 9.36 - 24930 / (T + 273) + dnno
    cond = Conditions(T_K=T + 273.15, log_fo2_bar=log_atm + np.log10(1.01325), X={"xTi": x_usp / 3})
    assert get("mt_FeTi_morgado2019").D(cond) == pytest.approx(printed, rel=0.015)


def test_costa_chakraborty_is_chakraborty_1997_at_its_reference_state():
    c = get("ol_FeMg_costa_chakraborty2004")
    ref = Conditions(T_K=1473.15, log_fo2_bar=-12.0, X={"XFe": 0.14}, axis="c")
    assert c.D(ref) == pytest.approx(get("ol_FeMg_chakraborty1997").D(ref.replace(axis=None)), rel=1e-4)
    # Fo78 adds 3 x 0.08 log units, one log unit of fO2 adds 1/6
    assert c.log10_D(ref.replace(X={"XFe": 0.22})) - c.log10_D(ref) == pytest.approx(0.24)
    assert c.log10_D(ref.replace(log_fo2_bar=-11.0)) - c.log10_D(ref) == pytest.approx(1 / 6)
    assert c.D(ref.replace(axis="a")) == pytest.approx(c.D(ref) / 6)


def test_mourey_composition_term_has_the_printed_sign():
    cond = Conditions(T_K=1456.15, P_Pa=60e6, log_fo2_bar=-8.37, X={"XFe": 0.1})
    tamed, mourey = get("ol_FeMg_dohmen_chakraborty2007_tamed"), get("ol_FeMg_mourey2023")
    assert mourey.D(cond) == pytest.approx(tamed.D(cond))
    shifted = cond.replace(X={"XFe": 0.2})
    assert mourey.log10_D(shifted) - tamed.log10_D(shifted) == pytest.approx(-0.6)


def test_mutch_regressions_read_columns_e_v_d_as_q_h_j():
    """At 1230 C, 0.8 GPa and QFM the TaMED regression lies within 0.4 log units of Dohmen & Chakraborty."""
    from diffusor.thermo.buffers import log_fo2_from_delta
    T, P = 1503.15, 0.8e9
    cond = Conditions(T_K=T, P_Pa=P, log_fo2_bar=log_fo2_from_delta("QFM", 0, T, P), X={"XFe": 0.12})
    for key in ("ol_FeMg_mutch2019_tamed", "ol_FeMg_mutch2019_global"):
        diff = get(key).log10_D(cond) - get("ol_FeMg_dohmen_chakraborty2007_tamed").log10_D(cond)
        assert abs(diff) < 0.4, key
    for c in study_laws.COEFFICIENTS:
        if c.covariance is not None:
            cov = np.asarray(c.covariance)
            assert np.allclose(cov, cov.T) and np.all(np.diag(cov) > 0)
            assert len(c.cov_order) == cov.shape[0] and set(c.cov_order) <= set(c.params)
    asio2 = get("ol_Ni_mutch2019_asio2")
    assert asio2.params["aSiO2"].value == 0.62


def test_dipra_laws_follow_table_1():
    cond = Conditions(T_K=1423.15, P_Pa=1e5, log_fo2_bar=-12.0, X={"XFe": 0.1}, angles_deg=(90, 90, 0))
    # Ca: the Coogan constants with fO2^0.3 instead of 0.31; equal at the reference 1e-12 bar
    assert get("ol_Ca_girona_costa2013").D(cond) == pytest.approx(get("ol_Ca_coogan2005").D(cond), rel=0.01)
    # Mn: the Fe-Mg TaMED law itself
    assert get("ol_Mn_girona_costa2013").D(cond.replace(angles_deg=None)) == pytest.approx(
        get("ol_FeMg_dohmen_chakraborty2007_tamed").D(cond.replace(angles_deg=None)))
    # Ni: 3.84e-9 exp(-220 kJ/RT) at fO2 = 1e-6 Pa, Fo90 and 1 bar
    ni = get("ol_Ni_girona_costa2013").D(Conditions(T_K=1423.15, log_fo2_bar=-11.0, X={"XFe": 0.1}))
    assert ni == pytest.approx(3.84e-9 * np.exp(-220e3 / (8.314462618 * 1423.15)), rel=1e-6)


def test_pyroxene_study_forms():
    T, lf = 1298.15, -10.0
    sato, dohmen = get("opx_FeMg_sato2022"), get("opx_FeMg_dohmen2016")
    c = Conditions(T_K=T, log_fo2_bar=lf, X={"XFe": 0.09}, axis="c")
    assert sato.D(c) == pytest.approx(dohmen.D(c))
    assert sato.D(c.replace(axis="b")) == pytest.approx(sato.D(c) / 3.5)
    assert sato.D(c.replace(X={"XFe": 0.3})) == sato.D(c)            # no composition term
    winslow, dias = get("opx_FeMg_winslow2026"), get("opx_FeMg_dias_dohmen2024")
    m = -2.711e4 / T + 23.5408
    w = Conditions(T_K=T, X={"XFe": 0.25}, axis="c")
    assert winslow.log10_D(w) - dias.log10_D(w) == pytest.approx(0.01 * m)
    ban, muller = get("cpx_FeMg_ban2026"), get("cpx_FeMg_muller2013")
    assert ban.D(Conditions(T_K=T, axis="c")) == pytest.approx(muller.D(Conditions(T_K=T)))
    assert ban.D(Conditions(T_K=T, axis="a")) == pytest.approx(muller.D(Conditions(T_K=T)) / 3.5)
    dw = get("cpx_FeMg_dimanov_wiedenbeck2006")
    assert dw.D(Conditions(T_K=1373.15, log_fo2_bar=-9.0)) > dw.D(Conditions(T_K=1373.15, log_fo2_bar=-12.0))


def test_druitt_sr_carries_the_printed_factor():
    c = Conditions(T_K=1150.0, X={"XAn": 0.5})
    assert get("plag_Sr_druitt2012").D(c) == pytest.approx(2.92 * get("plag_Sr_giletti_casserly1994").D(c), rel=1e-4)


def test_dfens_regressions_follow_the_2021_supporting_information():
    """Mutch et al. (2021): version 1 is the 2019 regression with more digits, Ni has a 10.7 anisotropy,
    version 2 drops the P/T term, and the plagioclase rows follow the calibrations they were fitted to."""
    cond = Conditions(T_K=1463.15, P_Pa=0.35e9, log_fo2_bar=-8.0, X={"XFe": 0.15, "XAn": 0.6}, axis="c")
    for sp in ("FeMg_{}_global", "FeMg_{}_tamed", "Ni_{}", "Mn_{}"):
        v1, v2 = get("ol_" + sp.format("mutch2021")), get("ol_" + sp.format("mutch2021") + "_v2")
        v19 = get("ol_" + sp.format("mutch2019"))
        assert abs(v1.log10_D(cond) - v19.log10_D(cond)) < 0.03 and abs(v2.log10_D(cond) - v1.log10_D(cond)) < 0.03
        assert "V" not in v2.params and v2.cov_order == ("A", "B", "C", "E", "D")
        np.testing.assert_allclose(v1.covariance, v19.covariance)
    ni = get("ol_Ni_mutch2021")
    assert ni.D(cond.replace(axis="a")) == pytest.approx(ni.D(cond) / 10.7)
    assert get("ol_FeMg_mutch2021_global").D(cond.replace(axis="b")) == pytest.approx(
        get("ol_FeMg_mutch2021_global").D(cond) / 6)
    plag = Conditions(T_K=1373.15, X={"XAn": 0.6})
    for dfens, calib in (("plag_Mg_mutch2021_vo14", "plag_Mg_vanorman2014"),
                         ("plag_Sr_mutch2021_cw", "plag_Sr_cherniak_watson1994"),
                         ("plag_Sr_mutch2021_gc", "plag_Sr_giletti_casserly1994"),
                         ("plag_Ba_mutch2021", "plag_Ba_cherniak2002")):
        assert abs(get(dfens).log10_D(plag) - get(calib).log10_D(plag)) < 0.1, dfens
    f13 = get("plag_Mg_mutch2021_f13")
    assert f13.D(plag) == pytest.approx(np.exp(-11.77 - 34080 / 1373.15))     # aSiO2 = 1 by default
    assert f13.D(plag, {"aSiO2": 0.6}) == pytest.approx(f13.D(plag) * 0.6 ** 2.931)
    from diffusor.coefficients.plagioclase import activity_A
    assert activity_A("Mg", 1463.15, "mutch2021") == -34.099


def test_behrens_water_in_melt():
    """Behrens et al. (2004): eq. 7 gives the 36 um2/s printed for dacite at 1608 K, D is proportional
    to C_H2Ot, and eq. 10 lies on the 1400 C line of Fig. 5 for dacite at 1 GPa (about -10.25)."""
    from diffusor.minerals import get_mineral
    dac = get("melt_H2O_behrens2004_dacite")
    one = Conditions(T_K=1608.0, P_Pa=1e9, X={"H2Ot": 1.0})
    assert dac.D(one) * 1e12 == pytest.approx(36.0, rel=0.02)
    assert dac.D(one.replace(X={"H2Ot": 4.0})) == pytest.approx(4 * dac.D(one))
    assert dac.D(one.replace(X={"H2Ot": 0.0})) == 0.0
    eq10 = get("melt_H2O_behrens2004_sio2")
    assert eq10.log10_D(Conditions(T_K=1673.15, P_Pa=1e9, X={"H2Ot": 1.0})) == pytest.approx(-10.25, abs=0.03)
    assert get_mineral("melt").composition_variable.key == "H2Ot"
    assert {c.species for c in study_laws.COEFFICIENTS if c.mineral == "melt"} == set()


def test_jollands_ti_in_quartz_and_vielzeuf_ca_in_garnet():
    """Jollands et al. (2020): log D = -8.3 - 311 kJ/(2.303 RT) along c. Vielzeuf et al. (2007) Table 5:
    D at 10 kbar of 6.96e-24, 1.56e-21 and 4.45e-20 m2/s at 700, 1000 and 1300 C."""
    q = get("qz_Ti_jollands2020")
    assert q.log10_D(Conditions(T_K=1473.15, axis="c")) == pytest.approx(-8.3 - 311e3 / (np.log(10) * 8.314462618 * 1473.15))
    assert get("qz_Ti_cherniak2007").log10_D(Conditions(T_K=1073.15, axis="c")) - q.log10_D(
        Conditions(T_K=1073.15, axis="c")) == pytest.approx(3.0, abs=0.3)      # "two to three orders" lower
    ca = get("grt_Ca_vielzeuf2007")
    for T_C, printed in ((700, 6.96e-24), (1000, 1.56e-21), (1300, 4.45e-20)):
        assert ca.D(Conditions(T_K=T_C + 273.15, P_Pa=1e9)) == pytest.approx(printed, rel=0.005)
