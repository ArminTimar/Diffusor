"""Each diffusion coefficient must reproduce a number printed in its source.

These are the tests that keep the registry honest: if a transcription is wrong,
the value published in the paper will not come back out.
"""
import numpy as np
import pytest

from diffusor.coefficients import Conditions, get, list_coefficients
from diffusor.references import REFERENCES


def test_every_coefficient_has_a_known_citation():
    for c in list_coefficients():
        assert c.citation in REFERENCES, f"{c.key} cites unknown key {c.citation}"
        for k in c.secondary_citations:
            assert k in REFERENCES, f"{c.key} cites unknown secondary key {k}"
        assert c.equation_text, f"{c.key} has no transcribed equation"
        assert c.verified_from, f"{c.key} does not say where its numbers came from"


def test_opx_dohmen2016_reproduces_run_OPXD_14():
    """Dohmen et al. (2016) Table 2, run OPXD_14: 950 C, log fO2 = -7 Pa, //c, Fs9.

    Published log D = -19.49 +/- 0.07 (sample 7D_20) and -19.59 +/- 0.07 (7D_21).
    """
    c = get("opx_FeMg_dohmen2016")
    cond = Conditions(T_K=950 + 273.15, log_fo2_bar=-12.0, X={"XFe": 0.09}, axis="c")
    assert np.log10(float(c.D(cond))) == pytest.approx(-19.49, abs=0.1)


def test_opx_dohmen2016_anisotropy_is_a_factor_of_3_5():
    c = get("opx_FeMg_dohmen2016")
    cond = Conditions(T_K=1223.15, log_fo2_bar=-12.0, X={"XFe": 0.09})
    assert float(c.D(cond.replace(axis="c")) / c.D(cond.replace(axis="a"))) == pytest.approx(3.5, rel=1e-9)
    assert float(c.D(cond.replace(axis="b")) / c.D(cond.replace(axis="c"))) == pytest.approx(1.0)


def test_opx_dohmen2016_abstract_form_matches_equation_1():
    """Abstract: D = 1.12e-6 (fO2[Pa])^0.053 exp(-308 kJ/mol / RT)."""
    from diffusor.constants import R_GAS
    c = get("opx_FeMg_dohmen2016")
    T, log_fo2_Pa = 1223.15, -7.0
    expected = 1.12e-6 * (10.0 ** log_fo2_Pa) ** 0.053 * np.exp(-308e3 / (R_GAS * T))
    cond = Conditions(T_K=T, log_fo2_bar=log_fo2_Pa - 5, X={"XFe": 0.09}, axis="c")
    assert float(c.D(cond)) == pytest.approx(expected, rel=0.01)


def test_opx_dohmen2016_fs1_abstract_form():
    """Abstract: D = 1.66e-4 exp(-377 kJ/mol / RT) for XFe = 0.01."""
    from diffusor.constants import R_GAS
    c = get("opx_FeMg_dohmen2016_fs1")
    T = 1273.15
    expected = 1.66e-4 * np.exp(-377e3 / (R_GAS * T))
    assert float(c.D(Conditions(T_K=T, axis="c"))) == pytest.approx(expected, rel=0.02)


def test_opx_composition_term_is_a_decade_per_unit_XFe():
    c = get("opx_FeMg_dohmen2016")
    base = Conditions(T_K=1223.15, log_fo2_bar=-12.0, X={"XFe": 0.09}, axis="c")
    hi = base.replace(X={"XFe": 1.09})
    assert float(c.D(hi) / c.D(base)) == pytest.approx(10.0, rel=1e-9)


def test_ganguly_tazzoli_matches_ostorero2022_equation_1():
    """Ostorero et al. (2022) eq. 1: log D [cm2/s] = -5.54 + 2.6 XFe - 12530/T."""
    c = get("opx_FeMg_ganguly_tazzoli1994")
    T, XFe = 1173.15, 0.30
    expected_cm2 = -5.54 + 2.6 * XFe - 12530.0 / T
    cond = Conditions(T_K=T, log_fo2_bar=-12.0, X={"XFe": XFe})
    got_m2 = np.log10(float(c.D(cond, {"use_fo2": 0.0})))
    assert got_m2 == pytest.approx(expected_cm2 - 4.0, abs=1e-9)


def test_cpx_muller2013_matches_the_abstract():
    from diffusor.constants import R_GAS
    c = get("cpx_FeMg_muller2013")
    T = 1273.15
    expected = 2.77e-7 * np.exp(-320.7e3 / (R_GAS * T))
    assert float(c.D(Conditions(T_K=T))) == pytest.approx(expected, rel=1e-9)


def test_cpx_dimanov_sautter_reproduces_petrone2016_table2():
    """Petrone et al. (2016) Table 2 footnote: D = 3.26e-20 at 1098 C, 1.20e-19 at 1150 C."""
    c = get("cpx_FeMg_dimanov_sautter2000")
    assert float(c.D(Conditions(T_K=1098 + 273.15))) == pytest.approx(3.26e-20, rel=0.005)
    assert float(c.D(Conditions(T_K=1150 + 273.15))) == pytest.approx(1.20e-19, rel=0.005)
    # Dimanov & Sautter (2000) p. 757: log D0 [cm2/s] = -0.02, abstract 0.955 cm2/s
    assert c.params["D0"].value == pytest.approx(0.955e-4, rel=1e-12)
    assert c.params["Q"].value == 406.0


def test_schwandt1998_table3_values():
    from diffusor.constants import R_GAS
    for axis, D0, Ea in (("a", 1.10e-4, 360e3), ("b", 6.93e-6, 339e3), ("c", 4.34e-9, 265e3)):
        c = get(f"opx_Mg_schwandt1998_{axis}")
        T = 1173.15
        assert float(c.D(Conditions(T_K=T))) == pytest.approx(D0 * np.exp(-Ea / (R_GAS * T)), rel=0.02)


def test_plag_Mg_vanorman2014_equation_4():
    from diffusor.constants import R_GAS
    c = get("plag_Mg_vanorman2014")
    T, xan = 1173.15, 0.67
    expected = np.exp(-6.06 - 7.96 * xan - 287e3 / (R_GAS * T))
    assert float(c.D(Conditions(T_K=T, X={"XAn": xan}))) == pytest.approx(expected, rel=1e-9)


def test_plag_Mg_costa2003_is_faster_than_vanorman_at_low_T_low_An():
    """Van Orman et al. (2014) p. 84: the Costa et al. (2003) expression
    over-predicts D, by a factor of about 10 for albite at 850 C."""
    v = get("plag_Mg_vanorman2014")
    k = get("plag_Mg_costa2003")
    cond = Conditions(T_K=850 + 273.15, X={"XAn": 0.0})
    ratio = float(k.D(cond) / v.D(cond))
    assert 3.0 < ratio < 30.0


def test_plag_Sr_giletti_casserly_matches_the_course_script():
    """DMG Short Course 2025 script: 8.3176e-5*exp(-276000/(T*8.314))*10^(-4.1*XAn)."""
    c = get("plag_Sr_giletti_casserly1994")
    T, xan = 1173.0, 0.5
    expected = 8.3176e-5 * np.exp(-276000.0 / (T * 8.314)) * 10 ** (-4.1 * xan)
    got = float(c.D(Conditions(T_K=T, X={"XAn": xan})))
    assert got == pytest.approx(expected, rel=1e-3)


def test_magnetite_Ti_reproduces_tomiya2013_shinmoedake():
    """Tomiya et al. (2013): Ti at 950 C, log fO2 = -11, X_Usp = 0.3 gives 4.3e-16 m2/s;
    6.9e-16 m2/s at 900 C."""
    c = get("mt_Ti_vanorman_crispin2010")
    d950 = float(c.D(Conditions(T_K=950 + 273.15, log_fo2_bar=-11.0, X={"xTi": 0.1})))
    d900 = float(c.D(Conditions(T_K=900 + 273.15, log_fo2_bar=-11.0, X={"xTi": 0.1})))
    assert d950 == pytest.approx(4.3e-16, rel=0.1)
    assert d900 == pytest.approx(6.9e-16, rel=0.1)


def test_magnetite_has_a_diffusion_minimum_with_temperature():
    """The competing vacancy and interstitial mechanisms make D non-monotonic in T
    (Van Orman & Crispin 2010; noted by Tomiya et al. 2013)."""
    c = get("mt_Ti_vanorman_crispin2010")
    T = np.linspace(1100, 1400, 40)
    D = np.array([float(c.D(Conditions(T_K=t, log_fo2_bar=-11.0, X={"xTi": 0.1}))) for t in T])
    assert np.argmin(D) not in (0, len(D) - 1), "expected an interior minimum in D(T)"


def test_magnetite_fo2_exponents_are_two_thirds():
    """The vacancy branch scales as fO2^(2/3); at oxidising conditions it dominates."""
    c = get("mt_Ti_vanorman_crispin2010")
    T = 1500.0
    d1 = float(c.D(Conditions(T_K=T, log_fo2_bar=-2.0, X={"xTi": 0.1})))
    d2 = float(c.D(Conditions(T_K=T, log_fo2_bar=-1.0, X={"xTi": 0.1})))
    assert np.log10(d2 / d1) == pytest.approx(2.0 / 3.0, rel=0.05)


def test_magnetite_feti_table11_values():
    for key, c0, a, b in (("mt_FeTi_freer_hauptman1978", -15.17, 13.3, 25870.0),
                          ("mt_FeTi_aragon1984", -22.71, 15.09, 19630.0)):
        c = get(key)
        T, xti = 1173.15, 0.15
        expected = np.exp(c0 + a * xti - b / T)
        assert float(c.D(Conditions(T_K=T, X={"xTi": xti}))) == pytest.approx(expected, rel=1e-9)


def test_olivine_anisotropy_is_a_factor_of_six():
    c = get("ol_FeMg_dohmen_chakraborty2007_tamed")
    cond = Conditions(T_K=1473.15, log_fo2_bar=-7.0, X={"XFe": 0.1})
    assert float(c.D(cond.replace(axis="c")) / c.D(cond.replace(axis="a"))) == pytest.approx(6.0)


def test_olivine_fo2_exponent_is_one_sixth():
    c = get("ol_FeMg_dohmen_chakraborty2007_tamed")
    base = Conditions(T_K=1473.15, log_fo2_bar=-7.0, X={"XFe": 0.1}, axis="c")
    hi = base.replace(log_fo2_bar=-4.0)
    assert np.log10(float(c.D(hi) / c.D(base))) == pytest.approx(3.0 / 6.0, rel=1e-9)


def test_unverified_entries_are_flagged_in_their_warnings():
    from dataclasses import replace
    cond = Conditions(T_K=1473.15, log_fo2_bar=-7.0, X={"XFe": 0.1})
    checked = get("ol_FeMg_dohmen_chakraborty2007_tamed")
    assert not any("NOT verified" in s for s in checked.check_conditions(cond))
    unchecked = replace(checked, verified=False, verified_from="a secondary summary")
    assert any("NOT verified" in s for s in unchecked.check_conditions(cond))


def test_every_entry_is_checked_against_its_primary_source():
    # all five former secondary-source entries were read against their PDFs on 1 October 2026
    assert [c.key for c in list_coefficients() if not c.verified] == []


def test_chakraborty1997_is_the_printed_fo86_fit():
    """Chakraborty (1997) p. 12,325: D0 = 5.38e-9 m2/s, Q = 226 kJ/mol, // [001], fO2 = 1e-12 bar."""
    from diffusor.constants import R_GAS
    c = get("ol_FeMg_chakraborty1997")
    T = 1473.15
    D = float(c.D(Conditions(T_K=T, axis="c")))
    assert D == pytest.approx(5.38e-9 * np.exp(-226e3 / (R_GAS * T)), rel=1e-9)
    with pytest.raises(ValueError):           # the paper measured [001] only
        c.D(Conditions(T_K=T, axis="a"))
    off = c.check_conditions(Conditions(T_K=T, log_fo2_bar=-8.0, axis="c"))
    assert any("log fO2" in w for w in off)
    on = c.check_conditions(Conditions(T_K=T, log_fo2_bar=-12.0, axis="c"))
    assert not any("log fO2" in w for w in on)


def test_brady_mccallister1983_eq5():
    """Brady & McCallister (1983) eq. 5: D = 3.89e-7 exp(-360.87 kJ/mol / RT) m2/s."""
    from diffusor.constants import R_GAS
    c = get("cpx_CaMg_brady1983")
    T = 1200 + 273.15
    assert float(c.D(Conditions(T_K=T))) == pytest.approx(3.89e-7 * np.exp(-360.87e3 / (R_GAS * T)), rel=1e-9)
    assert 86.25 * 4.184 == pytest.approx(360.87, abs=0.01)


def test_grove1984_reproduces_its_figure3_line():
    """Grove et al. (1984): D = 10.99 cm2/s exp(-123.4 kcal/mol / RT), ln D = -34.7 at 1400 C."""
    c = get("plag_NaSiCaAl_grove1984")
    lnD_cm2 = np.log(float(c.D(Conditions(T_K=1400 + 273.15))) * 1e4)
    assert lnD_cm2 == pytest.approx(-34.7, abs=0.1)
    lnD_cm2 = np.log(float(c.D(Conditions(T_K=1100 + 273.15))) * 1e4)
    assert lnD_cm2 == pytest.approx(-42.8, abs=0.1)


def test_out_of_range_conditions_produce_a_warning():
    c = get("cpx_FeMg_muller2013")
    w = c.check_conditions(Conditions(T_K=500.0))
    assert any("calibration range" in s for s in w)


def test_direction_cosines_must_close():
    from diffusor.minerals import direction_factor
    with pytest.raises(ValueError):
        direction_factor(1.0, 1.0, 1.0, 10.0, 10.0, 10.0)
    assert direction_factor(2.0, 4.0, 8.0, 90.0, 90.0, 0.0) == pytest.approx(8.0)


# Corrected on 1 October 2026 after a Crossref check; each value was wrong before.
@pytest.mark.parametrize("key, doi", [
    ("sato2022", "10.1016/j.jvolgeores.2022.107686"),
    ("polo_sanchez2023", "10.3389/feart.2023.1128083"),
    ("aggarwal_dieckmann2002", "10.1007/s00269-002-0284-0"),
    ("dohmen2017", "10.2138/rmg.2017.83.16"),
    ("grove1984", "10.1016/0016-7037(84)90391-0"),
    ("sneeringer1984", "10.1016/0016-7037(84)90329-6"),
    ("dimanov_sautter2000", "10.1127/0935-1221/2000/0012-0749"),
])
def test_corrected_reference_dois(key, doi):
    assert REFERENCES[key].doi == doi


def test_corrected_reference_records():
    sato = REFERENCES["sato2022"]
    assert sato.authors.startswith("Sato, M. and Ban, M. and Yuguchi, T. and Adachi, T.")
    assert (sato.volume, sato.pages) == ("432", "107686")
    polo = REFERENCES["polo_sanchez2023"]
    assert "Flaherty" in polo.authors and "Cluzel" not in polo.authors
    assert polo.pages == "1128083"
