"""Source-audit decisions for the accessory-mineral laws (ranges, labels and stated assumptions)."""
import numpy as np
import pytest

from diffusor.coefficients.registry import get
from diffusor.coefficients import isotopes as iso


def _celsius(c):
    return c.T_range.lo - 273.15, c.T_range.hi - 273.15


def test_zircon_ree_ranges_follow_the_data_of_each_element():
    # Tables 2 and 3 of Cherniak et al. (1997); the abstract's 1150-1400 C is the union of the three
    assert _celsius(get("zrn_Sm_cherniak1997")) == pytest.approx((1200, 1400))
    assert _celsius(get("zrn_Dy_cherniak1997")) == pytest.approx((1150, 1350))
    assert _celsius(get("zrn_Yb_cherniak1997")) == pytest.approx((1150, 1350))
    for sp in ("Sm", "Dy", "Yb"):
        c = get(f"zrn_{sp}_cherniak1997")
        assert "dryness not stated" in c.reference_state and "dry lattice" not in c.reference_state


def test_rutile_hf_ranges_use_the_part_common_to_abstract_and_table():
    # the abstract and Table 1 of Cherniak et al. (2007) assign 800-1000 and 750-1050 C to opposite directions
    for key in ("rt_Hf_cherniak2007_c", "rt_Hf_cherniak2007_a"):
        c = get(key)
        assert _celsius(c) == pytest.approx((800, 1000))
        assert any("part common to both readings" in n for n in c.calibration_notes)
    assert _celsius(get("rt_Zr_cherniak2007_c")) == pytest.approx((750, 1100))


def test_implanted_sm_is_a_tracer_type_exchange_law_and_says_why():
    c = get("ap_Sm_cherniak2000_implant")
    assert c.kind == "tracer" and "Diffusor's classification" in c.notes
    assert "Table 5" in c.uncertainty_note and "-6.20" in c.uncertainty_note
    # the label changes the warnings, not D
    from diffusor.coefficients.base import Conditions
    D = c.D(Conditions(1273.15))
    assert D == pytest.approx(10 ** np.log10(6.3e-7) * np.exp(-298e3 / (8.314462618 * 1273.15)))


def test_assumed_pressure_and_kind_are_labelled_as_assumptions():
    assumed = get("zrn_Sm_cherniak1997")
    assert "not stated by the authors" in assumed.P_range.unit and "pressure is not stated" in assumed.notes
    stated = get("qz_Ti_cherniak2007")
    assert "1 atm stated by the authors" in stated.P_range.unit
    assert "Diffusor's classification" in get("ap_Pb_cherniak1991").notes
    # the authors call the titanite and rutile coefficients chemical diffusion: no disclaimer
    assert "classification" not in get("ttn_Sr_cherniak1995").notes
    for c in (assumed, stated):
        assert (c.P_range.lo, c.P_range.hi) == (1.0e5, 101325.0)


def test_petry_note_quotes_the_printed_errors():
    note = get("ol_Ni_petry2004").uncertainty_note
    assert "216 +/- 6" in note and "+2.5/-1.5" in note and "give no error" not in note


def test_printed_alternatives_are_documented_where_values_were_kept():
    assert "2.95e-9" in " ".join(get("xtm_Pb_cherniak2006").calibration_notes)
    assert "-0.962" in " ".join(get("zrn_Pb_cherniak2001").calibration_notes)
    assert "563" in " ".join(get("mnz_Pb_cherniak2004").calibration_notes)
    assert "7.01e-8" in get("qz_Ti_cherniak2007").verified_from
    assert "parallel to (100)" in " ".join(get("ttn_Sr_cherniak1995").calibration_notes)
    assert "perpendicular to c" in " ".join(get("ap_Sr_cherniak1993").calibration_notes)
    # the equation values themselves are the abstract values
    assert get("xtm_Pb_cherniak2006").params["logD0"].value == pytest.approx(np.log10(3.0e-9))
    assert get("zrn_Pb_cherniak2001").params["logD0"].value == pytest.approx(np.log10(0.11))


def test_teng_beta_is_described_as_one_binary_with_borrowed_beta():
    for key in ("ol_Fe_teng2011", "ol_Mg_teng2011"):
        b = iso.BETAS[key]
        assert "independent" not in b.diffusion_model and "one coefficient" in b.diffusion_model
        assert b.assumed and "silicate melts" in b.notes
