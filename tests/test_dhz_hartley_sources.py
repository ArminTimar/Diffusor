"""Deer, Howie & Zussman (2013) citations and the Hartley et al. (2016) olivine example."""
import re

import pytest

import diffusor.datasets as ds
from diffusor.minerals.base import CompositionVariable
from diffusor.minerals.definitions import MINERALS, X_AN, X_FE, X_FO, X_OR
from diffusor.references import REFERENCES, get as get_reference


def test_second_edition_of_deer_howie_zussman_is_not_cited():
    assert "deer1992" not in REFERENCES
    assert get_reference("deer2013").year == 2013
    assert CompositionVariable("k", "l", "d").citation == "deer2013"
    for var in (X_FE, X_FO, X_AN, X_OR):
        assert var.citation == "deer2013"


@pytest.mark.parametrize("key, page", [
    ("olivine", "p. 5"), ("opx", "p. 102"), ("cpx", "p. 112"), ("plagioclase", "p. 292"),
    ("kfeldspar", "p. 253"), ("magnetite", "p. 402"), ("garnet", "p. 18"), ("quartz", "p. 311"),
    ("rutile", "p. 393"), ("titanite", "p. 15"), ("apatite", "p. 473"), ("zircon", "p. 12"),
    ("monazite", "p. 478")])
def test_every_dhz_statement_carries_its_printed_page(key, page):
    assert "DHZ 2013" in MINERALS[key].notes and page in MINERALS[key].notes


def test_xenotime_is_not_attributed_to_the_book():
    note = MINERALS["xenotime"].notes
    assert "no data sheet in Deer, Howie & Zussman (2013)" in note
    assert "Cherniak (2006)" in note


def test_quartz_note_gives_the_book_limits_without_inventing_a_pressure():
    note = MINERALS["quartz"].notes
    assert "573 C" in note and "870 C" in note and "no pressure" in note
    assert "1 bar" not in note


def test_anorthite_and_orthoclase_definitions_name_the_ternary_form_and_page():
    assert "Ca/(Ca+Na+K)" in X_AN.definition and "p. 489" in X_AN.definition
    assert "Ca/(Ca+Na) without K" in X_AN.definition
    assert "K/(K+Na+Ca)" in X_OR.definition and "p. 489" in X_OR.definition
    assert "p. 292" in MINERALS["plagioclase"].notes and "ternary" in MINERALS["plagioclase"].notes


def test_mineral_and_composition_text_has_no_semicolons():
    for m in MINERALS.values():
        assert ";" not in m.notes and ";" not in str(m.composition_variable), m.key


def test_hartley_example_states_what_was_checked():
    ex = ds.get("olivine_laki")
    s = ex.settings
    assert (s["T_C"], s["sigma_T_K"], s["buffer"], s["delta_buffer"], s["sigma_delta"]) == (
        1150.0, 30.0, "FMQ", -1.0, 0.5)
    text = ex.provenance + " ".join(ex.sources.values())
    assert "not checked" not in text.lower()
    assert "p. 61" in text
    assert re.search(r"1.5 kbar", text)
    # the law is the one the study cites (2007), not a later one
    assert s["coefficient"] == "ol_FeMg_dohmen_chakraborty2007_tamed"
    ref = get_reference("hartley2016")
    assert (ref.year, ref.volume, ref.pages, ref.doi) == (
        2016, "439", "58-70", "10.1016/j.epsl.2016.01.018")
