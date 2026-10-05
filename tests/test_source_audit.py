"""Rules from the October 2026 source audit, in which every law was compared with its source pages."""
import re

import numpy as np
import pytest

from diffusor import datasets as ds
from diffusor.coefficients import Conditions, get
from diffusor.coefficients.families import FAMILIES
from diffusor.coefficients.latex import (coefficient_latex, equation_lines, family_latex, num,
                                         render_png)
from diffusor.coefficients.registry import REGISTRY
from diffusor.constants import R_GAS


def test_correlated_parameters_are_never_sampled_independently_by_default():
    assert all(c.default_sampling_mode() != "independent" for c in REGISTRY.values())
    # the three Schwandt et al. (1998) laws, which spread log D by orders of magnitude when sampled so
    assert get("opx_Mg_schwandt1998_c").default_sampling_mode() == "none"


def test_every_sampled_scatter_states_its_basis():
    for c in REGISTRY.values():
        if c.sigma_logD is not None:
            assert re.match(r"(published|derived|assumed)", c.sigma_logD_basis), c.key
    assert get("opx_FeMg_dias2025").sigma_logD == 0.34
    assert get("opx_FeMg_ganguly_tazzoli1994").sigma_logD == 1.0
    assert get("ol_FeMg_chakraborty1997").sigma_logD_basis.startswith("assumed")


def test_every_law_and_family_has_a_rendering_latex_equation():
    for c in REGISTRY.values():
        for line in equation_lines(coefficient_latex(c)):
            assert render_png(line)[:4] == b"\x89PNG", c.key
    for f in FAMILIES.values():
        assert render_png(family_latex(f).split(r";\quad ")[0])[:4] == b"\x89PNG"


def test_equation_images_have_one_picture_pixel_per_screen_pixel():
    """Qt stretches a picture to the width and height it is given, which made the typeset
    laws soft on a scaled display. The picture is drawn at the screen's pixel ratio and sized
    in logical pixels, so the two agree to within the rounding of the size."""
    import base64
    from diffusor.coefficients.latex import equation_html, png_size
    for key in ("opx_FeMg_dias2025", "ol_FeMg_dohmen_chakraborty2007_tamed"):
        for ratio in (1.0, 1.25, 1.5, 2.0):
            html = equation_html(coefficient_latex(get(key)), max_width=640, ratio=ratio)
            tags = re.findall(r"base64,([^']+)' width='(\d+)' height='(\d+)'", html)
            assert tags
            for data, w, h in tags:
                pw, ph = png_size(base64.b64decode(data))
                assert abs(pw / ratio - int(w)) <= 0.5 and abs(ph / ratio - int(h)) <= 0.5, (key, ratio)
                assert int(w) <= 640, (key, ratio)


def test_long_equations_are_cut_at_operators_and_every_piece_renders():
    from diffusor.coefficients.latex import BASE_DPI, wrap_equation
    law = equation_lines(coefficient_latex(get("ol_FeMg_dohmen_chakraborty2007_tamed")))[0]
    pieces = wrap_equation(law, 600, BASE_DPI)
    assert len(pieces) >= 2
    for piece in pieces:
        assert render_png(piece)[:4] == b"\x89PNG"
    # nothing is cut inside parentheses, exponents or \left ... \right
    assert all(p.count("(") == p.count(")") for p in pieces)
    assert all(p.count(r"\left") == p.count(r"\right") for p in pieces)
    # a short line is left alone
    assert wrap_equation(r"D = 1", 600, BASE_DPI) == [r"D = 1"]


def test_latex_carries_the_parameter_values_the_code_uses():
    assert num(262914) == "262914" and num(2.77e-7) == r"2.77\times10^{-7}" and num(555.425) == "555.425"
    s = coefficient_latex(get("cpx_FeMg_muller2013"))
    assert r"2.77\times10^{-7}" in s and "320.7" in s
    s = coefficient_latex(get("plag_Mg_audetat2026"))
    assert "-2.99" in s and "262914" in s and "1.87" in s
    # the olivine Ca fO2 term belongs to every axis
    assert all("10^{-12}" in line for line in equation_lines(coefficient_latex(get("ol_Ca_coogan2005"))))


def test_dias_dohmen_2024_reference_composition_is_0_1():
    """Eqs 7 and 12: D(X_Fe) = D(X_Fe = 0.1) 10^(m (X_Fe - 0.1))."""
    c = get("opx_FeMg_dias_dohmen2024")
    T = 1273.15
    D = c.D(Conditions(T, X={"XFe": 0.1}, axis="c"))
    assert D == pytest.approx(3.8e-9 * np.exp(-261.07e3 / (R_GAS * T)))


def test_santorini_example_is_removed():
    assert "plag_santorini" not in {d.key for d in ds.DATASETS}
    assert not (ds.EXAMPLES_DIR / "plagioclase_santorini_druitt2012.csv").exists()


def test_secondary_only_laws_are_not_marked_verified():
    assert not get("mt_FeTi_aragon1984").verified
    assert "Costa et al. (2003) eq. 8" in get("plag_Mg_costa2003").verified_from
