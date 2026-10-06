"""Numbers that the README quotes in its verification table and registry counts.

The README states what Diffusor gives next to the published value. If a law changes, these
tests fail and the README table has to be recomputed from the code, not copied from an older
run. The published values are those printed in the sources; Diffusor's values are asserted
at the precision the README prints.
"""
import collections

import numpy as np
import pytest

from diffusor.coefficients import Conditions, get, list_coefficients
from diffusor.coefficients.families import get_family
from diffusor.coefficients.transport import ideal_ionic_matrix


def _tomiya(species, T_C):
    return float(get(f"mt_{species}_vanorman_crispin2010").D(
        Conditions(T_C + 273.15, log_fo2_bar=-11.0, X={"xTi": 0.1})))


@pytest.mark.parametrize("species, T_C, printed, shown", [
    ("Ti", 950, 4.3e-16, 4.15e-16), ("Ti", 900, 6.9e-16, 6.83e-16),
    ("Fe", 950, 4.4e-15, 4.17e-15), ("Fe", 900, 6.6e-15, 6.54e-15)])
def test_tomiya_rows_of_the_readme(species, T_C, printed, shown):
    got = _tomiya(species, T_C)
    assert got == pytest.approx(shown, rel=5e-3)           # the value the README shows
    assert abs(got / printed - 1.0) < 0.055                 # README: within 5.2 %


def test_dias_dohmen_2024_table1_runs_of_the_readme():
    """Per-run log D at X_Fe 0.09 from Table 1 (log D0 + m 0.09). Four runs agree with both laws
    within 0.10 log units. The 1000 C run lies 0.37 (2024) and 0.34 (2025) below them."""
    runs = [(1102, 3.7, -18.69), (1050, 3.0, -18.95), (1050, 2.9, -19.11), (950, 1.1, -19.72)]
    for law in ("opx_FeMg_dias_dohmen2024", "opx_FeMg_dias2025"):
        c = get(law)
        for T, m, logd0 in runs:
            cond = Conditions(T + 273.15, log_fo2_bar=-12.0, X={"XFe": 0.09}, axis="c")
            assert abs(float(c.log10_D(cond)) - (logd0 + 0.09 * m)) < 0.10, (law, T)
    cond = Conditions(1000 + 273.15, log_fo2_bar=-12.0, X={"XFe": 0.09}, axis="c")
    run = -19.74 + 0.09 * 2.4
    assert float(get("opx_FeMg_dias_dohmen2024").log10_D(cond)) - run == pytest.approx(0.37, abs=0.01)
    assert float(get("opx_FeMg_dias2025").log10_D(cond)) - run == pytest.approx(0.34, abs=0.01)


def test_petry_ni10_and_audetat_sr_gap_of_the_readme():
    cond = Conditions(1278.15, log_fo2_bar=-10.96, X={"XFe": 0.1}, axis="c")
    assert float(get("ol_Ni_petry2004").log10_D(cond)) == pytest.approx(-17.24, abs=0.005)
    cond = Conditions(1023.15, X={"XAn": 0.36})
    gap = float(get("plag_Sr_giletti_casserly1994").log10_D(cond)
                - get("plag_Sr_grocolas2025").log10_D(cond))
    assert gap == pytest.approx(2.77, abs=0.005)


def test_garnet_eq5_matrix_from_eq3_tracers_is_17_to_36_percent_low():
    """The README says the tracer inputs of the eq. 5 check are not the eq. 3 values, and that
    the eq. 3 values give a matrix 17 to 36 % below the printed one."""
    printed = np.array([[1.37e-12, -1.53e-13, -1.53e-13],
                        [-2.43e-13, 5.67e-13, -3.32e-14],
                        [-8.24e-13, -1.12e-13, 4.88e-13]])
    X = {"Fe": 0.61, "Mg": 0.18, "Mn": 0.20, "Ca": 0.01}
    d = get_family("grt_chakraborty_ganguly1992").tracer(1430.0 + 273.15, 41e8, X)
    D = np.array([d["Mn"], d["Mg"], d["Fe"], d["Ca"]], dtype=float) * 1e4          # cm2/s
    assert D[:3] == pytest.approx([1.52e-12, 4.46e-13, 4.52e-13], rel=5e-3)
    M = ideal_ionic_matrix(D, np.array([0.20, 0.18, 0.61, 0.01]), dependent=-1)
    low = 1.0 - M / printed
    assert low.min() == pytest.approx(0.17, abs=0.01) and low.max() == pytest.approx(0.36, abs=0.01)


def test_registry_counts_of_the_readme():
    cs = list_coefficients()
    assert len(cs) == 113 and sum(bool(c.verified) for c in cs) == 112
    assert dict(collections.Counter(c.mineral for c in cs)) == {
        "olivine": 10, "opx": 12, "cpx": 3, "plagioclase": 11, "kfeldspar": 7, "magnetite": 42,
        "garnet": 3, "quartz": 1, "rutile": 3, "titanite": 2, "apatite": 8, "zircon": 6,
        "monazite": 1, "xenotime": 4}
    with_cov = [c for c in cs if getattr(c, "cov_order", None)]
    scatter = [c for c in cs if c.sigma_logD is not None and not getattr(c, "cov_order", None)]
    neither = [c for c in cs if c.sigma_logD is None and not getattr(c, "cov_order", None)]
    assert (len(with_cov), len(scatter), len(neither)) == (2, 68, 43)
    stated = [c for c in scatter if not c.sigma_logD_basis.startswith("assumed")]
    assert (len(stated), len(scatter) - len(stated)) == (7, 61)


def test_published_study_rows_of_the_readme():
    """The Lynn, Gordeychik and Araya rows, from the saved validation results."""
    import json
    from diffusor.thermo import log_fo2_from_delta
    from diffusor.validation import ROOT
    lynn = json.loads((ROOT / "results.json").read_text(encoding="utf8"))
    ratio = [r["original_days"] / r["published_days"] for r in lynn]
    assert len(lynn) == 18 and all(r["within_published_interval"] for r in lynn)
    assert round(min(ratio), 2) == 0.75 and round(max(ratio), 2) == 1.06
    assert round(max(r["published_time_curve_rmse_Fo"] for r in lynn), 2) == 0.10
    ages = json.loads((ROOT / "gordeychik_results.json").read_text(encoding="utf8"))
    assert len(ages) == 32 and max(r["Fo_max_relative_error"] for r in ages) < 4e-16
    T = 966 + 273.15
    c = Conditions(T, P_Pa=1e5, log_fo2_bar=log_fo2_from_delta("NNO", 0, T, 1e5),
                   X={"XFe": .09}, axis="a")
    assert np.log10(get("opx_FeMg_dohmen2016").D(c)) == pytest.approx(-19.782, abs=5e-4)
    mourey = json.loads((ROOT / "mourey_results.json").read_text(encoding="utf8"))
    assert round(mourey["fits"][1]["days"]) == 406
