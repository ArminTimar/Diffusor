"""Independent source checks and preservation of scientific setup in the UI."""
import json
from pathlib import Path
import runpy

import numpy as np
import pandas as pd
import pytest

from diffusor.validation import ROOT, archived_model, records


def test_manifest_has_unique_profiles_with_valid_units_and_source_rows():
    rr=records()
    assert len({r['key'] for r in rr})==len(rr)
    assert len({r['file'] for r in rr})==len(rr)
    for r in rr:
        df=pd.read_csv(ROOT/r['file'])
        assert len(df)==r['n_points']>=3
        value=r.get('value_column','Fo_mol')
        assert np.isfinite(df[value]).all(),r['key']
        if value=='Fo_mol':
            assert df.Fo_mol.between(0,100).all(),r['key']
        assert np.isfinite(df.Distance_um).all()
        source=df.source_row if 'source_row' in df else df.source_cell
        assert source.is_unique,r['key']
        assert len(r['source_sha256'])==64
    mutch=next(r for r in rr if r['key']=='mutch2019_borg14_ol_c2_p1')
    df=pd.read_csv(ROOT/mutch['file'])
    # Independent cells C2/D2/Z2 of Data_S1, converted from mole fraction.
    assert df.Fo_mol.iloc[0]==pytest.approx(85.68603814)
    assert df.Fo_sigma.iloc[0]==pytest.approx(.11786643)
    assert df.Initial_Fo_mol.iloc[0]==pytest.approx(85.68603814)


def test_lynn_archived_model_excludes_unmodelled_rim_and_retains_conditions():
    r=next(r for r in records() if r['key']=='lynn2024_ol_1')
    model,df=archived_model(r)
    assert df.Distance_um.min()==10 # row 20 has measurements but no model
    assert model.conditions.P_Pa==42e6
    assert model.conditions.angles_deg==(93,29,117)
    assert model.initial.kind=='table'
    assert model.bc_left.value==82.4
    assert model.bc_right.value==88.5
    assert model._composition(88.5)==pytest.approx(.115)
    with pytest.raises(ValueError,match='two events'):
        archived_model(next(r for r in records() if r['key']=='lynn2024_ol_8'))


def test_lynn_ol19_refit_lands_inside_the_published_interval():
    # the shortest published time (14 +/- 4.2 d); the saved refit is 14.53 d
    from diffusor.constants import SEC_PER_DAY
    from diffusor.fitting import fit_time
    r=next(r for r in records() if r['key']=='lynn2024_ol_19')
    model,df=archived_model(r)
    fit=fit_time(model,df.Distance_um.to_numpy(),df.Fo_mol.to_numpy(),free_parameters=('t',),
                 t_min=SEC_PER_DAY*.001,t_max=SEC_PER_DAY*1e5)
    days=fit.t_seconds/SEC_PER_DAY
    assert days==pytest.approx(14.53,rel=1e-3)
    assert abs(days-r['setup']['published_days'])<=r['setup']['published_sigma_days']


def test_gordeychik_published_age_algebra_and_the_sm4_ni_formula():
    check=runpy.run_path(str(ROOT/'check_gordeychik.py'))['check']
    rr=json.loads((ROOT/'gordeychik_parameters.json').read_text(encoding='utf8'))
    results=[check(r) for r in rr]
    assert len(results)==32
    assert max(r['Fo_max_relative_error'] for r in results)<1e-12
    bad=[r for r in results if r['Ni_max_relative_error']>1e-6]
    assert len(bad)==5 and all(r['group']=='outer_core' for r in bad)
    assert results[0]['recomputed_Fo_days'][0]==pytest.approx(646.7225871921327)
    # an orientation is paired only when its geometric factor is the one the age row used,
    # so Diffusor's projected D_FeMg matches the spreadsheet's D x Apr within the constants
    for r in rr:
        if r['cosine_squared'] is not None:
            c=r['cosine_squared']
            assert c[0]/6+c[1]/6+c[2]==pytest.approx(r['Apr'],rel=1e-9),r['sample']
    paired=[v for r in results for v in r['registry_original_D_ratio']]
    assert len(paired)==50 and min(paired)>.996 and max(paired)<.998
    ol82=next(r for r in rr if r['sample']=='SHIV-08-05 17 Ol-8-2')
    assert ol82['cosine_squared'] is None and 'Ol-8-3' in ol82['orientation_note']


def test_lynn_gui_uses_identical_archived_model_and_resets_cooling():
    import os
    os.environ.setdefault('QT_QPA_PLATFORM','offscreen')
    from PySide6.QtWidgets import QApplication
    from diffusor.gui.main_window import MainWindow
    from diffusor.datasets import get
    from diffusor.dataio import ProfileSpec,read_table
    app=QApplication.instance() or QApplication([])
    window=MainWindow()
    try:
        window.chk_cooling.setChecked(True)
        d=get('lynn2024_ol_1')
        assert window._use_table(read_table(d.path),ProfileSpec(**d.spec),d.path,d,(None,False))
        model=window._model(d.settings['coefficient'])
        expected,df=archived_model(next(r for r in records() if r['key']==d.key))
        assert not window.chk_cooling.isChecked()
        assert window._free_parameters()==['t']
        assert model.initial.kind=='table'
        assert model.conditions.angles_deg==expected.conditions.angles_deg
        assert not model.boundaries_far
        np.testing.assert_allclose(model.initial.evaluate(model.x_grid),expected.initial.evaluate(expected.x_grid))
        assert model.bc_left.value==expected.bc_left.value
        assert model.comp_scale==expected.comp_scale
        # Changing to an ordinary example must remove the archived table option.
        d=get('olivine_laki')
        assert window._use_table(read_table(d.path),ProfileSpec(**d.spec),d.path,d,(None,False))
        assert not window._table_ic()
        assert window.cmb_ic.findData('table')==-1
        # The Ni study opens its measured Ni rather than its auxiliary Fo column.
        record=next(r for r in records() if r['study']=='ruprecht2013')
        window._load_validation_record(record)
        assert window.cmb_species.currentText()=='Ni'
        assert window.profile.spec.column_a=='Ni_ppm'
        assert window.dataset is None
        assert 'no verified model preset' in window.lbl_data.text()
        assert any('retained from your session' in label.text()
                   for label in window._prefill.values())
    finally:
        window.close()
        app.processEvents()


def test_iovine_refits_reproduce_tables_1_and_2():
    """Iovine et al. (2017), 930 C, Cherniak (2002) Ba law, sharp step with free plateaus."""
    result=json.loads((ROOT/'iovine_results.json').read_text(encoding='utf8'))
    assert result['n']==23 and result['coefficient']=='kfs_Ba_cherniak2002'
    assert all(.85<r['ratio']<1.3 for r in result['fits'])
    assert result['n_equal_after_rounding']>=18
    # the paper's Table 1 time is printed as a whole year; the two-step traverses stay data-only
    rr={r['key']:r for r in records() if r['study']=='iovine2017'}
    assert rr['iovine2017_a_cx6_grey_l1']['setup']['published_years']==6
    assert rr['iovine2017_b_cx1_xray']['setup'] is None
    assert rr['iovine2017_a_cx5_grey_l1']['setup'] is None
    assert len(rr)==43


def test_iovine_preset_fits_like_the_check_script():
    import os
    os.environ.setdefault('QT_QPA_PLATFORM','offscreen')
    from PySide6.QtWidgets import QApplication
    from diffusor.gui.main_window import MainWindow
    from diffusor.datasets import get
    from diffusor.dataio import ProfileSpec,read_table
    from diffusor.constants import SEC_PER_YEAR
    from diffusor.fitting import fit_time
    app=QApplication.instance() or QApplication([])
    window=MainWindow()
    try:
        d=get('iovine2017_d_cx6_grey_l1')
        assert d.mineral=='kfeldspar' and d.citation=='iovine2017'
        assert window._use_table(read_table(d.path),ProfileSpec(**d.spec),d.path,d,(None,False))
        assert window._checked_keys()==['kfs_Ba_cherniak2002']
        assert window._boundaries_far()
        assert set(window._free_parameters())=={'t','x0','C_left','C_right'}
        p=window.profile
        r=fit_time(window._model('kfs_Ba_cherniak2002'),p.x,p.C,p.sigma,window._free_parameters())
        saved=next(v for v in json.loads((ROOT/'iovine_results.json').read_text(encoding='utf8'))['fits']
                   if v['key']==d.key)
        assert r.t_seconds/SEC_PER_YEAR==pytest.approx(saved['diffusor_years'],rel=.01)
        assert r.t_seconds/SEC_PER_YEAR==pytest.approx(151,rel=.05)   # Table 1, D cx6_L1
        # a data-only pyroxene traverse opens as opx Fe-Mg from its oxides, without a preset
        record=next(r for r in records() if r['key']=='sato2022_okp_4_4_3_px_16_no1')
        window._load_validation_record(record)
        assert window.cmb_mineral.currentData()=='opx' and window.dataset is None
        assert window.profile.spec.column_a=='FeO_wt' and window.profile.spec.column_b=='MgO_wt'
    finally:
        window.close()
        app.processEvents()


def test_published_timescale_tables_recompute_from_their_own_laws():
    result=json.loads((ROOT/'timescale_tables_results.json').read_text(encoding='utf8'))
    for sheet,rows in (('Sr in sanidine',17),('Ba in sanidine',54),('Ti in quartz',151)):
        assert result['chamberlain2014'][sheet]['rows']==rows
        assert result['chamberlain2014'][sheet]['max_relative_difference']<1e-4
    petrone=result['petrone2018']
    assert petrone['rows_checked']==71 and petrone['within_5_percent_or_rounding']==62
