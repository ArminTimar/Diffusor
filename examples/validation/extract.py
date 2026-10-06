"""Extract downloaded primary data without digitising or inventing measurements.

Run from the repository root: python examples/validation/extract.py [source_dir]
Requires openpyxl; Ruth's legacy workbook additionally requires xlrd.
Original workbooks are read only. Row numbers and SHA256 hashes are retained.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import re
import sys

import openpyxl

OUT = Path(__file__).resolve().parent
SOURCE = Path(sys.argv[1]) if len(sys.argv) > 1 else Path('papers/Supplementaries and data')
RECORDS = []


def number(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def rows(path, sheet):
    w = openpyxl.load_workbook(path, read_only=True, data_only=True)
    try:
        return list(w[sheet].values)
    finally:
        w.close()


def sheets(path):
    w = openpyxl.load_workbook(path, read_only=True)
    result = w.sheetnames
    w.close()
    return result


def emit(study, sample, path, sheet, data, mapping, scale=1, setup=None, notes='', fo_scale=1):
    """data: (Excel row, row values); mapping: output name -> zero-based column."""
    selected = []
    for rownum, r in data:
        if not all(number(r[mapping[k]]) for k in ('Distance_um', 'Fo_mol')):
            continue
        d = {k: r[i] if number(r[i]) else '' for k, i in mapping.items()}
        d['Distance_um'] *= scale
        for k in ('Fo_mol','Initial_Fo_mol','Fo_sigma'):
            if k in d and number(d[k]):
                d[k] *= fo_scale
        d['source_row'] = rownum
        selected.append(d)
    if len(selected) < 3:
        return
    key = study + '_' + re.sub(r'[^a-z0-9]+', '_', str(sample).lower()).strip('_')
    dest = OUT / 'profiles' / (key + '.csv')
    dest.parent.mkdir(exist_ok=True)
    with dest.open('w', newline='', encoding='utf8') as f:
        writer = csv.DictWriter(f, fieldnames=list(mapping) + ['source_row'])
        writer.writeheader()
        writer.writerows(selected)
    RECORDS.append(dict(key=key, study=study, sample=str(sample), file='profiles/'+dest.name,
                        source_file=path.relative_to(SOURCE).as_posix(), source_sheet=sheet,
                        source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                        n_points=len(selected), setup=setup, notes=notes))


def grouped(study, path, sheet, id_column, mapping, start=1, scale=1, transform=str, notes='', fo_scale=1):
    groups = {}
    for i, r in enumerate(rows(path, sheet), 1):
        identity = ' '.join(str(r[c]) for c in id_column) if isinstance(id_column,tuple) else r[id_column]
        if i < start or not isinstance(identity, str):
            continue
        groups.setdefault(transform(identity), []).append((i, r))
    for name, data in groups.items():
        emit(study, name, path, sheet, data, mapping, scale=scale, notes=notes,fo_scale=fo_scale)


def extract():
    # All 19 archived Lynn traverses, including the two-event Ol 8.
    p = next((SOURCE/'lynn_2024_data').glob('*.xlsx'))
    for sheet in sheets(p):
        if not sheet.startswith('Ol '):
            continue
        rr = rows(p, sheet)
        setup = dict(T_C=rr[5][15], P_MPa=rr[5][16], log_fo2_bar=rr[5][17],
                     angles_deg=list(rr[9][15:18]), published_days=rr[13][15],
                     published_sigma_days=rr[13][16], coefficient='ol_FeMg_dohmen_chakraborty2007_tamed')
        if not all(number(setup[k]) for k in ('T_C','P_MPa','log_fo2_bar','published_days')):
            setup = None
        emit('lynn2024', sheet, p, sheet, enumerate(rr[3:],4),
             dict(Distance_um=1,Fo_mol=10,Initial_Fo_mol=12,Published_model_Fo_mol=13,
                  FeO_wt=3,MgO_wt=4,CaO_wt=5,NiO_wt=7), setup=setup,
             notes='Workbook inputs retained, including 42 MPa versus 45 MPa in the article. '
                   'Initial column is sampled on the EPMA positions. The sub-grid interface and original solver are not archived.')
    p = SOURCE/'mutch2019/Data_S1.xlsx'
    grouped('mutch2019',p,sheets(p)[0],0,
            dict(Distance_um=1,Fo_mol=2,Fo_sigma=3,Ni_ppm=8,Ni_sigma=9,Mn_ppm=10,
                 Mn_sigma=11,Initial_Fo_mol=25,Initial_Ni_ppm=26,Initial_Mn_ppm=27),start=2,fo_scale=100,
            notes='Source Fo, Fo_sigma and initial Fo mole fractions multiplied by 100 to mol%. '
                  'Joint Bayesian Fe-Mg/Ni/Mn inversion with study-specific D regressions; not a TaMED replication.')
    p = next((SOURCE/'mourey2023').glob('*.xlsx'))
    for sheet in sheets(p):
        if not sheet.startswith('KE62-'):
            continue
        rr = rows(p,sheet)
        # Individual sample sheets have different analytical routines/columns.
        header = [str(v).strip().lower() for v in rr[2]]
        find = lambda s: next(i for i,v in enumerate(header) if s in v)
        grouped('mourey2023',p,sheet,0,dict(Distance_um=find('reldist'),Fo_mol=find('fo'),
                FeO_wt=find('feo'),MgO_wt=find('mgo'),CaO_wt=find('cao'),NiO_wt=find('nio')),start=4,
                notes='Profile includes possible multiple events. Rim positions are constrained from Ca and images in the study; no numerical initial arrays in this workbook.')
    p = SOURCE/'sundermeyer2019/410_2019_1642_MOESM1_ESM.xlsx'
    for sheet in sheets(p):
        grouped('sundermeyer2020_reunion',p,sheet,0,
                dict(Distance_um=1,Fo_mol=17,FeO_wt=8,MgO_wt=10,CaO_wt=6,NiO_wt=14),start=7,
                notes='Whole measured traverse retained. The paper excludes the late outermost rim (<3 um) and fits only selected gradients; exact per-profile fit masks are not tabulated.')
    p = SOURCE/'sundermeyer2020/410_2020_1715_MOESM2_ESM.xlsx'
    for sheet in sheets(p):
        grouped('sundermeyer2020_eifel',p,sheet,1,
                dict(Distance_um=0,Fo_mol=25,FeO_wt=12,MgO_wt=10,CaO_wt=14,NiO_wt=18),start=7,
                transform=lambda x:re.sub(r'^Line\d+','',x),
                notes='Crystal names retain traverse suffixes. Growth zones and cooling/effective temperatures need to be matched before a published-time replication.')
    p = SOURCE/'gordeychik2018/41598_2018_30133_MOESM4_ESM.xlsx'
    grouped('gordeychik2018',p,'Table SM2-A',(3,4,5),
            dict(Distance_um=12,Fo_mol=25,FeO_wt=15,MgO_wt=14,CaO_wt=20,NiO_wt=18),start=3,scale=1000,
            notes='Source distances are mm, converted to um. Growth plus diffusion interpretation; retain whole traverse without assuming a sharp-step initial state.')
    p = SOURCE/'ruprecht_plank2023/41586_2013_BFnature12342_MOESM46_ESM.xlsx'
    for sheet in sheets(p):
        if not sheet.startswith('IZ-'):
            continue
        emit('ruprecht2013',sheet,p,sheet,enumerate(rows(p,sheet)[5:],6),
             dict(Distance_um=21,Fo_mol=39,Ni_ppm=36),
             notes='San Carlos-normalised LA-ICPMS values (columns V, AK, AN). '
                   'Study is 2013, despite the local folder name. Full composition/fO2-dependent Ni model is not the registry Fo90 fit.')
    import xlrd
    p = next((SOURCE/'ruth2018').glob('*.xls'))
    w = xlrd.open_workbook(str(p)); s = w.sheet_by_name('3 Data')
    groups = {}
    for i in range(6,s.nrows):
        r=s.row_values(i)
        if r[0]:groups.setdefault(r[0],[]).append((i+1,r))
    for name,data in groups.items():
        emit('ruth2018',name,p,s.name,data,dict(Distance_um=7,Fo_mol=18,FeO_wt=11,MgO_wt=13,CaO_wt=14,NiO_wt=16),
             notes='Measured rim-to-interior profiles; archived composition table does not include per-traverse initial conditions or EBSD angles.')
    (OUT/'manifest.json').write_text(json.dumps(RECORDS,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
    from collections import Counter
    print(dict(Counter(r['study'] for r in RECORDS)))


if __name__ == '__main__':
    extract()
