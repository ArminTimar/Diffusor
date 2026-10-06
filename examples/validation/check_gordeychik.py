"""Recompute published Fo ages and D-only sensitivity, not a raw-profile refit.

Uses the study's fitted widths/Dt and geometry unchanged. The alternative
coefficient is evaluated along the actual a/b/c orientation (its anisotropy
differs from the study's sixfold ratio). Also checks the published Ni equation.
"""
from pathlib import Path
import json
import numpy as np
from diffusor.coefficients import Conditions,get

ROOT=Path(__file__).resolve().parent

def check(r):
    out={k:r[k] for k in ('group','sample','part','row','source_file','published_Fo_days','published_Ni_days')}
    dt=r['Dt_mm2']*1e-6 if r['group']=='advanced_core' else (r['width_Fo_mm']*1e-3)**2/4
    ni_dt=.74*dt if r['group']=='advanced_core' else (r['width_Ni_mm']*1e-3)**2/4
    out['recomputed_Fo_days']=[dt/(d*r['Apr'])/86400 for d in r['D_Fo']]
    out['equation_Ni_days']=[ni_dt/(d*r['Apr'])/86400 for d in r['D_Ni']]
    out['Fo_max_relative_error']=max(abs(a/b-1) for a,b in zip(out['recomputed_Fo_days'],r['published_Fo_days']))
    out['Ni_max_relative_error']=max(abs(a/b-1) for a,b in zip(out['equation_Ni_days'],r['published_Ni_days']))
    out['alternative_Fo_days']=[]
    out['registry_original_D_ratio']=[]
    out['alternative_warnings']=[]
    if r['cosine_squared'] is None:
        out['alternative_to_published_ratio']=[]
        out['classification']=('Published fitted-width/Dt calculation; no matching a/b/c orientation row for alternative D. '
                               +(r.get('orientation_note') or ''))
        return out
    angles=tuple(np.degrees(np.arccos(np.sqrt(np.clip(r['cosine_squared'],0,1)))))
    for i,bound in enumerate(('high','low')):
        s=r['conditions'][bound]
        cond=Conditions(T_K=s['T_C']+273.15,P_Pa=s['P_Pa'],log_fo2_bar=np.log10(s['fo2_Pa'])-5,
                        X={'XFe':1-r['conditions']['Fo']/100},angles_deg=angles)
        old=get('ol_FeMg_dohmen_chakraborty2007_tamed')
        new=get('ol_FeMg_oeser2026')
        out['registry_original_D_ratio'].append(float(old.D(cond))/(r['D_Fo'][i]*r['Apr']))
        out['alternative_Fo_days'].append(dt/float(new.D(cond))/86400)
        out['alternative_warnings'].extend(new.check_conditions(cond))
    out['alternative_to_published_ratio']=[a/b for a,b in zip(out['alternative_Fo_days'],out['published_Fo_days'])]
    out['classification']='Published fitted-width/Dt calculation; not an independent profile fit. Oeser extrapolation, no precision claim.'
    return out

if __name__=='__main__':
    result=[check(r) for r in json.loads((ROOT/'gordeychik_parameters.json').read_text(encoding='utf8'))]
    (ROOT/'gordeychik_results.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
    print(len(result),'Fo age rows; largest relative discrepancy:',max(r['Fo_max_relative_error'] for r in result))
    print('Ni equation discrepancies:',[(r['group'],r['row'],r['Ni_max_relative_error']) for r in result if r['Ni_max_relative_error']>1e-6])
