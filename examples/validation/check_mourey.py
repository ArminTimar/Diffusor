"""Conditional Fe-Mg-only reconstruction of Mourey's single-event 3311_2_ol6.

The original interface was selected jointly with Ca and Ni. Its position is
not tabulated, so it is fitted here using Fo alone. This is not an exact replay.
"""
from pathlib import Path
from dataclasses import replace
import json
import numpy as np
import pandas as pd
from diffusor.coefficients import Conditions, get
from diffusor.constants import SEC_PER_DAY
from diffusor.fitting import DiffusionModel, fit_time
from diffusor.solvers import InitialCondition, dirichlet

ROOT = Path(__file__).resolve().parent


def compare():
    df = pd.read_csv(ROOT/'profiles/mourey2023_3311_ol6.csv')
    x, y = df.Distance_um.to_numpy(), df.Fo_mol.to_numpy()
    # MOESM1, Fe-Mg timescales, 3311_2_ol6: first-event temperature,
    # core/rim, EBSD angles and Fo80-82 scenario. Paper equations 2-4, p4.
    T = 1183 + 273.15
    conditions = Conditions(T_K=T, P_Pa=60e6,
        log_fo2_bar=8.912-25160/T,
        angles_deg=(24.7,65.3,89), X={'XFe':1-np.mean(y)/100})
    initial = InitialCondition('step', dict(x0=20., C_left=80.91, C_right=89.15))
    result = dict(study='mourey2023', sample='3311_2_ol6', published_days=395,
        published_plus_days=149, published_minus_days=106,
        classification='conditional Fo-only reconstruction, not exact replication',
        assumptions=[
            'Source T=1183 C, P=60 MPa and EBSD=(24.7,65.3,89). Equation 3 gives log fO2 = 8.912-25160/T "in Pa" for the QFM buffer named in the text. Read in bar it is -8.37 at 1183 C, 0.09 log units from Diffusor\'s QFM; read in Pa it is about 5 log units below QFM. Both readings are run.',
            'Source Fo core=89.15 and rim=80.91. Fixed compositions at both measured endpoints.',
            'All 25 measured Fo points, uniform weights, no beam correction, plane geometry.',
            'The source jointly located the initial rim using Ca, Ni and Fo. Here t and x0 are fitted to Fo only.',
            'Equation 2 prints the composition term as 3*(XMg-0.9); Dohmen & Chakraborty (2007, as corrected by their erratum) give 3*(XFe-0.1). The printed-forms branch keeps 3*(XMg-0.9) and the Pa reading.',
            'The main branch uses the Dohmen & Chakraborty (2007) TaMED law with the erratum and reads equation 3 in bar. These are Diffusor\'s readings, not a statement about the authors\' code.',
            'Alternative-law run changes only the coefficient. Oeser atmospheric-pressure and fO2 calibration is extrapolated.'],
        fits=[])
    canonical=get('ol_FeMg_dohmen_chakraborty2007_tamed')
    literal=replace(canonical, key='mourey2023_equation2_literal',
                    params={**canonical.params,'m':replace(canonical.params['m'],value=-3.)})
    cases=[('printed forms of equations 2-3',literal,
            replace(conditions,log_fo2_bar=conditions.log_fo2_bar-5)),
           ('canonical D, equation 3 interpreted in bar',canonical,conditions),
           ('Oeser D, same bar interpretation',get('ol_FeMg_oeser2026'),conditions)]
    for label,coefficient,cond in cases:
        model = DiffusionModel(coefficient, cond, initial,
            bc_left=dirichlet(80.91),bc_right=dirichlet(89.15),boundaries_far=False,
            comp_key='XFe',comp_scale=-.01,comp_offset=1,
            x_grid=np.linspace(x.min(),x.max(),401),n_nodes=401)
        fit = fit_time(model,x,y,free_parameters=('t','x0'),
                       t_min=.01*SEC_PER_DAY,t_max=1e5*SEC_PER_DAY)
        result['fits'].append(dict(label=label,coefficient=coefficient.key,days=fit.t_seconds/SEC_PER_DAY,
            success=bool(fit.success),warnings=fit.warnings,
            rmse_Fo=float(np.sqrt(np.mean((fit.C_model-y)**2))),
            initial_interface_um=float(fit.model.initial.params['x0'])))
        print(label,fit.t_seconds/SEC_PER_DAY,flush=True)
    result['alternative_to_original_ratio']=result['fits'][2]['days']/result['fits'][1]['days']
    return result


if __name__=='__main__':
    result=compare()
    (ROOT/'mourey_results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print(json.dumps(result,indent=2))
