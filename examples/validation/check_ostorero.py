"""Refit the existing Kizimen example with original and newer Opx coefficients.

The authors fitted BSE profiles. This is the archived EPMA surrogate, not an
exact replication. Conditions/window/beam match the bundled example preset.
"""
from pathlib import Path
import json
import numpy as np
from diffusor.datasets import get as dataset
from diffusor.coefficients import Conditions,get
from diffusor.dataio import ProfileSpec,build_profile,read_table
from diffusor.fitting import DiffusionModel,fit_time
from diffusor.solvers import dirichlet
from diffusor.solvers.initial import guess_step_from_data
from diffusor.thermo import log_fo2_from_delta
from diffusor.constants import SEC_PER_YEAR

if __name__=='__main__':
    d=dataset('opx_kizimen')
    p=build_profile(read_table(d.path),ProfileSpec(**d.spec))
    s=d.settings
    T=s['T_C']+273.15; P=s['P_MPa']*1e6
    cond=Conditions(T_K=T,P_Pa=P,log_fo2_bar=log_fo2_from_delta(s['buffer'],s['delta_buffer'],T,P),
                    X={'XFe':float(np.mean(p.C))},axis=s['axis'])
    initial=guess_step_from_data(p.x,p.C)
    result=dict(study='ostorero2022',sample='K9_L10C4',published_years=2.32,
                classification='EPMA surrogate for published BSE fit, not exact replication',fits=[])
    for key in (s['coefficient'],'opx_FeMg_dias2025'):
        model=DiffusionModel(get(key),cond,initial,comp_key='XFe',
                            bc_left=dirichlet(initial.params['C_left']),bc_right=dirichlet(initial.params['C_right']),
                            x_grid=np.linspace(p.x.min(),p.x.max(),301),beam_sigma_um=.5)
        fit=fit_time(model,p.x,p.C,p.sigma,free_parameters=('t','x0'))
        result['fits'].append(dict(coefficient=key,years=fit.t_seconds/SEC_PER_YEAR,
                                  success=bool(fit.success),warnings=fit.warnings))
        print(key,fit.t_seconds/SEC_PER_YEAR,flush=True)
    (Path(__file__).parent/'ostorero_results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
