"""Araya 2024 coefficient check, not a profile or age replication.

The paper gives constant log10 D=-19.78, 966 C and NNO. Its 'perpendicular
to c' orientation does not distinguish a from b. a reproduces the quoted D.
1 bar is a stated computational reference here, not a claimed magma pressure.
"""
from pathlib import Path
import json
import numpy as np
from diffusor.coefficients import Conditions,get
from diffusor.thermo import log_fo2_from_delta

if __name__=='__main__':
    T=966+273.15
    baseline=10**-19.78
    cond=Conditions(T_K=T,P_Pa=1e5,log_fo2_bar=log_fo2_from_delta('NNO',0,T,1e5),X={'XFe':.09},axis='a')
    old=get('opx_FeMg_dohmen2016')
    new=get('opx_FeMg_dias2025')
    result=dict(published_logD=-19.78,recomputed_logD=float(np.log10(old.D(cond))),
                source='Araya et al. 2024, sections 5.1.4 and 5.2, doi:10.1029/2023JB028558',
                assumptions='a axis reproduces the quoted D; 1 bar buffer reference; XFe=.09 disables the old composition correction. These are not fully specified per-profile inputs.',comparisons=[])
    for axis in ('a','b'):
        for xfe in (.09,.31,.39):
            c=cond.replace(axis=axis,X={'XFe':xfe})
            d=float(new.D(c))
            result['comparisons'].append(dict(axis=axis,XFe=xfe,alternative_logD=float(np.log10(d)),
                  time_multiplier_for_fixed_Dt=baseline/d,warnings=new.check_conditions(c)))
    result['scope']='Constant-D sensitivity only. XFe .31-.39 brackets the reported core Mg# 61-69. No new ages inferred without profiles/orientations and composition-dependent refits.'
    (Path(__file__).parent/'araya_results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print(result['recomputed_logD'])
