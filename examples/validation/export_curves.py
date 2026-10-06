"""Regenerate comparison curves at saved fitted times (no optimisation)."""
from pathlib import Path
import json
from diffusor.validation import archived_model,records

ROOT=Path(__file__).resolve().parent

if __name__=='__main__':
    lookup={r['key']:r for r in records()}
    dest=ROOT/'fits'; dest.mkdir(exist_ok=True)
    for result in json.loads((ROOT/'results.json').read_text(encoding='utf8')):
        record=lookup[result['key']]
        original,df=archived_model(record,n_nodes=result['n_nodes'])
        alternative,_=archived_model(record,coefficient='ol_FeMg_oeser2026',n_nodes=result['n_nodes'])
        x=df.Distance_um.to_numpy()
        df['Diffusor_original_fit_Fo']=original.profile(result['original_days']*86400,x)
        df['Diffusor_at_published_time_Fo']=original.profile(result['published_days']*86400,x)
        df['Diffusor_alternative_fit_Fo']=alternative.profile(result['alternative_days']*86400,x)
        df.to_csv(dest/(record['key']+'.csv'),index=False)
    print('Exported',len(list(dest.glob('*.csv'))),'curve comparisons')
