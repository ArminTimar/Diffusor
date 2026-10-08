"""Reproducible comparisons with archived measured profiles.

The archived initial curve is interpolated, never inferred from the observations.
These runs test agreement with the published model under its archived inputs;
they do not reproduce an unavailable MATLAB discretisation or uncertainty method.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .coefficients import Conditions, get
from .constants import SEC_PER_DAY
from .datasets import EXAMPLES_DIR
from .fitting import DiffusionModel, fit_time
from .solvers import InitialCondition, dirichlet

ROOT = EXAMPLES_DIR / 'validation'


def records():
    return json.loads((ROOT / 'manifest.json').read_text(encoding='utf8'))


def archived_model(record, *, coefficient=None, n_nodes=401):
    """Return a model and sorted table for a fully tabulated single-event case."""
    if record['study']=='lynn2024' and record['sample']=='Ol 8':
        raise ValueError('Ol 8 contains two events; a single-duration fit is not a replication.')
    s = record.get('setup')
    if not s or record['study'] != 'lynn2024':
        raise ValueError('No source-complete model setup for this profile; see the validation report.')
    df = pd.read_csv(ROOT / record['file']).dropna(
        subset=['Initial_Fo_mol','Published_model_Fo_mol']).sort_values('Distance_um')
    if not np.all(np.isfinite(df[['Distance_um','Fo_mol','Initial_Fo_mol']])):
        raise ValueError('Missing measured or initial profile values')
    x = df.Distance_um.to_numpy()
    if np.any(np.diff(x) <= 0):
        raise ValueError('Distances must be unique')
    initial = df.Initial_Fo_mol.to_numpy()
    cond = Conditions(T_K=s['T_C']+273.15,P_Pa=s['P_MPa']*1e6,
                      log_fo2_bar=s['log_fo2_bar'],angles_deg=tuple(s['angles_deg']),
                      X={'XFe':1-float(np.mean(initial))/100})
    model = DiffusionModel(get(coefficient or s['coefficient']),cond,
        InitialCondition('table',dict(x_table=x,C_table=initial)),
        bc_left=dirichlet(initial[0]),bc_right=dirichlet(initial[-1]),
        boundaries_far=False,comp_key='XFe',comp_scale=-.01,comp_offset=1,
        x_grid=np.linspace(x.min(),x.max(),n_nodes),n_nodes=n_nodes)
    return model,df


def compare(record, *, n_nodes=401):
    original,df = archived_model(record,n_nodes=n_nodes)
    x,y = df.Distance_um.to_numpy(),df.Fo_mol.to_numpy()
    archived = df.Published_model_Fo_mol.to_numpy()
    published = record['setup']['published_days']
    result = dict(key=record['key'],study=record['study'],sample=record['sample'],
                  published_days=published,published_sigma_days=record['setup']['published_sigma_days'],
                  n_nodes=n_nodes,classification='archived-input reconstruction')
    for label,model in [('original',original),('alternative',replace(original,coefficient=get('ol_FeMg_oeser2026')))]:
        # Same observations, weights, initial state, domain and boundaries.
        fit = fit_time(model,x,y,free_parameters=('t',),
                       t_min=SEC_PER_DAY*.001,t_max=SEC_PER_DAY*1e5)
        result[label+'_days']=fit.t_seconds/SEC_PER_DAY
        result[label+'_rmse_Fo']=float(np.sqrt(np.mean((fit.C_model-y)**2)))
        result[label+'_success']=bool(fit.success)
        result[label+'_warnings']=fit.warnings
        if label=='original':
            prediction=model.profile(published*SEC_PER_DAY,x)
            result['published_time_curve_rmse_Fo']=float(np.sqrt(np.mean((prediction-archived)**2)))
    result['original_to_published_ratio']=result['original_days']/published
    result['alternative_to_original_ratio']=result['alternative_days']/result['original_days']
    result['within_published_interval']=abs(result['original_days']-published)<=result['published_sigma_days']
    result['alternative_scope']='Oeser 2026 tracer-derived interdiffusion: sensitivity only; atmospheric-pressure, high-silica, ~1e-5 Pa fO2 calibration. No claim of greater precision.'
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--key',action='append',help='case key; repeat to select several')
    parser.add_argument('--nodes',type=int,default=401)
    parser.add_argument('--output',type=Path,default=ROOT/'results.json')
    args=parser.parse_args()
    selected=[r for r in records() if r['study']=='lynn2024' and r.get('setup') and r['sample']!='Ol 8']
    if args.key:
        selected=[r for r in selected if r['key'] in args.key]
        if len(selected)!=len(set(args.key)):
            parser.error('Unknown or unsupported case key')
    results=[]
    for record in selected:
        result=compare(record,n_nodes=args.nodes)
        results.append(result)
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(results,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
        print(f"{record['key']}: published {result['published_days']:.3g} d; "
              f"original {result['original_days']:.3g} d; alternative {result['alternative_days']:.3g} d",flush=True)


if __name__=='__main__':
    main()
