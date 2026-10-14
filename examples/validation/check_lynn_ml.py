"""Reconstruct the Lynn et al. (2024) Mauna Loa 2022 olivine timescales.

Paper setup (Nature Commun. 15:9451, Methods pp. 7-8): Fe-Mg interdiffusion with
the Dohmen & Chakraborty (2007) TaMED law along c, a and b six times slower,
projected on the EBSD angles; a step (Heaviside) initial condition at x_s, the
distance from the rim; rim and core held at Fo_r and Fo_c; Fo sigma 0.1 mol%.
The paper samples dt, T, p, fO2, x_s, Fo_r and Fo_c by MCMC (orientation fixed)
and reports posterior medians with equal-tailed 95% intervals (workbook sheet S6).

Diffusor holds T, P, log fO2 and the angles at the values printed on each
crystal's sheet (the centres of the paper's priors) and fits dt, x_s, Fo_r and
Fo_c by least squares to the points marked Valid = 1. The rim is at x = 0. This
is a best fit at fixed conditions, not a posterior median.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from diffusor.coefficients import Conditions, get
from diffusor.constants import SEC_PER_DAY
from diffusor.fitting import DiffusionModel, fit_time
from diffusor.solvers import InitialCondition, dirichlet

ROOT = Path(__file__).resolve().parent
LAW = 'ol_FeMg_dohmen_chakraborty2007_tamed'


def model_for(record, n_nodes=401):
    df = pd.read_csv(ROOT/record['file'])
    df = df[df.Valid == 1].sort_values('Distance_um')
    x, y = df.Distance_um.to_numpy(), df.Fo_mol.to_numpy()
    s, ts = record['source_conditions'], record['source_timescale']
    rim, core = float(ts['rim_Fo']), float(ts['core_Fo'])
    # start the step where the profile crosses half-way between rim and core
    half = 0.5 * (rim + core)
    x0 = float(np.interp(half, y, x)) if np.all(np.diff(y) > 0) else float(x[np.argmin(abs(y - half))])
    cond = Conditions(T_K=s['T_C'] + 273.15, P_Pa=s['P_MPa'] * 1e6, log_fo2_bar=s['log_fo2_bar'],
                      angles_deg=tuple(s['angles_deg']), X={'XFe': 1 - float(np.mean(y)) / 100})
    model = DiffusionModel(get(LAW), cond, InitialCondition('step', dict(x0=x0, C_left=rim, C_right=core)),
                           bc_left=dirichlet(rim), bc_right=dirichlet(core), boundaries_far=False,
                           comp_key='XFe', comp_scale=-.01, comp_offset=1,
                           x_grid=np.linspace(0.0, x.max(), n_nodes), n_nodes=n_nodes)
    return model, x, y


def refit(record):
    # 401 nodes resolve all but the shortest profile; a fit that ends on the edge of the time
    # window (a step sharper than the grid) is repeated on 801 nodes, which is converged
    # (ORI-4_ol18: 0.9121 d on 801, 0.9123 d on 1601 nodes)
    for n_nodes in (401, 801):
        model, x, y = model_for(record, n_nodes)
        fit = fit_time(model, x, y, sigma=np.full(y.size, 0.1), free_parameters=('t', 'x0', 'C_left', 'C_right'),
                       t_min=1e-3 * SEC_PER_DAY, t_max=1e5 * SEC_PER_DAY)
        if not any('edge of the search window' in w for w in fit.warnings):
            break
    ts = record['source_timescale']
    lo = ts['adjusted_min95_days'] + ts['day_erupted']
    hi = ts['adjusted_max95_days'] + ts['day_erupted']
    days = fit.t_seconds / SEC_PER_DAY
    return dict(key=record['key'], sample=record['sample'], location=ts['location'],
                published_median_days=ts['median_days'], published_95_days=[lo, hi],
                diffusor_days=days, ratio=days / ts['median_days'], within_95=bool(lo <= days <= hi),
                fitted=dict(x0_um=fit.free['x0'], Fo_rim=fit.free['C_left'], Fo_core=fit.free['C_right']),
                n_points=int(x.size), n_nodes=n_nodes, rmse_Fo=float(np.sqrt(np.mean((fit.C_model - y) ** 2))),
                success=bool(fit.success), warnings=fit.warnings)


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--key', action='append', help='recompute only these records and keep the saved others')
    args = parser.parse_args()
    records = [r for r in json.loads((ROOT/'manifest.json').read_text(encoding='utf8')) if r['study'] == 'lynn_ml2022']
    saved = {}
    if args.key and (ROOT/'lynn_ml_results.json').exists():
        saved = {v['key']: v for v in json.loads((ROOT/'lynn_ml_results.json').read_text(encoding='utf8'))['fits']}
    out = []
    for r in records:
        if args.key and r['key'] not in args.key and r['key'] in saved:
            out.append({**saved[r['key']], 'n_nodes': saved[r['key']].get('n_nodes', 401)})
            continue
        res = refit(r)
        out.append(res)
        print(f"{r['sample']:12s} published {res['published_median_days']:8.2f} d "
              f"[{res['published_95_days'][0]:7.2f}, {res['published_95_days'][1]:7.2f}]  "
              f"Diffusor {res['diffusor_days']:8.2f} d  ratio {res['ratio']:.2f}  {'in' if res['within_95'] else 'OUT'}",
              flush=True)
    ratios = np.array([r['ratio'] for r in out])
    result = dict(study='lynn_ml2022', coefficient=LAW,
                  classification='best fit at the sheet T, P, fO2 and angles; the paper reports MCMC posterior medians',
                  n=len(out), n_within_95=sum(r['within_95'] for r in out),
                  ratio_median=float(np.median(ratios)), ratio_min=float(ratios.min()), ratio_max=float(ratios.max()),
                  fits=out)
    (ROOT/'lynn_ml_results.json').write_text(json.dumps(result, indent=2, ensure_ascii=False) + '\n', encoding='utf8')
    print(f"{result['n_within_95']} of {len(out)} inside the 95% interval; ratio median {result['ratio_median']:.2f} "
          f"({result['ratio_min']:.2f}-{result['ratio_max']:.2f})")


if __name__ == '__main__':
    main()
