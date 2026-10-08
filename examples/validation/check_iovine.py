"""Refit Iovine et al. (2017) sanidine Ba traverses with the paper's own model.

Paper setup (Fig. 5e-f and p. 6-8): Ba diffusion with the Cherniak (2002) law,
D0 = 0.29 m2/s and Ea = 455 kJ/mol, no composition, fO2, pressure or orientation
term, at 930 C; a sharp initial step and the erfc solution with both plateaus
continuing (Fig. 5e). Diffusor fits the time, the step position and both
plateau values to every greyscale and X-ray traverse that Table 1 or 2 assigns
a time to. D is constant, so a linear grey-value or count-rate calibration
does not change the time and the values are fitted as recorded.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from diffusor.coefficients import Conditions, get
from diffusor.constants import SEC_PER_YEAR
from diffusor.fitting import DiffusionModel, fit_time
from diffusor.solvers import InitialCondition

ROOT = Path(__file__).resolve().parent


def model_for(record, df):
    s = record['setup']
    x = df.Distance_um.to_numpy()
    y = df[record['value_column']].to_numpy()
    order = np.argsort(x)
    x, y = x[order], y[order]
    # start the step at the steepest gradient between the end values
    i = int(np.argmax(np.abs(np.diff(y))))
    initial = InitialCondition('step', dict(x0=float(0.5*(x[i]+x[i+1])), C_left=float(y[0]), C_right=float(y[-1])))
    model = DiffusionModel(get(s['coefficient']), Conditions(T_K=s['T_C']+273.15), initial,
                           boundaries_far=True, n_nodes=401)
    return model, x, y


def refit(record):
    df = pd.read_csv(ROOT/record['file'])
    model, x, y = model_for(record, df)
    fit = fit_time(model, x, y, free_parameters=('t', 'x0', 'C_left', 'C_right'),
                   t_min=1e-3*SEC_PER_YEAR, t_max=1e5*SEC_PER_YEAR)
    years = fit.t_seconds/SEC_PER_YEAR
    s = record['setup']
    span = float(np.ptp(y))
    return dict(key=record['key'], sample=record['sample'], kind=record['value_column'],
                published_table=s['published_table'], published_label=f"{s['member']} {s['published_label']}",
                published_years=s['published_years'], diffusor_years=years,
                ratio=years/s['published_years'], n_points=len(x),
                rmse_relative_to_step=float(np.sqrt(np.mean((fit.C_model-y)**2))/span),
                success=bool(fit.success), warnings=fit.warnings)


def main():
    records = [r for r in json.loads((ROOT/'manifest.json').read_text(encoding='utf8'))
               if r['study'] == 'iovine2017' and r.get('setup')]
    out = [refit(r) for r in records]
    for r in out:
        print(f"{r['sample']:22s} {r['published_label']:16s} published {r['published_years']:5g} yr, "
              f"Diffusor {r['diffusor_years']:8.2f} yr ({r['ratio']:.2f})", flush=True)
    rounded = [r for r in out if round(r['diffusor_years']) == r['published_years']]
    result = dict(study='iovine2017', T_C=930.0, coefficient='kfs_Ba_cherniak2002',
                  classification='refit of the published greyscale and X-ray traverses with the paper\'s law, '
                                 'temperature, sharp initial step and free plateaus',
                  n=len(out), n_equal_after_rounding=len(rounded),
                  ratio_min=min(r['ratio'] for r in out), ratio_max=max(r['ratio'] for r in out),
                  ratio_median=float(np.median([r['ratio'] for r in out])), fits=out)
    (ROOT/'iovine_results.json').write_text(json.dumps(result, indent=2, ensure_ascii=False)+'\n', encoding='utf8')
    print(f"{len(rounded)} of {len(out)} round to the published year; ratios {result['ratio_min']:.2f}-{result['ratio_max']:.2f}")


if __name__ == '__main__':
    main()
