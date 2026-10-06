"""Reconstructed one-boundary Fe-Mg comparison for Sundermeyer et al. (2020).

Run from the repository root with .venv/Scripts/python.exe examples/validation/check_eifel.py.
This is explicitly a reconstruction/sensitivity, not an exact DIPRA replication.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from diffusor.coefficients import Conditions, get
from diffusor.constants import SEC_PER_DAY
from diffusor.fitting import DiffusionModel, fit_time
from diffusor.solvers import InitialCondition, dirichlet

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
PROFILE = ROOT / 'profiles' / 'sundermeyer2020_eifel_e41_4_1.csv'

# Published Table 1 age for the Fo/Mg-Fe fit; uncertainty columns in the paper
# are asymmetric (central value, lower, upper).
PUBLISHED = {'days': 291.0, 'minus_days': 213.0, 'plus_days': 103.0}
# The source paper gives no oxygen fugacity. The alternative law is calibrated
# near 1e-5 Pa, so use that common reference only to compare laws consistently.
COMPARISON_LOG_FO2_PA = -5.0
PRESSURE_PA = 2.0e8  # paper: minimum 2 kbar; similar estimate for Rothenberg.
FIT_MIN_UM = 30.0    # exclude the outer low-Fo rim/mantle gradient; fit inner boundary.
N_NODES = 121
REFINED_NODES = 241


def supplement_row():
    # Transcribed from Online Resource 4, row 21. Kept here so this shipped
    # check runs from the repository without the separately downloaded paper files.
    return {
        'source_row': 21,
        'crystal': 'E41-4-1',
        'T_C': 1159.0,
        'T_sigma_C': 23.0,
        'T_effective_C': 1101.0,
        'Dx_um': 1.49254,
        'Dt_s': 80000.0,
        'angles_deg_to_abc': [81.9, 9.2, 85.6],
        'Fo_core_mol_fraction': 0.838,
        'Fo_rim_mol_fraction': 0.859,
        'Fo_GB_mol_fraction': '-',
    }


def fit(coefficient_key, params, df, log_fo2_pa, n_nodes):
    # x increases from the more magnesian mantle toward the lower-Fo core in
    # this selected window; this orientation determines left/right step values.
    x = df['Distance_um'].to_numpy(float)
    fo = df['Fo_mol'].to_numpy(float) / 100.0
    initial = InitialCondition('step', {
        'x0': float(np.median(x)),
        'C_left': params['Fo_rim_mol_fraction'],
        'C_right': params['Fo_core_mol_fraction'],
    }, description='OR4 core/mantle Fo step; outer rim excluded')
    cond = Conditions(
        T_K=params['T_effective_C'] + 273.15,
        P_Pa=PRESSURE_PA,
        log_fo2_bar=log_fo2_pa - 5.0,
        X={'XFe': 1.0 - float(np.mean(fo))},
        angles_deg=tuple(params['angles_deg_to_abc']),
    )
    model = DiffusionModel(
        coefficient=get(coefficient_key), conditions=cond, initial=initial,
        comp_key='XFe', comp_scale=-1.0, comp_offset=1.0,
        bc_left=dirichlet(initial.params['C_left']),
        bc_right=dirichlet(initial.params['C_right']),
        boundaries_far=False,
        x_grid=np.linspace(float(x.min()), float(x.max()), n_nodes),
        n_nodes=n_nodes,
    )
    res = fit_time(model, x, fo, sigma=np.full_like(fo, 0.2 / 100.0),
                   free_parameters=('t', 'x0'),
                   t_min=SEC_PER_DAY * 0.01, t_max=SEC_PER_DAY * 5000)
    return {
        'coefficient': coefficient_key,
        'n_nodes': n_nodes,
        'days': float(res.t_seconds / SEC_PER_DAY),
        'fitted_interface_um': float(res.free.get('x0', initial.params['x0'])),
        'rmse_Fo_mol_percent': float(100 * np.sqrt(np.mean((res.C_model - fo) ** 2))),
        'success': bool(res.success),
        'solver': res.route,
        'warnings': list(res.warnings),
    }


def main():
    params = supplement_row()
    with PROFILE.open(newline='', encoding='utf8') as f:
        all_rows = list(csv.DictReader(f))
    df = [r for r in all_rows if float(r['Distance_um']) >= FIT_MIN_UM]
    import pandas as pd
    df = pd.DataFrame(df)
    fits_by_grid = {}
    for n_nodes in (N_NODES, REFINED_NODES):
        fits_by_grid[str(n_nodes)] = [
            fit('ol_FeMg_dohmen_chakraborty2007_tamed', params, df, COMPARISON_LOG_FO2_PA, n_nodes),
            fit('ol_FeMg_oeser2026', params, df, COMPARISON_LOG_FO2_PA, n_nodes),
        ]
    fits = fits_by_grid[str(REFINED_NODES)]
    for grid_fits in fits_by_grid.values():
        for f in grid_fits:
            f['ratio_to_published'] = f['days'] / PUBLISHED['days']
    convergence = []
    for base, refined in zip(fits_by_grid[str(N_NODES)], fits_by_grid[str(REFINED_NODES)]):
        convergence.append({
            'coefficient': base['coefficient'],
            'nodes_low': N_NODES,
            'nodes_high': REFINED_NODES,
            'days_low': base['days'],
            'days_high': refined['days'],
            'relative_time_change': refined['days'] / base['days'] - 1.0,
            'interface_low_um': base['fitted_interface_um'],
            'interface_high_um': refined['fitted_interface_um'],
        })
    result = {
        'study': 'Sundermeyer et al. (2020), Eifel; DOI 10.1007/s00410-020-01715-y',
        'case': 'E41-4-1',
        'classification': 'single-boundary reconstruction and D-law sensitivity; not an exact DIPRA replication',
        'published': {**PUBLISHED, 'source': 'paper Table 1, p. 16; Mg-Fe column'},
        'source_profile': 'profiles/sundermeyer2020_eifel_e41_4_1.csv',
        'fit_window': {'minimum_distance_um': FIT_MIN_UM, 'maximum_distance_um': max(float(v) for v in df.Distance_um),
                       'n_points': int(len(df)), 'rationale': 'exclude the outer low-Fo rim and fit the inner core-mantle boundary'},
        'source_model_parameters': {**params, 'source': 'Online Resource 4, row for E41-4-1'},
        'assumptions': {
            'pressure_Pa': PRESSURE_PA,
            'pressure_basis': 'Paper p. 12 gives 2 kbar for Laacher See and says a similar estimate is used for Rothenberg and Eppelsberg.',
            'temperature_C': params['T_effective_C'],
            'temperature_basis': 'Online Resource 4 effective temperature; isothermal case.',
            'fO2_log10_Pa': COMPARISON_LOG_FO2_PA,
            'fO2_basis': 'Not reported for this Eifel case; selected solely as the Oeser et al. nominal calibration reference (1e-5 Pa), for a same-condition law comparison.',
            'initial_condition': 'Sharp Fo step using OR4 core and rim values (0.838/0.859), with mantle on the lower-distance side and core on the higher-distance side; interface is fitted.',
            'geometry': '1-D plane; source paper says DIPRA uses finite 1-D numerical modeling but does not specify a geometry correction for this profile.',
            'boundaries': 'Dirichlet plateau values at the cropped window edges; this is an explicit reconstruction choice, not documented DIPRA configuration.',
            'fit_method': 'Diffusor fit_time; uniform 0.2 mol% Fo uncertainty as reported 2-sigma analytical precision; optimize time and interface position.',
        },
        'fits': fits,
        'grid_comparison': convergence,
        'grid_fits': fits_by_grid,
        'alternative_scope': 'Oeser et al. (2026) interdiffusion sensitivity only. It is outside its reported composition (San Carlos XFe=0.085) and pressure (atmospheric) calibration at this case; warnings are retained.',
        'interpretation': 'The reconstructed original-law fit is not expected to equal the published DIPRA age because fO2, exact fit window, full DIPRA initial/boundary setup, and original solver implementation are unavailable. Do not use this case as a GUI preset.',
    }
    (ROOT / 'eifel_results.json').write_text(json.dumps(result, indent=2, ensure_ascii=False) + '\n', encoding='utf8')
    print('Published:', PUBLISHED['days'], 'days')
    for n_nodes, grid_fits in fits_by_grid.items():
        for f in grid_fits:
            print(n_nodes, f['coefficient'], round(f['days'], 3), 'days;', 'RMSE', round(f['rmse_Fo_mol_percent'], 4), 'mol%')


if __name__ == '__main__':
    main()
