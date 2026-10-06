"""Build the reviewable report and figures from saved validation outputs."""
from collections import Counter
from pathlib import Path
import csv
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
load=lambda name:json.loads((ROOT/name).read_text(encoding='utf8'))


def main():
    records=load('manifest.json'); results=load('results.json'); studies=load('studies.json')
    counts=Counter(r['study'] for r in records)
    g=load('gordeychik_results.json'); opx=load('ostorero_results.json'); araya=load('araya_results.json')
    mourey=load('mourey_results.json')
    reunion=load('reunion_results.json'); eifel=load('eifel_results.json')
    published=np.array([r['published_days'] for r in results])
    fitted=np.array([r['original_days'] for r in results])
    alternative=np.array([r['alternative_days'] for r in results])
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(1,2,figsize=(11,4.8),layout='constrained')
    lo,hi=10,700
    xx=np.geomspace(lo,hi,100)
    axs[0].fill_between(xx,.7*xx,1.3*xx,color='#dde7e8',label='Published ±30% (approx.)')
    axs[0].plot(xx,xx,color='#7c898b',lw=1)
    axs[0].scatter(published,fitted,color='#1b6d77',s=32,zorder=3)
    for r in results:
        if r['sample'] in ('Ol 1','Ol 19','Ol 20'):
            axs[0].annotate(r['sample'],(r['published_days'],r['original_days']),xytext=(5,-12),textcoords='offset points')
    axs[0].set(xscale='log',yscale='log',xlim=(lo,hi),ylim=(lo,hi),xlabel='Published time (days)',ylabel='Diffusor, original D (days)',title='Lynn 2024: archived-input reconstruction')
    axs[0].legend(frameon=False,fontsize=8,loc='upper left')
    order=np.argsort(alternative/fitted)
    axs[1].axvline(1,color='#7c898b',lw=1)
    axs[1].scatter((alternative/fitted)[order],np.arange(len(results)),color='#b96336',s=30)
    axs[1].set(yticks=np.arange(len(results)),yticklabels=[results[i]['sample'] for i in order],
               xlabel='Time with Oeser 2026 / time with original D',title='Change D only: sensitivity, not greater precision')
    axs[1].grid(axis='x',alpha=.15)
    fig.savefig(ROOT/'comparison.png',dpi=180)
    fig.savefig(ROOT/'comparison.svg')
    plt.close(fig)
    fig,axs=plt.subplots(1,3,figsize=(12,3.7),layout='constrained')
    for ax,key in zip(axs,('lynn2024_ol_1','lynn2024_ol_19','lynn2024_ol_20')):
        df=pd.read_csv(ROOT/'fits'/f'{key}.csv')
        ax.plot(df.Distance_um,df.Initial_Fo_mol,color='#a0a0a0',ls=':',label='Archived initial')
        ax.scatter(df.Distance_um,df.Fo_mol,color='#303a40',s=15,label='Measured')
        ax.plot(df.Distance_um,df.Published_model_Fo_mol,color='#679869',lw=1.8,label='Published curve')
        ax.plot(df.Distance_um,df.Diffusor_original_fit_Fo,color='#1b6d77',label='Original-D refit')
        ax.plot(df.Distance_um,df.Diffusor_alternative_fit_Fo,color='#b96336',ls='--',label='Alternative-D refit')
        ax.set(title=key.replace('lynn2024_','').replace('_',' ').title(),xlabel='Distance (µm)',ylabel='Fo (mol%)')
    axs[0].legend(fontsize=7,frameon=False)
    fig.savefig(ROOT/'profiles.png',dpi=180)
    plt.close(fig)
    lines=['# Published-study validation — 7 October 2026','',
       f'{len(records)} measured profiles from eight downloaded studies are available in **File > Published validation library**. '
       'Eighteen Lynn profiles also have ready-to-run model presets in the usual example list. '
       'The remaining library entries are explicitly data-only. They require model setup before fitting.', '',
       '## What replicated','',
       f'**Lynn et al. (2024): all {len(results)} single-event reconstructions lie within the published uncertainty intervals.** '
       f'Refit/published time ratios range from {min(fitted/published):.3f} to {max(fitted/published):.3f}. '
       'This is agreement under the archived inputs, not bit-for-bit replication of the authors’ MATLAB discretisation or uncertainty calculation. '
       'Ol 8 has two events and is excluded from the single-duration comparison.','',
       f'**Changing only D changes the Lynn times by {(min(alternative/fitted)-1)*100:.1f}% to {(max(alternative/fitted)-1)*100:+.1f}%.** '
       'The alternative is Oeser, Dohmen & Weyer (2026), the registry’s interdiffusion law derived from its Fe and Mg tracer laws. '
       'It is not a tracer coefficient substituted directly into an exchange model. '
       'The calibration is at atmospheric pressure, approximately 10⁻⁵ Pa fO₂, high silica activity and 1100–1250 °C. '
       'The Lynn runs use 42 MPa and log fO₂ = −3.2 Pa, so the comparison extrapolates the calibration. '
       'Newer does not establish greater precision, and no narrower uncertainty interval is claimed.','',
       '![Published times and D sensitivity](comparison.png)','',
       '## Lynn setup and results','',
       'The workbook supplies T = 1200 °C, P = 42 MPa, log₁₀ fO₂ = −8.2 bar, three EBSD angles, measured Fo, Initial and Model columns, and published times. '
       'The article states 45 MPa; the workbook is used here. Both fits use the same data, uniform weights, plane geometry, fixed endpoint compositions, '
       'composition-dependent D and the archived initial values linearly interpolated between measurement positions. '
       'Fo is modelled in mol% with XFe = 1 − Fo/100. Only duration is fitted. '
       'Rows without an archived initial/model value are retained in the raw CSV but excluded from the model domain. '
       'No beam correction, additional free interface, inferred plateau or cooling path is added. '
       'The sub-grid sharp interface and original numerical grid are not published. '
       'The paper writes Fick’s second law as ∂C/∂t = D ∂²C/∂x² (eq. 1) and makes D depend on X_Fe (eq. 2); '
       'Diffusor solves ∂C/∂t = ∂/∂x(D ∂C/∂x) with its finite-volume solver on 401 nodes. '
       'The source model curves are therefore also compared directly. '
       'The workbook’s log₁₀ fO₂ = −8.2 bar is used as printed; the paper’s Methods name QFM+0.4, which Diffusor’s QFM buffer puts at −7.87 bar at 1200 °C and 42 MPa.','',
       '| Crystal | Published ± interval (d) | Original-D refit (d) | Alternative-D refit (d) | Alternative/original |',
       '|---|---:|---:|---:|---:|']
    for r in results:
        lines.append(f"| {r['sample']} | {r['published_days']:g} ± {r['published_sigma_days']:g} | {r['original_days']:.2f} | {r['alternative_days']:.2f} | {r['alternative_to_original_ratio']:.3f} |")
    lines+=['',f"The largest difference from the archived model curve at its published time is {max(r['published_time_curve_rmse_Fo'] for r in results):.3f} Fo mol% RMS. "
             'Ol 20 has the largest time discrepancy (about 25% shorter). It remains a discrepancy even though it is inside the reported interval. '
             'Doubling the grid to 801 nodes for Ol 1 and Ol 9 changes fitted times by less than 0.004%. '
             'This checks numerical convergence for those cases, not the unknown source initial-interface placement.','',
             '![Measured, initial and fitted profiles](profiles.png)','',
             'Full curves: [fits/](fits/). Machine-readable fit results and calibration warnings: [results.json](results.json). '
             'Grid checks: [convergence.json](convergence.json).','',
             '## Other executable checks','',
             '**Gordeychik et al. (2018).** All 32 published Fo age rows in Tables SM4-A, SM5-A and SM6-A reproduce to relative error below 4×10⁻¹⁶. '
             'These checks use the published fitted diffusion widths or Dt products and geometric factors. '
             'They verify the conversion into ages, not a new fit to the raw profiles or the jointly inferred growth/resorption geometry. '
             'The original registry D is about 0.3% below the spreadsheet D, largely reflecting literal constants (2.303 versus ln 10 and the gas constant). '
             'The reported width checks use the source D rather than hiding that distinction.','',
             'Equation SM4.4 of the supplement gives t = Δ²_Ni/(4D_Ni). In Table SM4-A, column O (higher D_Ni) computes this, '
             'while cells P4:P8 (lower D_Ni) multiply ΔFo × ΔNi. Diffusor reads P4:P8 as a slip in the spreadsheet formula and follows equation SM4.4, '
             'which gives values 7.4–7.6% lower in those five cells. Column P does not enter the table’s resulting time intervals, which take columns O and L. '
             'The spreadsheet formulas and the equation SM4.4 values are both kept in the results file.','']
    ratios=[v for r in g for v in r['alternative_to_published_ratio']]
    lines += [f'Where an unambiguous a/b/c orientation match exists, changing D to Oeser gives {min(ratios):.3f}–{max(ratios):.3f} times the published Fo age. '
              'The new anisotropy is projected separately; the old sixfold geometric factor is not reused. '
              'An orientation is paired with an age row only when the sample name and the geometric factor both agree with Table SM2-B. '
              f"{sum(1 for r in g if not r['alternative_to_published_ratio'])} of {len(g)} rows have no such match "
              '(for example, the age row for Ol-8-2 uses the factor 0.2539 that SM2-B lists for Ol-8-3). '
              'Those rows keep the age check but have no alternative-D result. '
              'At 0.6–1 GPa, these alternative coefficients are pressure extrapolations. '
              '[Source parameters and literal formulas](gordeychik_parameters.json), [results](gordeychik_results.json).','',
              f"**Ostorero et al. (2022).** The existing K9_L10C4 EPMA example gives {opx['fits'][0]['years']:.2f} yr with Ganguly & Tazzoli (1994) and {opx['fits'][1]['years']:.2f} yr with Dias et al. (2025), versus 2.32 yr published. "
              'This is an EPMA surrogate for the authors’ BSE fit, with the existing window, 850 °C, NNO+1.3, b-axis and beam correction held fixed. '
              'The newer-law result includes temperature/fO₂ extrapolation. [Results](ostorero_results.json).','',
              f"**Araya et al. (2024).** The stated log₁₀D = −19.78 is recovered as {araya['recomputed_logD']:.5f} with the original law at 966 °C, NNO, a-axis and no composition correction. "
              'The paper specifies perpendicular to c, which does not uniquely identify a or b. '
              'The 1-bar buffer reference is a computational assumption. Alternative-D factors are conditional on orientation and composition, not new ages. '
              '[Assumptions and D-only comparison](araya_results.json).','',
              '**Mourey et al. (2023), 3311_2_ol6.** A conditional Fo-only reconstruction gives '
              f"{mourey['fits'][1]['days']:.1f} days with the Dohmen & Chakraborty (2007) law as corrected by its erratum, versus 395 (+149/−106) days published, "
              f"and {mourey['fits'][2]['days']:.1f} days with Oeser D. "
              'Equation 3 of the paper gives fO₂ = 10^(8.912−25160/T) “in Pa”, for the QFM buffer named in the text. '
              'Read in bar, the expression is −8.37 at 1183 °C, 0.09 log units from Diffusor’s QFM (−8.45 bar at 60 MPa); read in Pa it lies about 5 log units below QFM. '
              'Diffusor reads it in bar. Equation 2 prints the composition term as `3(X_Mg−0.9)`; Dohmen & Chakraborty (2007, erratum) give `3(X_Fe−0.1)`, which Diffusor uses. '
              f"Taking both printed forms at face value gives {mourey['fits'][0]['days']:.0f} days in this reconstruction. "
              'The 406-day result rests on these two readings. It does not establish what the authors’ code used. '
              'The initial rim position is fitted to Fo alone here, whereas the paper used Ca and Ni as additional constraints. '
              '[All three branches and assumptions](mourey_results.json).','',
              '## Coverage and remaining requirements','',
              '| Study | Extracted profiles | What the downloaded material permits |','|---|---:|---|']
    section=['## Sundermeyer reconstructions: published ages not reproduced','',
        'These runs are executable sensitivity examples, not verified GUI presets. '
        'The original point selection and parts of the setup are unavailable. '
        'Réunion uses the available TaMED law as a proxy for the cited Chakraborty (2010) implementation. '
        'Eifel uses an explicitly assumed log fO₂ = −5 Pa because the downloaded source does not state it.','',
        '| Study / crystal | Published (d) | TaMED reconstruction (d) | Oeser comparison (d) |',
        '|---|---:|---:|---:|']
    for r in reunion:
        section.append(f"| Réunion / {r['sample']} | {r['published_days']:g} | {r['fits']['source_law']['fit_days']:.2f} | {r['fits']['alternative_law']['fit_days']:.2f} |")
    section.append(f"| Eifel / {eifel['case']} | {eifel['published']['days']:g} | {eifel['fits'][0]['days']:.2f} | {eifel['fits'][1]['days']:.2f} |")
    section += ['',
        'Réunion 150915-1-9 has a matching melt-inclusion temperature but a reconstructed eight-point window. '
        'Its TaMED time lies below the published interval. Extending its numerical domain changes the fit by only 0.28%, '
        'so the mismatch cannot be removed by that domain correction. '
        '150915-1-3 also needs a proxy temperature and gives a particularly poor reproduction. '
        'Eifel’s grid refinement changes the fitted time by less than 0.02%, but its unknown fO₂ and reconstructed window remain scientific limitations. '
        'These mismatches do not establish that the published ages are wrong or that Oeser is more accurate.','',
        '[Réunion methods and limitations](reunion_notes.md), [results](reunion_results.json). '
        '[Eifel methods and limitations](eifel_notes.md), [results](eifel_results.json).','']
    index=lines.index('## Coverage and remaining requirements')
    lines[index:index]=section
    for s in studies:
        link='https://doi.org/'+s.get('doi',s.get('data_doi',''))
        lines.append(f"| [{s['citation']}]({link}) | {counts[s['key']] if counts[s['key']] else '—'} | {s['status']} |")
    lines+=['','A workbook of core/rim compositions is not a measured diffusion traverse. '
            'A case requiring different dimensions, multiple events or a Bayesian joint inversion has not been labelled an exact replication merely because a scalar curve can be fitted. '
            'The supplements of Weller et al. (2026), Lynn et al. (2024, Mauna Loa) and Kahl et al. (2023) have not been downloaded yet.','',
            '## Reproduce','',
            'Run from the repository root in the project environment:','',
            '```text','python examples/validation/extract.py','python -m diffusor.validation',
            'python examples/validation/extract_gordeychik.py','python examples/validation/check_gordeychik.py',
            'python examples/validation/check_ostorero.py','python examples/validation/check_araya.py',
            'python examples/validation/check_mourey.py',
            'python examples/validation/check_reunion.py','python examples/validation/check_eifel.py',
            'python examples/validation/export_curves.py','python examples/validation/build_report.py','```','',
            'The extractors read the original supplementary workbooks, downloaded from the DOIs above into `papers/Supplementaries and data/<study>/`, '
            'and need openpyxl and xlrd (Ruth XLS only). '
            'Install the optional extraction dependency with `python -m pip install ".[validation]"`. '
            'Fitting and checks run offline from the bundled extracted files. '
            'For the 801-node checks use `python -m diffusor.validation --key lynn2024_ol_1 --key lynn2024_ol_9 --nodes 801 --output examples/validation/convergence.json`.','',
            'Each CSV retains `source_row`; the [manifest](manifest.json) stores source file, sheet, SHA256, sample and setup status. '
            'Original order and repeated distances are preserved. Review repeated positions before using a data-only file in a model. '
            'Mutch Fo and its uncertainty were converted from mole fraction to mol%; Gordeychik distances from mm to µm. '
            'No synthetic points or interpolated measurements were inserted. '
            '[Browse every extracted profile](catalogue.md). '
            '[Independent extraction audit](data_audit.md) checked all 11,758 source rows.','']
    lines += ['`tests/test_validation.py` rechecks the manifest, the Lynn model setup, the GUI preset and the Gordeychik arithmetic on every test run. '
              'Wheel contents and loading were verified in an isolated installation. '
              '[Mutch input requirements](mutch_notes.md) explain why its joint inversion has not been replaced by an assumed single-species fit.','']
    (ROOT/'report.md').write_text('\n'.join(lines),encoding='utf8')
    catalogue=['# Extracted measured profiles','','These are measurements. Only entries marked **preset** have a bundled model setup.','',
               '| Study / crystal | Points | Setup | CSV |','|---|---:|---|---|']
    for r in records:
        ready=r['study']=='lynn2024' and r.get('setup') and r['sample']!='Ol 8'
        catalogue.append(f"| {r['study']} / {r['sample']} | {r['n_points']} | {'preset' if ready else 'data only'} | [profile]({r['file']}) |")
    (ROOT/'catalogue.md').write_text('\n'.join(catalogue)+'\n',encoding='utf8')


if __name__=='__main__':
    main()
