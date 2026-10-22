"""Plot the complete audited source family. No individual data are accepted."""
from pathlib import Path
import datetime
import hashlib
import json
import math
import os
import re

from summarize_LN_specificity_source_aggregate_v1 import (
    ENDPOINTS, DIRECT, NULL, finite, validate_family,
)


ENDPOINT_LABELS = (
    "E0: PD1 / TIM3\nrelative-high CD8a-like",
    "E1: MKI67\nrelative-high CD8a-like",
    "E2: CCR7 / CD127\nrelative-high CD8a-like",
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def checked_source_report(task):
    """Require an audited actual exit and bind every plotted field to it."""
    source_path = task / 'private_LN_specificity_complete_source_v1.json'
    exit_path = task / 'private_LN_specificity_actual_exit_v1.json'
    report_path = task / 'private_LN_specificity_source_evidence_summary_v1.json'
    assert all(p.exists() for p in (source_path, exit_path, report_path)), 'ACTUAL_AUDITED_RESULTS_NOT_READY'
    source, exited, report = [json.loads(p.read_text()) for p in (source_path, exit_path, report_path)]
    frozen = json.loads((task / 'private_frozen_LN_specificity_execution_v1.json').read_text())
    assert exited['actual_exit_code'] == source['source_R_exit_code'] == 0
    assert exited['source_result_SHA256'] == report['actual_source_result_SHA256'] == digest(source_path)
    assert source['source_payload_RAM_released'] and source['no_dynamic_followup_nodes']
    assert source['source_object_SHA256'] == frozen['expected_file_SHA256']
    assert source['R_code_SHA256'] == frozen['R_code_SHA256']
    family = frozen['analysis_config']['family_contract_SHA256']
    assert source['frozen_secondary_family_SHA256'] == report['family_contract_SHA256'] == family
    summary = source['summary']
    assert summary['status'] == 'COMPLETE_FROZEN_SECONDARY_SPECIFICITY_AND_FULL_NUMERIC_AUDIT'
    production, audit = summary['production'], summary['audit']
    assert audit['status'] == 'PASS_FULL_SECONDARY_SPECIFICITY_NUMERIC_AUDIT'
    assert production['family_contract_SHA256'] == audit['family_contract_SHA256'] == family
    for rendered, native in (
        ('primary_rows', 'primary_rows'), ('scale_point_checks', 'sensitivities'),
        ('edge_margin_quality', 'edge_margin_quality'), ('coverage_reasons', 'coverage_reasons'),
        ('local_exchangeability_coverage', 'mobility'),
    ):
        assert report[rendered] == production[native], f'PLOTTED_FIELD_NOT_BOUND: {rendered}'
    validate_family(production)
    assert report['status'] == 'AUDITED_SOURCE_SECONDARY_EVIDENCE_SUMMARY'
    return report, report_path


def ratio(value):
    assert finite(value) and -700 < value < 700, 'NONFINITE_OR_UNPLOTTABLE_EFFECT'
    return math.exp(value)


def qtext(value):
    return f'{value:.3g}' if finite(value) else 'NA'


def ratio_axis(ax, log_values, label):
    """Keep the null visible; do not clip a point outside its bootstrap CI."""
    from matplotlib import ticker
    values = [0.0] + [v for v in log_values if finite(v)]
    low, high = min(values), max(values)
    span = max(high - low, .6)
    ax.set_xscale('log')
    left, right = ratio(low - .14 * span), ratio(high + .14 * span)
    ax.set_xlim(left, right)
    locator = (ticker.MaxNLocator(nbins=4) if right / left < 10 else
               ticker.LogLocator(base=10, subs=(1, 2, 5), numticks=6))
    ticks = sorted(set([1.0] + [float(x) for x in locator.tick_values(left, right) if left <= x <= right]))
    ax.set_xticks(ticks, [f'{x:.3g}' for x in ticks], fontsize=8)
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.axvline(1, color='#888888', ls='--', lw=.8)
    ax.set_xlabel(label, fontsize=9)
    ax.grid(axis='x', alpha=.15)


def plot_row(ax, row, y, values, include_interval=True, color='#2166ac'):
    point = row.get('median_effect')
    supported = row.get('status') == 'ESTIMABLE' and finite(row.get('primary_p'))
    marker = 'o' if supported else 'x'
    if finite(point):
        values.append(point)
        filled = supported and finite(row.get('BY')) and row['BY'] < .05
        ax.plot(ratio(point), y, marker=marker, ms=6, mec=color,
                mfc=color if filled else 'white', color=color)
        lower, upper = row.get('bootstrap_lower'), row.get('bootstrap_upper')
        if include_interval and finite(lower) and finite(upper):
            assert lower <= upper
            values.extend((lower, upper))
            ax.hlines(y, ratio(lower), ratio(upper), color=color, lw=1.4)
            ax.plot([ratio(lower), ratio(upper)], [y, y], '|', color=color, ms=7)
    else:
        ax.text(.5, y, 'not estimable', ha='center', va='center',
                transform=ax.get_yaxis_transform(), color='#777777', fontsize=8)
    suffix = '' if supported else '\n' + row['status'].replace('_', ' ')
    detail = f"n={row['cases']}; images={row['eligible_images']}\nBH={qtext(row.get('BH'))}; BY={qtext(row.get('BY'))}{suffix}"
    ax.text(.5, y + .29, detail, ha='center', va='top', fontsize=7.3,
            transform=ax.get_yaxis_transform())


def render_report(report, out, synthetic=False):
    assert report['status'] == ('SYNTHETIC_SOURCE_AGGREGATES_FOR_LAYOUT_ONLY' if synthetic
                                else 'AUDITED_SOURCE_SECONDARY_EVIDENCE_SUMMARY')
    records = validate_family({'primary_rows': report['primary_rows']})
    scales = {(r['endpoint'], r['grid']): r for r in report['scale_point_checks']}
    assert set(scales) == {(e, g) for e in ENDPOINTS for g in (2, 8)}
    assert len(report['scale_point_checks']) == 6
    for endpoint in ENDPOINTS:
        assert len({(records[endpoint, est]['eligible_images'], records[endpoint, est]['cases'])
                    for est in DIRECT}) == 1, 'DIRECT_SET_COUNTS_DIFFER'
        direct = records[endpoint, 'combined_MH']
        for row in (records[endpoint, NULL], scales[endpoint, 2], scales[endpoint, 8]):
            assert row['eligible_images'] <= direct['eligible_images'] and row['cases'] <= direct['cases'], 'NULL_SET_NOT_SUBSET_OF_DIRECT_COVERAGE'
    assert not out.exists(), 'PRESERVE_EXISTING_FIGURES_USE_NEW_VERSION'
    repository = Path(__file__).resolve().parents[1]
    assert out != repository and repository not in out.parents
    out.mkdir(parents=True)
    os.environ.setdefault('MPLCONFIGDIR', '/tmp/pancancer_LN_specificity_mpl_20261022')
    os.environ.setdefault('XDG_CACHE_HOME', '/tmp/pancancer_LN_specificity_fontcache_20261022')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9,
                         'axes.spines.top': False, 'axes.spines.right': False})
    prefix = 'SYNTHETIC LAYOUT ONLY — ' if synthetic else ''
    figure_paths = []

    def save(fig, stem):
        for extension in ('png', 'svg'):
            path = out / f'{stem}.{extension}'
            fig.savefig(path, dpi=180)
            figure_paths.append(path)
        plt.close(fig)

    fig, axes = plt.subplots(1, 4, figsize=(18, 7.2), sharey=True)
    titles = ('MKI67 remaining', 'MKI67 relative-high', 'Combined across\nMKI67 strata',
              'High / remaining\nodds-ratio ratio')
    for ax, estimand, title in zip(axes, DIRECT, titles):
        values = []
        for y, endpoint in enumerate(ENDPOINTS):
            plot_row(ax, records[endpoint, estimand], y, values)
        ratio_axis(ax, values, 'Odds-ratio ratio' if estimand == DIRECT[-1] else 'Odds ratio')
        ax.set_title(title, fontsize=11)
        ax.set_ylim(2.67, -.45)
    axes[0].set_yticks(range(3), ENDPOINT_LABELS)
    fig.suptitle(prefix + 'CD4-like regulatory candidate versus other CD4-like cells', fontsize=13)
    fig.text(.07, .035,
        'Each CD8a-like endpoint is contrasted with its own CD8a-like complement. MKI67 groups are relative marker groups, not validated cycling/quiescent cells.\n'
        'Points: median registered-patient-unit effect. Lines: 95% percentile patient bootstrap (2000 draws), conditional on fixed labels.\n'
        'BH/BY cover the complete 15-test secondary family. Filled circles: BY < 0.05; crosses: no primary inference. Missing tests remain visible.\n'
        'Only the direct high/remaining contrast tests a between-stratum difference. Different P values do not establish different effects.', fontsize=8.2)
    fig.tight_layout(rect=(.01, .20, .99, .94))
    save(fig, 'DIRECT_COMPARISONS')

    fig, axes = plt.subplots(1, 2, figsize=(12.8, 7.4), sharey=True,
                             gridspec_kw={'width_ratios': [1.05, 1.5]})
    values = []
    for y, endpoint in enumerate(ENDPOINTS):
        plot_row(axes[0], records[endpoint, NULL], y, values, color='#7b3294')
    ratio_axis(axes[0], values, 'Centered enrichment factor')
    axes[0].set_title('Primary 4 × 4 local grid\nConditional randomization test', fontsize=11)
    scale_colors = {2: '#1b9e77', 4: '#7b3294', 8: '#d95f02'}
    values = []
    for y, endpoint in enumerate(ENDPOINTS):
        for g, offset in ((2, -.18), (4, 0), (8, .18)):
            row = records[endpoint, NULL] if g == 4 else scales[endpoint, g]
            point = row.get('median_effect')
            if finite(point):
                values.append(point)
                axes[1].plot(ratio(point), y + offset, marker={2: 's', 4: 'o', 8: '^'}[g],
                             color=scale_colors[g], ms=5)
            # Each grid gets its own denominator; no connectors imply a matched comparison.
            axes[1].text(1.02, y + offset,
                f"{g}×{g}: n={row['cases']}, images={row['eligible_images']}" +
                ('' if row.get('status') == 'ESTIMABLE' else ' [no inference]'),
                transform=axes[1].get_yaxis_transform(), va='center', fontsize=7.8)
    ratio_axis(axes[1], values, 'Centered enrichment factor (point checks only)')
    axes[1].set_title('Grid sensitivity: separate eligible sets\nNo cross-grid test or attenuation percentage', fontsize=11)
    for ax in axes:
        ax.set_ylim(2.65, -.5)
    axes[0].set_yticks(range(3), ENDPOINT_LABELS)
    fig.suptitle(prefix + 'State pairing beyond the specified local-composition null', fontsize=13)
    fig.text(.07, .03,
        'The factor is exp(median patient-centered log effect), not the uncentered Mantel–Haenszel odds ratio.\n'
        'Shuffle blocks preserve local regulatory counts, MKI67 strata and CD8 degree; parent geometry and receiver states stay fixed.\n'
        'Fine-scale sender spatial autocorrelation is not fully preserved. Bootstrap intervals condition on the fitted null center.\n'
        'BH/BY cover all 15 secondary tests. Filled circles: BY < 0.05; crosses: no primary inference. Auxiliary grids have no P test.\n'
        'No claim of direct signalling, suppression, treatment resistance or clinical prediction follows from these tests.', fontsize=8.2)
    fig.tight_layout(rect=(.01, .22, .82, .93))
    save(fig, 'LOCAL_NULL_AND_SCALE')

    columns = ('Direct', 'Null 2×2', 'Null 4×4', 'Null 8×8')
    selected = [[records[e, 'combined_MH'], scales[e, 2], records[e, NULL], scales[e, 8]] for e in ENDPOINTS]
    fig, axes = plt.subplots(1, 3, figsize=(15.8, 7.4), gridspec_kw={'width_ratios': [1.2, 1.2, 1.05]})
    for ax, field, title in zip(axes[:2], ('eligible_images', 'cases'), ('Eligible images', 'Registered-patient units')):
        z = np.array([[r[field] for r in row] for row in selected], dtype=int)
        assert np.all(z >= 0)
        maximum = max(1, int(z.max()))
        im = ax.imshow(z, aspect='auto', cmap='Blues', vmin=0, vmax=maximum)
        for y in range(3):
            for x in range(4):
                ax.text(x, y, str(z[y, x]), ha='center', va='center',
                        color='white' if z[y, x] > maximum / 2 else 'black')
        ax.set_xticks(range(4), columns, rotation=25, ha='right')
        ax.set_yticks(range(3), ENDPOINT_LABELS if field == 'eligible_images' else ['E0', 'E1', 'E2'])
        ax.set_title(title, fontsize=11)
        fig.colorbar(im, ax=ax, fraction=.035, pad=.04)
    quality = {r['endpoint']: r for r in report['edge_margin_quality']}
    assert set(quality) == set(ENDPOINTS)
    zero_counts = [quality[e]['images_with_any_zero_edge_margin_stratum'] for e in ENDPOINTS]
    totals = [quality[e]['eligible_images'] for e in ENDPOINTS]
    assert all(0 <= n <= d for n, d in zip(zero_counts, totals))
    axes[2].barh(range(3), zero_counts, color='#bdbdbd')
    for y, (number, total) in enumerate(zip(zero_counts, totals)):
        axes[2].text(number + .04 * max(1, max(zero_counts)), y, f'{number} / {total}', va='center', fontsize=8)
    axes[2].set_xlim(0, max(1, max(zero_counts)) * 1.65)
    axes[2].set_yticks(range(3), ['E0', 'E1', 'E2'])
    axes[2].invert_yaxis()
    axes[2].set_title('Images with any zero\ncross-parent edge margin', fontsize=11)
    axes[2].set_xlabel('Images; continuity correction retained')
    mobility = {r['grid']: r for r in report['local_exchangeability_coverage']}
    assert set(mobility) == {2, 4, 8}
    mobility_text = ' | '.join(f"{g}×{g}: {mobility[g]['images_with_mobility']} / {mobility[g]['images_evaluated']} images" for g in (2, 4, 8))
    reasons = '; '.join(f'{k}: {v}' for k, v in sorted(report['coverage_reasons'].items()))
    fig.suptitle(prefix + 'Coverage and prespecified quality flags', fontsize=13)
    fig.text(.07, .04,
        'Exchangeability gate among evaluated images: ' + mobility_text + '\n'
        'Eligibility reason counts are image × endpoint opportunities, not unique images:\n' + reasons + '\n'
        'Missing eligibility is not a negative interaction. Equal counts do not prove that eligible image sets are identical.',
        fontsize=7.9, wrap=True)
    fig.tight_layout(rect=(.01, .24, .99, .94))
    save(fig, 'COVERAGE_AND_QUALITY')
    manifest = {
        'status': 'SYNTHETIC_LAYOUT_RENDERED_NOT_SCIENCE' if synthetic else 'ACTUAL_SOURCE_FIGURES_CREATED_PENDING_VISUAL_REVIEW',
        'observed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'family_contract_SHA256': report['family_contract_SHA256'],
        'actual_source_result_SHA256': None if synthetic else report['actual_source_result_SHA256'],
        'renderer_SHA256': digest(Path(__file__)),
        'figures': [{'name': p.name, 'SHA256': digest(p)} for p in figure_paths],
        'no_individual_arrays_accepted': True,
        'no_matched_set_attenuation_or_cross_grid_inference': True,
    }
    (out / 'private_figure_manifest_v1.json').write_text(json.dumps(manifest, indent=2) + '\n')
    return manifest


def main():
    task = Path(os.environ['PANCANCER_LN_STRUCTURE_TASK']).resolve()
    repository = Path(__file__).resolve().parents[1]
    assert task != repository and repository not in task.parents
    report, report_path = checked_source_report(task)
    version = os.environ.get('PANCANCER_LN_FIGURE_VERSION', 'v1')
    assert re.fullmatch(r'v[1-9][0-9]*', version)
    out = task / f'source_specificity_figures_{version}'
    manifest = render_report(report, out)
    manifest['source_evidence_summary_SHA256'] = digest(report_path)
    (out / 'private_figure_manifest_v1.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps({'status': manifest['status'], 'figure_files': len(manifest['figures']),
                      'family_contract_SHA256': manifest['family_contract_SHA256']}))


if __name__ == '__main__':
    main()
