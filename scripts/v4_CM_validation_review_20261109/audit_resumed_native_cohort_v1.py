"""Independent stdlib replay of all66 patient effects/P values/FDR; no identity claim."""
from pathlib import Path
from itertools import product
from statistics import median
from datetime import datetime, timezone
import csv
import json
import math
import sys

root = Path(sys.argv[1]).resolve()
out = root / 'outputs_resume_v1'
assert (out / 'terminal_PRIVATE.json').exists()
config = json.loads((root / 'config/cohort_PRIVATE.json').read_text())
cases = config['cases']
assert len(cases) == len({v['donor'] for v in cases}) == 6
assert len({v['receipt']['members'][0]['sha256'] for v in cases}) == 6


def table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def number(value):
    return float(value) if value not in (None, '') else math.nan


def exact_p(values):
    d = [v for v in values if math.isfinite(v) and abs(v) > 1e-12]
    if not d:
        return 1.
    absolute = list(map(abs, d))
    ranks = [1 + sum(a < v for a in absolute) + .5 * (sum(a == v for a in absolute) - 1)
             for v in absolute]
    observed = sum(rank for rank, value in zip(ranks, d) if value > 0)
    scores = [sum(rank for rank, include in zip(ranks, flags) if include)
              for flags in product([False, True], repeat=len(d))]
    return sum(score >= observed - 1e-12 for score in scores) / len(scores)


def compare(observed, expected):
    observed = number(observed)
    if math.isnan(expected):
        assert math.isnan(observed)
        return 0.
    assert math.isfinite(observed) and abs(observed - expected) < 1e-10, (observed, expected)
    return abs(observed - expected)


measurements = {}
conditional = []
for case in cases:
    unit = case['unit']
    folder = root / 'outputs_direct_v1/spatial_controls_v1' / unit
    if not (folder / 'complete.json').exists():
        folder = out / 'spatial_controls_v1' / unit
    receipt = json.loads((folder / 'complete.json').read_text())
    assert receipt['native_source_file_sha256'] == case['receipt']['members'][0]['sha256']
    assert receipt['coordinate_contract']['status'].startswith('PASS_')
    n = table(folder / 'neighborhood.tsv')
    r = table(folder / 'regional_ES.tsv')
    e = table(folder / 'conditional_LR_edges.tsv')
    assert len(n) == 12 and len(r) == len(e) == 27
    for row in n:
        radius = int(row['radius_um'])
        receiver = row['receiver']
        family = 'neighborhood_primary2' if radius == 50 and receiver in ('B', 'NK') else 'neighborhood_sensitivity10'
        key = family, f'{receiver}_{radius}um'
        effect = number(row.get('log_odds_enrichment'))
        if math.isfinite(effect):
            a, b = float(row['target_receiver']), float(row['other_receiver'])
            n1, n0 = float(row['n_target_nonEC']), float(row['n_other_nonEC'])
            independently = math.log((a + .5) / (n1 - a + .5)) - math.log((b + .5) / (n0 - b + .5))
            assert abs(independently - effect) < 1e-10
            measurements.setdefault(key, []).append([effect])
        else:
            measurements.setdefault(key, [])
        if radius == 50:
            conditional.append({'receiver': receiver,
                                'delta_conditional_lymph_fraction': number(row.get('target_conditional_lymph_fraction')) - number(row.get('other_conditional_lymph_fraction'))})
    for name, rows, status, columns in [
        ('regional_LR27', r, 'eligible_descriptive_effect', ['AB_minus_A', 'AB_minus_B']),
        ('conditional_edge_LR27', e, 'eligible_conditional_spatial_null', ['excess_edge_fraction'])]:
        for row in rows:
            key = name, f"{row['sender']}>{row['receiver']}:{row['ligand']}_{row['receptor']}"
            values = measurements.setdefault(key, [])
            if row['status'] == status:
                value = [number(row.get(column)) for column in columns]
                assert all(math.isfinite(v) for v in value)
                values.append(value)

saved = table(out / 'cohort_summary_PRIVATE/all66_patient_hypotheses.tsv')
assert len(saved) == len(measurements) == 66
errors, expected_p = [], {}
for row in saved:
    key = row['family'], row['hypothesis']
    values = measurements[key]
    assert int(row['n_patients']) == len(values)
    if row['family'] == 'regional_LR27':
        p = [exact_p([v[i] for v in values]) if len(values) >= 5 else math.nan for i in [0, 1]]
        for column, value in zip(['P_vs_A', 'P_vs_B'], p):
            errors.append(compare(row[column], value))
        expected_p[key] = max(p) if all(math.isfinite(x) for x in p) else math.nan
    else:
        expected_p[key] = exact_p([v[0] for v in values]) if len(values) >= 5 else math.nan
    errors.append(compare(row['P'], expected_p[key]))

families = {'neighborhood_primary2': 2, 'neighborhood_sensitivity10': 10,
            'regional_LR27': 27, 'conditional_edge_LR27': 27}
summary = []
for family, size in families.items():
    rows = [v for v in saved if v['family'] == family]
    assert len(rows) == size
    order = sorted(range(size), key=lambda i: expected_p[family, rows[i]['hypothesis']]
                   if math.isfinite(expected_p[family, rows[i]['hypothesis']]) else 1.)
    bh, running = {}, 1.
    for rank in range(size, 0, -1):
        i = order[rank - 1]
        p = expected_p[family, rows[i]['hypothesis']]
        running = min(running, (p if math.isfinite(p) else 1.) * size / rank)
        bh[i] = running
    harmonic = sum(1 / n for n in range(1, size + 1))
    for i, row in enumerate(rows):
        valid = math.isfinite(expected_p[family, row['hypothesis']])
        errors.append(compare(row['BH'], bh[i] if valid else math.nan))
        errors.append(compare(row['BY'], min(1., bh[i] * harmonic) if valid else math.nan))
    summary.append({'family': family, 'planned': size,
                    'estimable': sum(math.isfinite(number(v['P'])) for v in rows),
                    'BH05': sum(number(v['BH']) < .05 for v in rows),
                    'BY05': sum(number(v['BY']) < .05 for v in rows)})
conditional_summary = []
for receiver in ['B', 'T', 'NK']:
    values = [v['delta_conditional_lymph_fraction'] for v in conditional
              if v['receiver'] == receiver and math.isfinite(v['delta_conditional_lymph_fraction'])]
    conditional_summary.append({'receiver': receiver, 'n_sources': len(values),
                                'median_fraction_difference': median(values),
                                'positive_differences': sum(v > 0 for v in values),
                                'descriptive_only_no_new_P': True})
proof = {'UTC': datetime.now(timezone.utc).isoformat(),
         'status': 'PASS_independent6_source66_hypothesis_scope_exact_P_and_complete_FDR',
         'sources': 6, 'unique_donor_keys': 6, 'original_coordinate_contracts_checked': 6,
         'hypotheses': 66, 'maximum_numeric_error': max(errors), 'families': summary,
         'conditional_lymphocyte_fraction_at50um_descriptive': conditional_summary,
         'identity_visual_causal_or_clinical_acceptance_claim': False}
with (out / 'independent_numeric_audit.json').open('x') as handle:
    json.dump(proof, handle, indent=2)
print(json.dumps(proof), flush=True)
