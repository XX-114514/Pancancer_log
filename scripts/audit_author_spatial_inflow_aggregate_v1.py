"""Independent aggregate replay; no individual rows printed or new inference."""
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def adjusted(values, dependent=False):
    order = sorted(range(len(values)), key=lambda i: values[i])
    factor = sum(1 / (i + 1) for i in range(len(values))) if dependent else 1.
    result = [1.] * len(values)
    ceiling = 1.
    for rank in range(len(values), 0, -1):
        i = order[rank - 1]
        ceiling = min(ceiling, values[i] * len(values) * factor / rank)
        result[i] = ceiling
    return np.asarray(result)


def main():
    root = Path(os.environ['PANCANCER_SPATIAL_PRIVATE_OUTPUT'])
    terminal = root / 'private_terminal.json'
    t = json.loads(terminal.read_text())
    for name, entry in t['files'].items():
        p = root / name
        assert p.stat().st_size == entry['bytes'] and digest(p) == entry['sha256']
    results = pd.read_csv(root / 'private_complete7_conditional_source_inflow_results.tsv', sep='\t')
    fields = pd.read_csv(root / 'private_all_field_inflow_and_decompositions.tsv.gz', sep='\t')
    excluded = pd.read_csv(root / 'private_all_field_exclusions.tsv.gz', sep='\t')
    assert len(results) == 7 and not results.duplicated(['dataset_id', 'cancer']).any()
    assert not fields.duplicated(['sample_id', 'FOV', 'readout', 'scale', 'null']).any()
    assert not excluded.duplicated(['sample_id', 'FOV']).any()
    assert set(zip(excluded.sample_id, excluded.FOV)).isdisjoint(set(zip(fields.sample_id, fields.FOV)))
    assert (excluded[['T_cells', 'MPS_cells']].min(axis=1) < 20).all()
    terms = ['structural_opportunity', 'sender_mark', 'receiver_mark', 'joint_centered_marks']
    np.testing.assert_allclose(fields[terms].sum(axis=1), fields.inflow, atol=1e-11, rtol=1e-11)
    assert fields[['inflow', 'analytic_null_expected']].ge(-1e-12).all().all()
    primary = fields[(fields.readout == 'BACKGROUND_EXCESS') & (fields.scale == 'PRIMARY6NN_CAP2')]
    streams = np.load(root / 'private_complete_source_aggregate_null_streams.npz', allow_pickle=False)
    keys = set()
    for index, row in results.iterrows():
        part = primary[(primary.dataset_id == row.dataset_id) & (primary.cancer == row.cancer)]
        if row.measured_images == 0:
            assert part.empty and pd.isna(row.observed) and row.joint_P == 1
            continue
        assert part.sample_id.nunique() == row.measured_images
        assert len(part[['sample_id', 'FOV']].drop_duplicates()) == row.eligible_technical_FOVs
        probabilities = []
        for mode in ['DEPTH_AREA', 'DEPTH_AREA_XY']:
            section = part[part['null'] == mode]
            observed = section.groupby('sample_id').inflow.mean().mean()
            expectation = section.groupby('sample_id').analytic_null_expected.mean().mean()
            np.testing.assert_allclose([observed, expectation], [row.observed, row[mode + '_expected']], rtol=1e-12, atol=1e-12)
            key = 'stratum%02d_%s' % (index, mode)
            keys.add(key)
            values = streams[key]
            assert len(values) == 1999 and np.isfinite(values).all()
            # Python scalar enumeration independent of producer vectorized count.
            upper = (1 + sum(float(v) >= observed - 1e-12 for v in values)) / 2000
            lower = (1 + sum(float(v) <= observed + 1e-12 for v in values)) / 2000
            probability = min(1., 2 * min(upper, lower)) if np.std(values) > 1e-14 else 1.
            np.testing.assert_allclose([probability, np.std(values, ddof=1)], [row[mode + '_P'], row[mode + '_null_SD']], atol=1e-12)
            probabilities.append(probability)
        assert row.joint_P == max(probabilities)
    assert set(streams.files) == keys
    for column, dependent in [('BH7', False), ('BY7', True)]:
        np.testing.assert_allclose(results[column], adjusted(results.joint_P.tolist(), dependent), atol=1e-12)
    # Conditional diagnostic descriptions keep every readout/scale/null, no P.
    summaries = []
    for labels, section in fields.groupby(['dataset_id', 'cancer', 'readout', 'scale', 'null']):
        values = section.groupby('sample_id')[['inflow', 'analytic_null_expected'] + terms].mean().mean()
        summaries.append(dict(zip(['dataset_id', 'cancer', 'readout', 'scale', 'null'], labels), **values.to_dict()))
    pd.DataFrame(summaries).to_csv(root / 'private_complete_source_point_sensitivities.tsv', sep='\t', index=False)
    report = dict(status='ACCEPTED_ALL7_SOURCE_P_BH_BY_HIERARCHICAL_AGGREGATION_AND_ALL_FIELD_DECOMPOSITIONS',
                  reviewed_UTC=datetime.now(timezone.utc).isoformat(), producer_terminal_sha256=digest(terminal),
                  source_positions=7, aggregate_null_values_replayed=sum(len(streams[k]) for k in keys),
                  field_decomposition_positions=len(fields), new_P_family=False,
                  raw_or_graph_reconstruction_scope='Producer full-native parity and explicit-edge checks reviewed; this audit replays stored aggregates, not a second raw/graph rebuild.',
                  original_identity_or_inference_changed=False, scientific_outputs_public=False)
    path = root / 'private_independent_aggregate_acceptance.json'
    assert not path.exists()
    path.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(dict(status=report['status'], audit_sha256=digest(path), aggregate_null_values_replayed=report['aggregate_null_values_replayed'], field_decomposition_positions=len(fields))))


if __name__ == '__main__':
    main()
