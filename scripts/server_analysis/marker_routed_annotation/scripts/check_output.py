"""Read back route/model isolation and source-ID completeness."""
import argparse
import json
from pathlib import Path
import sys
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from marker_policy import counterexamples, routes, EXPLICIT_NONIMMUNE

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('output', type=Path)
p.add_argument('--malignancy', type=Path, required=True)
a = p.parse_args()
summary = json.loads((a.output / 'summary.json').read_text())
d = pd.read_csv(a.output / 'candidate_labels.tsv.gz', sep='\t')
m = pd.read_csv(a.malignancy, sep='\t', usecols=['cell_id', 'malignancy_call']).set_index('cell_id')
assert d.cell_id.is_unique and len(d) == summary['output_cells']
if not summary['pilot']:
    assert set(d.cell_id) == set(m.index)
assert (d.malignancy_call.to_numpy() == m.reindex(d.cell_id).malignancy_call.to_numpy()).all()
assert (d.prediction_route == routes(d.marker_parent, d.malignancy_call)).all()
for route, models in [('immune_celltypist', ['Immune_All_High', 'Immune_All_Low']),
                      ('vascular_celltypist', ['Adult_Human_Vascular'])]:
    mask = d.prediction_route == route
    for model in models:
        if mask.any():
            assert d[model + '_candidate'].notna().equals(mask)
            for suffix in ['_max_probability', '_margin']:
                assert d.loc[mask, model + suffix].between(0, 1).all()
        else:
            assert model + '_candidate' not in d or d[model + '_candidate'].isna().all()
sr = d.prediction_route == 'singler_hpca'
if sr.any():
    assert d.hpca_main_raw.notna().equals(sr)
assert not d.reference_status.eq('pending').any()
if 'Immune_All_Low_candidate' in d:
    conflict = d.Immune_All_Low_candidate.isin(EXPLICIT_NONIMMUNE)
    assert (d.loc[conflict, 'reference_status'] == 'reference_parent_conflict_or_unresolved').all()
assert (d.loc[d.malignancy_call == 'malignant', 'prediction_route'] == 'excluded_malignant').all()
assert sum(summary['route_counts'].values()) == len(d)
print(json.dumps({'source': summary['source'], 'rows': len(d),
                  'route_isolation': True, 'source_ids': True,
                  'malignancy_preserved': True, 'policy_counterexamples': counterexamples()}))
