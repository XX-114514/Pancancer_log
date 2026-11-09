# Annotation v2 (2026-11-09, user decision): stage 2 of the ten-seed consensus for the four tissue pools, which the 2026-11-08
# continuation had no seed files for. Applies the same consensus() and eligibility rule as consensus_labels_v1.py to the R2 reviewed
# clusters (keep_tissue is the action that marks tissue identity here).
# Usage: consensus_tissue_pools_v1.py CELL_LABELS_R3 REVIEWED_R2 SEED_DIR OUT_DIR
import sys, json, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / 'annotation_v2_resume_20261108'))
import pandas as pd
from consensus_labels_v1 import consensus

lab, rev, seeds, out = sys.argv[1], sys.argv[2], pathlib.Path(sys.argv[3]), pathlib.Path(sys.argv[4])
out.mkdir(exist_ok=False)
L = pd.read_csv(lab, sep='\t', dtype=str, keep_default_na=False).set_index('cell_id')
R2 = pd.read_csv(rev, sep='\t', dtype=str, keep_default_na=False)
sums = []
for pool in ['Kupffer', 'Microglia', 'Alveolar_macrophage', 'Sinusoidal_EC']:
    rv = R2[R2.pool == pool].set_index('cluster')
    s = pd.read_csv(seeds / f'{pool}.seeds.tsv.gz', sep='\t', dtype=str).set_index('cell_id'); assert s.index.is_unique
    original = s.run_r1.map(rv.r2_label); actions = s.run_r1.map(rv.action); assert original.notna().all()
    sc = [f'seed_{n}' for n in range(2026110800, 2026110810)]; assert not s[sc].isna().any().any()
    label, stability, ties, pairs, maps, members = consensus(original.values, s[sc].values)
    ll = L.loc[s.index]; low = (ll.r1_action == 'keep_flag_low_depth') | (actions == 'keep_flag_low_depth')
    retained = (ll.major_label == 'r3_' + pool) & actions.str.startswith('keep')
    allowed = set(rv.loc[rv.action.isin(['keep', 'keep_tissue']), 'r2_label']) - set(rv.loc[~rv.action.isin(['keep', 'keep_tissue']), 'r2_label'])
    safe = {n: all(x in allowed for x in v) for n, v in members.items()}
    label_s = pd.Series(label, index=s.index)
    elig = retained & ~low & actions.isin(['keep', 'keep_tissue']) & (stability >= .7) & ~ties & label_s.map(safe)
    d = pd.DataFrame({'cell_id': s.index, 'pool': pool, 'cluster': s.run_r1.values, 'review_label': original.values, 'review_action': actions.values,
                      'consensus_label': label, 'stability': stability, 'mode_tie': ties, 'retained': retained.values, 'low_depth': low.values, 'include_main_subtype': elig.values})
    d.to_csv(out / f'{pool}.cells.tsv.gz', sep='\t', index=False, compression='gzip')
    pd.DataFrame(pairs).to_csv(out / f'{pool}.pairs.tsv', sep='\t', index=False)
    (out / f'{pool}.components.json').write_text(json.dumps(members, indent=2))
    d.groupby(['consensus_label', 'review_action', 'include_main_subtype']).agg(cells=('cell_id', 'size'), stability_median=('stability', 'median')).round(3).reset_index().to_csv(out / f'{pool}.summary.tsv', sep='\t', index=False)
    r = d[d.retained]
    sums.append(dict(pool=pool, object_cells=len(d), retained=int(retained.sum()), retained_stable=int((r.stability >= .7).sum()), retained_stability_median=float(r.stability.median()),
                     keep_tissue_cells=int((d.review_action == 'keep_tissue').sum()), main_subtype_cells=int(elig.sum()), components=len(members), merged_pairs=sum(p['merge'] for p in pairs),
                     tied_majorities=sum(m['tied_majority'] for m in maps), modal_ties=int(ties.sum())))
    print(json.dumps(sums[-1]), flush=True)
pd.DataFrame(sums).to_csv(out / 'summary.tsv', sep='\t', index=False)
(out / 'DONE.json').write_text(json.dumps({'status': 'computed_not_final_annotation', 'n_seeds': 10, 'merge_threshold': .25, 'stable_threshold': .7}))
