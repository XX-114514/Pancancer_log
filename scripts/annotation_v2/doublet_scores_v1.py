# Annotation v2 step A0 (2026-11-06): one uniform doublet score per cell for every unit (sample), used only as
# cluster-level supporting evidence (methods/tme-annotation-v2-lodi-workflow-20261105.md). Scrublet on raw counts of
# ALL cells of the unit (malignant included: heterotypic doublets form with any cell in the droplet run).
# Expected rate scales with the unit's recovered cells: 0.008 per 1,000 cells (10x user guide, the figure Lodi cite),
# clipped to [0.01, 0.15]. Many sources were author-filtered, so the absolute rate is an upper bound; scores are
# compared between clusters of the same object, never thresholded per cell for removal.
# Usage: doublet_scores_v1.py SHARD.h5ad OUT.tsv.gz SEED
import sys, json
import numpy as np, pandas as pd, anndata as ad, scrublet as scr

shard, out, seed = sys.argv[1], sys.argv[2], int(sys.argv[3])
a = ad.read_h5ad(shard)
C = (a.layers['counts'] if 'counts' in a.layers else a.X).tocsr()
rows = []
for u, idx in a.obs.groupby('unit_id', observed=True).indices.items():
    idx = np.asarray(idx); n = len(idx)
    rate = float(np.clip(0.008 * n / 1000, 0.01, 0.15))
    rec = {'cell_id': a.obs_names[idx], 'unit_id': u, 'expected_rate': rate, 'doublet_score': np.nan, 'scrublet_predicted': np.nan, 'scrublet_threshold': np.nan}
    if n >= 100:
        try:
            s = scr.Scrublet(C[idx], expected_doublet_rate=rate, random_state=seed)
            sc_, pred = s.scrub_doublets(n_prin_comps=min(30, n - 1), verbose=False)
            rec.update(doublet_score=sc_, scrublet_predicted=pred if pred is not None else np.nan,
                       scrublet_threshold=getattr(s, 'threshold_', np.nan))
        except Exception as e:
            print(json.dumps({'unit': u, 'error': str(e)[:200]}), flush=True)
    rows.append(pd.DataFrame(rec))
D = pd.concat(rows)
D.to_csv(out, sep='\t', index=False, compression='gzip')
print(json.dumps({'cells': len(D), 'units': int(D.unit_id.nunique()), 'scored': int(D.doublet_score.notna().sum()),
                  'median_score': float(D.doublet_score.median()), 'predicted_frac': float(pd.to_numeric(D.scrublet_predicted, errors='coerce').mean())}))
