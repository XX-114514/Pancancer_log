"""Build V7 within-major node frequency inputs for CoVarNet (primary strict branch and legacy all-display branch)."""
import json, lzma, pickle, sys, hashlib
from pathlib import Path
import numpy as np, pandas as pd

run = Path(sys.argv[1]); codes_dir = run / 'inputs/v7_codes'
meta = json.load(open(codes_dir / 'v7_codes_meta.json'))
blob = pickle.loads(lzma.decompress((codes_dir / 'v7_codes.pkl.xz').read_bytes()))
c = {k: np.frombuffer(b, dtype=dt) for k, (dt, b) in blob.items()}
n = meta['n']; assert all(len(v) == n for v in c.values())
cat = meta['categories']
def lab(k): return np.asarray(cat[k], dtype=object)[c[k]]

sample = lab('sampleID'); major = lab('major_celltype_v7'); node = lab('final_annotation_with_malignancy_v7')
branches = {
    'primary_strict': (c['strict_major_eligible_v7'] == 1) & (c['exclude_from_primary_major_v2'] != 1) & (c['strict_downstream_subtype_eligible_v7'] == 1),
    'legacy_all_display': np.ones(n, bool),
}
MIN_SAMPLE_CELLS, MIN_NODE_SAMPLES = 100, 10
for name, keep in branches.items():
    out = run / 'inputs' / name; out.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame({'sampleID': sample[keep], 'majorCluster': major[keep], 'subCluster': node[keep]})
    per_sample = df.groupby('sampleID').size()
    samples = per_sample.index[per_sample >= MIN_SAMPLE_CELLS]
    df = df[df.sampleID.isin(samples)]
    counts = df.groupby(['subCluster', 'majorCluster', 'sampleID']).size().unstack('sampleID', fill_value=0)
    assert not counts.index.get_level_values('subCluster').duplicated().any(), 'node spans two majors'
    major_tot = df.groupby(['majorCluster', 'sampleID']).size().unstack('sampleID', fill_value=0)
    denom = major_tot.reindex(counts.index.get_level_values('majorCluster')).to_numpy()
    freq = np.divide(counts.to_numpy(), denom, out=np.zeros(counts.shape), where=denom > 0)  # missing major -> 0 (reference algorithm)
    freq = pd.DataFrame(freq, index=counts.index, columns=counts.columns)
    present = (counts > 0).sum(axis=1)
    kept = (present >= MIN_NODE_SAMPLES) & (freq.max(axis=1) > freq.min(axis=1))
    freq = freq[kept.to_numpy()].reset_index(level='majorCluster')
    freq.index.name = None
    freq.to_csv(out / 'frequency.tsv', sep='\t')
    counts.to_csv(out / 'node_sample_counts.tsv.gz', sep='\t')
    sm = pd.DataFrame({'sampleID': freq.columns[1:]})
    sm['source'] = sm.sampleID.str.split(r'::|\||__', regex=True).str[0]
    sm['n_cells'] = per_sample.reindex(sm.sampleID).to_numpy()
    sm.to_csv(out / 'metadata.tsv', sep='\t', index=False)
    summ = dict(branch=name, cells_in_filter=int(keep.sum()), cells_in_kept_samples=int(len(df)), samples=int(len(samples)),
                nodes_observed=int(len(counts)), nodes_kept=int(kept.sum()), sources=int(sm.source.nunique()),
                dropped_nodes=sorted(map(str, counts.index.get_level_values('subCluster')[~kept.to_numpy()])),
                min_sample_cells=MIN_SAMPLE_CELLS, min_node_samples=MIN_NODE_SAMPLES,
                v7_checkpoint_sha256=meta['source_sha256'], v7_index_sha256=meta['index_sha256_sorted_by_global_position'],
                frequency_sha256=hashlib.sha256((out / 'frequency.tsv').read_bytes()).hexdigest())
    json.dump(summ, open(out / 'input_summary.json', 'w'), indent=1)
    print(name, {k: v for k, v in summ.items() if k != 'dropped_nodes'})
