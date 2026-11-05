"""Sample a V7-labelled atlas reference (raw counts, panel-union genes) from the canonical zarr."""
import json, lzma, pickle, sys, glob
from pathlib import Path
import numpy as np, pandas as pd, zarr, anndata as ad, scipy.sparse as sp

run = Path(sys.argv[1]); Z = sys.argv[2]; xen = sys.argv[3]
PER_LABEL, PER_SAMPLE, HELDOUT_FRAC, SEED = 4000, 200, 0.2, 20261103
rng = np.random.default_rng(SEED)
meta = json.load(open(run / 'inputs/v7_codes/v7_codes_meta.json'))
blob = pickle.loads(lzma.decompress((run / 'inputs/v7_codes/v7_codes.pkl.xz').read_bytes()))
c = {k: np.frombuffer(b, dtype=dt) for k, (dt, b) in blob.items()}
lab = lambda k: np.asarray(meta['categories'][k], dtype=object)[c[k]]
strict = (c['strict_major_eligible_v7'] == 1) & (c['exclude_from_primary_major_v2'] != 1) & (c['strict_downstream_subtype_eligible_v7'] == 1)
major = lab('major_celltype_v7'); fine = lab('final_annotation_v7'); sample = lab('sampleID')
label = np.where(strict, fine, np.char.add('Other_', major.astype(str))).astype(object)
obs = pd.DataFrame({'label': label, 'sampleID': sample, 'major': major, 'fine_v7': fine, 'strict': strict,
                    'node_v7': lab('final_annotation_with_malignancy_v7'), 'malignancy_call_v7': lab('malignancy_call_v7')})
obs['source'] = obs.sampleID.str.split(r'::|\||__', regex=True).str[0]
# held-out samples are disjoint from training samples
samples = obs.sampleID.unique(); held = set(rng.choice(samples, int(len(samples) * HELDOUT_FRAC), replace=False))
obs['split'] = np.where(obs.sampleID.isin(held), 'heldout', 'train')
pick = []
for (lb, sp_), g in obs.groupby(['label', 'split'], sort=True):
    g = g.groupby('sampleID', group_keys=False).apply(lambda x: x.sample(min(len(x), PER_SAMPLE), random_state=SEED))
    cap = PER_LABEL if sp_ == 'train' else PER_LABEL // 4
    pick.append(g.sample(min(len(g), cap), random_state=SEED).index.to_numpy())
rows = np.sort(np.concatenate(pick))
var = zarr.open_group(Z + '/var', mode='r')['_index'][:].astype(str)
genes = set()
for f in glob.glob(xen + '/*.counts.h5ad'):
    genes |= set(ad.read_h5ad(f, backed='r').var_names)
gidx = np.array(sorted(i for i, g in enumerate(var) if g in genes)); remap = -np.ones(len(var), np.int64); remap[gidx] = np.arange(len(gidx))
X = zarr.open_group(Z + '/X', mode='r'); indptr = X['indptr'][:]
data_out, ind_out, ptr = [], [], [0]
for i, r in enumerate(rows):
    a, b = indptr[r], indptr[r + 1]
    j = remap[X['indices'][a:b]]; d = X['data'][a:b]; k = j >= 0
    ind_out.append(j[k]); data_out.append(d[k]); ptr.append(ptr[-1] + int(k.sum()))
    if i % 20000 == 0: print(i, len(rows), flush=True)
M = sp.csr_matrix((np.concatenate(data_out).astype(np.float32), np.concatenate(ind_out), np.array(ptr)), shape=(len(rows), len(gidx)))
o = obs.iloc[rows].copy(); o.index = [f'gp{r}' for r in rows]
ref = ad.AnnData(M, obs=o, var=pd.DataFrame(index=var[gidx]))
ref.uns['provenance'] = dict(v7_checkpoint_sha256=meta['source_sha256'], per_label=PER_LABEL, per_sample=PER_SAMPLE, heldout_frac=HELDOUT_FRAC, seed=SEED, zarr=Z)
(run / 'reference').mkdir(exist_ok=True); ref.write_h5ad(run / 'reference/v7_spatial_reference_counts.h5ad', compression='gzip')
print(ref); print(ref.obs.groupby(['label', 'split']).size().unstack().to_string())
