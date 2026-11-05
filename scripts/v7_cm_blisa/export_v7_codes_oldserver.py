"""Run on the old server: export compact V7 per-cell label codes (global_position order) from the frozen V7 checkpoint.
Usage: python export_v7_codes_oldserver.py <final_annotation_v7.pkl.gz> <out_dir>"""
import sys, os, json, hashlib, lzma, pickle
import numpy as np, pandas as pd
p, out = sys.argv[1], sys.argv[2]; os.makedirs(out, exist_ok=True)
df = pd.read_pickle(p).sort_values('global_position')
assert (df.global_position.to_numpy() == np.arange(len(df))).all()
cat_cols = ['sampleID', 'major_celltype_v7', 'final_annotation_v7', 'final_annotation_with_malignancy_v7', 'malignancy_call_v7',
            'malignancy_tier_v7', 'annotation_confidence_v7', 'major_celltype_final']
bool_cols = ['strict_major_eligible_v7', 'exclude_from_primary_major_v2', 'strict_downstream_subtype_eligible_v7', 'strict_subtype_eligible_v7',
             'is_malignant_primary_v7', 'is_malignant_high_confidence_v7', 'strict_malignant_downstream_eligible_v7']
blob, cats = {}, {}
for c in cat_cols:
    s = df[c].astype('category'); v = s.cat.codes.to_numpy(); v = v.astype(np.int16 if v.max() >= 127 else np.int8)
    blob[c] = (str(v.dtype), v.tobytes()); cats[c] = [str(x) for x in s.cat.categories]
for c in bool_cols:  # nullable booleans: -1 = NA
    s = df[c]; v = np.where(s.isna(), -1, s.fillna(False).astype(bool).astype(np.int8)).astype(np.int8)
    blob[c] = (str(v.dtype), v.tobytes())
open(f'{out}/v7_codes.pkl.xz', 'wb').write(lzma.compress(pickle.dumps(blob, protocol=4), preset=9 | lzma.PRESET_EXTREME))
meta = dict(n=len(df), categories=cats, source_sha256=hashlib.sha256(open(p, 'rb').read()).hexdigest(),
            index_sha256_sorted_by_global_position=hashlib.sha256('\n'.join(df.index.astype(str)).encode()).hexdigest())
json.dump(meta, open(f'{out}/v7_codes_meta.json', 'w')); print(meta['n'], meta['source_sha256'])
