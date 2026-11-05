"""Tumour-tissue-only V7 primary inputs: keep samples whose 20261010 sample_class is primary/metastatic tumour; re-apply node filters."""
import sys, json
from pathlib import Path
import pandas as pd
run, sample_meta = Path(sys.argv[1]), sys.argv[2]
src = run / 'inputs/primary_strict'; out = run / 'inputs/primary_strict_tumour_only'; out.mkdir(exist_ok=True)
fq = pd.read_csv(src / 'frequency.tsv', sep='\t', index_col=0); meta = pd.read_csv(src / 'metadata.tsv', sep='\t')
counts = pd.read_csv(src / 'node_sample_counts.tsv.gz', sep='\t', index_col=[0, 1]); counts.index = counts.index.get_level_values(0)
s = pd.read_csv(sample_meta, sep='\t'); s['sid'] = s.sampleID.str.split('::', n=1).str[1]
cls = s.set_index('sid').sample_class.reindex(meta.sampleID)
keep = meta.sampleID[cls.isin(['primary_or_unspecified_tumor', 'metastatic_tumor']).to_numpy()].tolist()
f = fq[['majorCluster'] + keep]; num = f.drop(columns='majorCluster')
ok = ((counts.loc[f.index, keep] > 0).sum(1) >= 10) & (num.max(1) > num.min(1))
f[ok.values].to_csv(out / 'frequency.tsv', sep='\t'); meta[meta.sampleID.isin(keep)].to_csv(out / 'metadata.tsv', sep='\t', index=False)
summ = dict(samples=len(keep), sources=int(meta[meta.sampleID.isin(keep)].source.nunique()), nodes_kept=int(ok.sum()),
            nodes_dropped=list(f.index[~ok.values]), class_counts=cls.value_counts(dropna=False).to_dict(), sample_class_source='20261010 expression downstream sample_metadata.tsv')
json.dump(summ, open(out / 'input_summary.json', 'w'), indent=1); print({k: v for k, v in summ.items()})
