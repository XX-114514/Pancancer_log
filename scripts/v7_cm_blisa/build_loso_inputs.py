"""Leave-one-source-out frequency inputs for the V7 primary branch (nodes must keep >=10 observed samples and nonzero span)."""
import sys, json
from pathlib import Path
import pandas as pd
run = Path(sys.argv[1]); src = run / 'inputs/primary_strict'
fq = pd.read_csv(src / 'frequency.tsv', sep='\t', index_col=0); meta = pd.read_csv(src / 'metadata.tsv', sep='\t')
counts = pd.read_csv(src / 'node_sample_counts.tsv.gz', sep='\t', index_col=[0, 1])
counts.index = counts.index.get_level_values(0)
rows = []
for s in sorted(meta.source.unique()):
    keep = meta.sampleID[meta.source != s].tolist(); f = fq[['majorCluster'] + keep]
    num = f.drop(columns='majorCluster')
    ok = ((counts.loc[f.index, keep] > 0).sum(1) >= 10) & (num.max(1) > num.min(1))
    d = run / f'inputs/loso/{s}'; d.mkdir(parents=True, exist_ok=True)
    f[ok.values].to_csv(d / 'frequency.tsv', sep='\t'); meta[meta.source != s].to_csv(d / 'metadata.tsv', sep='\t', index=False)
    rows.append(dict(source=s, samples_left_out=int((meta.source == s).sum()), nodes_kept=int(ok.sum()), nodes_dropped=';'.join(f.index[~ok.values])))
pd.DataFrame(rows).to_csv(run / 'inputs/loso/loso_inputs.tsv', sep='\t', index=False); print(pd.DataFrame(rows).describe())
