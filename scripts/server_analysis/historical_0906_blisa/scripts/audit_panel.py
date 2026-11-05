import os
from pathlib import Path
import csv, json
import h5py
root=Path(__file__).resolve().parents[1]
axes=list(csv.DictReader((root/'tables/focus_axes_definition.tsv').open(),delimiter='\t'))
paths=[(Path(os.environ['XENIUM_PILOT_ROOT']) / 'processed/PDAC_3_counts.h5ad')]
paths+=sorted((Path(os.environ['SPATIAL_ROOT']) / 'processed/h5ad').glob('*.counts.h5ad'))
rows=[]
for p in paths:
    with h5py.File(p) as f:
        key=f['var'].attrs.get('_index','_index')
        genes=set(f['var'][key].asstr()[:])
        obskeys=list(f['obs'])
        for a in axes:
            needed=set(a['ligand'].split('_')+a['receptor'].split('_'))
            missing=sorted(needed-genes)
            rows.append(dict(file=str(p),axis=a['axis'],cm=a['cm'],n_genes=len(genes),missing_genes=','.join(missing),gene_eligible=not missing,obs_columns=','.join(obskeys)))
with (root/'tables/panel_axis_coverage.tsv').open('w') as out:
    w=csv.DictWriter(out,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
print('Files checked:',len(paths))
for p in paths:
    a=[r for r in rows if r['file']==str(p) and r['gene_eligible']]
    if len(a)>=8 or p==paths[0]:print(p.name,len(a),'of 12 axes measurable at gene-list level')
