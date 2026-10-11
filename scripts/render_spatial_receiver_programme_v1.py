from pathlib import Path
import os, json, hashlib
from datetime import datetime, timezone
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
out = root / 'immune_same_name_family_source_review_v1/spatial_receiver_programme_exposure_v1'
d = pd.read_csv(out/'all70_primary_programme_exposure.tsv',sep='\t')
mapping = pd.read_csv(out/'all75_record_mapping.tsv',sep='\t')
assert len(d) == 70
f = out/'figures_v1';f.mkdir(exist_ok=False)
rows = d[['programme_interface_id','ligand','receptor','sender_family','receiver_family']].drop_duplicates().sort_values('programme_interface_id')
files = sorted(d.sample_id.unique()); labs = ['GSE287472\nauthor unavailable','LUAD\nTMA1','LUAD\nTMA2','Mesothelioma\nTMA1','Mesothelioma\nTMA2']
matrix = np.full((14,5),np.nan)
for i,row in enumerate(rows.itertuples()):
    for j,sid in enumerate(files):
        r = d[(d.programme_interface_id == row.programme_interface_id)&(d.sample_id == sid)].iloc[0]
        if r.status == 'ESTIMABLE_CONDITIONAL':matrix[i,j] = r.beta_per_unit_L_fraction
fig,ax = plt.subplots(figsize=(16,12)); cmap = plt.get_cmap('RdBu_r').copy();cmap.set_bad('#e2e2e2')
im = ax.imshow(np.ma.masked_invalid(matrix),cmap=cmap,norm=TwoSlopeNorm(vmin=-1,vcenter=0,vmax=1),aspect='auto')
labels = []
for row in rows.itertuples():
    programmes = sorted(mapping.loc[mapping.programme_interface_id == row.programme_interface_id,'program'].unique())
    assert len(programmes) == 1
    labels.append(row.programme_interface_id+' '+row.ligand+' → '+row.receptor+' | '+programmes[0]+'\n'+row.sender_family.replace('_',' ')+' → '+row.receiver_family.replace('_',' '))
ax.set_yticks(range(14),labels,fontsize=9);ax.set_xticks(range(5),labs,fontsize=10)
codes = {'AUTHOR_NOT_RETRIEVED':'A','AUTHOR_B_CATEGORY_UNAVAILABLE':'B','INSUFFICIENT_ELIGIBLE_RECEIVERS':'N','INSUFFICIENT_PROGRAMME_PANEL':'P','INSUFFICIENT_RESIDUAL_DF':'D'}
for i,row in enumerate(rows.itertuples()):
    for j,sid in enumerate(files):
        r = d[(d.programme_interface_id == row.programme_interface_id)&(d.sample_id == sid)].iloc[0]
        if r.status == 'ESTIMABLE_CONDITIONAL':txt = f'β={r.beta_per_unit_L_fraction:+.3f}\nn={int(r.eligible_receivers)}'
        else:txt = 'NA · '+('L' if r.status.startswith('MISSING_LIGAND') else codes[r.status])
        ax.text(j,i,txt,ha='center',va='center',fontsize=9,color='white' if np.isfinite(matrix[i,j]) and abs(matrix[i,j])>.68 else 'black')
ax.set_title('Receiver programme vs ligand-positive sender-neighbour fraction\n22 estimable / 70 fixed models; 0 nominal p<0.05; 0 primary BH q<0.05',fontsize=14)
cb = fig.colorbar(im,ax=ax,pad=.02);cb.set_label('Adjusted β per 0→1 ligand-positive fraction (colour clipped at ±1)')
fig.text(.5,.043,'NA reasons: A author labels unavailable; B author B category unavailable; L ligand gene absent;\nN <40 eligible receivers; P programme panel coverage insufficient; D residual degrees of freedom insufficient.',ha='center',fontsize=10)
fig.text(.5,.012,'Within R+ identities; full native RNA, author/FOV, family membership, area/density and sender-family degree adjusted.\n10,000 ligand-label permutations; full 70-test BH. FOV/TMA are not independent patients; generic programmes are not ligand-specific target evidence.',ha='center',fontsize=9)
fig.subplots_adjust(left=.40,bottom=.145,right=.88,top=.93)
for ext in ['png','pdf']:fig.savefig(f/('all70_receiver_programme_exposure.'+ext),dpi=160)
plt.close(fig)
d.drop(columns=['covariate_names']).to_csv(f/'figure_source.tsv',sep='\t',index=False,na_rep='NA')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
(f/'render_manifest.json').write_text(json.dumps(dict(rendered_UTC=datetime.now(timezone.utc).isoformat(),input_sha256=sha(out/'all70_primary_programme_exposure.tsv'),renderer_sha256=sha(Path(__file__)),source_rows=70,covariate_level_names_omitted=True,outputs={q.name:sha(q) for q in f.iterdir() if q.is_file()},visual_inspection='PENDING',scientific_release=False),indent=2)+'\n')
print(json.dumps(dict(status='RENDERED_PENDING_PIXEL_REVIEW',source_rows=70)))
