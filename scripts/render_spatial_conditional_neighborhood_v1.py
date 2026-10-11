from pathlib import Path
import os,json,hashlib
from datetime import datetime,timezone
import numpy as np,pandas as pd
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=root/'immune_same_name_family_source_review_v1/spatial_contact_conditional_neighborhood_v1';d=pd.read_csv(out/'all45_interface_results.tsv',sep='\t');assert len(d)==45
f=out/'figures_v3';f.mkdir(exist_ok=False)
axes=d[['interface_id','ligand','receptor','sender_family','receiver_family']].drop_duplicates().sort_values('interface_id');files=sorted(d.sample_id.unique());labs=['GSE287472\nauthor unavailable','LUAD\nTMA1','LUAD\nTMA2','Mesothelioma\nTMA1','Mesothelioma\nTMA2']
m=np.full((9,5),np.nan)
for i,ax in enumerate(axes.itertuples()):
 for j,sid in enumerate(files):
  r=d[(d.interface_id==ax.interface_id)&(d.sample_id==sid)].iloc[0]
  if np.isfinite(r.enrichment_ratio):m[i,j]=np.clip(np.log2(max(r.enrichment_ratio,.25)),-2,2)
fig,chart=plt.subplots(figsize=(13,8));cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#dddddd');im=chart.imshow(np.ma.masked_invalid(m),cmap=cmap,norm=TwoSlopeNorm(vmin=-2,vcenter=0,vmax=2),aspect='auto')
chart.set_xticks(range(5),labs);chart.set_yticks(range(9),[r.ligand+' – '+r.receptor+'\n'+r.sender_family.replace('_',' ')+' → '+r.receiver_family.replace('_',' ') for r in axes.itertuples()],fontsize=9)
for i,ax in enumerate(axes.itertuples()):
 for j,sid in enumerate(files):
  r=d[(d.interface_id==ax.interface_id)&(d.sample_id==sid)].iloc[0]
  if r.status!='ESTIMABLE_CONDITIONAL':text='NA\n'+('author unavailable' if r.status=='AUTHOR_NOT_RETRIEVED' else 'B category\nunavailable' if r.status=='AUTHOR_B_CATEGORY_UNAVAILABLE' else 'gene absent' if r.status.startswith('MISSING') else 'null invariant')
  else:text=f'{int(r.observed_pairs)} / {r.null_mean:.1f}\nq={r.BH_full45:.3f}'
  chart.text(j,i,text,ha='center',va='center',fontsize=8,color='white' if np.isfinite(m[i,j]) and abs(m[i,j])>1.35 else 'black')
chart.set_title('Fixed spatial RNA interfaces: observed / conditional-null mean\n29 estimable of 45; 0 primary BH q<0.05; NA is not zero',fontsize=12)
bar=fig.colorbar(im,ax=chart,pad=.02);bar.set_label('log2 observed/null ratio (colour clipped at ±2)')
fig.text(.5,.016,'Within-FOV area-normalized finite neighbourhood; no self-pairs. Conditional on family location, RNA depth and density.\nTMA/FOV ≠ independent patient; centroid proximity ≠ protein binding. All 45 primary tests retained.',ha='center',fontsize=9)
fig.subplots_adjust(left=.28,bottom=.13,right=.89,top=.9)
for ext in ['png','pdf']:fig.savefig(f/('all45_conditional_neighborhood.'+ext),dpi=160)
plt.close(fig);d.to_csv(f/'figure_source.tsv',sep='\t',index=False,na_rep='NA')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
meta=dict(rendered_UTC=datetime.now(timezone.utc).isoformat(),input_sha256=sha(out/'all45_interface_results.tsv'),renderer_sha256=sha(Path(__file__)),source_rows=45,outputs={q.name:sha(q) for q in f.iterdir() if q.is_file()},visual_inspection='PENDING',scientific_release=False)
(f/'render_manifest.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(dict(status='RENDERED_PENDING_PIXEL_REVIEW',source_rows=45)))
