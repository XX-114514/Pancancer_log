from pathlib import Path
import os,json,hashlib
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);O=R/'immune_same_name_family_source_review_v1/spatial_contact_native_preflight_v3';F=O/'figures_v1';assert not F.exists();F.mkdir()
a=json.loads((O/'independent_audit.json').read_text());assert a['status'].startswith('PASS_')
d=pd.read_csv(O/'all_contact_record_FOV_counts.tsv',sep='\t');t=pd.read_csv(O/'all_contact_record_testability.tsv',sep='\t');x=d.merge(t[['model_key','sample_id','ligand','receptor','sender_family','receiver_family']],on=['model_key','sample_id'],validate='many_to_one');x['interface']=x.ligand+'–'+x.receptor;x['min_endpoint']=np.minimum(x.sender_positive,x.receiver_positive)
# Duplicate discovery/programme rows must agree; a contact interface is counted once.
for _,g in x.groupby(['sample_id','interface','arm','fov']):assert g.min_endpoint.nunique(dropna=False)==1
x=x.drop_duplicates(['sample_id','interface','arm','fov']);y=x.groupby(['sample_id','interface','arm'],dropna=False).min_endpoint.max().reset_index()
interfaces=sorted(x.interface.unique());samples=list(t.sample_id.drop_duplicates());short=['LUAD study2','LUAD TMA1','LUAD TMA2','MESO TMA1','MESO TMA2'];arms=['qc_detected','qc_background','qc_background_pure','qc_background_author'];titles=['RNA > 0 + QC','Above background + QC','Above background + pure RNA lineage','Above background + author lineage']
fig,axs=plt.subplots(1,4,figsize=(17,5.4),sharey=True,layout='constrained');source=[]
for ax,arm,title in zip(axs,arms,titles):
 m=np.full((len(interfaces),len(samples)),np.nan)
 for i,lr in enumerate(interfaces):
  for j,sid in enumerate(samples):
   v=y.loc[(y['interface']==lr)&(y.sample_id==sid)&(y.arm==arm),'min_endpoint'];assert len(v)==1;m[i,j]=v.iloc[0];source.append(dict(interface=lr,sample_id=sid,arm=arm,max_FOV_min_endpoint=m[i,j]))
 cmap=plt.get_cmap('viridis').copy();cmap.set_bad('#eeeeee');im=ax.imshow(np.ma.masked_invalid(np.log1p(m)),cmap=cmap,aspect='auto',vmin=0,vmax=np.log1p(max(20,float(x.min_endpoint.max()))))
 ax.set_title(title,fontsize=10);ax.set_xticks(range(len(samples)),short,rotation=60,ha='right',fontsize=9);ax.set_yticks(range(len(interfaces)),interfaces,fontsize=10)
 for i in range(len(interfaces)):
  for j in range(len(samples)):
   v=m[i,j];label='NA' if np.isnan(v) else str(int(v));ax.text(j,i,label,ha='center',va='center',fontsize=8,color='black' if np.isnan(v) else 'white')
fig.suptitle('Fixed 9 LR interfaces: maximum per-FOV minimum of sender/receiver positive cells\nSame-cell overlap is counted in endpoints, not a contact; no proximity test performed',fontsize=12)
fig.savefig(F/'all9_interfaces_native_spatial_coverage.png',dpi=160);fig.savefig(F/'all9_interfaces_native_spatial_coverage.pdf');plt.close(fig)
pd.DataFrame(source).to_csv(F/'figure_source.tsv',sep='\t',index=False,na_rep='NA');(F/'render_manifest.json').write_text(json.dumps(dict(renderer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),producer_audit_sha256=hashlib.sha256((O/'independent_audit.json').read_bytes()).hexdigest(),plotted_interfaces=9,assays=5,arms=4,plot_is_descriptive=True,colour_scale='log1p counts, same log1p scale across all panels; printed values are cell counts',outputs={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in F.iterdir() if p.is_file()}),indent=2)+'\n');print('rendered fixed9x5x4 coverage with printed count denominators')
