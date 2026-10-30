"""Full-member native marker review panels, separate from DEG discovery."""
import argparse,datetime,hashlib,json,sys,time
from pathlib import Path
import numpy as np
import pandas as pd
def load(p):return json.loads(Path(p).read_text())
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
 return h.hexdigest()
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def main(R):
 cfg=load(R/'config/spec_v1.json');P=Path(cfg['marker_producer']);A=Path(cfg['marker_auditor']);accepted=A/'outputs/all_native_profiles_acceptance_v1.json'
 while not accepted.exists():
  if (A/'provenance/supervisor_terminal_v1.json').exists():assert load(A/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
  time.sleep(5)
 accept=load(accepted);verified={r['source']:r['source_profile_SHA256'] for r in accept['source_checks']};records=load(P/'outputs/native_profiles_index_v1.json');contract=load(P/'outputs/marker_input_contract_v1.json');ng=contract['n_clusters'];N=contract['n_fit'];J=Path(cfg['joined']);sys.path.insert(0,str(J/'scripts'))
 from native_panel_aggregate_v1 import LINEAGE_PANELS,STATE_PANELS,IDENTITY_CANDIDATE_PANELS
 from extract_full_source_local_marker_evidence_v1 import storage_gate
 storage_gate(R,cfg)
 panels=dict(LINEAGE_PANELS);states={'T_core':LINEAGE_PANELS['T'],'NK_related':LINEAGE_PANELS['NK_related'],'CD4_CD8_review':['CD4','CD8A','CD8B'],**STATE_PANELS,**IDENTITY_CANDIDATE_PANELS};genes=list(dict.fromkeys(g for panel in [panels,states] for gs in panel.values() for g in gs));by={g:i for i,g in enumerate(genes)};sums=np.zeros((ng,len(genes)));detected=np.zeros_like(sums);denom=np.zeros_like(sums);source_support=np.zeros_like(sums);n=np.zeros(ng,np.int64)
 for r in records:
  p=P/'outputs'/r['path'];assert sha(p)==r['SHA256']==verified[r['source']]
  with np.load(p,allow_pickle=False) as f:
   names=f['native_gene_names'].tolist();index={g:i for i,g in enumerate(names)};ns=f['n_cells'];n+=ns
   for g in genes:
    if g not in index:continue
    j=by[g];k=index[g];sums[:,j]+=f['logCP10K_sum'][:,k];detected[:,j]+=f['detected'][:,k];denom[:,j]+=ns;source_support[:,j]+=(ns>0)
 assert n.sum()==N and (detected<=denom).all();mean=np.full_like(sums,np.nan);fraction=np.full_like(sums,np.nan);np.divide(sums,denom,out=mean,where=denom>0);np.divide(detected,denom,out=fraction,where=denom>0);table=[]
 for k in range(ng):
  for j,g in enumerate(genes):table.append(dict(cluster=k,gene=g,mean_log1p_CP10K=mean[k,j],detection_fraction=fraction[k,j],n_measured_cells=int(denom[k,j]),n_measured_sources=int(source_support[k,j]),n_all_cluster_cells=int(n[k]),measurement_status='measured' if denom[k,j]>0 else 'not_measured'))
 pd.DataFrame(table).to_csv(R/'outputs/native_marker_panel_values_v1.tsv',sep='\t',index=False)
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 from matplotlib.cm import ScalarMappable
 from matplotlib.colors import Normalize
 outputs=[]
 for name,panel in [('lineage',panels),('T_NK_state_and_candidate',states)]:
  selected=list(dict.fromkeys(g for gs in panel.values() for g in gs));ix=[by[g] for g in selected];m=mean[:,ix];d=fraction[:,ix];scaled=np.zeros_like(m)
  for j in range(len(selected)):
   valid=np.isfinite(m[:,j]);lo=m[valid,j].min() if valid.any() else 0;hi=m[valid,j].max() if valid.any() else 0
   scaled[valid,j]=(m[valid,j]-lo)/(hi-lo) if hi>lo else 0;scaled[~valid,j]=np.nan
  fig,ax=plt.subplots(figsize=(max(12,len(selected)*.25+3),max(9,ng*.24+2.8)));xx,yy=np.meshgrid(np.arange(len(selected)),np.arange(ng));ok=np.isfinite(m);ax.scatter(xx[ok],yy[ok],s=110*d[ok],c=scaled[ok],cmap='viridis',vmin=0,vmax=1,edgecolors='none');ax.scatter(xx[~ok],yy[~ok],s=12,c='#999999',marker='x');ax.set_xticks(np.arange(len(selected)),selected,rotation=90,fontsize=8);ax.set_yticks(np.arange(ng),[str(k)+'  (n='+format(int(n[k]),',')+')' for k in range(ng)],fontsize=7);ax.invert_yaxis();ax.set_xlim(-.7,len(selected)-.3);ax.set_ylim(ng-.3,-.7);ax.set_ylabel('Numeric cluster | all fitted members');ax.grid(axis='y',color='#eeeeee',lw=.5);ax.set_axisbelow(True)
  handles=[ax.scatter([],[],s=110*f,color='#777777',label=str(int(f*100))+'%') for f in [.1,.25,.5,1]];ax.legend(handles=handles,title='Detection among measured cells',loc='upper center',bbox_to_anchor=(.5,-.18),ncol=4,frameon=False,fontsize=8);fig.colorbar(ScalarMappable(norm=Normalize(0,1),cmap='viridis'),ax=ax,pad=.015,fraction=.025,label='Gene-wise relative mean log1p(CP10K)');ax.set_title('Current TME | full native '+name.replace('_',' ')+' marker review\nPrior panels for biological review; not new DEG calls or final identities',fontsize=13);fig.tight_layout();p=R/'outputs'/('new_cluster_'+name+'_native_dotplot_v1.png');fig.savefig(p,dpi=170,bbox_inches='tight');plt.close(fig);outputs.append(dict(name=p.name,bytes=p.stat().st_size,SHA256=sha(p)))
 save(R/'outputs/native_marker_panel_review_manifest_v1.json',dict(n_fit=N,n_clusters=ng,n_marker_genes=len(genes),all_full_member_native_profiles_independently_accepted=True,prespecified_panels_not_new_DEG_or_independent_identity_validation=True,mean_expression_color_scaled_within_gene_for_display_only=True,point_size_detection_among_measured_cells=True,unmeasured_not_filled_as_zero=True,noncommon_marker_coverage_reported=True,source_profile_acceptance_SHA256=sha(accepted),figures=outputs))
 save(R/'provenance/terminal_v1.json',dict(status='FULL_MEMBER_NATIVE_LINEAGE_STATE_MARKER_REVIEW_PANELS_COMPLETE_NOT_FINAL_ANNOTATION',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),n_fit=N,files=outputs))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
