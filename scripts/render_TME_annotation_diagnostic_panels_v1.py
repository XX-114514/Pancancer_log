"""Static marker/decision-boundary diagnostics, without changing cell labels."""
import argparse,collections,datetime,hashlib,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from TME_all_lineage_cluster_policy_v1 import BROAD

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')

GENES=['CD3D','CD3E','TRAC','TRBC1','TRBC2','TRDC','TRGC1','TRGC2','NKG7','GNLY','CLEC4C','IL3RA','TCF4','LILRA4','SNAP25','SYT1','RBFOX3','STMN2','PF4','PPBP','GP9','HBB','HBA1','AHSP','GYPA','RGS5','CSPG4','PDGFRB','ACTA2','TAGLN','MYH11','COL1A1','DCN','LUM','SOX10','PLP1','OLIG1','OLIG2','MBP','MOG','GFAP','AQP4','MPZ','PRX','PMP22','PMEL','MLANA','MITF','TYR','DCT','TYRP1']

def main(root):
 cfg=load(root/'config/spec_v1.json');P=Path(cfg['producer']);A=Path(cfg['auditor']);D=Path(cfg['delivery'])
 for p in [P,A,D]:assert load(p/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 nodes=load(P/'outputs/cluster_hierarchy_named_catalog_v1.json')['nodes'];leaves=[r for r in nodes if r['final_leaf']];types=sorted({r['cluster_RNA_proposal']['celltype'] for r in leaves},key=lambda t:(t=='unresolved',t));source_support=collections.defaultdict(list);source_mean=collections.defaultdict(list);source_detect=collections.defaultdict(list)
 for entry in load(P/'outputs/cluster_hierarchy_index_v1.json'):
  mp=P/'outputs'/entry['metadata'];assert sha(mp)==entry['metadata_SHA256'];m=load(mp);assert sha(m['statistics']['path'])==m['statistics']['SHA256'];by={g:j for j,g in enumerate(m['native_measured_genes'])}
  with np.load(m['statistics']['path'],allow_pickle=False) as f:n=f['n'];log=f['log_sum'];det=f['detected']
  for typ in types:
   rr=[r for r in leaves if r['source']==entry['source'] and r['cluster_RNA_proposal']['celltype']==typ]
   if not rr:continue
   total=sum(r['n_records'] for r in rr);support=[]
   for family in BROAD:
    measurable=[r for r in rr if r['cluster_RNA_proposal']['native_gate_fractions'][family] is not None];den=sum(r['n_records'] for r in measurable);support.append(sum(r['cluster_RNA_proposal']['native_gate_fractions'][family]*r['n_records'] for r in measurable)/den if den else np.nan)
   source_support[typ].append(support);ids=[r['source_node_code'] for r in rr];assert int(n[ids].sum())==total;mean=np.full(len(GENES),np.nan);pct=mean.copy()
   for j,g in enumerate(GENES):
    if g in by:mean[j]=log[ids,by[g]].sum()/total;pct[j]=det[ids,by[g]].sum()/total
   source_mean[typ].append(mean);source_detect[typ].append(pct)
 plt.rcParams.update({'font.size':9,'pdf.fonttype':42});O=root/'outputs';files=[]
 with PdfPages(O/'annotation_boundary_diagnostics.pdf') as pdf:
  x=np.asarray([np.nanmean(source_support[t],axis=0) for t in types]);fig,ax=plt.subplots(figsize=(13,9));im=ax.imshow(np.ma.masked_invalid(x),aspect='auto',cmap='YlOrRd',vmin=0,vmax=.8);ax.set_xticks(range(len(BROAD)),BROAD,rotation=65,ha='right');ax.set_yticks(range(len(types)),types);ax.set_title('Same-cell coherent RNA modules used by the classifier\nEqual source weight; unmeasured modules masked; this is not independent accuracy');fig.colorbar(im,ax=ax,label='Fraction of cells with the complete module');fig.tight_layout();p=O/'classifier_core_support.png';fig.savefig(p,dpi=150,bbox_inches='tight');pdf.savefig(fig,bbox_inches='tight');plt.close(fig);files.append(p)
  mean=np.asarray([np.nanmean(source_mean[t],axis=0) for t in types]);pct=np.asarray([np.nanmean(source_detect[t],axis=0) for t in types]);z=(mean-np.nanmean(mean,axis=0))/np.maximum(np.nanstd(mean,axis=0),.05);xx,yy=np.meshgrid(range(len(GENES)),range(len(types)));ok=np.isfinite(z)&np.isfinite(pct);fig,ax=plt.subplots(figsize=(19,10));sc=ax.scatter(xx[ok],yy[ok],c=z[ok],s=120*pct[ok],cmap='RdBu_r',vmin=-2,vmax=2,edgecolors='#666666',linewidths=.2);ax.scatter(xx[~ok],yy[~ok],marker='x',c='#bbbbbb',s=10);ax.set_xticks(range(len(GENES)),GENES,rotation=90);ax.set_yticks(range(len(types)),types);ax.invert_yaxis();ax.set_title('Boundary markers: pDC / neuronal / platelet / erythroid / mural-fibroblast / melanocytic-glial\nFull native group profiles; color = gene-wise standardized mean, size = detection fraction');fig.colorbar(sc,ax=ax,shrink=.5);fig.tight_layout();p=O/'boundary_marker_dotplot.png';fig.savefig(p,dpi=150,bbox_inches='tight');pdf.savefig(fig,bbox_inches='tight');plt.close(fig);files.append(p)
  R=Path(cfg['overlap_producer']);RA=Path(cfg['overlap_auditor'])
  if (RA/'provenance/supervisor_terminal_v1.json').exists() and load(RA/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0:
   assert load(RA/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(R/'provenance/terminal_v1.json');rr=load(R/'outputs/unresolved_cluster_core_overlap_review_v1.json');joint=collections.defaultdict(float);den=collections.defaultdict(float)
   for r in rr:
    for p in r['pair_patterns']:key=(p['core_a'],p['core_b']);joint[key]+=p['n_both'];den[key]+=min(p['n_core_a'],p['n_core_b'])
   mat=np.full((len(BROAD),len(BROAD)),np.nan)
   for (a,b),v in den.items():
    if v:mat[BROAD.index(a),BROAD.index(b)]=mat[BROAD.index(b),BROAD.index(a)]=joint[a,b]/v
   fig,ax=plt.subplots(figsize=(11,10));im=ax.imshow(np.ma.masked_invalid(mat),vmin=0,vmax=1,cmap='magma');ax.set_xticks(range(len(BROAD)),BROAD,rotation=65,ha='right');ax.set_yticks(range(len(BROAD)),BROAD);ax.set_title('Core coexpression inside unresolved groups\nJoint positive / smaller marginal core count; masked = no supported pair');fig.colorbar(im,ax=ax,label='Same-cell core overlap');fig.tight_layout();p=O/'unresolved_core_overlap_heatmap.png';fig.savefig(p,dpi=150,bbox_inches='tight');pdf.savefig(fig,bbox_inches='tight');plt.close(fig);files.append(p)
 files.append(O/'annotation_boundary_diagnostics.pdf');save(O/'diagnostic_panel_manifest_v1.json',dict(files=[dict(name=p.name,bytes=p.stat().st_size,SHA256=sha(p)) for p in files],all_panels_are_diagnostics_not_relabeling=True,no_new_identity_or_CM=True));save(root/'provenance/terminal_v1.json',dict(status='ANNOTATION_BOUNDARY_MARKER_AND_CLASSIFIER_PANELS_RENDERED',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),n_files=len(files),manifest_SHA256=sha(O/'diagnostic_panel_manifest_v1.json'),all_labels_counts_and_states_unchanged=True))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
 try:main(root)
 except BaseException:
  import traceback
  save(root/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
