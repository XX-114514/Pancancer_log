"""Fixed diagnostic panels on all accepted discovery cells; not population prevalence."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,gc
import numpy as np,pandas as pd,anndata as ad
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);P=Path('/tmp/all18_discovery_panel_coexpression_plan_v1.json');plan=json.loads(P.read_text());assert sha(__file__)==plan['producer_sha256'];A=R/'all_lineage_state_native_input_audit_20261014_v2/accepted_input_inventory.json';assert sha(A)==plan['accepted_input_inventory_sha256'];iv=json.loads(A.read_text());panels=plan['panels'];O=V/'all18_discovery_panel_coexpression_20261014_v1';assert not O.exists();O.mkdir();rows=[];joint=[];coverage=[];total=0;maxdiff=0
for item in iv['files']:
 assert sha(item['path'])==item['sha256'];a=ad.read_h5ad(item['path']);major=item['major'];assert a.obs.major.eq(major).all();total+=a.n_obs;tot=a.obs.full_RNA_total.to_numpy();values={};ok={}
 for name,genes in panels.items():
  ix=a.var_names.get_indexer(genes);present=ix>=0;coverage.append(dict(major=major,panel=name,required_genes=len(genes),measured_genes=int(present.sum()),complete_gene_axis=bool(present.all()),missing_genes=';'.join(np.array(genes)[~present])));x=a.X[:,ix[present]];det=np.asarray((x>0).sum(axis=1)).ravel();fraction=np.asarray(x.sum(axis=1)).ravel()/tot;ok[name]=bool(present.all());values[name]=(det,fraction,(det>=3)&(fraction>=.01))
  # Independent gene-by-gene accumulation for every cell in the diagnostic panels.
  dd=np.zeros(a.n_obs,np.int64);ss=np.zeros(a.n_obs)
  for col in ix[present]:
   v=a.X[:,int(col)].toarray().ravel();dd+=(v>0);ss+=v
  assert np.array_equal(det,dd) and np.allclose(fraction,ss/tot,atol=1e-12);maxdiff=max(maxdiff,float(np.max(np.abs(fraction-ss/tot),initial=0)))
 for keys,positions in a.obs.groupby(['gse_id','node_id','subCluster'],observed=True,sort=True).indices.items():
  idx=np.asarray(positions);sid,node,label=keys
  for name,(det,frac,strong) in values.items():rows.append(dict(gse_id=sid,major=major,node_id=node,subCluster=label,cells=len(idx),panel=name,all_panel_genes_measured=ok[name],cells_detecting_atleast3_panel_genes=int((det[idx]>=3).sum()) if ok[name] else np.nan,cells_atleast3genes_and1pct_native_RNA=int(strong[idx].sum()) if ok[name] else np.nan,panel_native_RNA_fraction_median=float(np.median(frac[idx])) if ok[name] else np.nan,diagnostic_marker_support_not_state_positive_rate=True,balanced_discovery_not_cohort_prevalence=True))
  for left,right in plan['joint_panels']:
   complete=ok[left] and ok[right];n=int((values[left][2][idx]&values[right][2][idx]).sum()) if complete else np.nan;joint.append(dict(gse_id=sid,major=major,node_id=node,subCluster=label,cells=len(idx),panel_A=left,panel_B=right,complete_gene_axes=complete,joint3genes_and1pct_each_cells=n,joint_fraction=n/len(idx) if complete else np.nan,not_contamination_or_doublet_diagnosis=True))
 print(json.dumps(dict(major=major,cells=a.n_obs,phase='PANELS_CHECKED',UTC=datetime.now(timezone.utc).isoformat())),flush=True);del a;gc.collect()
assert total==100766
for name,tab in [('all_major_panel_gene_coverage.tsv',coverage),('all_source_node_panel_support.tsv',rows),('all_source_node_joint_panel_support.tsv',joint)]:pd.DataFrame(tab).to_csv(O/name,sep='\t',index=False)
out=dict(status='COMPLETE_ALL100766_DISCOVERY_PANEL_DIAGNOSTICS_NOT_STATE_LABELS',observed_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(P),accepted_input_inventory_sha256=sha(A),cells=total,majors=18,panels=len(panels),source_node_panel_rows=len(rows),joint_panel_rows=len(joint),all_cell_panel_quantities_independently_accumulated=True,max_fraction_difference=maxdiff,full_cohort_prevalence=False,formal_state=False,original_labels_changed=False,output_sha256={p.name:sha(p) for p in O.glob('*.tsv')});(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
