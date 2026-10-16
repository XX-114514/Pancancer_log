"""Freeze all-source identity-internal, non-discovery-cell spatial proxy calibration."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,shutil,importlib.metadata
import numpy as np,pandas as pd,anndata as ad
from scipy.sparse import csr_matrix
import pyucell as uc

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);planp=Path('/tmp/state_spatial_proxy_calibration_plan_v1.json');plan=json.loads(planp.read_text());assert sha(__file__)==plan['producer_sha256'];assert importlib.metadata.version('pyucell')=='0.8.0';assert shutil.disk_usage(V).free>500*2**30
F=R/'all_CM_state_native_input_preflight_20261014_v1';iv=json.loads((F/'inventory.json').read_text());assert sha(F/'inventory.json')==plan['scope_inventory_sha256'];up=F/'private_all_frozen_sample_node_units.tsv.gz';rp=F/'private_canonical_row_selection.npz';assert sha(up)==iv['private_unit_table_sha256'] and sha(rp)==iv['private_row_selection_sha256'];units=pd.read_csv(up,sep='\t',keep_default_na=False);z=np.load(rp,allow_pickle=False);rows=z['global_rows'];uid=z['unit_index'];assert len(rows)==4532933 and np.all(np.diff(rows)>0)
D=V/'all_lineage_state_candidate_dictionary_20261014_v4/frozen_programme_dictionary.json';assert sha(D)==plan['dictionary_sha256'];records=json.loads(D.read_text())['programmes'];assert len(records)==814
P=V/'all_source_native_UCell_projection_20261014_v2';assert sha(P/'inventory.json')==plan['projection_inventory_sha256'];pi=json.loads((P/'inventory.json').read_text());assert pi['cells']==4532933
M=R/'all_lineage_state_discovery_row_selection_20261014_v1/private_discovery_metadata.tsv.gz';training=pd.read_csv(M,sep='\t',usecols=['core_position','gse_id','major']);assert len(training)==100766;excluded=np.isin(rows,training.core_position.to_numpy());pool=np.flatnonzero(~excluded);order=pool[np.argsort(uid[pool],kind='stable')];offset=np.r_[0,np.cumsum(np.bincount(uid[pool],minlength=len(units)))];poolcounts=np.diff(offset);assert len(order)==4432167;source_train=set(zip(training.gse_id,training.major));rng=np.random.default_rng(plan['seed']);selected=[];coverage=[]
for (sid,major,node,stratum),g in units.groupby(['gse_id','major','node_id','stratum'],sort=True):
 indices=g.unit_index.to_numpy(int);eligible=indices[poolcounts[indices]>0];bg=units.iloc[eligible].groupby('effective_resampling_block',sort=True);blocks=list(bg.groups);choose=sorted(rng.choice(blocks,min(plan['max_blocks_per_source_node_stratum'],len(blocks)),replace=False).tolist()) if blocks else [];count=0;reviewed=0
 for block in choose:
  ii=bg.get_group(block).unit_index.to_numpy(int);cs=poolcounts[ii];cum=np.r_[0,np.cumsum(cs)];take=rng.choice(int(cum[-1]),min(plan['max_cells_per_block'],int(cum[-1])),replace=False)
  for k in take:
   j=int(np.searchsorted(cum,k,side='right')-1);selected.append(int(order[offset[ii[j]]+k-cum[j]]));count+=1
  reviewed+=int(units.iloc[ii].unit_kind.eq('reviewed_patient').all())
 coverage.append(dict(gse_id=sid,major=major,node_id=node,stratum=stratum,original_cells=int(g.n_cells.sum()),non_discovery_cell_pool=int(poolcounts[indices].sum()),available_analysis_blocks=len(blocks),selected_analysis_blocks=len(choose),selected_reviewed_patient_blocks=reviewed,selected_cells=count,major_source_used_for_discovery=(sid,major) in source_train,status='FIXED_NON_DISCOVERY_CELLS' if count else 'NO_NON_DISCOVERY_CELLS_NOT_BIOLOGICAL_NEGATIVE'))
selected=np.array(sorted(selected),int);assert len(np.unique(selected))==len(selected) and not excluded[selected].any();assert len(selected)<=plan['maximum_selected_cells'];O=V/'state_spatial_proxy_calibration_inputs_20261016_v1';assert not O.exists();O.mkdir(mode=0o700);np.savez_compressed(O/'private_fixed_calibration_rows.npz',global_rows=rows[selected],unit_index=uid[selected],full_selection_index=selected);pd.DataFrame(coverage).to_csv(O/'all_source_node_stratum_calibration_coverage.tsv',sep='\t',index=False)
S=V/'frozen_state_spatial_panel_preflight_20261016_v2';assert sha(S/'inventory.json')==plan['spatial_preflight_inventory_sha256'];si=json.loads((S/'inventory.json').read_text());axespath=S/'all14_spatial_axes.tsv';assert sha(axespath)==si['output_sha256'][axespath.name];axes=pd.read_csv(axespath,sep='\t');manifest=R/'immune_same_name_family_source_review_v1/spatial_native_input_repair_v2/analysis_input_manifest.tsv';assert sha(manifest)==si['spatial_manifest_sha256'];mf=pd.read_csv(manifest,sep='\t');panels=[];panelcounts=[]
for r in axes.itertuples():
 row=mf[mf.sample_id.eq(r.sample_id)];assert len(row)==1;p=Path(row.iloc[0].analysis_path);assert sha(p)==r.input_sha256;a=ad.read_h5ad(p,backed='r');assert a.n_vars==r.features;lookup={};preferred=[]
 for j,name in enumerate(a.var_names):
  names={str(name)};symbol=str(a.var.iloc[j]['gene_symbol']) if 'gene_symbol' in a.var else str(name);preferred.append(symbol if symbol not in ['','nan','None'] else str(name))
  names.add(symbol)
  for gene in names:
   if gene and gene not in ['nan','None']:lookup.setdefault(gene,set()).add(j)
 genes=sorted(g for g,js in lookup.items() if len(js)==1);axis_sha=hashlib.sha256(('\n'.join(genes)+'\n').encode()).hexdigest();panels.append(dict(dataset_id=r.dataset_id,sample_id=r.sample_id,features=int(r.features),unambiguous_symbols=genes,symbol_to_feature_index={g:next(iter(lookup[g])) for g in genes},preferred_symbol_by_feature=preferred,gene_axis_sha256=axis_sha,input_sha256=r.input_sha256));panelcounts.append(dict(dataset_id=r.dataset_id,sample_id=r.sample_id,features=int(r.features),unambiguous_symbols=len(genes),gene_axis_sha256=axis_sha));a.file.close()
(O/'fixed_all14_panel_gene_axes.json').write_text(json.dumps(panels,indent=2)+'\n');pd.DataFrame(panelcounts).to_csv(O/'all14_panel_axis_identity.tsv',sep='\t',index=False)
# Tiny official zero/tie/rank-cap checks before any biological-cell proxy scoring.
checks=[]
for cap in [5,8]:
 x=np.array([[0,0,0,0,0,0,0,0],[5,5,1,0,0,0,0,0],[1,2,3,4,5,6,7,8]],float);a=ad.AnnData(csr_matrix(x),var=pd.DataFrame(index=['A','B','C','D','E','F','G','H']));uc.compute_ucell_scores(a,signatures={'two':['A','B'],'missing':['A','Z']},max_rank=cap,ties_method='average',missing_genes='skip',n_jobs=1,device='cpu',suffix='_UCell');manual=[]
 for v in x:
  ranks=np.full(8,cap,float);idx=np.flatnonzero(v>0)
  for j in idx:ranks[j]=min(cap,1+np.count_nonzero(v[idx]>v[j])+(np.count_nonzero(v[idx]==v[j])-1)/2)
  manual.append([1-(ranks[:2].sum()-3)/(2*cap-3),1-(ranks[0]-1)/(cap-1)])
 actual=a.obs[['two_UCell','missing_UCell']].to_numpy();err=float(np.max(abs(actual-np.array(manual))));assert err<1e-6;assert (actual[0]==0).all();checks.append(dict(max_rank=cap,rows=3,signatures=2,max_manual_error=err,zero_is_zero=True,missing_skipped_not_filled=True))
(O/'official_UCell_small_panel_checks.json').write_text(json.dumps(checks,indent=2)+'\n');out=dict(status='FROZEN_ALL52_SOURCE_IDENTITY_INTERNAL_NONDISCOVERY_ROW_SELECTION_AND14_PANEL_AXES_NOT_CALIBRATED_YET',observed_UTC=datetime.now(timezone.utc).isoformat(),plan_sha256=sha(planp),producer_sha256=sha(__file__),dictionary_sha256=sha(D),projection_inventory_sha256=sha(P/'inventory.json'),spatial_preflight_inventory_sha256=sha(S/'inventory.json'),discovery_metadata_sha256=sha(M),frozen_cells=len(selected),sources=int(units.iloc[np.unique(uid[selected])].gse_id.nunique()),majors=int(units.iloc[np.unique(uid[selected])].major.nunique()),source_node_strata=len(coverage),strata_with_non_discovery_cells=sum(v['selected_cells']>0 for v in coverage),panels=14,unique_symbol_axes=len(set(v['gene_axis_sha256'] for v in panels)),original_full_signature_gate_unchanged=True,formal_state_or_spatial_biological_validation=False,files={p.name:sha(p) for p in O.iterdir() if p.is_file()});(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='files'}))
