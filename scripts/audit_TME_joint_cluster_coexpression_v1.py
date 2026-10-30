"""Independent native-row Boolean replay of all new joint-cluster modules."""
import argparse,collections,datetime,hashlib,json,sys,time
from pathlib import Path
import numpy as np
import pandas as pd
from TME_all_lineage_cluster_policy_v1 import GATES,GATE_ORDER,FINE_BY_PARENT,annotate_cluster
from TME_assay_aware_T_anchor_policy_v1 import ACCESSORY,SIGNAL,CHAIN
from audit_TME_full_cluster_markers_v1 import read_rows
def load(p):return json.loads(Path(p).read_text())
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
 return h.hexdigest()
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def independent_flags(detected,genes):
 by={g:i for i,g in enumerate(genes)};n=len(detected);flags={};assessed={}
 def anyof(module):
  present=[by[g] for g in module if g in by]
  return np.logical_or.reduce(detected[:,present],axis=1) if present else np.zeros(n,bool)
 for name,modules in GATES.items():
  assessed[name]=all(any(g in by for g in module) for module in modules);flags[name]=np.logical_and.reduce([anyof(module) for module in modules])
 mode='original_chain_and_CD3_signaling_native_module'
 if not assessed['T']:
  if any(g in by for g in CHAIN):mode='incomplete_other_T_module_not_rescued'
  elif not (all(g in by for g in ['CD3D','CD3E']) and any(g in by for g in ACCESSORY) and sum(g in by for g in SIGNAL)>=2):mode='chain_unmeasured_and_alternative_anchor_unassessable'
  else:
   assessed['T']=True;flags['T']=anyof(['CD3D'])&anyof(['CD3E'])&anyof(ACCESSORY)&(np.sum([anyof([g]) for g in SIGNAL],axis=0)>=2);mode='T_chain_unmeasured_alternate_CD3_accessory_two_signal_RNA_evidence'
 for parent,children in FINE_BY_PARENT.items():
  for name in children:assessed[name]=assessed[name] and assessed[parent];flags[name]&=flags[parent]
 flags['CD4_T']&=~flags['Treg']
 return np.stack([flags[g] for g in GATE_ORDER],axis=1),assessed,mode
def main(R):
 cfg=load(R/'config/spec_v1.json');P=Path(cfg['producer']);terminal=P/'provenance/supervisor_terminal_v1.json'
 while not terminal.exists():time.sleep(5)
 assert load(terminal)['actual_child_exit_code']==0;pc=load(P/'config/spec_v1.json');H=Path(pc['producer']);G=Path(pc['graph']);hc=load(H/'config/spec_v1.json');J=Path(hc['joined']);sys.path.insert(0,str(J/'scripts'))
 from native_archive_collection_v1 import NativeArchiveCollection
 from anndata._io.specs import read_elem
 from extract_full_source_local_marker_evidence_v1 import storage_gate
 storage_gate(R,cfg);m=load(P/'outputs/joint_cluster_coexpression_manifest_v1.json');assert sha(P/'outputs/joint_cluster_coexpression_manifest_v1.json')==load(P/'provenance/terminal_v1.json')['manifest_SHA256'];assert m['gate_order']==GATE_ORDER
 for r in m['outputs']:assert sha(P/'outputs'/r['path'])==r['SHA256']
 for name,digest in m['policy_SHA256'].items():assert sha(R/'scripts'/name)==digest
 with np.load(H/'outputs/private_cell_axes_v1.npz',allow_pickle=False) as f:axes={k:f[k] for k in f.files}
 assert sha(H/'outputs/private_cell_axes_v1.npz')==m['input_axes_SHA256'] and sha(G/'outputs/private_full_graph_cluster_axes_v1.tsv.gz')==m['graph_axes_SHA256'];labels=pd.read_csv(G/'outputs/private_full_graph_cluster_axes_v1.tsv.gz',sep='\t',usecols=['leiden_0_5']).iloc[:,0].to_numpy(np.int32);ng=m['n_clusters'];N=m['n_fit'];assert len(labels)==N
 native=NativeArchiveCollection(load(J/'outputs/joined_native_atlas_delivery_v1.json')['archive_manifest']);lineage=load(H/'outputs/source_input_lineage_v1.json');source_meta={r['source']:r for r in load(P/'outputs/private_source_measurement_modes_v1.json')};proposals={(r['source'],r['cluster']):r for r in load(P/'outputs/private_source_cluster_coherent_proposals_v1.json')};bits=np.load(P/'outputs/private_all_cell_coherent_lineage_bitmask_v1.npy',mmap_mode='r');assert bits.shape==(N,) and bits.dtype==np.uint32 and not (bits>>len(GATE_ORDER)).any();totals=np.load(H/'outputs/private_analysis_RNA_totals_v1.npy',mmap_mode='r');counts=np.zeros((ng,len(GATE_ORDER)),np.int64);denom=np.zeros_like(counts);observed=np.bincount(labels,minlength=ng);sources_done=[]
 for code,r in enumerate(lineage):
  source=r['source'];wh=np.flatnonzero(axes['source_code']==code);z=native.store(source);var=read_elem(z['var']);lookup=native.feature_lookup(source);meta=source_meta[source];wanted={g for modules in GATES.values() for module in modules for g in module}|set(ACCESSORY+SIGNAL+CHAIN+['CD3D','CD3E']);genes=[];columns=[]
  for j,idx in enumerate(var.index):
   name=str(var.iloc[j]['gene_symbol']) if 'gene_symbol' in var else str(idx)
   if name in wanted and lookup.get(name,set())=={j}:genes.append(name);columns.append(j)
  assert genes==meta['measured_native_panel_gene_names'] and columns==meta['native_panel_physical_columns'];delta=z['layers/ambient_removed_counts'] if 'ambient_removed_counts' in z['layers'] else None;expected=np.zeros((ng,len(GATE_ORDER)),np.int64);n=np.bincount(labels[wh],minlength=ng);units=native.frozen_expression_units(source)[axes['source_native_position'][wh]].astype(str);missing=np.isin(np.char.lower(np.char.strip(units)),['','nan','none','na','n/a','unknown','__missing__','__source_value_missing__'])
  for lo in range(0,len(wh),1024):
   ix=wh[lo:lo+1024];raw=read_rows(z[r['counts_key']],axes['native_row'][ix],True);np.testing.assert_array_equal(np.asarray(raw.sum(axis=1)).ravel(),axes['native_total_RNA'][ix]);x=raw if delta is None else (raw-read_rows(delta,axes['native_row'][ix],False)).tocsr();assert np.isfinite(x.data).all() and x.data.min(initial=0)>=-1e-8;x.data[x.data<0]=0;x.eliminate_zeros();np.testing.assert_allclose(np.asarray(x.sum(axis=1)).ravel(),totals[ix],rtol=1e-10,atol=1e-7);f,assessed,mode=independent_flags(x[:,columns].toarray()>0,genes);packed=np.zeros(len(ix),np.uint32)
   for j,key in enumerate(GATE_ORDER):np.add.at(expected[:,j],labels[ix],f[:,j].astype(np.int64));packed+=f[:,j].astype(np.uint32)*(2**j)
   np.testing.assert_array_equal(bits[ix],packed)
  assert assessed==meta['assessed'] and mode==meta['T_evidence_mode'];assert meta['all_original_marker_alternatives_measured']=={g:all(all(x in genes for x in module) for module in modules) for g,modules in GATES.items()}
  for k in np.flatnonzero(n):
   row=proposals[source,str(k)];assert row['n_cells']==int(n[k]) and row['positive_gate_counts']==dict(zip(GATE_ORDER,expected[k].tolist()));assert row['source_cluster_RNA_proposal']==annotate_cluster(int(n[k]),dict(zip(GATE_ORDER,expected[k].tolist())),assessed);sel=labels[wh]==k;known=units[sel][~missing[sel]];unique,values=np.unique(known,return_counts=True);assert row['n_declared_expression_units_not_verified_patients_or_libraries']==len(unique) and row['n_missing_declared_expression_unit']==int(missing[sel].sum());actual=float(values.max()/values.sum()) if len(values) else None;assert row['largest_known_expression_unit_fraction']==actual
  counts+=expected;denom+=n[:,None]*np.asarray([assessed[g] for g in GATE_ORDER]);sources_done.append(source);native._fields.clear();save(R/'checkpoints'/(str(time.time_ns())+'.json'),dict(stage='independent_full_native_same_cell_modules_accepted',observed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),source=source,n_fit=len(wh)))
 df=pd.read_csv(P/'outputs/joint_cluster_coherent_module_coverage_v1.tsv',sep='\t',dtype={'cluster':str,'gate':str});assert len(df)==ng*len(GATE_ORDER)
 for k in range(ng):
  entries=[r for (source,cluster),r in proposals.items() if cluster==str(k)]
  for j,key in enumerate(GATE_ORDER):
   q=df[(df.cluster==str(k))&(df.gate==key)];assert len(q)==1;q=q.iloc[0];assert q.n_cluster_cells==observed[k] and q.n_module_assessable_cells==denom[k,j] and q.n_positive_coherent_cells==counts[k,j];np.testing.assert_allclose(q.positive_fraction_among_assessable,counts[k,j]/denom[k,j] if denom[k,j] else np.nan,equal_nan=True,atol=1e-12);eligible=[r for r in entries if r['assessed'][key]];adequate=[r for r in eligible if r['n_cells']>=30];assert q.n_assessable_sources==len(eligible) and q.n_adequate_sources==len(adequate)
   for subset,col in [(eligible,'source_equal_mean_fraction'),(adequate,'adequate_source_equal_mean_fraction')]:np.testing.assert_allclose(q[col],np.mean([r['positive_gate_counts'][key]/r['n_cells'] for r in subset]) if subset else np.nan,equal_nan=True,atol=1e-12)
 save(R/'outputs/accepted_joint_cluster_coexpression_v1.json',dict(n_fit=N,n_sources=len(sources_done),n_clusters=ng,all_native_rows_integer_counts_and_adopted_analysis_totals_replayed=True,all_cell_coexpression_bits_independently_reproduced=True,all_source_cluster_counts_measurability_modes_and_RNA_proposals_replayed=True,all_global_module_denominators_and_source_equal_fractions_replayed=True,anonymous_expression_unit_support_replayed_not_patient_validation=True,producer_terminal_SHA256=sha(P/'provenance/terminal_v1.json'),no_final_identity_or_malignancy_promoted=True));storage_gate(R,cfg)
 save(R/'provenance/terminal_v1.json',dict(status='ALL_NEW_JOINT_CLUSTER_NATIVE_SAME_CELL_MODULES_MEASUREMENT_AND_UNIT_SUPPORT_INDEPENDENTLY_REPLAYED_REQUIRES_BIOLOGICAL_REVIEW',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),n_fit=N,n_sources=len(sources_done)))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
