"""Native same-cell lineage evidence on the accepted joint TME cluster axis.

Reuse frozen all-lineage rules; do not assign global final identities, resolve
malignancy or interpret unmeasured modules as biological negative cells.
"""
import argparse,collections,datetime,hashlib,json,sys,time
from pathlib import Path
import numpy as np
import pandas as pd
from TME_all_lineage_cluster_policy_v1 import GATES,GATE_ORDER,coherent_flags,annotate_cluster
from TME_assay_aware_T_anchor_policy_v1 import ACCESSORY,SIGNAL,CHAIN,assay_aware_flags
from audit_TME_standard_Harmony_v2 import gather_nonnegative_float_delta

def load(p):return json.loads(Path(p).read_text())
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
 return h.hexdigest()
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def checkpoint(R,stage,**kw):save(R/'checkpoints'/(str(time.time_ns())+'.json'),dict(stage=stage,observed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),**kw))
def main(R):
 start=time.monotonic();cfg=load(R/'config/spec_v1.json');H=Path(cfg['producer']);G=Path(cfg['graph']);A=Path(cfg['graph_auditor']);assert load(A/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0;accepted=load(A/'outputs/accepted_same_graph_delivery_v1.json');assert accepted['graph_terminal_SHA256']==sha(G/'provenance/terminal_v1.json')
 for name,digest in cfg['policy_SHA256'].items():assert sha(R/'scripts'/name)==digest
 import scipy.sparse as sp
 toy=sp.csr_matrix(np.asarray([[1,0,0],[0,1,1],[1,1,1]],float));f,a=coherent_flags(toy,['CD3D','TRAC','LCK']);assert not f[:2,GATE_ORDER.index('T')].any() and f[2,GATE_ORDER.index('T')]
 g=['CD3D','CD3E','CD7','LCK','LAT'];f,a,mode=assay_aware_flags(sp.csr_matrix(np.ones((1,5))),g);assert a['T'] and f[0,0] and not a['gamma_delta_T_RNA_candidate'] and 'alternate' in mode
 f,a,mode=assay_aware_flags(sp.csr_matrix([[1,1,1,1,1,0]]),g+['TRAC']);assert a['T'] and not f[0,0] and mode=='original_chain_and_CD3_signaling_native_module'
 save(R/'provenance/coexpression_and_measured_negative_chain_smoke_v1.json',dict(same_cell_not_disjoint_cell_expression_pass=True,unmeasured_chain_alternate_mode_explicit_pass=True,measured_negative_chain_not_rescued_pass=True,gamma_delta_missing_chain_not_assigned=True,synthetic_not_biological_validation=True))
 hc=load(H/'config/spec_v1.json');J=Path(hc['joined']);sys.path.insert(0,str(J/'scripts'))
 from native_archive_collection_v1 import NativeArchiveCollection
 from native_annotation_review_v1 import native_measured_genes
 from native_panel_aggregate_v1 import gather_csr_rows
 from extract_full_source_local_marker_evidence_v1 import storage_gate
 import psutil
 storage_gate(R,cfg);assert psutil.virtual_memory().available>=cfg['minimum_available_RAM_bytes']
 native=NativeArchiveCollection(load(J/'outputs/joined_native_atlas_delivery_v1.json')['archive_manifest']);lineage=load(H/'outputs/source_input_lineage_v1.json');scope=load(H/'outputs/scope_and_common_genes_v1.json');N=scope['n_fit_records'];ng=accepted['clusters']['leiden_0_5'];assert len(GATE_ORDER)<=32
 with np.load(H/'outputs/private_cell_axes_v1.npz',allow_pickle=False) as f:axes={k:f[k] for k in f.files}
 graph_axes=G/'outputs/private_full_graph_cluster_axes_v1.tsv.gz';assert sha(graph_axes)==next(r['SHA256'] for r in accepted['all_output_file_SHA256'] if r['name']==graph_axes.name);labels=pd.read_csv(graph_axes,sep='\t',usecols=['leiden_0_5']).iloc[:,0].to_numpy(np.int32);assert len(labels)==N
 totalRNA=np.load(H/'outputs/private_analysis_RNA_totals_v1.npy',mmap_mode='r');flags_path=R/'outputs/private_all_cell_coherent_lineage_bitmask_v1.npy';bits=np.lib.format.open_memmap(flags_path,mode='w+',dtype=np.uint32,shape=(N,));observed=np.zeros((ng,len(GATE_ORDER)),np.int64);measured=np.zeros_like(observed);positive=np.zeros_like(observed);source_rows=[];source_meta=[]
 wanted=set(ACCESSORY+SIGNAL+CHAIN+['CD3D','CD3E'])|{g for modules in GATES.values() for m in modules for g in m}
 for code,r in enumerate(lineage):
  source=r['source'];assert source==scope['sources'][code];wh=np.flatnonzero(axes['source_code']==code);assert len(wh)==r['n_fit'];z=native.store(source);names,columns=native_measured_genes(native,source);take_gene=[i for i,g in enumerate(names) if g in wanted];genes=[names[i] for i in take_gene];physical=np.asarray(columns)[take_gene];delta=z['layers/ambient_removed_counts'] if 'ambient_removed_counts' in z['layers'] else None;count=np.zeros((ng,len(GATE_ORDER)),np.int64);n=np.bincount(labels[wh],minlength=ng)
  units=native.frozen_expression_units(source)[axes['source_native_position'][wh]].astype(str);missing=np.isin(np.char.lower(np.char.strip(units)),['','nan','none','na','n/a','unknown','__missing__','__source_value_missing__'])
  assessed=None;mode=None
  for lo in range(0,len(wh),cfg['read_block_size']):
   ix=wh[lo:lo+cfg['read_block_size']];raw=gather_csr_rows(z[r['counts_key']],axes['native_row'][ix],max_span=8192).astype(np.float64).tocsr();np.testing.assert_array_equal(np.asarray(raw.sum(axis=1)).ravel(),axes['native_total_RNA'][ix]);assert np.isfinite(raw.data).all() and (raw.data>=0).all() and np.equal(raw.data,np.floor(raw.data)).all();x=raw
   if delta is not None:x=(raw-gather_nonnegative_float_delta(delta,axes['native_row'][ix])).tocsr()
   x.eliminate_zeros();assert np.isfinite(x.data).all() and x.data.min(initial=0)>=-1e-8;x.data[x.data<0]=0;x.eliminate_zeros();np.testing.assert_allclose(np.asarray(x.sum(axis=1)).ravel(),totalRNA[ix],rtol=1e-10,atol=1e-7);panel=x[:,physical].tocsr();flags,assessed=coherent_flags(panel,genes);_,_,mode=assay_aware_flags(panel,genes);encoded=np.zeros(len(ix),np.uint32)
   for j,key in enumerate(GATE_ORDER):
    count[:,j]+=np.bincount(labels[ix],weights=flags[:,j],minlength=ng).astype(np.int64);encoded|=flags[:,j].astype(np.uint32)<<j
   bits[ix]=encoded
  assess=np.asarray([assessed[g] for g in GATE_ORDER],bool);measured+=n[:,None]*assess;observed+=n[:,None];positive+=count
  complete={g:all(all(x in genes for x in module) for module in modules) for g,modules in GATES.items()}
  for k in np.flatnonzero(n):
   sel=labels[wh]==k;known=units[sel][~missing[sel]];ct=collections.Counter(known.tolist());source_rows.append(dict(source=source,cluster=str(k),n_cells=int(n[k]),positive_gate_counts=dict(zip(GATE_ORDER,count[k].tolist())),assessed=assessed,all_original_marker_alternatives_measured=complete,T_evidence_mode=mode,source_cluster_RNA_proposal=annotate_cluster(int(n[k]),dict(zip(GATE_ORDER,count[k].tolist())),assessed),n_declared_expression_units_not_verified_patients_or_libraries=len(ct),n_missing_declared_expression_unit=int(missing[sel].sum()),largest_known_expression_unit_fraction=float(max(ct.values())/len(known)) if ct else None))
  source_meta.append(dict(source=source,n_fit=len(wh),measured_native_panel_gene_names=genes,native_panel_physical_columns=physical.tolist(),assessed=assessed,all_original_marker_alternatives_measured=complete,T_evidence_mode=mode,counts_key=r['counts_key'],adopted_ambient_delta_applied=delta is not None))
  native._fields.clear();bits.flush();storage_gate(R,cfg);checkpoint(R,'native_same_cell_cluster_modules_complete',source=source,n_fit=len(wh))
 assert int(np.bincount(labels,minlength=ng).sum())==N;assert np.all(positive<=measured)
 aggregate=[]
 for k in range(ng):
  entries=[r for r in source_rows if r['cluster']==str(k)]
  for j,key in enumerate(GATE_ORDER):
   eligible=[r for r in entries if r['assessed'][key]];adequate=[r for r in eligible if r['n_cells']>=30]
   aggregate.append(dict(cluster=str(k),gate=key,n_cluster_cells=int(observed[k,j]),n_module_assessable_cells=int(measured[k,j]),n_positive_coherent_cells=int(positive[k,j]),positive_fraction_among_assessable=float(positive[k,j]/measured[k,j]) if measured[k,j] else None,n_assessable_sources=len(eligible),n_adequate_sources=len(adequate),source_equal_mean_fraction=float(np.mean([r['positive_gate_counts'][key]/r['n_cells'] for r in eligible])) if eligible else None,adequate_source_equal_mean_fraction=float(np.mean([r['positive_gate_counts'][key]/r['n_cells'] for r in adequate])) if adequate else None,not_a_global_final_celltype_or_probability=True))
 bits.flush();del bits;pd.DataFrame(aggregate).to_csv(R/'outputs/joint_cluster_coherent_module_coverage_v1.tsv',sep='\t',index=False);save(R/'outputs/private_source_cluster_coherent_proposals_v1.json',source_rows);save(R/'outputs/private_source_measurement_modes_v1.json',source_meta)
 outputs=[dict(path=p.name,bytes=p.stat().st_size,SHA256=sha(p)) for p in (R/'outputs').iterdir() if p.is_file()];save(R/'outputs/joint_cluster_coexpression_manifest_v1.json',dict(n_fit=N,n_sources=len(lineage),n_clusters=ng,gate_order=GATE_ORDER,outputs=outputs,policy_SHA256=cfg['policy_SHA256'],input_axes_SHA256=sha(H/'outputs/private_cell_axes_v1.npz'),graph_axes_SHA256=sha(graph_axes),analysis_layer='inherited_adopted_analysis_counts_no_rounding',missing_module_denominators_preserved=True,existing_policy_is_exploratory_RNA_proposal_not_truth=True,no_global_final_identity_or_malignancy_promoted=True,expression_units_not_verified_patients=True));storage_gate(R,cfg)
 save(R/'provenance/terminal_v1.json',dict(status='FULL_NATIVE_SAME_CELL_LINEAGE_MODULES_ON_NEW_JOINT_CLUSTERS_COMPLETE_REQUIRES_INDEPENDENT_READBACK_AND_BIOLOGICAL_REVIEW',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),n_fit=N,n_sources=len(lineage),n_clusters=ng,wall_seconds=time.monotonic()-start,manifest_SHA256=sha(R/'outputs/joint_cluster_coexpression_manifest_v1.json')))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
