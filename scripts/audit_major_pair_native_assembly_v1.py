"""Independent fine-unit row sums, all major RNA vectors and readiness records."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,gc
import numpy as np,pandas as pd
from scipy.sparse import load_npz

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for block in iter(lambda:f.read(8388608),b''):h.update(block)
 return h.hexdigest()
def main():
 R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);P=R/'immune_same_name_family_source_review_v1';U=AR/'20261010_uniform_expression_annotation_v1';B=AR/'20261010_expression_annotation_downstream_v1';O=P/'major_pair_native_assembly_v1';planpath=Path(os.environ['PANCANCER_MAJOR_NATIVE_PLAN']);plan=json.loads(planpath.read_text());assert sha(__file__)==plan['auditor_sha256'];inv=json.loads((O/'inventory.json').read_text());assert inv['plan_sha256']==sha(planpath) and not (O/'independent_audit.json').exists()
 for n,h in inv['outputs'].items():assert sha(O/n)==h
 original=pd.read_csv(U/'all_node_number_catalog.tsv',sep='\t');current=pd.read_csv(B/'downstream_review_v2/number_to_annotation_dictionary.tsv',sep='\t').set_index('node_id');scope=pd.read_csv(P/'patient_block_provenance_overlay_v1/private_all1290_patient_provenance_overlay.tsv',sep='\t',keep_default_na=False).set_index('sampleID');counts=pd.read_csv(B/'inputs/sample_major_subcluster_counts.tsv',sep='\t');counts=counts[counts.sampleID.isin(scope.index)];n=counts.groupby(['sampleID','majorCluster']).n_cells.sum();groups=pd.read_csv(O/'all276_major_direction_contexts.tsv',sep='\t');assert len(groups)==276 and groups.group_id.is_unique;resource=pd.read_csv(Path(os.environ['PANCANCER_FROZEN_LR_RESOURCE']),sep='\t');lrgenes=set(v for pair in resource.itertuples(index=False,name=None) for gene in pair for v in gene.split('_'));breast=pd.read_csv(P/'GSE161529_official_metadata_recon_v2/private_all69_official_metadata_exact_mapping.tsv',sep='\t',dtype=str,keep_default_na=False).set_index('sampleID');reports=[];nunits=0;ngroups=0;ncases=0
 for sid,gg in groups.groupby('gse_id'):
  root=U/'by_gse'/sid;units=pd.read_csv(root/'units.tsv',sep='\t');genes=pd.read_csv(root/'genes.tsv',sep='\t').gene.to_numpy(dtype=str);meta=pd.read_csv(O/'private_source_major_arrays'/(sid+'_units.tsv'),sep='\t',keep_default_na=False);axes=dict(np.load(O/'private_source_major_arrays'/(sid+'_axes.npz'),allow_pickle=False));mean=load_npz(O/'private_source_major_arrays'/(sid+'_meanlog.npz')).tocsr();pct=load_npz(O/'private_source_major_arrays'/(sid+'_pct_positive.npz')).tocsr();assert np.array_equal(axes['genes'],genes);logs=load_npz(root/'sum_log1p_CP10K.npz');positive=load_npz(root/'positive_counts.npz');keys={}
  for row in units.itertuples():
   if row.candidate_eligible and row.sample_uid in scope.index:
    major=current.loc[original.iloc[row.node_index].node_id,'majorCluster'];keys.setdefault((row.sample_uid,major),[]).append(row.unit_row)
  keyindex={}
  for row in meta.itertuples():
   key=(row.sampleID,row.majorCluster);take=keys[key];cells=int(units.iloc[take].n_cells.sum());assert cells==row.n_cells==int(n.loc[key]);lm=logs[take].sum(axis=0)/cells;pm=positive[take].sum(axis=0)/cells;assert np.allclose(np.asarray(lm).ravel(),mean[row.unit_row].toarray().ravel(),rtol=1e-10,atol=1e-10);assert np.allclose(np.asarray(pm).ravel(),pct[row.unit_row].toarray().ravel(),rtol=1e-10,atol=1e-10);nonlr=~pd.Index(genes).isin(lrgenes);assert np.isclose(np.asarray(pm).ravel()[nonlr].sum(),axes['mean_detected_nonLR_genes'][row.unit_row]);assert np.isclose(np.asarray(lm).ravel()[nonlr].sum(),axes['mean_sum_log_nonLR_genes'][row.unit_row]);keyindex[key]=row.unit_row;nunits+=1
  for row in gg.itertuples():
   f=dict(np.load(O/'private_group_arrays'/(row.group_id+'.npz'),allow_pickle=False));samples=f['samples'];mm=scope.loc[samples];assert mm.analysis_eligible.all() and mm.stratum.eq(row.stratum).all() and len(samples)==row.n_endpoint_libraries and mm.effective_resampling_block.nunique()==row.n_analysis_blocks;assert np.array_equal(f['patient_blocks'],mm.effective_resampling_block.to_numpy(dtype=str));assert np.array_equal(f['reviewed_patient'],mm.unit_kind.eq('reviewed_patient'));
   if sid=='GSE161529':
    bm=breast.loc[samples];assert np.array_equal(f['corrected_material_for_next_inference'],bm.official_material_role);tumor=bm.current_class.eq('primary_or_unspecified_tumor').to_numpy();expected=np.where(tumor,('GSE161529:tumor_author:'+bm.author_patient_token).to_numpy(dtype=str),mm.effective_resampling_block.to_numpy(dtype=str));assert np.array_equal(f['corrected_blocks_for_next_inference'],expected)
   else:assert np.array_equal(f['corrected_blocks_for_next_inference'],f['patient_blocks'])
   si=[keyindex[s,row.source] for s in samples];ti=[keyindex[s,row.target] for s in samples];assert np.array_equal(f['source_n'],meta.iloc[si].n_cells) and np.array_equal(f['target_n'],meta.iloc[ti].n_cells);assert (f['source_n']>=20).all() and (f['target_n']>=20).all();gi={g:i for i,g in enumerate(genes)};L=mean[si].toarray();T=mean[ti].toarray();LP=pct[si].toarray();TP=pct[ti].toarray()
   for j,(lig,rec) in enumerate(resource.itertuples(index=False,name=None)):
    measured=all(v in gi for v in lig.split('_')+rec.split('_'));assert measured==bool(f['native_measured'][j]);colsL=[gi[v] for v in lig.split('_')] if measured else [];colsR=[gi[v] for v in rec.split('_')] if measured else []
    for array,expected in [('ligand_meanlog',L[:,colsL]),('receptor_meanlog',T[:,colsR]),('ligand_min_positive_fraction',LP[:,colsL]),('receptor_min_positive_fraction',TP[:,colsR])]:
     if not measured:assert np.isnan(f[array][:,j]).all()
     else:assert np.allclose(np.minimum.reduce(expected.T),f[array][:,j],rtol=1e-10,atol=1e-10)
   ngroups+=1
  reports.append(dict(gse_id=sid,major_units=len(meta),groups=len(gg),status='PASS_SEPARATE_FINE_UNIT_ROW_SUMS_AND_ALL_COMPLEX_VALUES'));print(json.dumps(reports[-1]),flush=True);del logs,positive,mean,pct,L,T,LP,TP;gc.collect()
 # Complete case key and readiness check, streamed chunks, no individual exports.
 index=groups.set_index('group_id');cache={};group_order={g:i for i,g in enumerate(groups.group_id)};seen=np.zeros((276,4624),bool)
 for d in pd.read_csv(O/'all1276224_native_LR_case_readiness.tsv.gz',sep='\t',chunksize=50000,keep_default_na=False):
  for gid,g in d.groupby('group_id',sort=False):
   if gid not in cache:cache={gid:dict(np.load(O/'private_group_arrays'/(gid+'.npz'),allow_pickle=False))}
   f=cache[gid]
   for row in g.itertuples():
    j=int(row.case_id.rsplit('LR',1)[1])-1;assert row.case_id.startswith(gid+'_LR') and not seen[group_order[gid],j];seen[group_order[gid],j]=True;assert (row.ligand,row.receptor)==tuple(resource.iloc[j]);ok=bool(f['native_measured'][j]);assert row.native_all_subunits_measured==ok
    if ok:assert int(row.n_positive_ligand_libraries)==int((f['ligand_meanlog'][:,j]>0).sum()) and int(row.n_positive_receptor_libraries)==int((f['receptor_meanlog'][:,j]>0).sum());assert row.ligand_variable==bool(np.ptp(f['ligand_meanlog'][:,j])>1e-12) and row.receptor_variable==bool(np.ptp(f['receptor_meanlog'][:,j])>1e-12)
    else:assert not row.readiness_only and row.n_positive_ligand_libraries=='' and row.n_positive_receptor_libraries==''
    assert row.readiness_only==bool(ok and row.ligand_variable and row.receptor_variable);ncases+=1
 assert ngroups==276 and ncases==1276224 and seen.all();out=dict(status='PASS_ALL276_MAJOR_DIRECTIONS_1276224_CASES_AND_NATIVE_FINE_UNIT_AGGREGATION',observed_UTC=datetime.now(timezone.utc).isoformat(),auditor_sha256=sha(__file__),producer_inventory_sha256=sha(O/'inventory.json'),native_major_units_independently_checked=nunits,directions=ngroups,case_records=ncases,source_reports=reports,physical_binding_or_induced_response_validation=False,new_p=False,fresh_raw_cell_recount_this_audit=False);(O/'independent_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
if __name__=='__main__':main()
