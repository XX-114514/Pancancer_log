"""V2 accepts losslessly integral decimal count serialization; science unchanged."""
import argparse,json
from pathlib import Path
from datetime import datetime,timezone
import h5py,numpy as np,pandas as pd
from scipy.sparse import csr_matrix
from scipy.spatial import KDTree
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource

def main(a):
 resource(a.root,8*1024**2);assert not a.output.exists();a.output.mkdir(mode=0o700)
 term=json.loads((a.producer/'terminal_public.json').read_text());assert term['n_slots']==83232
 for name,checksum in term['files'].items():assert sha(a.producer/name)==checksum
 plan=json.loads((a.producer/'private_bound_plan_before_molecular_counts.json').read_text())
 refpath=a.root/'LR_binding_resource_review_v1/all4624_reference_evidence_status.tsv';assert sha(refpath)==plan['reference_sha256']
 reference=pd.read_csv(refpath,sep='\t',keep_default_na=False);assert len(reference)==4624
 result=pd.read_csv(a.producer/'private_all83232_spatial_molecular_opportunities.tsv.gz',sep='\t',keep_default_na=False);assert len(result)==83232 and not result.duplicated(['source_index','direction','ligand','receptor']).any()
 coord=a.root/'GSE299886_original_coordinate_sidecars_20261018_v1';assert sha(coord/'private_terminal.json')==plan['coordinate_terminal_sha256']
 reports=json.loads((coord/'private_all_source_coordinate_reconstruction_receipt.json').read_text())
 schemas=json.loads((a.root/'CM02_existing_targeted_spatial_schema_20261018_v2/private_all_eligible_targeted_input_schema.json').read_text())
 state=Path(json.loads((a.root/'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root'])
 ax=state/'state_spatial_proxy_calibration_inputs_20261016_v1/fixed_all14_panel_gene_axes.json';assert sha(ax)==plan['axes_sha256'];axes=json.loads(ax.read_text())
 fp=a.root.parents[1]/'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json';assert sha(fp)==plan['fixed_panels_sha256'];panels=json.loads(fp.read_text())['panels']
 checked=0;raw_values=0;differences=[]
 for number,r in enumerate(reports):
  schema=next(x for x in schemas if x['sample_id']==r['sample_id']);mapping=next(x for x in axes if x['sample_id']==r['sample_id'])['symbol_to_feature_index'];matrix=Path(schema['analysis_file']);assert sha(matrix)==r['matrix_SHA256']
  z=np.load(coord/r['output'],allow_pickle=False);xy=z['xy_source'];n=len(xy);det=np.empty((n,500),bool)
  with h5py.File(matrix,'r') as h:
   # Independent CSR reconstruction, different chunk size and sparse operator.
   for start in range(0,n,3072):
    stop=min(n,start+3072);p=h['X/indptr'][start:stop+1];lo,hi=int(p[0]),int(p[-1])
    b=csr_matrix((h['X/data'][lo:hi],h['X/indices'][lo:hi],p-lo),shape=(stop-start,500));b.sum_duplicates()
    v=b.toarray();assert np.isfinite(v).all() and (v>=0).all();det[start:stop]=(v>0);raw_values+=v.size
  captures={k:sum(det[:,mapping[g]].astype(np.int8) for g in panels[k])>=2 for k in ['pan_T','regulatory_associated','cytotoxic']}
  t,rr,cc=(captures[k] for k in ['pan_T','regulatory_associated','cytotoxic'])
  ameasured={'HLA-DRA','CD74'}<=set(mapping)
  masks={'R':t&rr&~cc,'C':t&cc&~rr,'A':(~t)&det[:,mapping['HLA-DRA']]&det[:,mapping['CD74']] if ameasured else np.zeros(n,bool)}
  indices={k:np.where(m)[0] for k,m in masks.items()};complexes=set(reference.ligand)|set(reference.receptor)
  complete={k:set(k.split('_'))<=set(mapping) for k in complexes}
  positive={}
  for group,ix in indices.items():
   positive[group]={}
   for c in complexes:
    if complete[c]:
     selected=np.ones(len(ix),bool)
     for g in c.split('_'):selected &= det[ix,mapping[g]]
     positive[group][c]=selected
  unit=np.median(KDTree(xy).query(xy,k=2)[0][:,1]);radius=5*unit
  for left,right in [('R','C'),('C','R'),('R','A'),('A','R'),('C','A'),('A','C')]:
   ix,iy=indices[left],indices[right];ns,nt=len(ix),len(iy)
   balls=KDTree(xy[ix]).query_ball_tree(KDTree(xy[iy]),radius) if ns and nt else [[] for _ in ix]
   ei=np.repeat(np.arange(ns,dtype=np.int32),[len(v) for v in balls]);ej=np.array([j for v in balls for j in v],np.int32);ne=len(ei)
   if ne:assert np.all(np.linalg.norm(xy[ix[ei]]-xy[iy[ej]],axis=1)<=radius*(1+1e-12)) and not np.any(ix[ei]==iy[ej])
   saved=result[(result.source_index==number)&result.direction.eq(left+'_to_'+right)].set_index(['ligand','receptor']);assert len(saved)==4624
   for entry in reference.itertuples(index=False):
    s=saved.loc[(entry.ligand,entry.receptor)];missing=sorted(set(entry.ligand.split('_')+entry.receptor.split('_'))-set(mapping));measured=not missing
    assert bool(s.subunits_measured)==measured and s.missing_subunits==';'.join(missing)
    assert int(s.sender_cells)==ns and int(s.receiver_cells)==nt and int(s.directed_group_edges)==ne
    assert s.reference_status==entry.reference_status and int(s.n_reference_supported)==entry.n_supported
    lp=int(positive[left][entry.ligand].sum()) if complete[entry.ligand] else None
    rp=int(positive[right][entry.receptor].sum()) if complete[entry.receptor] else None
    for value,column in [(lp,'ligand_all_subunits_positive'),(rp,'receptor_all_subunits_positive')]:
     assert (s[column]=='' if value is None else int(float(s[column]))==value)
    if measured:
     observed=int(np.logical_and(positive[left][entry.ligand][ei],positive[right][entry.receptor][ej]).sum())
     assert int(float(s.ligand_to_complete_receptor_edges))==observed
     rate=observed/ne if ne else None;global_rate=lp/ns*rp/nt if ns and nt else None
     ratio=rate/global_rate if rate is not None and global_rate else None
     for value,column in [(rate,'edge_positive_fraction'),(global_rate,'global_group_pair_positive_fraction'),(ratio,'edge_over_global_pair_rate')]:
      if value is None:assert s[column]==''
      else:
       err=abs(float(s[column])-value);differences.append(err);assert np.isclose(float(s[column]),value,rtol=1e-10,atol=1e-12)
    else:assert s.ligand_to_complete_receptor_edges==''
    reasons=[]
    if not measured:reasons.append('unmeasured_subunits')
    if 'A' in [left,right] and not ameasured:reasons.append('unmeasured_HLAII_group')
    if min(ns,nt)<30:reasons.append('insufficient_group_cells')
    if measured and min(lp,rp)<10:reasons.append('insufficient_joint_positive_cells')
    if not ne:reasons.append('no_directed_neighbor_edges')
    assert s.reason==';'.join(reasons) and bool(s.measurement_qualified)==(not reasons)
    assert bool(s.curated_reference_rejected)==(entry.n_curators_rejected>0 and entry.n_supported==0)
    checked+=1
  del det,positive
  print(json.dumps(dict(phase='all_source_resource_positions_replayed',source_index=number,checked=checked)),flush=True)
 assert checked==83232
 receipt=dict(status='ALL83232_SPATIAL_MOLECULAR_MEASUREMENTS_INDEPENDENTLY_REPLAYED',finished_utc=datetime.now(timezone.utc).isoformat(),n_slots=checked,n_raw_values_re_read=raw_values,maximum_fraction_error=max(differences,default=0),producer_terminal_sha256=sha(a.producer/'terminal_public.json'),auditor_sha256=sha(__file__),no_new_P_or_FDR=True,scope='Independent original RNA/group/joint-subunit/neighbor-edge/all-resource replay; not independent biological validation or a conditional spatial enrichment test.')
 (a.output/'independent_audit_public.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))

if __name__=='__main__':
 p=argparse.ArgumentParser()
 for name in ['root','producer','output']:p.add_argument('--'+name,type=Path,required=True)
 main(p.parse_args())
