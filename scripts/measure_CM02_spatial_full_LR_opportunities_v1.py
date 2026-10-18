"""All-resource spatial RNA opportunities; descriptive, no new P or labels."""
import argparse, hashlib, json
from pathlib import Path
from datetime import datetime, timezone
import h5py
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from repair_GSE299886_original_coordinate_sidecars_v1 import dense_rows
from qualify_GSE299886_fixed_panels_negative_controls_v1 import verify

DIRECTIONS=[('R','C'),('C','R'),('R','A'),('A','R'),('C','A'),('A','C')]

def main(a):
 root=a.root;budget=resource(root,64*1024**2);assert not a.output.exists();a.output.mkdir(mode=0o700)
 contract=json.loads(a.plan.read_text());assert contract['slots_expected']==83232 and not contract['new_P_or_FDR']
 refpath=root/'LR_binding_resource_review_v1/all4624_reference_evidence_status.tsv'
 ref=pd.read_csv(refpath,sep='\t',keep_default_na=False);assert len(ref)==4624 and not ref.duplicated(['ligand','receptor']).any()
 fp=root.parents[1]/'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json';panels=json.loads(fp.read_text())['panels']
 coord=root/'GSE299886_original_coordinate_sidecars_20261018_v1';verify(coord)
 reports=json.loads((coord/'private_all_source_coordinate_reconstruction_receipt.json').read_text());assert len(reports)==3
 schema=json.loads((root/'CM02_existing_targeted_spatial_schema_20261018_v2/private_all_eligible_targeted_input_schema.json').read_text())
 state=Path(json.loads((root/'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root'])
 axespath=state/'state_spatial_proxy_calibration_inputs_20261016_v1/fixed_all14_panel_gene_axes.json';axes=json.loads(axespath.read_text())
 bound=dict(contract=contract,reference_sha256=sha(refpath),coordinate_terminal_sha256=sha(coord/'private_terminal.json'),fixed_panels_sha256=sha(fp),axes_sha256=sha(axespath),producer_sha256=sha(__file__),budget=budget,started_utc=datetime.now(timezone.utc).isoformat())
 (a.output/'private_bound_plan_before_molecular_counts.json').write_text(json.dumps(bound,indent=2)+'\n')
 rows=[];source=[]
 for number,rec in enumerate(reports):
  sid=rec['sample_id'];s=next(x for x in schema if x['sample_id']==sid);axis=next(x for x in axes if x['sample_id']==sid)
  matrix=Path(s['analysis_file']);before=(matrix.stat().st_size,matrix.stat().st_mtime_ns);assert sha(matrix)==rec['matrix_SHA256']==axis['input_sha256']
  mapping=axis['symbol_to_feature_index'];assert all(g in mapping for k in ['pan_T','regulatory_associated','cytotoxic'] for g in panels[k])
  z=np.load(coord/rec['output'],allow_pickle=False);xy=z['xy_source'];n=len(xy);det=np.zeros((n,500),bool)
  with h5py.File(matrix,'r') as h:
   for start in range(0,n,2048):
    stop=min(n,start+2048);v=dense_rows(h['X'],start,stop,500);assert np.isfinite(v).all() and (v>=0).all();det[start:stop]=v>=1
  T=det[:,[mapping[g] for g in panels['pan_T']]].sum(1)>=2
  R=det[:,[mapping[g] for g in panels['regulatory_associated']]].sum(1)>=2
  C=det[:,[mapping[g] for g in panels['cytotoxic']]].sum(1)>=2
  a_measured=all(g in mapping for g in ['HLA-DRA','CD74'])
  A=(~T)&det[:,mapping['HLA-DRA']]&det[:,mapping['CD74']] if a_measured else np.zeros(n,bool)
  ids={'R':np.flatnonzero(T&R&~C),'C':np.flatnonzero(T&C&~R),'A':np.flatnonzero(A)}
  complexes=set(ref.ligand)|set(ref.receptor);complete={c:all(g in mapping for g in c.split('_')) for c in complexes}
  positive={k:{c:np.all(det[np.ix_(ix,[mapping[g] for g in c.split('_')])],axis=1) for c in sorted(complexes) if complete[c]} for k,ix in ids.items()}
  nn=cKDTree(xy).query(xy,k=2,workers=1)[0][:,1];unit=float(np.median(nn));assert np.isfinite(unit) and unit>0
  summary=dict(source_index=number,n_cells=n,T_capture_cells=int(T.sum()),R_exclusive_cells=len(ids['R']),C_exclusive_cells=len(ids['C']),double_positive_excluded=int((T&R&C).sum()),HLAII_nonT_cells=len(ids['A']),HLAII_markers_measured=a_measured,unit_source_coordinates=unit,physical_units_unverified=True)
  for left,right in DIRECTIONS:
   ix,iy=ids[left],ids[right];ns,nt=len(ix),len(iy)
   if ns and nt:
    distance=cKDTree(xy[ix]).sparse_distance_matrix(cKDTree(xy[iy]),5*unit,output_type='coo_matrix');adj=csr_matrix((np.ones(len(distance.data),np.int32),(distance.row,distance.col)),shape=(ns,nt));adj.sum_duplicates()
   else:adj=csr_matrix((ns,nt),dtype=np.int32)
   assert (adj.data==1).all();edge_count=int(adj.sum());hits={}
   receptors=sorted(c for c in set(ref.receptor) if complete[c])
   for start in range(0,len(receptors),16):
    cc=receptors[start:start+16];rp=np.column_stack([positive[right][c] for c in cc]).astype(np.int32) if nt else np.zeros((0,len(cc)),np.int32)
    near=adj@rp
    for j,c in enumerate(cc):
     for l in ref.loc[ref.receptor.eq(c),'ligand']:
      if complete[l]:hits[l,c]=int(near[:,j][positive[left][l]].sum(dtype=np.int64))
   for r in ref.itertuples(index=False):
    missing=sorted(set(r.ligand.split('_')+r.receptor.split('_'))-set(mapping));measured=not missing
    lp=int(positive[left][r.ligand].sum()) if complete[r.ligand] else None;rp=int(positive[right][r.receptor].sum()) if complete[r.receptor] else None
    observed=hits.get((r.ligand,r.receptor));opportunity=observed/edge_count if observed is not None and edge_count else None
    global_rate=lp/ns*rp/nt if measured and ns and nt else None
    reasons=[]
    if not measured:reasons.append('unmeasured_subunits')
    if 'A' in [left,right] and not a_measured:reasons.append('unmeasured_HLAII_group')
    if min(ns,nt)<30:reasons.append('insufficient_group_cells')
    if measured and min(lp,rp)<10:reasons.append('insufficient_joint_positive_cells')
    if not edge_count:reasons.append('no_directed_neighbor_edges')
    rows.append(dict(source_index=number,direction=left+'_to_'+right,ligand=r.ligand,receptor=r.receptor,reference_status=r.reference_status,n_reference_supported=r.n_supported,
     sender_cells=ns,receiver_cells=nt,subunits_measured=measured,missing_subunits=';'.join(missing),ligand_all_subunits_positive=lp,receptor_all_subunits_positive=rp,
     directed_group_edges=edge_count,ligand_to_complete_receptor_edges=observed,edge_positive_fraction=opportunity,global_group_pair_positive_fraction=global_rate,
     edge_over_global_pair_rate=opportunity/global_rate if opportunity is not None and global_rate else None,measurement_qualified=not reasons,reason=';'.join(reasons),
     curated_reference_rejected=(r.n_curators_rejected>0 and r.n_supported==0),no_new_P=True,no_binding_or_signaling_claim=True))
  assert before==(matrix.stat().st_size,matrix.stat().st_mtime_ns) and sha(matrix)==rec['matrix_SHA256']
  summary['matrix_sha256']=rec['matrix_SHA256'];source.append(summary);del det,positive
  print(json.dumps(dict(phase='assay_complete',source_index=number,n_resource_slots=len(rows))),flush=True)
 result=pd.DataFrame(rows);assert len(result)==83232
 result.to_csv(a.output/'private_all83232_spatial_molecular_opportunities.tsv.gz',sep='\t',index=False)
 (a.output/'private_all3_source_measurement_summary.json').write_text(json.dumps(source,indent=2)+'\n')
 size=sum(p.stat().st_size for p in a.output.iterdir() if p.is_file());assert size<64*1024**2
 terminal=dict(status='ALL83232_MOLECULAR_MEASUREMENT_SLOTS_COMPLETE_INDEPENDENT_AUDIT_PENDING',finished_utc=datetime.now(timezone.utc).isoformat(),n_slots=len(result),no_P_or_FDR=True,no_old_families_reopened=True,output_bytes=size,files={p.name:sha(p) for p in a.output.iterdir() if p.is_file()})
 (a.output/'terminal_public.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps({k:v for k,v in terminal.items() if k!='files'}),flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser()
 for name in ['root','plan','output']:p.add_argument('--'+name,type=Path,required=True)
 main(p.parse_args())
