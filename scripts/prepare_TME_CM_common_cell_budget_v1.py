"""Identical cell draws for all annotation schemes on a frozen CM panel.

The budget is applied within each material and broad parent, retaining every
sampled fallback cell in the denominator. No scheme-specific resampling occurs.
"""
import argparse,datetime,hashlib,json,sys
from pathlib import Path
import numpy as np,pandas as pd

def load(p):return json.loads(Path(p).read_text())
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
 return h.hexdigest()
def bound(b):
 p=Path(b['path']);assert sha(p)==b['SHA256'];return p
def binding(p):return dict(path=str(p),SHA256=sha(p))
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False)
def main(R):
 cfg=load(R/'config/spec_v1.json');sys.path.insert(0,str(Path(cfg['joined'])/'scripts'));from extract_full_source_local_marker_evidence_v1 import storage_gate
 storage_gate(R,cfg);P=Path(cfg['preparation']);common=np.load(P/'outputs/private_common_assessed_cell_mask_v1.npy');C=Path(load(P/'config/spec_v1.json')['unit_closure']);unit=np.load(bound(load(C/'outputs/accepted_fit_qualification_unit_bridge_v3.json')['unit_codes']));nu=int(unit.max())+1
 items=[];sample_ids=None;groups=None;used_units=None;shared_parent=None
 for s in cfg['schemes']:
  root=Path(s['root']);assert load(root/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
  ac=load(root/'config/spec_v1.json')['input_acceptance'];a=load(bound(ac));membership=load(bound(a['annotation_acceptance']));p=np.load(bound(membership['parent_codes']));n=np.load(bound(membership['node_codes']));train=np.load(bound(a['training_mask_binding']));meta=pd.read_csv(bound(a['metadata']),sep='\t');fq=pd.read_csv(bound(a['frequency']),sep='\t',index_col=0)
  assert fq.columns[1:].tolist()==meta.sampleID.tolist()
  parents=[v['code'] for v in membership['parents'] if v['name'] in a['requested_parent_panel']['parents']];mask=common&train&np.isin(unit,meta.qualification_unit_code.to_numpy())&np.isin(p,parents)
  if sample_ids is None:
   sample_ids=meta.sampleID.tolist();used_units=meta.qualification_unit_code.to_numpy();groups=mask.copy();shared_parent=p.copy()
  assert sample_ids==meta.sampleID.tolist();np.testing.assert_array_equal(mask,groups);np.testing.assert_array_equal(p[mask],shared_parent[mask])
  nodemap={v['name']:v['code'] for v in membership['nodes']};rows=np.array([nodemap[v] for v in fq.index],int)
  items.append(dict(id=s['id'],root=str(root),acceptance=ac,membership=membership,node=n,rows=rows,nodes=fq.index.tolist(),parents=fq.iloc[:,0].tolist()))
 parent_codes=np.unique(shared_parent[groups]);chunks=[]
 for u in used_units:
  for p in parent_codes:
   rows=np.flatnonzero(groups&(unit==u)&(shared_parent==p));assert len(rows)>=cfg['cells_per_parent_unit'];chunks.append((int(u),int(p),rows))
 rng=np.random.default_rng(cfg['seed']);repeats=[];scalar_probes=[]
 for repeat in range(cfg['n_repeats']):
  selected=np.concatenate([rng.choice(rows,cfg['cells_per_parent_unit'],replace=False) for u,p,rows in chunks]);assert len(selected)==len(chunks)*cfg['cells_per_parent_unit'] and len(np.unique(selected))==len(selected)
  draw=R/'outputs'/('private_repeat_'+str(repeat)+'_same_cell_rows_v1.npy');np.save(draw,selected);files=[]
  for item in items:
   n=item['node'];count=np.bincount(n[selected].astype(np.int64)*nu+unit[selected],minlength=len(item['membership']['nodes'])*nu).reshape(len(item['membership']['nodes']),nu)
   parent=np.bincount(shared_parent[selected].astype(np.int64)*nu+unit[selected],minlength=len(item['membership']['parents'])*nu).reshape(len(item['membership']['parents']),nu)
   for p in parent_codes:
    child=np.array([v['parent_code']==p for v in item['membership']['nodes']]);np.testing.assert_array_equal(count[child].sum(0),parent[p]);assert (parent[p,used_units]==cfg['cells_per_parent_unit']).all()
   pids=np.array([item['membership']['nodes'][j]['parent_code'] for j in item['rows']]);numerator=count[np.ix_(item['rows'],used_units)];denominator=parent[np.ix_(pids,used_units)];frequency=numerator/denominator;assert np.isfinite(frequency).all() and (frequency>=0).all() and (frequency<=1).all()
   fp=R/'outputs'/('private_repeat_'+str(repeat)+'_'+item['id']+'_frequency_v1.tsv.gz');frame=pd.DataFrame(frequency,index=item['nodes'],columns=sample_ids);frame.insert(0,'majorCluster',item['parents']);frame.to_csv(fp,sep='\t',index_label='subCluster',compression='gzip')
   if repeat==0:
    for positive in [True,False]:
     at=np.argwhere((numerator>0) if positive else (numerator==0))
     if len(at):
      j,k=at[0];u=used_units[k];nd=item['rows'][j];pa=pids[j];ss=selected[unit[selected]==u];num=np.count_nonzero(n[ss]==nd);den=np.count_nonzero(shared_parent[ss]==pa);assert frequency[j,k]==num/den;scalar_probes.append(dict(scheme=item['id'],case='observed_positive' if positive else 'observed_zero',direct_cell_recount_pass=True))
   files.append(dict(scheme=item['id'],frequency=binding(fp)))
  repeats.append(dict(repeat=repeat,n_selected_cells=len(selected),same_cell_draw=binding(draw),schemes=files));storage_gate(R,cfg)
 op=R/'outputs/accepted_common_cell_budget_inputs_v1.json';save(op,dict(schema='TME_CM_common_cell_budget_inputs_v1',schemes=[{k:i[k] for k in ['id','root','acceptance']} for i in items],n_analysis_units=len(used_units),n_parents=len(parent_codes),n_repeats=cfg['n_repeats'],cells_per_parent_unit=cfg['cells_per_parent_unit'],repeats=repeats,sample_ids=sample_ids,direct_scalar_recounts=scalar_probes,identical_draws_across_annotation_schemes=True,all_fallbacks_retained_in_complete_parent_denominator=True,subgroup_observed_zeros_not_filtered_or_filled_from_missing=True,annotation_models_not_refit=True,author_minmax_reestimated_in_each_budget_refit=True,clinical_outcomes_used=False))
 storage_gate(R,cfg);save(R/'provenance/terminal_v1.json',dict(status='SAME_CELL_BUDGET_CM_INPUTS_PREPARED_AND_RECOUNTED',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),acceptance_SHA256=sha(op)))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
