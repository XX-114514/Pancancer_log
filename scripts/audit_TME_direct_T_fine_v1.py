"""Independent full-scope replay, rebuilding fine cues from native quantities."""
import argparse,collections,datetime,gc,hashlib,json,sys,time
from pathlib import Path
import numpy as np

TYPES=['T','Treg','CD4_T','CD8_T','T_gamma_delta_RNA_candidate']
REASONS=['outside_coarse_T_target','direct_native_and_unique_reference_fine_RNA','gamma_RNA_candidate_requires_orthogonal_review','no_direct_fine_intersection','multiple_direct_fine_candidates','coarse_T_native_support_not_qualified','zero_native_RNA']
def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
 with Path(p).open('x') as f:json.dump(d,f,indent=2,allow_nan=False);f.write('\n')

def main(root):
 start=time.monotonic();P=Path(load(root/'config/spec_v1.json')['producer']);assert load(P/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0;pc=load(P/'config/spec_v1.json')
 for b in pc['bindings'].values():assert sha(b['path'])==b['SHA256']
 accepted=load(pc['working_acceptance']);assert sha(accepted['producer_manifest']['path'])==accepted['producer_manifest']['SHA256'];w=load(accepted['producer_manifest']['path']);wi={r['source']:r for r in load(w['index']['path'])};base=Path(w['index']['path']).parent
 F=Path(pc['candidate_producer']);fi={r['source']:r for r in load(F/'outputs/full_native_cell_identity_candidate_index_v1.json')};E=Path(pc['native_evidence']);ei={r['source']:r for r in load(E/'outputs/TME_full_native_evidence_index_v1.json')}
 ref=Path(load(F/'config/full_cell_candidate_spec_v1.json')['reference']);ri={r['source']:r for r in load(ref/'outputs/full_archive_reference_index_v1.json')};models=load(ref/'outputs/reference_models_v1.json');sys.path.insert(0,str(F/'scripts'));from scoped_low_reference_identity_v1 import scoped_low_identity
 refmap=np.asarray([sum(1<<j for j,n in enumerate(TYPES[1:]) if n in scoped_low_identity(label)[1]) for label in models['immune_low']['labels']],np.uint16)
 ip=P/'outputs/direct_T_fine_proposal_index_v1.json';assert sha(ip)==load(P/'provenance/terminal_v1.json')['index_SHA256'];pi=load(ip);assert {r['source'] for r in pi}==set(wi)==set(ei)
 totals=collections.Counter();labels=collections.Counter();why=collections.Counter();n_native=0
 for entry in pi:
  source=entry['source'];mp=P/'outputs'/entry['metadata'];assert sha(mp)==entry['metadata_SHA256'];m=load(mp);pa=P/'outputs'/m['arrays'];assert sha(pa)==m['arrays_SHA256']
  with np.load(pa,allow_pickle=False) as f:out={k:f[k] for k in f.files}
  wr=wi[source];wm=load(base/wr['metadata']);wa=base/wr['arrays'];assert sha(base/wr['metadata'])==wr['metadata_SHA256'] and sha(wa)==wr['arrays_SHA256']
  with np.load(wa,allow_pickle=False) as f:q={k:f[k] for k in ['native_row','source_native_position','consensus_major_code','consensus_celltype_code','native_family_bitmask','annotation_review_code','native_total_RNA','annotation_group_code']}
  er=ei[source];ep=E/'outputs'/er['metadata'];assert sha(ep)==er['metadata_SHA256'];em=load(ep);cp=em['native_context'];assert sha(cp['path'])==cp['SHA256']
  with np.load(cp['path'],allow_pickle=False) as f:x=f['native_gene_values'];mask=f['measured'];np.testing.assert_array_equal(f['native_row'],q['native_row'])
  by={g:j for j,g in enumerate(em['context_genes'])};zeros=np.zeros(len(x),bool)
  def positive(g):return x[:,by[g]]>0 if g in by and mask[by[g]] else zeros
  Tpanel=np.column_stack([positive(g) for g in ['CD3D','CD3E','TRAC','CD247']]).sum(axis=1)>=2
  cues=[positive('FOXP3')&(positive('IL2RA')|positive('CTLA4')),positive('CD4')&Tpanel,
        positive('CD8A')&positive('CD8B')&(positive('CD3D')|positive('CD3E')),positive('TRDC')&(positive('TRGC1')|positive('TRGC2'))&positive('CD3D')]
  own=sum((v.astype(np.uint16)*(1<<j) for j,v in enumerate(cues)),np.zeros(len(x),np.uint16))
  fr=fi[source];fm=F/'outputs'/fr['metadata'];fa=F/'outputs'/fr['arrays'];assert sha(fm)==fr['metadata_sha256'] and sha(fa)==fr['arrays_sha256'];assert load(fm)['fine_candidate_order'][:4]==TYPES[1:]
  rr=ri[source];rp=ref/'outputs'/rr['metadata'];ra=ref/'outputs'/rr['arrays'];assert sha(rp)==rr['metadata_sha256'] and sha(ra)==rr['arrays_sha256']
  pos=q['source_native_position']
  with np.load(ra,allow_pickle=False) as f:
   np.testing.assert_array_equal(q['native_row'],f['native_row'][pos]);low=f['immune_low__unique_code'][pos];reference=np.zeros(len(x),np.uint16);valid=low>=0;reference[valid]=refmap[low[valid]];assert not np.any(valid&~f['immune_reference_scope'][pos])
  with np.load(fa,allow_pickle=False) as f:
   np.testing.assert_array_equal(q['native_row'],f['native_row'][pos]);np.testing.assert_array_equal(own,f['own_fine_bitmask'][pos]&15);np.testing.assert_array_equal(reference,f['reference_fine_bitmask'][pos]&15)
   np.testing.assert_array_equal(out['own_fine_bitmask'],f['own_fine_bitmask'][pos]);np.testing.assert_array_equal(out['reference_fine_bitmask'],f['reference_fine_bitmask'][pos])
  for name in ['native_row','source_native_position']:np.testing.assert_array_equal(out[name],q[name])
  major=np.asarray(wm['candidate_dictionaries']['consensus_major'],object)[q['consensus_major_code']];fine=np.asarray(wm['candidate_dictionaries']['consensus_celltype'],object)[q['consensus_celltype_code']];target=(major=='T')&(fine=='T');support=(q['native_family_bitmask']==1)|(q['annotation_review_code']==6)
  lookup=np.zeros(16,np.uint8);pruned=np.zeros(16,np.uint8);status=np.zeros(16,np.uint8)
  for v in range(16):
   u=v&~2 if v&1 else v;pruned[v]=u;hits=[j+1 for j in range(4) if u&(1<<j)];lookup[v]=hits[0] if len(hits)==1 else 0;status[v]=4 if len(hits)>1 else 3
  raw=own&reference;direct=pruned[raw];proposal=np.where(target&support&(q['native_total_RNA']>0),lookup[raw],0).astype(np.uint8);reason=np.zeros(len(x),np.uint8);reason[target]=status[raw[target]];reason[target&~support]=5;reason[target&(q['native_total_RNA']<=0)]=6;reason[proposal>0]=1;reason[proposal==4]=2
  for name,value in [('proposed_fine_code',proposal),('fine_review_code',reason),('direct_intersection_bitmask',direct),('eligible_coarse_T',target)]:np.testing.assert_array_equal(out[name],value)
  changed=proposal>0;cc=dict(collections.Counter(np.asarray(TYPES,object)[proposal[changed]]));rrc=dict(collections.Counter(REASONS[int(v)] for v in reason[target]));assert cc==m['new_fine_counts'] and rrc==m['target_review_counts'];assert m['fine_dictionary']==TYPES and m['review_dictionary']==REASONS
  assert m['working_arrays_SHA256']==sha(wa) and m['working_metadata_SHA256']==sha(base/wr['metadata']) and m['frozen_fine_metadata_SHA256']==sha(fm) and m['frozen_fine_arrays_SHA256']==sha(fa)
  affected=[]
  for code in np.unique(q['annotation_group_code'][changed&(q['annotation_group_code']>=0)]):
   sel=q['annotation_group_code']==code;after=fine[sel].copy();after[changed[sel]]=np.asarray(TYPES,object)[proposal[sel][changed[sel]]];affected.append({'previous_group_code':int(code),'n_previous_group':int(sel.sum()),'new_fine_counts':dict(collections.Counter(after))})
  assert affected==m['affected_previous_groups']
  nc=int(changed.sum());nu=int((changed&(q['annotation_group_code']<0)).sum());assert nc==m['n_fine_changes']==entry['n_fine_changes'] and int(target.sum())==m['n_target_coarse_T'] and len(x)==m['n_scope_records'] and nu==m['n_changed_without_fitted_group']
  totals.update(n_sources=1,n_scope_records=len(x),n_target_coarse_T=int(target.sum()),n_fine_changes=nc,n_affected_previous_groups=len(affected),n_changed_without_fitted_group=nu);labels.update(cc);why.update(rrc);n_native+=len(x)*4
  save(root/'outputs'/(source+'_independent_direct_fine_replay_v1.json'),{'source':source,'n_scope_records':len(x),'own_native_fine_cues_and_scoped_unique_reference_rebuilt':True,'all_proposal_fields_exactly_replayed':True});del q,out,x;gc.collect()
 pt=load(P/'provenance/terminal_v1.json');assert dict(totals)==pt['totals'] and dict(labels)==pt['new_fine_counts'];sp=P/'outputs/direct_T_fine_proposal_summary_v1.json';assert sha(sp)==pt['summary_SHA256'];s=load(sp);assert s['target_review_counts']==dict(why)
 save(root/'provenance/terminal_v1.json',{'status':'FULL_SCOPE_DIRECT_FINE_POLICY_NATIVE_CUES_AND_REFERENCE_REPLAY_ACCEPTED','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'totals':dict(totals),'new_fine_counts':dict(labels),
  'producer_terminal_SHA256':sha(P/'provenance/terminal_v1.json'),'producer_index_SHA256':sha(ip),'own_fine_cues_rebuilt_from_independently_audited_native_quantities':True,'new_production_policy_not_imported':True,'not_orthogonal_biological_identity_accuracy':True,'wall_seconds':time.monotonic()-start})

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);r=Path(p.parse_args().output)
 try:main(r)
 except BaseException:
  import traceback
  save(r/'provenance/failure_v1.json',{'traceback':traceback.format_exc()});raise
