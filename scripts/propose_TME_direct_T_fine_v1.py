"""Full current-scope fine-cue replay; keep current labels until release audit."""
import argparse,collections,datetime,gc,hashlib,json,subprocess,sys,time
from pathlib import Path
import numpy as np
from TME_direct_T_fine_policy_v1 import TYPES,FROZEN_BITS,REASONS,direct_fine

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
 with Path(p).open('x') as f:json.dump(d,f,indent=2,allow_nan=False);f.write('\n')

def main(root):
 start=time.monotonic();cfg=load(root/'config/spec_v1.json')
 for b in cfg['bindings'].values():assert sha(b['path'])==b['SHA256']
 test=subprocess.run([sys.executable,'-B',str(root/'tests/test_TME_direct_T_fine_policy_v1.py')],capture_output=True,text=True,timeout=60);(root/'provenance/test_stderr.txt').write_text(test.stderr)
 assert test.returncode==0 and 'Ran 12 tests' in test.stderr and 'skipped' not in test.stderr.lower()
 J=Path(cfg['joined']);sys.path.insert(0,str(J/'scripts'));from extract_full_source_local_marker_evidence_v1 import storage_gate
 accepted=load(cfg['working_acceptance']);assert sha(accepted['producer_manifest']['path'])==accepted['producer_manifest']['SHA256'];m=load(accepted['producer_manifest']['path']);assert sha(m['index']['path'])==m['index']['SHA256'];base=Path(m['index']['path']).parent
 F=Path(cfg['candidate_producer']);fp=F/'outputs/full_native_cell_identity_candidate_index_v1.json';assert sha(fp)==load(F/'provenance/terminal_v1.json')['index_sha256'];fi={r['source']:r for r in load(fp)}
 records=[];counts=collections.Counter();labels=collections.Counter();why=collections.Counter()
 for item in load(m['index']['path']):
  storage_gate(root,cfg);source=item['source'];wp=base/item['metadata'];wa=base/item['arrays'];assert sha(wp)==item['metadata_SHA256'] and sha(wa)==item['arrays_SHA256'];wm=load(wp)
  with np.load(wa,allow_pickle=False) as f:q={k:f[k] for k in ['native_row','source_native_position','consensus_major_code','consensus_celltype_code','annotation_review_code','native_family_bitmask','native_total_RNA','annotation_group_code','fit_pool_candidate_not_final_QC']}
  r=fi[source];fm=F/'outputs'/r['metadata'];fa=F/'outputs'/r['arrays'];assert sha(fm)==r['metadata_sha256'] and sha(fa)==r['arrays_sha256'];assert load(fm)['fine_candidate_order'][:4]==FROZEN_BITS
  with np.load(fa,allow_pickle=False) as f:
   p=q['source_native_position'];np.testing.assert_array_equal(q['native_row'],f['native_row'][p]);own=f['own_fine_bitmask'][p];ref=f['reference_fine_bitmask'][p]
  major=np.asarray(wm['candidate_dictionaries']['consensus_major'],object)[q['consensus_major_code']];fine=np.asarray(wm['candidate_dictionaries']['consensus_celltype'],object)[q['consensus_celltype_code']]
  result=direct_fine(major,fine,own,ref,q['native_family_bitmask'],q['annotation_review_code'],q['native_total_RNA']>0);changed=result['proposed_fine_code']>0;old_new=collections.Counter(zip(fine[changed],np.asarray(TYPES,object)[result['proposed_fine_code'][changed]]));affected=[]
  for group in np.unique(q['annotation_group_code'][changed&(q['annotation_group_code']>=0)]):
   sel=q['annotation_group_code']==group;after=fine[sel].copy();after[changed[sel]]=np.asarray(TYPES,object)[result['proposed_fine_code'][sel][changed[sel]]]
   affected.append({'previous_group_code':int(group),'n_previous_group':int(sel.sum()),'new_fine_counts':dict(collections.Counter(after))})
  out=root/'outputs'/(source+'_direct_T_fine_proposals_v1.npz');arrays=dict(result,native_row=q['native_row'],source_native_position=q['source_native_position'],own_fine_bitmask=own,reference_fine_bitmask=ref)
  with out.open('xb') as f:np.savez_compressed(f,**arrays)
  with np.load(out,allow_pickle=False) as f:
   for k,v in arrays.items():np.testing.assert_array_equal(f[k],v)
  cc=dict(collections.Counter(np.asarray(TYPES,object)[result['proposed_fine_code'][changed]]));rr=dict(collections.Counter(REASONS[int(k)] for k in result['fine_review_code'][result['eligible_coarse_T']]))
  meta={'source':source,'arrays':out.name,'arrays_SHA256':sha(out),'n_scope_records':len(fine),'n_target_coarse_T':int(result['eligible_coarse_T'].sum()),'n_fine_changes':int(changed.sum()),'n_changed_without_fitted_group':int((changed&(q['annotation_group_code']<0)).sum()),
   'new_fine_counts':cc,'target_review_counts':rr,'affected_previous_groups':affected,'working_metadata_SHA256':sha(wp),'working_arrays_SHA256':sha(wa),'frozen_fine_metadata_SHA256':sha(fm),'frozen_fine_arrays_SHA256':sha(fa),
   'fine_dictionary':TYPES,'review_dictionary':REASONS,'old_to_new_fine_counts':[{'old':str(a),'new':str(b),'n':int(n)} for (a,b),n in sorted(old_new.items())],'current_labels_and_states_not_changed':True,'not_independent_identity_truth':True}
  mp=root/'outputs'/(source+'_direct_T_fine_proposals_v1.json');save(mp,meta);records.append({'source':source,'metadata':mp.name,'metadata_SHA256':sha(mp),'n_fine_changes':int(changed.sum())})
  counts.update(n_sources=1,n_scope_records=len(fine),n_target_coarse_T=meta['n_target_coarse_T'],n_fine_changes=meta['n_fine_changes'],n_affected_previous_groups=len(affected),n_changed_without_fitted_group=meta['n_changed_without_fitted_group']);labels.update(cc);why.update(rr);del q,arrays,result;gc.collect()
 ip=root/'outputs/direct_T_fine_proposal_index_v1.json';save(ip,records);sp=root/'outputs/direct_T_fine_proposal_summary_v1.json';save(sp,{'totals':dict(counts),'new_fine_counts':dict(labels),'target_review_counts':dict(why),'labels_not_applied':True})
 save(root/'provenance/terminal_v1.json',{'status':'FULL_CURRENT_SCOPE_DIRECT_T_FINE_CUES_REPLAYED_INDEPENDENT_AUDIT_REQUIRED','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'totals':dict(counts),'new_fine_counts':dict(labels),'index_SHA256':sha(ip),'summary_SHA256':sha(sp),'tests_actual':12,'tests_skipped':0,'current_labels_QC_states_and_CM_unchanged':True,'wall_seconds':time.monotonic()-start})

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);r=Path(p.parse_args().output)
 try:main(r)
 except BaseException:
  import traceback
  save(r/'provenance/failure_v1.json',{'traceback':traceback.format_exc()});raise
