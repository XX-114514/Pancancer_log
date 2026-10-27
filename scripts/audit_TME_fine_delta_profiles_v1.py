"""Independent native CSR replay for every changed fine-label exact group."""
import argparse,collections,datetime,gc,hashlib,json,re,sys,time
from pathlib import Path
import numpy as np

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
 with Path(p).open('x') as f:json.dump(d,f,indent=2,allow_nan=False);f.write('\n')
def membership(source,ids):
 assert len(ids)==len(set(ids));h=hashlib.sha256()
 for v in ['atlas_exact_membership_v1',source,*sorted(ids)]:
  b=v.encode();h.update(len(b).to_bytes(8,'big'));h.update(b)
 return h.hexdigest()

def main(root):
 start=time.monotonic();S=Path(load(root/'config/spec_v1.json')['producer']);assert load(S/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 sc=load(S/'config/spec_v1.json');P=Path(sc['proposal']);pc=load(P/'config/spec_v1.json');E=Path(pc['native_evidence']);oldS=Path(sc['previous_changed_profiles']);J=Path(pc['joined']);NA=Path(pc['native_audit'])
 for stage in [P,Path(sc['proposal_audit']),E,NA,oldS,Path(sc['previous_changed_audit'])]:assert load(stage/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 sys.path.insert(0,str(J/'scripts'));sys.path.insert(0,str(NA/'scripts'))
 from native_archive_collection_v1 import NativeArchiveCollection
 from audit_full_TME_native_evidence_v1 import gather
 native=NativeArchiveCollection(load(J/'outputs/joined_native_atlas_delivery_v1.json')['archive_manifest'])
 accepted=load(pc['working_acceptance']);assert sha(accepted['producer_manifest']['path'])==accepted['producer_manifest']['SHA256'];w=load(accepted['producer_manifest']['path']);base=Path(w['index']['path']).parent
 wi={r['source']:r for r in load(w['index']['path'])};catalog=load(w['group_catalog']['path'])['groups']
 ei={r['source']:r for r in load(E/'outputs/TME_full_native_evidence_index_v1.json')};oi={r['source']:r for r in load(oldS/'outputs/TCR_changed_membership_native_index_v1.json')};pi={r['source']:r for r in load(P/'outputs/direct_T_fine_proposal_index_v1.json')}
 ip=S/'outputs/fine_delta_native_profile_index_v1.json';assert sha(ip)==load(S/'provenance/terminal_v1.json')['index_SHA256'];totals=collections.Counter();n_direct=0
 for item in load(ip):
  ds=item['source'];mp=S/'outputs'/item['metadata'];assert sha(mp)==item['metadata_SHA256'];m=load(mp);sp=m['statistics'];assert sha(sp['path'])==sp['SHA256']
  with np.load(sp['path'],allow_pickle=False) as f:actual={k:f[k] for k in f.files}
  wr=wi[ds];wm=load(base/wr['metadata']);wa=base/wr['arrays'];assert sha(wa)==wr['arrays_SHA256']
  with np.load(wa,allow_pickle=False) as f:q={k:f[k] for k in ['native_row','source_native_position','annotation_group_code','native_total_RNA','consensus_celltype_code']}
  pm=P/'outputs'/pi[ds]['metadata'];assert sha(pm)==pi[ds]['metadata_SHA256'];p=load(pm);pp=P/'outputs'/p['arrays'];assert sha(pp)==p['arrays_SHA256']
  with np.load(pp,allow_pickle=False) as f:proposal=f['proposed_fine_code'];np.testing.assert_array_equal(f['native_row'],q['native_row'])
  oldfine=np.asarray(wm['candidate_dictionaries']['consensus_celltype'],object)[q['consensus_celltype_code']];after=oldfine.copy();after[proposal>0]=np.asarray(p['fine_dictionary'],object)[proposal[proposal>0]]
  er=ei[ds];em=E/'outputs'/er['metadata'];assert sha(em)==er['metadata_SHA256'];e=load(em);genes=e['native_measured_gene_names'];cols=np.asarray(e['native_columns']);assert genes==m['native_measured_genes'];np.testing.assert_array_equal(actual['native_columns'],cols)
  available={};buffers=[]
  for meta,groups,stats,rowfield,keyfield in [(e,e['groups'],e['group_statistics'],'source_annotation_group_code','stable_group_sha256')]:
   assert sha(stats['path'])==stats['SHA256']
   with np.load(stats['path'],allow_pickle=False) as f:v={k:f[k] for k in ['n','cp_sum','log_sum','detected','native_columns']}
   np.testing.assert_array_equal(v['native_columns'],cols);buffers.append(v)
   for g in groups:available[g[keyfield]]=(v,int(g[rowfield]))
  if ds in oi:
   r=oi[ds];op=oldS/'outputs'/r['metadata'];assert sha(op)==r['metadata_SHA256'];o=load(op);assert o['native_measured_genes']==genes and sha(o['statistics']['path'])==o['statistics']['SHA256']
   with np.load(o['statistics']['path'],allow_pickle=False) as f:v={k:f[k] for k in ['n','cp_sum','log_sum','detected','native_columns']}
   np.testing.assert_array_equal(v['native_columns'],cols);buffers.append(v)
   for g in o['groups']:available[g['stable_membership_SHA256']]=(v,g['measurement_subset_code'])
  parents=sorted(set(map(int,q['annotation_group_code'][(proposal>0)&(q['annotation_group_code']>=0)])));expected_groups=[];codes=np.full(len(after),-1,np.int32);current={g['source_annotation_group_code']:g for g in catalog if g['source']==ds};rawN=0;reused=0
  for parent in parents:
   domain=q['annotation_group_code']==parent;labels=sorted(set(after[domain]));reuse=len(labels)==1
   if not reuse:rawN+=int(domain.sum())
   for label in labels:
    sel=np.flatnonzero(domain&(after==label));codes[sel]=len(expected_groups);expected_groups.append((parent,label,sel,reuse));reused+=int(reuse)
  take=np.flatnonzero(codes>=0)
  for k,v in [('scope_position',take),('measurement_subset_code',codes[take]),('native_row',q['native_row'][take])]:np.testing.assert_array_equal(actual[k],v)
  n=len(expected_groups);assert len(m['groups'])==n;expected={'n':np.zeros(n,np.int64),**{k:np.zeros((n,len(cols))) for k in ['cp_sum','log_sum','detected']}}
  csr=native.store(ds)[native.legacy_config['counts_key'] if native.is_legacy(ds) else native.primary.partitions[ds]['counts_key']]
  for lo in range(0,len(take),1024):
   use=take[lo:lo+1024];raw=gather(csr,q['native_row'][use]);rna=np.asarray(raw.sum(axis=1)).ravel();np.testing.assert_array_equal(rna,q['native_total_RNA'][use]);x=raw[:,cols].tocsr();assert (rna>0).all()
   for j in np.unique(codes[use]):
    sub=x[codes[use]==j];den=rna[codes[use]==j];cp=sub.multiply((10000/den)[:,None]).tocsr();logged=cp.copy();logged.data=np.log1p(logged.data);expected['n'][j]+=len(den)
    for k,values in [('cp_sum',cp),('log_sum',logged),('detected',sub>0)]:expected[k][j]+=np.asarray(values.sum(axis=0)).ravel()
  for k in ['n','detected']:np.testing.assert_array_equal(actual[k],expected[k])
  for k in ['cp_sum','log_sum']:np.testing.assert_allclose(actual[k],expected[k],rtol=2e-10,atol=1e-7)
  pools={'T_NK':{'n':0,'cp':np.zeros(len(cols)),'det':np.zeros(len(cols))},'all':{'n':0,'cp':np.zeros(len(cols)),'det':np.zeros(len(cols))}}
  for oldg in current.values():
   a,j=available[oldg['stable_group_sha256']];assert a['n'][j]==oldg['n_scope_group_rows']
   choices=['all']+(['T_NK'] if oldg['major_distribution'] and set(oldg['major_distribution']).issubset({'T','NK','T_NK'}) else [])
   for name in choices:
    pools[name]['n']+=int(a['n'][j]);pools[name]['cp']+=a['cp_sum'][j];pools[name]['det']+=a['detected'][j]
  ids=native.native_ids(ds)
  for j,(parent,label,sel,reuse) in enumerate(expected_groups):
   g=m['groups'][j];oldg=current[parent];a,k=available[oldg['stable_group_sha256']];children=[u for u,c in enumerate(expected_groups) if c[0]==parent]
   for name in ['n','detected']:np.testing.assert_array_equal(actual[name][children].sum(axis=0),a[name][k])
   for name in ['cp_sum','log_sum']:np.testing.assert_allclose(actual[name][children].sum(axis=0),a[name][k],rtol=2e-10,atol=1e-7)
   assert g['measurement_subset_code']==j and g['previous_group_code']==parent and g['previous_group_number']==oldg['global_subcluster_number'] and g['previous_group_SHA256']==oldg['stable_group_sha256']
   assert g['fine_celltype']==label and g['n_exact_subset']==len(sel) and g['exact_membership_reused']==reuse
   assert g['stable_membership_SHA256']==membership(ds,ids[q['source_native_position'][sel]].astype(str).tolist()) and g['native_rows_SHA256']==hashlib.sha256(np.asarray(q['native_row'][sel],'<i8').tobytes()).hexdigest()
   if reuse:
    ref=g['accepted_profile_reference'];assert ref['stable_membership_SHA256']==oldg['stable_group_sha256'] and sha(ref['statistics_path'])==ref['statistics_SHA256'] and sha(ref['metadata_path'])==ref['metadata_SHA256']
    with np.load(ref['statistics_path'],allow_pickle=False) as f:
     for name in ['n','detected']:np.testing.assert_array_equal(f[name][ref['profile_row']],actual[name][j])
     for name in ['cp_sum','log_sum']:np.testing.assert_allclose(f[name][ref['profile_row']],actual[name][j],rtol=2e-10,atol=1e-7)
   else:assert g['accepted_profile_reference'] is None
   pool=pools['T_NK'];count=pool['n']-len(sel);scope='other_current_source_T_NK_fitted_group_cells'
   if count==0:pool=pools['all'];count=pool['n']-len(sel);scope='other_current_source_fitted_group_cells_no_T_NK_comparator'
   assert count>=0;own=actual['cp_sum'][j]/len(sel);det=actual['detected'][j]/len(sel)
   if count:
    diff=pool['cp']-actual['cp_sum'][j];dd=pool['det']-actual['detected'][j];assert diff.min()>=-1e-7 and dd.min()>=0;other=np.maximum(diff,0)/count;otherdet=dd/count;effect=np.log2((own+.1)/(other+.1))
   else:other=otherdet=effect=None
   eligible=np.asarray([not re.match(r'^(MT-|RPL|RPS|ENSG|IG[HKL][VDJC])',g) for g in genes])&(det>=.1)&(own>0)
   if count:
    rank=effect*np.maximum(det-otherdet,0);order=sorted(np.flatnonzero(eligible&(effect>0)&(det>otherdet)),key=lambda k:(-rank[k],genes[k]))[:20];kind='own_subset_expression_specificity_not_inferential_DEG'
   else:order=[];kind='own_subset_expression_without_comparator_not_inferential_DEG'
   if not order:
    order=sorted(np.flatnonzero(eligible),key=lambda k:(-own[k]*det[k],genes[k]))[:20]
    if count:kind='own_subset_expression_without_positive_specificity_not_inferential_DEG'
   candidates=g['own_full_native_gene_marker_candidates_not_DEG'];assert [v['gene'] for v in candidates]==[genes[k] for k in order]
   for k,c in zip(order,candidates):
    assert c['native_mean_CP10K']==float(own[k]) and c['native_detection_fraction']==float(det[k]) and c['marker_evidence_type']==kind
    if count:
     np.testing.assert_allclose(c['descriptive_log2_mean_ratio'],float(effect[k]),rtol=1e-12,atol=1e-12);assert c['comparator_detection_fraction']==float(otherdet[k])
    else:assert c['descriptive_log2_mean_ratio'] is None and c['comparator_detection_fraction'] is None
   assert g['comparator_scope']==scope and g['n_comparator_cells']==count
  assert m['proposal_metadata_SHA256']==sha(pm) and m['working_arrays_SHA256']==sha(wa) and m['source_complete_profile_metadata_SHA256']==sha(em)
  counts=dict(n_sources=1,n_parent_groups=len(parents),n_exact_subsets=n,n_raw_records_read=rawN,n_exact_memberships_reused=reused,n_fine_changes=int((proposal>0).sum()),n_unfitted_fine_changes=int(((proposal>0)&(q['annotation_group_code']<0)).sum()))
  for k,v in counts.items():
   if k in m:assert m[k]==v
  totals.update(counts);n_direct+=len(take);save(root/'outputs'/(ds+'_independent_fine_delta_native_audit_v1.json'),dict(source=ds,n_exact_subsets=n,n_native_records_directly_reread=len(take),accepted=True));native._fields.clear();del q,actual,expected,available,buffers;gc.collect()
 assert dict(totals)==load(S/'provenance/terminal_v1.json')['totals']
 save(root/'provenance/terminal_v1.json',dict(status='ALL_FINE_DELTA_NATIVE_PROFILES_MEMBERSHIPS_REUSE_AND_MARKERS_ACCEPTED',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),totals=dict(totals),n_native_records_directly_reread=n_direct,producer_terminal_SHA256=sha(S/'provenance/terminal_v1.json'),new_production_measurement_resolver_and_ranking_not_imported=True,not_orthogonal_identity_accuracy=True,wall_seconds=time.monotonic()-start))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);r=Path(p.parse_args().output)
 try:main(r)
 except BaseException:
  import traceback
  save(r/'provenance/failure_v1.json',{'traceback':traceback.format_exc()});raise
