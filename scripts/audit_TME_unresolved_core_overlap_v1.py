"""Full independent CSR replay of native single- and dual-core counts."""
import argparse,collections,datetime,hashlib,json,sys,time
from pathlib import Path
import numpy as np
from TME_all_lineage_cluster_policy_v1 import GATES,BROAD
from measure_TME_unresolved_core_overlap_v1 import pair_pattern

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')

def main(root):
 cfg=load(root/'config/spec_v1.json');R=Path(cfg['producer']);EA=Path(cfg['independent_native_auditor'])
 while not (R/'provenance/supervisor_terminal_v1.json').exists():time.sleep(3)
 assert load(R/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0;rt=load(R/'provenance/terminal_v1.json');rc=load(R/'config/spec_v1.json');P=Path(rc['producer']);pc=load(P/'config/spec_v1.json');assert rt['producer_terminal_SHA256']==sha(P/'provenance/terminal_v1.json')
 for name,digest in rt['files'].items():assert sha(R/'outputs'/name)==digest
 index=load(R/'outputs/unresolved_core_overlap_index_v1.json');reviews=load(R/'outputs/unresolved_cluster_core_overlap_review_v1.json');summary=load(R/'outputs/unresolved_core_overlap_summary_v1.json');targets=[r for r in load(P/'outputs/cluster_hierarchy_named_catalog_v1.json')['nodes'] if r['final_leaf'] and r['cluster_RNA_proposal']['major']=='unresolved'];review_by={r['stable_membership_SHA256']:r for r in reviews};assert len(review_by)==len(reviews)==len(targets)
 source_index={r['source']:r for r in index['sources']};assert set(source_index)=={r['source'] for r in targets};pairs=[(a,b) for i,a in enumerate(BROAD) for b in BROAD[i+1:]];assert index['broad_core_order']==BROAD and index['pair_order']==[list(v) for v in pairs]
 accepted=load(Path(pc['working_delivery'])/'outputs/accepted_TME_scoped_working_annotation_v1.json');w=load(accepted['producer_manifest']['path']);wb=Path(w['index']['path']).parent;wi={r['source']:r for r in load(w['index']['path'])};J=Path(rc['joined']);jm=load(J/'outputs/joined_native_atlas_delivery_v1.json');sys.path.insert(0,str(J/'scripts'));sys.path.insert(0,str(EA/'scripts'))
 from native_archive_collection_v1 import NativeArchiveCollection
 from audit_full_TME_native_evidence_v1 import gather
 native=NativeArchiveCollection(jm['archive_manifest']);totals=collections.Counter();pair_totals=collections.defaultdict(collections.Counter)
 for entry in load(P/'outputs/cluster_hierarchy_index_v1.json'):
  source=entry['source'];selected=[r for r in targets if r['source']==source]
  if not selected:continue
  m=load(P/'outputs'/entry['metadata']);genes=m['native_measured_genes'];by={g:j for j,g in enumerate(genes)};assert sha(m['arrays']['path'])==m['arrays']['SHA256'] and sha(m['statistics']['path'])==m['statistics']['SHA256']
  with np.load(m['arrays']['path'],allow_pickle=False) as f:pool=f['source_scope_position'];leaf=f['final_leaf_node'];rows=f['native_row']
  ids=np.asarray([r['source_node_code'] for r in selected]);lookup={int(n):j for j,n in enumerate(ids)};take=np.asarray([j for j,n in enumerate(leaf) if int(n) in lookup],np.int64);codes=np.asarray([lookup[int(n)] for n in leaf[take]],np.int32)
  wr=wi[source];assert sha(wb/wr['arrays'])==wr['arrays_SHA256']
  with np.load(wb/wr['arrays'],allow_pickle=False) as f:np.testing.assert_array_equal(rows[take],f['native_row'][pool[take]]);rna=f['native_total_RNA'][pool[take]]
  with np.load(m['statistics']['path'],allow_pickle=False) as f:cols=f['native_columns']
  n=np.zeros(len(ids),np.int64);singles=np.zeros((len(ids),len(BROAD)),np.int64);joint=np.zeros((len(ids),len(pairs)),np.int64);csr=native.store(source)[native.legacy_config['counts_key'] if native.is_legacy(source) else native.primary.partitions[source]['counts_key']]
  for lo in range(0,len(take),1024):
   pos=np.arange(lo,min(lo+1024,len(take)));raw=gather(csr,rows[take[pos]]).astype(np.float64).tocsr();raw.sum_duplicates();raw.eliminate_zeros();assert np.isfinite(raw.data).all() and (raw.data>=0).all() and np.equal(raw.data,np.floor(raw.data)).all();np.testing.assert_array_equal(np.asarray(raw.sum(axis=1)).ravel(),rna[pos]);x=raw[:,cols].tocsr();flags=[]
   def any_of(names):
    c=[by[g] for g in names if g in by];return np.asarray((x[:,c]>0).sum(axis=1)).ravel()>0 if c else np.zeros(len(pos),bool)
   for name in BROAD:
    flag=np.ones(len(pos),bool)
    for module in GATES[name]:flag&=any_of(module)
    if name=='T' and not any(g in by for g in ['TRAC','TRBC1','TRBC2','TRDC']) and all(g in by for g in ['CD3D','CD3E']) and any(g in by for g in ['CD2','CD7','TRAT1']) and sum(g in by for g in ['CD247','LCK','LAT','ZAP70','BCL11B'])>=2:
     sig=np.zeros(len(pos),np.int16)
     for g in ['CD247','LCK','LAT','ZAP70','BCL11B']:sig+=any_of([g])
     flag=any_of(['CD3D'])&any_of(['CD3E'])&any_of(['CD2','CD7','TRAT1'])&(sig>=2)
    flags.append(flag)
   f=np.column_stack(flags);pair_flags=np.column_stack([f[:,BROAD.index(a)]&f[:,BROAD.index(b)] for a,b in pairs]);np.add.at(n,codes[pos],1);np.add.at(singles,codes[pos],f.astype(np.int64));np.add.at(joint,codes[pos],pair_flags.astype(np.int64))
  output=source_index[source];assert sha(output['path'])==output['SHA256'] and output['n_records_read']==len(take)
  with np.load(output['path'],allow_pickle=False) as f:np.testing.assert_array_equal(f['source_leaf_node_code'],ids);np.testing.assert_array_equal(f['n'],n);np.testing.assert_array_equal(f['single_core_counts'],singles);np.testing.assert_array_equal(f['pair_core_counts'],joint)
  for j,r in enumerate(selected):
   report=review_by[r['stable_membership_SHA256']];active=r['cluster_RNA_proposal']['coherent_supported_families'];pp=[];assert n[j]==r['n_records']
   for k,(a,b) in enumerate(pairs):
    if a not in active or b not in active:continue
    aa=int(singles[j,BROAD.index(a)]);bb=int(singles[j,BROAD.index(b)]);ab=int(joint[j,k]);kind=pair_pattern(int(n[j]),aa,bb,ab,rc['minimum_core_population'],rc['separated_maximum_overlap'],rc['coexpressed_minimum_overlap']);pp.append(dict(core_a=a,core_b=b,n_core_a=aa,n_core_b=bb,n_both=ab,pattern=kind));pair_totals[a+'__'+b].update(n_groups=1,n_records_in_groups=int(n[j]),n_both_core_positive_records=ab)
   assert report['pair_patterns']==pp and report['active_core_families']==active;separate=any(p['pattern']=='separated_core_populations' for p in pp);coexp=any(p['pattern']=='same_cell_core_coexpression' for p in pp);kind='mixed_separation_and_coexpression' if separate and coexp else 'separated_core_populations' if separate else 'same_cell_core_coexpression' if coexp else 'intermediate_overlap' if pp else 'no_multiple_active_cores';deeper=separate and n[j]>=rc['minimum_deeper_parent_records'];assert report['group_pattern']==kind and report['further_clustering_candidate']==deeper;totals.update(n_groups=1,n_records=int(n[j]),n_deeper_clustering_candidate_groups=int(deeper),n_deeper_clustering_candidate_records=int(n[j]) if deeper else 0);totals[kind+'_groups']+=1;totals[kind+'_records']+=int(n[j])
  native._fields.clear()
 assert dict(totals)==rt['totals']==summary['totals'];assert {k:dict(v) for k,v in pair_totals.items()}==summary['pair_counts'];assert totals['n_records']==sum(r['n_records'] for r in targets);save(root/'provenance/terminal_v1.json',dict(status='ALL_UNRESOLVED_NATIVE_CORE_INTERSECTION_COUNTS_AND_RULE_REPLAY_ACCEPTED_NOT_BIOLOGICAL_DIAGNOSIS',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),producer_terminal_SHA256=sha(R/'provenance/terminal_v1.json'),totals=dict(totals),independent_native_CSR_and_boolean_counts=True,current_working_labels_and_clusters_unchanged=True,new_CM_started=False))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
 try:main(root)
 except BaseException:
  import traceback
  save(root/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
