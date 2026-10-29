"""Distinguish separate core-positive populations from same-cell coexpression.

Read all unresolved fitted cluster members in the accepted full tree. This is
a descriptive diagnostic, never a doublet/progenitor/EMT diagnosis or relabel.
"""
import argparse,collections,datetime,hashlib,json,subprocess,sys,time
from pathlib import Path
import numpy as np
from TME_all_lineage_cluster_policy_v1 import BROAD,GATE_ORDER,coherent_flags

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')

def pair_pattern(n,a,b,both,minimum_population=30,separated_maximum_overlap=.1,coexpressed_minimum_overlap=.5):
 if min(n,a,b,both)<0 or max(a,b)>n or not max(0,a+b-n)<=both<=min(a,b):raise ValueError('Inconsistent pair counts')
 if min(a,b)==0:return 'no_dual_positive_core_support'
 ratio=both/min(a,b)
 if ratio<=separated_maximum_overlap:return 'separated_core_populations' if min(a-both,b-both)>=minimum_population else 'separated_small_populations'
 if ratio>=coexpressed_minimum_overlap:return 'same_cell_core_coexpression'
 return 'intermediate_core_overlap'

def main(root):
 from TME_native_measurement_core_v1 import canonical_native
 start=time.monotonic();cfg=load(root/'config/spec_v1.json');P=Path(cfg['producer']);A=Path(cfg['auditor']);D=Path(cfg['delivery']);J=Path(cfg['joined']);pc=load(P/'config/spec_v1.json')
 for p in [P,A,D]:assert load(p/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 assert load(A/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(P/'provenance/terminal_v1.json');assert sha(P/'outputs/cluster_hierarchy_named_catalog_v1.json')==load(P/'provenance/terminal_v1.json')['catalog_SHA256']
 test=subprocess.run([sys.executable,'-B',str(root/'tests/test_TME_unresolved_core_overlap_v1.py')],capture_output=True,text=True,timeout=60);(root/'provenance/test_stderr.txt').write_text(test.stderr);assert test.returncode==0 and 'Ran 4 tests' in test.stderr and 'skipped' not in test.stderr
 catalog=load(P/'outputs/cluster_hierarchy_named_catalog_v1.json')['nodes'];targets=[r for r in catalog if r['final_leaf'] and r['cluster_RNA_proposal']['major']=='unresolved'];target_by={(r['source'],r['source_node_code']):r for r in targets};baseline=load(Path(pc['working_delivery'])/'outputs/accepted_TME_scoped_working_annotation_v1.json');w=load(baseline['producer_manifest']['path']);wb=Path(w['index']['path']).parent;wi={r['source']:r for r in load(w['index']['path'])};jm=load(J/'outputs/joined_native_atlas_delivery_v1.json');sys.path.insert(0,str(J/'scripts'))
 from native_archive_collection_v1 import NativeArchiveCollection
 from native_panel_aggregate_v1 import gather_csr_rows
 native=NativeArchiveCollection(jm['archive_manifest']);pairs=[(a,b) for j,a in enumerate(BROAD) for b in BROAD[j+1:]];single_columns=[GATE_ORDER.index(k) for k in BROAD];pair_columns=[(BROAD.index(a),BROAD.index(b)) for a,b in pairs];outputs=[];reviews=[];summary=collections.Counter();pair_summary=collections.defaultdict(collections.Counter)
 for e in load(P/'outputs/cluster_hierarchy_index_v1.json'):
  source=e['source'];selected=[r for r in targets if r['source']==source]
  if not selected:continue
  mp=P/'outputs'/e['metadata'];assert sha(mp)==e['metadata_SHA256'];m=load(mp);assert sha(m['arrays']['path'])==m['arrays']['SHA256'] and sha(m['statistics']['path'])==m['statistics']['SHA256']
  with np.load(m['arrays']['path'],allow_pickle=False) as f:pool=f['source_scope_position'];rows=f['native_row'];leaf=f['final_leaf_node']
  ids=np.array([r['source_node_code'] for r in selected],np.int64);take=np.flatnonzero(np.isin(leaf,ids));lookup={int(n):j for j,n in enumerate(ids)};codes=np.asarray([lookup[int(n)] for n in leaf[take]],np.int32)
  with np.load(m['statistics']['path'],allow_pickle=False) as f:cols=f['native_columns']
  wr=wi[source];assert sha(wb/wr['arrays'])==wr['arrays_SHA256']
  with np.load(wb/wr['arrays'],allow_pickle=False) as f:rna=f['native_total_RNA'][pool[take]];np.testing.assert_array_equal(rows[take],f['native_row'][pool[take]])
  counts=np.zeros(len(ids),np.int64);singles=np.zeros((len(ids),len(BROAD)),np.int64);both=np.zeros((len(ids),len(pairs)),np.int64);csr=native.store(source)[native.legacy_config['counts_key'] if native.is_legacy(source) else native.primary.partitions[source]['counts_key']]
  for lo in range(0,len(take),cfg['read_block_size']):
   p=np.arange(lo,min(lo+cfg['read_block_size'],len(take)));raw=canonical_native(gather_csr_rows(csr,rows[take[p]],max_span=8192));np.testing.assert_array_equal(np.asarray(raw.sum(axis=1)).ravel(),rna[p]);flags,assessed=coherent_flags(raw[:,cols].tocsr(),m['native_measured_genes']);v=flags[:,single_columns];cc=codes[p];counts+=np.bincount(cc,minlength=len(ids))
   for j in range(len(BROAD)):singles[:,j]+=np.bincount(cc,weights=v[:,j],minlength=len(ids)).astype(np.int64)
   for j,(a,b) in enumerate(pair_columns):both[:,j]+=np.bincount(cc,weights=v[:,a]&v[:,b],minlength=len(ids)).astype(np.int64)
  audit=load(A/'outputs'/(source+'_cluster_hierarchy_audit_v1.json'));dominance={r['source_node_code']:r for r in audit['anonymous_unit_dominance']}
  for j,node in enumerate(ids):
   r=target_by[source,int(node)];n=r['n_records'];assert counts[j]==n;call=r['cluster_RNA_proposal'];active=call['coherent_supported_families'];pp=[]
   for k,family in enumerate(BROAD):
    fraction=call['native_gate_fractions'][family]
    if fraction is None:assert singles[j,k]==0
    else:assert int(round(fraction*n))==singles[j,k]
   for k,(a,b) in enumerate(pairs):
    if a not in active or b not in active:continue
    aa=int(singles[j,BROAD.index(a)]);bb=int(singles[j,BROAD.index(b)]);ab=int(both[j,k]);pattern=pair_pattern(n,aa,bb,ab,cfg['minimum_core_population'],cfg['separated_maximum_overlap'],cfg['coexpressed_minimum_overlap']);pp.append(dict(core_a=a,core_b=b,n_core_a=aa,n_core_b=bb,n_both=ab,pattern=pattern));pair_summary[a+'__'+b].update(n_groups=1,n_records_in_groups=n,n_both_core_positive_records=ab)
   separate=any(p['pattern']=='separated_core_populations' for p in pp);coexpression=any(p['pattern']=='same_cell_core_coexpression' for p in pp);kind='mixed_separation_and_coexpression' if separate and coexpression else 'separated_core_populations' if separate else 'same_cell_core_coexpression' if coexpression else 'intermediate_overlap' if pp else 'no_multiple_active_cores';deeper=separate and n>=cfg['minimum_deeper_parent_records'];dom=dominance[int(node)];reviews.append(dict(source=source,subcluster=r['display_label'],stable_membership_SHA256=r['stable_membership_SHA256'],source_node_code=int(node),n_records=n,active_core_families=active,pair_patterns=pp,group_pattern=kind,further_clustering_candidate=deeper,n_anonymous_units_at_least_20_cells=dom['n_declared_units_with_at_least_20_cells'],not_independent_patient_support=True,not_doublet_or_EMT_or_progenitor_diagnosis=True));summary.update(n_groups=1,n_records=n,n_deeper_clustering_candidate_groups=int(deeper),n_deeper_clustering_candidate_records=n if deeper else 0);summary[kind+'_groups']+=1;summary[kind+'_records']+=n
  op=root/'outputs'/(source+'_unresolved_core_overlap_v1.npz')
  with op.open('xb') as f:np.savez_compressed(f,source_leaf_node_code=ids,n=counts,single_core_counts=singles,pair_core_counts=both)
  outputs.append(dict(source=source,path=str(op),SHA256=sha(op),n_records_read=len(take)));native._fields.clear()
 assert summary['n_records']==sum(r['n_records'] for r in targets) and summary['n_groups']==len(targets);save(root/'outputs/unresolved_core_overlap_index_v1.json',dict(broad_core_order=BROAD,pair_order=pairs,sources=outputs));save(root/'outputs/unresolved_cluster_core_overlap_review_v1.json',reviews);save(root/'outputs/unresolved_core_overlap_summary_v1.json',dict(totals=dict(summary),pair_counts={k:dict(v) for k,v in pair_summary.items()},pair_record_counts_can_overlap_and_are_not_summed_as_unique_cells=True,all_current_unresolved_fit_records_measured=True,current_identities_and_clusters_unchanged=True,proposed_deeper_clustering_not_executed=True))
 save(root/'provenance/terminal_v1.json',dict(status='ALL_UNRESOLVED_CLUSTER_NATIVE_CORE_OVERLAP_MEASURED_REQUIRES_REVIEW',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),totals=dict(summary),producer_terminal_SHA256=sha(P/'provenance/terminal_v1.json'),files={p.name:sha(p) for p in (root/'outputs').iterdir()},no_label_QC_or_cluster_changes=True,new_CM_started=False,wall_seconds=time.monotonic()-start))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
 try:main(root)
 except BaseException:
  import traceback
  save(root/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
