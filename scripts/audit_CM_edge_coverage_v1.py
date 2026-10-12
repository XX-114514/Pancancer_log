"""Independently replay original CM coverage via row-wise counters, no producer import."""
from pathlib import Path
from datetime import datetime,timezone
from collections import Counter,defaultdict
import os,json,csv,hashlib

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return list(csv.DictReader(Path(p).open(),delimiter='\t'))
def main():
 start=datetime.now(timezone.utc).isoformat();R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);B=AR/'20261010_expression_annotation_downstream_v1';P=R/'immune_same_name_family_source_review_v1';O=P/'CM_edge_coverage_v1';assert not (O/'independent_audit.json').exists();plan=json.loads(Path(os.environ['PANCANCER_CM_COVERAGE_PLAN']).read_text());assert sha(__file__)==plan['code_sha256'][Path(__file__).name];inv=json.loads((O/'inventory.json').read_text());assert inv['plan_sha256']==sha(os.environ['PANCANCER_CM_COVERAGE_PLAN']);assert inv['input_sha256']==plan['input_sha256']
 for n,h in inv['outputs'].items():assert sha(O/n)==h
 paths=dict(edges=B/'downstream_review_v2/corrected_CM_edges_source_and_marker_evidence.tsv',counts=B/'inputs/sample_major_subcluster_counts.tsv',sample_metadata=B/'inputs/sample_metadata.tsv',W=B/'covarnet/tables/nmf_W_subcluster_by_cm.tsv',H=B/'covarnet/tables/nmf_H_cm_by_sample.tsv',F=B/'covarnet/tables/frequency_normalized.tsv',catalogue=B/'downstream_review_v2/number_to_annotation_dictionary.tsv',liana_status=B/'liana/sample_status.tsv',liana_audit=B/'liana/independent_audit.json',provenance_inventory=P/'patient_block_provenance_overlay_v1/inventory.json',provenance_audit=P/'patient_block_provenance_overlay_v1/independent_audit.json',provenance_scope=P/'patient_block_provenance_overlay_v1/private_all1290_patient_provenance_overlay.tsv')
 for k,p in paths.items():assert sha(p)==plan['input_sha256'][k]['sha256'],k
 u={r['sampleID']:r for r in read(P/'patient_block_provenance_overlay_v1/private_all1290_patient_provenance_overlay.tsv')};assert len(u)==1290;count=Counter();major=Counter();total=Counter()
 for r in read(B/'inputs/sample_major_subcluster_counts.tsv'):
  s=r['sampleID'];n=int(r['n_cells']);assert n>=0
  if s in u:count[(s,r['subCluster'])]+=n;major[(s,r['majorCluster'])]+=n;total[s]+=n
 meta={r['sampleID']:r for r in read(B/'inputs/sample_metadata.tsv')};assert all(total[s]==int(meta[s]['n_cells']) for s in u);cat={r['subCluster']:r['majorCluster'] for r in read(B/'downstream_review_v2/number_to_annotation_dictionary.tsv')};lstat={r['sampleID']:r['status'] for r in read(B/'liana/sample_status.tsv')};groups=defaultdict(list)
 for s,r in u.items():
  assert r['effective_resampling_block']==r['analysis_block'];assert r['analysis_eligible'] in ['True','False']
  coarse='||'.join([r['gse_id'],r['cancerType'],r['sample_class']]);groups[('reviewed_context',r['stratum'])].append(s);groups[('coarse_GSE_cancer_class_diagnostic',coarse)].append(s)
 edges=read(B/'downstream_review_v2/corrected_CM_edges_source_and_marker_evidence.tsv');pairs=sorted({tuple(sorted([r['subCluster1'],r['subCluster2']])) for r in edges});assert len(edges)==120 and len(pairs)==60
 expectedmap=set()
 for r in edges:
  a,b=sorted([r['subCluster1'],r['subCluster2']]);pid='PAIR'+str(pairs.index((a,b))+1).zfill(3)
  expectedmap.update([(r['cm'],pid,a,b),(r['cm'],pid,b,a)])
 mapping=read(O/'all240_directed_CM_edge_mapping.tsv');assert len(mapping)==240 and {(r['cm'],r['pair_id'],r['source'],r['target']) for r in mapping}==expectedmap
 rows=read(O/'all60_pair_context_coverage.tsv');assert len(rows)==60*len(groups);seen=set();checked=0;repeated=0;failed20=0
 for r in rows:
  mode,st=r['scope_mode'],r['stratum'];pid=r['pair_id'];i=int(pid[4:])-1;a,b=pairs[i];assert r['source']==a and r['target']==b and r['source_major']==cat[a] and r['target_major']==cat[b];key=(mode,st,pid);assert key not in seen;seen.add(key);ids=groups[(mode,st)];base=u[ids[0]]
  for label in ['gse_id','cancerType','sample_class']:assert r[label]==base[label] and all(u[s][label]==base[label] for s in ids)
  eligible=[s for s in ids if u[s]['analysis_eligible']=='True'];bm=[s for s in eligible if major[s,cat[a]]>=20 and major[s,cat[b]]>=20];bn=[s for s in eligible if count[s,a]>=20 and count[s,b]>=20];paired=[s for s in bn if lstat[s]=='success'];review=[s for s in paired if u[s]['unit_kind']=='reviewed_patient'];blocks={u[s]['effective_resampling_block'] for s in paired};n=len(paired);nb=len(blocks)
  expected=dict(n_scope_libraries=len(ids),n_context_eligible=len(eligible),n_both_major_ge20=len(bm),n_both_nodes_ge20=len(bn),n_endpoint_eligible_LIANA_success=n,n_analysis_blocks=nb,n_reviewed_patient_libraries=len(review),n_reviewed_source_patient_blocks=len({u[s]['effective_resampling_block'] for s in review}),n_sample_fallback_libraries=n-len(review),n_repeated_library_rows_above_block_count=n-nb,n_endpoint_eligible_LIANA_not_success=len(bn)-n,n_source_major_not_captured=sum(major[s,cat[a]]==0 for s in eligible),n_target_major_not_captured=sum(major[s,cat[b]]==0 for s in eligible),n_source_node_not_captured=sum(count[s,a]==0 for s in eligible),n_target_node_not_captured=sum(count[s,b]==0 for s in eligible),n_source_cells_on_paired=sum(count[s,a] for s in paired),n_target_cells_on_paired=sum(count[s,b] for s in paired))
  for k,v in expected.items():assert int(r[k])==v,(key,k)
  vari=(len({count[s,a] for s in paired})>1 and len({count[s,b] for s in paired})>1) if n else False
  for k,v in dict(node_counts_both_variable=vari,coverage_ge10libraries_ge10blocks=n>=10 and nb>=10,coverage_ge20libraries_ge10blocks=n>=20 and nb>=10).items():assert r[k]==str(v),(key,k)
  checked+=len(expected)+3;repeated+=int(n>nb);failed20+=int(n<20 or nb<10)
 # Independently check every public aggregate without treating duplicated CM maps as evidence.
 for mode in ['reviewed_context','coarse_GSE_cancer_class_diagnostic']:
  g=[r for r in rows if r['scope_mode']==mode];ge={k:[r for r in g if r['coverage_ge'+str(k)+'libraries_ge10blocks']=='True'] for k in [10,20]};summary=dict(contexts=len({r['stratum'] for r in g}),pair_context_records=len(g))
  for k,sel in ge.items():
   summary.update({f'pair_contexts_ge{k}':len(sel),f'sources_ge{k}':len({r['gse_id'] for r in sel}),f'pairs_ge{k}':len({r['pair_id'] for r in sel})});ps=defaultdict(set);pc=defaultdict(set)
   for r in sel:ps[r['pair_id']].add(r['gse_id']);pc[r['pair_id']].add(r['cancerType'])
   summary[f'pairs_ge{k}_in_atleast2GSE']=sum(len(v)>=2 for v in ps.values());summary[f'pairs_ge{k}_in_atleast2cancer_labels']=sum(len(v)>=2 for v in pc.values())
  assert summary==inv['summary'][mode]
 out=dict(status='PASS_ALL60_PAIR_CONTEXT_COUNT_BLOCK_STATUS_AND240_MAPPING_ROW_REPLAY',started_UTC=start,finished_UTC=datetime.now(timezone.utc).isoformat(),auditor_sha256=sha(__file__),producer_inventory_sha256=sha(O/'inventory.json'),coverage_rows_replayed=len(rows),mapping_rows_replayed=240,integer_and_boolean_fields_checked=checked,rows_with_repeated_library_blocks=repeated,rows_below20_library10_block_gate=failed20,method='Independent CSV row counters and explicit sample sets; no producer import/pivot or sample identities exported',new_LR_tests=False,independent_biological_validation=False)
 (O/'independent_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':main()
