"""Group-level cross-version, marker-role and measurement-coverage diagnostics.

This stage records candidate interpretation flags; it never diagnoses ambient
RNA, excludes cells, or silently changes accepted cluster identities/markers.
All output stays private unless a separate aggregate export is authorized.
"""
from pathlib import Path
import argparse,collections,datetime,hashlib,json
import numpy as np

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main(root):
 spec=load(root/'config/spec_v1.json');R=Path(spec['producer']);A=Path(spec['auditor']);D=Path(spec['delivery'])
 for r in [R,A,D]:assert load(r/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 assert sha(R/'outputs/cluster_hierarchy_named_catalog_v1.json')==load(R/'provenance/terminal_v1.json')['catalog_SHA256']
 cfg=load(R/'config/spec_v1.json');baseline=load(Path(cfg['working_delivery'])/'outputs/accepted_TME_scoped_working_annotation_v1.json');w=load(baseline['producer_manifest']['path']);base=Path(w['index']['path']).parent;wi={r['source']:r for r in load(w['index']['path'])};nodes=load(R/'outputs/cluster_hierarchy_named_catalog_v1.json')['nodes'];by={(r['source'],r['source_node_code']):r for r in nodes};oldgroups={(r['source'],r['source_annotation_group_code']):r for r in load(w['group_catalog']['path'])['groups']};aliases={'Tcell':'T','NKcell':'NK','gamma_delta_T_RNA_candidate':'T_gamma_delta_RNA_candidate'};roles={'erythroid_associated_leading_marker':set('HBA1 HBA2 HBB HBD HBG1 HBG2 HBM HBZ AHSP ALAS2 GYPA SLC4A1'.split()),'platelet_associated_leading_marker':set('PF4 PPBP GP9 GP1BA GP1BB ITGA2B'.split()),'heatshock_stress_or_processing_marker':set('HSPA1A HSPA1B DNAJB1 HSPH1 HSP90AA1 HSPA6'.split()),'immediate_early_program_marker':set('FOS FOSB JUN JUNB JUND DUSP1 EGR1 ATF3'.split()),'broad_RNA_content_marker':set('MALAT1 NEAT1 ACTB GAPDH TMSB4X TMSB10 EEF1A1'.split())};cross=[];membership_cross=[];flags=[];native_gene_sets={};sources_by_type=collections.defaultdict(set);support=[]
 assert sha(baseline['producer_manifest']['path'])==baseline['producer_manifest']['SHA256']
 assert sha(w['index']['path'])==w['index']['SHA256'] and sha(w['group_catalog']['path'])==w['group_catalog']['SHA256']
 for r in load(R/'outputs/cluster_hierarchy_index_v1.json'):
  source=r['source'];mp=R/'outputs'/r['metadata'];assert sha(mp)==r['metadata_SHA256'];m=load(mp);native_gene_sets[source]=set(m['native_measured_genes']);assert sha(m['arrays']['path'])==m['arrays']['SHA256'];wm=load(base/wi[source]['metadata'])
  with np.load(m['arrays']['path'],allow_pickle=False) as f:pool=f['source_scope_position'];leaf=f['final_leaf_node']
  with np.load(base/wi[source]['arrays'],allow_pickle=False) as f:old=np.asarray(wm['candidate_dictionaries']['consensus_celltype'],object)[f['consensus_celltype_code'][pool]];oldcode=f['annotation_group_code'][pool]
  new=np.asarray([by[source,int(n)]['cluster_RNA_proposal']['celltype'] for n in leaf],object)
  for (a,b),n in sorted(collections.Counter(zip(map(str,old),map(str,new))).items()):cross.append(dict(source=source,baseline_v5_celltype=a,cluster_RNA_proposal=b,canonical_cluster_proposal=aliases.get(b,b),n_records=n,counts_are_cluster_membership_proposals_not_biological_accuracy=True))
  pair=collections.Counter(zip(map(int,oldcode),map(int,leaf)));old_to_new=collections.defaultdict(set);new_to_old=collections.defaultdict(set);oldpool=collections.Counter(map(int,oldcode))
  for oc,nc in pair:old_to_new[oc].add(nc);new_to_old[nc].add(oc)
  for (oc,nc),n in sorted(pair.items()):
   og=oldgroups[source,oc];ng=by[source,nc];split=len(old_to_new[oc])>1;merged=len(new_to_old[nc])>1;outside=og['n_scope_group_rows']-oldpool[oc];assert outside>=0
   if not split and not merged and outside==0:assert og['stable_group_sha256']==ng['stable_membership_SHA256'] and og['global_subcluster_number']==ng['global_display_number']
   relation='many_to_many' if split and merged else 'split' if split else 'merge' if merged else ('partial_previous_group' if outside else 'same_membership_view')
   membership_cross.append(dict(source=source,old_v5_group_SHA256=og['stable_group_sha256'],old_v5_display_number=og['global_subcluster_number'],old_v5_display=og['subcluster'],new_cluster_group_SHA256=ng['stable_membership_SHA256'],new_cluster_display_number=ng['global_display_number'],new_cluster_display=ng['display_label'],n_records_in_pair=n,n_old_group_in_first_round_pool=oldpool[oc],n_old_group_outside_first_round_pool=outside,relation_within_first_round_pool=relation))
  ap=A/'outputs'/(source+'_cluster_hierarchy_audit_v1.json');aud=load(ap);assert aud['producer_metadata_SHA256']==sha(mp);dom={r['source_node_code']:r for r in aud['anonymous_unit_dominance']}
  for node in np.unique(leaf):
   g=by[source,int(node)];call=g['cluster_RNA_proposal'];sources_by_type[call['celltype']].add(source)
   rr=[role for role,genes in roles.items() if (g['marker'] or '').upper() in genes]
   if rr:flags.append(dict(source=source,subcluster=g['display_label'],celltype_proposal=call['celltype'],marker=g['marker'],n_records=g['n_records'],review_roles=rr,own_ranked_marker_candidate_names=[x['gene'] for x in g['marker_candidates']],does_not_prove_contamination_or_dissociation=True,no_automatic_cell_exclusion_or_relabel=True))
   support.append(dict(source=source,subcluster=g['display_label'],celltype_proposal=call['celltype'],n_records=g['n_records'],n_declared_units_not_verified_patients=g['n_declared_expression_units_not_verified_patients_or_libraries'],**{k:v for k,v in dom[int(node)].items() if k!='source_node_code'}))
 coverage=[]
 for typ,src in sorted(sources_by_type.items()):
  common=set.intersection(*(native_gene_sets[s] for s in sorted(src)));coverage.append(dict(cluster_proposal=typ,n_sources=len(src),n_common_unique_native_measured_genes=len(common),not_a_shared_taxonomy_model_fit=True))
 allcommon=set.intersection(*native_gene_sets.values());summary={k:dict(n_clusters=sum(k in x['review_roles'] for x in flags),n_records=sum(x['n_records'] for x in flags if k in x['review_roles'])) for k in roles}
 assert sum(r['n_records_in_pair'] for r in membership_cross)==load(R/'provenance/terminal_v1.json')['totals']['n_pool_records']
 out=dict(observed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),status='AGGREGATE_MARKER_ROLE_CROSSWALK_AND_SHARED_FEATURE_COVERAGE_REVIEW_NOT_NEW_ANNOTATION',crosswalk=cross,version_membership_crosswalk=membership_cross,canonical_token_aliases=aliases,marker_role_review=flags,anonymous_unit_support=support,shared_native_gene_coverage_by_cluster_proposal=coverage,n_common_unique_native_genes_all_sources=len(allcommon),marker_role_counts=summary,no_individual_identifiers_expression_values_or_parameters=True,current_v5_cluster_proposals_and_markers_unchanged=True);op=root/'outputs/TME_cluster_crosswalk_marker_role_and_coverage_review_v1.json'
 with op.open('x') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
 with (root/'provenance/terminal_v1.json').open('x') as f:json.dump(dict(status='CLUSTER_ROLE_CROSSWALK_AND_COVERAGE_REVIEW_READY_NO_RELABELING',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),output_SHA256=sha(op),n_first_round_records=sum(r['n_records_in_pair'] for r in membership_cross),current_working_labels_and_cluster_assignments_unchanged=True),f,indent=2);f.write('\n')

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
 try:main(root)
 except BaseException:
  import traceback
  with (root/'provenance/failure_v1.json').open('x') as f:json.dump(dict(traceback=traceback.format_exc()),f,indent=2)
  raise
