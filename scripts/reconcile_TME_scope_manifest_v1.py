"""Reconcile a historical unit list with an audited, cell-aware source gate.

The generated registry is preparation-only. It does not authorize module fitting,
infer patients from library names, or recover tissue weights from sorted cells.
The evidence review and its stable parent relationships must already exist.
"""
import argparse,collections,csv,datetime,hashlib,json
from pathlib import Path

CORE_FIELDS=('source','qualification_unit_id','qualification_unit_key_sha256',
             'source_qualification','n_archive_rows','native_parent_unit_sha256s')
STATUS=('来源资格通过','待核实','排除')

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
    with Path(p).open('x') as f:json.dump(d,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
def table(p,rows):
    if not rows:raise ValueError('Empty output requires a separately specified schema')
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with Path(p).open('x',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader()
        w.writerows({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v
                     for k,v in r.items()} for r in rows)
def stats(rows):
    return {'work_units_not_patients':len(rows),'archive_cell_records':sum(int(r['n_archive_rows']) for r in rows),
            'sources':len({r['source'] for r in rows})}

def reconcile(old_path,review_path,live_path,out):
    out=Path(out);out.mkdir(exist_ok=True,parents=True)
    old=list(csv.DictReader(Path(old_path).open(encoding='utf-8-sig'),delimiter='\t'))
    by_old={r['单位稳定SHA']:r for r in old}
    if len(old)!=len(by_old):raise ValueError('Duplicate historical parent identity')
    review=load(review_path);units=review['units'];live=load(live_path)
    if len({r['qualification_unit_key_sha256'] for r in units})!=len(units):
        raise ValueError('Duplicate qualification work-unit identity')
    # Reconstruct the actual minimal server authority without reusing display order.
    core=[{k:(r['TME_preanalysis_source_gate_v7'] if k=='source_qualification' else r[k])
           for k in CORE_FIELDS} for r in units]
    core=sorted(core,key=lambda r:r['qualification_unit_key_sha256'])
    digest=hashlib.sha256(json.dumps(core,ensure_ascii=False,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    if digest!=live['core_semantic_SHA256']:
        raise ValueError('Current source review differs from live server authority')
    if live['actual_producer_exit']['actual_exit_code']!=0:
        raise ValueError('Source gate has no accepted actual process exit')
    qualified=[];held=[];excluded=[];allrows=[];parent=collections.defaultdict(list)
    for r in units:
        if r['TME_preanalysis_source_gate_v7'] not in STATUS:raise ValueError('Unknown decision')
        if len(r['native_parent_unit_sha256s'])!=1:
            raise ValueError('Multi-parent work units need an explicit cell allocation contract')
        parent[r['native_parent_unit_sha256s'][0]].append(r)
        row=dict(r)
        row['source_qualification_historical_v6']=r['source_qualification']
        row['source_qualification']=r['TME_preanalysis_source_gate_v7']
        row['current_source_gate_version']='v7'
        row['CM_input_status']='preparation_only_final_identity_QC_specimen_duplicate_gates_pending'
        row['included_in_historical_v2_source_list']=r['native_parent_unit_sha256s'][0] in by_old
        allrows.append(row)
        {'来源资格通过':qualified,'待核实':held,'排除':excluded}[row['source_qualification']].append(row)
    cross=[]
    for key,r in by_old.items():
        children=parent.get(key,[])
        if not children or sum(u['n_archive_rows'] for u in children)!=int(r['细胞记录数']):
            raise ValueError('Historical parent is missing, duplicated, or incompletely allocated')
        if any(u['source']!=r['来源ID'] for u in children):raise ValueError('Source identity differs')
        for u in children:
            cross.append({'source':u['source'],'historical_v2_parent_key_SHA256':key,
                          'current_qualification_unit_key_SHA256':u['qualification_unit_key_sha256'],
                          'historical_v2_decision':'来源资格通过',
                          'current_v7_decision':u['TME_preanalysis_source_gate_v7'],
                          'n_archive_rows':u['n_archive_rows'],
                          'reason':u['v7_additional_preanalysis_hold_reason'] or u['primary_reason']})
    added=[r for r in qualified if not r['included_in_historical_v2_source_list']]
    transition=[]
    for source in sorted({r['source'] for r in cross}):
        for status in STATUS:
            rr=[r for r in cross if r['source']==source and r['current_v7_decision']==status]
            if rr:transition.append({'source':source,'current_v7_decision':status,
                'historical_v2_parent_units':len({r['historical_v2_parent_key_SHA256'] for r in rr}),
                'current_work_units':len(rr),'archive_cell_records':sum(r['n_archive_rows'] for r in rr),
                'reasons':' | '.join(sorted({r['reason'] for r in rr}))})
    summaries=[]
    for source in sorted({r['source'] for r in units}):
        rr=[r for r in allrows if r['source']==source];q=[r for r in rr if r['source_qualification']==STATUS[0]]
        s={'source':source,'cancer_registration':' / '.join(sorted({r['cancer_type'] for r in rr})),
           'archive_cell_records':sum(r['n_archive_rows'] for r in rr),'all_work_units':len(rr)}
        for status,label in zip(STATUS,['prepare','hold','exclude_current_TME_purpose']):
            z=[r for r in rr if r['source_qualification']==status]
            s[label+'_units']=len(z);s[label+'_archive_records']=sum(r['n_archive_rows'] for r in z)
        s.update(material_categories=' / '.join(sorted({r['material'] for r in q})),
                 RNA_modalities=' / '.join(sorted({r['RNA_modality'] for r in q})),
                 sorting_evidence=' | '.join(sorted({r['experimental_sorting'] for r in q})),
                 patient_block_basis_missing_units=sum(not r['patient_block_basis_available'] for r in q),
                 fresh_frozen_unknown_units=sum(r['fresh_frozen']=='待核实' for r in q),
                 treatment_phase_unknown_units=sum(r['treatment_phase'].startswith('待核实') for r in q),
                 final_CM_input_status='pending_final_identity_QC_and_actual_specimen_duplicate_review')
        summaries.append(s)
    for name,rows in [('TME_source_qualified_units_v7.tsv',qualified),('TME_source_held_units_v7.tsv',held),
                      ('TME_source_excluded_current_purpose_units_v7.tsv',excluded),
                      ('all_source_work_units_current_v7.tsv',allrows),('v2_to_current_v7_unit_crosswalk.tsv',cross),
                      ('newly_qualified_since_v2_units_v7.tsv',added),('v2_to_current_v7_source_summary.tsv',transition),
                      ('TME_all65_source_scope_summary_v7.tsv',summaries)]:table(out/name,rows)
    counts={status:stats([r for r in allrows if r['source_qualification']==status]) for status in STATUS}
    for status,s in counts.items():
        if s['work_units_not_patients']!=live['counts'][status]['units'] or s['archive_cell_records']!=live['counts'][status]['archive_rows']:
            raise ValueError('Live scope totals differ')
    files={p.name:sha(p) for p in out.glob('*.tsv')}
    summary={'schema':'TME_source_scope_reconciliation_v1','created_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'historical_v2':{'work_units_not_patients':len(old),'archive_cell_records':sum(int(r['细胞记录数']) for r in old),'sources':len({r['来源ID'] for r in old})},
        'historical_v2_current_disposition':{status:{'work_units':sum(r['current_v7_decision']==status for r in cross),'archive_cell_records':sum(r['n_archive_rows'] for r in cross if r['current_v7_decision']==status)} for status in STATUS},
        'newly_qualified_since_v2':stats(added),'current_counts':counts,
        'all_archive_records':sum(r['n_archive_rows'] for r in units),'full_archive_preserved':True,
        'source_gate_version':'v7','source_manifest_SHA256':live['source_manifest_SHA256'],
        'key_sorted_core_semantic_SHA256':digest,'files':files,
        'source_scoped_parent_keys_used_not_display_names_or_GSE_only':True,'cell_aware_overrides_required_not_native_parent_whitelist':True,
        'final_CM_ready':False,'new_CM_fitted':False,'native_expression_read':False,'old_CM_results_unchanged':True}
    save(out/'private_scope_reconciliation_acceptance_v1.json',summary)
    contract={'schema':'current_TME_module_analysis_objects_v1','created_UTC':summary['created_UTC'],
        'scope':'unselected_tumor_tissue_scRNA_preparation','eligibility_version':'v7',
        'source_view_manifest_SHA256':live['source_manifest_SHA256'],'joined_native_atlas_manifest_SHA256':live['joined_manifest_SHA256'],
        'key_sorted_source_core_SHA256':digest,'source_qualified_objects':counts[STATUS[0]],
        'source_qualified_unit_table':{'path':'TME_source_qualified_units_v7.tsv','SHA256':files['TME_source_qualified_units_v7.tsv']},
        'files':files,'historical_v2_list_is_superseded_for_current_preparation':True,
        'required_cell_level_resolution':'Use the bound source view including inherited per-cell overrides; parent-unit whitelists are unsafe for split author samples.',
        'analysis_units_are_not_independent_patients':True,
        'material_and_phase':'Primary, metastatic, involved lymph-node and unknown tumor site retain separate strata; normal, blood, product, enriched and snRNA are separate uses.',
        'composition':'Unresolved identities stay in the reported denominator; local group numbers do not establish shared biological types.',
        'source_or_cancer_dependence':'Report supported contexts; dependence is not an automatic rejection of a module.',
        'module_rank':'To be evaluated from current data, not inherited as nine modules.',
        'execution_status':'PREPARATION_ONLY_NOT_READY_FOR_MODULE_FIT',
        'pending_gates':['final_common_identity_and_QC','actual_specimen_and_technical_repeat_resolution','source_alias_duplicate_review','composition_coverage_and_analysis_unit_freeze'],
        'clinical_use_pending':'patient_phase_and_endpoint_review','old_dated_CM_versions_remain_historical':True}
    save(out/'current_TME_module_analysis_objects_v1.json',contract)
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--old',required=True);p.add_argument('--review',required=True)
    p.add_argument('--live',required=True);p.add_argument('--output',required=True);a=p.parse_args()
    s=reconcile(a.old,a.review,a.live,a.output)
    print(json.dumps({k:v for k,v in s.items() if k not in ['files']},ensure_ascii=False))
