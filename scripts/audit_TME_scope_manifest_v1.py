"""Read-back audit of historical-to-current TME source registries.

No producer helpers are imported. The audit checks exact membership, source
decisions, split-parent accounting and the non-executable preparation contract.
"""
import argparse,csv,datetime,hashlib,json
from pathlib import Path

def load(p):return json.loads(Path(p).read_text())
def rows(p):return list(csv.DictReader(Path(p).open(encoding='utf-8-sig'),delimiter='\t'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main(a):
    out=Path(a.output);d=Path(a.delivery);review=load(a.review)['units'];old=rows(a.old);live=load(a.live)
    accepted=load(d/'private_scope_reconciliation_acceptance_v1.json')
    contract=load(d/'current_TME_module_analysis_objects_v1.json')
    for name,h in accepted['files'].items():assert sha(d/name)==h,(name,'hash')
    authoritative={u['qualification_unit_key_sha256']:u for u in review}
    assert len(authoritative)==len(review)
    allrows=rows(d/'all_source_work_units_current_v7.tsv')
    actual={u['qualification_unit_key_sha256']:u for u in allrows}
    assert len(actual)==len(allrows) and set(actual)==set(authoritative)
    for key,u in authoritative.items():
        r=actual[key]
        assert r['source']==u['source'] and int(r['n_archive_rows'])==u['n_archive_rows']
        assert r['source_qualification']==u['TME_preanalysis_source_gate_v7']
        assert r['source_qualification_historical_v6']==u['source_qualification']
        assert json.loads(r['native_parent_unit_sha256s'])==u['native_parent_unit_sha256s']
        assert r['CM_input_status']=='preparation_only_final_identity_QC_specimen_duplicate_gates_pending'
    partitions={}
    for status,name in [('来源资格通过','TME_source_qualified_units_v7.tsv'),('待核实','TME_source_held_units_v7.tsv'),('排除','TME_source_excluded_current_purpose_units_v7.tsv')]:
        found=rows(d/name);keys=[r['qualification_unit_key_sha256'] for r in found]
        expected={k for k,u in authoritative.items() if u['TME_preanalysis_source_gate_v7']==status}
        assert len(keys)==len(set(keys)) and set(keys)==expected
        # Separate tables have different optional evidence columns; absent and
        # empty serialized fields agree, while all nonempty values must match.
        assert all(all(r.get(f,'')==actual[r['qualification_unit_key_sha256']].get(f,'')
                       for f in set(r)|set(actual[r['qualification_unit_key_sha256']])) for r in found)
        partitions[status]=set(keys)
        assert sum(int(r['n_archive_rows']) for r in found)==live['counts'][status]['archive_rows']
    assert not(partitions['来源资格通过']&partitions['待核实'] or partitions['来源资格通过']&partitions['排除'] or partitions['待核实']&partitions['排除'])
    oldby={r['单位稳定SHA']:r for r in old};assert len(oldby)==len(old)
    links=rows(d/'v2_to_current_v7_unit_crosswalk.tsv');seen=set()
    for r in links:
        k=r['current_qualification_unit_key_SHA256'];u=authoritative[k];parent=r['historical_v2_parent_key_SHA256']
        assert k not in seen;seen.add(k)
        assert parent in u['native_parent_unit_sha256s'] and parent in oldby
        assert r['source']==u['source']==oldby[parent]['来源ID']
        assert int(r['n_archive_rows'])==u['n_archive_rows'] and r['current_v7_decision']==u['TME_preanalysis_source_gate_v7']
    for key,r in oldby.items():
        rr=[v for v in links if v['historical_v2_parent_key_SHA256']==key]
        assert rr and sum(int(v['n_archive_rows']) for v in rr)==int(r['细胞记录数'])
    expected_added={k for k in partitions['来源资格通过'] if not set(authoritative[k]['native_parent_unit_sha256s'])&set(oldby)}
    added=rows(d/'newly_qualified_since_v2_units_v7.tsv')
    assert len(added)==len(expected_added) and {r['qualification_unit_key_sha256'] for r in added}==expected_added
    sources=rows(d/'TME_all65_source_scope_summary_v7.tsv')
    assert len(sources)==len({u['source'] for u in review})
    for r in sources:
        for status,label in [('来源资格通过','prepare'),('待核实','hold'),('排除','exclude_current_TME_purpose')]:
            z=[u for u in review if u['source']==r['source'] and u['TME_preanalysis_source_gate_v7']==status]
            assert int(r[label+'_units'])==len(z) and int(r[label+'_archive_records'])==sum(u['n_archive_rows'] for u in z)
    assert contract['source_view_manifest_SHA256']==live['source_manifest_SHA256']
    assert contract['execution_status']=='PREPARATION_ONLY_NOT_READY_FOR_MODULE_FIT'
    assert accepted['final_CM_ready'] is False and accepted['new_CM_fitted'] is False
    assert 'overrides' in contract['required_cell_level_resolution']
    assert contract['files']==accepted['files']
    report={'schema':'independent_TME_scope_manifest_audit_v1','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'status':'COMPLETE_MEMBERSHIP_DECISION_PARENT_ALLOCATION_AND_PREPARATION_CONTRACT_ACCEPTED',
        'all_work_units_checked':len(actual),'historical_parent_units_checked':len(oldby),
        'all_archive_records':sum(int(r['n_archive_rows']) for r in allrows),
        'all_delivered_tables_hash_checked':True,'all_source_partition_memberships_exact':True,
        'historical_list_not_modified':True,'final_CM_not_authorized_by_source_gate':True,
        'producer_helpers_imported':False,'native_expression_read':False,
        'contract_SHA256':sha(d/'current_TME_module_analysis_objects_v1.json'),
        'source_view_manifest_SHA256':live['source_manifest_SHA256']}
    with out.open('x') as f:json.dump(report,f,ensure_ascii=False,indent=2);f.write('\n')
    print(json.dumps(report))

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for k in ['old','review','delivery','live','output']:p.add_argument('--'+k,required=True)
    main(p.parse_args())
