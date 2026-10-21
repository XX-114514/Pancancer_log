"""Independently accept complete source-only slots, hashes and denominators.

This is numerical/contract acceptance, not biological or clinical validation.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import math


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):return json.loads(path.read_text())


def check_interval(v,maximum):
    n=v['estimable_images'];assert isinstance(n,int) and 0<=n<=maximum
    assert 0<=v['positive_images']+v['negative_images']<=n
    if n==0:
        assert v['median'] is None and v['IQR'] is None and v['leave_one_image_out_median_range'] is None
    else:
        median=v['median'];q=v['IQR'];assert all(math.isfinite(x) for x in [median]+q)
        assert q[0]<=median<=q[1]
        r=v['leave_one_image_out_median_range']
        if n==1:assert r is None
        else:assert len(r)==2 and all(math.isfinite(x) for x in r) and r[0]<=median<=r[1]


def accept(task,code):
    task=Path(task);code=Path(code)
    plan=read(task/'private_frozen_measurement_v1.json');mt=read(task/'private_measurement_terminal_v1.json')
    m=read(task/'private_complete_measurement_v1.json');c=read(task/'private_frozen_description_v1.json')
    d=read(task/'private_complete_description_v1.json');t=read(task/'private_terminal_v1.json')
    qa=read(task/'private_description_synthetic_acceptance_v2.json')
    assert t['source_exit_code']==0 and t['source_RAM_released'] and t['no_raw_or_individual_files']
    assert t['status']=='COMPLETE_DESCRIPTION_RAM_RELEASED' and d['status']=='COMPLETE_AUTHOR_TYPE_SPATIAL_DESCRIPTION_NO_INFERENCE'
    assert mt['frozen_contract_SHA256']==sha(task/'private_frozen_measurement_v1.json')==t['frozen_contract_SHA256']
    assert c['measurement_terminal_SHA256']==sha(task/'private_measurement_terminal_v1.json')
    assert c['combined_source_RAM_array_SHA256']==m['combined_source_RAM_array_SHA256']
    assert t['description_contract_SHA256']==sha(task/'private_frozen_description_v1.json')
    assert qa['code_SHA256']==c['code_SHA256']==sha(code/c['code_filename'])
    for name,digest in plan['code_SHA256'].items():assert sha(code/name)==digest
    n=m['images'];assert n==len(plan['sources'])==t['all_files']==d['image_count']
    assert m['rows']==plan['expected_total_rows'] and m['coordinates_finite']-m['duplicate_xy_rows_excluded']==m['coordinates_usable']
    assert m['usable_CD4']==m['usable_CD4_Treg']+m['usable_CD4_nonTreg']
    counts=m['concordant_fine_T_counts']
    assert sum(v for k,v in counts.items() if k.startswith('CD4_'))==m['usable_CD4']
    assert sum(v for k,v in counts.items() if k.startswith('CD8_'))==m['usable_CD8']
    assert counts.get('CD4_Treg',0)==m['usable_CD4_Treg']
    assert m['all_exact_coordinate_and_fine_parent_independent_replays_accepted']
    for field in m['fields'].values():assert sum(field['qualified_counts'].values())+field['unqualified_rows']+field['missing_rows']==m['rows']
    expected={(label,k) for label in c['all_CD8_labels'] for k in c['k_values']}
    slots=d['all_spatial_slots'];assert len(slots)==len(expected)
    assert {(s['CD8_author_label'],s['k']) for s in slots}==expected
    for s in slots:
        assert s['primary']==(s['k']==10) and s['ambiguous_CD8_Treg']==(s['CD8_author_label']=='CD8_Treg')
        for key,reasons in [('global_native_graph','missing_reasons'),('same_author_region_MH','region_missing_reasons')]:
            check_interval(s[key],n)
            assert s[key]['estimable_images']+sum(s[reasons].values())==n
        assert s['same_author_region_MH']['estimable_images']<=s['global_native_graph']['estimable_images']
        assert s['global_native_graph']['estimable_images']<=m['images_meeting_future_20each_fourgroup_gate'].get(s['CD8_author_label'],0)
    regions=d['all_regional_composition_slots'];expected={(r,k) for r in c['all_qualified_region_labels'] for k in ['CD4_Treg']+c['all_CD8_labels']}
    assert len(regions)==len(expected) and {(r['author_region_label'],r['author_type']) for r in regions}==expected
    for r in regions:
        check_interval(r['relative_to_same_parent_elsewhere'],n)
        assert r['relative_to_same_parent_elsewhere']['estimable_images']+r['unavailable_images']==n
    assert d['independent_flat_edge_replay_all_estimable_designs'] and d['direct_dense_neighbor_checks_first_middle_last_CD4_per_image_all_k']
    files=['private_frozen_measurement_v1.json','private_measurement_terminal_v1.json','private_complete_measurement_v1.json',
           'private_frozen_description_v1.json','private_complete_description_v1.json','private_terminal_v1.json']
    result=dict(status='COMPLETE_SOURCE_NUMERICAL_AND_CONTRACT_ACCEPTANCE_NOT_BIOLOGICAL_CONFIRMATION',observed_UTC=datetime.now(timezone.utc).isoformat(),
        files_SHA256={f:sha(task/f) for f in files},source_exit0_and_RAM_released=True,complete_slots_and_missing_denominators=True,
        original_and_handoff_code_SHA_rechecked=True,source_independent_geometry_hierarchy_and_edge_replay_accepted=True,
        no_independent_patient_phase_endpoint_RNA_state_or_original_CM_exact_replication_claim=True,visual_QA_still_pending=True)
    path=task/'private_numeric_acceptance_v1.json';assert not path.exists();path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'])


if __name__=='__main__':
    import sys
    accept(sys.argv[1],Path(__file__).parent)
