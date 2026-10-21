"""Reuse complete pre-existing source aggregate capture, never individual arrays.

The frozen contract supplies biological groups. No new expression access,
patient linkage, hypothesis tests, weights or matrix output.
"""
from pathlib import Path
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import statistics

METRICS=('cells_any_gene','cells_ge2_genes','cells_ge2_genes_each_value_ge2')


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(path.read_text())


def interval(values):
    values=sorted(values)
    if not values:return dict(files=0,median=None,minimum=None,maximum=None)
    return dict(files=len(values),median=statistics.median(values),minimum=min(values),maximum=max(values))


def audit(task):
    task=Path(task);contract=task/'private_frozen_existing_capture_audit_v1.json';plan=read(contract)
    assert sha(Path(__file__))==plan['code_SHA256']
    folder=Path(plan['old_capture_folder']);terminal=read(folder/'private_terminal.json')
    assert sha(folder/'private_terminal.json')==plan['old_terminal_SHA256']
    assert terminal['status']=='COMPLETE_DECLARED_CELLBIN_CAPTURE_SCOPE_NO_SPATIAL_EFFECTS'
    spatial=read(Path(plan['accepted_metadata_summary']));assert sha(Path(plan['accepted_metadata_summary']))==plan['accepted_metadata_SHA256']
    inputs=read(Path(plan['accepted_spatial_frozen_measurement']))
    expected={v['SHA256'] for v in inputs['sources']};seen=set();rows=0;group_records={};all_names=Counter()
    for name,info in terminal['files'].items():
        if not name.endswith('_private_qualification.json'):continue
        path=folder/name;assert sha(path)==info['sha256'] and path.stat().st_size==info['bytes']
        source=read(path);key=source['public_file_sha256'];assert key in expected and key not in seen;seen.add(key)
        rows+=source['observations'];fine=[r for r in source['fixed_panel_capture'] if r['annotation_field']=='anno2']
        panels=sorted({r['panel'] for r in fine});assert panels==plan['all_original_panels']
        index={(r['label'],r['panel']):r for r in fine};nlabel={}
        for label,panel in index:
            rec=index[label,panel];assert rec['axis_complete'] and rec['unique_measured_genes']==rec['required_genes']
            assert label not in nlabel or nlabel[label]==rec['n'];nlabel[label]=rec['n']
            a,b,c=[rec[k] for k in METRICS];assert 0<=c<=b<=a<=rec['n']
        all_names.update(nlabel)
        for group,labels in plan['groups'].items():
            for panel in panels:
                selected=[index[label,panel] for label in labels if (label,panel) in index]
                n=sum(r['n'] for r in selected);counts={metric:sum(r[metric] for r in selected) for metric in METRICS}
                group_records.setdefault((group,panel),[]).append(dict(n=n,counts=counts))
    assert seen==expected and len(seen)==spatial['images'] and rows==spatial['rows']
    assert all(all_names[name]==count for name,count in spatial['concordant_fine_T_counts'].items())
    results=[]
    for (group,panel),records in sorted(group_records.items()):
        n=sum(v['n'] for v in records);counts={metric:sum(v['counts'][metric] for v in records) for metric in METRICS}
        # Independently flatten all existing qualified rows into a Counter tally.
        tally=Counter()
        for r in records:
            tally['n']+=r['n']
            for metric,value in r['counts'].items():tally[metric]+=value
        assert tally['n']==n and all(tally[k]==v for k,v in counts.items())
        results.append(dict(author_group=group,original_panel=panel,observations=n,
            marker_capture_counts=counts,source_pooled_observation_fractions={k:v/n if n else None for k,v in counts.items()},
            equal_file_capture_fraction_summary={metric:interval([v['counts'][metric]/v['n'] for v in records if v['n']]) for metric in METRICS},
            no_pooled_cell_median_reconstructed_from_group_medians=True))
    contrasts=[]
    for first,second in plan['all_group_contrasts']:
        for panel in plan['all_original_panels']:
            left=group_records[first,panel];right=group_records[second,panel];assert len(left)==len(right)
            for metric in METRICS:
                values=[a['counts'][metric]/a['n']-b['counts'][metric]/b['n'] for a,b in zip(left,right) if a['n'] and b['n']]
                contrasts.append(dict(first=first,second=second,panel=panel,metric=metric,
                    paired_same_file_fraction_difference=interval(values),positive_files=sum(v>0 for v in values),negative_files=sum(v<0 for v in values)))
    result=dict(status='COMPLETE_EXISTING_NATIVE_CAPTURE_REUSE_NO_NEW_EXPRESSION_OR_TESTS',completed_UTC=datetime.now(timezone.utc).isoformat(),
        source='GSE328481',source_files=len(seen),source_observations=rows,all_group_panel_capture=results,all_paired_capture_descriptions=contrasts,
        source_group_labels_and_all_native_input_SHA_reconciled=True,all_capture_counts_replayed_independently=True,
        no_new_expression_patient_clinical_or_cell_arrays=True,no_marker_specificity_or_formal_identity_truth_claim=True,
        panel_CXCL13_or_gammadelta_marker_validation_not_present=True,original_panel_definition_not_changed=True)
    out=task/'private_complete_existing_RNA_capture_audit_v1.json';assert not out.exists();out.write_text(json.dumps(result,indent=2)+'\n')
    report=dict(status=result['status'],completed_UTC=result['completed_UTC'],contract_SHA256=sha(contract),results_SHA256=sha(out),
                source_terminal_and_all_file_SHA_checked=True,raw_category_or_object_identifiers_never_emitted=True)
    out=task/'private_audit_terminal_v1.json';assert not out.exists();out.write_text(json.dumps(report,indent=2)+'\n')
    print(report['status'])


if __name__=='__main__':
    import sys
    audit(sys.argv[1])
