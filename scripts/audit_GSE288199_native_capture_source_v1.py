"""Independent source-summary conservation and marker-parent invariants.

Does not recreate individual arrays. Case quantiles are checked for support
and bounds, not advertised as an independent per-case quantile replay.
"""
from datetime import datetime, timezone
from itertools import product
from pathlib import Path
import argparse
import hashlib
import json
import math

PARENTS=['dual_regulatory_evidence','costate_regulatory_sensitivity']
QC=['all_released_count_columns','human_RNA_QC500_250_MT20']
DEPTH=['all_capture','low_capture','mid_capture','high_capture']
SORT=['all_source_materials','CD3_CD45_enriched','CD45_enriched']
PHASES=['pre','post']
STATES=['R_TNFRSF9','R_LAG3','R_cycling','E_dysfunction','E_cycling','E_memory']
METRICS=['T_support','regulatory_only','CD8_only','dual_ambiguous','neither_parent']+STATES+['B_cycling','epithelial_cycling']
GENES=['CD3D','CD3E','TRAC','CD4','CD8A','FOXP3','IL2RA','CTLA4','TNFRSF9','LAG3','PDCD1','TOX','HAVCR2','MKI67','TOP2A','UBE2C','TCF7','IL7R','MS4A1','CD79A','EPCAM','SOX2','CXCL13']


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def rates(values):
    if values is not None:
        assert len(values)==5 and values==sorted(values)
        assert all(math.isfinite(v) and 0<=v<=1 for v in values)


def main():
    p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('output',type=Path)
    args=p.parse_args();src=args.source;out=args.output
    assert not out.exists();out.mkdir(mode=0o700)
    t=json.loads((src/'private_terminal.json').read_text())
    assert t['status']=='COMPLETE_ALL60_NATIVE_RNA_OBSERVED_STATE_CAPTURE_NO_EFFECT_TESTS'
    for name,r in t['files'].items():assert sha(src/name)==r['SHA256']
    rows=json.loads((src/'private_all96_native_state_capture_opportunities.json').read_text())
    budget=json.loads((src/'private_complete_header_budget_gate.json').read_text())
    stream=json.loads((src/'private_complete_stream_conservation.json').read_text())
    fields=['phase','sort_scope','QC','parent','capture_scope']
    keys=[tuple(r[k] for k in fields) for r in rows]
    expected=set(product(PHASES,SORT,QC,PARENTS,DEPTH))
    assert len(rows)==96 and len(set(keys))==96 and set(keys)==expected
    lookup=dict(zip(keys,rows));checks=0
    for r in rows:
        n=r['coded_source_units'];with_cells=r['coded_source_units_with_cells']
        assert 0<=with_cells<=n and 0<=r['paired_parent_units_both_atleast20']<=with_cells
        m=r['metric_positive_observations']
        if n<5:
            assert r['cells'] is None and m is None
        else:
            assert m is not None and set(m)==set(METRICS)
            assert all(isinstance(v,int) and 0<=v<=r['cells'] for v in m.values())
            assert sum(m[k] for k in ['regulatory_only','CD8_only','dual_ambiguous','neither_parent'])==r['cells']
            assert m['regulatory_only']+m['CD8_only']+m['dual_ambiguous']<=m['T_support']
            assert m['B_cycling']<=r['cells']-m['T_support'] and m['epithelial_cycling']<=r['cells']-m['T_support']
        assert set(r['state_capture'])==set(STATES)
        for state,s in r['state_capture'].items():
            parent='regulatory_only' if state.startswith('R_') else 'CD8_only'
            supported=s['coded_source_units_parent_atleast20']
            assert 0<=supported<=with_cells
            rates(s['coded_unit_balanced_fraction_min_q25_median_q75_max'])
            if supported<5:assert s['coded_unit_balanced_fraction_min_q25_median_q75_max'] is None
            else:assert s['coded_unit_balanced_fraction_min_q25_median_q75_max'] is not None
            if m is None:assert s['parent_cell_observations'] is None and s['state_positive_observations'] is None
            else:
                assert s['parent_cell_observations']==m[parent] and s['state_positive_observations']==m[state]
                assert m[state]<=m[parent] and supported*20<=m[parent]
        for parent in ['regulatory_only','CD8_only']:
            markers=r['parent_marker_capture'][parent];assert set(markers)==set(GENES)
            for g,v in markers.items():
                rates(v['coded_unit_balanced_detection_min_q25_median_q75_max'])
                if m is None:assert v['positive_observations'] is None
                else:assert 0<=v['positive_observations']<=m[parent]
            if m is not None:
                positive={g:v['positive_observations'] for g,v in markers.items()}
                mandatory='FOXP3' if parent=='regulatory_only' else 'CD8A'
                assert positive[mandatory]==m[parent]
                values=markers[mandatory]['coded_unit_balanced_detection_min_q25_median_q75_max']
                if values is not None:assert all(v==1 for v in values)
                assert positive['CD3D']+positive['CD3E']+positive['TRAC']>=m[parent]
                if parent=='regulatory_only':
                    assert positive['TNFRSF9']==m['R_TNFRSF9'] and positive['LAG3']==m['R_LAG3']
                    assert m['R_cycling']<=positive['MKI67'] and m['R_cycling']<=positive['TOP2A']+positive['UBE2C']
                    if r['parent']==PARENTS[0]:assert positive['IL2RA']==m[parent]
                    else:assert positive['IL2RA']+positive['CTLA4']>=m[parent]
                else:
                    assert max(positive['TCF7'],positive['IL7R'])<=m['E_memory']<=positive['TCF7']+positive['IL7R']
                    assert m['E_dysfunction']<=positive['PDCD1'] and m['E_dysfunction']<=positive['TOX']+positive['HAVCR2']
                    assert m['E_cycling']<=positive['MKI67'] and m['E_cycling']<=positive['TOP2A']+positive['UBE2C']
        checks+=1
    for phase,sort,qc,parent in product(PHASES,SORT,QC,PARENTS):
        whole=lookup[phase,sort,qc,parent,DEPTH[0]]
        children=[lookup[phase,sort,qc,parent,d] for d in DEPTH[1:]]
        assert all(c['coded_source_units']==whole['coded_source_units'] for c in children)
        if whole['metric_positive_observations'] is not None:
            assert sum(c['cells'] for c in children)==whole['cells']
            for key in METRICS:assert sum(c['metric_positive_observations'][key] for c in children)==whole['metric_positive_observations'][key]
            for p in ['regulatory_only','CD8_only']:
                for g in GENES:
                    assert sum(c['parent_marker_capture'][p][g]['positive_observations'] for c in children)==whole['parent_marker_capture'][p][g]['positive_observations']
    for phase,qc,parent,depth in product(PHASES,QC,PARENTS,DEPTH):
        whole=lookup[phase,SORT[0],qc,parent,depth]
        children=[lookup[phase,s,qc,parent,depth] for s in SORT[1:]]
        assert sum(c['coded_source_units'] for c in children)==whole['coded_source_units']
        assert sum(c['coded_source_units_with_cells'] for c in children)==whole['coded_source_units_with_cells']
        assert sum(c['paired_parent_units_both_atleast20'] for c in children)==whole['paired_parent_units_both_atleast20']
        for state in STATES:
            assert sum(c['state_capture'][state]['coded_source_units_parent_atleast20'] for c in children)==whole['state_capture'][state]['coded_source_units_parent_atleast20']
        if whole['metric_positive_observations'] is not None:
            for child in children:
                if child['metric_positive_observations'] is not None:
                    assert child['cells']<=whole['cells']
                    assert all(child['metric_positive_observations'][k]<=whole['metric_positive_observations'][k] for k in METRICS)
    for phase,sort,parent in product(PHASES,SORT,PARENTS):
        before=lookup[phase,sort,QC[0],parent,DEPTH[0]];after=lookup[phase,sort,QC[1],parent,DEPTH[0]]
        if before['metric_positive_observations'] is not None:
            assert after['cells']<=before['cells']
            assert all(after['metric_positive_observations'][k]<=before['metric_positive_observations'][k] for k in METRICS)
    for phase,sort,qc,depth in product(PHASES,SORT,QC,DEPTH):
        first=lookup[phase,sort,qc,PARENTS[0],depth];second=lookup[phase,sort,qc,PARENTS[1],depth]
        if first['metric_positive_observations'] is not None:
            a=first['metric_positive_observations'];b=second['metric_positive_observations']
            assert first['cells']==second['cells'] and a['T_support']==b['T_support']
            assert b['regulatory_only']>=a['regulatory_only'] and b['dual_ambiguous']>=a['dual_ambiguous']
            assert b['CD8_only']<=a['CD8_only'] and b['neither_parent']<=a['neither_parent']
            assert all(b[k]>=a[k] for k in STATES[:3]) and all(b[k]<=a[k] for k in STATES[3:])
    total=sum(lookup[phase,SORT[0],QC[0],PARENTS[0],DEPTH[0]]['cells'] for phase in PHASES)
    assert total==budget['declared_columns']
    assert stream['observed_entries']==budget['declared_nonzero_entries']
    assert stream['source_header_manifest_SHA256']==budget['source_header_manifest_SHA256']
    result=dict(status='PASS_COMPLETE_SOURCE96_PARTITIONS_PARENT_MARKER_QC_SUPPORT_AND_HEADER_CONSERVATION',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(Path(__file__)),
        original_terminal_SHA256=sha(src/'private_terminal.json'),source_opportunities=checks,
        all_QC_all_capture_subsets_checked=True,all_within_QC_depth_partitions_checked=True,
        no_false_QC_subset_requirement_for_recomputed_depth_tertiles=True,
        all_parent_overlap_directions_and_marker_required_genes_checked=True,
        all_header_columns_nnz_source_totals_exact=True,small_source_distributions_suppressed=True,
        case_quantiles_support_and_bounds_only_not_independent_case_array_replay=True,
        raw_numeric_replay_scope='Producer every raw chunk sparse versus coordinate arithmetic and all cell scalar recipes; source aggregate invariants checked here.',
        no_new_clinical_or_source_state_pair_effects=True)
    target=out/'private_source_aggregate_audit.json';target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],receipt_SHA256=sha(target))))


if __name__=='__main__':main()
