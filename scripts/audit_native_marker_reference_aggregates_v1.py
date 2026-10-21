"""Audit all native marker source summaries against accepted reference counts.

Uses source aggregates only. It cannot independently reconstruct case-wise
rates without individual data, and makes that limitation explicit.
"""
from collections import Counter
from datetime import datetime,timezone
from itertools import product
from pathlib import Path
import argparse
import hashlib
import json


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main(a):
    terminal=json.loads((a.native/'private_terminal.json').read_text())
    prior=json.loads((a.reference/'private_terminal.json').read_text())
    for root,t in [(a.native,terminal),(a.reference,prior)]:
        for n,v in t['files'].items():
            if (root/n).exists():assert sha(root/n)==v['SHA256']
    rows=json.loads((a.native/'private_all1080_native_marker_reference_opportunities.json').read_text())
    background=json.loads((a.native/'private_all_source_major_capture_background.json').read_text())
    reference=json.loads((a.reference/'private_unique_fine_state_crosscounts.json').read_text())
    denominators=json.loads((a.reference/'private_all_node_phase_reference_denominators.json').read_text())
    contexts={'full_current_source':None,'CM02_regulatory_naming_anchor':'39697d077b6cb8bf',
              'CM02_cytotoxic_naming_anchor':'c3f761a9da61f578'}
    keys=['phase','context','capture_scope','regulatory_parent','reference_label']
    indexed={tuple(r[k] for k in keys):r for r in rows}
    assert len(rows)==len(indexed)==1080
    axes=[sorted({r[k] for r in rows}) for k in keys]
    assert set(indexed)==set(product(*axes))
    native_depth_names=['low_native_capture','mid_native_capture','high_native_capture']
    for phase,context,parent,label in product(axes[0],axes[1],axes[3],axes[4]):
        row=indexed[phase,context,'all_native_filtered',parent,label]
        node=contexts[context]
        matched=lambda r:r['phase']==phase and (node is None or r['current_node']==node)
        if label=='no_unique_fine_reference':
            expected=sum(r['no_fine_reference_cells'] for r in denominators if matched(r))
        elif label=='overlapping_fine_reference':
            expected=sum(r['overlapping_reference_cells'] for r in denominators if matched(r))
        else:
            expected=sum(r['matched_cells'] for r in reference if matched(r) and r['author_label']==label)
        assert row['cell_observations']==expected
        parts=[indexed[phase,context,d,parent,label] for d in native_depth_names]
        assert sum(r['cell_observations'] for r in parts)==expected
        for name,v in row['measurements'].items():
            assert sum(r['measurements'][name]['positive_cell_observations'] for r in parts)==v['positive_cell_observations']
    for row in rows:
        n=row['cell_observations'];m=row['measurements']
        assert sum(m[k]['positive_cell_observations'] for k in ['regulatory_only','CD8_only','dual_ambiguous','neither_parent'])==n
        for k,v in m.items():
            positive=v['positive_cell_observations'];assert 0<=positive<=n
            assert v['pooled_cell_observation_fraction']==(positive/n if n else None)
            q=v['case_balanced_rate_min_q25_median_q75_max']
            if row['reviewed_cases_atleast20_reference_cells']>=5:
                assert q is not None and len(q)==5 and q==sorted(q) and 0<=q[0]<=q[-1]<=1
            else:assert q is None
            if k.startswith('R_'):assert positive<=m['regulatory_only']['positive_cell_observations']
            if k.startswith('E_'):assert positive<=m['CD8_only']['positive_cell_observations']
    for phase,parent in product(axes[0],axes[3]):
        group=[r for r in rows if r['phase']==phase and r['regulatory_parent']==parent and r['context']=='full_current_source' and r['capture_scope']=='all_native_filtered']
        bg=[r for r in background if r['phase']==phase and r['regulatory_parent']==parent]
        assert sum(r['cell_observations'] for r in group)==sum(r['cell_observations'] for r in bg)
        for metric in group[0]['measurements']:
            assert sum(r['measurements'][metric]['positive_cell_observations'] for r in group)==sum(r['observed_marker_counts'][metric] for r in bg)
    result=dict(status='ALL1080_NATIVE_REFERENCE_SOURCE_DENOMINATORS_DEPTH_AND_MARKER_PARTITIONS_REPLAYED',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(__file__),
        native_terminal_SHA256=sha(a.native/'private_terminal.json'),prior_reference_terminal_SHA256=sha(a.reference/'private_terminal.json'),
        all_context_phase_parent_reference_denominators_match_previous_accepted_crosscounts=True,
        all_depth_partitions_parent_partitions_marker_nesting_and_major_sums_checked=True,
        all_fixed_opportunities_missing_scopes_and_rate_bounds_preserved=True,
        case_rate_quantiles_bounds_only_not_independent_per_case_replay=True,
        native_all_cell_scalar_and_raw_CSR_proofs_accepted_separately=True,
        author_reference_same_study_not_independent_biological_truth=True,
        no_patient_cell_rows_expression_arrays_or_new_effects_read=True)
    assert not a.output.exists();a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for n in ['native','reference','output']:p.add_argument('--'+n,type=Path,required=True)
    main(p.parse_args())
