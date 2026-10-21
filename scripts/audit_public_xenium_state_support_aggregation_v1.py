"""Independent Counter aggregation of all frozen support opportunities.

Reuses the qualified input decoder and all-cell scalar checks. Reconstructs
section eligibility and source-case grouping separately; no expression reread.
"""
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import os

import xenium_observed_state_context_v2 as context


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    root=Path(os.environ['PANCANCER_OBSERVED_SUPPORT']);terminal=context.verified(root)
    assert terminal['all17_measured_sections_all_cells_four_modes_scalar_checked']
    rows=json.loads((root/'private_all18_source_slot_marker_support.json').read_text())
    opportunities=json.loads((root/'private_all36_source_state_support_opportunities.json').read_text())
    pooled=json.loads((root/'private_all4_source_marker_observation_counts.json').read_text())
    anatomy=json.loads((root/'private_source_anatomical_scope.json').read_text())
    ctx=context.load_context()
    assert len(rows)==18 and [r['ordinal'] for r in rows]==list(range(18))
    metadata=[p['source_metadata'] for p in ctx['pairs']]
    selected=[i for i,m in enumerate(metadata) if m['diagnosis']=='HNSCC']
    assert all(rows[i]['excluded']==(i not in selected) for i in range(18))
    assert len(selected)==17 and len({metadata[i]['case'] for i in selected})==10
    for mi,target in enumerate(pooled):
        counts=Counter();rc=Counter();ec=Counter();rp=Counter();ep=Counter()
        for k in selected:
            m=rows[k]['modes'][mi];counts.update(m['counts'])
            assert m['all_cell_observations_independent_scalar_checked']
            assert rows[k]['raw_CSC_marker_columns_and_human_totals_independent_bincount_checked']
            for j,v in enumerate(m['R_observed_state_counts']):rc[j]+=v
            for j,v in enumerate(m['E_observed_state_counts']):ec[j]+=v
            for j,v in enumerate(m['R_joint_state_patterns']):rp[j]+=v
            for j,v in enumerate(m['E_joint_state_patterns']):ep[j]+=v
        assert dict(counts)==target['cell_observation_slots_across_sections_not_patient_prevalence']
        assert [rc[j] for j in range(3)]==target['R_observed_counts']
        assert [ec[j] for j in range(3)]==target['E_observed_counts']
        assert [rp[j] for j in range(8)]==target['R_patterns']
        assert [ep[j] for j in range(8)]==target['E_patterns']
        assert sum(rp.values())==counts['regulatory'] and sum(ep.values())==counts['CD8']
        for pi in range(9):
            i,j=divmod(pi,3);eligible=[];statuses=Counter()
            for k in selected:
                m=rows[k]['modes'][mi];nr=m['counts']['regulatory'];ne=m['counts']['CD8']
                cr=m['R_observed_state_counts'][i];ce=m['E_observed_state_counts'][j]
                if nr<20 or ne<20:s='INSUFFICIENT_PARENT_OBSERVATIONS'
                elif cr<5 or ce<5:s='INSUFFICIENT_OBSERVED_CHILD_SUPPORT'
                elif cr==nr or ce==ne:s='NO_WITHIN_PARENT_STATE_VARIATION'
                else:s='OBSERVED_STATE_COUNT_SCOPE_READY';eligible.append(k)
                assert m['opportunities'][pi]['status']==s
                statuses[s]+=1
            target_op=opportunities[mi*9+pi]
            assert target_op['section_status_counts']==dict(statuses)
            assert target_op['supporting_sections']==len(eligible)
            assert target_op['supporting_source_cases']==len({metadata[k]['case'] for k in eligible})
            for phase in ['primary','recurrent']:
                assert target_op['supporting_'+phase+'_source_cases']==len({metadata[k]['case'] for k in eligible if metadata[k]['recurrence']==phase})
    assert anatomy['case_counts_by_anatomical_site']=={s:len({metadata[k]['case'] for k in selected if metadata[k]['anatomical_site']==s}) for s in {metadata[k]['anatomical_site'] for k in selected}}
    out=Path(os.environ['PANCANCER_RAM_OUTPUT']);assert out.parent.parent==Path('/dev/shm') and not out.exists();out.mkdir(mode=0o700)
    result=dict(status='ALL36_SOURCE_CASE_SUPPORT_AND_ALL4_OBSERVATION_AGGREGATES_INDEPENDENT_COUNTER_REPLAYED',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(__file__),
        context_SHA256=sha(context.__file__),support_terminal_SHA256=sha(root/'private_terminal.json'),
        full36_opportunities_and_section_missing_reasons_recomputed=True,
        no_cross_section_endpoint_union_or_technical_replicate_as_patient=True,
        all4_marker_counts_and_joint_patterns_reconstructed=True,
        all18_payload_source_and_all_cell_scalar_checks_previously_accepted=True,
        source_anatomical_categories_recomputed=True,
        same_decoder_not_claimed_independent_source_metadata_parser=True,
        observational_RNA_pattern_not_cell_identity_or_function=True,
        common_cycling_markers_require_cross_lineage_specificity_and_RNA_bridge=True,
        no_new_spatial_or_clinical_effects=True,no_case_or_cell_rows_saved=True)
    (out/'private_independent_aggregation_acceptance.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':main()
