"""Validate source-only phase summaries; suppress small source-phase values.

Preserve original private outputs. No clinical or correlation recalculation.
"""
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import argparse
import hashlib
import itertools
import json
import math


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('output',type=Path)
    args=p.parse_args();src=args.source;out=args.output
    assert not out.exists();out.mkdir(mode=0o700)
    terminal=json.loads((src/'private_terminal.json').read_text())
    assert terminal['status']=='COMPLETE_OUTCOME_FREE_PHASE_ACTIVITY_STATE_CORRESPONDENCE'
    for name in ['private_all48_phase_activity_state_correspondences.json','private_source_phase_projection_profiles.json','private_original_model_parity.json']:
        assert sha(src/name)==terminal['files'][name]['SHA256']
    rows=json.loads((src/'private_all48_phase_activity_state_correspondences.json').read_text())
    fields=['parent','context','phase','state']
    expected=set(itertools.product(['dual_regulatory_evidence','costate_regulatory_sensitivity'],
        ['full_current_source','CM02_network_member_union'],['pre','post'],
        ['R_TNFRSF9','R_LAG3','R_cycling','E_dysfunction','E_cycling','E_memory']))
    assert len(rows)==48 and {tuple(r[k] for k in fields) for r in rows}==expected
    for r in rows:
        assert r['source']=='GSE286827'
        assert r['reviewed_cases']+r['parent_or_T_support_excluded']==r['projected_phase_cases']
        if r['status']=='DESCRIPTIVE_SOURCE_PHASE_CORRESPONDENCE':
            assert r['reviewed_cases']>=20 and r['bootstrap_manual_tie_rank_replay']
            assert r['no_P_or_FDR'] and math.isfinite(r['spearman_rho']) and abs(r['spearman_rho'])<=1+1e-12
            assert r['bootstrap_CI'] is None or (r['bootstrap_valid']>=1800 and len(r['bootstrap_CI'])==2 and r['bootstrap_CI'][0]<=r['bootstrap_CI'][1])
            assert r['leave_one_case_same_sign']<=r['leave_one_case_valid']<=r['reviewed_cases']
            assert r['reviewed_cases']-1>=20 or r['leave_one_case_valid']==0
    phase=json.loads((src/'private_source_phase_projection_profiles.json').read_text())
    assert {r['phase'] for r in phase}=={'pre','post','post2'}
    suppressed=0
    for r in phase:
        if r['reviewed_cases']<5:
            for k in ['activity_min_q25_median_q75_max','relative_projection_residual_min_q25_median_q75_max']:
                r[k]=None
            r['status']='SMALL_SOURCE_PHASE_DISTRIBUTION_SUPPRESSED'
            suppressed+=1
        else:r['status']='DESCRIPTIVE_SOURCE_PHASE_PROJECTION'
    def save(n,v):(out/n).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
    save('private_all48_phase_activity_state_correspondences.json',rows)
    save('private_source_phase_projection_profiles.json',phase)
    save('private_source_output_qualification.json',dict(
        status='PASS_COMPLETE_SOURCE_SCOPE_AND_SMALL_PHASE_SUPPRESSION',
        observed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(Path(__file__)),
        original_terminal_SHA256=sha(src/'private_terminal.json'),
        all48_unique_frozen_opportunities=True,status_counts=dict(Counter(r['status'] for r in rows)),
        all_source_denominators_and_interval_LOO_bounds_checked=True,
        source_phase_distribution_minimum_cases=5,suppressed_source_phase_distributions=suppressed,
        original_private_history_preserved=True,no_effects_or_bootstraps_recalculated=True,
        files={p.name:dict(SHA256=sha(p),bytes=p.stat().st_size) for p in out.iterdir()}))
    print(json.dumps(dict(status='PASS_COMPLETE_SOURCE_SCOPE_AND_SMALL_PHASE_SUPPRESSION',
        receipt_SHA256=sha(out/'private_source_output_qualification.json'))))


if __name__=='__main__':main()
