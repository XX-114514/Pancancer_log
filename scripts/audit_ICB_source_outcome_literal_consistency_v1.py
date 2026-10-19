"""Metadata fingerprints and formatting-only outcome diagnostics, no recoding."""
import hashlib
import json
import os
import shutil
from datetime import datetime,timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata._io.specs import read_elem


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert out.parent.parent==Path('/dev/shm') and not out.exists() and shutil.disk_usage('/dev/shm').free>512*1024**2
    parent=Path('/dev/shm/pancancer-ICB-sample-scope-20261019-v2/private_output')
    assert sha(parent/'private_terminal.json')=='6e5d98bdb3c0b0cd780ca1e14a672222ee15a55776612021b8d143b38ff1049d'
    terminal=json.loads((parent/'private_terminal.json').read_text())
    for name,d in terminal['files'].items():assert sha(parent/name)==d['sha256']
    plan=json.loads((parent/'plan_before_new_original_sample_unit_counts.json').read_text())
    proof=json.loads((parent/'private_current_source_metadata_and_selector_proofs.json').read_text())
    expected=pd.read_csv(parent/'private_complete16_source_population_unit_and_capture_scope.tsv',sep='\t').set_index(['source','population'])
    contract=root.parents[2]/'Project_v5/runs/20261013_immunotherapy_atlas_v1/config/source_contract_v2_private.json'
    assert sha(contract)==plan['source_contract_sha256']
    records=[r for r in json.loads(contract.read_text())['records'] if r['dataset'].startswith('ICB_')]
    paths={r['path'] for r in records};assert len(paths)==1
    missing=set(plan['missing_tokens'])
    out.mkdir(mode=0o700);save=lambda name,value:(out/name).write_text(json.dumps(value,indent=2)+'\n')
    save('plan_before_format_only_outcome_diagnostics.json',dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),parent_terminal_sha256=sha(parent/'private_terminal.json'),scope='POST_SCOPE_DIAGNOSTIC: exact all16 current metadata fingerprints and raw outcome literal counts; stripping/casefold only, no semantic synonym mapping or clinical reclassification. Missing/source-overlap/closed-family limits remain. No RNA, outcomes association or P.'))
    with h5py.File(next(iter(paths)),'r') as h:
        data=pd.DataFrame({name:pd.Series(read_elem(h['obs'][name])).astype(str) for name in plan['metadata_fields']})
    rows=[];checks={(r['source'],r['population']):r for r in proof['checks']}
    for record in records:
        full=np.flatnonzero(data.Study_name.eq(record['dataset'][4:]).to_numpy())
        retained=np.load(record['row_positions_ref'],allow_pickle=False)
        for population,positions in [('contract_retained',retained),('all_available_same_study_in_cache',full)]:
            f=data.iloc[positions].reset_index(drop=True)
            fingerprint=hashlib.sha256(pd.util.hash_pandas_object(f,index=False).to_numpy().tobytes()).hexdigest()
            assert fingerprint==checks[(record['dataset'],population)]['all_selected_metadata_fields_sha256']
            pre=f[f.pre_post.eq('Pre')&~f.PMID_donor_id.str.strip().str.casefold().isin(missing)]
            raw=0;formatted=0;units=0
            for _,g in pre.groupby('PMID_donor_id',sort=True):
                values=g.outcome[~g.outcome.str.strip().str.casefold().isin(missing)]
                raw+=int(values.nunique()>1)
                formatted+=int(values.str.strip().str.casefold().nunique()>1)
                units+=1
            assert raw==expected.loc[(record['dataset'],population),'outcome_within_donor_phase_conflicts']
            assert units==expected.loc[(record['dataset'],population),'baseline_author_identifier_tokens']
            assert formatted<=raw
            rows.append(dict(source=record['dataset'],population=population,baseline_author_identifier_tokens=units,raw_literal_conflict_units=raw,case_whitespace_normalized_conflict_units=formatted,conflicts_explained_by_case_or_whitespace=raw-formatted,clinical_endpoint_or_patient_identity_accepted=False,nontrivial_literal_difference_not_automatically_clinical_contradiction=True))
    assert len(rows)==16
    pd.DataFrame(rows).to_csv(out/'private_complete16_format_only_outcome_consistency.tsv',sep='\t',index=False)
    report=dict(status='ACCEPTED_ALL16_CURRENT_METADATA_FINGERPRINTS_AND_FORMAT_ONLY_OUTCOME_DIAGNOSTICS',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),parent_terminal_sha256=sha(parent/'private_terminal.json'),current_metadata_fingerprints_checked=16,original_capture_numbers_already_reproduced=40,no_semantic_response_recoding=True,no_new_clinical_inference_or_RNA=True,scientific_outputs_public=False,files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()})
    assert sum(d['bytes'] for d in report['files'].values())<512*1024
    save('private_terminal.json',report);print(json.dumps(dict(status=report['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':
    main()
