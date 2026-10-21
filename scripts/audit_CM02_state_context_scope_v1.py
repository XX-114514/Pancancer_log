"""Independent Counter/Fraction replay of source-level author-state scope."""
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction
from itertools import product
from pathlib import Path
import csv
import hashlib
import json
import os

REG=['4-1BB+ Treg','HLA-DR+ Treg','LAG3+ Treg','Prolif Treg']
OTHER=['Tn/Tmem','IFN-stim CD4T','GZMK+T','Tmem','Teff','Tstr','Tprolif']


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);parent=Path(os.environ['PANCANCER_PARENT_OUTPUT'])
    out=Path(os.environ['PANCANCER_RAM_OUTPUT']);assert out.parent.parent==Path('/dev/shm') and not out.exists()
    terminal=json.loads((parent/'private_terminal.json').read_text())
    for n,v in terminal['files'].items():assert sha(parent/n)==v['SHA256'] and (parent/n).stat().st_size==v['bytes']
    with (root/'HNSCC_author_state_count_layer_preflight_20261016_v1/private_all57_original_author_counts.tsv').open() as f:
        counts={r['library_id']:r for r in csv.DictReader(f,delimiter='\t')}
    cols=['library_id','patient_registry','timepoint_registry','neo_type','HPV/p16']
    with (root/'HNSCC_author_clinical_table_user_20261014_v1/private_key_and_response_audit_v1/all57_private_author_table_join.tsv').open() as f:
        metadata=[{k:r[k] for k in cols} for r in csv.DictReader(f,delimiter='\t')]
    actual=json.loads((parent/'private_all56_state_pair_opportunities.json').read_text())
    indexed={(r['phase'],r['reg_author_state'],r['other_author_T_state']):r for r in actual}
    assert len(indexed)==len(actual)==56 and set(indexed)==set(product(['pre','post'],REG,OTHER))
    phase_rows={r['phase']:r for r in json.loads((parent/'private_all2_phase_scope.json').read_text())}
    for phase in ['pre','post']:
        md=[r for r in metadata if r['timepoint_registry']==phase]
        duplicates=Counter(r['patient_registry'] for r in md)
        raw=[];eligible=[]
        for m in md:
            c={k:int(counts[m['library_id']][k]) for k in REG+OTHER}
            nr=sum(c[k] for k in REG);ne=sum(c[k] for k in OTHER)
            raw.append((m,c,nr,ne))
            if nr>=20 and ne>=20 and m['neo_type'] in ['D','D+T'] and m['HPV/p16'] in ['positive','negative'] and duplicates[m['patient_registry']]==1:
                eligible.append((m,c,nr,ne))
        s=phase_rows[phase]
        assert s['author_materials']==len(md) and s['unique_author_patients']==len(duplicates)
        assert s['duplicate_phase_patients']==sum(n>1 for n in duplicates.values())
        assert s['both_disjoint_parents_atleast20']==sum(nr>=20 and ne>=20 for _,_,nr,ne in raw)
        assert s['eligible_for_count_structure']==len(eligible)
        assert s['missing_or_unsupported_treatment']==sum(m['neo_type'] not in ['D','D+T'] for m in md)
        assert s['missing_or_unsupported_HPV']==sum(m['HPV/p16'] not in ['positive','negative'] for m in md)
        for r in s['by_treatment']:assert r['author_patients']==sum(m['neo_type']==r['arm'] for m,_,_,_ in eligible)
        for a,b in product(REG,OTHER):
            r=indexed[phase,a,b];va={Fraction(c[a],nr) for _,c,nr,_ in eligible};vb={Fraction(c[b],ne) for _,c,_,ne in eligible}
            assert r['eligible_author_patients']==len(eligible)
            assert r['reg_positive_patients']==sum(c[a]>0 for _,c,_,_ in eligible)
            assert r['other_positive_patients']==sum(c[b]>0 for _,c,_,_ in eligible)
            assert r['reg_unique_fractions']==len(va) and r['other_unique_fractions']==len(vb)
            ready=len(eligible)>=20 and len(va)>1 and len(vb)>1
            assert (r['status']=='COUNT_SCOPE_READY_NOT_STATE_OR_CLINICAL_VALIDATION')==ready
            assert r['zero_child_counts_retained'] and r['other_T_not_automatically_CD8_or_cytotoxic'] and r['exact_CM02_per_cell_membership_not_established']
    result=dict(status='ACCEPTED_ALL56_AUTHOR_STATE_OPPORTUNITIES_AND_TWO_PHASE_DENOMINATORS',completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(__file__),parent_terminal_SHA256=sha(parent/'private_terminal.json'),independent_prior_count_table_Counter_and_exact_rational_variation=True,all56_keys_and_coverage_statuses_replayed=True,no_new_outcome_expression_or_association_tests=True,annotation_field_source_read_not_independently_repeated=True,individual_rows_exported=False)
    out.mkdir(mode=0o700);(out/'private_independent_acceptance.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
