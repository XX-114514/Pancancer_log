"""Independent field-total replay from all saved metadata terms, no raw reread."""
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from gamma_author_text_semantics_v2 import KINDS, label_kind


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def as_integer(value):
    number=float(value)
    assert np.isfinite(number) and number>=0 and number==np.floor(number)
    return int(number)


def main():
    root=Path(os.environ['PANCANCER_FINE_GAMMA_PRIVATE_OUTPUT'])
    terminal=root/'private_terminal.json'
    assert sha(terminal)=='117add0af9b1af35cd9194e30df037436cfeadb6b3d90a424aff8c93572d1b82'
    t=json.loads(terminal.read_text())
    for n,d in t['files'].items():
        assert sha(root/n)==d['sha256'] and (root/n).stat().st_size==d['bytes']
    fields=pd.read_csv(root/'private_complete21_annotation_field_qualification.tsv',sep='\t',keep_default_na=False)
    terms=pd.read_csv(root/'private_complete_source_field_author_label_counts.tsv.gz',sep='\t',keep_default_na=False)
    proof=json.loads((root/'private_all21_canonical_selector_and_source_proofs.json').read_text())
    assert len(proof)==21 and fields.record_index.nunique()==21
    assert not fields.duplicated(['record_index','field']).any()
    assert not terms.duplicated(['record_index','field','author_label']).any()
    checked=0
    for row in fields.itertuples():
        counts=terms[(terms.record_index==row.record_index)&(terms.field==row.field)]
        if row.status=='ANNOTATION_VALUES_COUNTED_NOT_REFERENCE_ACCEPTANCE':
            assert len(counts)==as_integer(row.distinct_values) and counts.n_cells.sum()==row.selected_rows
            sums={k:0 for k in KINDS}
            for value in counts.itertuples():
                kind=KINDS[label_kind(value.author_label)]
                assert kind==value.semantic_kind
                sums[kind]+=int(value.n_cells)
            assert sums[KINDS[0]]==as_integer(row.explicit_gamma_cells)
            assert sums[KINDS[1]]==as_integer(row.delta_subtype_text_cells)
            assert sums[KINDS[2]]==as_integer(row.gamma_marker_text_cells)
            assert sums[KINDS[8]]==as_integer(row.mixed_T_NK_text_cells)
            checked+=1
        else:
            assert counts.empty and row.status in ['NO_MATCHING_ANNOTATION_FIELD','HIGH_CARDINALITY_FIELD_NOT_TREATED_AS_IDENTITY_REFERENCE']
    for p in proof:
        actual=fields[fields.record_index==p['record_index']]
        assert len(actual)==max(1,p['field_count'])
        assert actual.selected_rows.eq(p['selected_rows']).all()
        assert p['all_obs_schema_equal_prior'] and p['source_stat_unchanged']
    report=dict(status='ACCEPTED_ALL21_RECORDS_AND_ALL_SAVED_FINE_ANNOTATION_TERM_TOTALS',
                reviewed_UTC=datetime.now(timezone.utc).isoformat(),producer_terminal_sha256=sha(terminal),
                record_positions=21,field_positions=len(fields),measured_field_positions=checked,
                statuses=fields.status.value_counts().to_dict(),source_or_field_overlap_not_independent=True,
                audit_scope='All saved terms and field metadata, not a second canonical source metadata read or receptor/identity accuracy validation.',
                formal_identity_or_clinical_inference_changed=False,scientific_outputs_public=False)
    p=root/'private_independent_field_aggregate_acceptance.json';assert not p.exists();p.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],audit_sha256=sha(p),field_positions=len(fields),measured_field_positions=checked)))


if __name__=='__main__':
    main()
