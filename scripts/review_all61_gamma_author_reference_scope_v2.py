"""Rebuild source-level semantics from all saved author terms; no new RNA."""
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from gamma_author_text_semantics_v2 import KINDS, label_kind, fixture


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    root=Path(os.environ['PANCANCER_GAMMA_AUTHOR_PRIVATE_OUTPUT'])
    terminal=root/'private_terminal.json'
    assert sha(terminal)=='017deb6e9448c22daa22eb929fa1b9c4c1b1782d375acaf785465c5c8f6a4a40'
    t=json.loads(terminal.read_text())
    for n,d in t['files'].items():
        assert sha(root/n)==d['sha256'] and (root/n).stat().st_size==d['bytes']
    fixture()
    terms=pd.read_csv(root/'private_complete_source_author_label_semantics.tsv.gz',sep='\t',keep_default_na=False)
    old=pd.read_csv(root/'private_complete61_source_author_reference_scope.tsv',sep='\t',keep_default_na=False)
    assert len(old)==61 and not terms.duplicated(['arm','dataset','author_label']).any()
    terms['reviewed_semantic_kind']=[KINDS[label_kind(v)] for v in terms.author_label]
    table=terms.groupby(['arm','dataset','reviewed_semantic_kind']).n_cells.sum().unstack(fill_value=0).reindex(columns=KINDS,fill_value=0)
    assert len(table)==61
    legacy=old.set_index(['arm','dataset']).reindex(table.index)
    np.testing.assert_array_equal(table.sum(axis=1),legacy.canonical_cells)
    # The revised NK/T interpretation must leave every gamma reference count exact.
    np.testing.assert_array_equal(table[KINDS[0]],legacy[KINDS[0]])
    for key,section in terms.groupby(['arm','dataset'],sort=True):
        independent={name:0 for name in KINDS}
        for row in section.itertuples():
            independent[row.reviewed_semantic_kind]+=int(row.n_cells)
        assert independent==table.loc[key].to_dict()
    table=table.reset_index()
    for column in ['canonical_cells','libraries','libraries_with_explicit_gamma','explicit_gamma_in_frozen_lymphoid','explicit_gamma_outside_frozen_lymphoid']:
        table[column]=legacy[column].to_numpy()
    table.to_csv(root/'private_complete61_reviewed_source_author_reference_scope_v2.tsv',sep='\t',index=False)
    terms.to_csv(root/'private_all_source_author_terms_semantics_v2.tsv.gz',sep='\t',index=False,compression='gzip')
    report=dict(status='ACCEPTED_ALL61_SOURCE_METADATA_TOTALS_AND_REVIEWED_T_NK_SEMANTICS_NO_NEW_EXPRESSION',
                reviewed_UTC=datetime.now(timezone.utc).isoformat(),producer_terminal_sha256=sha(terminal),
                source_arm_positions=61,term_positions=len(terms),distinct_author_terms=terms.author_label.nunique(),
                gamma_reference_counts_unchanged=True,
                missing_reference_not_cell_absence=True,
                NK_text_correction='Original v1 NK-text category can include T/NK mixed labels; revised v2 separates these and recognizes stand-alone T text. No pure-negative-reference inference from v1.',
                scope='Saved whole source-author term counts checked, not an independent second all-cell metadata read or a receptor identity audit.',
                source_expression_or_clinical_effects_read=False,formal_identity_changed=False,scientific_outputs_public=False)
    p=root/'private_reviewed_source_reference_acceptance_v2.json';assert not p.exists();p.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],audit_sha256=sha(p),source_positions=61,term_positions=len(terms))))


if __name__=='__main__':
    main()
