"""Source-level author/frozen role cross-tabs; no per-cell or learned-gene export."""
import json
import re
from datetime import datetime,timezone

import numpy as np
import pandas as pd

from NSCLC_external_lineage_inputs_v2 import load,cell_arms,matched_quota,sha,VIEWS


def author_disagrees(frozen, parsed, text):
    if frozen < 0:
        return False
    if parsed >= 0:
        return frozen != parsed
    # Explicit CD4 is contrary to a CD8 endpoint, while it is not proof of Treg.
    s=re.sub(r'[_-]+',' ',str(text).casefold())
    return frozen==0 and bool(re.search(r'\bcd4(?:\b|t)',s)) and not bool(re.search(r'\bcd8(?:\b|t)',s))


def main():
    data=load();root=data['root'];matched_quota(root)
    prior=root/'NSCLC_native_MFA_external_lineage_qualification_20261016_v2/private_independent_lineage_case_audit.json'
    assert sha(prior)=='b07bf0d9b682e8b6cc2ef72192a0ff110a722ba99bd7a931c5e7745c7640699a'
    out=root/'NSCLC_native_MFA_external_author_role_crosscheck_20261016_v1'
    assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),
        upstream_audit_sha256=sha(prior),input_sha256=data['hashes'],
        question='For EVERY source, how do frozen candidate pool cells relate to retained source-author labels?',
        rules='Explicit different CD8/Treg/Macrophage/B role is a disagreement flag, not automatic truth. Explicit CD4-only text contrary to CD8 is also flagged. Broad CD4 is compatible but insufficient to prove Treg. Missing/broad/unparsed labels remain unresolved, never a negative lineage.',
        scope='All7202990 selected metadata records and all61 source-arm entries; original172/atlas270 frozen annotations unchanged.',
        purpose='Review semantic compatibility before external expression; do not select a different channel based on downstream effects.',
        no_expression_no_P=True,no_annotations_changed=True)
    (out/'plan_before_author_frozen_role_crosscheck.json').write_text(json.dumps(plan,indent=2)+'\n')
    records=[];seen=0
    for arm in cell_arms(data):
        source=arm['source'];author=arm['author'];frozen=arm['roles'][0]
        ns=len(source.categories);na=len(author.categories)+1
        packed=(source.codes.astype(np.int64)*5+(frozen+1))*na+(author.codes+1)
        code,counts=np.unique(packed,return_counts=True)
        for value,count in zip(code,counts):
            src=int(value//(5*na));role=int((value//na)%5)-1;ac=int(value%na)-1
            text=str(author.categories[ac]) if ac>=0 else ''
            parsed=int(arm['author_mapping'][ac]) if ac>=0 else -1
            records.append(dict(arm=arm['arm'],dataset=str(source.categories[src]),
                frozen_role=VIEWS[role] if role>=0 else 'UNRESOLVED',author_label=text,
                author_parsed_role=VIEWS[parsed] if parsed>=0 else 'AMBIGUOUS' if parsed==-2 else 'UNRESOLVED',
                author_role_disagrees=author_disagrees(role,parsed,text),
                explicit_same_role=bool(role>=0 and parsed==role),cells=int(count)))
        seen+=len(source)
    tab=pd.DataFrame(records)
    assert seen==7202990 and tab.cells.sum()==seen
    tab.to_csv(out/'private_all_source_frozen_author_label_crosstabs.tsv.gz',sep='\t',index=False,compression='gzip')
    summary=[]
    for (arm,source),g in data['lib'].groupby(['arm','dataset'],sort=True):
        st=tab[(tab.arm==arm)&(tab.dataset==source)]
        assert st.cells.sum()==g.annotated_scope_cells.sum()
        for role in VIEWS:
            t=st[st.frozen_role==role];n=int(t.cells.sum())
            bad=int(t.loc[t.author_role_disagrees,'cells'].sum());same=int(t.loc[t.explicit_same_role,'cells'].sum())
            summary.append(dict(arm=arm,dataset=source,frozen_role=role,candidate_cells=n,
                source_author_explicit_same_cells=same,source_author_disagree_cells=bad,
                source_author_unresolved_or_other_broad_cells=n-same-bad,
                any_explicit_disagreement=bool(bad),annotation_replaced=False,external_model_accepted=False))
    assert len(summary)==244
    pd.DataFrame(summary).to_csv(out/'private_all244_source_view_author_concordance.tsv',sep='\t',index=False)
    terminal=dict(status='COMPLETE_ALL61_SOURCE_ROLE_AUTHOR_CROSSTABS_NO_IDENTITY_REPLACEMENT',
        finished_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),
        metadata_cells=seen,files={p.name:sha(p) for p in out.iterdir() if p.is_file()},new_expression_or_P=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'],finished_UTC=terminal['finished_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
