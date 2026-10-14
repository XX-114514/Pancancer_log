"""Inspect frozen state candidates' metadata only, no learned gene export."""
import hashlib,json,os
from datetime import datetime,timezone
from pathlib import Path
from collections import Counter


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    directory=root/'state_programme_dictionary_20261014_v2'
    path=directory/'frozen_programme_dictionary.json'
    inventory=json.loads((directory/'inventory.json').read_text())
    assert hashlib.sha256(path.read_bytes()).hexdigest()==inventory['dictionary_sha256']
    dictionary=json.loads(path.read_text());records=dictionary['programmes']
    out=root/'frozen_T_sender_options_metadata_20261014_v1';assert not out.exists();out.mkdir()
    assert len(records)==98 and len({r['programme_id'] for r in records})==98
    meta=[]
    fields=['programme_id','method','allowed_majors','stability_main_candidate','foreign_RNA_sensitive','basic_expression_review_flag']
    for record in records:
        row={field:record.get(field) for field in fields}
        row['n_primary_signature_genes']=len(record['primary_genes'])
        meta.append(row)
    T=[r for r in meta if 'T cell' in r['allowed_majors']]
    report=dict(status='FROZEN_DICTIONARY_META_ONLY_AUDITED_NO_NEW_SCORES',finished_UTC=datetime.now(timezone.utc).isoformat(),
      dictionary_sha256=inventory['dictionary_sha256'],fields=fields,T_programmes=len(T),T_method_counts=dict(Counter(r['method'] for r in T)),
      T_stability_main_candidates=sum(bool(r['stability_main_candidate']) for r in T),
      T_foreign_RNA_flags=sum(bool(r['foreign_RNA_sensitive']) for r in T),
      T_expression_review_flags=sum(bool(r['basic_expression_review_flag']) for r in T),T_candidate_metadata=T,
      no_learned_genes_exported=True,no_formal_states=True)
    (out/'private_candidate_metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
    (out/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()
