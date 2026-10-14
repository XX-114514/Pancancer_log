"""Read grouping/flag consistency before interpreting any association."""
import json,os
from datetime import datetime,timezone
from pathlib import Path
import pandas as pd


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    p=root/'HNSCC_frozen_T_sender_associations_20261014_v1/private_signature_dedup_contract.tsv'
    candidates=pd.read_csv(p,sep='\t')
    rows=[]
    for column in ['observed_signature_group','orthogonal_group']:
        conflicts=[]
        for group,frame in candidates.groupby(column):
            if frame.primary_admissible.any() and (frame.foreign_RNA_sensitive.any() or frame.basic_expression_review_flag.any()):
                conflicts.append(dict(group=group,programme_ids=sorted(frame.programme_id)))
        rows.append(dict(grouping=column,unique_groups=candidates[column].nunique(),admissible_flag_conflicts=conflicts))
    report=dict(status='SIGNATURE_FLAG_CONSISTENCY_CHECKED_WITHOUT_EFFECTS',observed_UTC=datetime.now(timezone.utc).isoformat(),groupings=rows)
    out=root/'HNSCC_T_sender_signature_flag_audit_20261014_v1';assert not out.exists();out.mkdir()
    (out/'private_flag_audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()
