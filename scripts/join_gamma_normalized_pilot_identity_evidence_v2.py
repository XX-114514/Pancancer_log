"""Versioned metadata sidecar only; preserve every original identity column."""
import hashlib
import io
import json
import os
import shutil
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert out.parent.parent==Path('/dev/shm') and not out.exists() and shutil.disk_usage('/dev/shm').free>512*1024**2
    parent=Path('/dev/shm/pancancer-gamma-normalized-pilot-20261019-v2/private_output')
    audit=Path('/dev/shm/pancancer-gamma-normalized-audit-20261019-v2/private_output')
    assert sha(parent/'private_terminal.json')=='ef64454f3d5e9867f904c5cf50599a421670910a6bf6c39e1cb904a6abca5209'
    assert sha(audit/'private_terminal.json')=='666e90720403e7b764c3e3dd21cea7ac76eb54b97fbcd6d942412fcbf4d2fa0f'
    join=root/'gamma_reference_candidate_AB_evidence_join_20261019_v1'
    assert sha(join/'private_terminal.json')=='2d19ac4fa9dfdbb0f246d71e5cd9e8a5c72ae863fbc87ea84208ed821d60749f'
    for folder in [parent,join]:
        for name,d in json.loads((folder/'private_terminal.json').read_text())['files'].items():
            assert sha(folder/name)==d['sha256']
    old=root/'NSCLC_gamma_scGate_fixed_pilot_20261019_v1/private_full_input_and_output_transfer_proofs.json'
    new=parent/'private_full_input_and_output_transfer_proofs.json'
    assert json.loads(old.read_text())['binary_input_sha256']==json.loads(new.read_text())['binary_input_sha256']
    out.mkdir(mode=0o700);save=lambda name,value:(out/name).write_text(json.dumps(value,indent=2)+'\n')
    save('plan_before_normalized_metadata_sidecar.json',dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),original_join_terminal_sha256=sha(join/'private_terminal.json'),normalized_audit_sha256=sha(audit/'private_terminal.json'),codes={'0':'NOT_ASSESSED','1':'Pure','2':'Impure'},original_union_rows_and_columns_preserved=True,original_literal_HiTME_failure_slot_not_overwritten=True,no_new_expression_or_clinical_test=True,final_identity_adopted=False,per_cell_outputs_remain_server_tmpfs=True))
    original=dict(np.load(join/'private_all_original_union_rows_with_AB_evidence.npz',allow_pickle=False))
    pilot=dict(np.load(parent/'private_complete_pilot_candidate_classification_sidecar.npz',allow_pickle=False))
    rows=original['native_row'];pr=pilot['native_row'];names=pilot['models'].tolist();classes=pilot['classes']
    assert np.all(np.diff(rows)>0) and np.all(np.diff(pr)>0) and len(pr)==20000
    pos=np.searchsorted(rows,pr);assert (pos<len(rows)).all();np.testing.assert_array_equal(rows[pos],pr)
    normalized=classes[:,names.index('HiTME_full_symbol_normalized')]
    codes=np.zeros(len(rows),np.uint8);codes[pos]=np.where(normalized=='Pure',1,2)
    joined=dict(original);joined['normalized_HiTME_pilot_class_code']=codes
    blob=io.BytesIO();np.savez_compressed(blob,**joined)
    p=out/'private_all_original_union_with_normalized_HiTME_pilot_evidence.npz';p.write_bytes(blob.getvalue())
    reread=dict(np.load(p,allow_pickle=False))
    for name,values in original.items():np.testing.assert_array_equal(reread[name],values)
    lookup={int(row):(1 if label=='Pure' else 2) for row,label in zip(pr,normalized)}
    for i,row in enumerate(rows):assert int(reread['normalized_HiTME_pilot_class_code'][i])==lookup.get(int(row),0)
    assert np.count_nonzero(codes)==20000
    summary=[]
    for reference in [0,1,2]:
        for raw in [False,True]:
            for ab in [0,1,2]:
                for call in [0,1,2]:
                    take=(original['author_reference_code']==reference)&((original['raw_rule_bits']!=0)==raw)&(original['AB_annotation_status']==ab)&(codes==call)
                    summary.append(dict(source='GSE243013',author_reference_code=reference,any_old_raw_rule=raw,AB_annotation_status=ab,normalized_HiTME_pilot_class_code=call,cells=int(take.sum()),scope='ORIGINAL_EVIDENCE_UNION_NOT_ALL_NONREFERENCE_CELLS'))
    table=pd.DataFrame(summary);assert len(table)==54 and table.cells.sum()==len(rows)
    table.to_csv(out/'private_complete54_normalized_evidence_descriptions.tsv',sep='\t',index=False)
    same=normalized==classes[:,names.index('HiTME_native_subset')]
    save('private_normalized_model_equality_and_evidence_contract.json',dict(source='GSE243013',same_pilot_cells=20000,full_normalized_vs_native_subset_identical_calls=bool(same.all()),classification_disagreements=int((~same).sum()),unchanged_old_same_cell_binary_input=True,not_assessed_code_not_negative=True,old_literal_failure_slot_preserved=True,no_AB_record_not_gamma_truth=True,complete_AB_not_productivity_truth=True,formal_identity_changed=False))
    t=dict(status='COMPLETE_ALL_ORIGINAL_EVIDENCE_ROWS_PLUS_NORMALIZED_PILOT_FIELD_NO_IDENTITY_ADOPTION',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),original_columns_exact=True,new_field_every_row_dictionary_checked=True,pilot_metadata_rows=20000,source_description_positions=54,formal_identity_or_old_inference_changed=False,scientific_outputs_public=False,files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()})
    assert sum(d['bytes'] for d in t['files'].values())<1024**2
    save('private_terminal.json',t);print(json.dumps(dict(status=t['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':
    main()
