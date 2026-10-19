"""Dependent all-call/transfer audit and complete source metadata descriptions."""
import hashlib
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    control=Path('/dev/shm/pancancer-gamma-normalized-pilot-20261019-v2')
    parent=control/'private_output'
    out=Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert out.parent.parent==Path('/dev/shm') and not out.exists()
    out.mkdir(mode=0o700)
    save=lambda name,value:(out/name).write_text(json.dumps(value,indent=2)+'\n')
    launch=json.loads((control/'launch_receipt.json').read_text())
    assert launch['source_sha256']=='3383c57a5c946946c1e359ce809d290a94b19c71256b3d96f9f90597ff04d95b'
    save('plan_before_new_classification_metadata_review.json',dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),parent_source_sha256=launch['source_sha256'],fixed_pilot_cells=20000,model_slots=4,all_author_types=51,all_AB_annotation_codes=[0,1,2],all_NK_overlap=[False,True],all_old_rule_overlap=[False,True],scope='Full source descriptions and exact prior-anchor/anonymous-transfer replay. No new RNA/graph rebuild, outcomes, P/FDR or identity adoption. No-record AB is not negative and complete same-row AB is not receptor productivity truth.'))
    deadline=time.monotonic()+7200
    while not (parent/'private_terminal.json').exists():
        path=Path('/proc')/str(launch['pid'])/'stat'
        assert path.exists() and time.monotonic()<deadline
        p=path.read_text().rsplit(')',1)[1].split()
        assert p[19]==str(launch['proc_start_ticks']) and p[0]!='Z'
        save('private_progress.json',dict(stage='WAITING_FOR_EXACT_ORIGINAL_NORMALIZED_PILOT',observed_UTC=datetime.now(timezone.utc).isoformat()))
        time.sleep(15)
    t=json.loads((parent/'private_terminal.json').read_text())
    assert t['status']=='COMPLETE_SAME_CELL_SYMBOL_NORMALIZED_FULL_SCGATE_AND_THREE_ANCHORS_REVIEW_PENDING'
    assert t['source_sha256']==launch['source_sha256']
    for name,d in t['files'].items():
        assert sha(parent/name)==d['sha256'] and (parent/name).stat().st_size==d['bytes']
    original=root/'NSCLC_gamma_scGate_fixed_pilot_20261019_v1'
    previous=dict(np.load(original/'private_complete_pilot_candidate_classification_sidecar.npz',allow_pickle=False))
    got=dict(np.load(parent/'private_complete_pilot_candidate_classification_sidecar.npz',allow_pickle=False))
    rows=got['native_row'];names=got['models'].tolist();classes=got['classes']
    assert names==['HiTME_full_symbol_normalized','PBMC_exact','HiTME_native_subset','NK_HiTME_exact']
    assert classes.shape==(20000,4) and np.isin(classes,['Pure','Impure']).all()
    expected=np.sort(np.random.default_rng(2026101901).choice(1254749,20000,replace=False))
    np.testing.assert_array_equal(rows,expected)
    np.testing.assert_array_equal(rows,previous['native_row'])
    for name in names[1:]:
        np.testing.assert_array_equal(classes[:,names.index(name)],previous['classes'][:,previous['models'].tolist().index(name)])
    # Reconstruct the exact anonymous R table from the recorded per-call labels.
    table=pd.DataFrame(classes,columns=names);table.insert(0,'position',np.arange(len(rows)))
    stream=table.to_csv(index=False,sep='\t',lineterminator='\n').encode()
    transfer=json.loads((parent/'private_full_input_and_output_transfer_proofs.json').read_text())
    assert hashlib.sha256(stream).hexdigest()==transfer['classification_TSV_sha256']
    assert transfer['raw_entries']==2010550708 and transfer['all_native_RNA_and_detection_exact']
    assert transfer['all_selected_old_hard_rule_bits_exact'] and transfer['three_unchanged_model_classifications_exactly_reproduced']
    assert transfer['anonymous_R_cells'] and not transfer['clinical_or_author_labels_supplied_to_R']
    join=root/'gamma_reference_candidate_AB_evidence_join_20261019_v1'
    assert sha(join/'private_terminal.json')=='2d19ac4fa9dfdbb0f246d71e5cd9e8a5c72ae863fbc87ea84208ed821d60749f'
    jt=json.loads((join/'private_terminal.json').read_text())
    side=join/'private_all_original_union_rows_with_AB_evidence.npz'
    assert sha(side)==jt['files'][side.name]['sha256']
    evidence=dict(np.load(side,allow_pickle=False))
    pos=np.searchsorted(evidence['native_row'],rows)
    assert (pos<len(evidence['native_row'])).all()
    np.testing.assert_array_equal(evidence['native_row'][pos],rows)
    status=evidence['AB_annotation_status'][pos]
    assert np.isin(status,[0,1,2]).all()
    oldrule=evidence['raw_rule_bits'][pos]!=0
    native=root/'NSCLC_full_native_gene_units_20261016_v1'
    unit=np.load(native/'private_cell_unit_index.npy',allow_pickle=False)
    u=pd.read_csv(native/'private_all_native_units_and_quality.tsv',sep='\t',keep_default_na=False).sort_values('unit_index')
    metadata=u.iloc[unit[rows]].reset_index(drop=True)
    labels=u[['major_cell_type','sub_cell_type']].drop_duplicates().sort_values(['major_cell_type','sub_cell_type'])
    assert len(labels)==51
    reported=pd.read_csv(parent/'private_complete_author_type_model_classification_summary.tsv',sep='\t',keep_default_na=False)
    full=[];joint=[]
    nk=classes[:,3]=='Pure'
    for kind in labels.itertuples():
        at=(metadata.sub_cell_type==kind.sub_cell_type).to_numpy()
        assert metadata.loc[at,'major_cell_type'].eq(kind.major_cell_type).all()
        for j,name in enumerate(names):
            pure=classes[:,j]=='Pure'
            row=dict(source='GSE243013',model=name,major_cell_type=kind.major_cell_type,sub_cell_type=kind.sub_cell_type,pilot_cells=int(at.sum()),pure_cells=int((at&pure).sum()),status='ASSESSED_FIXED_SOURCE_CALLS_NOT_IDENTITY_TRUTH' if at.any() else 'NO_PILOT_CELL_FOR_AUTHOR_REFERENCE')
            full.append(row)
            saved=reported[(reported.model==name)&(reported.sub_cell_type==kind.sub_cell_type)]
            if at.any():
                assert len(saved)==1 and int(saved.pilot_cells.iloc[0])==row['pilot_cells'] and int(saved.pure_cells.iloc[0])==row['pure_cells']
            else:assert saved.empty
            for ab in [0,1,2]:
                for conflict in [False,True]:
                    for raw in [False,True]:
                        selected=at&(status==ab)&(nk==conflict)&(oldrule==raw)
                        joint.append(dict(source='GSE243013',model=name,major_cell_type=kind.major_cell_type,sub_cell_type=kind.sub_cell_type,AB_annotation_status=ab,NK_model_pure=conflict,any_original_raw_rule=raw,pilot_cells=int(selected.sum()),model_pure_cells=int((selected&pure).sum())))
    full=pd.DataFrame(full);joint=pd.DataFrame(joint)
    assert len(full)==204 and len(joint)==2448
    assert full.groupby('model').pilot_cells.sum().eq(20000).all()
    assert joint.groupby('model').pilot_cells.sum().eq(20000).all()
    np.testing.assert_array_equal(full.groupby('model').pure_cells.sum().reindex(names),np.count_nonzero(classes=='Pure',axis=0))
    full.to_csv(out/'private_complete51_by4_author_reference_counts.tsv',sep='\t',index=False)
    joint.to_csv(out/'private_complete2448_author_AB_NK_rule_model_descriptions.tsv.gz',sep='\t',index=False,compression='gzip')
    receipt=dict(status='ACCEPTED_ALL80000_MODEL_CALLS_THREE_ANCHORS_ANONYMOUS_TRANSFER_AND_COMPLETE_SOURCE_DESCRIPTIONS',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),producer_terminal_sha256=sha(parent/'private_terminal.json'),model_call_positions=80000,original_anchor_call_positions=60000,complete_author_model_positions=204,complete_joint_description_positions=2448,independent_full_raw_or_graph_rebuild=False,TCR_productivity_or_external_identity_truth=False,new_clinical_or_P_FDR=False,formal_identity_changed=False,scientific_outputs_public=False,files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()})
    assert sum(v['bytes'] for v in receipt['files'].values())<2*1024**2
    save('private_terminal.json',receipt)
    print(json.dumps(dict(status=receipt['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':
    main()
