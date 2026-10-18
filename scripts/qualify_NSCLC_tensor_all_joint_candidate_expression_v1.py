"""Post-result coverage for every joint candidate; does not retest or refit."""
import json
import os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,2*1024**2);out=root/'NSCLC_tensor_joint_candidate_expression_qualification_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    joint=root/'NSCLC_tensor_joint_clinical_contributions_20261019_v1';original=root/'NSCLC_communication_tensor_family_20261019_v1'
    save=lambda n,d:(out/n).write_text(json.dumps(d,indent=2)+'\n')
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,scope='Every pre-fixed8000candidate row, both clinical contrasts and both rankings. Eligibility description only, no per-axis P/FDR, no main score/missingness or candidate changes.',threshold='In an original heldout author record both original endpoint types have >=20cells; all ligand components each detected in >=10percent of sender cells and all receptor components each detected in >=10percent of receiver cells. Use original raw positive-cell-count cache and exact native component columns.',interpretation='Nonzero mean score is not this expression-proportion support. Record-count support is not design-estimability or independent-patient confirmation. Full original main inference unchanged. No phenotype used to choose candidates or expression cutoff; counts separately by original binary label only describe coverage.',storage='All component-detection and complete LR-pair support arrays RAM only; private8000row coverage sidecar plus full40slot summaries,2MiB output cap and standing500GB+1GiB reserve.')
    save('plan_before_all_candidate_detection_qualification.json',plan)
    t=json.loads((joint/'private_terminal.json').read_text());assert t['status'].startswith('ACCEPTED_ALL50597760') and not t['main_inference_changed'] and not t['new_molecular_P_values']
    for name,d in t['files'].items():assert sha(joint/name)==d['sha256']
    axes=dict(np.load(original/'private_tensor_axes.npz'));cases=axes['heldout_ids'].tolist();states=axes['states'].tolist();available=axes['availability'][129:];vocab=json.loads((original/'private_selected_LR_vocabulary.json').read_text());L=len(vocab);S=len(states);N=len(cases);assert (N,L,S)==(87,876,38)
    raw=root/'NSCLC_full_native_gene_units_20261016_v1';units=pd.read_csv(raw/'private_all_native_units_and_quality.tsv',sep='\t');assert units.unit_index.tolist()==list(range(len(units)));lookup={(r.sampleID,r.sub_cell_type):(r.unit_index,r.n_cells) for r in units.itertuples()};cols=axes['native_gene_columns'];reverse={int(c):i for i,c in enumerate(cols)};detect=np.zeros((N,S,len(cols)),bool);checked=np.zeros((N,S),bool);counts=np.load(raw/'positive_cell_counts.npy',mmap_mode='r')
    for i,case in enumerate(cases):
        for j,state in enumerate(states):
            row=lookup.get((case,state))
            if row is not None:
                unit,n=row;checked[i,j]=n>=20
                if n>=20:
                    values=counts[unit,cols];assert np.all(values>=0) and np.all(values<=n);detect[i,j]=values*10>=n
    np.testing.assert_array_equal(available,checked);del counts
    ligand=np.zeros((N,L,S),bool);receptor=np.zeros_like(ligand)
    for j,d in enumerate(vocab):
        ligand[:,j]=detect[:,:,[reverse[c] for c in d['native_ligand_columns']]].all(2)&available;receptor[:,j]=detect[:,:,[reverse[c] for c in d['native_receptor_columns']]].all(2)&available
    sample_x=None
    registry=pd.read_csv(original/'private_all20_factor_registry.tsv',sep='\t')
    for slot in registry.model_id:
        d=dict(np.load(original/('private_'+slot+'_clinical_inputs.npz')));assert d['author_record_ids'].tolist()==cases
        if sample_x is None:sample_x=d['x']
        else:np.testing.assert_array_equal(d['x'],sample_x)
    supported=np.einsum('nls,nlr->lsr',ligand.astype(np.int16),receptor.astype(np.int16),optimize=True)
    bylabel={v:np.einsum('nls,nlr->lsr',ligand[sample_x==v].astype(np.int16),receptor[sample_x==v].astype(np.int16),optimize=True) for v in [0,1]}
    np.testing.assert_array_equal(supported,bylabel[0]+bylabel[1]);observed=available.astype(int).T@available.astype(int);assert np.all(supported<=observed[None])
    table=pd.read_csv(joint/'private_all8000_joint_axis_candidate_terms.tsv.gz',sep='\t',keep_default_na=False);assert len(table)==8000;li={d['canonical_resource_LR_id']:i for i,d in enumerate(vocab)};si={s:i for i,s in enumerate(states)};records=[]
    for number,row in enumerate(table.itertuples()):
        l=li[row.canonical_resource_LR_id];s=si[row.sender_original_author_type];r=si[row.receiver_original_author_type];assert row.observed_author_records==observed[s,r] and supported[l,s,r]<=row.nonzero_expression_score_author_records
        records.append(dict(candidate_table_row=number,model_id=row.model_id,scenario=row.scenario,ranking=row.ranking,candidate_rank=row.candidate_rank,canonical_resource_LR_id=row.canonical_resource_LR_id,sender_original_author_type=row.sender_original_author_type,receiver_original_author_type=row.receiver_original_author_type,observed_author_records=int(observed[s,r]),all_complex_components_10percent_supported_author_records=int(supported[l,s,r]),supported_nonMPR_author_records=int(bylabel[0][l,s,r]),supported_MPR_author_records=int(bylabel[1][l,s,r]),status='DESCRIPTIVE_EXPRESSION_COVERAGE_NOT_AXIS_SIGNIFICANCE_OR_DESIGN_ESTIMABILITY'))
    result=pd.DataFrame(records)
    # Replay every distinct candidate directly from original integer counts.
    # This bypasses the full-array complex/minimum and einsum calculation.
    direct=np.load(raw/'positive_cell_counts.npy',mmap_mode='r');unique=result.drop_duplicates(['canonical_resource_LR_id','sender_original_author_type','receiver_original_author_type']);checked_axes=0
    for row in unique.itertuples():
        d=vocab[li[row.canonical_resource_LR_id]];flags=[]
        for case in cases:
            send=lookup.get((case,row.sender_original_author_type));receive=lookup.get((case,row.receiver_original_author_type));ok=False
            if send is not None and receive is not None and send[1]>=20 and receive[1]>=20:
                ok=all(int(direct[send[0],col])*10>=send[1] for col in d['native_ligand_columns']) and all(int(direct[receive[0],col])*10>=receive[1] for col in d['native_receptor_columns'])
            flags.append(ok)
        flags=np.asarray(flags,bool);assert int(flags.sum())==row.all_complex_components_10percent_supported_author_records;assert int(flags[sample_x==0].sum())==row.supported_nonMPR_author_records;assert int(flags[sample_x==1].sum())==row.supported_MPR_author_records;checked_axes+=1
    del direct
    result.to_csv(out/'private_all8000_joint_candidate_detection_coverage.tsv.gz',sep='\t',index=False,compression='gzip');summary=[]
    for keys,frame in result.groupby(['model_id','scenario','ranking'],sort=False):
        s=frame.all_complex_components_10percent_supported_author_records
        summary.append(dict(model_id=keys[0],scenario=keys[1],ranking=keys[2],candidate_rows=len(frame),support_under10=int((s<10).sum()),support10to19=int(((s>=10)&(s<20)).sum()),support20to42=int(((s>=20)&(s<43)).sum()),support43plus=int((s>=43).sum()),minimum_supported_author_records=int(s.min()),maximum_supported_author_records=int(s.max())))
    assert len(summary)==80;pd.DataFrame(summary).to_csv(out/'private_all80_candidate_ranking_coverage_summaries.tsv',sep='\t',index=False)
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()};assert sum(d['bytes'] for d in files.values())<2*1024**2
    receipt=dict(status='COMPLETE_ALL8000_JOINT_CANDIDATE_EXACT_NATIVE_DETECTION_COVERAGE_NO_NEW_INFERENCE',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),joint_terminal_sha256=sha(joint/'private_terminal.json'),files=files,candidate_rows=len(result),ranking_summaries=len(summary),unique_axes_direct_original_integer_count_replay=checked_axes,main_inference_changed=False,new_molecular_P_values=False,scientific_outputs_public=False)
    save('private_terminal.json',receipt);print(json.dumps(dict(status=receipt['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
