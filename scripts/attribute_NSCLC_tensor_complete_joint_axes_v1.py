"""Full joint LR x sender x receiver descriptive contrast partitions in RAM only."""
import json
import os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC_immune_RNA_MCP_cohort_math_v1 import design,estimate
from attribute_NSCLC_tensor_all20_clinical_contrasts_v1 import source
from NSCLC_tensor_contribution_math_v1 import log_activity_multiplier
from NSCLC_tensor_joint_contribution_math_v1 import gram_joint,qr_joint,fixture


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,8*1024**2);out=root/'NSCLC_tensor_joint_clinical_contributions_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    marginal=root/'NSCLC_tensor_clinical_contrast_attribution_20261019_v1'
    assert sha(root/'NSCLC_tensor_clinical_contrast_attribution_audit_20261019_v1/private_terminal.json')=='edf548e9cc0d3b24158725d78f34df203f779bfc999fadf3ad1a9a96ceedddca'
    save=lambda name,data:(out/name).write_text(json.dumps(data,indent=2)+'\n')
    save('plan_before_complete_joint_attribution.json',dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),math_sha256=sha(Path(__file__).with_name('NSCLC_tensor_joint_contribution_math_v1.py')),resource_guard=guard,scope='All original20 programmes, both original adjusted and unadjusted contrasts, all876LR x38x38 ordered author pairs. Post-result descriptive only; no new molecular P/FDR, no main analysis or endpoint changes.',method='Each joint term Y_lsr * [D_A inv(D_A.T D_A)]_lsr,k on original NNLS active face. Multiply original log-activity scaling and clinical contrast row. Complete Gram and independent explicit QR joint arrays recomputed and compared in RAM over every position, not just candidate terms. Full original LR and pair margins and effect sums anchored to accepted marginal analysis.',ranking='For all40 programme-scenario slots, both top100 absolute contributions and top100 contributions aligned to original total contrast, canonical flattened index breaks ties. Both positive and negative terms kept. Coverage and positive expression-score counts described, no single-axis hypothesis test.',storage='Complete joint arrays only RAM; bounded private ranking/diagnostics, original sample identifiers and learned features never public.8MiB maximum output and standing500GB+1GiB reserve.'))
    save('private_joint_synthetic_gate.json',fixture())
    src,axes,cases,states,available,vocab,lig,rec=source(root);stored=dict(np.load(marginal/'private_all20_complete_contribution_arrays.npz'));slots=stored['model_ids'].tolist();scenarios=stored['scenarios'].tolist();L=len(vocab);S=len(states);rms=float(axes['global_RMS']);registry=pd.read_csv(src/'private_all20_factor_registry.tsv',sep='\t').set_index('model_id');rows=np.zeros((2,20,len(cases)))
    for j,slot in enumerate(slots):
        d=dict(np.load(src/('private_'+slot+'_clinical_inputs.npz')));assert d['author_record_ids'].tolist()==cases;m,status=design(d['x'],d['y'],d['C'],d['names'],'clinical');assert m is not None;g,effect=estimate(m);rows[0,j]=g['inverse'][-1];x=d['x'];rows[1,j]=np.where(x==1,1/(x==1).sum(),-1/(x==0).sum());np.testing.assert_allclose(rows[:,j]@d['y'],stored['expected_contrast'][:,j],atol=1e-10,rtol=1e-10)
    coverage=available.astype(int).T@available.astype(int);positive=np.zeros((L,S,S),np.int16)
    for i in range(len(cases)):
        a=available[i];positive+=((lig[i]>0)[:,:,None]&(rec[i]>0)[:,None,:]&a[None,:,None]&a[None,None,:]).astype(np.int16)
    candidates=[];diagnostics=[];rank_receipts=[];maximum=0.
    for rank in range(2,7):
        bundle=dict(np.load(src/('private_selected_K%d_factors_and_projection.npz'%rank)));B,U,V=[np.asarray(bundle[k],float) for k in ['LR','sender','receiver']];indices=[slots.index('TC%02dF%02d'%(rank,j+1)) for j in range(rank)];med=np.array([registry.loc[slots[j]].positive_training_coefficient_median for j in indices]);sd=np.array([registry.loc[slots[j]].training_log_activity_sd for j in indices]);joint=np.zeros((2,rank,L,S,S));audit=np.zeros_like(joint)
        for i in range(len(cases)):
            a=available[i];Y=np.log1p(lig[i].astype(float)[:,:,None]*rec[i].astype(float)[:,None,:]).astype(np.float32).astype(float)/rms;Y[:,~a,:]=0;Y[:,:,~a]=0;c=bundle['heldout_coefficient'][i];g=gram_joint(Y,B,U,V,a,c);q=qr_joint(Y,B,U,V,a,c);np.testing.assert_allclose(g,q,atol=1e-7,rtol=1e-7);maximum=max(maximum,float(np.max(np.abs(g-q))));gamma=log_activity_multiplier(c,med,sd)
            for j,index in enumerate(indices):
                for sc in range(2):weight=rows[sc,index,i]*gamma[j];joint[sc,j]+=weight*g[j];audit[sc,j]+=weight*q[j]
            if i%10==0:save('private_progress.json',dict(phase='FULL_JOINT_GRAM_AND_EXPLICIT_QR_REPLAY',rank=rank,author_records_done=i+1,observed_UTC=datetime.now(timezone.utc).isoformat()))
        np.testing.assert_allclose(joint,audit,atol=1e-7,rtol=1e-7);np.testing.assert_allclose(joint.sum((3,4)),stored['LR'][:,indices],atol=1e-7,rtol=1e-7);np.testing.assert_allclose(joint.sum(2),stored['cellpair'][:,indices],atol=1e-7,rtol=1e-7);np.testing.assert_allclose(joint.sum((2,3,4)),stored['expected_contrast'][:,indices],atol=1e-7,rtol=1e-7)
        for sc,scenario in enumerate(scenarios):
            for j,index in enumerate(indices):
                slot=slots[index];values=joint[sc,j].ravel();effect=stored['expected_contrast'][sc,index];absolute=np.abs(values);mass=absolute/absolute.sum();aligned=np.sign(effect)*values;order_abs=np.argsort(-absolute,kind='stable')[:100];order_aligned=np.argsort(-aligned,kind='stable')[:100]
                diagnostics.append(dict(model_id=slot,scenario=scenario,total_effect=effect,sum_all_joint=float(values.sum()),joint_positions=len(values),effective_absolute_joint_terms=float(1/(mass@mass)),top1_absolute_mass=float(mass[order_abs[0]]),top100_absolute_mass=float(mass[order_abs].sum()),maximum_Gram_QR_difference=float(np.max(np.abs(joint[sc,j]-audit[sc,j])))))
                for kind,order in [('absolute',order_abs),('aligned_to_programme_effect',order_aligned)]:
                    for position,flat in enumerate(order):
                        l,s,r=np.unravel_index(flat,(L,S,S));d=vocab[l];candidates.append(dict(model_id=slot,scenario=scenario,ranking=kind,candidate_rank=position+1,canonical_resource_LR_id=d['canonical_resource_LR_id'],ligand_complex='_'.join(d['ligand_components']),receptor_complex='_'.join(d['receptor_components']),sender_original_author_type=states[s],receiver_original_author_type=states[r],observed_author_records=int(coverage[s,r]),nonzero_expression_score_author_records=int(positive[l,s,r]),clinical_contrast_contribution=float(values[flat]),absolute_mass_share=float(mass[flat]),status='POST_RESULT_JOINT_DESCRIPTIVE_CANDIDATE_NO_MOLECULAR_P_OR_FDR'))
        rank_receipts.append(dict(rank=rank,joint_positions_compared=int(joint.size),maximum_difference=float(np.max(np.abs(joint-audit)))))
        del joint,audit
    pd.DataFrame(candidates).to_csv(out/'private_all8000_joint_axis_candidate_terms.tsv.gz',sep='\t',index=False,compression='gzip');pd.DataFrame(diagnostics).to_csv(out/'private_all40_complete_joint_partition_diagnostics.tsv',sep='\t',index=False)
    assert len(candidates)==8000 and len(diagnostics)==40 and sum(r['joint_positions_compared'] for r in rank_receipts)==50597760
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()};assert sum(r['bytes'] for r in files.values())<8*1024**2
    terminal=dict(status='ACCEPTED_ALL50597760_JOINT_POSITIONS_GRAM_VS_EXPLICIT_QR_AND_ALL_MARGINS_8000_PRIVATE_CANDIDATE_ROWS',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),math_sha256=sha(Path(__file__).with_name('NSCLC_tensor_joint_contribution_math_v1.py')),original_marginal_audit_sha256='edf548e9cc0d3b24158725d78f34df203f779bfc999fadf3ad1a9a96ceedddca',rank_receipts=rank_receipts,maximum_per_record_joint_Gram_QR_difference=maximum,files=files,independent_implementations_executed_in_same_job=True,separate_process_joint_audit=False,main_inference_changed=False,new_molecular_P_values=False,scientific_outputs_public=False)
    save('private_terminal.json',terminal);print(json.dumps(dict(status=terminal['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
