"""Explicit observed-template QR replay of every LR and senderreceiver term."""
import json
import math
import os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC_immune_RNA_MCP_cohort_audit_math_v1 import prepare,evaluate
from attribute_NSCLC_tensor_all20_clinical_contrasts_v1 import source


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,2*1024**2);out=root/'NSCLC_tensor_clinical_contrast_attribution_audit_20261019_v1';assert not out.exists();out.mkdir(mode=0o700);src=root/'NSCLC_tensor_clinical_contrast_attribution_20261019_v1';t=json.loads((src/'private_terminal.json').read_text());assert t['source_sha256']=='e6cc0693d8d5db8176c1fbc97ef608d107574cca72828a1606908e1e95a27e3f' and not t['all_primary_results_changed'] and not t['new_molecular_P_values']
    for n,r in t['files'].items():assert sha(src/n)==r['sha256']
    original,axes,cases,states,available,vocab,lig,rec=source(root);stored=dict(np.load(src/'private_all20_complete_contribution_arrays.npz'));slots=stored['model_ids'].tolist();assert len(slots)==20 and stored['scenarios'].tolist()==['original_main','unadjusted'];assert stored['states'].tolist()==states and stored['LR_ids'].tolist()==[d['canonical_resource_LR_id'] for d in vocab];L=len(vocab);S=len(states);rms=float(axes['global_RMS']);registry=pd.read_csv(original/'private_all20_factor_registry.tsv',sep='\t').set_index('model_id');LR=np.zeros((2,20,L));pairs=np.zeros((2,20,S,S));rows=np.zeros((2,20,87));expected=np.zeros((2,20));maximum=0.
    for index,slot in enumerate(slots):
        a=dict(np.load(original/('private_'+slot+'_clinical_inputs.npz')));assert a['author_record_ids'].tolist()==cases;m,status=prepare(a['x'],a['y'],a['C'],a['names'],'clinical');assert m is not None;g=evaluate(m);rows[0,index]=g['row'];x=a['x'];D=np.c_[np.ones(len(x)),x];q,R=np.linalg.qr(D,mode='reduced');target=np.array([0.,1.]);rows[1,index]=np.linalg.solve(R.T,target)@q.T;expected[:,index]=rows[:,index]@a['y']
    np.testing.assert_allclose(expected,stored['expected_contrast'],atol=1e-10,rtol=1e-10)
    for rank in range(2,7):
        f=dict(np.load(original/('private_selected_K%d_factors_and_projection.npz'%rank)));indices=[slots.index('TC%02dF%02d'%(rank,j+1)) for j in range(rank)];med=np.array([registry.loc[slots[j]].positive_training_coefficient_median for j in indices]);sd=np.array([registry.loc[slots[j]].training_log_activity_sd for j in indices])
        for i in range(87):
            a=available[i];ix=np.flatnonzero(a);Y=np.log1p(lig[i].astype(float)[:,:,None]*rec[i].astype(float)[:,None,:]).astype(np.float32).astype(float)/rms;Y=Y[:,a][:,:,a];coef=f['heldout_coefficient'][i];active=np.flatnonzero(coef>0);lp=np.zeros((L,rank));pp=np.zeros((S,S,rank))
            if len(active):
                D=np.einsum('lk,sk,rk->lsrk',f['LR'][:,active],f['sender'][a][:,active],f['receiver'][a][:,active]).reshape(-1,len(active));Q,R=np.linalg.qr(D,mode='reduced');operator=np.linalg.solve(R,Q.T);decomposed=operator.reshape(len(active),L,len(ix),len(ix))*Y[None,:,:,:];lp[:,active]=decomposed.sum((2,3)).T;local=decomposed.sum(1).transpose(1,2,0)
                for si,s in enumerate(ix):
                    for ti,receiver in enumerate(ix):pp[s,receiver,active]=local[si,ti]
            np.testing.assert_allclose(lp.sum(0),coef,atol=1e-7,rtol=1e-7);np.testing.assert_allclose(pp.sum((0,1)),coef,atol=1e-7,rtol=1e-7);maximum=max(maximum,float(np.max(np.abs(lp.sum(0)-coef))))
            for j,index in enumerate(indices):
                gamma=math.log1p(coef[j]/med[j])/(coef[j]*sd[j]) if coef[j]>0 else 1/(med[j]*sd[j])
                for sc in range(2):LR[sc,index]+=rows[sc,index,i]*gamma*lp[:,j];pairs[sc,index]+=rows[sc,index,i]*gamma*pp[:,:,j]
            if i%10==0:(out/'private_progress.json').write_text(json.dumps(dict(phase='EXPLICIT_ACTIVE_FACE_QR_ALL_LR_AND_PAIR_PARTITIONS',rank=rank,author_records_done=i+1,observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n')
    np.testing.assert_allclose(LR,stored['LR'],atol=1e-7,rtol=1e-7);np.testing.assert_allclose(pairs,stored['cellpair'],atol=1e-7,rtol=1e-7);np.testing.assert_allclose(LR.sum(2),expected,atol=1e-7,rtol=1e-7);np.testing.assert_allclose(pairs.sum((2,3)),expected,atol=1e-7,rtol=1e-7);np.testing.assert_array_equal(stored['cellpair_observed_contexts'],available.astype(int).T@available.astype(int))
    lrtable=pd.read_csv(src/'private_all35040_LR_clinical_contribution_terms.tsv.gz',sep='\t',keep_default_na=False);ptable=pd.read_csv(src/'private_all57760_cellpair_clinical_contribution_terms.tsv.gz',sep='\t',keep_default_na=False);families=pd.read_csv(src/'private_all_exact_complex_family_contributions.tsv.gz',sep='\t',keep_default_na=False);summary=pd.read_csv(src/'private_all40_complete_partition_diagnostics.tsv',sep='\t').set_index(['scenario','model_id']);assert len(lrtable)==35040 and len(ptable)==57760 and len(summary)==40
    for sc,scenario in enumerate(stored['scenarios'].tolist()):
        for index,slot in enumerate(slots):
            l=lrtable[(lrtable.scenario==scenario)&(lrtable.model_id==slot)];p=ptable[(ptable.scenario==scenario)&(ptable.model_id==slot)];assert len(l)==L and len(p)==S*S and l.canonical_resource_LR_id.tolist()==stored['LR_ids'].tolist();assert l.ligand_complex.tolist()==['_'.join(d['ligand_components']) for d in vocab] and l.receptor_complex.tolist()==['_'.join(d['receptor_components']) for d in vocab];assert p.sender_original_author_type.tolist()==[a for a in states for b in states] and p.receiver_original_author_type.tolist()==[b for a in states for b in states]
            np.testing.assert_allclose(l.clinical_contrast_contribution,LR[sc,index],atol=1e-7,rtol=1e-7);np.testing.assert_allclose(p.clinical_contrast_contribution,pairs[sc,index].ravel(),atol=1e-7,rtol=1e-7);np.testing.assert_array_equal(p.observed_author_records,stored['cellpair_observed_contexts'].ravel());effect=expected[sc,index];mass=np.abs(LR[sc,index]);mass/=mass.sum();aligned=np.maximum(np.sign(effect)*LR[sc,index],0);aligned/=aligned.sum();np.testing.assert_allclose(l.absolute_mass_share,mass,atol=1e-8);np.testing.assert_allclose(l.aligned_mass_share,aligned,atol=1e-8);np.testing.assert_allclose(l.signed_net_share,LR[sc,index]/effect,atol=1e-6,rtol=1e-7);np.testing.assert_allclose(p.signed_net_share,pairs[sc,index].ravel()/effect,atol=1e-6,rtol=1e-7)
            r=summary.loc[(scenario,slot)];np.testing.assert_allclose([r.expected_contrast,r.full_LR_sum,r.full_pair_sum,r.top1_absolute_contribution_mass,r.top5_absolute_contribution_mass,r.effective_absolute_LR_count],[effect,LR[sc,index].sum(),pairs[sc,index].sum(),mass.max(),np.sort(mass)[-5:].sum(),1/(mass@mass)],atol=1e-6,rtol=1e-7)
            for kind,col in [('ligand','ligand_complex'),('receptor','receptor_complex')]:
                grouped=l.groupby(col).clinical_contrast_contribution.agg(['sum','count']);ff=families[(families.scenario==scenario)&(families.model_id==slot)&(families.family_kind==kind)].set_index('exact_complex');assert set(ff.index)==set(grouped.index);np.testing.assert_array_equal(ff.loc[grouped.index].LR_terms,grouped['count']);np.testing.assert_allclose(ff.loc[grouped.index].clinical_contrast_contribution,grouped['sum'],atol=1e-10,rtol=1e-10)
    receipt=dict(status='ACCEPTED_ALL20_EXPLICIT_QR_35040_LR_57760_PAIR_AND_COMPLETE_COMPLEX_FAMILY_CONTRAST_PARTITIONS',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),producer_terminal_sha256=sha(src/'private_terminal.json'),resource_guard=guard,LR_terms=35040,cellpair_terms=57760,max_explicit_QR_raw_coefficient_error=maximum,all_original_inference_changed=False,new_molecular_P_values=False,scientific_outputs_public=False)
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'status':receipt['status'],'terminal_sha256':sha(out/'private_terminal.json')}))


if __name__=='__main__':main()
