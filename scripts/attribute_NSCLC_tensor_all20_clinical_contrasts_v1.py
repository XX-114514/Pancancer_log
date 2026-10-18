"""All20 complete LR/pair decomposition of original contrasts, no new P values."""
import json
import os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC_immune_RNA_MCP_cohort_math_v1 import design,estimate
from NSCLC_tensor_contribution_math_v1 import projection_parts,log_activity_multiplier


def fixture():
    rng=np.random.default_rng(2026101929);L,S,K,N=5,4,2,60;B=rng.uniform(.1,1,(L,K));B[-1]=0;U=rng.uniform(.1,1,(S,K));V=rng.uniform(.1,1,(S,K));wanted=rng.uniform(.1,1,(N,K));wanted[0]=0;wanted[1,1]=0;med=np.array([.5,.6]);sd=np.array([.7,.8]);lr=[];pairs=[];c=[];maximum=0.
    for i in range(N):
        available=np.ones(S,bool);available[i%S]=False;D=np.einsum('lk,sk,rk->lsrk',B,U,V);Y=np.einsum('lsrk,k->lsr',D,wanted[i]);Y[:,~available]=0;Y[:,:,~available]=0;mask=available[:,None]&available[None,:];matrix=D[:,mask,:].reshape(-1,K);y=Y[:,mask].ravel();coef,_=nnls(matrix,y);a,b=projection_parts(Y,B,U,V,available,coef);mult=log_activity_multiplier(coef,med,sd);lr.append(a*mult);pairs.append(b*mult);c.append(coef);maximum=max(maximum,float(np.max(np.abs(a.sum(0)-coef))))
    lr=np.asarray(lr);pairs=np.asarray(pairs);c=np.asarray(c);activity=np.log1p(c/med)/sd;np.testing.assert_allclose(lr.sum(1),activity,atol=1e-9);np.testing.assert_allclose(pairs.sum((1,2)),activity,atol=1e-9)
    x=np.r_[np.zeros(N//2),np.ones(N//2)];cov=rng.normal(size=(N,2))
    for j in range(K):
        y=activity[:,j]-.123;m,status=design(x,y,cov,['a','b'],'clinical');assert m is not None;g,beta=estimate(m);row=g['inverse'][-1];np.testing.assert_allclose(np.sum(np.einsum('n,nl->l',row,lr[:,:,j])),beta,atol=1e-10);np.testing.assert_allclose(np.sum(np.einsum('n,nsr->sr',row,pairs[:,:,:,j])),beta,atol=1e-10)
    return dict(status='PASS_SYNTHETIC_TRUE_ZERO_MISSING_ACTIVE_BOUNDARY_LOG_AND_COMPLETE_CONTRAST_ADDITIVITY',cases=N,components=K,maximum_raw_coefficient_partition_error=maximum)


def source(root):
    src=root/'NSCLC_communication_tensor_family_20261019_v1';axes=dict(np.load(src/'private_tensor_axes.npz'));cases=axes['heldout_ids'].tolist();states=axes['states'].tolist();available=axes['availability'][129:];vocab=json.loads((src/'private_selected_LR_vocabulary.json').read_text());native=root/'NSCLC_full_native_gene_units_20261016_v1';u=pd.read_csv(native/'private_all_native_units_and_quality.tsv',sep='\t');assert u.unit_index.tolist()==list(range(len(u)));lookup={(r.sampleID,r.sub_cell_type):r.unit_index for r in u.itertuples()};columns=axes['native_gene_columns'];reverse={int(c):j for j,c in enumerate(columns)};M=np.zeros((87,38,len(columns)));a=np.zeros((87,38),bool);mm=np.load(native/'mean_logCP10K.npy',mmap_mode='r')
    for i,case in enumerate(cases):
        for j,state in enumerate(states):
            unit=lookup.get((case,state))
            if unit is not None:M[i,j]=mm[unit,columns];a[i,j]=u.iloc[unit].n_cells>=20
    np.testing.assert_array_equal(a,available);del mm;lig=np.zeros((87,len(vocab),38),np.float32);rec=np.zeros_like(lig)
    for l,d in enumerate(vocab):lig[:,l]=M[:,:,[reverse[c] for c in d['native_ligand_columns']]].min(2);rec[:,l]=M[:,:,[reverse[c] for c in d['native_receptor_columns']]].min(2)
    return src,axes,cases,states,available,vocab,lig,rec


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,16*1024**2);out=root/'NSCLC_tensor_clinical_contrast_attribution_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    assert sha(root/'NSCLC_communication_tensor_family_audit_20261019_v2/private_terminal.json')=='e83389e75443b9168cd9e3d56a63980c14d42021ddfe1faedb3d3371c6d37101'
    save=lambda n,d:(out/n).write_text(json.dumps(d,indent=2)+'\n')
    plan={'created_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'math_sha256':sha(Path(__file__).with_name('NSCLC_tensor_contribution_math_v1.py')),'resource_guard':guard,'scope':'Every20original programme, all876LRterms and38x38ordered senderreceiverpairs, both originaladjusted andunadjusted clinicalcontrast. No programme/axis selected by newP. Results are post-result descriptive partition for prospective hypothesis generation, not binding/causal effects or additional molecular discoveries.','formula':'Observed NNLS activeface A: G_A=D_A.T@D_A; perLR partialrawcoefficient c_L=solve(G_A,D_L,A.T@y_L); analogously percellpair. Signedpartials sumtooriginalnonnegative coefficient. Multiplypartial by log1p(c/m)/(c*discoverySD), continuouszero limit1/(m*SD). Originalcontrast rowr=Dclinical_pinv[-1]; LRclinicalpart=sum_i r_i*scaledpartial_i. Intercept makes discoverycenterconstant cancel. FullLR andfullpair partitions must eachsum tooriginalclinicaleffect.','interpretation':'Signedterms canopposeeachother duejointdeconvolution; absoluteor alignedsharesare descriptive,not proportionsofcausality. Sharedligand/receptorexactcomplexfamilies aggregatefullterms; noindependent evidence counting. No molecular P/FDR inherited fromwholeprogramme. No rank/seed/model/endpointretraining or maininferencechange.','gate':'Syntheticobservedtruezero/missing/NNLSboundary/logtransform/contrastadditivity before realterms. Independentauditmustconstructexplicitobservedtemplates and activefaceQR operators for allcases/ranks, allterms/sums andlabels.','storage':'Only aggregateLR/pair effectarrays and completeprivatetables persisted; no cellmatrix/fulltensor/percaseattributionarray.16MiBcap and standingreserve maintained.'}
    save('plan_before_all20_descriptive_contribution_terms.json',plan);save('private_synthetic_partition_gate.json',fixture());src,axes,cases,states,available,vocab,lig,rec=source(root);L=len(vocab);S=len(states);rms=float(axes['global_RMS']);results=pd.read_csv(src/'private_all20_heldout_clinical_results.tsv',sep='\t').set_index('model_id');sens=pd.read_csv(src/'private_all60_background_sensitivity_effects.tsv',sep='\t').set_index(['model_id','sensitivity']);registry=pd.read_csv(src/'private_all20_factor_registry.tsv',sep='\t').set_index('model_id');slots=results.index.tolist();lr_terms=np.zeros((2,20,L));pair_terms=np.zeros((2,20,S,S));expected=np.zeros((2,20));operators=np.zeros((2,20,87));scenarios=['original_main','unadjusted'];max_raw_error=0.
    for index,slot in enumerate(slots):
        d=dict(np.load(src/('private_'+slot+'_clinical_inputs.npz')));assert d['author_record_ids'].tolist()==cases;m,status=design(d['x'],d['y'],d['C'],d['names'],'clinical');assert m is not None;g,effect=estimate(m);operators[0,index]=g['inverse'][-1];x=d['x'];operators[1,index]=np.where(x==1,1/(x==1).sum(),-1/(x==0).sum());expected[:,index]=[effect,sens.loc[(slot,'unadjusted')].effect];np.testing.assert_allclose(operators[:,index].sum(1),0,atol=1e-10);np.testing.assert_allclose(operators[:,index]@d['y'],expected[:,index],atol=1e-10,rtol=1e-10)
    for rank in range(2,7):
        bundle=dict(np.load(src/('private_selected_K%d_factors_and_projection.npz'%rank)));B,U,V=[bundle[k] for k in ['LR','sender','receiver']];indices=[slots.index('TC%02dF%02d'%(rank,j+1)) for j in range(rank)];med=np.array([registry.loc[slots[j]].positive_training_coefficient_median for j in indices]);sd=np.array([registry.loc[slots[j]].training_log_activity_sd for j in indices])
        for i in range(87):
            Y=np.log1p(lig[i].astype(float)[:,:,None]*rec[i].astype(float)[:,None,:]).astype(np.float32).astype(float)/rms;a=available[i];Y[:,~a,:]=0;Y[:,:,~a]=0;c=bundle['heldout_coefficient'][i];lp,pp=projection_parts(Y,B,U,V,a,c);max_raw_error=max(max_raw_error,float(np.max(np.abs(lp.sum(0)-c))));gamma=log_activity_multiplier(c,med,sd)
            for j,index in enumerate(indices):
                for scenario in range(2):
                    weight=operators[scenario,index,i]*gamma[j];lr_terms[scenario,index]+=weight*lp[:,j];pair_terms[scenario,index]+=weight*pp[:,:,j]
            if i%10==0:save('private_progress.json',dict(phase='ALL20_FULL_LR_AND_CELLPAIR_CONTRAST_PARTITION',rank=rank,author_records_done=i+1,observed_UTC=datetime.now(timezone.utc).isoformat()))
    np.testing.assert_allclose(lr_terms.sum(2),expected,atol=1e-7,rtol=1e-7);np.testing.assert_allclose(pair_terms.sum((2,3)),expected,atol=1e-7,rtol=1e-7);np.savez_compressed(out/'private_all20_complete_contribution_arrays.npz',LR=lr_terms,cellpair=pair_terms,expected_contrast=expected,scenarios=np.asarray(scenarios),model_ids=np.asarray(slots),LR_ids=np.asarray([d['canonical_resource_LR_id'] for d in vocab]),states=np.asarray(states),cellpair_observed_contexts=available.astype(int).T@available.astype(int))
    lr_rows=[];pair_rows=[];summary=[];families=[]
    for scenario,name in enumerate(scenarios):
        for index,slot in enumerate(slots):
            effect=expected[scenario,index];l=lr_terms[scenario,index];p=pair_terms[scenario,index];mass=np.abs(l);mass/=mass.sum();aligned=np.maximum(np.sign(effect)*l,0);aligned/=aligned.sum()
            summary.append(dict(model_id=slot,scenario=name,expected_contrast=effect,full_LR_sum=float(l.sum()),full_pair_sum=float(p.sum()),top1_absolute_contribution_mass=float(mass.max()),top5_absolute_contribution_mass=float(np.sort(mass)[-5:].sum()),effective_absolute_LR_count=float(1/np.sum(mass**2)),interpretation='descriptive_signed_partition_no_molecular_P_or_causal_claim'))
            for j,d in enumerate(vocab):lr_rows.append(dict(model_id=slot,scenario=name,canonical_resource_LR_id=d['canonical_resource_LR_id'],ligand_complex='_'.join(d['ligand_components']),receptor_complex='_'.join(d['receptor_components']),clinical_contrast_contribution=l[j],signed_net_share=l[j]/effect if effect else np.nan,absolute_mass_share=mass[j],aligned_mass_share=aligned[j]))
            for s in range(S):
                for t in range(S):pair_rows.append(dict(model_id=slot,scenario=name,sender_original_author_type=states[s],receiver_original_author_type=states[t],observed_author_records=int(np.sum(available[:,s]&available[:,t])),clinical_contrast_contribution=p[s,t],signed_net_share=p[s,t]/effect if effect else np.nan))
            for kind,key in [('ligand','ligand_components'),('receptor','receptor_components')]:
                groups={}
                for d,value in zip(vocab,l):groups.setdefault('_'.join(d[key]),[]).append(float(value))
                for group,values in groups.items():families.append(dict(model_id=slot,scenario=name,family_kind=kind,exact_complex=group,LR_terms=len(values),clinical_contrast_contribution=sum(values)))
    pd.DataFrame(lr_rows).to_csv(out/'private_all35040_LR_clinical_contribution_terms.tsv.gz',sep='\t',index=False,compression='gzip');pd.DataFrame(pair_rows).to_csv(out/'private_all57760_cellpair_clinical_contribution_terms.tsv.gz',sep='\t',index=False,compression='gzip');pd.DataFrame(families).to_csv(out/'private_all_exact_complex_family_contributions.tsv.gz',sep='\t',index=False,compression='gzip');pd.DataFrame(summary).to_csv(out/'private_all40_complete_partition_diagnostics.tsv',sep='\t',index=False)
    files={p.name:{'sha256':sha(p),'bytes':p.stat().st_size} for p in out.iterdir() if p.is_file()};assert sum(r['bytes'] for r in files.values())<16*1024**2;receipt=dict(status='COMPLETE_ALL20_FULL_LR_PAIR_AND_COMPLEX_FAMILY_CLINICAL_CONTRAST_PARTITIONS_INDEPENDENT_AUDIT_PENDING',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),math_sha256=sha(Path(__file__).with_name('NSCLC_tensor_contribution_math_v1.py')),files=files,LR_terms=35040,cellpair_terms=57760,maximum_raw_partition_error=max_raw_error,all_primary_results_changed=False,new_molecular_P_values=False,scientific_outputs_public=False)
    save('private_terminal.json',receipt);print(json.dumps({'status':receipt['status'],'terminal_sha256':sha(out/'private_terminal.json')}))


if __name__=='__main__':main()
