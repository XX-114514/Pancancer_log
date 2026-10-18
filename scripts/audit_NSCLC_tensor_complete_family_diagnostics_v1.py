"""Independent source scores/QR/full wild-t replay for all post-result diagnostics."""
import json
import os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import t as student_t
from scipy.cluster.hierarchy import linkage,fcluster
from scipy.spatial.distance import squareform
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC_immune_RNA_MCP_cohort_audit_math_v1 import prepare,evaluate
from NSCLC_tensor_followup_math_v1 import seed


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,2*1024**2);src=root/'NSCLC_tensor_complete_family_diagnostics_20261019_v1';out=root/'NSCLC_tensor_complete_family_diagnostics_audit_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    t=json.loads((src/'private_terminal.json').read_text());assert t['source_sha256']=='f4674b3b06cb3b5f64fce1623a233f2ca629ac029e736ad77bdcd6695a9bd6c1' and not t['original_main_results_changed'] and not t['original_bootstrap_CI_gate_changed']
    for n,r in t['files'].items():assert sha(src/n)==r['sha256']
    original=root/'NSCLC_communication_tensor_family_20261019_v1';read=lambda p:pd.read_csv(p,sep='\t');base=read(original/'private_all20_heldout_clinical_results.tsv').set_index('model_id');slots=base.index.tolist();axes=dict(np.load(original/'private_tensor_axes.npz'));cases=axes['heldout_ids'].tolist();states=axes['states'].tolist();available=axes['availability'][129:];columns=axes['native_gene_columns'];rev={int(c):i for i,c in enumerate(columns)};rms=float(axes['global_RMS']);native=root/'NSCLC_full_native_gene_units_20261016_v1';u=read(native/'private_all_native_units_and_quality.tsv');assert u.unit_index.tolist()==list(range(len(u)));index={(r.sampleID,r.sub_cell_type):r.unit_index for r in u.itertuples()};mm=np.load(native/'mean_logCP10K.npy',mmap_mode='r');M=np.zeros((87,38,len(columns)));quality=np.zeros((87,38,2));nn=np.zeros((87,38),int)
    for i,case in enumerate(cases):
        for j,state in enumerate(states):
            unit=index.get((case,state))
            if unit is not None:M[i,j]=mm[unit,columns];quality[i,j]=[u.iloc[unit].mean_logRNA,u.iloc[unit].mean_loggenes];nn[i,j]=u.iloc[unit].n_cells
    np.testing.assert_array_equal(nn>=20,available);del mm;vocab=json.loads((original/'private_selected_LR_vocabulary.json').read_text());total=np.zeros(87);mask=available[:,:,None]&available[:,None,:]
    for d in vocab:
        L=M[:,:,[rev[c] for c in d['native_ligand_columns']]].min(2).astype(np.float32);R=M[:,:,[rev[c] for c in d['native_receptor_columns']]].min(2).astype(np.float32);score=np.log1p(L.astype(float)[:,:,None]*R.astype(float)[:,None,:]).astype(np.float32);total+=np.sum(score.astype(float)*mask,axis=(1,2))
    burden=np.log1p(total/(len(vocab)*available.sum(1)**2)/rms);capture=available.sum(1).astype(float);stored=dict(np.load(src/'private_all87_source_background_inputs.npz'));assert stored['author_record_ids'].tolist()==cases;np.testing.assert_allclose(burden,stored['log1p_mean_observed_LRscore'],atol=1e-12,rtol=1e-12);np.testing.assert_array_equal(quality,stored['unit_meanlogRNA_loggenes']);np.testing.assert_array_equal(capture,stored['captured_author_types'])
    fs={};vectors=[[],[],[],[]]
    for slot in slots:
        k=int(slot[2:4]);j=int(slot[-2:])-1;data=dict(np.load(original/('private_selected_K%d_factors_and_projection.npz'%k)));f=[data[key][:,j] for key in ['context','LR','sender','receiver']];fs[slot]=f
        for a,v in enumerate(f):vectors[a].append(v/np.linalg.norm(v))
    cosine=[np.stack(v)@np.stack(v).T for v in vectors];mean=sum(cosine)/4;distance=np.maximum.reduce([1-mean,.5*(1-np.minimum.reduce(cosine[1:])),np.zeros((20,20))]);distance=(distance+distance.T)/2;np.fill_diagonal(distance,0);pairs=read(src/'private_all190_cross_rank_factor_pairs.tsv');assert len(pairs)==190
    for r in pairs.itertuples():
        i,j=slots.index(r.left),slots.index(r.right);np.testing.assert_allclose([r.context_cosine,r.LR_cosine,r.sender_cosine,r.receiver_cosine,r.mean4_cosine,r.distance],[*[c[i,j] for c in cosine],mean[i,j],distance[i,j]],atol=1e-12,rtol=1e-12)
    groups=fcluster(linkage(squareform(distance),method='complete'),t=.1,criterion='distance');members=sorted([sorted(slots[i] for i in range(20) if groups[i]==g) for g in set(groups)]);mapping={slot:'CPG%02d'%(i+1) for i,g in enumerate(members) for slot in g};roles=read(src/'private_all20_original_author_role_and_redundancy.tsv').set_index('model_id');meta=read(root/'NSCLC_communication_tensor_inputs_20261019_v1/private_all51_author_type_coverage.tsv').set_index('original_author_type')
    categories=[]
    for label in states:
        major=meta.loc[label,'major'].lower();categories.append('Treg' if label.startswith('CD4T_Treg_') else 'CD8' if label.startswith('CD8T_') else 'Myeloid' if 'myeloid' in major else 'B' if major.startswith('b') else 'Other_T_NK')
    for slot in slots:
        assert roles.loc[slot].descriptive_redundancy_group==mapping[slot]
        for kind,a in [('sender',fs[slot][2]),('receiver',fs[slot][3])]:
            totals={c:sum(v for g,v in zip(categories,a) if g==c) for c in set(categories)};key=max(totals,key=totals.get);category=key if totals[key]/sum(totals.values())>=.5 else 'Mixed';assert roles.loc[slot,kind+'_coarse_category']==category and roles.loc[slot,kind+'_highest_loading_original_author_label']==states[np.argmax(a)]
    intervals=read(src/'private_all20_supplementary_fixed_design_intervals.tsv').set_index('model_id');sensitivity=read(src/'private_all60_total_signal_and_endpoint_RNA_sensitivities.tsv').set_index(['model_id','sensitivity']);diagnostics={r['model_id']:r for r in json.loads((src/'private_all20_original_case_bootstrap_diagnostics.json').read_text())};pvalues=[];maxwild=0.
    for slot in slots:
        data=dict(np.load(original/('private_'+slot+'_clinical_inputs.npz')));x,y,C,names=data['x'],data['y'],data['C'],data['names'].tolist();model,status=prepare(x,y,C,names,'clinical');assert model is not None;g=evaluate(model);Q=model['q'];h=model['h'];fitted=Q@(Q.T@y);residual=(y-fitted)/(1-h);rng=np.random.default_rng(seed(slot,'full_model_wild_t'));values=np.empty(19999)
        for start in range(0,19999,1000):
            signs=np.where(rng.random((min(1000,19999-start),len(y)))<.5,-1.,1.).T;Y=fitted[:,None]+residual[:,None]*signs;error=Y-Q@(Q.T@Y);se=np.sqrt(np.sum((g['row'][:,None]*error/(1-h[:,None]))**2,axis=0));values[start:start+Y.shape[1]]=(g['row']@Y-g['beta'])/se
        old=np.load(src/('private_'+slot+'_fixed_design_wild_t.npz'))['values'];maxwild=max(maxwild,float(np.max(np.abs(values-old))));np.testing.assert_allclose(values,old,atol=1e-8,rtol=1e-8);qlo,qhi=np.quantile(values,[.025,.975]);ci=[g['beta']-qhi*g['se'],g['beta']-qlo*g['se']];critical=student_t.ppf(.975,model['df']);P=2*student_t.sf(abs(g['t']),model['df']);pvalues.append(P);r=intervals.loc[slot]
        assert r.original_percentile_CI_status=='UNAVAILABLE_ORIGINAL_VALID_DRAW_GATE_UNMET';np.testing.assert_allclose([r.effect,r.original_wild_P,r.original_BH20,r.original_BY20],[base.loc[slot].effect,base.loc[slot].wild_P,base.loc[slot].BH20,base.loc[slot].BY20],atol=1e-10,rtol=1e-10);np.testing.assert_allclose([r.supplementary_fixed_design_wild_t_CI_low,r.supplementary_fixed_design_wild_t_CI_high,r.supplementary_analytic_HC3_t_CI_low,r.supplementary_analytic_HC3_t_CI_high,r.secondary_analytic_HC3_t_P],[*ci,g['beta']-critical*g['se'],g['beta']+critical*g['se'],P],atol=1e-8,rtol=1e-8)
        U,V=fs[slot][2:];sq=np.array([[sum(U[j]*quality[i,j,q] for j in range(38) if available[i,j])/sum(U[available[i]]) for q in range(2)] for i in range(87)]);rq=np.array([[sum(V[j]*quality[i,j,q] for j in range(38) if available[i,j])/sum(V[available[i]]) for q in range(2)] for i in range(87)]);quality_extra=np.c_[sq,rq];total_extra=np.c_[burden,capture];quality_names=['sender_factor_meanlogRNA','sender_factor_meanloggenes','receiver_factor_meanlogRNA','receiver_factor_meanloggenes']
        for label,extra,en in [('total_LR_and_capture',total_extra,['total_LRscore','captured_type_count']),('endpoint_RNA_quality',quality_extra,quality_names),('both',np.c_[total_extra,quality_extra],['total_LRscore','captured_type_count']+quality_names)]:
            mm,status=prepare(x,y,np.c_[C,extra],names+en,'clinical');rr=sensitivity.loc[(slot,label)];assert rr.status==status
            if mm is not None:
                e=evaluate(mm);critical=student_t.ppf(.975,mm['df']);np.testing.assert_allclose([rr.effect,rr.HC3_t_CI_low,rr.HC3_t_CI_high],[e['beta'],e['beta']-critical*e['se'],e['beta']+critical*e['se']],atol=1e-8,rtol=1e-8);assert rr.nuisance_rank==mm['nuisance_rank']
        streams=dict(np.load(original/('private_'+slot+'_all_inference_streams.npz')));d=diagnostics[slot];keys,counts=np.unique(streams['bootstrap_status'],return_counts=True);assert d['original_bootstrap_reasons']==dict(zip(keys,counts.astype(int).tolist()));loo=streams['leave_one_author_effect'];assert d['leave_one_valid']==np.isfinite(loo).sum() and d['leave_one_same_sign']==(np.sign(loo[np.isfinite(loo)])==np.sign(g['beta'])).sum();np.testing.assert_allclose([d['leave_one_min'],d['leave_one_max']],[np.nanmin(loo),np.nanmax(loo)],atol=1e-12)
        Q0=model['Q'];h0=np.sum(Q0**2,axis=1);np.testing.assert_allclose(d['original_nuisance_max_leverage'],h0.max(),atol=1e-12);assert d['nuisance_records_over_half_leverage']==(h0>.5).sum()
        for r in d['minority_control_diagnostics']:
            flag=C[:,names.index(r['control'])];unique,count=np.unique(flag,return_counts=True);minority=flag==unique[np.argmin(count)];assert r['original_minority_records']==minority.sum();bad=streams['bootstrap_indices'][streams['bootstrap_status']!='ESTIMABLE'];size=minority[bad].sum(1);assert r['failed_draws_with_minority_absent']==(size==0).sum() and r['failed_draws_with_single_minority_record']==(size==1).sum()
        (out/'private_progress.json').write_text(json.dumps(dict(phase='ALL20_INDEPENDENT_POST_RESULT_REPLAY',model_id=slot,observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n')
    from statsmodels.stats.multitest import multipletests
    np.testing.assert_allclose(intervals.secondary_analytic_HC3_BH20,multipletests(pvalues,method='fdr_bh')[1],atol=1e-12);np.testing.assert_allclose(intervals.secondary_analytic_HC3_BY20,multipletests(pvalues,method='fdr_by')[1],atol=1e-12)
    mm,status=prepare(x,burden,C[:,:-2],names[:-2],'clinical');e=evaluate(mm);critical=student_t.ppf(.975,mm['df']);g=json.loads((src/'private_total_observed_LR_score_description.json').read_text());assert g['status']==status;np.testing.assert_allclose([g['effect'],g['HC3_t_CI_low'],g['HC3_t_CI_high']],[e['beta'],e['beta']-critical*e['se'],e['beta']+critical*e['se']],atol=1e-8,rtol=1e-8)
    receipt=dict(status='ACCEPTED_ALL20_SUPPLEMENTARY_WILD_T_INTERVALS_ALL190_PAIRS_60_BACKGROUNDS_AND_BOOTSTRAP_FAILURE_DIAGNOSTICS',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),producer_terminal_sha256=sha(src/'private_terminal.json'),resource_guard=guard,all_fixed_design_wild_statistics=399980,max_independent_wild_t_error=maxwild,original_main_or_case_CI_gate_changed=False,scientific_outputs_public=False)
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'status':receipt['status'],'terminal_sha256':sha(out/'private_terminal.json')}))


if __name__=='__main__':main()
