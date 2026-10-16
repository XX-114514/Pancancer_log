"""Complete audited clinical landscape and exact-case cross-layer interpretation."""
import hashlib,json,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def truth(value):return str(value).lower() in {'true','1'}


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);source=r/'HNSCC_function_change_clinical_results_20261016_v1';audit=r/'HNSCC_function_change_clinical_audit_20261016_v1/private_terminal.json';rc=json.loads((r/'HNSCC_receiver_followup_job_control_20261016_v2/audit_HNSCC_function_change_clinical_results_v1/launch_receipt.json').read_text());deadline=time.monotonic()+43200
    while not audit.exists():
        p=Path('/proc')/str(rc['pid'])/'stat';assert p.exists();a=p.read_text().rsplit(')',1)[1].split();assert a[0]!='Z' and a[19]==str(rc['proc_start_ticks']) and time.monotonic()<deadline;time.sleep(5)
    proof=json.loads(audit.read_text());assert proof['status'].startswith('PASS_ALL_UNIQUE_MODELS');terminal=json.loads((source/'private_terminal.json').read_text());assert sha(source/'private_terminal.json')==proof['producer_terminal_sha256']
    for name,h in terminal['files'].items():assert sha(source/name)==h
    pre=r/'HNSCC_function_change_clinical_preflight_20261016_v1';original=pd.read_csv(pre/'private_all8000_clinical_scope.tsv',sep='\t',keep_default_na=False).set_index('slot');table=pd.read_csv(source/'private_all8000_longitudinal_function_clinical_results.tsv',sep='\t',keep_default_na=False,na_values=['NA']).set_index('slot');assert set(table.index)==set(original.index)
    for column in original.columns:pd.testing.assert_series_equal(table.loc[original.index,column],original[column],check_dtype=False,check_names=False)
    native=r/'HNSCC_native_multicell_functions_20261016_v1';masks=pd.read_csv(native/'private_all16_mask_axis.tsv',sep='\t');refs=pd.read_csv(native/'private_all50_reference_axis.tsv',sep='\t');mi={v:i for i,v in enumerate(masks.mask_id)};ri={v:i for i,v in enumerate(refs.reference)};archive=np.load(native/'private_all_multicell_function_scores_and_QA.npz',allow_pickle=False);function={k:archive[k] for k in archive.files};cases=pd.read_csv(r/'HNSCC_paired_multicell_function_changes_20261016_v2/private_all_slot_exact_patient_pair_differences.tsv.gz',sep='\t',dtype={'patient_block':str});groups={k:g.set_index('patient_block') for k,g in cases.groupby(['mask_id','reference','pool'],sort=False)}
    turn=pd.read_csv(pre/'private_all_patient_pool_fine_composition_turnover.tsv',sep='\t',dtype={'patient_block':str}).set_index(['patient_block','pool']);coord=pd.read_csv(r/'HNSCC_complete_multicell_function_review_20261016_v3/private_all4800_same_patient_change_coordination.tsv',sep='\t').set_index(['mask_id','reference','pool_A','pool_B']);designs={};model_rows=pd.read_csv(pre/'private_unique_numeric_designs.tsv',sep='\t')
    for row in model_rows.itertuples():
        a=np.load(pre/'private_designs'/(row.model_id+'.npz'),allow_pickle=False);d=a['D'];y=a['Y'];z=d[:,:-1];x=d[:,-1];ex=x-z@np.linalg.lstsq(z,x,rcond=None)[0];ey=y-z@np.linalg.lstsq(z,y,rcond=None)[0];raw=float(np.corrcoef(x,y)[0,1]);partial=float(np.corrcoef(ex,ey)[0,1]);designs[row.model_id]=dict(patients=a['patients'].astype(str).tolist(),raw_clinical_rank_rho=raw,quality_baseline_treatment_HPV_partial_rho=partial)
    out=r/'HNSCC_function_change_clinical_complete_review_20261016_v1';assert not out.exists();out.mkdir(mode=0o700);enriched=[]
    for slot,row in table.iterrows():
        item=row.to_dict();item['slot']=slot;endpoints=['A'] if row.family=='single3200' else ['A','B'];joint_CI=True;joint_LOO=True;people=None
        for label in endpoints:
            identifier=row[label+'_model'];model=designs.get(identifier)
            if model is None:joint_CI=False;joint_LOO=False;continue
            current=model['patients']
            if people is None:people=current
            else:assert people==current
            pool=row['pool_'+label];g=groups[(row.mask_id,row.reference,pool)];assert set(current)<=set(g.index);g=g.loc[current];i=mi[row.mask_id];j=ri[row.reference]
            item[label+'_raw_clinical_rank_rho']=model['raw_clinical_rank_rho'];item[label+'_complete_model_partial_rank_rho']=model['quality_baseline_treatment_HPV_partial_rho'];item[label+'_median_pool_composition_TV']=float(np.median([turn.loc[(person,pool),'composition_total_variation'] for person in current]));item[label+'_median_change_functional_RNA']=float(g.rounded_difference.median())
            for key in ['target_mean','nontarget_mean','expressed_target_genes']:item[label+'_same_patient_median_change_'+key]=float(np.median(function[key][i,g.post_pool_index.to_numpy(int),j]-function[key][i,g.pre_pool_index.to_numpy(int),j]))
            lo=float(row[label+'_CI_low']);hi=float(row[label+'_CI_high']);joint_CI &= np.isfinite([lo,hi]).all() and lo*hi>0;joint_LOO &= truth(row[label+'_LOO_same_sign']) and int(row[label+'_LOO_valid'])==len(current)
        item['all_endpoint_CI_exclude_zero']=bool(joint_CI);item['all_endpoint_LOO_valid_same_sign']=bool(joint_LOO);item['same_exact_patient_sets_confirmed']=people is not None and len(people)==int(row.clinical_complete_patients);item['new_direct_communication_or_independent_validation_accepted']=False
        if row.family=='pair4800':
            expected=sorted(set(groups[(row.mask_id,row.reference,row.pool_A)].index)&set(groups[(row.mask_id,row.reference,row.pool_B)].index));same_cases=people==expected;item['historical_coordination_same_exact_patients']=same_cases;item['clinical_endpoint_betas_same_sign']=bool(np.sign(row.A_beta)==np.sign(row.B_beta))
            if same_cases:
                h=coord.loc[(row.mask_id,row.reference,row.pool_A,row.pool_B)];assert h.same_patients==len(people)
                for field in ['raw_change_rho','quality12_adjusted_change_rho','quality12_status','B_function_only_adjusted_rho','B_only_status','quality12_plus_B_rho','quality12_plus_B_status']:item['matched_coordination_'+field]=h[field]
        if not str(row.association_status).startswith('ESTIMATED_'):category='INSUFFICIENT_CLINICAL_INFORMATION'
        elif row.BH_q>=.05:category='NO_FULL_FAMILY_SHARED_CLINICAL_SUPPORT' if len(endpoints)==2 else 'NO_FULL_FAMILY_SINGLE_POOL_CLINICAL_SUPPORT'
        elif row.BY_q<.05 and row.MC_upper_BH_q<.05 and joint_CI and joint_LOO:category='SOURCE_LIMITED_SHARED_CLINICAL_CANDIDATE' if len(endpoints)==2 else 'SOURCE_LIMITED_SINGLE_POOL_CLINICAL_CANDIDATE'
        else:category='EXPLORATORY_CLINICAL_ASSOCIATION_WITH_WEAK_OR_SENSITIVE_DIAGNOSTICS'
        item['evidence_category']=category;enriched.append(item)
    evidence=pd.DataFrame(enriched);assert len(evidence)==8000 and evidence.slot.is_unique;evidence.to_csv(out/'private_all8000_exact_case_multilevel_clinical_evidence.tsv',sep='\t',index=False,na_rep='NA');summaries={}
    for family,keys in [('pair4800',['reference','pool_A','pool_B']),('single3200',['reference','pool_A'])]:
        rows=[]
        for key,g in evidence[evidence.family.eq(family)].groupby(keys,sort=True):
            if not isinstance(key,tuple):key=(key,)
            assert len(g)==16;item=dict(zip(keys,key));item.update(mask_count=16,primary_BH_masks=int((g.BH_q<.05).sum()),BY_masks=int((g.BY_q<.05).sum()),MCupper_BH_masks=int((g.MC_upper_BH_q<.05).sum()),robust_source_masks=int(g.evidence_category.str.contains('SOURCE_LIMITED').sum()),min_BH_q=float(g.BH_q.min()),max_BH_q=float(g.BH_q.max()),min_patients=int(g.clinical_complete_patients.min()),max_patients=int(g.clinical_complete_patients.max()),independent_mask_replications=False)
            for label in ['A'] if family=='single3200' else ['A','B']:
                for field in ['beta','raw_clinical_rank_rho','complete_model_partial_rank_rho','median_pool_composition_TV']:item[label+'_median_'+field]=float(g[label+'_'+field].median())
            rows.append(item)
        summary=pd.DataFrame(rows);assert len(summary)==(300 if family=='pair4800' else 200);summary.to_csv(out/('private_all300_pair_function_summaries.tsv' if family=='pair4800' else 'private_all200_single_function_summaries.tsv'),sep='\t',index=False);summaries[family]=summary
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    pair=summaries['pair4800'];pairs=sorted(set(zip(pair.pool_A,pair.pool_B)));labels=lambda p:p.replace('candidate::','').replace('extension::','').replace('_plus_',' + ');order=refs.reference.tolist();names=[s.replace('HALLMARK_','').replace('_',' ') for s in order];matrix=[];ynames=[]
    for A,B in pairs:
        g=pair[pair.pool_A.eq(A)&pair.pool_B.eq(B)].set_index('reference').reindex(order)
        for label,pool in [('A',A),('B',B)]:matrix.append(g[label+'_median_complete_model_partial_rank_rho'].to_numpy());ynames.append(labels(A)+' / '+labels(B)+' | '+labels(pool))
    fig,ax=plt.subplots(figsize=(20,9),constrained_layout=True);im=ax.imshow(np.array(matrix),aspect='auto',vmin=-1,vmax=1,cmap='coolwarm');ax.set_yticks(range(12),ynames,fontsize=8);ax.set_xticks(range(50),names,rotation=70,ha='right',fontsize=7);ax.set_title('Complete HNSCC clinical landscape: rank residual correlations; each endpoint own full controls\nMedian across16 correlated exclusion masks; descriptive colour, inference uses audited HC3 wild tests');fig.colorbar(im,ax=ax,shrink=.6);fig.savefig(out/'private_all50_function_pair_endpoint_clinical_map.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(3,1,figsize=(20,12),constrained_layout=True)
    for ax,column,title in zip(axes,['primary_BH_masks','BY_masks','MCupper_BH_masks'],['Primary4800-family BH','BY dependence sensitivity','Simultaneous Monte Carlo upper-bound BH']):
        matrix=np.array([pair[pair.pool_A.eq(a)&pair.pool_B.eq(b)].set_index('reference').reindex(order)[column].to_numpy() for a,b in pairs]);im=ax.imshow(matrix,aspect='auto',vmin=0,vmax=16,cmap='Blues');ax.set_yticks(range(6),[labels(a)+' / '+labels(b) for a,b in pairs],fontsize=8);ax.set_xticks(range(50),names,rotation=70,ha='right',fontsize=7);ax.set_title(title+': number of16 sensitivity masks, not independent replications');fig.colorbar(im,ax=ax,shrink=.6)
    fig.savefig(out/'private_all50_function_pair_full_family_support.png',dpi=150);plt.close(fig)
    counts={family:{'slots':len(g),'estimated':int(g.association_status.str.startswith('ESTIMATED_').sum()),'primary_BH':int((g.BH_q<.05).sum()),'BY':int((g.BY_q<.05).sum()),'MCupper_BH':int((g.MC_upper_BH_q<.05).sum()),'source_candidates':int(g.evidence_category.str.contains('SOURCE_LIMITED').sum())} for family,g in evidence.groupby('family')}
    report='# Complete longitudinal functional clinical review\n\n'+json.dumps(counts,indent=2)+'\n\nAll8000 original scope fields, exact patient-key sets and all endpoints retained. Primary4800 and secondary3200 are separate families; gene exclusion masks are sensitivity analyses, not independent cohorts. Clinical association is exploratory, conditional on each endpoint own six quality changes, baseline function, treatment and HPV. Both endpoint marginal associations on the same patients do not establish unique incremental information or direct communication.\n\nComposition turnover is retained. The Treg pool has one frozen fine node, making its TV structurally zero; that is not evidence of unchanged cell state. CD8, myeloid and B/Plasma pool functional means can reflect subtype mixture and within-node changes. No state-only causal interpretation. Conditional bootstrap, small samples and unspecified age/site/processing effects limit transfer. Source programme clinical/LR families and original172/2026 October CM remain unchanged. All data, figures and clinical effects server private/current authorized chat only.\n'
    (out/'PRIVATE_COMPLETE_FUNCTION_CLINICAL_REVIEW.md').write_text(report);result=dict(status='COMPLETE_ALL8000_EXACT_SCOPE_CASE_EVIDENCE_ALL300_200_SUMMARIES_AND_COMPLETE_CLINICAL_MAPS_REVIEW_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),clinical_audit_sha256=sha(audit),files={p.name:sha(p) for p in out.iterdir() if p.is_file()},new_P_or_clinical_prediction_or_direct_communication=False)
    (out/'private_terminal.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
