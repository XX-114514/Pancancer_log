"""Complete same-context/endpoints evidence synthesis after full independent audit."""
import hashlib,json,os,time
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);result=r/'NSCLC_paired_programme_clinical_increment_results_20261016_v1';proof=result/'private_independent_audit.json'
    job=r/'HNSCC_receiver_followup_job_control_20261016_v2/audit_NSCLC_paired_clinical_increment_results_v1';receipt=json.loads((job/'launch_receipt.json').read_text());deadline=time.monotonic()+21600
    while not proof.exists():
        p=Path('/proc')/str(receipt['pid'])/'stat';assert p.exists();s=p.read_text().rsplit(')',1)[1].split()
        assert s[19]==str(receipt['proc_start_ticks']) and s[0]!='Z' and time.monotonic()<deadline;time.sleep(10)
    assert json.loads(proof.read_text())['status'].startswith('PASS_ALL572')
    t=pd.read_csv(result/'private_all572_paired_clinical_results.tsv',sep='\t');cpath=r/'NSCLC_two_end_coordination_results_20261016_v1/private_all572_two_end_results.tsv';c=pd.read_csv(cpath,sep='\t')
    clinical_path=r/'NSCLC_frozen_signature_clinical_results_20261016_v1/private_all_planned_signature_clinical_results.tsv';clinical=pd.read_csv(clinical_path,sep='\t')
    for sub in ['NSCLC_two_end_coordination_results_20261016_v1','NSCLC_frozen_signature_clinical_results_20261016_v1']:
        assert json.loads((r/sub/'private_independent_audit.json').read_text())['status'].startswith('PASS_')
    keys=['record','context','program','author_subtype_A','author_subtype_B'];assert t[keys].equals(c[keys]) and len(t)==572
    out=r/'NSCLC_paired_clinical_complete_evidence_review_20261016_v1';assert not out.exists();out.mkdir(mode=0o700)
    metadata=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),input_audit_sha256=sha(proof),coordination_source_sha256=sha(cpath),single_end_source_sha256=sha(clinical_path),
                  classification='Primary conjunction q<.05 plus both original bootstrap CIs excluding0 and both LOO direction-stable; then annotate existing same-case coordination, nonshared directions, MC and BY. These are not independent sources.',
                  display='Every planned row retained. Top12 lowest primary q for compact display only; all estimated context/programme forests retained. No new P or endpoint/family selection.',
                  limits='Post-treatment clinical association and multicellular coordination, not directed communication, synergistic causal action, baseline prediction or proof of an original CM axis.')
    (out/'private_review_plan.json').write_text(json.dumps(metadata,indent=2)+'\n')
    extras=c[['record','eligible_paired_cases','partial_rank_rho','nonshared_gene_partial_rho','BH_q_full572']].rename(columns={'partial_rank_rho':'old_pathology_controlled_coordination_rho','nonshared_gene_partial_rho':'old_nonshared_coordination_rho','BH_q_full572':'old_coordination_BH_q'})
    combined=t.merge(extras,on='record',validate='one_to_one');combined['old_single_end_eligibility_differs_possible']=True
    cr=clinical.set_index(['context','program','type','sub_cell_type']);rows=[]
    for row in combined.to_dict('records'):
        for typ,key in [('typeA','A'),('typeB','B')]:
            orig=cr.loc[(row['context'],row['program'],typ,row['author_subtype_'+key])]
            row['old_'+key+'_single_end_cases']=int(orig.eligible_cases);row['old_'+key+'_single_end_BH_q']=orig.BH_q_full_planned_family
        estimated=pd.notna(row['pair_IUT_P']);row['exact_case_join_verified']=False
        if estimated:
            a=np.load(r/'NSCLC_paired_programme_clinical_increment_preflight_20261016_v1/private_designs'/(row['record']+'.npz'),allow_pickle=False)
            b=np.load(r/'NSCLC_two_end_coordination_preflight_20261016_v1/private_pair_designs'/(row['record']+'.npz'),allow_pickle=False)
            assert np.array_equal(a['cases'],b['cases']) and np.array_equal(a['raw_A'],b['X']) and np.array_equal(a['raw_B'],b['Y'])
            assert np.allclose(a['C_raw'],np.delete(b['covariates'],8,axis=1)) and np.array_equal(a['Y'],b['response'])
            row['exact_case_join_verified']=True
        stable=estimated and row['BH_q_pair572']<.05 and bool(row['A_LOO_same_sign']) and bool(row['B_LOO_same_sign'])
        stable=stable and all(row[k+'_CI_low']>0 or row[k+'_CI_high']<0 for k in ['A','B'])
        row['both_conditional_clinical_supported']=bool(stable)
        row['coordination_supported_same_cases']=bool(pd.notna(row['old_coordination_BH_q']) and row['old_coordination_BH_q']<.05)
        row['conditional_sign_reversal_vs_matched_single_end']=bool(estimated and any(np.sign(row[k+'_beta'])!=np.sign(row[k+'_matched_single_end_beta']) for k in ['A','B']))
        row['nonshared_conditional_directions_retained']=bool(estimated and row['nonshared_conditional_status']=='DESCRIPTIVE_ONLY_NO_NEW_P' and all(np.sign(row['nonshared_'+k+'_beta'])==np.sign(row[k+'_beta']) for k in ['A','B']))
        row['nonshared_coordination_direction_retained']=bool(row['coordination_supported_same_cases'] and pd.notna(row['old_nonshared_coordination_rho']) and np.sign(row['old_nonshared_coordination_rho'])==np.sign(row['old_pathology_controlled_coordination_rho']))
        row['evidence_class']='INPUT_INSUFFICIENT' if not estimated else 'TWO_END_CONDITIONAL_CLINICAL_AND_COORDINATION' if stable and row['coordination_supported_same_cases'] else 'TWO_END_CONDITIONAL_CLINICAL_ONLY' if stable else 'COORDINATED_PROGRAMME_WITHOUT_VERIFIED_BOTH_END_INCREMENT' if row['coordination_supported_same_cases'] else 'NO_VERIFIED_PAIRED_CLINICAL_INCREMENT'
        rows.append(row)
    table=pd.DataFrame(rows);assert len(table)==572 and table.record.is_unique
    table.to_csv(out/'private_all572_complete_same_context_evidence.tsv',sep='\t',index=False,na_rep='NA')
    estimated=table[table.pair_IUT_P.notna()];display=estimated.sort_values(['BH_q_pair572','pair_IUT_P','record']).head(12)
    display.to_csv(out/'private_top12_display_only.tsv',sep='\t',index=False,na_rep='NA')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8})
    for (context,program),f in estimated.groupby(['context','program']):
        f=f.sort_values(['author_subtype_A','author_subtype_B']);y=np.arange(len(f));fig,axes=plt.subplots(1,2,figsize=(15,max(5,len(f)*.24+1.8)),sharey=True)
        for j,key in enumerate(['A','B']):
            ax=axes[j]
            for i,row in enumerate(f.itertuples()):
                color='#206b93' if row.both_conditional_clinical_supported else '#888888';lo=getattr(row,key+'_CI_low');hi=getattr(row,key+'_CI_high');point=getattr(row,key+'_beta')
                if np.isfinite(lo) and np.isfinite(hi):ax.plot([lo,hi],[i,i],lw=.8,color=color)
                ax.scatter(point,i,s=14,color=color)
            ax.axvline(0,color='#555555',lw=.7,ls='--');ax.set(title=key+' programme conditional on the other end',xlabel='Conditional MPR difference per original matched-case programme SD')
        axes[0].set_yticks(y,[f'{x.record} | {x.author_subtype_A} / {x.author_subtype_B}' for x in f.itertuples()],fontsize=7);axes[0].invert_yaxis()
        fig.suptitle(f'Context {context}, {program}: all estimable pairs; one post-treatment NSCLC source',fontsize=11)
        fig.subplots_adjust(left=.34,right=.98,bottom=.065,top=.94,wspace=.15);fig.savefig(out/f'private_complete_forest_C{context}_{program}.png',dpi=130);plt.close(fig)
    columns=['record','context','program','author_subtype_A','author_subtype_B','eligible_cases','A_beta','A_CI_low','A_CI_high','B_beta','B_CI_low','B_CI_high','BH_q_pair572','BY_q_pair572','MC_upper_BH_q_pair572','old_pathology_controlled_coordination_rho','old_coordination_BH_q','old_nonshared_coordination_rho','conditional_sign_reversal_vs_matched_single_end','nonshared_conditional_directions_retained','evidence_class']
    report=['# Complete paired clinical programme evidence','',
            'All572 planned pairs retained; every estimable pair has exactly the same cases, outcomes, raw programmes and clinical/technical covariates as its original coordination record. Single-end historical P values used more available cases and are not relabelled as matched-case tests.','',
            '## Evidence classes','',table.evidence_class.value_counts().to_string(),'','## Compact display (lowest primary conjunction q; not a new inferential family)','',display[columns].to_string(index=False),'',
            '## All estimable records','',estimated[columns].to_string(index=False),'','## Limits and next action','',
            'BOTH conditional information is a strict incremental question. Highly coordinated programmes may share a clinically meaningful signal without both retaining unique conditional information. Failure of this conjunction does not erase original single-end clinical or coordination results.','',
            'Conditional sign reversals may be suppression/contrast patterns, not co-activation or directional signalling. Nonshared effects are descriptive only; original MPR association is post-treatment, not baseline prediction. These context/programme/cell-pair records are not independent cohorts or separately validated original CM axes.','',
            'Prioritise source-limited programmes with coherent exact-endpoint clinical and nonshared coordination support; require independent patient/cohort and suitable spatial/protein/function evidence before a direct mechanism claim. No threshold or endpoint replacement within this family.']
    (out/'PRIVATE_COMPLETE_PAIRED_CLINICAL_REVIEW.md').write_text('\n'.join(report)+'\n')
    result=dict(status='COMPLETE_ALL572_SAME_CONTEXT_CASE_EVIDENCE_AND_ALL_ESTIMABLE_FORESTS_DISPLAY_AUDIT_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),
                files={p.name:sha(p) for p in out.iterdir() if p.is_file()},independent_cohort_or_direct_communication_claim=False)
    (out/'private_terminal.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
