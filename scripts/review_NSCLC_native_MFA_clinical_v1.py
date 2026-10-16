"""All40 forests/all60 evidence and complete descriptive sensitivity review."""
import json
from datetime import datetime,timezone

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from NSCLC_native_MFA_clinical_inputs_v1 import read,matched_quota,sha,VIEWS


def main():
    root,pre,plan,scope,pairs,arrays=read();src=root/'NSCLC_native_MFA_clinical_results_20261016_v1';audit=src/'private_independent_clinical_result_audit.json'
    assert json.loads(audit.read_text())['status'].startswith('PASS_ALL40_')
    out=root/'NSCLC_native_MFA_clinical_review_20261016_v1';assert not out.exists();matched_quota(root);out.mkdir(mode=0o700)
    single=pd.read_csv(src/'private_all40_native_factor_clinical_endpoints.tsv',sep='\t')
    joint=pd.read_csv(src/'private_all60_native_factor_clinical_pairs.tsv',sep='\t')
    cr=root/'NSCLC_native_MFA_holdout_20261016_v2';ct=json.loads((cr/'private_terminal.json').read_text())
    cf=cr/'private_all60_internal_holdout_coordination.tsv';assert sha(cf)==ct['files'][cf.name]
    old=pd.read_csv(cf,sep='\t');keys=['factor','view_A','view_B']
    old=old[keys+['block_rank_rho','CI_low','CI_high','BH_q','BY_q','author_cases']].rename(columns={k:'coordination_'+k for k in ['block_rank_rho','CI_low','CI_high','BH_q','BY_q','author_cases']})
    evidence=joint.merge(old,on=keys,how='left',validate='one_to_one');assert len(evidence)==60 and evidence.coordination_block_rank_rho.notna().all()
    evidence['coordination_and_clinical_candidate']=evidence.robust_internal_coassociation&evidence.coordination_BY_q.lt(.05)&evidence.coordination_CI_low.gt(0)
    evidence.to_csv(out/'private_all60_coordination_and_clinical_evidence.tsv',sep='\t',index=False)
    factors=sorted(single.factor.unique());summaries=[]
    for factor in factors:
        s=single[single.factor==factor];p=evidence[evidence.factor==factor]
        summaries.append(dict(factor=factor,single_nominal=int(s.P.lt(.05).sum()),single_BH=int(s.BH_q.lt(.05).sum()),single_BY=int(s.BY_q.lt(.05).sum()),pair_BH=int(p.BH_q.lt(.05).sum()),pair_BY=int(p.BY_q.lt(.05).sum()),pair_min_BH=float(p.BH_q.min()),robust_coassociation_pairs=int(p.robust_internal_coassociation.sum()),coordination_and_clinical_pairs=int(p.coordination_and_clinical_candidate.sum()),
            min_valid_bootstrap=int(s.valid_bootstraps.min()),max_valid_bootstrap=int(s.valid_bootstraps.max()),
            clinical_beta_by_view='; '.join('%s:%.3f'%(r.view,r.beta) for r in s.itertuples())))
    summary=pd.DataFrame(summaries);summary.to_csv(out/'private_all10_native_factor_clinical_summary.tsv',sep='\t',index=False)
    labels={'author_CD8':'CD8','author_Treg':'Treg','author_Macrophage':'Macrophage','author_B_Plasma':'B / Plasma'}
    fig,axes=plt.subplots(1,4,figsize=(19,7),sharey=True)
    for ax,view in zip(axes,VIEWS):
        data=single[single.view==view].set_index('factor').loc[factors]
        ax.axvline(0,color='#777777',lw=.8)
        for j,r in enumerate(data.itertuples()):
            color='#9d1427' if r.BH_q<.05 else '#265580'
            if np.isfinite(r.CI_low) and np.isfinite(r.CI_high):ax.plot([r.CI_low,r.CI_high],[j,j],color=color,lw=1.5)
            ax.plot(r.beta,j,'o',color=color,ms=5)
        ax.set_title(labels[view]);ax.set_yticks(range(10),[x.rsplit('_',1)[1] for x in factors]);ax.set_xlabel('Response1 - response0\nDiscovery score SD units');ax.grid(axis='x',alpha=.18)
    axes[0].invert_yaxis()
    fig.suptitle('All40 native factor endpoint clinical associations / original internal holdout\nHC3 wild tests; 95% blocked case bootstrap CI. Red=auxiliary endpoint40 BH<.05; no direct interaction/prediction claim.',fontsize=11)
    fig.tight_layout(rect=[0,0,1,.9]);fig.savefig(out/'private_all40_clinical_forest.png',dpi=150,bbox_inches='tight');plt.close(fig)
    columns=list(dict.fromkeys(zip(joint.view_A,joint.view_B)));colnames=[labels[a]+' - '+labels[b] for a,b in columns]
    fig,axes=plt.subplots(1,2,figsize=(17,9))
    rho=np.array([[evidence[(evidence.factor==f)&(evidence.view_A==a)&(evidence.view_B==b)].iloc[0].coordination_block_rank_rho for a,b in columns] for f in factors])
    qvals=np.array([[evidence[(evidence.factor==f)&(evidence.view_A==a)&(evidence.view_B==b)].iloc[0].BH_q for a,b in columns] for f in factors])
    im=axes[0].imshow(rho,aspect='auto',cmap='RdBu_r',vmin=-1,vmax=1);fig.colorbar(im,ax=axes[0],fraction=.04);axes[0].set_title('Original within-block coordination /63cases')
    im=axes[1].imshow(-np.log10(np.maximum(qvals,1e-12)),aspect='auto',cmap='Blues',vmin=0,vmax=3);fig.colorbar(im,ax=axes[1],fraction=.04);axes[1].set_title('Clinical conjunction -log10(BH q) /62cases')
    for ax in axes:
        ax.set_yticks(range(10),[x.rsplit('_',1)[1] for x in factors]);ax.set_xticks(range(6),colnames,rotation=40,ha='right',fontsize=8)
    for i in range(10):
        for j in range(6):axes[0].text(j,i,'%.2f'%rho[i,j],ha='center',va='center',fontsize=8);axes[1].text(j,i,'%.2g'%qvals[i,j],ha='center',va='center',fontsize=8)
    fig.suptitle('All10 factors x all6 pairs: expression coordination and clinical coassociation are separate evidence\nClinical cell labels show BH q; score color scales differ. Complete families retained.',fontsize=12)
    fig.tight_layout(rect=[0,0,1,.92]);fig.savefig(out/'private_all60_multilevel_evidence.png',dpi=150,bbox_inches='tight');plt.close(fig)
    fields=['beta','mixture_residual_beta','agent_controlled_beta','discovery_descriptive_beta','all_frozen_four_view_cases_descriptive_beta']
    fig,axes=plt.subplots(2,2,figsize=(13,12))
    for ax,view in zip(axes.ravel(),VIEWS):
        values=single[single.view==view].set_index('factor').loc[factors,fields].to_numpy()
        im=ax.imshow(values,aspect='auto',cmap='RdBu_r',vmin=-1.5,vmax=1.5)
        ax.set_title(labels[view]);ax.set_yticks(range(10),[x.rsplit('_',1)[1] for x in factors]);ax.set_xticks(range(5),['Holdout primary','Holdout mixture residual','Holdout + agent blocks','Discovery descriptive','All complete descriptive'],rotation=25,ha='right',fontsize=8)
        for i in range(10):
            for j in range(5):ax.text(j,i,'%.2f'%values[i,j] if np.isfinite(values[i,j]) else 'NA',ha='center',va='center',fontsize=8)
    fig.suptitle('Complete descriptive coefficient sensitivity / identical discovery score units\nHoldout62, discovery140, complete202; no additional P or alternative primary family. Color clipped to+/-1.5.',fontsize=11)
    fig.tight_layout(rect=[0,0,.94,.92]);fig.colorbar(im,ax=axes.ravel().tolist(),fraction=.025,pad=.04);fig.savefig(out/'private_all40_clinical_sensitivity.png',dpi=150,bbox_inches='tight');plt.close(fig)
    report=['# Complete native factor clinical review','', 'Input/full null calibration/complete endpoint random streams independently audited.',
        'All10 factors/40 endpoints/60 pairs retained. Internal relative to the new native model, not an independent project cohort or pretreatment prediction.',
        'Pair conjunction concerns both marginal adjusted clinical associations, not a unique conditional contribution or direct LR mechanism.',
        '', '```text',summary.to_string(index=False),'```','',
        'If no robust coordinated clinical pair, close this native clinical family; preserve meaningful expression coordination and all uncertainty. Do not refit, change outcomes or add favorable subsets to rescue P. Next external/space qualification must use frozen programme and comparable lineages without exporting learned genes.']
    (out/'PRIVATE_COMPLETE_NATIVE_CLINICAL_REVIEW.md').write_text('\n'.join(report)+'\n')
    terminal=dict(status='COMPLETE_ALL40_ALL60_CLINICAL_AGGREGATE_REVIEW_VISUAL_REVIEW_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),
        audit_sha256=sha(audit),files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()},direct_interaction_or_independent_prediction_confirmed=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'],finished_UTC=terminal['finished_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
