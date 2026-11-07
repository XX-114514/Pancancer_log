"""Complete original CM family, both prose/code comparators; absent pairs stay NA."""
import json, itertools
from datetime import datetime, timezone
import numpy as np, pandas as pd
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests

r=PROJECT_ROOT/'runs/20261107_V4_CoVarNet_communication_v1';out=r/'results/summary';out.mkdir(exist_ok=False)
members=pd.read_csv(r/'audit/network_members.tsv',sep='\t');allcm=['CM%02d'%i for i in range(1,15)]
cache={}
def read(tag):
    if tag in cache:return cache[tag]
    p=r/'results'/tag/'audited_pair_counts.tsv'
    if not p.exists():return None
    a=pd.read_csv(p,sep='\t');cache[tag]=a.set_index(a.source+'|'+a.target).n_candidate_CCI;return cache[tag]
globalpair=read('global');g=pd.read_csv(r/'results/global/audited_pair_counts.tsv',sep='\t');cmglobal=[];comparison=[];tests=[]
for cm in allcm:
    nodes=members.loc[members.cm.eq(cm),'subCluster'].tolist()
    for other in allcm:
        othernodes=members.loc[members.cm.eq(other),'subCluster'].tolist();pairs=g[g.source.isin(nodes)&g.target.isin(othernodes)]
        cmglobal.append({'source_CM':cm,'target_CM':other,'n_evaluable_identity_pairs':len(pairs),
            'mean_candidate_CCI':pairs.n_candidate_CCI.mean() if len(nodes)*len(othernodes)>0 else np.nan,
            'membership_counting':'inclusive overlapping membership; cells/pairs are not independent patients'})
    for sender,receiver in itertools.product(nodes,repeat=2):
        pair=sender+'|'+receiver;high=read(cm+'_High');low=read(cm+'_Low');other_counts=[]
        for oc in allcm:
            if oc==cm:continue
            values=read(oc+'_High');other_counts.append(values.get(pair,np.nan) if values is not None else np.nan)
        hv=high.get(pair,np.nan) if high is not None else np.nan;lv=low.get(pair,np.nan) if low is not None else np.nan
        # Code-style comparator over available CMT pairs; missing CMTs never silently zero-filled.
        median=float(np.nanmedian(other_counts)) if np.isfinite(other_counts).any() else np.nan
        comparison.append({'cm':cm,'source':sender,'target':receiver,'self_pair':sender==receiver,'High':hv,'Low_pooled':lv,'Low_other_CMT_median':median,
            'n_other_CMT_pairs_available':int(np.isfinite(other_counts).sum()),'n_other_CMTs_requested':13,
            'missing_pairs_not_zero_filled':True})
pairframe=pd.DataFrame(comparison);pairframe.to_csv(out/'all_original_CM_CCI_pair_comparisons.tsv',sep='\t',index=False)
pd.DataFrame(cmglobal).to_csv(out/'global_CM_by_CM_mean_candidate_CCI.tsv',sep='\t',index=False)
for cm in allcm:
    selected=pairframe[pairframe.cm.eq(cm)&~pairframe.self_pair]
    for column in ['Low_pooled','Low_other_CMT_median']:
        eligible=selected[['High',column]].dropna();diff=eligible.High-eligible[column]
        p=float(wilcoxon(eligible.High,eligible[column],alternative='two-sided',zero_method='wilcox',method='approx',correction=True).pvalue) if len(diff) and (diff!=0).any() else 1. if len(diff) else np.nan
        tests.append({'cm':cm,'comparator':column,'n_paired_nonself_identity_pairs':len(eligible),
            'mean_High':eligible.High.mean(),'mean_Low':eligible[column].mean(),'median_delta_CCI':diff.median(),
            'p_paired_identity_pair_Wilcoxon_descriptive':p,'status':'EXPLORATORY_DEPENDENT_CELLPAIRS_SPATIAL_GATE_PENDING' if len(eligible) else 'EMPTY_OR_UNASSESSABLE_NOT_ZERO',
            'n_independent_patients':'not_certified'})
stats=pd.DataFrame(tests);stats['BH_all14_CM_diagnostic']=np.nan
for _,g in stats.groupby('comparator'):
    q=multipletests(g.p_paired_identity_pair_Wilcoxon_descriptive.fillna(1),method='fdr_bh')[1];stats.loc[g.index,'BH_all14_CM_diagnostic']=q
stats.loc[stats.p_paired_identity_pair_Wilcoxon_descriptive.isna(),'BH_all14_CM_diagnostic']=np.nan
stats.to_csv(out/'all14_CM_CCI_exploratory_tests.tsv',sep='\t',index=False)
terminal=json.loads((r/'launch/all_communication_terminal.json').read_text())
receipt={'status':'FULL_ORIGINAL14_CM_COMMUNICATION_SUMMARIES_WRITTEN_PENDING_COHORT_AUDIT',
    'UTC':datetime.now(timezone.utc).isoformat(),'n_original_CM':14,'n_CM_pair_cells':len(cmglobal),'n_identity_pair_rows':len(comparison),
    'n_CPDB_groups_requested':len(terminal['groups']),'CellChat_exit_code':terminal['CellChat_exit_code'],
    'cytokine_response_exit_code':terminal['cytokine_response_exit_code'],
    'full_original_empty_networks_in_denominator':True,'patient_contact_mechanism_or_publication_validation_claim':False}
(out/'complete.json').open('x').write(json.dumps(receipt,indent=2));print(json.dumps(receipt))
