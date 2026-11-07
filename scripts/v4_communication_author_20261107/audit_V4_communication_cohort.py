"""All requested counts, comparisons, enrichment arithmetic and alternative model readback."""
import json
from datetime import datetime, timezone
import numpy as np,pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

r=PROJECT_ROOT/'runs/20261107_V4_CoVarNet_communication_v1'
terminal=json.loads((r/'launch/all_communication_terminal.json').read_text())
assert len(terminal['groups'])==25 and len({z['tag'] for z in terminal['groups']})==25
for z in terminal['groups']:
    assert z['fit_exit_code']==z['audit_exit_code']==0
    proof=json.loads((r/'results'/z['tag']/'independent_audit.json').read_text());assert proof['status'].startswith(('PASS_','UNASSESSABLE_'))
pair=pd.read_csv(r/'results/summary/all_original_CM_CCI_pair_comparisons.tsv',sep='\t')
tests=pd.read_csv(r/'results/summary/all14_CM_CCI_exploratory_tests.tsv',sep='\t');assert len(tests)==28
for row in pair.itertuples():
    for group,col in [('High','High'),('Low','Low_pooled')]:
        p=r/'results'/(row.cm+'_'+group)/'audited_pair_counts.tsv'
        if p.exists():
            table=pd.read_csv(p,sep='\t');selected=table[table.source.eq(row.source)&table.target.eq(row.target)]
            assert len(selected)==1 or pd.isna(getattr(row,col))
            if len(selected):assert float(selected.n_candidate_CCI.iloc[0])==getattr(row,col)
assert tests[tests.cm.isin(['CM01','CM02','CM05'])].n_paired_nonself_identity_pairs.eq(0).all()
response_status='FAILED_OR_UNAVAILABLE'
if terminal['cytokine_response_exit_code']==0:
    response=pd.read_csv(r/'results/cytokine/response_full_family.tsv',sep='\t');tested=response[response.status.eq('tested')]
    expected=hypergeom.sf(tested.overlap-1,tested.universe_size,tested.signature_size,tested.DE_size)
    assert np.allclose(expected,tested.p,atol=1e-12)
    for _,g in response.groupby('branch'):
        q=multipletests(g.p.fillna(1),method='fdr_bh')[1];valid=g.status.eq('tested');assert np.allclose(q[valid],g.loc[valid,'q_full_family'])
    edges=pd.read_csv(r/'results/cytokine/production_response_network.tsv',sep='\t')
    assert edges.mean_expression.gt(.1).all() and edges.response_q.lt(.05).all() and edges.producer.ne(edges.receiver).all()
    response_status='PASS_hypergeometric_complete_family_BH_and_subunit_threshold_readback'
cellchat_status='FAILED_OR_UNAVAILABLE'
if terminal['CellChat_exit_code']==0:
    p=r/'results/CellChat_global';edges=pd.read_csv(p/'candidate_LR_by_identity.tsv',sep='\t');counts=pd.read_csv(p/'identity_pair_counts.tsv',sep='\t')
    assert not counts.duplicated(['source','target']).any()
    actual=edges.groupby(['source','target']).size()
    for row in counts.itertuples():assert actual.get((row.source,row.target),0)==row.n_candidate_CCI
    cellchat_status='PASS_CellChat_saved_candidate_rows_and_pair_counts_readback'
receipt={'status':'PASS_all_original_CM_CPDB_family_and_available_secondary_outputs_independent_readback',
    'UTC':datetime.now(timezone.utc).isoformat(),'n_original_CM':14,'n_CPDB_comparison_groups':25,
    'n_CM_comparator_rows':len(tests),'n_cellpair_comparison_rows':len(pair),
    'response_status':response_status,'CellChat_status':cellchat_status,
    'all_complete_family_science_success':terminal['CellChat_exit_code']==terminal['cytokine_response_exit_code']==0,
    'spatial_gate_accepted_by_this_receipt':False,'patient_independence_or_intervention_or_mechanism_claim':False}
(r/'audit/FINAL_COMMUNICATION_COHORT_AUDIT.json').open('x').write(json.dumps(receipt,indent=2));print(json.dumps(receipt))
