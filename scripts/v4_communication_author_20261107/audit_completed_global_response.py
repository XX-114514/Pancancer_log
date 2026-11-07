"""Independent real-output coverage, expression/DE, enrichment and network readback."""
import json, zipfile, hashlib
from datetime import datetime, timezone
import numpy as np, pandas as pd, anndata as ad
from scipy.stats import ttest_ind_from_stats, hypergeom
from statsmodels.stats.multitest import multipletests

r=PROJECT_ROOT/'runs/20261107_V4_CoVarNet_communication_v1';out=r/'audit/completed_global_response_v1';out.mkdir(exist_ok=False)
proof=json.loads((r/'results/global/independent_audit.json').read_text());assert proof['status'].startswith('PASS_') and proof['iterations']==1000
pair=pd.read_csv(r/'results/global/audited_pair_counts.tsv',sep='\t');members=pd.read_csv(r/'audit/network_members.tsv',sep='\t')
assert len(pair)==45*45 and not pair.duplicated(['source','target']).any()
allcm=['CM%02d'%i for i in range(1,15)];cmrows=[]
for cm in allcm:
    nodes=set(members.loc[members.cm.eq(cm),'subCluster'])
    for other in allcm:
        othernodes=set(members.loc[members.cm.eq(other),'subCluster']);sel=pair[pair.source.isin(nodes)&pair.target.isin(othernodes)]
        cmrows.append({'source_CM':cm,'target_CM':other,'n_identity_pairs':len(sel),'mean_candidate_CCI':sel.n_candidate_CCI.mean() if len(nodes)*len(othernodes) else np.nan,'meaning':'global expression-specificity candidates, inclusive memberships; not CM environment or interaction validation'})
pd.DataFrame(cmrows).to_csv(out/'global_all14_CM_by_CM_candidate_CCI.tsv',sep='\t',index=False)
with zipfile.ZipFile(r/'sources/cellphonedb_v4.1.0.zip') as z:
    tables={name:pd.read_csv(z.open(name)) for name in ['gene_table.csv','protein_table.csv','complex_composition_table.csv','interaction_table.csv']}
join=tables['gene_table.csv'].merge(tables['protein_table.csv'],left_on='protein_id',right_on='id_protein')
protein_genes=join.groupby('protein_multidata_id').hgnc_symbol.apply(lambda x:set(x.dropna())).to_dict()
components={k:[k] for k in protein_genes}
for k,g in tables['complex_composition_table.csv'].groupby('complex_multidata_id'):components[k]=g.protein_multidata_id.astype(int).tolist()
a=ad.read_h5ad(r/'inputs/full_common_PRIVATE.h5ad');measured=set(a.var_names)
retained=pd.read_csv(next((r/'results/global').glob('statistical_analysis_pvalues_*.txt')),sep='\t');ids=set(retained.id_cp_interaction)
coverage=[]
for row in tables['interaction_table.csv'].itertuples():
    subunits=components.get(row.multidata_1_id,[])+components.get(row.multidata_2_id,[])
    alternatives=[protein_genes.get(k,set()) for k in subunits];allnames=set().union(*alternatives) if alternatives else set()
    good=bool(subunits and all(g and g&measured for g in alternatives))
    coverage.append({'id_cp_interaction':row.id_cp_interaction,'n_required_protein_components':len(subunits),
        'n_components_with_measured_HGNC':sum(bool(g&measured) for g in alternatives),
        'all_components_have_a_measured_symbol':good,'all_named_genes_measured_strict':bool(allnames and allnames<=measured),
        'n_components_with_multiple_possible_symbols':sum(len(g)>1 for g in alternatives),
        'required_gene_alternatives':';'.join('|'.join(sorted(g)) for g in alternatives),
        'retained_in_formal_global_CPDB':row.id_cp_interaction in ids})
coverage=pd.DataFrame(coverage);coverage.to_csv(out/'full_database_component_measurement.tsv',sep='\t',index=False)
response=pd.read_csv(r/'results/cytokine/response_full_family.tsv',sep='\t');tested=response[response.status.eq('tested')]
expected=hypergeom.sf(tested.overlap-1,tested.universe_size,tested.signature_size,tested.DE_size)
hypererr=float(np.max(np.abs(expected-tested.p)));assert hypererr<1e-12
for _,g in response.groupby('branch'):
    q=multipletests(g.p.fillna(1),method='fdr_bh')[1];valid=g.status.eq('tested');assert np.allclose(q[valid],g.loc[valid,'q_full_family'],atol=1e-12)
eligibility=pd.read_csv(r/'results/cytokine/CM_member_signature_eligibility.tsv',sep='\t')
maxp=0.;maxfc=0.;nDE=0;audited=[]
for i,row in eligibility.iterrows():
    if row.status!='READY_CELL_LEVEL_ONLY':continue
    b=a[a.obs.subCluster.eq(row.node)].copy();high=b.obs.CMT.eq(row.cm).to_numpy();x=b.X.astype(float)
    means=[];stds=[]
    for mask in [high,~high]:
        mat=x[mask];mean=np.asarray(mat.mean(0)).ravel();sq=mat.copy();sq.data**=2
        var=np.maximum((np.asarray(sq.mean(0)).ravel()-mean**2)*len(mask[mask])/(len(mask[mask])-1),0)
        means.append(mean);stds.append(np.sqrt(var))
    _,p=ttest_ind_from_stats(means[0],stds[0],int(high.sum()),means[1],stds[1],int((~high).sum()),equal_var=False)
    p=np.nan_to_num(p,nan=1);fc=np.log2((np.expm1(means[0])+1e-9)/(np.expm1(means[1])+1e-9))
    de=pd.read_csv(r/'results/cytokine'/(row.cm+'_'+str(i+1)+'_DE_PRIVATE.tsv.gz'),sep='\t').set_index('names').reindex(a.var_names)
    ep=float(np.max(np.abs(p-de.pvals)));ef=float(np.max(np.abs(fc-de.logfoldchanges)))
    assert ep<1e-6 and ef<1e-4,(row.cm,row.node,ep,ef)
    maxp=max(maxp,ep);maxfc=max(maxfc,ef);nDE+=len(de)
    audited.append({'cm':row.cm,'node':row.node,'n_genes':len(de),'n_high_cells':int(high.sum()),'n_low_cells':int((~high).sum()),'max_P_error':ep,'max_log2FC_error':ef})
pd.DataFrame(audited).to_csv(out/'direct_Scanpy_DE_independent_expression_audit.tsv',sep='\t',index=False)
edges=pd.read_csv(r/'results/cytokine/production_response_network.tsv',sep='\t')
for e in edges.itertuples():
    assert e.producer!=e.receiver and e.mean_expression>.1 and e.response_q<.05
    mask=a.obs.subCluster.eq(e.producer)&a.obs.CMT.eq(e.cm);mean=float(a[mask,e.subunit].X.mean())
    assert np.isclose(mean,e.mean_expression,atol=1e-5)
response.groupby('branch').agg(n_requested_tests=('p','size'),n_tested=('p','count'),n_member_q_lt05=('q_within_CM_member',lambda x:int(x.lt(.05).sum())),n_full_family_q_lt05=('q_full_family',lambda x:int(x.lt(.05).sum()))).reset_index().to_csv(out/'response_full_family_count_summary.tsv',sep='\t',index=False)
for name in ['CM_member_signature_eligibility.tsv','CM_member_DE_summary.tsv','response_full_family.tsv','CM_minimum_member_q_descriptive.tsv','production_response_network.tsv','cytokine_producer_symbol_mapping.tsv']:
    (out/name).write_bytes((r/'results/cytokine'/name).read_bytes())
result={'status':'PASS_independent_formal_global_and_direct_DE_enrichment_and_producer_expression_readback',
    'UTC':datetime.now(timezone.utc).isoformat(),'n_global_cells':38711,'n_global_labels':45,'iterations':1000,
    'n_candidate_cellpair_LR_events':int(pair.n_candidate_CCI.sum()),'n_database_interactions':len(coverage),
    'n_full_component_measurable_interactions':int(coverage.all_components_have_a_measured_symbol.sum()),
    'n_CPDB_retained_interactions':len(ids),'n_retained_with_incomplete_component_measurement':int((coverage.retained_in_formal_global_CPDB&~coverage.all_components_have_a_measured_symbol).sum()),
    'n_CM_member_DE_rows_audited':nDE,'n_CM_member_DE_models_audited':len(audited),'maximum_DE_P_error':maxp,'maximum_DE_log2FC_error':maxfc,
    'hypergeometric_maximum_error':hypererr,'n_production_response_subunit_edges_checked':len(edges),
    'full25_CMT_family_complete':False,'independent_patient_or_contact_or_causal_or_spatial_CM_validation_claim':False}
(out/'complete.json').open('x').write(json.dumps(result,indent=2));print(json.dumps(result))
