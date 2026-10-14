"""Freeze all90 cNMF + eight stable Spectra candidates before projection."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib
import numpy as np,pandas as pd

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def frame(p):
    with np.load(p,allow_pickle=True) as z:return pd.DataFrame(z['data'],index=z['index'],columns=z['columns'])
def markers(x,genes):
    ids=np.argsort(-x,kind='stable');return genes[ids[x[ids]>0][:50]].tolist()
def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);P=Path(os.environ['PANCANCER_STATE_DICTIONARY_PLAN']);plan=json.loads(P.read_text());assert sha(__file__)==plan['producer_sha256'];C=R/'state_method_comparison_20261014_v2';S=R/'Spectra_multiseed_diagnostics_audit_20261014_v1';assert sha(C/'inventory.json')==plan['comparison_inventory_sha256'];assert sha(S/'all173_candidates_interpretation_safe.tsv')==plan['Spectra_catalogue_sha256']
    for parent in ['state_method_comparison_audit_20261014_v1','state_lineage_coexpression_audit_20261014_v1']:assert json.loads((R/parent/'independent_audit.json').read_text())['status'].startswith('PASS_')
    cn=pd.read_csv(C/'all90_cnmf_candidates.tsv',sep='\t');sp=pd.read_csv(S/'all173_candidates_interpretation_safe.tsv',sep='\t');panels=json.loads((Path('/tmp/state_lineage_coexpression_plan_v2.json')).read_text())['panels'];records=[]
    for row in cn.itertuples():
        name='T_cell' if row.major=='T cell' else 'Macrophage';b=R/'cNMF_state_contrast_20261014_v1'/name/'cnmf/cnmf_tmp';coef=frame(b/f'cnmf.gene_spectra_score.k_{row.K}.dt_0_5.df.npz');tpm=frame(b/f'cnmf.gene_spectra_tpm.k_{row.K}.dt_0_5.df.npz');genes=coef.columns.to_numpy(str)
        records.append(dict(programme_id=row.candidate_id,method='cNMF',allowed_majors=[row.major],K=row.K,component=row.component,primary_genes=markers(np.maximum(coef.loc[row.component].values,0),genes),TPM_sensitivity_genes=markers(tpm.loc[row.component].values,genes),stability_main_candidate=bool(row.stable_repetition_rule and row.crossK_supported),source_r2=float(row.source_r2),source_node_r2=float(row.source_node_r2),source_node_increment=float(row.source_node_r2-row.source_r2),RNA_partial_rank_after_source_node=float(row.RNA_partial_rank_after_source_node),detected_partial_rank_after_source_node=float(row.detected_partial_rank_after_source_node),gene_signature_origin='top50_positive_fullgene_TPM_regression_coefficients',crossK_independent_state_count=False))
    models=[]
    for seed in [20261014,20261015,20261016]:
        with np.load(R/'Spectra_state_discovery_20261014_v1/private_seed_models'/('seed'+str(seed))/'factors_and_discovery_scores.npz',allow_pickle=False) as z:models.append(dict(z))
    for row in sp[sp.stable_by_frozen_rule].itertuples():
        ix=[row.anchor_factor,row.seed2_factor,row.seed3_factor];weights=[]
        for model,j in zip(models,ix):
            v=model['factors'][j].astype(float);weights.append(v/np.linalg.norm(v));assert np.array_equal(model['genes'],models[0]['genes'])
        v=np.median(np.stack(weights),axis=0);records.append(dict(programme_id=row.candidate_id,method='Spectra',allowed_majors=['T cell','Macrophage'] if row.scope=='global' else [row.scope],primary_genes=markers(v,models[0]['genes']),TPM_sensitivity_genes=[],stability_main_candidate=True,scope=row.scope,matched_factor_indices=ix,source_r2=float(row.max_source_r2),source_node_r2=float(row.max_source_node_r2),source_node_increment=None,complexity_flag=bool(row.complexity_flag),best_prior_descriptor=row.best_prior_descriptor,prior_descriptor_is_confirmed_state=False,gene_signature_origin='top50_median_l2_normalized_weights_across_all3_matched_seeds'))
    for r in records:
        assert len(r['primary_genes'])==50 and len(set(r['primary_genes']))==50
        overlap={k:len(set(r['primary_genes'])&set(v)) for k,v in panels.items()};r['diagnostic_panel_overlaps']=overlap;r['descriptive_programmes']=[k for k,n in overlap.items() if n>=3];r['foreign_RNA_sensitive']=any(overlap[k]>=3 for k in ['epithelial_RNA','stromal_RNA','mural_RNA','neural_RNA']);r['ribosomal_top50_count']=sum(g.startswith(('RPL','RPS')) for g in r['primary_genes']);r['mitochondrial_top50_count']=sum(g.startswith('MT-') for g in r['primary_genes']);r['basic_expression_review_flag']=r['ribosomal_top50_count']>=10 or r['mitochondrial_top50_count']>=3;r['formal_state_name_assigned']=False;r['clinical_LR_p_used']=False;r['state_naming_status']='PENDING_NATIVE_CELL_SUPPORT_AND_IDENTITY_WITHIN_NODE_REVIEW';r['cross_method_agreement_assumed']=False
    assert len(records)==98 and sum(r['stability_main_candidate'] for r in records)==65
    O=R/'state_programme_dictionary_20261014_v1';assert not O.exists();O.mkdir();data=dict(schema='frozen_state_programme_dictionary_v1',frozen_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(P),input_comparison_inventory_sha256=sha(C/'inventory.json'),programme_count=98,stability_main_candidates=65,not_unique_states=True,programmes=records,formal_state_labels=False,UCell_projected=False);(O/'frozen_programme_dictionary.json').write_text(json.dumps(data,indent=2)+'\n')
    rows=[]
    for r in records:rows.append({k:(';'.join(v) if isinstance(v,list) else json.dumps(v) if isinstance(v,dict) else v) for k,v in r.items()})
    pd.DataFrame(rows).to_csv(O/'all98_programme_catalogue.tsv',sep='\t',index=False)
    inv=dict(status='FROZEN_ALL98_CANDIDATE_SIGNATURES_NOT_STATE_ACCEPTANCE',producer_sha256=sha(__file__),frozen_UTC=data['frozen_UTC'],programmes=98,main_stability_candidates=65,dictionary_sha256=sha(O/'frozen_programme_dictionary.json'),catalogue_sha256=sha(O/'all98_programme_catalogue.tsv'),formal_state_labels=False);(O/'inventory.json').write_text(json.dumps(inv,indent=2)+'\n');print(json.dumps(inv))
if __name__=='__main__':main()
