"""Independent native detection and tie-aware manual prior-AP replay."""
import hashlib,json,os
from datetime import datetime,timezone
from pathlib import Path
import anndata as ad
import numpy as np,pandas as pd,pyreadr
from scipy.stats import rankdata
from diagnose_frozen_MCP_two_end_functions_v1 import hnscc_inputs


def manual_AP(y,prediction):
    order=np.argsort(-prediction,kind='stable');sorted_values=prediction[order]
    cuts=np.flatnonzero(np.r_[sorted_values[:-1]!=sorted_values[1:],True])
    hit=np.cumsum(y[order])[cuts];recall=hit/y.sum();precision=hit/(cuts+1)
    return float(np.diff(np.r_[0,recall])@precision)


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);project=Path(os.environ['PANCANCER_PROJECT_ROOT'])
    out=root/'HNSCC_receiver_prior_ligand_ranking_20261016_v1';audit_path=out/'private_independent_audit.json';assert not audit_path.exists()
    terminal=json.loads((out/'private_terminal.json').read_text());assert terminal['status']=='HNSCC_FIXED_RECEIVER_ALL_PRIOR_LIGANDS_RANKED_NO_ASSOCIATION_TESTS'
    ranking_path=out/'private_all_ligand_prior_rankings.tsv'
    assert hashlib.sha256(ranking_path.read_bytes()).hexdigest()==terminal['prior_ranking_sha256']
    plan=json.loads((out/'plan_before_native_readout_or_ranking.json').read_text())
    prior_path=root/'immune_same_name_family_source_review_v1/fullgene_receiver_target_readiness_v1/verified_prior_inputs_v2'
    objects=pyreadr.read_r(str(prior_path/'ligand_target_matrix_nsga2r_final.rds'));assert len(objects)==1
    target=next(iter(objects.values()));target.index=target.index.astype(str);target.columns=target.columns.astype(str)
    lr=next(iter(pyreadr.read_r(str(prior_path/'lr_network_human_21122021.rds')).values()))
    background=pd.read_csv(out/'private_native_prior_background_and_programme.tsv',sep='\t')
    assert background.gene.is_unique and background.native_column.is_unique
    g,source_rows,cells,var,eligible,measured=hnscc_inputs(root,project)
    source_receiver=pd.read_csv(root/'HNSCC_within_macro_subtype_diagnostics_20261014_v1/private_patient_subtype_scores.tsv',sep='\t')
    receiver=source_receiver[source_receiver.subtype_candidate.eq('Macro_C1QC')].set_index('patient_id')
    keep=cells.timepoint.eq('pre')&cells.patient_id.isin(receiver.index)&cells.library_id.isin(eligible.library_id)&cells.major.eq('Myeloid')&cells.subtype.eq('Macro_C1QC')&cells.evidence.eq('multigene_subtype_candidate')
    positions=np.flatnonzero(keep.to_numpy());assert len(positions)==terminal['native_receiver_cells']
    hits=np.zeros(len(var),np.int64);matrix=ad.io.sparse_dataset(g['layers/counts'])
    for start in range(0,len(positions),1000):
        x=matrix[source_rows[positions[start:start+1000]],:].tocsr();x.eliminate_zeros();assert x.has_canonical_format
        # Independently count native nonzero physical-column occurrences.
        hits+=np.bincount(x.indices,minlength=len(var))
    columns=background.native_column.to_numpy(int)
    assert np.array_equal(hits[columns],background.receiver_detected_cells)
    assert measured[columns].all() and (hits[columns]*10>=len(positions)).all()
    weights=pd.read_csv(root/'HNSCC_frozen_MCP_gene_axis_20261014_v2/Macnode18_private_frozen_gene_columns.tsv',sep='\t')
    positive=set(weights.loc[weights.weight.gt(0)&weights.native_column.ge(0),'native_column'].astype(int))
    negative=set(weights.loc[weights.weight.lt(0)&weights.native_column.ge(0),'native_column'].astype(int))
    assert not (positive&negative)
    expected=np.array(['positive' if c in positive else 'negative' if c in negative else 'background' for c in columns])
    assert np.array_equal(expected,background.programme_sign)
    rankings=pd.read_csv(ranking_path,sep='\t');B=target.loc[background.gene]
    assert not rankings.duplicated(['programme_sign','ligand']).any()
    for sign,group in rankings.groupby('programme_sign'):
        y=background.programme_sign.eq(sign).to_numpy(int);assert y.sum()>=10
        assert len(group)==len(target.columns) and set(group.ligand)==set(target.columns)
        for record in group.itertuples():
            v=B[record.ligand].to_numpy(float)
            ap=manual_AP(y,v);assert np.isclose(ap,record.average_precision,atol=1e-12,rtol=1e-12)
            assert np.isclose(ap-y.mean(),record.AP_excess_over_prevalence,atol=1e-12,rtol=1e-12)
            ranks=rankdata(v);n=int(y.sum());other=len(y)-n
            auc=(ranks[y==1].sum()-n*(n+1)/2)/(n*other)
            assert np.isclose(auc,record.AUROC,atol=1e-12,rtol=1e-12)
            if v.std()>0:
                a=v-v.mean();b=y-y.mean();pearson=float(a@b/(np.linalg.norm(a)*np.linalg.norm(b)))
                assert np.isclose(pearson,record.Pearson,atol=1e-12,rtol=1e-12)
        ordered=group.sort_values(['AP_excess_over_prevalence','ligand'],ascending=[False,True])
        assert np.array_equal(ordered.all_prior_rank,np.arange(1,len(ordered)+1))
        chosen=ordered[ordered.receiver_prior_compatible&ordered.AP_excess_over_prevalence.gt(0)].head(plan['top_ligands_per_sign'])
        assert set(group.loc[group.prospective_top20,'ligand'])==set(chosen.ligand)
    # Prior receptor strings must be literal network edges; no fuzzy symbol joins.
    compatibility=pd.read_csv(out/'private_all_ligand_receptor_compatibility.tsv',sep='\t',keep_default_na=False)
    edges=set(zip(lr['from'].astype(str),lr['to'].astype(str)))
    for record in compatibility.itertuples():
        for receptor in record.receiver_expressed_receptors.split('|') if record.receiver_expressed_receptors else []:
            assert (record.ligand,receptor) in edges
    report=dict(status='PASS_NATIVE_RECEIVER_DETECTION_ALL_PRIOR_AP_AUC_PEARSON_AND_FROZEN_TOP20',
        finished_UTC=datetime.now(timezone.utc).isoformat(),ranked_rows=len(rankings),
        prior_ranking_sha256=terminal['prior_ranking_sha256'],no_sender_receiver_effect_tests=True,
        metrics_are_custom_AP_not_official_NicheNet_activity=True,actual_complex_binding_not_verified=True)
    audit_path.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()
