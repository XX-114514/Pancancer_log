"""Independent signed-gene arithmetic, all clinical/function effects and LOO."""
import hashlib,json,os
from datetime import datetime,timezone
from pathlib import Path
import anndata as ad
import numpy as np,pandas as pd
from diagnose_frozen_MCP_two_end_functions_v1 import hnscc_inputs


def ranked(x):
    x=np.asarray(x,float)
    return np.array([np.count_nonzero(x<v)+(np.count_nonzero(x==v)+1)/2 for v in x],float)


def clinical_effect(x,y,C,cov):
    rx=ranked(x);ry=ranked(y)
    if min(rx.std(),ry.std())<1e-12:return np.nan
    rx=(rx-rx.mean())/rx.std();ry=(ry-ry.mean())/ry.std()
    extra=np.column_stack([ranked(cov[:,i]) for i in range(cov.shape[1])]) if cov.shape[1] else cov
    design=np.column_stack((C,extra,rx))
    if len(x)-design.shape[1]<10 or np.linalg.matrix_rank(design)!=design.shape[1]:return np.nan
    Q,R=np.linalg.qr(design,mode='reduced')
    return float(np.linalg.solve(R,Q.T@ry)[-1])


def function_effect(x,y,cov):
    N=np.column_stack((np.ones(len(x)),*[ranked(cov[:,i]) for i in range(cov.shape[1])]))
    if len(x)-N.shape[1]<10 or np.linalg.matrix_rank(N)!=N.shape[1]:return np.nan
    Q,_=np.linalg.qr(N,mode='reduced');a=ranked(x);b=ranked(y);a-=Q@(Q.T@a);b-=Q@(Q.T@b)
    if min(a.std(),b.std())<1e-12:return np.nan
    return float(a@b/(np.linalg.norm(a)*np.linalg.norm(b)))


def same(a,b):
    if not np.isfinite(a):assert pd.isna(b)
    else:assert np.isclose(a,b,atol=1e-8,rtol=2e-7),(a,b)


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);project=Path(os.environ['PANCANCER_PROJECT_ROOT'])
    out=root/'HNSCC_frozen_receiver_component_attribution_20261016_v1';audit=out/'private_independent_audit.json';assert not audit.exists()
    source=root/'HNSCC_within_macro_subtype_diagnostics_20261014_v1';crc=root/'CRC_fine_pair_DIALOGUE_20261014_v1'
    m=pd.read_csv(out/'private_patient_signed_and_quality_components.tsv',sep='\t',index_col=[0,1]).sort_index()
    original=pd.read_csv(source/'private_patient_subtype_scores.tsv',sep='\t').set_index(['library_id','subtype_candidate']).loc[m.index]
    selected=pd.read_csv(source/'private_selected_cells.tsv',sep='\t')
    g,source_rows,cells,var,eligible,measured=hnscc_inputs(root,project)
    where=cells.index.get_indexer(selected.cell);assert (where>=0).all()
    raw=ad.io.sparse_dataset(g['layers/counts'])[source_rows[where],:].tocsr();raw.eliminate_zeros()
    total=np.asarray(raw.sum(1)).ravel()
    weights=pd.read_csv(root/'HNSCC_frozen_MCP_gene_axis_20261014_v2/Macnode18_private_frozen_gene_columns.tsv',sep='\t')
    w=weights[weights.native_column.ge(0)]
    values=raw[:,w.native_column.to_numpy(int)].toarray();positive=np.zeros(len(where));negative=np.zeros(len(where))
    for j,r in enumerate(w.itertuples()):
        contribution=float(r.weight)*(np.log1p(values[:,j]*10000/total)-float(r.mean))/float(r.sd)
        if r.weight>0:positive+=contribution
        elif r.weight<0:negative+=contribution
    native=cells.iloc[where]
    expected=pd.DataFrame(dict(library_id=native.library_id.to_numpy(),subtype_candidate=selected.subtype_candidate.to_numpy(),positive=positive,negative=negative)).groupby(['library_id','subtype_candidate']).mean().loc[m.index]
    assert np.allclose(expected.positive,m.positive_weight_component,atol=1e-9,rtol=1e-9)
    assert np.allclose(expected.negative,m.negative_weight_component,atol=1e-9,rtol=1e-9)
    q=np.load(crc/'private_inputs/Macnode18_frozen_projection.npz')
    coefficients=pd.read_csv(crc/'private_transfer_fixed_weights_v1/Macnode18_MCP1_quality.tsv',sep='\t').set_index('term').value
    quality=(original[['mean_logRNA','mean_log_genes']].to_numpy()-q['quality_mean'])/q['quality_sd']
    correction=-coefficients['(Intercept)']-quality[:,0]*coefficients.cellQ-quality[:,1]*coefficients.geneQ
    assert np.allclose(correction,m.transferred_quality_correction,atol=1e-9,rtol=1e-9)
    assert np.allclose(m.positive_weight_component+m.negative_weight_component+correction,original.score,atol=1e-9,rtol=1e-9)
    clinical=pd.read_csv(out/'all45_private_component_clinical_descriptors.tsv',sep='\t')
    functions=pd.read_csv(out/'all1500_private_component_function_descriptors.tsv',sep='\t')
    means=pd.read_csv(root/'HNSCC_full50_orthogonal_function_diagnostics_20261014_v1/private_patient_subtype_function_means.tsv',sep='\t',index_col=[0,1])
    author=pd.read_csv(root/'HNSCC_author_clinical_table_user_20261014_v1/private_key_and_response_audit_v1/all57_private_author_table_join.tsv',sep='\t')
    author=author[author.time.eq('Pre')].set_index('patient_registry')
    cache={}
    for subtype in m.index.get_level_values(1).unique():
        md=original.xs(subtype,level=1).sort_index();components=m.xs(subtype,level=1).loc[md.index]
        fn=means.xs(subtype,level=1).loc[md.index];cm=author.loc[md.patient_id]
        C=np.column_stack((np.ones(len(md)),cm.neo_type.eq('D+T').to_numpy(float),cm['HPV/p16'].eq('positive').to_numpy(float)))
        qc=np.column_stack((md.mean_logRNA,md.mean_log_genes,np.log1p(md.native_subtype_cells)))
        cache[subtype]=(components,fn,cm,C,qc)
    for row in clinical.itertuples():
        components,fn,cm,C,qc=cache[row.subtype_candidate];x=components[row.component].to_numpy();y=cm['Tumor.regression'].to_numpy(float)
        if row.variant=='arm_HPV':cov=np.empty((len(x),0))
        elif row.variant=='arm_HPV_native_quality_count':cov=qc
        else:
            male=cm.Gender.str.strip().str.lower().map({'male':1,'m':1,'female':0,'f':0}).to_numpy(float)
            C=np.column_stack((C,male));cov=cm.Age.to_numpy(float)[:,None]
        same(clinical_effect(x,y,C,cov),row.beta);leave=[]
        for i in range(len(x)):
            k=np.arange(len(x))!=i;v=clinical_effect(x[k],y[k],C[k],cov[k])
            if np.isfinite(v):leave.append(v)
        assert len(leave)==row.LOO_valid
        if leave:same(min(leave),row.LOO_min);same(max(leave),row.LOO_max)
    for row in functions.itertuples():
        components,fn,cm,C,qc=cache[row.subtype_candidate];x=components[row.component].to_numpy();y=fn[row.reference].to_numpy()
        cov=qc if row.variant=='native_quality_count' else np.empty((len(x),0))
        same(function_effect(x,y,cov),row.rho);leave=[]
        for i in range(len(x)):
            k=np.arange(len(x))!=i;v=function_effect(x[k],y[k],cov[k])
            if np.isfinite(v):leave.append(v)
        assert len(leave)==row.LOO_valid
        if leave:same(min(leave),row.LOO_min);same(max(leave),row.LOO_max)
    assert len(clinical)==45 and len(functions)==1500
    report=dict(status='PASS_SIGNED_GENE_ARITHMETIC_ORIGINAL_SCORE_ALL1545_DESCRIPTORS_AND_COMPLETE_LOO',
        finished_UTC=datetime.now(timezone.utc).isoformat(),component_table_sha256=hashlib.sha256((out/'private_patient_signed_and_quality_components.tsv').read_bytes()).hexdigest(),
        new_P_BH_or_primary_model_changes=False)
    audit.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()
