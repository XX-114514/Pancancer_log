"""Post-result description of both anchor programmes, all views/directions/sets."""
import gzip
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.stats import hypergeom

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


DIAGNOSTICS={
    'sex_linked_sentinels':['XIST','TSIX','JPX','DDX3Y','UTY','KDM5D','EIF1AY','ZFY','NLGN4Y','USP9Y','TMSB4Y'],
    'immediate_early_stress_associated':['FOS','FOSB','JUN','JUNB','JUND','ATF3','EGR1','DUSP1','HSPA1A','HSPA1B','DNAJB1','IER2','IER3','NR4A1','BTG2'],
    'erythroid_associated':['HBA1','HBA2','HBB','HBD','AHSP','ALAS2'],
    'epithelial_associated':['EPCAM','KRT8','KRT18','KRT19','KRT7','KRT5','KRT14','KRT17','MUC1'],
    'collagen_ECM_associated':['COL1A1','COL1A2','COL3A1','DCN','LUM','COL6A1','COL6A2']}


def qvalues(p):
    p=np.asarray(p,float);order=np.argsort(p,kind='stable');n=len(p)
    q=np.minimum.accumulate((p[order]*n/np.arange(1,n+1))[::-1])[::-1]
    result=np.empty(n);result[order]=np.minimum(q,1)
    return result


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget=resource(root,8*1024**2)
    audit=root/'immune_RNA_MCP_stability_coordination_audit_20261018_v1/private_terminal.json'
    assert sha(audit)=='27efdf31dfee1a3937de6790a03947ff13ab7af594f8a74ddd3e5ffeeaf11a6a'
    refs_path=Path(__file__).parent/'hallmark50_public_reference_v3'
    provenance=json.loads((refs_path/'reference_provenance.json').read_text());assert len(provenance['records'])==50
    refs={}
    for row in provenance['records']:
        p=refs_path/(row['name']+'.json');assert sha(p)==row['sha256']
        refs[row['name']]=set(str(g).strip().upper() for g in json.loads(p.read_text())[row['name']]['geneSymbols'])
    projection=root/'immune_RNA_MCP_frozen_projection_20261018_v1'
    terminal=json.loads((projection/'private_terminal.json').read_text())
    for name,digest in terminal['files'].items():assert sha(projection/name)==digest
    inputs=root/'immune_RNA_MCP_native_expression_inputs_20261018_v1/private_inputs'
    weights=pd.read_csv(projection/'private_R_projection/private_frozen_gene_weights.tsv',sep='\t')
    cases=pd.read_csv(projection/'private_R_projection/private_all_case_mean_scores.tsv',sep='\t',dtype={'sample':str})
    anchor='K2_seed20261018_source_adjusted'
    weights=weights[weights.fit_id.eq(anchor)].copy();cases=cases[cases.fit_id.eq(anchor)].copy()
    output=root/'immune_RNA_MCP_functional_description_20261018_v1'
    assert not output.exists();output.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),budget=budget,
        numerical_acceptance_sha256=sha(audit),public_reference_manifest_sha256=sha(refs_path/'reference_provenance.json'),
        result_followup='All6 anchor source tests have already been reviewed. This is descriptive interpretation of both programmes/all3 views, never new confirmatory evidence, gene re-selection or replacement of primary models.',
        annotation='All2 programmes x3 views x2 signed directions x50 public Hallmarks=600 slots. Exact candidate5000 per-view background, all positive NNLS signed coefficients, no top-N selection. Hypergeometric/BH600 describe annotation of selected genes, not clinical/causal validation.',
        diagnostics=DIAGNOSTICS,
        diagnostic_limits='Small public sentinel groups are not validated sex/purity/stress estimators. Their absence does not exclude confounding; stress/erythroid/epithelial/ECM RNA can be biological or technical. No automatic cell/identity or mechanism label.',
        contribution='For each sentinel group separately, recompute its signed RNA component using frozen discovery mean/SD and the same discovery-only nuisance projection. Group plus complement equals the original case score. Within-source covariance(group,total)/variance(total) is additive for that partition and may be outside[0,1]; groups overlap and cannot be added across categories. It is not causal explained variance.',
        privacy='All learned genes/weights, case components and scientific tables remain private. No new communication/clinical P family or dataset fit.')
    (output/'plan_before_complete_functional_description.json').write_text(json.dumps(plan,indent=2)+'\n')
    annotation=[];diagnostic=[];component_records=[];checks=[]
    controls=['cellQ','geneQ','sourceCC','fraction','logN']
    for view in ['R_FOXP3_capture','E_effector_capture','Myeloid']:
        parts={};genes=None
        for part in ['discovery','TNBC_holdout']:
            prefix=inputs/(view+'_'+part)
            md=pd.read_csv(str(prefix)+'_metadata.tsv',sep='\t',dtype={'cell':str,'sample':str})
            g=Path(str(prefix)+'_genes.tsv').read_text().splitlines()
            assert genes is None or genes==g;genes=g
            with gzip.open(str(prefix)+'_tpm.mtx.gz','rb') as handle:mat=mmread(handle).tocsr()
            assert mat.shape==(len(g),len(md));parts[part]=(md,mat)
        background=set(g.upper() for g in genes);assert len(background)==len(genes)==5000
        for program in ['MCP1','MCP2']:
            w=weights[weights.view.eq(view)&weights.program.eq(program)].copy()
            assert len(w)>0 and not w.gene.duplicated().any()
            w['symbol']=w.gene.str.upper();assert set(w.symbol)<=background
            for direction,sign in [('positive',1),('negative',-1)]:
                active=set(w.loc[w.weight*sign>0,'symbol'])
                for name,ref in sorted(refs.items()):
                    measured=ref&background;hit=active&measured;M=len(background);N=len(measured);n=len(active);k=len(hit)
                    status='ANNOTATION_AVAILABLE' if N and n else 'NO_EFFECTIVE_REFERENCE_OR_SIGNED_GENES'
                    p=float(hypergeom.sf(k-1,M,N,n)) if N and n else np.nan
                    annotation.append(dict(view=view,program=program,direction=direction,reference=name,
                        candidate_background=M,reference_in_background=N,signed_active_genes=n,overlap=k,
                        enrichment_fold=k*M/(n*N) if N and n else np.nan,annotation_p=p,status=status))
            index=pd.Index(genes).get_indexer(w.gene);assert (index>=0).all()
            md,mat=parts['discovery'];D=np.c_[np.ones(len(md)),md[controls].to_numpy(float)]
            assert np.linalg.matrix_rank(D)==6
            mu=w.discovery_mean.to_numpy();scale=w.discovery_sd.to_numpy();value=w.weight.to_numpy()/scale
            raw=np.asarray(value@mat[index]).ravel()-value@mu
            coef=np.linalg.lstsq(D,raw,rcond=None)[0]
            mdh,mh=parts['TNBC_holdout'];Dh=np.c_[np.ones(len(mdh)),mdh[['cellQ','geneQ','fraction','logN']].to_numpy(float)]
            ci=[0,1,2,4,5]
            total={'discovery':raw-D@coef,'TNBC_holdout':np.asarray(value@mh[index]).ravel()-value@mu-Dh@coef[ci]}
            for part,sc in total.items():
                meta=parts[part][0]
                got=pd.DataFrame(dict(sample=meta['sample'].to_numpy(),score=sc)).groupby('sample').score.mean()
                saved=cases[cases.view.eq(view)&cases.program.eq(program)&cases.partition.eq(part)].set_index('sample').score
                assert set(got.index)==set(saved.index)
                error=float(np.max(np.abs(got-saved.reindex(got.index))))
                assert error<1e-8
                checks.append(dict(view=view,program=program,partition=part,max_case_score_error=error))
            for label,symbols in DIAGNOSTICS.items():
                take=w.symbol.isin(symbols).to_numpy();v=value[take];i=index[take];m=mu[take]
                raw_group=np.asarray(v@mat[i]).ravel()-v@m if take.any() else np.zeros(len(md))
                coefficient=np.linalg.lstsq(D,raw_group,rcond=None)[0]
                groups={'discovery':raw_group-D@coefficient,
                    'TNBC_holdout':(np.asarray(v@mh[i]).ravel()-v@m if take.any() else np.zeros(len(mdh)))-Dh@coefficient[ci]}
                diagnostic.append(dict(view=view,program=program,diagnostic_class=label,
                    active_sentinel_genes=int(take.sum()),sentinels_in_candidate_background=len(set(symbols)&background),
                    absolute_coefficient_fraction=float(np.abs(w.weight.to_numpy()[take]).sum()/np.abs(w.weight.to_numpy()).sum())))
                for part,group in groups.items():
                    meta=parts[part][0]
                    frame=pd.DataFrame(dict(sample=meta['sample'].to_numpy(),dataset=meta.dataset.to_numpy(),group=group,total=total[part]))
                    mean=frame.groupby(['sample','dataset'])[['group','total']].mean().reset_index()
                    mean['complement']=mean.total-mean.group
                    for source,g in mean.groupby('dataset'):
                        a=g.group.to_numpy();b=g.total.to_numpy();rest=g.complement.to_numpy()
                        a=a-a.mean();b=b-b.mean();rest=rest-rest.mean();variance=float(b@b)
                        fraction=float(a@b/variance) if variance>1e-20 else np.nan
                        other=float(rest@b/variance) if variance>1e-20 else np.nan
                        if variance>1e-20:assert abs(fraction+other-1)<1e-10
                        corr=float(a@b/np.sqrt((a@a)*variance)) if (a@a)>1e-20 and variance>1e-20 else np.nan
                        component_records.append(dict(view=view,program=program,partition=part,dataset=source,
                            diagnostic_class=label,case_or_sample_units=len(g),component_total_correlation=corr,
                            group_total_covariance_fraction=fraction,complement_total_covariance_fraction=other,
                            interpretation='Descriptive partition; neither causal fraction nor independent validation'))
    table=pd.DataFrame(annotation);assert len(table)==600
    table['annotation_BH600']=qvalues(table.annotation_p.fillna(1));table.loc[table.annotation_p.isna(),'annotation_BH600']=np.nan
    table.to_csv(output/'private_all600_Hallmark_signed_annotations.tsv',sep='\t',index=False)
    pd.DataFrame(diagnostic).to_csv(output/'private_all30_sentinel_gene_content_descriptions.tsv',sep='\t',index=False)
    pd.DataFrame(component_records).to_csv(output/'private_all90_source_sentinel_RNA_component_descriptions.tsv',sep='\t',index=False)
    pd.DataFrame(checks).to_csv(output/'private_all12_original_case_score_reconstruction.tsv',sep='\t',index=False)
    assert len(diagnostic)==30 and len(component_records)==90 and len(checks)==12
    files={p.name:sha(p) for p in sorted(output.iterdir()) if p.is_file()}
    (output/'private_terminal.json').write_text(json.dumps(dict(status='COMPLETE_ALL_ANCHOR_FUNCTIONAL_AND_NUISANCE_DESCRIPTIONS_REVIEW_PENDING',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),files=files,
        annotation_slots=600,diagnostic_content_slots=30,source_component_slots=90,
        no_new_clinical_or_communication_test=True),indent=2)+'\n')
    print('ALL_FUNCTIONAL_AND_SENTINEL_DESCRIPTIONS_DONE_REVIEW_PENDING')


if __name__=='__main__':main()
