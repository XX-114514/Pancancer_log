"""Recompute all functional annotations and source RNA diagnostic partitions."""
import gzip
import json
import os
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.special import gammaln,logsumexp

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def close(a,b):
    assert (pd.isna(a) and pd.isna(b)) or np.isclose(a,b,rtol=1e-8,atol=1e-10)


def exact_tail(M,N,n,k):
    lo=max(k,n-(M-N),0);hi=min(N,n);x=np.arange(lo,hi+1)
    choose=lambda a,b:gammaln(a+1)-gammaln(b+1)-gammaln(a-b+1)
    return min(1.,float(np.exp(logsumexp(choose(N,x)+choose(M-N,n-x)-choose(M,n)))))


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget=resource(root,4*1024**2)
    source=root/'immune_RNA_MCP_functional_description_20261018_v1'
    terminal=source/'private_terminal.json';record=json.loads(terminal.read_text())
    assert record['status']=='COMPLETE_ALL_ANCHOR_FUNCTIONAL_AND_NUISANCE_DESCRIPTIONS_REVIEW_PENDING'
    for name,digest in record['files'].items():assert sha(source/name)==digest
    plan=json.loads((source/'plan_before_complete_functional_description.json').read_text())
    refs_path=Path(__file__).parent/'hallmark50_public_reference_v3'
    assert sha(refs_path/'reference_provenance.json')==plan['public_reference_manifest_sha256']
    refs={}
    for r in json.loads((refs_path/'reference_provenance.json').read_text())['records']:
        p=refs_path/(r['name']+'.json');assert sha(p)==r['sha256']
        refs[r['name']]=set(g.upper() for g in json.loads(p.read_text())[r['name']]['geneSymbols'])
    output=root/'immune_RNA_MCP_functional_description_audit_20261018_v1'
    assert not output.exists();output.mkdir(mode=0o700)
    (output/'plan_before_all_functional_replay.json').write_text(json.dumps(dict(
        created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),budget=budget,
        producer_terminal_sha256=sha(terminal),scope='All600 signed Hallmark/background overlaps and independent combinatorial tail sums/BH. All30 gene-content and90 source RNA partitions; dense standardized gene arithmetic and independent QR nuisance projection. No effect refit, signature selection or formal identity.'),indent=2)+'\n')
    projection=root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection'
    weights=pd.read_csv(projection/'private_frozen_gene_weights.tsv',sep='\t')
    weights=weights[weights.fit_id.eq('K2_seed20261018_source_adjusted')].copy()
    weights['symbol']=weights.gene.str.upper()
    total_case=pd.read_csv(projection/'private_all_case_mean_scores.tsv',sep='\t',dtype={'sample':str})
    total_case=total_case[total_case.fit_id.eq('K2_seed20261018_source_adjusted')]
    ann=pd.read_csv(source/'private_all600_Hallmark_signed_annotations.tsv',sep='\t')
    content=pd.read_csv(source/'private_all30_sentinel_gene_content_descriptions.tsv',sep='\t')
    components=pd.read_csv(source/'private_all90_source_sentinel_RNA_component_descriptions.tsv',sep='\t')
    assert len(ann)==600 and len(content)==30 and len(components)==90
    inputs=root/'immune_RNA_MCP_native_expression_inputs_20261018_v1/private_inputs'
    maximum_error=0.;source_rows=0
    for view in ['R_FOXP3_capture','E_effector_capture','Myeloid']:
        genes=Path(str(inputs/(view+'_discovery'))+'_genes.tsv').read_text().splitlines()
        background=set(g.upper() for g in genes);assert len(background)==5000
        parts={}
        for part in ['discovery','TNBC_holdout']:
            prefix=inputs/(view+'_'+part)
            md=pd.read_csv(str(prefix)+'_metadata.tsv',sep='\t',dtype={'sample':str})
            with gzip.open(str(prefix)+'_tpm.mtx.gz','rb') as h:matrix=mmread(h).tocsr()
            parts[part]=(md,matrix)
        md,matrix=parts['discovery'];controls=['cellQ','geneQ','sourceCC','fraction','logN']
        D=np.column_stack([np.ones(len(md)),md[controls]])
        Q,R=np.linalg.qr(D,mode='reduced');assert np.linalg.matrix_rank(R)==6
        mh,_=parts['TNBC_holdout'];Dh=np.column_stack([np.ones(len(mh)),mh[['cellQ','geneQ','fraction','logN']]])
        for program in ['MCP1','MCP2']:
            w=weights[weights.view.eq(view)&weights.program.eq(program)]
            ix=pd.Index(genes).get_indexer(w.gene);assert (ix>=0).all()
            for direction,sign in [('positive',1),('negative',-1)]:
                active=set(w.loc[w.weight*sign>0,'symbol'])
                a=ann[ann.view.eq(view)&ann.program.eq(program)&ann.direction.eq(direction)]
                assert len(a)==50 and set(a.reference)==set(refs)
                for row in a.itertuples(index=False):
                    measured=refs[row.reference]&background;M=len(background);N=len(measured);n=len(active);k=len(active&measured)
                    assert (row.candidate_background,row.reference_in_background,row.signed_active_genes,row.overlap)==(M,N,n,k)
                    p=exact_tail(M,N,n,k) if N and n else np.nan
                    if N and n:assert np.isclose(p,row.annotation_p,rtol=1e-8,atol=1e-300)
                    else:assert pd.isna(row.annotation_p)
                    close(row.enrichment_fold,k*M/(n*N) if N and n else np.nan)
            for label,symbols in plan['diagnostics'].items():
                take=w.symbol.isin(symbols).to_numpy();raw={}
                for part,(meta,mat) in parts.items():
                    if take.any():
                        Z=(mat[ix[take]].toarray()-w.discovery_mean.to_numpy()[take,None])/w.discovery_sd.to_numpy()[take,None]
                        raw[part]=w.weight.to_numpy()[take]@Z
                    else:raw[part]=np.zeros(len(meta))
                beta=np.linalg.solve(R,Q.T@raw['discovery'])
                values={'discovery':raw['discovery']-Q@(Q.T@raw['discovery']),
                    'TNBC_holdout':raw['TNBC_holdout']-Dh@beta[[0,1,2,4,5]]}
                c=content[content.view.eq(view)&content.program.eq(program)&content.diagnostic_class.eq(label)].iloc[0]
                assert c.active_sentinel_genes==take.sum() and c.sentinels_in_candidate_background==len(set(symbols)&background)
                close(c.absolute_coefficient_fraction,np.abs(w.weight.to_numpy()[take]).sum()/np.abs(w.weight.to_numpy()).sum())
                for part,v in values.items():
                    meta=parts[part][0];sample,inv=np.unique(meta['sample'],return_inverse=True)
                    mean=np.bincount(inv,weights=v)/np.bincount(inv)
                    # Every sample has one source; map rather than positional source order.
                    source_map=dict(zip(meta['sample'],meta.dataset));dataset=np.array([source_map[s] for s in sample])
                    baseline=total_case[total_case.view.eq(view)&total_case.program.eq(program)&total_case.partition.eq(part)].set_index('sample').score
                    b=baseline.loc[sample].to_numpy()
                    for dataset_name in np.unique(dataset):
                        keep=dataset==dataset_name;group=mean[keep];total=b[keep];rest=total-group
                        variance=np.var(total,ddof=1)
                        cov=np.cov(group,total,ddof=1)[0,1];other=np.cov(rest,total,ddof=1)[0,1]
                        total_ss=variance*(len(group)-1);group_ss=np.var(group,ddof=1)*(len(group)-1)
                        fraction=cov/variance if total_ss>1e-20 else np.nan
                        complement=other/variance if total_ss>1e-20 else np.nan
                        correlation=np.corrcoef(group,total)[0,1] if group_ss>1e-20 and total_ss>1e-20 else np.nan
                        row=components[components.view.eq(view)&components.program.eq(program)&components.partition.eq(part)&components.dataset.eq(dataset_name)&components.diagnostic_class.eq(label)]
                        assert len(row)==1;row=row.iloc[0];assert row.case_or_sample_units==keep.sum()
                        for value,saved in [(fraction,row.group_total_covariance_fraction),(complement,row.complement_total_covariance_fraction),(correlation,row.component_total_correlation)]:
                            close(value,saved)
                            if np.isfinite(value):maximum_error=max(maximum_error,abs(value-saved))
                        source_rows+=1
    p=ann.annotation_p.fillna(1).to_numpy();m=len(p)
    q=np.array([min(1.,min(m*b/np.count_nonzero(p<=b) for b in p[p>=a])) for a in p])
    q[ann.annotation_p.isna()]=np.nan
    assert np.allclose(q,ann.annotation_BH600,equal_nan=True,rtol=1e-8,atol=1e-300)
    assert source_rows==90
    (output/'private_terminal.json').write_text(json.dumps(dict(
        status='ACCEPTED_ALL600_FUNCTIONAL_30_CONTENT_90_SOURCE_RNA_DESCRIPTIONS',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),
        producer_terminal_sha256=sha(terminal),annotation_slots=600,content_slots=30,source_component_slots=90,
        maximum_source_component_difference=maximum_error,descriptive_not_clinical_or_causal_validation=True),indent=2)+'\n')
    print('ALL_FUNCTIONAL_AND_NUISANCE_DESCRIPTIONS_INDEPENDENTLY_ACCEPTED')


if __name__=='__main__':main()
