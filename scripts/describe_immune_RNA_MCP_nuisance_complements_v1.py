"""All fixed-sentinel complement point sensitivities, never new primary P tests."""
import gzip
import json
import os
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.stats import spearmanr

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from test_immune_RNA_MCP_source_coordination_v1 import fit_rank_partial
from audit_immune_RNA_MCP_stability_and_coordination_v1 import independent_fit


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget=resource(root,8*1024**2)
    audit=root/'immune_RNA_MCP_functional_description_audit_20261018_v1/private_terminal.json'
    assert sha(audit)=='e96be1d1a2630885557cec7443b931723a8f162f76dbff2b920f92b638f37b42'
    function=root/'immune_RNA_MCP_functional_description_20261018_v1'
    proof=json.loads(audit.read_text());assert sha(function/'private_terminal.json')==proof['producer_terminal_sha256']
    prior=json.loads((function/'plan_before_complete_functional_description.json').read_text())
    projection=root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection'
    weights=pd.read_csv(projection/'private_frozen_gene_weights.tsv',sep='\t')
    weights=weights[weights.fit_id.eq('K2_seed20261018_source_adjusted')]
    baseline=pd.read_csv(projection/'private_all_case_mean_scores.tsv',sep='\t',dtype={'sample':str})
    baseline=baseline[baseline.fit_id.eq('K2_seed20261018_source_adjusted')]
    contract=json.loads((root/'immune_RNA_MCP_postfit_contract_20261018_v1/plan_before_learned_outputs_or_projection.json').read_text())
    output=root/'immune_RNA_MCP_nuisance_complements_20261018_v1'
    assert not output.exists();output.mkdir(mode=0o700)
    (output/'plan_before_all_sentinel_complements.json').write_text(json.dumps(dict(
        created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),budget=budget,
        prior_function_acceptance_sha256=sha(audit),diagnostics=prior['diagnostics'],
        reason='Post-result immediate-early/stress-associated RNA contribution in one retained-source view. Examine every previously fixed sentinel class and both programmes/all view pairs/all sources.',
        scope='6pairs x5classes x3sources=90 point-only sensitivities, plus all90 original counterparts. Remove the fixed gene-class component from both endpoints using original coefficients, discovery scaling and nuisance fit; no genes re-learned.',
        inference='No new permutation P, CI or BH, no replacing primary6 tests. Descriptive sensitivity cannot prove absence of all stress/sex/ambient effects or causal communication.',
        absent_gene='A gene outside learned coefficients contributes zero to that fixed score; this does not impute missing expression or rule out indirect background effects.'),indent=2)+'\n')
    cases=[];inputs=root/'immune_RNA_MCP_native_expression_inputs_20261018_v1/private_inputs'
    controls=['cellQ','geneQ','sourceCC','fraction','logN']
    for view in contract['view_to_official_name']:
        parts={};genes=None
        for part in ['discovery','TNBC_holdout']:
            p=inputs/(view+'_'+part)
            md=pd.read_csv(str(p)+'_metadata.tsv',sep='\t',dtype={'sample':str})
            g=Path(str(p)+'_genes.tsv').read_text().splitlines();assert genes is None or genes==g;genes=g
            with gzip.open(str(p)+'_tpm.mtx.gz','rb') as h:mat=mmread(h).tocsr()
            parts[part]=(md,mat)
        md,mat=parts['discovery'];D=np.c_[np.ones(len(md)),md[controls].to_numpy(float)]
        mh,matrix_h=parts['TNBC_holdout'];Dh=np.c_[np.ones(len(mh)),mh[['cellQ','geneQ','fraction','logN']].to_numpy(float)]
        Q,R=np.linalg.qr(D,mode='reduced')
        for program in ['MCP1','MCP2']:
            w=weights[weights.view.eq(view)&weights.program.eq(program)].copy()
            for label,symbols in prior['diagnostics'].items():
                selected=w.gene.str.upper().isin(symbols).to_numpy()
                v=w.weight.to_numpy()[selected]/w.discovery_sd.to_numpy()[selected]
                ix=pd.Index(genes).get_indexer(w.gene[selected]);assert (ix>=0).all()
                mu=w.discovery_mean.to_numpy()[selected]
                raw=np.asarray(v@mat[ix]).ravel()-v@mu if selected.any() else np.zeros(len(md))
                beta=np.linalg.solve(R,Q.T@raw)
                pieces={'discovery':raw-Q@(Q.T@raw),
                    'TNBC_holdout':(np.asarray(v@matrix_h[ix]).ravel()-v@mu if selected.any() else np.zeros(len(mh)))-Dh@beta[[0,1,2,4,5]]}
                for part,piece in pieces.items():
                    metadata=parts[part][0]
                    frame=metadata[['sample','dataset','fraction','logN','cellQ','geneQ']].copy();frame['piece']=piece
                    mean=frame.groupby(['sample','dataset'],sort=True).mean().reset_index()
                    base=baseline[baseline.view.eq(view)&baseline.program.eq(program)&baseline.partition.eq(part)][['sample','score']]
                    mean=mean.merge(base,on='sample',validate='one_to_one')
                    mean['complement']=mean.score-mean.piece;mean['view']=view;mean['program']=program;mean['diagnostic_class']=label
                    cases.extend(mean.to_dict('records'))
    frame=pd.DataFrame(cases);assert len(frame)==1980
    frame.to_csv(output/'private_all_case_fixed_complements.tsv',sep='\t',index=False)
    rows=[]
    for slot,item in enumerate(contract['primary_source_association']['slots']):
        for label in prior['diagnostics']:
            a=frame[frame.view.eq(item['view_A'])&frame.program.eq(item['program'])&frame.diagnostic_class.eq(label)]
            b=frame[frame.view.eq(item['view_B'])&frame.program.eq(item['program'])&frame.diagnostic_class.eq(label)]
            pair=a.merge(b,on=['sample','dataset','program','diagnostic_class'],suffixes=('_A','_B'),validate='one_to_one')
            names=[c+'_'+s for s in ['A','B'] for c in ['fraction','logN','cellQ','geneQ']]
            for source,g in pair.groupby('dataset',sort=True):
                C=g[names].to_numpy(float);row=dict(slot=slot,program=item['program'],view_A=item['view_A'],view_B=item['view_B'],
                    dataset=source,diagnostic_class=label,case_or_sample_units=len(g),point_sensitivity_no_new_P_or_CI=True)
                for kind in ['score','complement']:
                    x=g[kind+'_A'].to_numpy();y=g[kind+'_B'].to_numpy()
                    fit,reason=fit_rank_partial(x,y,C);ref,rr=independent_fit(x,y,C)
                    assert reason==rr
                    effect=fit['effect'] if fit is not None else np.nan
                    if fit is not None:assert abs(effect-ref['effect'])<1e-9
                    row[kind+'_unadjusted_rho']=float(spearmanr(x,y).statistic) if len(g)>=3 and min(np.ptp(x),np.ptp(y))>0 else np.nan
                    row[kind+'_adjusted_rho']=effect;row[kind+'_status']=reason
                rows.append(row)
    result=pd.DataFrame(rows);assert len(result)==90
    result.to_csv(output/'private_all90_nuisance_complement_point_sensitivities.tsv',sep='\t',index=False)
    files={p.name:sha(p) for p in sorted(output.iterdir()) if p.is_file()}
    (output/'private_terminal.json').write_text(json.dumps(dict(
        status='COMPLETE_ALL90_FIXED_NUISANCE_COMPLEMENT_POINTS_TWO_NUMERIC_SOLVERS_AGREE',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),files=files,
        new_P_CI_BH_or_model_selection=False,original_six_tests_unchanged=True),indent=2)+'\n')
    print('ALL_FIXED_NUISANCE_COMPLEMENT_POINT_SENSITIVITIES_DONE')


if __name__=='__main__':main()
