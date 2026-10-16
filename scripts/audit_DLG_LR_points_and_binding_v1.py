"""Independent all-estimable OLS point/sensitivity replay and full scope/BH.

This audit does not claim replay of every original permutation/bootstrap draw.
Curated binding failures remain in the original statistical family and q values.
"""
import hashlib,json,os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import rankdata
from statsmodels.stats.multitest import multipletests
from threadpoolctl import threadpool_limits


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def corr(x,y,D=None):
    a=rankdata(x).astype(float);b=rankdata(y).astype(float)
    a-=a.mean();b-=b.mean()
    if D is not None:
        a-=D@np.linalg.lstsq(D,a,rcond=None)[0];b-=D@np.linalg.lstsq(D,b,rcond=None)[0]
    return float(a@b/(np.linalg.norm(a)*np.linalg.norm(b)))


def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    V=Path(json.loads((R/'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root'])
    S=V/'state_after_DIALOGUE_LR_response_20261016_v1';iv=json.loads((S/'inventory.json').read_text())
    for name,h in iv['files'].items():assert sha(S/name)==h,name
    t=pd.read_csv(S/'all_estimable_exploratory_LR_receiver_associations.tsv',sep='\t')
    counts=pd.read_csv(S/'all258944_hypothesis_scope_counts.tsv',sep='\t')
    assert counts.records.sum()==258944 and len(t)==iv['estimable_tests']
    assert not t.duplicated(['context','source','target','MCP','ligand','receptor']).any()
    expected=multipletests(t.permutation_p,method='fdr_bh')[1];assert np.allclose(expected,t.BH_q,atol=1e-12)
    Q=V/'all_source_native_UCell_unit_sensitivity_20261014_v1';qi=json.loads((Q/'inventory.json').read_text());quality={}
    for info in qi['private_outputs']:
        p=Q/info['relative_file'];assert sha(p)==info['sha256'];z=np.load(p,allow_pickle=False)
        quality.update({int(k):(float(r),float(g)) for k,r,g in zip(z['unit_index'],z['mean_full_RNA'],z['mean_full_detected_genes'])})
    E=R/'immune_same_name_family_source_review_v1/CM_edge_LR_extraction_v2';ei=json.loads((E/'inventory.json').read_text())
    alignment=pd.read_csv(V/'state_after_DIALOGUE_all_LR_preflight_20261016_v1/all28_cache_alignment_checks.tsv',sep='\t')
    component_root=V/'state_after_DIALOGUE_gene_components_20261016_v1';ci=json.loads((component_root/'inventory.json').read_text())
    p=component_root/'private_all_block_gene_components.tsv.gz';assert sha(p)==ci['files'][p.name]
    components=pd.read_csv(p,sep='\t')
    out=R/'DLG_LR_point_binding_review_20261016_v1';assert not out.exists();out.mkdir()
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),
        original_result_sha256=sha(S/'all_estimable_exploratory_LR_receiver_associations.tsv'),
        scope='All original estimable records, original full scope and BH, original interface-gene masks and three sensitivity indicators.',
        method='Separate rank vectors regressed by full least-squares nuisance, not producer QR residual implementation.',
        randomization='Original p and bootstrap distributions are not regenerated here; CI availability and total-family BH checked, report this scope precisely.',
        binding='Cross-reference every corrected-significant record to existing constituent evidence without modifying p/q or filling unsupported references.')
    (out/'plan_before_independent_point_review.json').write_text(json.dumps(plan,indent=2)+'\n')
    rows=[]
    with threadpool_limits(limits=1):
        for ar in alignment.itertuples():
            partroot=V/'state_after_DIALOGUE_inputs_20261016_v3'/f'context_{int(ar.context):04d}'
            original=json.loads((partroot/'private_input_inventory.json').read_text())
            part=pd.read_csv(partroot/'private_block_library_partition.tsv',sep='\t')
            reverse=ar.source==original['endpoints'][1];typ='typeA' if reverse else 'typeB'
            path=E/'private_group_arrays'/(ar.old_group+'.npz');assert sha(path)==ei['outputs'][str(path.relative_to(E))]
            z=np.load(path,allow_pickle=False);indices=pd.Index(z['samples'].astype(str)).get_indexer(part.sample_uid.astype(str));assert (indices>=0).all()
            assert list(z['patient_blocks'][indices])==list(part.effective_resampling_block)
            source_units=part['unit_index_B' if reverse else 'unit_index_A'];target_units=part['unit_index_A' if reverse else 'unit_index_B']
            a=np.asarray([quality[int(u)] for u in source_units]);b=np.asarray([quality[int(u)] for u in target_units])
            cov=np.column_stack((np.log1p(z['source_n'][indices]),np.log1p(z['target_n'][indices]),z['source_within_major'][indices],z['target_within_major'][indices],np.log1p(a[:,0]),np.log1p(b[:,0]),np.log1p(a[:,1]),np.log1p(b[:,1])))
            lookup={(str(l),str(r)):i for i,(l,r) in enumerate(zip(z['ligands'],z['receptors']))}
            for mcp in ['MCP1','MCP2']:
                selected=t[t.context.eq(ar.context)&t.source.eq(ar.source)&t.target.eq(ar.target)&t.MCP.eq(mcp)]
                expected_n=counts[counts.context.eq(ar.context)&counts.source.eq(ar.source)&counts.target.eq(ar.target)&counts.MCP.eq(mcp)&counts.status.eq('ESTIMABLE_EXPLORATORY_FULL_BLOCK_ASSOCIATION')].records.sum()
                assert len(selected)==expected_n
                if not len(selected):continue
                c=components[components.context.eq(ar.context)&components.type.eq(typ)&components.program.eq(mcp)]
                W=c.groupby('gene').weight.first();B=c.pivot(index='sample',columns='gene',values='contribution').reindex(part.effective_resampling_block)
                w=W.reindex(B.columns).to_numpy();full=B.to_numpy()@w
                for row in selected.itertuples():
                    j=lookup[row.ligand,row.receptor];scores=z['scores'][indices,j,:];mask=np.isfinite(scores[:,0]);assert mask.sum()==row.blocks and z['native_measured'][j]
                    remove=set(row.ligand.split('_')+row.receptor.split('_'))&set(B.columns)
                    assert remove==set(str(row.masked_genes).split(';')) if pd.notna(row.masked_genes) else not remove
                    keep=~B.columns.isin(remove);Y=(B.loc[:,keep].to_numpy()@w[keep])[mask];X=np.log1p(scores[mask,0])
                    cv=cov[mask];D=np.column_stack((np.ones(len(X)),(cv-cv.mean(axis=0))/cv.std(axis=0)))
                    assert np.linalg.matrix_rank(D)==9 and len(X)-10>=10
                    computed=[corr(X,Y),corr(X,Y,D),corr(X,full[mask],D),corr(np.log1p(scores[mask,1]),Y,D),corr(-scores[mask,2],Y,D),corr(-scores[mask,3],Y,D)]
                    saved=[row.raw_masked_rho,row.adjusted_masked_rho,row.adjusted_unmasked_rho,row.lr_means_adjusted_rho,row.negative_specificity_rank_adjusted_rho,row.negative_magnitude_rank_adjusted_rho]
                    assert np.allclose(computed,saved,atol=1e-8,rtol=1e-8,equal_nan=True),(row.context,row.ligand,row.receptor)
                    tr=part.partition.eq('discovery').to_numpy()[mask];te=part.partition.eq('holdout').to_numpy()[mask]
                    if tr.sum()>=20 and te.sum()>=10 and np.linalg.matrix_rank(D[tr])==9:
                        xh=X[te]-D[te]@np.linalg.lstsq(D[tr],X[tr],rcond=None)[0];yh=Y[te]-D[te]@np.linalg.lstsq(D[tr],Y[tr],rcond=None)[0]
                        assert np.isclose(corr(xh,yh),row.internal_holdout_rho,atol=1e-8)
                    else:assert pd.isna(row.internal_holdout_rho)
                    rows.append(dict(context=row.context,source=row.source,target=row.target,MCP=mcp,ligand=row.ligand,receptor=row.receptor,max_point_error=float(np.max(np.abs(np.asarray(computed)-saved)))))
    assert len(rows)==len(t)
    pd.DataFrame(rows).to_csv(out/'all4012_private_point_replay.tsv',sep='\t',index=False)
    reference_path=R/'LR_binding_resource_review_v1/all4624_reference_evidence_status.tsv';reference=pd.read_csv(reference_path,sep='\t')
    significant=t[t.BH_q.lt(.05)].copy();joined=significant.merge(reference,on=['ligand','receptor'],how='left',validate='many_to_one')
    assert joined.reference_status.notna().all()
    joined['constituent_reference_eligible']=joined.reference_status.isin(['simplex_reference_supported','all_constituents_supported_complex_assembly_unverified'])
    joined.to_csv(out/'all14_private_binding_and_point_evidence.tsv',sep='\t',index=False,na_rep='NA')
    summary=dict(status='PASS_ALL4012_POINT_SENSITIVITY_HOLDOUT_SCOPE_AND_BH_REFERENCE_CROSSCHECK_NOT_RANDOMIZATION_REPLAY',
        finished_UTC=datetime.now(timezone.utc).isoformat(),original_hypotheses=258944,estimated_records=len(rows),corrected_records=len(joined),
        maximum_point_error=max(r['max_point_error'] for r in rows),reference_eligible_records=int(joined.constituent_reference_eligible.sum()),
        binding_rejected_records=int((~joined.constituent_reference_eligible).sum()),
        original_result_sha256=plan['original_result_sha256'],binding_reference_sha256=sha(reference_path),
        evidence_sha256=sha(out/'all14_private_binding_and_point_evidence.tsv'),
        original_p_q_or_candidate_family_changed=False,randomization_draws_replayed=False,direct_binding_in_current_tumor_verified=False)
    (out/'private_terminal.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))


if __name__=='__main__':main()
