"""Replay every native programme score with an independent linear expansion."""
import hashlib,json,os,time,math
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);O=R/'NSCLC_context_native_programme_projection_20261016_v1'
    launch=json.loads((R/'HNSCC_receiver_followup_job_control_20261016_v2/project_NSCLC_context_programmes_native_v1/launch_receipt.json').read_text());deadline=time.monotonic()+7200
    while not (O/'private_terminal.json').exists():
        p=Path('/proc')/str(launch['pid'])/'stat';assert p.exists()
        s=p.read_text().rsplit(')',1)[1].split();assert s[19]==str(launch['proc_start_ticks']) and s[0]!='Z'
        assert time.monotonic()<deadline;time.sleep(10)
    out=O/'private_independent_audit.json';assert not out.exists()
    report=json.loads((O/'private_terminal.json').read_text())
    for name,h in report['files'].items():assert sha(O/name)==h
    assert json.loads((R/'NSCLC_full_native_gene_unit_audit_20261016_v1/private_terminal.json').read_text())['status'].startswith('PASS_ALL_OUTPUTS_')
    V=Path(json.loads((R/'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root'])
    F=R/'NSCLC_context_frozen_parameters_20261016_v1';N=R/'NSCLC_full_native_gene_units_20261016_v1'
    params=pd.read_csv(F/'private_all53_frozen_gene_parameters.tsv',sep='\t')
    q=pd.read_csv(F/'private_all28_training_quality_transforms.tsv',sep='\t').set_index(['context','type'])
    u=pd.read_csv(N/'private_all_native_units_and_quality.tsv',sep='\t',keep_default_na=False).set_index('unit_index',drop=False)
    scores=pd.read_csv(O/'private_native_author_subtype_programme_scores.tsv.gz',sep='\t')
    scope=pd.read_csv(O/'private_all56_projection_scope.tsv',sep='\t')
    original=pd.read_csv(V/'state_after_DIALOGUE_frozen_projection_20261016_v1/public_slot_scope.tsv',sep='\t').set_index(['context','type','program'])
    domain=pd.read_csv(O/'private_all_projected_slots_training_domain.tsv',sep='\t').set_index(['context','type','program','covariate'])
    old=pd.read_csv(R/'all_CM_state_native_input_preflight_20261014_v1/private_all_frozen_sample_node_units.tsv.gz',sep='\t',keep_default_na=False).set_index('unit_index',drop=False)
    old_major=old.groupby(['sample_uid','major']).n_cells.sum()
    genes=pd.Index(pd.read_csv(N/'private_native_gene_axis.tsv',sep='\t').gene);M=np.load(N/'mean_logCP10K.npy',mmap_mode='r')
    allT=u.sub_cell_type.str.startswith(('CD4T_','CD8T_','T_gdT_'));allMac=u.sub_cell_type.str.startswith('Mφ_')
    totals={'T cell':u[allT].groupby('sampleID').n_cells.sum(),'Macrophage':u[allMac].groupby('sampleID').n_cells.sum()}
    maximum=0.;count=0
    for row in scope.itertuples():
        key=(row.context,row.type,row.program);orig=original.loc[key];qt=q.loc[(row.context,row.type)];endpoint=qt.original_endpoint
        f=scores[scores.context.eq(row.context)&scores.type.eq(row.type)&scores.program.eq(row.program)]
        family='Macrophage' if endpoint.startswith('Macrophage_') else 'CD8' if endpoint.startswith('CD8_T_cell_') else 'Treg' if endpoint.startswith('Regulatory_T_cell_') else None
        assert row.original_endpoint==endpoint and row.original_status==orig.status
        if orig.status!='FROZEN_RECONSTRUCTION_PASS':assert row.status==orig.status and f.empty;continue
        if family is None:assert row.status=='NO_CURATED_COMPATIBLE_NATIVE_LINEAGE' and f.empty;continue
        p=params[params.context.eq(row.context)&params.type.eq(row.type)&params.program.eq(row.program)]
        ix=genes.get_indexer(p.gene);take=ix>=0;fraction=math.fsum(abs(w) for w,ok in zip(p.weight,take) if ok)/math.fsum(abs(w) for w in p.weight)
        assert abs(fraction-row.native_weight_fraction)<1e-12 and bool(row.complete_native)==bool(take.all())
        if fraction<.95:assert row.status=='NATIVE_WEIGHT_BELOW_FROZEN95_RULE' and f.empty;continue
        match=allMac if family=='Macrophage' else u.sub_cell_type.str.startswith('CD8T_') if family=='CD8' else u.sub_cell_type.str.startswith('CD4T_Treg_')
        expected=u[match].sort_values('unit_index');f=f.sort_values('unit_index');assert f.unit_index.tolist()==expected.unit_index.tolist()
        assert f.sampleID.tolist()==expected.sampleID.tolist() and f.sub_cell_type.tolist()==expected.sub_cell_type.tolist()
        assert np.array_equal(f.n_cells,expected.n_cells) and f.original_endpoint.eq(endpoint).all() and f.source_family.eq(family).all()
        major='Macrophage' if family=='Macrophage' else 'T cell';n=expected.n_cells.to_numpy()
        C=np.column_stack([(expected.mean_logRNA.to_numpy()-qt.mean_logRNA)/qt.sd_logRNA,(expected.mean_loggenes.to_numpy()-qt.mean_loggenes)/qt.sd_loggenes,
            n/totals[major].reindex(expected.sampleID).to_numpy(),np.log1p(n)])
        assert np.allclose(C,f[['cellQ','geneQ','fraction','logN']],atol=1e-12,rtol=1e-12)
        p=p.loc[take];coef=p.weight.to_numpy()/p.gene_sd.to_numpy()
        constant=math.fsum(float(w)*float(mu)/float(sd) for w,mu,sd in zip(p.weight,p.gene_mean,p.gene_sd))
        raw=np.asarray(M[np.ix_(expected.unit_index.to_numpy(int),ix[take])])@coef-constant
        names=['correction_intercept','correction_cellQ','correction_geneQ','correction_fraction','correction_logN']
        correction_co=np.array([math.fsum(float(w)*float(b) for w,b in zip(p.weight,p[name])) for name in names])
        correction=correction_co[0]+C@correction_co[1:];replay=raw-correction
        for name,actual in [('raw_weighted_score',raw),('frozen_quality_composition_correction',correction),('score',replay)]:
            assert np.allclose(f[name],actual,atol=1e-7,rtol=1e-9),key
            maximum=max(maximum,float(np.max(np.abs(f[name].to_numpy()-actual))))
        assert np.array_equal(f.later_inference_cell_count_eligible,n>=20)
        root=V/'state_after_DIALOGUE_inputs_20261016_v3'/f'context_{row.context:04d}'
        md=pd.read_csv(root/(row.type+'_discovery_metadata.tsv'),sep='\t');tr=md.groupby('unit_index')[['cellQ','geneQ']].mean();a=old.loc[tr.index]
        tr['fraction']=a.n_cells.to_numpy()/np.array([old_major.loc[(s,m)] for s,m in zip(a.sample_uid,a.major)])
        tr['logN']=np.log1p(a.n_cells.to_numpy())
        flags=[]
        for j,name in enumerate(['cellQ','geneQ','fraction','logN']):
            lo=tr[name].min();hi=tr[name].max();outside=(C[:,j]<lo)|(C[:,j]>hi);flags.append(outside)
            assert np.array_equal(outside,f[name+'_outside_training_unit_range'])
            dr=domain.loc[key+(name,)];assert np.isclose(lo,dr.discovery_unit_minimum) and np.isclose(hi,dr.discovery_unit_maximum)
            assert outside.sum()==dr.native_units_outside and np.sum(outside&(n>=20))==dr.native_ge20_outside
        assert np.array_equal(np.any(flags,axis=0),f.any_covariate_outside_training_unit_range)
        count+=len(f)
    assert count==len(scores)==report['native_unit_score_records'] and len(scope)==56
    a=dict(status='PASS_ALL56_SCOPE_ALL_NATIVE_AUTHOR_LINEAGE_SCORES_AND_TRAINING_DOMAIN_FLAGS',finished_UTC=datetime.now(timezone.utc).isoformat(),
        source_terminal_sha256=sha(O/'private_terminal.json'),audit_code_sha256=sha(__file__),score_records=count,maximum_expanded_linear_replay_error=maximum,
        clinical_outcomes_used=False,fine_identity_mapping_validated=False)
    out.write_text(json.dumps(a,indent=2)+'\n');print(json.dumps(a))


if __name__=='__main__':main()
