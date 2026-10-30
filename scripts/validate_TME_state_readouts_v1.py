"""Development-only identity-conditioned programme readouts and split halves.

Correlations are descriptive measurement diagnostics, not patient-independent
tests or biological-function validation. Held-out score arrays are not opened.
"""
import argparse,datetime,hashlib,json,sys,time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata

def load(p):return json.loads(Path(p).read_text())
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
    return h.hexdigest()
def save(p,v):
    with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False)
def binding(p):return dict(path=str(p),SHA256=sha(p))
def bound(b):
    assert sha(b['path'])==b['SHA256'];return Path(b['path'])
def column_spearman(a,b):
    a=np.asarray(a,float);b=np.asarray(b,float)
    if b.ndim==1:b=np.broadcast_to(b[:,None],a.shape)
    assert a.ndim==2 and a.shape==b.shape and np.isfinite(a).all() and np.isfinite(b).all()
    ar=rankdata(a,axis=0,method='average');br=rankdata(b,axis=0,method='average')
    ar-=ar.mean(axis=0);br-=br.mean(axis=0);den=np.sqrt((ar*ar).sum(axis=0)*(br*br).sum(axis=0))
    return np.divide((ar*br).sum(axis=0),den,out=np.full(a.shape[1],np.nan),where=den>0)
def split_means(values,seed):
    assert len(values)>=2
    order=np.random.default_rng(np.random.SeedSequence(seed)).permutation(len(values));n=len(values)//2
    return values[order[:n]].mean(axis=0,dtype=np.float64),values[order[n:]].mean(axis=0,dtype=np.float64),n,len(values)-n

def main(R):
    cfg=load(R/'config/spec_v1.json');start=time.monotonic();A=Path(cfg['score_auditor'])
    test=Path(cfg['math_test_stage']);assert load(test/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    test_accept=load(test/'outputs/readout_math_acceptance_v1.json');assert test_accept['passed'] and test_accept['n_tests']==4
    import ast
    def math_nodes(path):
        tree=ast.parse(Path(path).read_text())
        return {n.name:ast.dump(n,include_attributes=False) for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in {'column_spearman','split_means'}}
    assert math_nodes(test/'scripts/validate_TME_state_readouts_v1.py')==math_nodes(__file__)
    while not (A/'provenance/supervisor_terminal_v1.json').exists():
        assert time.monotonic()-start<cfg['dependency_wait_seconds'];time.sleep(5)
    assert load(A/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    sys.path.insert(0,str(Path(cfg['joined'])/'scripts'))
    from extract_full_source_local_marker_evidence_v1 import storage_gate
    storage_gate(R,cfg);ap=A/'outputs/accepted_native_program_score_execution_v1.json'
    assert sha(ap)==load(A/'provenance/terminal_v1.json')['acceptance_SHA256'];audit=load(ap)
    score_index=load(bound(audit['producer_index']));P=Path(cfg['score_producer']);pc=load(P/'config/spec_v1.json')
    S=Path(cfg['applicability']);L=Path(cfg['source_roles'])
    for root in [S,L]:assert load(root/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    sap=S/'outputs/accepted_state_program_applicability_view_v1.json';app=load(sap)
    assert sha(sap)==load(S/'provenance/terminal_v1.json')['acceptance_SHA256']
    assert app['scope_policy_mode']=='resource_hierarchy_primary'
    lap=L/'outputs/accepted_state_source_roles_v1.json';roles=load(lap)
    assert sha(lap)==load(L/'provenance/terminal_v1.json')['acceptance_SHA256']
    assert app['common_cell_mask_SHA256']==roles['common_cell_mask_SHA256']==score_index['common_cell_mask_SHA256']
    role={r['source']:r['role'] for r in roles['source_roles']}
    development=np.load(bound(roles['development_mask']));ctx=np.load(bound(app['context_codes']))
    prep=Path(pc['preparation']);prep_spec=load(prep/'config/spec_v1.json')
    design=Path(prep_spec['design']);closure=Path(prep_spec['unit_closure'])
    unit_axis=np.load(bound(load(closure/'outputs/accepted_fit_qualification_unit_bridge_v3.json')['unit_codes']))
    common=np.load(prep/'outputs/private_common_assessed_cell_mask_v1.npy')
    unit_eligible=np.load(design/'outputs/private_all_qualified_tumour_unit_mask_v1.npy').copy()
    common_counts=np.bincount(unit_axis[common],minlength=len(unit_eligible))
    unit_eligible&=common_counts>=cfg['minimum_total_cells']
    applicability=np.load(bound(app['context_program_applicability']));nc=app['n_contexts']
    names=score_index['program_ids'];assert names==[p['program_id'] for p in app['programs']]
    totals=np.load(Path(pc['model'])/'outputs/private_analysis_RNA_totals_v1.npy',mmap_mode='r')
    output=[];diagnostics=[];n_cells=0;n_groups=0
    for source_number,record in enumerate(score_index['sources']):
        source=record['source']
        if role[source]!='development':continue
        with np.load(bound(record['cell_axis']),allow_pickle=False) as f:
            rows=f['fit_row'];unit=f['qualification_unit_code'];detected=f['detected_native_genes']
        assert development[rows].all();context=ctx[rows];assert (context>=0).all()
        y=np.load(bound(record['scores']),mmap_mode='r')
        measured=np.asarray([p['native_coverage_pass'] for p in record['coverage_and_interpretation']],bool)
        keys,inverse=np.unique(unit.astype(np.int64)*nc+context,return_inverse=True)
        counts=np.bincount(inverse,minlength=len(keys));order=np.argsort(inverse,kind='stable');starts=np.r_[0,np.cumsum(counts)]
        means=np.full((len(keys),len(names)),np.nan,np.float32)
        halves_a=np.full((cfg['n_splits'],len(keys),len(names)),np.nan,np.float32);halves_b=np.full_like(halves_a,np.nan)
        mean_rna=np.empty(len(keys),float);mean_detected=np.empty(len(keys),float)
        for g,key in enumerate(keys):
            positions=order[starts[g]:starts[g+1]];c=int(key%nc);u=int(key//nc)
            mean_rna[g]=np.log1p(totals[rows[positions]]).mean();mean_detected[g]=detected[positions].mean()
            if len(positions)<cfg['minimum_readout_cells']:continue
            js=np.flatnonzero(measured&applicability[c]);values=np.asarray(y[positions][:,js])
            assert np.isfinite(values).all();means[g,js]=values.mean(axis=0,dtype=np.float64)
            if len(positions)<cfg['minimum_split_cells']:continue
            for split in range(cfg['n_splits']):
                a,b,na,nb=split_means(values,[cfg['seed'],source_number,u,c,split])
                halves_a[split,g,js]=a;halves_b[split,g,js]=b
                if split==0:np.testing.assert_allclose((na*a+nb*b)/(na+nb),means[g,js],rtol=1e-5,atol=2e-7)
        for c in np.unique(keys%nc):
            group=np.flatnonzero((keys%nc==c)&(counts>=cfg['minimum_split_cells'])&unit_eligible[keys//nc])
            js=np.flatnonzero(measured&applicability[c])
            if len(group)<cfg['minimum_reliability_units']:
                for j in js:diagnostics.append(dict(source=source,context_code=int(c),program_id=names[j],status='insufficient_within_source_context_units',n_units=len(group),reliable_for_exploratory_state_patterns=False))
                continue
            correlations=np.stack([column_spearman(halves_a[s,group][:,js],halves_b[s,group][:,js]) for s in range(cfg['n_splits'])])
            depth_rna=column_spearman(means[group][:,js],mean_rna[group])
            depth_detected=column_spearman(means[group][:,js],mean_detected[group])
            for k,j in enumerate(js):
                values=correlations[:,k];finite=values[np.isfinite(values)]
                median=float(np.median(finite)) if len(finite) else None
                positive=float(np.mean(finite>0)) if len(finite) else None
                reliable=median is not None and median>=cfg['minimum_split_median_Spearman'] and positive>=cfg['minimum_positive_split_fraction']
                diagnostics.append(dict(source=source,context_code=int(c),program_id=names[j],
                    status='estimable' if len(finite) else 'no_between_unit_rank_variation',n_units=len(group),
                    n_finite_splits=len(finite),split_half_Spearman_median=median,positive_split_fraction=positive,
                    mean_score_RNA_Spearman=float(depth_rna[k]) if np.isfinite(depth_rna[k]) else None,
                    mean_score_detected_gene_Spearman=float(depth_detected[k]) if np.isfinite(depth_detected[k]) else None,
                    reliable_for_exploratory_state_patterns=reliable,
                    depth_association_is_diagnostic_not_automatic_exclusion=True))
        out=R/'outputs'/('private_'+source+'_development_state_readouts_v1.npz')
        np.savez_compressed(out,unit_code=keys//nc,context_code=keys%nc,n_cells=counts,
            unit_eligible_for_module_analysis=unit_eligible[keys//nc],
            means=means,split_means_A=halves_a,split_means_B=halves_b,
            mean_log1p_RNA=mean_rna,mean_detected_native_genes=mean_detected)
        # Measurement and identity mask apply to every emitted row.
        for g,key in enumerate(keys):
            ok=measured&applicability[int(key%nc)]&(counts[g]>=cfg['minimum_readout_cells'])
            assert np.isfinite(means[g,ok]).all() and np.isnan(means[g,~ok]).all()
        output.append(dict(source=source,n_cells=len(rows),n_unit_context_groups=len(keys),readouts=binding(out)))
        n_cells+=len(rows);n_groups+=len(keys);storage_gate(R,cfg)
        save(R/'checkpoints'/('source_'+str(source_number)+'_complete_v1.json'),dict(source=source,n_complete_cells=n_cells,n_complete_sources=len(output)))
    assert n_cells==roles['n_development_cells']
    table=R/'outputs/private_development_state_readout_diagnostics_v1.tsv';pd.DataFrame(diagnostics).to_csv(table,sep='\t',index=False)
    # Gene-set overlap is a resource property, not a correlation discovered in
    # held-out expression. Keep it separate from split and depth diagnostics.
    catalog=load(bound(pc['catalog']));genes={p['program_id']:set(p['genes']) for p in catalog['programs']};overlap=[]
    for i,a in enumerate(names):
        for b in names[i+1:]:
            common=len(genes[a]&genes[b]);union=len(genes[a]|genes[b])
            overlap.append(dict(program_a=a,program_b=b,n_shared_genes=common,Jaccard=common/union))
    op=R/'outputs/private_fixed_program_gene_overlap_v1.tsv';pd.DataFrame(overlap).to_csv(op,sep='\t',index=False)
    index=R/'outputs/accepted_development_state_readouts_v1.json'
    save(index,dict(schema='accepted_TME_development_state_readout_diagnostics_v1',
        n_cells=n_cells,n_sources=len(output),n_programs=len(names),program_ids=names,n_unit_context_groups=n_groups,
        readout_sources=output,diagnostics=binding(table),signature_overlap=binding(op),
        source_roles=binding(lap),applicability=binding(sap),score_numerical_audit=binding(ap),
        exact_tested_math_functions_reused=True,math_counterexample_acceptance=binding(test/'outputs/readout_math_acceptance_v1.json'),
        minimum_total_cells_for_reliability_units=cfg['minimum_total_cells'],
        n_descriptively_reliable_source_context_programs=sum(r['reliable_for_exploratory_state_patterns'] for r in diagnostics),
        n_lockbox_score_arrays_read=0,clinical_outcomes_used=False,
        split_halves_measure_readout_precision_not_independent_patients=True,
        depth_adjustment_thinning_and_feature_redundancy_review_not_replaced=True,
        final_state_feature_freeze_and_module_admission_still_pending=True))
    storage_gate(R,cfg);save(R/'provenance/terminal_v1.json',dict(status='DEVELOPMENT_ONLY_STATE_READOUT_SPLIT_HALF_AND_DEPTH_DIAGNOSTICS_COMPLETE',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),acceptance_SHA256=sha(index)))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
    try:main(R)
    except BaseException:
        import traceback
        save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
