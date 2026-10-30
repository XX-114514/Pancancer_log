"""Official UCell on the new TME common domain and native gene universes.

This produces diagnostic expression-program scores, not biological state or
module admission. Reference applicability and reliable within-identity readout
are separate gates. No learned parameters or clinical outcomes are consumed.
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
def bound(b):
    assert sha(b['path'])==b['SHA256'];return Path(b['path'])
def binding(p):return dict(path=str(p),SHA256=sha(p))

def scalar_scores(x,genes,signatures,max_rank):
    """Independent sorted positive-rank equation, not the package scorer."""
    lookup={g:i for i,g in enumerate(genes)};out=np.empty((x.shape[0],len(signatures)))
    positions=[]
    for sig in signatures.values():
        assert sig and all(g in lookup for g in sig);positions.append([lookup[g] for g in sig])
    for i,vec in enumerate(x.toarray().astype(np.float32)):
        positive=np.flatnonzero(vec>0);ranks=np.full(len(vec),max_rank,float)
        ranks[positive]=np.minimum(rankdata(-vec[positive],method='average'),max_rank)
        for j,pos in enumerate(positions):
            n=len(pos);base=n*(n+1)/2
            out[i,j]=max(0.,1-(sum(ranks[g] for g in pos)-base)/(n*max_rank-base))
    return out

def main(R):
    cfg=load(R/'config/spec_v1.json');started=time.monotonic()
    J=Path(cfg['joined']);H=Path(cfg['model']);P=Path(cfg['preparation']);C=Path(cfg['coverage']);S=Path(cfg['scope'])
    sys.path.insert(0,str(J/'scripts'));sys.path.insert(0,str(Path(cfg['vendor_parent'])))
    from native_archive_collection_v1 import NativeArchiveCollection
    from native_annotation_review_v1 import native_measured_genes
    from native_panel_aggregate_v1 import gather_csr_rows
    from audit_TME_standard_Harmony_v2 import gather_nonnegative_float_delta
    from extract_full_source_local_marker_evidence_v1 import storage_gate
    import anndata,pyucell as uc
    assert uc.__version__=='0.8.0' and Path(uc.__file__).resolve().is_relative_to(Path(cfg['vendor_parent']).resolve());storage_gate(R,cfg)
    for root in cfg.get('prior_accepted_stages',[]):
        p=Path(root)/'provenance/supervisor_terminal_v1.json'
        while not p.exists():
            assert time.monotonic()-started<cfg['dependency_wait_seconds'];time.sleep(5)
        assert load(p)['actual_child_exit_code']==0,'Prior annotation/module or smoke acceptance failed'
    for root in [P,C,S]:assert load(root/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    coverage_p=C/'outputs/state_native_coverage_preparation_v1.json'
    assert sha(coverage_p)==load(C/'provenance/terminal_v1.json')['manifest_SHA256']
    cov=load(coverage_p);catalog=load(bound(cfg['catalog']));curation=load(bound(cfg['curation']))
    assert cfg['catalog']['SHA256']==cov['catalog_SHA256']==curation['scored_catalog_SHA256']
    assert cfg['curation']['SHA256']==cov['curation_SHA256']
    assert not curation['all_module_admission'];programs=catalog['programs'];names=[x['program_id'] for x in programs]
    assert len(names)==len(set(names))==cov['n_programs']
    src_cov=load(C/'outputs/private_state_program_source_native_coverage_v1.json')
    by_cov={(x['source'],x['program_id']):x for x in src_cov}
    scope=load(S/'outputs/accepted_shared_compartment_scope_v1.json')
    assert sha(S/'outputs/accepted_shared_compartment_scope_v1.json')==load(S/'provenance/terminal_v1.json')['acceptance_SHA256']
    compartments=np.load(bound(scope['codes']));selected=np.load(P/'outputs/private_common_assessed_cell_mask_v1.npy')
    auditor=Path(cfg['input_auditor']);assert load(auditor/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    iap=auditor/'outputs/accepted_comparable_CM_inputs_v1.json';assert sha(iap)==load(auditor/'provenance/terminal_v1.json')['acceptance_SHA256']
    one=load(bound(load(iap)['inputs'][0]['acceptance']))
    assert sha(P/'outputs/private_common_assessed_cell_mask_v1.npy')==one['common_assessed_cell_mask_SHA256']
    membership=load(bound(one['annotation_acceptance']))
    assert membership['scope_acceptance']['SHA256']==sha(S/'outputs/accepted_shared_compartment_scope_v1.json')
    assert selected.dtype.kind=='b' and len(selected)==scope['n_fit'];assert (compartments[selected]>=0).all()
    assert scope['input_axis_SHA256']==sha(H/'outputs/private_cell_axes_v1.npz')
    with np.load(H/'outputs/private_cell_axes_v1.npz',allow_pickle=False) as f:axes={k:f[k] for k in f.files}
    src_names=load(H/'outputs/scope_and_common_genes_v1.json')['sources']
    common_genes=load(H/'outputs/scope_and_common_genes_v1.json')['genes']
    totals=np.load(H/'outputs/private_analysis_RNA_totals_v1.npy',mmap_mode='r')
    native=NativeArchiveCollection(load(J/'outputs/joined_native_atlas_delivery_v1.json')['archive_manifest'])
    lineage={x['source']:x for x in load(H/'outputs/source_input_lineage_v1.json')}
    spec=load(P/'config/spec_v1.json');U=Path(spec['unit_closure'])
    units=np.load(bound(load(U/'outputs/accepted_fit_qualification_unit_bridge_v3.json')['unit_codes']))
    ncomp=len(scope['compartments']);rng=np.random.default_rng(cfg['seed'])
    source_bindings=[];ncomplete=0;all_probes=[]
    vendor_files=sorted(Path(cfg['vendor_parent']).glob('pyucell/*.py'))
    save(R/'provenance/frozen_official_pyucell_source_bindings_v1.json',[binding(p) for p in vendor_files])
    for source_code,source in enumerate(src_names):
        if cfg['mode']=='smoke' and source not in cfg['smoke_sources']:continue
        ids=np.flatnonzero(selected&(axes['source_code']==source_code))
        if cfg['mode']=='smoke':
            take=[]
            for p in np.unique(compartments[ids]):
                x=ids[compartments[ids]==p];take.extend(rng.choice(x,min(len(x),cfg['smoke_cells_per_compartment']),replace=False).tolist())
            ids=np.sort(take)
        if not len(ids):continue
        native_names,physical_columns=native_measured_genes(native,source)
        assert len(native_names)==len(set(native_names)) and set(common_genes)<=set(native_names)
        native_names=list(native_names);native_set=set(native_names);ph=np.asarray(physical_columns)
        sig={};score_cols=[];reasons=[]
        for j,pr in enumerate(programs):
            c=by_cov[source,pr['program_id']];actual=sorted(set(pr['genes'])&native_set)
            assert actual==c['native_measured_genes']
            valid=len(actual)>=cfg['minimum_genes'] and len(actual)/len(set(pr['genes']))>=cfg['minimum_fraction']
            assert valid==c['native_coverage_pass']
            if valid:
                assert len(actual)<cfg['max_rank'];sig[pr['program_id']]=actual;score_cols.append(j)
            reasons.append(dict(program_id=pr['program_id'],native_coverage_pass=valid,
                n_measured_genes=len(actual),n_original_genes=len(set(pr['genes'])),
                biological_applicability_and_module_admission_pending=True))
        score_path=R/'outputs'/('private_'+source+'_native_UCell_v1.npy')
        assert sig,'No genuinely measurable signatures; report this scope before scoring'
        scores=np.lib.format.open_memmap(score_path,mode='w+',dtype='float32',shape=(len(ids),len(names)));scores[:]=np.nan
        detected=np.empty(len(ids),np.int32);z=native.store(source);key=lineage[source]['counts_key']
        for lo in range(0,len(ids),cfg['read_block_size']):
            hi=min(lo+cfg['read_block_size'],len(ids));fit_rows=ids[lo:hi];nr=axes['native_row'][fit_rows]
            x=gather_csr_rows(z[key],nr,max_span=8192).astype(np.float64).tocsr()
            if 'ambient_removed_counts' in z['layers']:x=(x-gather_nonnegative_float_delta(z['layers/ambient_removed_counts'],nr)).tocsr()
            x.eliminate_zeros();assert np.isfinite(x.data).all() and x.data.min(initial=0)>=-1e-8
            x.data[x.data<0]=0;x.eliminate_zeros()
            np.testing.assert_allclose(np.asarray(x.sum(axis=1)).ravel(),totals[fit_rows],rtol=1e-10,atol=1e-7)
            x=x[:,ph].astype(np.float32).tocsr();detected[lo:hi]=np.diff(x.indptr)
            a=anndata.AnnData(x,var=pd.DataFrame(index=native_names))
            uc.compute_ucell_scores(a,signatures=sig,max_rank=cfg['max_rank'],missing_genes='skip',ties_method='average',chunk_size=cfg['score_chunk_size'],n_jobs=cfg['threads'],device='cpu')
            y=a.obs[[n+'_UCell' for n in sig]].to_numpy(np.float32)
            assert np.isfinite(y).all() and y.min()>=-1e-6 and y.max()<=1.000001
            scores[lo:hi,score_cols]=y
            if lo==0:
                probe_n=min(4,x.shape[0]);manual=scalar_scores(x[:probe_n],native_names,sig,cfg['max_rank'])
                np.testing.assert_allclose(manual,y[:probe_n],rtol=1e-5,atol=3e-6)
                # UCell ranks are cell-wise; the same selected input is also
                # evaluated with a different package chunk layout.
                b=anndata.AnnData(x[:probe_n].copy(),var=pd.DataFrame(index=native_names))
                uc.compute_ucell_scores(b,signatures=sig,max_rank=cfg['max_rank'],missing_genes='skip',ties_method='average',chunk_size=1,n_jobs=1,device='cpu')
                np.testing.assert_array_equal(b.obs[[n+'_UCell' for n in sig]].to_numpy(np.float32),y[:probe_n])
                all_probes.append(dict(source=source,n_manual_scores=manual.size,
                    max_manual_error=float(np.max(np.abs(manual-y[:probe_n]))),package_chunk_invariance=True))
            (R/'checkpoints/current_scoring_progress_v1.json').write_text(json.dumps(dict(source=source,processed_source_cells=hi,n_source_cells=len(ids),n_completed_cells=ncomplete)))
        scores.flush();rows_p=R/'outputs'/('private_'+source+'_score_cell_axis_v1.npz')
        np.savez(rows_p,fit_row=ids,qualification_unit_code=units[ids],compartment_code=compartments[ids],detected_native_genes=detected)
        assert np.array_equal(np.load(rows_p)['fit_row'],ids)
        available=np.zeros(len(names),bool);available[score_cols]=True
        assert np.isnan(scores[:,~available]).all() and np.isfinite(scores[:,available]).all()
        source_bindings.append(dict(source=source,n_cells=len(ids),n_native_ranking_genes=len(native_names),
            native_gene_axis_SHA256=hashlib.sha256('\0'.join(native_names).encode()).hexdigest(),
            n_coverage_pass_programs=len(sig),coverage_and_interpretation=reasons,
            scores=binding(score_path),cell_axis=binding(rows_p),rank_universe='source_native_uniquely_measured_gene_symbols'))
        ncomplete+=len(ids);storage_gate(R,cfg);native._fields.clear()
        save(R/'checkpoints'/('source_'+str(source_code)+'_complete_v1.json'),dict(source=source,n_complete_cells=ncomplete,n_sources_complete=len(source_bindings)))
    assert ncomplete>0
    if cfg['mode']=='full':assert ncomplete==int(selected.sum())
    output=R/'outputs/native_UCell_expression_program_score_index_v1.json'
    save(output,dict(schema='TME_native_UCell_expression_program_scores_v1',mode=cfg['mode'],
        n_cells=ncomplete,n_sources=len(source_bindings),n_programs=len(names),program_ids=names,
        catalog=cfg['catalog'],curation=cfg['curation'],input_axes_SHA256=scope['input_axis_SHA256'],
        common_cell_mask_SHA256=sha(P/'outputs/private_common_assessed_cell_mask_v1.npy'),
        scoring_method='official_pyucell_0.8.0_no_smoothing',max_rank=cfg['max_rank'],
        native_missing_signature_genes_skipped_after_explicit_coverage_gate=True,
        unavailable_program_values_are_NaN_not_zero=True,
        raw_scores_require_separate_biological_parent_applicability_mask=True,
        original_article_curation_functional_and_module_admission_not_granted=True,
        no_old_TNK_scores_copied=True,sources=source_bindings,scalar_and_chunk_probes=all_probes,
        clinical_outcomes_used=False,wall_seconds=time.monotonic()-started))
    save(R/'provenance/terminal_v1.json',dict(status='OFFICIAL_NATIVE_UCELL_EXPRESSION_PROGRAM_SCORES_COMPLETE_APPLICABILITY_AND_READOUT_GATES_SEPARATE',
        completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),mode=cfg['mode'],n_cells=ncomplete,n_sources=len(source_bindings),index_SHA256=sha(output)))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
    try:main(R)
    except BaseException:
        import traceback
        save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
