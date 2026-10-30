"""Independent row, measurement, range and scalar-equation checks of UCell.

Numerical acceptance does not grant cell-state applicability, function or
module admission. Every source is checked, but not every value is rescored.
"""
import argparse,datetime,hashlib,json,sys,time
from pathlib import Path
import numpy as np
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

def main(R):
    cfg=load(R/'config/spec_v1.json');P=Path(cfg['producer']);began=time.monotonic()
    terminal=P/'provenance/supervisor_terminal_v1.json'
    while not terminal.exists():
        assert time.monotonic()-began<cfg['dependency_wait_seconds'];time.sleep(5)
    assert load(terminal)['actual_child_exit_code']==0
    spec=load(P/'config/spec_v1.json');J=Path(spec['joined']);H=Path(spec['model'])
    sys.path.insert(0,str(J/'scripts'))
    from native_archive_collection_v1 import NativeArchiveCollection
    from native_annotation_review_v1 import native_measured_genes
    from native_panel_aggregate_v1 import gather_csr_rows
    from audit_TME_standard_Harmony_v2 import gather_nonnegative_float_delta
    from extract_full_source_local_marker_evidence_v1 import storage_gate
    storage_gate(R,cfg)
    ip=P/'outputs/native_UCell_expression_program_score_index_v1.json'
    assert sha(ip)==load(P/'provenance/terminal_v1.json')['index_SHA256'];m=load(ip)
    assert m['mode']==spec['mode']=='full' and m['raw_scores_require_separate_biological_parent_applicability_mask']
    catalog=load(bound(spec['catalog']));names=[r['program_id'] for r in catalog['programs']]
    assert names==m['program_ids'] and len(names)==m['n_programs']
    program_genes={r['program_id']:set(r['genes']) for r in catalog['programs']}
    prep=Path(spec['preparation']);selected=np.load(prep/'outputs/private_common_assessed_cell_mask_v1.npy')
    assert sha(prep/'outputs/private_common_assessed_cell_mask_v1.npy')==m['common_cell_mask_SHA256']
    pscope=Path(spec['scope'])/'outputs/accepted_shared_compartment_scope_v1.json'
    scope=load(pscope);compartment=np.load(bound(scope['codes']))
    with np.load(H/'outputs/private_cell_axes_v1.npz',allow_pickle=False) as f:axis={k:f[k] for k in f.files}
    assert sha(H/'outputs/private_cell_axes_v1.npz')==m['input_axes_SHA256']
    src=load(H/'outputs/scope_and_common_genes_v1.json')['sources']
    unit_root=Path(load(prep/'config/spec_v1.json')['unit_closure'])
    unit=np.load(bound(load(unit_root/'outputs/accepted_fit_qualification_unit_bridge_v3.json')['unit_codes']))
    totals=np.load(H/'outputs/private_analysis_RNA_totals_v1.npy',mmap_mode='r')
    lineage={r['source']:r for r in load(H/'outputs/source_input_lineage_v1.json')}
    native=NativeArchiveCollection(load(J/'outputs/joined_native_atlas_delivery_v1.json')['archive_manifest'])
    expected_sources=[s for i,s in enumerate(src) if np.any(selected&(axis['source_code']==i))]
    assert [r['source'] for r in m['sources']]==expected_sources
    checked=0;manual=0;maximum_error=0.;zero_scores=0;rows=[]
    for record in m['sources']:
        source=record['source'];code=src.index(source)
        expected=np.flatnonzero(selected&(axis['source_code']==code))
        with np.load(bound(record['cell_axis']),allow_pickle=False) as f:
            np.testing.assert_array_equal(f['fit_row'],expected)
            np.testing.assert_array_equal(f['qualification_unit_code'],unit[expected])
            np.testing.assert_array_equal(f['compartment_code'],compartment[expected])
            assert (f['detected_native_genes']>=0).all()
        y=np.load(bound(record['scores']),mmap_mode='r')
        assert y.dtype==np.float32 and y.shape==(len(expected),len(names))
        genes,columns=native_measured_genes(native,source);genes=list(genes)
        assert len(genes)==len(set(genes))==record['n_native_ranking_genes']
        assert hashlib.sha256('\0'.join(genes).encode()).hexdigest()==record['native_gene_axis_SHA256']
        actual=[sorted(program_genes[n]&set(genes)) for n in names]
        available=np.array([len(g)>=spec['minimum_genes'] and len(g)/len(program_genes[n])>=spec['minimum_fraction'] for n,g in zip(names,actual)])
        np.testing.assert_array_equal(available,[r['native_coverage_pass'] for r in record['coverage_and_interpretation']])
        assert np.isnan(y[:,~available]).all()
        for lo in range(0,len(y),8192):
            v=y[lo:lo+8192,available]
            assert np.isfinite(v).all() and v.min()>=-1e-6 and v.max()<=1.000001
            checked+=v.size;zero_scores+=int(np.count_nonzero(v==0))
        positions=np.unique(np.linspace(0,len(expected)-1,min(len(expected),cfg['scalar_cells_per_source']),dtype=int))
        fit=expected[positions];nr=axis['native_row'][fit];z=native.store(source)
        x=gather_csr_rows(z[lineage[source]['counts_key']],nr,max_span=8192).astype(np.float64).tocsr()
        if 'ambient_removed_counts' in z['layers']:x=(x-gather_nonnegative_float_delta(z['layers/ambient_removed_counts'],nr)).tocsr()
        x.eliminate_zeros();assert x.data.min(initial=0)>=-1e-8;x.data[x.data<0]=0;x.eliminate_zeros()
        np.testing.assert_allclose(np.asarray(x.sum(axis=1)).ravel(),totals[fit],rtol=1e-10,atol=1e-7)
        values=x[:,np.asarray(columns)].toarray().astype(np.float32);lookup={g:i for i,g in enumerate(genes)}
        for i,vec in enumerate(values):
            positive=np.flatnonzero(vec>0);ranks=np.full(len(vec),spec['max_rank'],float)
            ranks[positive]=np.minimum(rankdata(-vec[positive],method='average'),spec['max_rank'])
            for j in np.flatnonzero(available):
                n=len(actual[j]);baseline=n*(n+1)/2
                expected_value=max(0.,1-(sum(ranks[lookup[g]] for g in actual[j])-baseline)/(n*spec['max_rank']-baseline))
                error=abs(expected_value-float(y[positions[i],j]));maximum_error=max(maximum_error,error)
                assert error<3e-6;manual+=1
        rows.append(dict(source=source,n_cells=len(expected),n_coverage_pass_programs=int(available.sum()),n_scalar_values=len(positions)*int(available.sum())))
        save(R/'checkpoints'/('source_'+str(code)+'_accepted_v1.json'),rows[-1]);native._fields.clear()
    assert sum(r['n_cells'] for r in rows)==m['n_cells']==int(selected.sum())
    out=R/'outputs/accepted_native_program_score_execution_v1.json'
    save(out,dict(schema='accepted_TME_native_UCell_numerical_execution_v1',
        n_cells=m['n_cells'],n_sources=len(rows),n_programs=len(names),n_values_range_and_measurement_checked=checked,
        n_explicit_observed_zero_scores=zero_scores,n_independent_scalar_scores=manual,max_scalar_error=maximum_error,
        producer_index=dict(path=str(ip),SHA256=sha(ip)),
        producer_actual_exit_zero=True,all_cell_rows_and_unit_and_compartment_axes_replayed=True,
        all_values_independently_rescored=False,all_sources_scalar_probed=True,
        applicability_readout_reliability_and_state_module_admission_not_granted=True,
        completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat()))
    storage_gate(R,cfg);save(R/'provenance/terminal_v1.json',dict(status='ALL_SOURCE_NATIVE_UCELL_NUMERICAL_AUDIT_COMPLETE_BIOLOGICAL_READOUT_GATES_SEPARATE',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),acceptance_SHA256=sha(out)))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
    try:main(R)
    except BaseException:
        import traceback
        save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
