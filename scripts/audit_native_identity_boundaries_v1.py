"""Independently replay a fixed native-RNA boundary probe, including all counts.

This validates measurement and sampling, not identity truth or a doublet call.
It neither imports the producer nor uses its CSR gathering/statistic helpers.
"""
import argparse,collections,datetime,gc,hashlib,json,resource,sys,time
from pathlib import Path
import numpy as np

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):Path(p).write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False)+'\n')

def spread(choices,limit):
    n=len(choices);k=min(n,limit)
    if k<=1:return choices[:k]
    indexes=(np.arange(k,dtype=np.int64)*(n-1))//(k-1)
    return choices[indexes]

def qualification(s,view,core):
    ds=s['source'];q=np.empty(len(s['rows']),np.int8)
    if ds in view['new_source_overrides']:
        entry=view['new_source_overrides'][ds];assert sha(entry['path'])==entry['sha256']
        with np.load(entry['path'],allow_pickle=False) as f:np.testing.assert_array_equal(f['native_row'],s['rows']);q[:]=f['source_qualification_code']
    else:
        previous=load(view['base_v4_manifest']['path']);origin=load(previous['bindings']['author_bridge_registry']['path'])['source']
        if ds==origin:
            assert sha(previous['override_arrays_path'])==previous['override_arrays_sha256']
            with np.load(previous['override_arrays_path'],allow_pickle=False) as f:np.testing.assert_array_equal(f['native_row'],s['rows']);q[:]=f['source_qualification_code']
        else:
            by={u['qualification_unit_key_sha256']:u for u in core}
            for j,u in enumerate(s['um']['unit_records']):q[s['ua']['native_expression_unit_code']==j]={'排除':0,'来源资格通过':1,'待核实':2}[by[u['native_unit_key_sha256']]['source_qualification']]
    return q

def main(root):
    start=time.monotonic();cfg=load(root/'config/spec_v1.json');R=Path(cfg['producer'])
    t=load(R/'provenance/terminal_v1.json');s=load(R/'provenance/supervisor_terminal_v1.json');assert s['actual_child_exit_code']==0 and not s['wall_timeout']
    assert sha(R/'outputs/native_identity_boundary_probe_index_v1.json')==t['index_SHA256'] and sha(R/'outputs/native_identity_boundary_probe_summary_v1.json')==t['summary_SHA256']
    config=load(R/'config/spec_v1.json');J=Path(config['joined']);sys.path.insert(0,str(J/'scripts'))
    for b in config['bindings'].values():assert sha(b['path'])==b['sha256']
    from joined_native_atlas_delivery_v1 import JoinedNativeAtlas
    atlas=JoinedNativeAtlas(J/'outputs/joined_native_atlas_delivery_v1.json');view=load(config['source_view']);core=load(view['work_units']['path'])['units'];source_n=collections.Counter()
    assert sha(view['work_units']['path'])==view['work_units']['sha256'] and sha(view['base_v4_manifest']['path'])==view['base_v4_manifest']['sha256']
    for u in core:
        if u['source_qualification']=='来源资格通过':source_n[u['source']]+=u['n_archive_rows']
    ordered=sorted(source_n,key=lambda ds:(source_n[ds],ds));expected_sources=[ordered[i] for i in spread(np.arange(len(ordered)),8)]
    index=load(R/'outputs/native_identity_boundary_probe_index_v1.json');summary=load(R/'outputs/native_identity_boundary_probe_summary_v1.json');by_source={v['source']:v['results'] for v in summary['source_results']}
    assert [e['source'] for e in index]==expected_sources==config['selected_sources_by_size_spaced_not_phenotype']
    total=collections.Counter();source_results=[];zero_cases=missing_cases=0
    for entry in index:
        ds=entry['source'];assert sha(entry['metadata'])==entry['SHA256'];m=load(entry['metadata']);p=Path(entry['metadata']).parent/m['arrays'];assert sha(p)==m['arrays_SHA256']
        with np.load(p,allow_pickle=False) as f:d={k:f[k] for k in f.files}
        a=atlas.source(ds);ca=a['ca'];native=ca['native_family_bitmask'];ref=ca['reference_family_bitmask'];major=np.asarray(a['cm']['candidate_dictionaries']['consensus_major'],object)[ca['consensus_major_code']];eligible=(qualification(a,view,core)==1)&ca['fit_pool_candidate']
        frames=[eligible&(native==5),eligible&(native==132),eligible&(native==160),eligible&(major=='T_NK'),eligible&(native==1)&(ref==1)&(major=='T'),eligible&(native==4)&(ref==4)&(major=='Myeloid'),eligible&(native==128)&(ref==128)&(major=='Epithelial')]
        positions=[];codes=[]
        for j,mask in enumerate(frames):
            frame=np.flatnonzero(mask);pick=spread(frame,32);positions.extend(pick.tolist());codes.extend([j]*len(pick));assert len(frame)==m['strata'][j]['n_source_qualified_fit_records_in_frame'] and len(pick)==m['strata'][j]['n_sampled']
        np.testing.assert_array_equal(d['source_native_position'],positions);np.testing.assert_array_equal(d['stratum_code'],codes);np.testing.assert_array_equal(d['native_row'],a['rows'][positions]);assert len(positions)==len(set(positions))
        np.testing.assert_array_equal(d['native_family_bitmask'],native[positions]);np.testing.assert_array_equal(d['reference_family_bitmask'],ref[positions])
        genes=m['genes'];lookup=atlas.native.feature_lookup(ds);physical=[next(iter(lookup[g])) if len(lookup.get(g,set()))==1 else None for g in genes];measured=np.asarray([v is not None for v in physical]);np.testing.assert_array_equal(d['measured'],measured)
        key=atlas.native.legacy_config['counts_key'] if atlas.native.is_legacy(ds) else atlas.native.primary.partitions[ds]['counts_key'];csr=atlas.native.store(ds)[key];shape=tuple(csr.attrs['shape']);expected=np.zeros_like(d['native_gene_values']);totals=np.empty(len(positions));nnz=0
        # Direct per-row indptr/data/indices, without production CSR helpers.
        for i,row in enumerate(d['native_row']):
            lo,hi=int(csr['indptr'][int(row)]),int(csr['indptr'][int(row)+1]);values=np.asarray(csr['data'][lo:hi],dtype=np.float64);columns=np.asarray(csr['indices'][lo:hi],dtype=np.int64);assert np.isfinite(values).all() and (values>=0).all();nnz+=len(values);totals[i]=values.sum()
            per_gene=np.bincount(columns,weights=values,minlength=shape[1])
            for j,col in enumerate(physical):
                if col is not None:expected[i,j]=per_gene[col]
        np.testing.assert_array_equal(expected,d['native_gene_values']);np.testing.assert_allclose(totals,d['native_total_RNA'],rtol=1e-10,atol=1e-6);np.testing.assert_allclose(totals,a['ea']['native_total_RNA'][positions],rtol=1e-10,atol=1e-6)
        assert (totals>0).all();result=[]
        for j,stratum in enumerate(m['strata']):
            mask=d['stratum_code']==j;out={**stratum,'native_total_RNA_quantiles':np.quantile(totals[mask],[0,.25,.5,.75,1]).tolist() if mask.any() else None,'panels':{}}
            for name,gs in m['panels'].items():
                cols=[genes.index(g) for g in gs if measured[genes.index(g)]];x=expected[:,cols];n1=np.count_nonzero(x>0,axis=1);n2=np.count_nonzero(x>=2,axis=1);sums=x.sum(axis=1)
                for suffix,values in [('n_genes_detected',n1),('n_genes_value_ge2',n2),('count_slot_sum',sums),('fraction_of_native_RNA',sums/totals)]:np.testing.assert_allclose(d[name+'__'+suffix],values,rtol=1e-12,atol=1e-12)
                out['panels'][name]={'measured_genes':[g for g in gs if measured[genes.index(g)]],'n_measured':len(cols),'n_requested':len(gs),'fraction_with_at_least2_detected_genes':float(np.mean(n1[mask]>=2)) if mask.any() and len(cols)>=2 else None,'fraction_with_at_least2_genes_each_value_ge2':float(np.mean(n2[mask]>=2)) if mask.any() and len(cols)>=2 else None,'median_native_RNA_fraction':float(np.median((sums/totals)[mask])) if mask.any() and cols else None}
            result.append(out)
        assert result==m['summaries']==by_source[ds];zero_cases+=int(np.any(expected[:,measured]==0));missing_cases+=int((~measured).any())
        total.update(n_sources=1,n_sampled_rows=len(positions),n_native_gene_values_recomputed=int(measured.sum())*len(positions),n_native_nnz_read=nnz);source_results.append({'source':ds,'n_complete_original_rows_recomputed':len(positions),'n_native_gene_values_recomputed':int(measured.sum())*len(positions)})
        save(root/'checkpoints'/(ds+'.json'),{'source':ds,'n_completed_sources':len(source_results)});atlas._cache.clear();atlas.native._fields.clear();del a,ca,d,expected;gc.collect()
    assert total['n_sources']==t['counts']['n_sources'] and total['n_sampled_rows']==t['counts']['n_sampled_rows'];assert zero_cases>0
    out={'status':'ALL_BOUNDARY_PROBE_ROWS_NATIVE_COUNTS_MEASUREMENTS_SAMPLING_AND_SUMMARIES_INDEPENDENTLY_ACCEPTED','totals':dict(total),'sources':source_results,'n_sources_with_measured_zero_cases':zero_cases,'n_sources_with_unmeasured_gene_cases':missing_cases,'missingness_not_invented_if_no_missing_gene_in_this_probe':True,'producer_functions_or_CSR_helpers_imported':False,'native_read_scope':'every row of bounded probe only, not complete archive','not_population_error_rate_or_biological_identity_validation':True,'no_labels_or_QC_changed':True,'producer_terminal_SHA256':sha(R/'provenance/terminal_v1.json')};save(root/'outputs/independent_native_identity_boundary_probe_audit_v1.json',out)
    save(root/'provenance/terminal_v1.json',{'status':out['status'],'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'output_SHA256':sha(root/'outputs/independent_native_identity_boundary_probe_audit_v1.json'),'wall_seconds':time.monotonic()-start,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except BaseException as e:
        import traceback
        save(root/'provenance/failure_v1.json',{'error_class':type(e).__name__,'error':str(e),'traceback':traceback.format_exc()});raise
