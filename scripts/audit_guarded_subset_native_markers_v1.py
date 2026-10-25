"""Independent exact-membership, native-statistic and marker audit.

Checks every saved subset and one complete native group per source; additionally
replays all native rows of a bounded source. This is computational verification,
not independent biological confirmation of the candidate identities.
"""
import argparse,datetime,gc,json,resource,subprocess,sys,time
from pathlib import Path
import numpy as np
from annotation_representation_pilot_v1 import save
from audit_full_source_local_graphs_v1 import index_at,direct_native_rows
from audit_full_source_local_marker_evidence_v1 import verify_markers,replay_complete_source
from extract_full_source_local_marker_evidence_v1 import storage_gate
from full_context_atlas_v1 import sha256_file
from native_archive_collection_v1 import NativeArchiveCollection,row_axis_sha256,text_axis_sha256
from summarize_full_native_annotation_evidence_v1 import read_source


def main(root):
    began=time.monotonic();spec=json.loads((root/'config/guarded_marker_audit_spec_v1.json').read_text());producer=Path(spec['producer']);pc=producer/'config/native_marker_evidence_spec_v1.json'
    if sha256_file(pc)!=spec['producer_config_sha256']:raise ValueError('Guarded marker producer configuration changed')
    cfg=json.loads(pc.read_text());request=json.loads((producer/'config/supervisor_request_v1.json').read_text())
    for name,digest in request['script_SHA256'].items():
        if sha256_file(producer/name)!=digest:raise ValueError('Guarded marker deployed source changed')
    for b in cfg['bindings'].values():
        if sha256_file(b['path'])!=b['sha256']:raise ValueError('Guarded marker bound input changed')
    r=subprocess.run([sys.executable,'-B',str(root/'tests/test_native_marker_unit_support_v1.py')],capture_output=True,text=True,timeout=120)
    (root/'provenance/test_stderr.txt').write_text(r.stderr)
    if r.returncode or 'Ran 3 tests' not in r.stderr or 'skipped' in r.stderr.lower():raise ValueError('Actual independent marker unit numerical tests failed')
    save(root/'provenance/numerical_tests_v1.json',{'actual_exit_code':0,'n_tests':3,'no_skips':True})
    full=spec['scope']=='full_completed_archive'
    if full:records=list(index_at(producer,'full_source_local_marker_evidence_index_v1.json').values())
    else:records=spec['completed_sources']
    candidates,graphs,reference=map(Path,[cfg['candidates'],cfg['graphs'],cfg['reference']])
    ci=index_at(candidates,'full_native_cell_identity_candidate_index_v1.json');gi=index_at(graphs,'full_native_source_local_graph_index_v1.json');ri=index_at(reference,'full_archive_reference_index_v1.json');reader=NativeArchiveCollection(cfg['archive_manifest'])
    observed={r['source'] for r in records}
    if not observed or len(observed)!=len(records) or not observed<=set(ci) or (full and observed!=set(ci)):raise ValueError('Guarded marker audit source scope incomplete or duplicated')
    reports=[]
    for record in records:
        storage_gate(root,spec);ds=record['source'];mp=producer/'outputs'/record['metadata']
        if sha256_file(mp)!=record['metadata_sha256']:raise ValueError('Own-subset marker metadata changed')
        meta=json.loads(mp.read_text());cm,ca=read_source(candidates,ci[ds]);gm,ga=read_source(graphs,gi[ds]);_,rv=read_source(reference,ri[ds]);sm,subset=read_source(producer,meta['source_graph_record']);rows=reader.rows(ds);n=len(rows)
        for a in [ca,ga,rv,subset]:np.testing.assert_array_equal(a['native_row'],rows)
        if meta['native_ids_sha256']!=text_axis_sha256(reader.native_ids(ds)) or meta['native_rows_sha256']!=row_axis_sha256(rows):raise ValueError('Own marker native identity changed')
        fit=ga['fit_pool_candidate'];parent=ga['source_local_group_code'];types=ca['consensus_celltype_code'];names=cm['dictionaries']['consensus_celltype']
        np.testing.assert_array_equal(ca['parent_source_local_group_code'],parent);np.testing.assert_array_equal(ca['fit_pool_candidate'],fit)
        pairs=sorted(set(zip(parent[fit].tolist(),types[fit].tolist())));lookup={pair:j for j,pair in enumerate(pairs)}
        expected=np.full(n,-1,np.int32)
        for i in np.flatnonzero(fit):expected[i]=lookup[(int(parent[i]),int(types[i]))]
        np.testing.assert_array_equal(subset['source_local_group_code'],expected)
        for key in ['fit_pool_candidate','fine_graph_fitted','provisional_family_code']:np.testing.assert_array_equal(subset[key],ga[key])
        if len(sm['guarded_subsets'])!=len(pairs) or sm['partition_kind']!='exact_frozen_parent_by_guarded_celltype_candidate' or sm['guarded_subsets']!=meta['guarded_subsets']:raise ValueError('Subset partition definitions differ')
        if meta['n_archive_rows']!=n or meta['n_fitted_rows']!=int(fit.sum()) or meta['n_unfitted_archive_rows_retained']!=int((~fit).sum()):raise ValueError('Own subset denominator changed')
        if not meta['all_subset_markers_recomputed_from_exact_membership'] or meta['parent_marker_copied_into_subset'] or meta['new_final_identity'] or meta['inferential_DEG_computed'] or meta['parent_group_markers_borrowed'] or meta['patient_independence_claimed']:raise ValueError('Own-subset evidence scope changed')
        original={g['source_local_group_code']:g for g in cm['parent_group_guarded_counts']};split=[]
        for group,((p,t),desc) in enumerate(zip(pairs,sm['guarded_subsets'])):
            member=expected==group;whole=parent==p;identical=np.array_equal(member,whole)
            if desc['guarded_subset_code']!=group or desc['parent_source_local_group_code']!=p or desc['guarded_celltype_candidate']!=names[t] or desc['n_subset_fit_rows']!=int(member.sum()) or desc['n_parent_fit_rows']!=int(whole.sum()) or desc['membership_identical_to_parent']!=identical or desc['native_subset_rows_sha256']!=row_axis_sha256(rows[member]) or desc['stable_parent_group_sha256']!=original[p]['stable_parent_group_sha256']:raise ValueError('Subset parent/membership crosswalk differs')
            split.append(not identical)
        stats=[]
        for name,digest in [('group_native_statistics','group_native_statistics_sha256'),('marker_unit_arrays','marker_unit_arrays_sha256')]:
            p=producer/'outputs'/meta[name]
            if sha256_file(p)!=meta[digest]:raise ValueError('Own marker statistic bytes changed')
            with np.load(p,allow_pickle=False) as f:stats.append({k:f[k] for k in f.files})
        stat,units=stats
        np.testing.assert_array_equal(stat['n'],np.bincount(expected[fit],minlength=len(pairs)))
        genes=meta['native_measured_gene_names'];features=reader.feature_lookup(ds)
        if len(stat['native_columns'])!=len(genes) or any(features.get(g,set())!={int(col)} for g,col in zip(genes,stat['native_columns'])):raise ValueError('Own marker measured gene axes differ')
        count=verify_markers(meta,stat,units)
        # Independent native reads and dense sums for one complete exact subset
        # per source, prioritizing changed membership within a fixed size cap.
        choices=[j for j,size in enumerate(stat['n']) if 0<size<=spec['maximum_group_replay_rows']]
        group=min(choices,key=lambda j:(not split[j],-int(stat['n'][j]),j)) if choices else int(np.argmin(stat['n']))
        selected=np.flatnonzero(expected==group)
        store=reader.store(ds);key=reader.legacy_config['counts_key'] if reader.is_legacy(ds) else reader.primary.partitions[ds]['counts_key']
        direct={name:np.zeros(len(genes),float) for name in ['cp_sum','log_sum','detected']}
        for low in range(0,len(selected),512):
            take=selected[low:low+512];raw=direct_native_rows(store[key],rows[take]);total=np.asarray(raw.sum(axis=1,dtype=np.float64)).ravel()
            np.testing.assert_array_equal(total,rv['native_total_RNA'][take]);x=raw[:,stat['native_columns']].astype(float).toarray();cp=x*(10000./total[:,None])
            direct['cp_sum']+=cp.sum(axis=0);direct['log_sum']+=np.log1p(cp).sum(axis=0);direct['detected']+=(x>0).sum(axis=0)
        error=0.
        for field,value in direct.items():
            np.testing.assert_allclose(value,stat[field][group],rtol=2e-10,atol=2e-7);error=max(error,float(np.max(np.abs(value-stat[field][group]))))
        reports.append({'source':ds,'n_archive_rows':n,'n_guarded_subsets':len(pairs),'n_markers_arithmetic_checked':count,'n_split_subsets':sum(split),'complete_group_native_replay_code':group,'n_complete_group_rows_replayed':len(selected),'maximum_complete_group_stat_absolute_error':error,'source_marker_metadata_sha256':record['metadata_sha256']})
        save(root/'checkpoints'/(ds+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_complete_sources':len(reports),**reports[-1]});reader._fields.clear();del ca,ga,rv,subset,meta,stat,units,raw,x,cp;gc.collect()
    bounded=[r for r in records if r['n_fitted_rows']<=spec['maximum_full_replay_rows']]
    if not bounded:raise ValueError('No bounded complete native source for marker replay')
    chosen=min(bounded,key=lambda r:(r['n_fitted_rows'],r['source']));ds=chosen['source'];meta=json.loads((producer/'outputs'/chosen['metadata']).read_text());_,graph=read_source(producer,meta['source_graph_record']);_,rv=read_source(reference,ri[ds])
    with np.load(producer/'outputs'/meta['group_native_statistics'],allow_pickle=False) as f:stat={k:f[k] for k in f.files}
    with np.load(producer/'outputs'/meta['marker_unit_arrays'],allow_pickle=False) as f:units={k:f[k] for k in f.files}
    replay=replay_complete_source(reader,ds,graph,rv,meta,stat,units)
    out=root/'outputs/guarded_subset_native_marker_audit_v1.json';save(out,{'scope':spec['scope'],'sources':reports,'complete_source_native_replay':replay,'all_subsets_memberships_and_marker_arithmetic_checked':True,'all_native_counts_recomputed':False,'independent_biological_validation':False,'new_final_identity':False})
    terminal={'status':'ALL_GUARDED_SUBSET_MARKER_MEMBERSHIP_AND_ARITHMETIC_AUDITED' if full else 'GUARDED_SUBSET_MARKER_COMPLETED_SOURCE_PREFLIGHT_ONLY','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_sources':len(reports),'n_guarded_subsets':sum(r['n_guarded_subsets'] for r in reports),'n_markers_checked':sum(r['n_markers_arithmetic_checked'] for r in reports),'output_sha256':sha256_file(out),'scope':spec['scope'],'wall_seconds':time.monotonic()-began,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,'full_biological_identity_accepted':False}
    if full:terminal['producer_terminal_sha256']=sha256_file(producer/'provenance/terminal_v1.json')
    save(root/'provenance/terminal_v1.json',terminal)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)});raise
