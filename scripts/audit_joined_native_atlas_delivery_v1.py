"""Audit delivery bindings, every exact group and every model-domain join.

All candidate-cell axes already have full independent upstream audits. This
checks newly joined metadata at every group and retained exclusion category,
and bounded independent native CSR values per source. It is not a new
biological truth set or an exhaustive reread of every count.
"""
import argparse,collections,datetime,gc,json,resource,subprocess,sys,time
from pathlib import Path
import numpy as np
from annotation_representation_pilot_v1 import save
from audit_full_source_local_graphs_v1 import direct_native_rows
from extract_full_source_local_marker_evidence_v1 import storage_gate
from full_context_atlas_v1 import sha256_file
from joined_native_atlas_delivery_v1 import JoinedNativeAtlas,completed


def assert_carried_axis_equal(actual,expected):
    """Preserve missing objects without turning them into the string 'nan'."""
    a,b=np.asarray(actual),np.asarray(expected)
    if a.shape!=b.shape:raise AssertionError('Carried metadata axes differ')
    for x,y in zip(a.ravel(),b.ravel()):
        if x is y:continue
        if isinstance(x,(float,np.floating)) and isinstance(y,(float,np.floating)) and np.isnan(x) and np.isnan(y):continue
        same=x==y
        if not isinstance(same,(bool,np.bool_)) or not same:raise AssertionError('Carried metadata values or missingness differ')


def main(root):
    start=time.monotonic();cfg=json.loads((root/'config/joined_delivery_audit_spec_v1.json').read_text());producer=Path(cfg['producer']);t=completed(producer)
    for b in cfg['bindings'].values():
        if sha256_file(b['path'])!=b['sha256']:raise ValueError('Joined audit bound evidence changed')
    test=subprocess.run([sys.executable,'-B',str(root/'tests/test_joined_metadata_missingness_audit_v1.py')],capture_output=True,text=True,timeout=120)
    (root/'provenance/test_stderr.txt').write_text(test.stderr)
    if test.returncode or 'Ran 4 tests' not in test.stderr or 'skipped' in test.stderr.lower():raise ValueError('Carried metadata missingness counterexamples failed')
    save(root/'provenance/numerical_tests_v1.json',{'actual_exit_code':0,'n_tests':4,'no_skips':True})
    mp=producer/'outputs/joined_native_atlas_delivery_v1.json'
    if sha256_file(mp)!=t['manifest_sha256']:raise ValueError('Delivery manifest changed')
    reader=JoinedNativeAtlas(mp);records=[];total=collections.Counter()
    for ds in reader.sources:
        storage_gate(root,cfg);s=reader.source(ds);ca,ea=s['ca'],s['ea'];n=len(s['rows']);group=ca['global_subcluster_number'];_,representative=np.unique(group,return_index=True)
        # Include the exclusions and each coarse/fine evidence category even
        # when they have no source-local fitted subgroup.
        more=[]
        for field in ['consensus_major_code','consensus_celltype_code','coarse_evidence_code','fine_evidence_code','fit_pool_candidate']:
            _,p=np.unique(ca[field],return_index=True);more.extend(p.tolist())
        p=np.unique(np.r_[representative,more,0,n-1]).astype(np.int64);checked=0
        local={v['source_local_group_code']:v for v in s['lm']['groups']}
        for low in range(0,len(p),4096):
            q=p[low:low+4096][::-1];out=reader.metadata(ds,q);v=out['columns'];np.testing.assert_array_equal(v['native_cell_id'],s['ids'][q])
            np.testing.assert_array_equal(v['native_row'],s['rows'][q]);np.testing.assert_array_equal(v['global_subcluster_number'],ca['global_subcluster_number'][q])
            author=reader.native.author_labels(ds);assert_carried_axis_equal(v['author_label'],author['values'][q])
            for field in ['consensus_major','consensus_celltype','coarse_evidence','fine_evidence']:
                expected=np.asarray(s['cm']['candidate_dictionaries'][field],object)[ca[field+'_code'][q]];np.testing.assert_array_equal(v[field],expected)
            for j,position in enumerate(q):
                code=int(ca['global_subcluster_number'][position]);parent=int(ca['parent_source_local_group_code'][position])
                expected_major=local[parent]['local_major_proposal'] if parent in local else 'unresolved_no_fitted_local_group'
                if v['local_major'][j]!=expected_major:raise ValueError('Local major was silently replaced by consensus')
                if code<0:
                    if any(v[f][j] is not None for f in ['global_subcluster_id','marker','stable_group_sha256']):raise ValueError('Archive exclusion acquired a fabricated subgroup')
                else:
                    g=reader.groups[code]
                    if g['source']!=ds or g['source_subset_code']!=ca['guarded_source_subset_code'][position]:raise ValueError('Global group source/subaxis differs')
                    for field in ['global_subcluster_id','subcluster','marker','marker_origin','marker_evidence_type','stable_group_sha256']:
                        if v[field][j]!=g[field]:raise ValueError('Delivery group annotation provenance differs')
            for j,state in enumerate(s['em']['state_order']):np.testing.assert_array_equal(v['RNA_state__'+state],ea['state_native_logCP10K'][q,j])
            unit=s['ua']['native_expression_unit_code'][q]
            np.testing.assert_array_equal(v['native_expression_unit'],s['ua']['expression_unit_dictionary'][unit])
            for j,u in enumerate(unit):
                old=s['um']['unit_records'][u]['inherited_frozen_unit_review'] or {}
                if v['patient_review_status'][j]!=old.get('patient_review_status','unknown'):raise ValueError('Reviewed patient status changed')
                for k,view in enumerate(s['um']['view_order']):
                    expected=s['um']['view_status_decoder'][str(s['ua']['unit_view_status_code'][u,k])]
                    if v['view_status__'+view][j]!=expected:raise ValueError('Scientific gate was promoted by delivery')
            checked+=len(q)
        native=reader.native;genes=['TRAC','CD3D','LYZ','LST1','EPCAM','COL1A1','__UNMEASURED_SENTINEL__'];probe=np.asarray([n-1,0,n//2,n-1],np.int64);counts=reader.read_counts(ds,probe,genes);direct=direct_native_rows(native.store(ds)[counts['count_origin']],s['rows'][probe]);lookup=native.feature_lookup(ds)
        for j,g in enumerate(genes):
            cols=lookup.get(g,set());measured=len(cols)==1
            if bool(counts['measured'][j])!=measured:raise ValueError('Source measurement mask differs')
            expected=direct[:,next(iter(cols))].toarray().ravel() if measured else np.zeros(len(probe));np.testing.assert_array_equal(counts['counts'][:,j].toarray().ravel(),expected)
        np.testing.assert_array_equal(counts['native_total_counts'],np.asarray(direct.sum(axis=1,dtype=np.float64)).ravel())
        model_count=0;missing_domain_tests=0;empty_count=0;patch_count=0;raw_only_count=0
        for key,entry in reader.pi.items():
            if key[0]!=ds:continue
            rec=entry['base_projection_record'];model_count+=1
            if not rec['n_model_domain_rows']:
                v=reader.projection(ds,key[1],probe)
                if v['model_domain_available'].any() or v['arrays']:raise ValueError('Empty lineage domain was zero-filled')
                empty_count+=1;continue
            path=Path(reader.projections['base_projection'])/'outputs'/rec['arrays']
            with np.load(path,allow_pickle=False) as f:
                domain=f['source_native_position'];positions=domain[np.unique([0,len(domain)//2,len(domain)-1])];absent=np.flatnonzero(~np.isin(np.arange(n),domain))[:1];q=np.r_[positions,absent,positions[:1]][::-1];v=reader.projection(ds,key[1],q);expected=np.isin(q,domain)
                np.testing.assert_array_equal(v['model_domain_available'],expected)
                for field in ['raw_PCA','candidate_Harmony','reference_fit_exposure','adaptation_eligible','adaptation_included']:
                    np.testing.assert_array_equal(v['arrays'][field].mask[~expected],True)
                    idx=np.searchsorted(domain,q[expected]);np.testing.assert_array_equal(v['arrays'][field][expected],f[field][idx])
                if not np.isfinite(f['candidate_Harmony']).any():raw_only_count+=1
                missing_domain_tests+=int(bool(len(absent)))
            for patch in entry['exact_neighbor_overlays']:
                with np.load(patch['repair_arrays'],allow_pickle=False) as f:a={k:f[k] for k in f.files}
                v=reader.projection(ds,key[1],a['source_native_position'][::-1]);basis=patch['basis'];np.testing.assert_allclose(v['arrays'][basis+'_reference_distance'],a['exact_nearest_other_reference_distance'][::-1],rtol=0,atol=0);np.testing.assert_array_equal(v['arrays'][basis+'_outside_reference_radius_q99'],a['outside_original_reference_radius_q99'][::-1]);patch_count+=1
        row={'source':ds,'n_archive_rows':n,'n_exact_groups_checked':int((np.unique(group)>=0).sum()),'n_metadata_representatives_checked':checked,'n_model_domains_checked':model_count,'n_empty_model_domains_checked':empty_count,'n_absent_lineage_row_tests':missing_domain_tests,'n_raw_only_domains_retained':raw_only_count,'n_exact_overlay_domains_checked':patch_count,'n_independent_count_probe_rows':len(probe)};records.append(row);total.update({k:v for k,v in row.items() if k!='source'})
        save(root/'checkpoints'/(ds+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_complete_sources':len(records),**row});reader._cache.clear();native._fields.clear();native._native_arrays.clear();del ca,ea,s,counts,direct;gc.collect()
    if total['n_archive_rows']!=t['n_archive_rows'] or total['n_exact_overlay_domains_checked']!=t['n_exact_overlay_domains_exercised']:raise ValueError('Complete joined scope differs')
    out=root/'outputs/joined_native_atlas_delivery_audit_v1.json';save(out,{'sources':records,'totals':dict(total),'upstream_full_row_audits_bound_not_replaced':True,'new_join_audit_scope':'every_group_and_model_domain_with_category_and_bounded_direct_native_probes','new_authoritative_taxonomy_or_clinical_eligibility':False})
    save(root/'provenance/terminal_v1.json',{'status':'ALL_SOURCE_READONLY_DELIVERY_GROUP_MODEL_AND_NATIVE_PROBE_JOINS_AUDITED','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_sources':len(records),**dict(total),'output_sha256':sha256_file(out),'producer_terminal_sha256':sha256_file(producer/'provenance/terminal_v1.json'),'new_authoritative_taxonomy':False,'standalone_Zarr':False,'wall_seconds':time.monotonic()-start,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except Exception as e:save(root/'provenance/failure_v1.json',{'error_class':type(e).__name__,'error':str(e)});raise
