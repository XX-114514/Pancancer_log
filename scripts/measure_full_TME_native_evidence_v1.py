"""One bounded native pass for all qualified rows and exact-group markers.

All uniquely measured native genes contribute group statistics. A fixed small
context panel retains per-cell amounts for subsequent identity review. Neither
the current annotation nor source/QC/module eligibility changes in this stage.
"""
import argparse,collections,datetime,gc,hashlib,json,resource,subprocess,sys,time
from pathlib import Path
import numpy as np
import psutil
from TME_native_measurement_core_v1 import CONTEXT_GENES,CONTEXT_PANELS,canonical_native,context_values,update_group_stats,expression_candidates

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
    with Path(p).open('x') as f:json.dump(d,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
def family(major):
    if major in ['T','NK','T_NK']:return 'T_NK'
    if major in ['B','Plasma','B_Plasma']:return 'B_Plasma'
    return major

def main(root):
    begin=time.monotonic();cfg=load(root/'config/spec_v1.json')
    for b in cfg['bindings'].values():assert sha(b['path'])==b['SHA256']
    test=subprocess.run([sys.executable,'-B',str(root/'tests/test_TME_native_measurement_core_v1.py')],capture_output=True,text=True,timeout=120)
    (root/'provenance/test_stderr_v1.txt').write_text(test.stderr)
    assert test.returncode==0 and 'Ran 4 tests' in test.stderr and 'skipped' not in test.stderr.lower()
    sys.path.insert(0,str(Path(cfg['joined'])/'scripts'))
    from native_archive_collection_v1 import NativeArchiveCollection,text_axis_sha256
    from native_annotation_review_v1 import native_measured_genes
    from native_panel_aggregate_v1 import gather_csr_rows
    from extract_full_source_local_marker_evidence_v1 import storage_gate
    accepted=load(cfg['working_acceptance']);manifest=load(accepted['producer_manifest']['path'])
    assert sha(accepted['producer_manifest']['path'])==accepted['producer_manifest']['SHA256']
    assert sha(manifest['index']['path'])==manifest['index']['SHA256'];assert sha(manifest['group_catalog']['path'])==manifest['group_catalog']['SHA256']
    source_index=load(manifest['index']['path']);groups=load(manifest['group_catalog']['path'])['groups'];base=Path(manifest['index']['path']).parent
    archive=load(Path(cfg['joined'])/'outputs/joined_native_atlas_delivery_v1.json')['archive_manifest'];reader=NativeArchiveCollection(archive)
    records=[];totals=collections.Counter()
    for item in source_index:
        storage_gate(root,cfg);started=time.monotonic();ds=item['source'];mp=base/item['metadata'];ap=base/item['arrays'];assert sha(mp)==item['metadata_SHA256'] and sha(ap)==item['arrays_SHA256'];m=load(mp)
        with np.load(ap,allow_pickle=False) as f:
            q={k:f[k] for k in ['source_native_position','native_row','annotation_group_code','native_total_RNA','consensus_major_code','consensus_celltype_code','annotation_review_code','fit_pool_candidate_not_final_QC']}
        positions=q['source_native_position'];rows=q['native_row'];n=len(rows);np.testing.assert_array_equal(reader.rows(ds)[positions],rows);assert text_axis_sha256(reader.native_ids(ds))==m['native_ids_source_SHA256']
        source_groups=sorted([g for g in groups if g['source']==ds],key=lambda g:g['source_annotation_group_code']);ng=len(source_groups);codes=q['annotation_group_code']
        assert [g['source_annotation_group_code'] for g in source_groups]==list(range(ng))
        genes,columns=native_measured_genes(reader,ds);p=len(genes);lookup={g:int(c) for g,c in zip(genes,columns)};physical=[lookup.get(g) for g in CONTEXT_GENES]
        stat_bytes=ng*p*8*3+n*len(CONTEXT_GENES)*8
        if stat_bytes+cfg['RAM_reserve_bytes']>psutil.virtual_memory().available:raise MemoryError('Fresh native-evidence RAM gate failed')
        stat={'n':np.zeros(ng,np.int64),**{k:np.zeros((ng,p),np.float64) for k in ['cp_sum','log_sum','detected']}}
        context=np.full((n,len(CONTEXT_GENES)),np.nan,np.float64);measured=np.asarray([j is not None for j in physical],bool)
        key=reader.legacy_config['counts_key'] if reader.is_legacy(ds) else reader.primary.partitions[ds]['counts_key'];csr=reader.store(ds)[key]
        nnz_read=0;canonical_nnz=0
        for low in range(0,n,cfg['read_block_size']):
            high=min(n,low+cfg['read_block_size']);raw=gather_csr_rows(csr,rows[low:high],max_span=cfg['max_native_row_span']);nnz_read+=raw.nnz
            raw=canonical_native(raw);canonical_nnz+=raw.nnz;rna=np.asarray(raw.sum(axis=1,dtype=np.float64)).ravel();np.testing.assert_array_equal(rna,q['native_total_RNA'][low:high])
            values,mask=context_values(raw,physical);np.testing.assert_array_equal(mask,measured);context[low:high]=values
            update_group_stats(stat,raw[:,columns].tocsr(),rna,codes[low:high])
            if low//cfg['read_block_size']%64==0:
                storage_gate(root,cfg)
                (root/'checkpoints/current_source_progress_v1.json').write_text(json.dumps({'observed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':ds,'n_complete_scope_rows_in_source':high,'n_source_scope_rows':n,'n_completed_sources':len(records)})+'\n')
        np.testing.assert_array_equal(stat['n'],np.bincount(codes[codes>=0],minlength=ng));assert stat['n'].sum()==q['fit_pool_candidate_not_final_QC'].sum()
        source_families=[]
        for group in source_groups:
            kinds=list(group['major_distribution']);source_families.append(family(kinds[0]) if len(kinds)==1 else 'mixed_major_review')
        source_families=np.asarray(source_families,object);marker_records=[]
        for j,group in enumerate(source_groups):
            compare=source_families==source_families[j];compare[j]=False;scope='other_current_TME_source_groups_in_same_analysis_family';count=int(stat['n'][compare].sum())
            if count==0:
                compare=np.arange(ng)!=j;count=int(stat['n'][compare].sum());scope='other_current_TME_source_groups_no_within_family_comparator'
            own_cp=stat['cp_sum'][j]/stat['n'][j];own_det=stat['detected'][j]/stat['n'][j]
            other_cp=stat['cp_sum'][compare].sum(axis=0)/count if count else None;other_det=stat['detected'][compare].sum(axis=0)/count if count else None
            candidates=expression_candidates(genes,own_cp,own_det,other_cp,other_det,cfg['minimum_detection'],cfg['marker_limit'])
            marker_records.append({'source_annotation_group_code':j,'global_subcluster_number':group['global_subcluster_number'],'stable_group_sha256':group['stable_group_sha256'],
                'n_exact_group_cells':int(stat['n'][j]),'current_analysis_family':str(source_families[j]),'comparator_scope':scope,'n_comparator_cells':count,
                'own_native_expression_marker_candidates':candidates,'previous_marker':group['marker'],'previous_marker_evidence_type':group['marker_evidence_type'],
                'no_new_marker_selected_if_candidates_empty':not bool(candidates),'marker_candidates_not_replicate_aware_DEG':True,
                'original_group_and_cell_identity_not_changed':True,'expression_unit_and_patient_support_not_newly_inferred':True})
        sp=root/'outputs'/(ds+'_TME_full_gene_group_statistics_v1.npz')
        with sp.open('xb') as f:np.savez_compressed(f,**stat,native_columns=columns)
        cp=root/'outputs'/(ds+'_TME_native_context_amounts_v1.npz')
        with cp.open('xb') as f:np.savez_compressed(f,source_native_position=positions,native_row=rows,native_gene_values=context,measured=measured,native_total_RNA=q['native_total_RNA'],annotation_group_code=codes)
        for path,arrays in [(sp,stat),(cp,{'source_native_position':positions,'native_row':rows,'native_gene_values':context,'measured':measured,'native_total_RNA':q['native_total_RNA'],'annotation_group_code':codes})]:
            with np.load(path,allow_pickle=False) as f:
                for k,v in arrays.items():np.testing.assert_array_equal(f[k],v)
        output={'source':ds,'n_current_scope_records':n,'n_fit_group_records':int(stat['n'].sum()),'n_exact_groups':ng,'native_ids_source_SHA256':m['native_ids_source_SHA256'],
            'scope_annotation_metadata_SHA256':item['metadata_SHA256'],'scope_annotation_arrays_SHA256':item['arrays_SHA256'],
            'native_measured_gene_names':genes,'all_uniquely_measured_native_genes_included':True,'native_columns':columns.tolist(),
            'context_genes':CONTEXT_GENES,'context_panels':CONTEXT_PANELS,'context_gene_measured':measured.tolist(),
            'group_statistics':{'path':str(sp),'SHA256':sha(sp)},'native_context':{'path':str(cp),'SHA256':sha(cp)},
            'groups':marker_records,'native_nnz_read':int(nnz_read),'canonical_native_nnz':int(canonical_nnz),
            'duplicate_sparse_slots_coalesced_in_RAM_before_log_not_written_to_archive':True,'native_count_origin':key,
            'all_context_missing_features_are_NaN_with_explicit_mask':True,'no_inferential_DEG_or_new_labels_or_QC':True,
            'unit_support_for_new_marker_candidates_not_yet_estimated':True,'wall_seconds':time.monotonic()-started}
        op=root/'outputs'/(ds+'_TME_native_evidence_v1.json');save(op,output)
        records.append({'source':ds,'metadata':op.name,'metadata_SHA256':sha(op),'n_current_scope_records':n,'n_exact_groups':ng,'group_statistics_SHA256':sha(sp),'native_context_SHA256':sha(cp)})
        totals.update(n_sources=1,n_current_scope_records=n,n_fit_group_records=int(stat['n'].sum()),n_exact_groups=ng,n_native_nnz_read=int(nnz_read),n_groups_with_own_marker_candidates=sum(bool(g['own_native_expression_marker_candidates']) for g in marker_records))
        save(root/'checkpoints'/('source_'+ds+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'totals':dict(totals),'source_wall_seconds':output['wall_seconds']})
        del q,stat,context,raw,output;reader._fields.clear();gc.collect()
    assert totals['n_current_scope_records']==sum(r['n_source_qualified_rows'] for r in source_index) and totals['n_exact_groups']==len(groups)
    ip=root/'outputs/TME_full_native_evidence_index_v1.json';save(ip,records)
    save(root/'provenance/terminal_v1.json',{'status':'COMPLETE_TME_SCOPE_NATIVE_CONTEXT_AMOUNTS_AND_FULL_GENE_EXACT_GROUP_EXPRESSION_MEASURED_INDEPENDENT_AUDIT_REQUIRED',
        'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'totals':dict(totals),'index_SHA256':sha(ip),
        'n_context_genes':len(CONTEXT_GENES),'marker_ranking_is_descriptive_not_DEG':True,'no_labels_or_QC_or_CM_changed':True,
        'output_bytes':sum(p.stat().st_size for p in root.rglob('*') if p.is_file()),'wall_seconds':time.monotonic()-begin,
        'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);r=Path(p.parse_args().output)
    try:main(r)
    except BaseException:
        import traceback
        save(r/'provenance/failure_v1.json',{'traceback':traceback.format_exc()});raise
