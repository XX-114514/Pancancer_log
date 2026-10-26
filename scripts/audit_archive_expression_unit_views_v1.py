"""Full native unit/block axes and inherited review audit without expression."""
import argparse,collections,datetime,gc,hashlib,json,resource,time
from pathlib import Path
import numpy as np
import pandas as pd
from annotation_representation_pilot_v1 import save
from audit_full_source_local_graphs_v1 import index_at
from extract_full_source_local_marker_evidence_v1 import storage_gate
from full_context_atlas_v1 import sha256_file
from native_archive_collection_v1 import NativeArchiveCollection,text_axis_sha256
from summarize_full_native_annotation_evidence_v1 import read_source


def main(root):
    began=time.monotonic();cfg=json.loads((root/'config/unit_view_audit_spec_v1.json').read_text());producer=Path(cfg['producer']);index=index_at(producer,'archive_unit_view_index_v1.json');terminal=json.loads((producer/'provenance/terminal_v1.json').read_text());spec=json.loads((producer/'config/archive_unit_view_spec_v1.json').read_text())
    for b in [*cfg['bindings'].values(),*spec['bindings'].values()]:
        if sha256_file(b['path'])!=b['sha256']:raise ValueError('Unit view audit input binding changed')
    reader=NativeArchiveCollection(spec['archive_manifest']);candidate=Path(spec['candidates']);ci=index_at(candidate,'full_native_cell_identity_candidate_index_v1.json');frozen=pd.read_csv(reader.legacy_config['sample_metadata_path'],sep='\t',dtype=str,keep_default_na=False)
    if not frozen.sampleID.is_unique:raise ValueError('Frozen unit metadata not unique')
    old=frozen.set_index('sampleID').to_dict('index');reports=[];view_totals=collections.defaultdict(collections.Counter);reviewed_total=0
    if set(index)!=set(ci) or set(index)!=set(reader.archive_partitions):raise ValueError('Unit view audit all-source scope differs')
    for ds,record in sorted(index.items()):
        storage_gate(root,cfg);mp=producer/'outputs'/record['metadata'];ap=producer/'outputs'/record['arrays']
        if sha256_file(mp)!=record['metadata_sha256'] or sha256_file(ap)!=record['arrays_sha256']:raise ValueError('Unit view output changed')
        m=json.loads(mp.read_text())
        with np.load(ap,allow_pickle=False) as f:a={k:f[k] for k in f.files}
        _,ca=read_source(candidate,ci[ds]);rows=reader.rows(ds);native=reader.frozen_expression_units(ds).astype(str);names=a['expression_unit_dictionary'];unit=a['native_expression_unit_code'];fit=ca['fit_pool_candidate']
        np.testing.assert_array_equal(a['native_row'],rows);np.testing.assert_array_equal(a['fit_pool_candidate'],fit)
        if len(names)!=len(set(names)) or len(unit)!=len(rows) or np.any(unit<0) or np.any(unit>=len(names)):raise ValueError('Unit axis code invalid')
        np.testing.assert_array_equal(names[unit],native)
        if m['native_ids_sha256']!=text_axis_sha256(reader.native_ids(ds)) or m['native_unit_axis_sha256']!=text_axis_sha256(native):raise ValueError('Native cell/unit axis digest differs')
        if a['unit_view_status_code'].shape!=(len(names),len(m['view_order'])) or not np.isin(a['unit_view_status_code'],[0,1,2]).all():raise ValueError('Unit view status array shape/domain invalid')
        if len(m['unit_records'])!=len(names) or record['n_units']!=len(names):raise ValueError('Unit table denominator differs')
        decoded_blocks=m['analysis_block_dictionary'];block_lookup={(v['kind'],v['original_key']):j for j,v in enumerate(decoded_blocks)}
        if len(block_lookup)!=len(decoded_blocks):raise ValueError('Duplicate analysis block dictionary entry')
        expected_block=np.full(len(names),-1,np.int32);unit_counts=np.bincount(unit,minlength=len(names));fit_counts=np.bincount(unit[fit],minlength=len(names));source_reviewed=0
        for j,r in enumerate(m['unit_records']):
            name=str(names[j]);bound=old.get(name) if reader.is_legacy(ds) else None
            if bound is not None and bound['gse_id']!=ds:raise ValueError('Frozen unit joined across sources')
            reviewed=bound is not None and bound['patient_review_status']=='reviewed_identity' and bound['unit_kind']=='reviewed_patient' and bool(bound['analysis_unit']);source_reviewed+=reviewed
            usable=bound is not None or name.strip().lower() not in {'','nan','none','na','n/a','unknown','<na>','unreviewed','not_available'}
            kind='inherited_reviewed_patient' if reviewed else 'native_expression_unit_fallback' if usable else 'unresolved_unit_key_archive_bucket';key=bound['analysis_unit'] if reviewed else name
            if r['source']!=ds or r['native_expression_unit']!=name or r['n_archive_rows']!=int(unit_counts[j]) or r['n_fit_pool_candidate_rows']!=int(fit_counts[j]) or r['analysis_block_kind']!=kind or r['analysis_block_original_key']!=key or r['inherited_reviewed_patient_block_preserved']!=reviewed or r['native_unit_key_usable_or_frozen_reviewed']!=usable:raise ValueError('Native unit/count/reviewed-block mapping differs')
            digest=hashlib.sha256(json.dumps([ds,name],ensure_ascii=False).encode()).hexdigest()
            if r['native_unit_key_sha256']!=digest:raise ValueError('Stable unit namespace key differs')
            expected=None if bound is None else {k:bound[k] for k in r['inherited_frozen_unit_review']}
            if r['inherited_frozen_unit_review']!=expected:raise ValueError('Frozen review fields changed')
            if (kind,key) not in block_lookup:raise ValueError('Reviewed block dictionary omits expected unit')
            expected_block[j]=block_lookup[(kind,key)]
            technical=r['technical_review']
            if technical['patient_identity']!=('accepted' if reviewed else 'unknown') or technical['expression_unit_join']!=('accepted' if usable else 'unknown'):raise ValueError('Patient/unit qualification promoted or erased')
            # This producer supplies no new clinical/assay/source-alias or final
            # identity evidence. Their unresolved gates must stay explicit.
            for field in ['annotation_release','prior_QC_scope','assay_scope','source_alias_scope','treatment_phase','clinical_endpoint','patient_sample_library_join','repeated_sampling_plan']:
                if technical[field]!='unknown':raise ValueError('Unavailable scientific review gate silently promoted')
            for k,view in enumerate(m['view_order']):
                q=r['qualified_views'][view];status=m['view_status_decoder'][str(int(a['unit_view_status_code'][j,k]))]
                if q['status']!=status or q['eligible']!=(status=='accepted'):raise ValueError('Unit status array/table differs')
                if view=='all_archive':
                    if status!='accepted' or q['reasons']:raise ValueError('Archive member erased by scientific view gate')
                elif status=='accepted' or not q['reasons']:raise ValueError('No scientific view can bypass unresolved required review')
                if view=='whole_tissue_composition' and technical['capture_class'] in ['immune_enriched','T_enriched_nuclei','selected_nonmalignant_subset','cell_product','mixed_sorted_fractions'] and status!='rejected':raise ValueError('Selected capture presented as whole tissue')
                view_totals[view][status+'_units']+=1;view_totals[view][status+'_archive_rows']+=int(unit_counts[j])
        np.testing.assert_array_equal(a['analysis_block_code'],expected_block[unit])
        if m['new_patient_identity_clinical_eligibility_or_independence_assigned']:raise ValueError('Unit table changed independence status')
        reviewed_total+=source_reviewed;reports.append({'source':ds,'n_archive_rows':len(rows),'n_units_not_patients':len(names),'n_preserved_inherited_reviewed_patient_units':source_reviewed})
        save(root/'checkpoints'/(ds+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_complete_sources':len(reports),**reports[-1]});reader._fields.clear();del a,ca,native;gc.collect()
    summary_path=producer/'outputs/archive_unit_view_coverage_v1.json'
    if sha256_file(summary_path)!=terminal['summary_sha256']:raise ValueError('Unit view coverage summary changed')
    summary=json.loads(summary_path.read_text())
    if summary['views']!={k:dict(v) for k,v in view_totals.items()} or summary['totals']['n_inherited_reviewed_patient_units']!=reviewed_total or sum(r['n_archive_rows'] for r in reports)!=terminal['n_archive_rows']:raise ValueError('Whole unit-view coverage totals differ')
    output=root/'outputs/archive_unit_view_audit_v1.json';save(output,{'sources':reports,'n_preserved_reviewed_patient_units_not_patients':reviewed_total,'all_scientific_unknown_gates_retained':True,'new_clinical_qualification':False})
    save(root/'provenance/terminal_v1.json',{'status':'ALL_ARCHIVE_NATIVE_UNIT_AND_INHERITED_REVIEW_BLOCK_AXES_AUDITED_SCIENTIFIC_VIEW_GATES_REMAIN_EXPLICIT','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_sources':len(reports),'n_archive_rows':sum(r['n_archive_rows'] for r in reports),'n_units_not_patients':sum(r['n_units_not_patients'] for r in reports),'output_sha256':sha256_file(output),'producer_terminal_sha256':sha256_file(producer/'provenance/terminal_v1.json'),'native_expression_read':False,'new_clinical_qualification':False,'wall_seconds':time.monotonic()-began,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)});raise
