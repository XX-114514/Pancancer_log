"""Independently replay every local identity-review evidence join.

This audits joins and rule application; it does not establish biological truth
or promote source-local proposals to final common cell identities.
"""
import argparse
import collections
import datetime
import gc
import hashlib
import json
from pathlib import Path
import resource
import time

import numpy as np
from annotation_representation_pilot_v1 import save
from audit_full_source_local_graphs_v1 import index_at
from audit_native_author_label_evidence_v1 import author_context
from full_context_atlas_v1 import sha256_file
from native_archive_collection_v1 import NativeArchiveCollection,row_axis_sha256,text_axis_sha256
from native_annotation_review_v1 import classify_group
from native_lineage_routing_v1 import FAMILIES
from summarize_full_native_annotation_evidence_v1 import read_source


def counts(values):return {str(k):int(v) for k,v in collections.Counter(np.asarray(values).astype(str)).items()}


def main(root):
    began=time.monotonic();spec=json.loads((root/'config/identity_review_audit_spec_v1.json').read_text());producer=Path(spec['producer'])
    cp=producer/'config/local_identity_review_spec_v1.json'
    if sha256_file(cp)!=spec['producer_config_sha256']:raise ValueError('Identity-review config changed')
    cfg=json.loads(cp.read_text());request=json.loads((producer/'config/supervisor_request_v1.json').read_text())
    for name,h in request['script_SHA256'].items():
        if sha256_file(producer/name)!=h:raise ValueError('Identity-review deployed code changed')
    for binding in cfg['bindings'].values():
        if sha256_file(binding['path'])!=binding['sha256']:raise ValueError('Identity-review input binding changed')
    index=index_at(producer,'full_local_identity_review_index_v1.json');graphroot=Path(cfg['graphs']);refroot=Path(cfg['reference']);suproot=Path(cfg['supplemental']);markerroot=Path(cfg['markers'])
    graphs=index_at(graphroot,'full_native_source_local_graph_index_v1.json');refs=index_at(refroot,'full_archive_reference_index_v1.json');supp=index_at(suproot,'supplemental_native_index_v1.json')
    markers=index_at(markerroot,'full_source_local_marker_evidence_index_v1.json');reader=NativeArchiveCollection(cfg['archive_manifest'])
    if any(set(i)!=set(reader.archive_partitions) for i in [index,graphs,refs,supp,markers]):raise ValueError('Full identity review source scope differs')
    models=json.loads((refroot/'outputs/reference_models_v1.json').read_text());reports=[];total_states=0
    for ds in sorted(index):
        record=index[ds];op=producer/'outputs'/record['output']
        if sha256_file(op)!=record['output_sha256']:raise ValueError('Identity review output changed')
        out=json.loads(op.read_text());gm,ga=read_source(graphroot,graphs[ds]);rm,rv=read_source(refroot,refs[ds]);sm,sv=read_source(suproot,supp[ds])
        mp=markerroot/'outputs'/markers[ds]['metadata']
        if sha256_file(mp)!=markers[ds]['metadata_sha256']:raise ValueError('Identity marker binding differs')
        mm=json.loads(mp.read_text());sp=markerroot/'outputs'/mm['group_native_statistics']
        if sha256_file(sp)!=mm['group_native_statistics_sha256']:raise ValueError('Identity marker statistic binding differs')
        with np.load(sp,allow_pickle=False) as f:stat={k:f[k] for k in f.files}
        rows=reader.rows(ds);selected=ga['fit_pool_candidate'];codes=ga['source_local_group_code'];ng=gm['n_source_local_groups']
        if out['source']!=ds or out['native_ids_sha256']!=text_axis_sha256(reader.native_ids(ds)) or out['n_archive_rows']!=len(rows) or out['n_fit_rows']!=int(selected.sum()) or out['n_unfit_archive_rows_retained']!=int((~selected).sum()) or out['n_local_groups']!=ng:raise ValueError('Identity source denominator differs')
        if out['source_graph_record']!=graphs[ds] or out['source_marker_record']['metadata_sha256']!=markers[ds]['metadata_sha256']:raise ValueError('Identity upstream snapshot differs')
        for side in [ga,rv,sv]:np.testing.assert_array_equal(rows,side['native_row'])
        field=rm['author_field'];author=reader.field(ds,field).astype(str) if field else np.full(len(rows),'',dtype='<U1')
        names,inv=np.unique(author,return_inverse=True);context=np.asarray([author_context(v) for v in names])[inv]
        labels=reader.frozen_labels(ds);units=reader.frozen_expression_units(ds).astype(str)
        if len(out['groups'])!=ng or len({g['stable_local_group_sha256'] for g in out['groups']})!=ng:raise ValueError('Identity groups missing or duplicate')
        genes=mm['native_measured_gene_names'];gene_index={g:j for j,g in enumerate(genes)}
        for code,group in enumerate(out['groups']):
            member=(codes==code)&selected;n=int(member.sum());members=rows[member]
            if group['source_local_group_code']!=code or group['n_fit_cells']!=n or n!=stat['n'][code]:raise ValueError('Identity group membership differs')
            stable=hashlib.sha256((cfg['local_label_version']+'\0'+ds+'\0'+row_axis_sha256(members)).encode()).hexdigest()
            expected_id=ds+'::'+cfg['local_label_version']+'::'+str(code+1)
            if group['stable_local_group_sha256']!=stable or group['local_subcluster_id']!=expected_id or group['native_members_row_sha256']!=row_axis_sha256(members):raise ValueError('Local stable identity differs')
            if group['n_expression_units_not_assumed_patients']!=len(set(units[member])):raise ValueError('Identity expression-unit coverage differs')
            if group['native_author_label_counts']!=counts(author[member]) or group['author_coarse_context_counts']!=counts(context[member]):raise ValueError('Identity original author join differs')
            if group['inherited_major_counts']!=counts(labels['major'][member]) or group['inherited_fine_counts']!=counts(labels['fine'][member]):raise ValueError('Identity inherited-label crosswalk differs')
            family=stat['provisional_family_code'][code];expected_route=FAMILIES[family] if family>=0 else 'unresolved_context'
            if group['provisional_routing_family']!=expected_route:raise ValueError('Identity provisional route differs')
            panel={name:float(np.count_nonzero(rv['native_panel_bitmask'][member]&(1<<j))/n) for j,name in enumerate(rm['native_panel_bit_order'])}
            family_fraction={name:float(np.count_nonzero(rv['native_family_bitmask'][member]&(1<<j))/n) for j,name in enumerate(rm['native_family_bit_order'])}
            supplemental={name:float(sv['cue_support'][member,j].mean()) for j,name in enumerate(sm['panel_order'])}
            detection={g:float(stat['detected'][code,gene_index[g]]/n) for g in rm['review_genes'] if g in gene_index}
            if panel!=group['native_panel_fraction'] or family_fraction!=group['native_family_fraction'] or supplemental!=group['supplemental_native_multigene_fraction'] or detection!=group['native_review_gene_detection']:raise ValueError('Identity native evidence join differs')
            if group['supplemental_lineage_review_required']!=any(v>=.2 for v in supplemental.values()):raise ValueError('Supplemental lineage review scope differs')
            for model in rm['source_reference_models']:
                positive=rv[model+'__n_positive_labels'][member];unique=rv[model+'__unique_code'][member];valid=unique>=0
                scope=rv['immune_reference_scope'][member] if model.startswith('immune_') else rv['organ_reference_scope'][member]
                expected=counts(np.asarray(models[model]['labels'])[unique[valid]])
                observed=group['reference_proposals'][model]
                if observed['unique_label_counts']!=expected or observed['n_scored']!=int((positive>=0).sum()) or observed['n_in_scope']!=int(scope.sum()) or observed['n_scored_multiple']!=int((positive>1).sum()) or observed['n_scored_no_positive']!=int((positive==0).sum()):raise ValueError('Identity scoped-reference join differs')
            proposal=classify_group(panel,detection,group['reference_proposals']['immune_low']['unique_label_counts'],n)
            for actual,expected in [('local_major_proposal','local_major'),('local_celltype_proposal','local_celltype'),('proposal_evidence_status','annotation_evidence_status'),('proposal_confidence','confidence'),('native_supported_families_in_legacy_rule','native_supported_families')]:
                if group[actual]!=proposal[expected]:raise ValueError('Identity pilot rule application differs')
            markerinfo=mm['groups'][code];marker=markerinfo['own_group_markers'][0]['gene'] if markerinfo['own_group_markers'] else None
            display=proposal['local_celltype']+'_'+(marker if marker else 'marker_unresolved')+'_'+str(code+1)
            if group['display_marker_proposal']!=marker or group['display_subcluster_proposal']!=display or group['marker_origin']!='own_source_local_group' or group['marker_evidence_type']!='own_group_expression_specificity_not_inferential_DEG' or group['native_marker_evidence_metadata_sha256']!=markers[ds]['metadata_sha256']:raise ValueError('Identity display-marker origin differs')
            for j,name in enumerate(rm['state_order']):
                values=rv['state_native_logCP10K'][member,j];values=values[np.isfinite(values)].astype(np.float64);observed=group['descriptive_RNA_states'][name]
                if observed['n_measured']!=len(values) or observed['measurement_coverage']!=rm['state_measurement_coverage'][name] or observed['functional_identity_inferred']:raise ValueError('Identity state measurement scope differs')
                if len(values):np.testing.assert_allclose(observed['mean_native_logCP10K'],values.mean(),rtol=1e-12,atol=1e-12)
                elif observed['mean_native_logCP10K'] is not None:raise ValueError('Missing RNA state assigned numeric value')
                total_states+=1
            if group['new_authoritative_identity'] or group['malignancy_from_lineage_inferred'] or not group['cross_source_common_identity_not_assigned']:raise ValueError('Identity review output scope promoted')
        reports.append({'source':ds,'n_archive_rows':len(rows),'n_fit_rows':int(selected.sum()),'n_groups':ng,'output_sha256':record['output_sha256']})
        save(root/'checkpoints'/(ds+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':ds,'completed_sources':len(reports)})
        reader._fields.clear();del stat,ga,rv,sv,author,context,labels,units,out,mm;gc.collect()
    output=root/'outputs/full_identity_review_join_audit_v1.json';save(output,{'sources':reports,'n_group_state_summaries_replayed':total_states,
      'cell_identity_biological_truth_established':False,'not_new_native_count_pass':True,'new_final_identity':False})
    save(root/'provenance/terminal_v1.json',{'status':'FULL_IDENTITY_REVIEW_JOIN_AUDIT_COMPLETE_NOT_FINAL_ANNOTATION','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
      'n_sources':len(reports),'n_archive_rows':sum(r['n_archive_rows'] for r in reports),'n_local_groups':sum(r['n_groups'] for r in reports),
      'n_group_state_summaries_replayed':total_states,'output_sha256':sha256_file(output),'producer_terminal_sha256':sha256_file(producer/'provenance/terminal_v1.json'),
      'wall_seconds':time.monotonic()-began,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,'new_final_identity':False})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)});raise
