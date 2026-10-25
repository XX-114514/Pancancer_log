"""All-source reference exposure, rare reserves and native feature-axis replay."""
import argparse
import datetime
import gc
import json
from pathlib import Path
import resource
import time

import numpy as np
from annotation_representation_pilot_v1 import save
from audit_full_source_local_graphs_v1 import index_at
from full_context_atlas_v1 import sha256_file
from native_archive_collection_v1 import NativeArchiveCollection,text_axis_sha256
from native_lineage_routing_v1 import FAMILIES,local_family
from summarize_full_native_annotation_evidence_v1 import read_source


def main(root):
    began=time.monotonic();spec=json.loads((root/'config/reference_input_audit_spec_v1.json').read_text());producer=Path(spec['producer'])
    cp=producer/'config/common_reference_input_spec_v1.json'
    if sha256_file(cp)!=spec['producer_config_sha256']:raise ValueError('Reference input configuration changed')
    cfg=json.loads(cp.read_text());request=json.loads((producer/'config/supervisor_request_v1.json').read_text())
    for name,h in request['script_SHA256'].items():
        if sha256_file(producer/name)!=h:raise ValueError('Reference input deployed code changed')
    for binding in cfg['bindings'].values():
        if sha256_file(binding['path'])!=binding['sha256']:raise ValueError('Reference input upstream binding changed')
    index=index_at(producer,'full_common_reference_input_index_v1.json');candidate=Path(cfg['candidate']);ci=index_at(candidate,'full_native_cell_identity_candidate_index_v1.json');native_inputs=Path(cfg['native_local_inputs']);fi=index_at(native_inputs,'full_native_local_discovery_input_index_v1.json')
    reader=NativeArchiveCollection(cfg['archive_manifest']);sources=sorted(index)
    if set(sources)!=set(reader.archive_partitions) or set(sources)!=set(ci) or set(sources)!=set(fi):raise ValueError('Reference audit full scope differs')
    allocation=json.loads((producer/'outputs/reference_pool_allocation_v1.json').read_text())
    if allocation['sources']!=sources or allocation['pools']!=['global_overview',*FAMILIES]:raise ValueError('Reference allocation axes differ')
    featurepath=producer/'outputs/common_measured_model_feature_candidates_v1.json';featuremeta=json.loads(featurepath.read_text());genes=featuremeta['genes']
    if genes!=sorted(set(genes)) or len(genes)<cfg['minimum_common_candidates'] or featuremeta['feature_selection_fitted']:raise ValueError('Common feature candidates invalid')
    fp=Path(cfg['feature_coverage'])/'outputs/full_archive_feature_coverage_v1.json';expected_common=set(json.loads(fp.read_text())['common_nonheldout_feature_candidates'])
    local_features={}
    for ds in sources:
        mp=native_inputs/'outputs'/fi[ds]['metadata']
        if sha256_file(mp)!=fi[ds]['metadata_sha256']:raise ValueError('Source feature metadata changed')
        m=json.loads(mp.read_text());p=native_inputs/'outputs'/m['feature_candidates']
        if sha256_file(p)!=m['feature_candidates_sha256']:raise ValueError('Source feature candidates changed')
        f=json.loads(p.read_text());expected_common&=set(f['genes']);local_features[ds]=f
    if genes!=sorted(expected_common):raise ValueError('Common feature intersection differs')
    modepath=Path(cfg['modality_overlay'])/'outputs/source_modality_overlay_v1.json';modes={r['source']:r for r in json.loads(modepath.read_text())['sources']}
    reports=[];totals={name:0 for name in allocation['pools']}
    for si,ds in enumerate(sources):
        meta,side=read_source(producer,index[ds]);cm,ca=read_source(candidate,ci[ds]);rows=reader.rows(ds);n=len(rows)
        np.testing.assert_array_equal(rows,side['native_row']);np.testing.assert_array_equal(rows,ca['native_row'])
        if meta['native_ids_sha256']!=text_axis_sha256(reader.native_ids(ds)) or meta['candidate_record']!=ci[ds] or meta['source_modality_context']!=modes[ds]:raise ValueError('Reference identity/modality/candidate binding differs')
        eligible=ca['fit_pool_candidate']&ca['positive_native_RNA'];np.testing.assert_array_equal(side['fit_pool_candidate'],eligible)
        np.testing.assert_array_equal(side['parent_source_local_group_code'],ca['parent_source_local_group_code'])
        major=np.asarray(cm['dictionaries']['consensus_major'])[ca['consensus_major_code']];fine=np.asarray(cm['dictionaries']['consensus_celltype'])[ca['consensus_celltype_code']]
        expected_route=np.full(n,-1,np.int16)
        for name in np.unique(major):
            family=local_family(str(name))
            if family is not None:expected_route[major==name]=FAMILIES.index(family)
        np.testing.assert_array_equal(side['guarded_route_code'],expected_route)
        gamma=(ca['own_fine_bitmask']&(1<<cm['fine_candidate_order'].index('T_gamma_delta_RNA_candidate')))!=0
        # Integer triples independently encode the local/fine/gamma stratum.
        fine_names,fine_code=np.unique(fine,return_inverse=True)
        parent=ca['parent_source_local_group_code'];strata=(parent.astype(np.int64)+1)*(2*len(fine_names))+2*fine_code+gamma
        units=reader.frozen_expression_units(ds).astype(str)
        if np.any(side['global_reference_fit']&~eligible) or np.any(side['lineage_reference_fit']&(~eligible|(expected_route<0))):raise ValueError('Reference selection outside eligible route')
        for name in allocation['pools']:
            member=eligible if name=='global_overview' else eligible&(expected_route==FAMILIES.index(name))
            selected=side['global_reference_fit'] if name=='global_overview' else side['lineage_reference_fit']&(expected_route==FAMILIES.index(name))
            _,all_counts=np.unique(strata[member],return_counts=True);labels,selected_counts=np.unique(strata[selected],return_counts=True);all_labels=np.unique(strata[member]);observed=dict(zip(labels,selected_counts))
            reserve=int(np.minimum(all_counts,cfg['rare_reserve']).sum());count=int(selected.sum());capacity=int(member.sum());summary=meta['reference_pools'][name]
            if capacity!=allocation['capacities'][name][si] or reserve!=allocation['minimum_reserves'][name][si] or count!=allocation['allocated'][name][si]:raise ValueError('Reference allocation denominator differs')
            if any(observed.get(label,0)<min(number,cfg['rare_reserve']) for label,number in zip(all_labels,all_counts)):raise ValueError('Reference local rare reserve not preserved')
            expected={'n_fit_pool_rows':capacity,'n_selected_reference_rows':count,'n_local_guarded_gamma_strata':len(all_labels),'n_selected_local_guarded_gamma_strata':len(labels),'minimum_reserved_rows':reserve,'n_selected_expression_units':len(np.unique(units[selected]))}
            if summary!=expected:raise ValueError('Reference pool coverage summary differs')
            totals[name]+=count
        featurefile=producer/'outputs'/meta['native_feature_columns']
        if sha256_file(featurefile)!=meta['native_feature_columns_sha256'] or meta['common_feature_list_sha256']!=sha256_file(featurepath):raise ValueError('Reference native feature mapping changed')
        with np.load(featurefile,allow_pickle=False) as f:
            if f.files!=['native_columns']:raise ValueError('Unexpected native feature mapping keys')
            columns=f['native_columns']
        local=local_features[ds];mapping=dict(zip(local['genes'],local['native_columns']));np.testing.assert_array_equal(columns,[mapping[g] for g in genes]);lookup=reader.feature_lookup(ds)
        if len(set(columns))!=len(columns) or any(lookup.get(g,set())!={int(c)} for g,c in zip(genes,columns)) or set(columns)&set(local['heldout_physical_alias_columns_excluded']):raise ValueError('Reference physical feature/alias axis differs')
        if not meta['local_fine_graphs_and_identities_not_replaced'] or not meta['modality_and_anatomy_transfer_not_accepted_by_sampling'] or meta['new_final_identity']:raise ValueError('Reference input scope promoted')
        reports.append({'source':ds,'n_archive_rows':n,'n_global_reference_rows':int(side['global_reference_fit'].sum()),'n_lineage_reference_rows':int(side['lineage_reference_fit'].sum()),'arrays_sha256':index[ds]['arrays_sha256'],'metadata_sha256':index[ds]['metadata_sha256']})
        save(root/'checkpoints'/(ds+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':ds,'completed_sources':len(reports)})
        reader._fields.clear();del side,ca,eligible,major,fine,expected_route,gamma,strata,units;gc.collect()
    for name,value in totals.items():
        cap=cfg['global_reference_cap'] if name=='global_overview' else cfg['lineage_reference_cap']
        if value!=min(cap,sum(allocation['capacities'][name])):raise ValueError('Reference total cap differs')
    output=root/'outputs/common_reference_input_audit_v1.json';save(output,{'sources':reports,'reference_pool_totals':totals,'n_common_feature_candidates':len(genes),'all_archive_axes_exposure_rare_strata_and_feature_columns_replayed':True,'coordinates_or_final_labels_accepted':False})
    save(root/'provenance/terminal_v1.json',{'status':'COMMON_REFERENCE_ROWS_FEATURE_AXES_AUDIT_COMPLETE_NO_MODEL_FITTED','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
      'n_sources':len(reports),'n_archive_rows':sum(r['n_archive_rows'] for r in reports),'n_global_reference_rows':totals['global_overview'],
      'n_lineage_reference_rows':sum(totals[name] for name in FAMILIES),'n_common_feature_candidates':len(genes),'output_sha256':sha256_file(output),
      'producer_terminal_sha256':sha256_file(producer/'provenance/terminal_v1.json'),'wall_seconds':time.monotonic()-began,
      'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,'new_final_identity':False})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)});raise
