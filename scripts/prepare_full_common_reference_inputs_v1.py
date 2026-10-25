"""Freeze common-coordinate reference rows without replacing local identities.

Includes an all-fit-row global overview pool and known-lineage pools. Rare
local/guarded/gamma strata are reserved; clinical outcomes are never read.
"""
import argparse
import datetime
import gc
import json
from pathlib import Path
import resource
import subprocess
import sys
import time

import numpy as np
from annotation_pilot_sampling_v1 import stable_seed
from annotation_representation_pilot_v1 import save
from audit_full_source_local_graphs_v1 import index_at
from common_reference_sampling_v1 import reserve_then_balance,stratum_reserve,select_reference
from extract_full_source_local_marker_evidence_v1 import storage_gate
from full_context_atlas_v1 import sha256_file
from native_archive_collection_v1 import NativeArchiveCollection,text_axis_sha256
from native_lineage_routing_v1 import FAMILIES,local_family
from summarize_full_native_annotation_evidence_v1 import read_source


def source_pool(reader,root,record):
    meta,side=read_source(root,record);ds=record['source'];rows=reader.rows(ds)
    np.testing.assert_array_equal(rows,side['native_row'])
    if meta['native_ids_sha256']!=text_axis_sha256(reader.native_ids(ds)):raise ValueError('Reference pool native identity changed')
    eligible=side['fit_pool_candidate']&side['positive_native_RNA'];major=np.asarray(meta['dictionaries']['consensus_major'])[side['consensus_major_code']]
    fine=np.asarray(meta['dictionaries']['consensus_celltype'])[side['consensus_celltype_code']]
    route=np.full(len(rows),-1,np.int16)
    for value in np.unique(major):
        family=local_family(str(value))
        if family is not None:route[major==value]=FAMILIES.index(family)
    gamma=(side['own_fine_bitmask']&(1<<meta['fine_candidate_order'].index('T_gamma_delta_RNA_candidate')))!=0
    strata=np.char.add(np.char.add(np.char.add(side['parent_source_local_group_code'].astype(str),'::'),fine.astype(str)),np.where(gamma,'::own_gamma_RNA','::no_own_gamma_RNA'))
    units=reader.frozen_expression_units(ds).astype(str)
    return meta,side,eligible,route,strata,units


def save_arrays(path,**arrays):
    with path.open('xb') as f:np.savez_compressed(f,**arrays)
    with np.load(path,allow_pickle=False) as f:
        if set(f.files)!=set(arrays):raise ValueError('Reference input serialized keys differ')
        for k,v in arrays.items():np.testing.assert_array_equal(v,f[k])


def main(root):
    began=time.monotonic();cfg=json.loads((root/'config/common_reference_input_spec_v1.json').read_text())
    for binding in cfg['bindings'].values():
        if sha256_file(binding['path'])!=binding['sha256']:raise ValueError('Common reference input binding changed')
    test=subprocess.run([sys.executable,'-B',str(root/'tests/test_common_reference_sampling_v1.py')],capture_output=True,text=True,timeout=120)
    (root/'provenance/test_stderr.txt').write_text(test.stderr)
    if test.returncode or 'Ran 3 tests' not in test.stderr or 'skipped' in test.stderr.lower():raise ValueError('Reference sampler counterexamples failed')
    save(root/'provenance/numerical_tests_v1.json',{'actual_exit_code':0,'n_tests':3,'no_skips':True})
    candidate=Path(cfg['candidate']);audit=Path(cfg['candidate_audit'])
    for stage in [candidate,audit]:
        sup=json.loads((stage/'provenance/supervisor_terminal_v1.json').read_text())
        if sup['actual_child_exit_code'] or sup['wall_timeout']:raise ValueError('Completed candidate and independent audit required')
    at=json.loads((audit/'provenance/terminal_v1.json').read_text())
    if at['producer_terminal_sha256']!=sha256_file(candidate/'provenance/terminal_v1.json'):raise ValueError('Candidate audit belongs to different producer')
    index=index_at(candidate,'full_native_cell_identity_candidate_index_v1.json');inputs=Path(cfg['native_local_inputs']);fi=index_at(inputs,'full_native_local_discovery_input_index_v1.json')
    reader=NativeArchiveCollection(cfg['archive_manifest']);sources=sorted(index)
    if set(sources)!=set(reader.archive_partitions) or set(fi)!=set(sources):raise ValueError('Common reference source scopes differ')
    fp=Path(cfg['feature_coverage']);ft=json.loads((fp/'provenance/terminal_v1.json').read_text());fpath=fp/'outputs/full_archive_feature_coverage_v1.json'
    if sha256_file(fpath)!=ft['output_sha256']:raise ValueError('Full feature coverage changed')
    feature_spec=json.loads(fpath.read_text());common=set(feature_spec['common_nonheldout_feature_candidates']);source_feature_files={}
    for ds in sources:
        mp=inputs/'outputs'/fi[ds]['metadata']
        if sha256_file(mp)!=fi[ds]['metadata_sha256']:raise ValueError('Native feature input metadata changed')
        fm=json.loads(mp.read_text());path=inputs/'outputs'/fm['feature_candidates']
        if sha256_file(path)!=fm['feature_candidates_sha256']:raise ValueError('Native feature candidate list changed')
        f=json.loads(path.read_text());common&=set(f['genes']);source_feature_files[ds]=(path,fm['feature_candidates_sha256'])
    genes=sorted(common)
    if len(genes)<cfg['minimum_common_candidates']:raise ValueError('Too few shared measured alias-excluded candidates')
    features=root/'outputs/common_measured_model_feature_candidates_v1.json'
    save(features,{'genes':genes,'n_sources':len(sources),'all_sources_measured_unambiguous':True,'heldout_RNA_and_physical_aliases_excluded':True,
      'feature_selection_fitted':False,'original_native_all_gene_archive_unchanged':True,'feature_coverage_source_sha256':ft['output_sha256']})
    modality=Path(cfg['modality_overlay']);mt=json.loads((modality/'provenance/terminal_v1.json').read_text());mfile=modality/'outputs/source_modality_overlay_v1.json'
    if sha256_file(mfile)!=mt['output_sha256']:raise ValueError('Source modality overlay changed')
    modes={r['source']:r for r in json.loads(mfile.read_text())['sources']}
    if set(modes)!=set(sources):raise ValueError('Modality overlay sources differ')
    pools=['global_overview',*FAMILIES];capacities={name:[] for name in pools};reserves={name:[] for name in pools}
    for ds in sources:
        meta,side,eligible,route,strata,units=source_pool(reader,candidate,index[ds])
        for name in pools:
            member=eligible if name=='global_overview' else eligible&(route==FAMILIES.index(name))
            capacities[name].append(int(member.sum()));reserves[name].append(stratum_reserve(strata,member,cfg['rare_reserve']))
        save(root/'checkpoints'/('inventory_'+ds+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':ds,'stage':'reference_pool_inventory'})
        reader._fields.clear();del side,eligible,route,strata,units;gc.collect()
    allocations={name:reserve_then_balance(capacities[name],reserves[name],cfg['global_reference_cap'] if name=='global_overview' else cfg['lineage_reference_cap'],stable_seed(cfg['seed'],'all_sources',name)).tolist() for name in pools}
    save(root/'outputs/reference_pool_allocation_v1.json',{'sources':sources,'pools':pools,'capacities':capacities,'minimum_reserves':reserves,'allocated':allocations,
      'expression_units_not_assumed_independent_patients':True,'coarse_candidate_labels_are_sampling_strata_not_validation_truth':True,
      'clinical_or_CM_selection':False,'all_source_transductive_reference_exposure':True})
    output_index=[]
    for si,ds in enumerate(sources):
        storage_gate(root,cfg);meta,side,eligible,route,strata,units=source_pool(reader,candidate,index[ds]);n=len(eligible)
        selected_global=np.zeros(n,bool);selected_lineage=np.zeros(n,bool);reports={}
        for name in pools:
            member=eligible if name=='global_overview' else eligible&(route==FAMILIES.index(name))
            selection=select_reference(units,strata,member,allocations[name][si],stable_seed(cfg['seed'],ds,name),cfg['rare_reserve'])
            if name=='global_overview':selected_global[selection]=True
            else:selected_lineage[selection]=True
            reports[name]={'n_fit_pool_rows':int(member.sum()),'n_selected_reference_rows':len(selection),
              'n_local_guarded_gamma_strata':int(len(np.unique(strata[member]))),
              'n_selected_local_guarded_gamma_strata':int(len(np.unique(strata[selection]))),
              'minimum_reserved_rows':reserves[name][si],'n_selected_expression_units':int(len(np.unique(units[selection])))}
        if np.any(selected_lineage&(route<0)) or np.any((selected_global|selected_lineage)&~eligible):raise ValueError('Reference exposure outside eligible route')
        source_feature_path,source_feature_sha=source_feature_files[ds];f=json.loads(source_feature_path.read_text())
        if sha256_file(source_feature_path)!=source_feature_sha:raise ValueError('Source feature file changed during sampling')
        column_by_gene=dict(zip(f['genes'],f['native_columns']));columns=np.asarray([column_by_gene[g] for g in genes],np.int64);lookup=reader.feature_lookup(ds)
        if len(set(columns))!=len(columns) or any(lookup.get(g,set())!={int(c)} for g,c in zip(genes,columns)):raise ValueError('Common native feature axis differs')
        if set(columns)&set(f['heldout_physical_alias_columns_excluded']):raise ValueError('Heldout physical alias entered common feature axis')
        colfile=root/'outputs'/(ds+'_common_native_feature_columns_v1.npz');save_arrays(colfile,native_columns=columns)
        arrayfile=root/'outputs'/(ds+'_common_reference_rows_v1.npz');save_arrays(arrayfile,native_row=side['native_row'],fit_pool_candidate=eligible,
          guarded_route_code=route,parent_source_local_group_code=side['parent_source_local_group_code'],global_reference_fit=selected_global,lineage_reference_fit=selected_lineage)
        info={'source':ds,'n_archive_rows':n,'native_ids_sha256':meta['native_ids_sha256'],'candidate_record':index[ds],
          'source_modality_context':modes[ds],'native_feature_columns':colfile.name,'native_feature_columns_sha256':sha256_file(colfile),
          'common_feature_list':features.name,'common_feature_list_sha256':sha256_file(features),'guarded_route_dictionary':FAMILIES,'unresolved_route_code':-1,
          'reference_pools':reports,'reference_fit_exposure_is_not_cell_identity_accuracy':True,'local_fine_graphs_and_identities_not_replaced':True,
          'modality_and_anatomy_transfer_not_accepted_by_sampling':True,'source_alias_and_patient_qualification_not_promoted':True,
          'global_and_lineage_fit_exposures_can_overlap_not_independent_models':True,'native_counts_read':False,'new_final_identity':False}
        mp=arrayfile.with_suffix('.json');save(mp,info)
        output_index.append({'source':ds,'n_archive_rows':n,'arrays':arrayfile.name,'arrays_sha256':sha256_file(arrayfile),'metadata':mp.name,'metadata_sha256':sha256_file(mp)})
        save(root/'checkpoints'/('selected_'+ds+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':ds,'stage':'common_reference_rows_frozen','completed_sources':len(output_index)})
        reader._fields.clear();del side,eligible,route,strata,units;gc.collect()
    output=root/'outputs/full_common_reference_input_index_v1.json';save(output,output_index)
    save(root/'provenance/terminal_v1.json',{'status':'FULL_COMMON_REFERENCE_ROWS_AND_MEASURED_FEATURE_INPUTS_PREPARED_NO_COORDINATES_FITTED',
      'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_sources':len(output_index),'n_archive_rows':sum(r['n_archive_rows'] for r in output_index),
      'n_common_alias_excluded_feature_candidates':len(genes),'n_global_reference_rows':sum(allocations['global_overview']),
      'n_lineage_reference_rows':sum(sum(allocations[name]) for name in FAMILIES),'index_sha256':sha256_file(output),
      'candidate_auditor_terminal_sha256':sha256_file(audit/'provenance/terminal_v1.json'),'wall_seconds':time.monotonic()-began,
      'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,'native_counts_read':False,'new_final_identity':False})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)});raise
