"""True source-by-lineage fine discovery and reviewed-route joint discovery.

Shared features and engine stay fixed. Unknown pools remain source-local;
small branches retain coarse memberships without pretending a fine fit exists.
"""
import argparse
import collections
import datetime
import gzip
import json
from pathlib import Path
import resource
import subprocess
import sys
import time

import numpy as np

from annotation_representation_pilot_v1 import pca_fit,harmonize,save
from compare_annotation_granularity_v1 import graph_grid
from annotation_granularity_metrics_v1 import granularity,unit_out_state,gamma_group_retention
from paired_granularity_review_v1 import paired_scores
from full_context_atlas_v1 import FullContextAtlas,sha256_file
from load_annotation_pilot_native_v1 import load_pilot
from native_lineage_routing_v1 import FAMILIES,UNRESOLVED


def fit_branches(x,source,route,positive,seed,parameters,resolutions,checkpoint=None,fit=None):
    source,route,positive=np.asarray(source),np.asarray(route),np.asarray(positive,bool)
    if x.shape[0]!=len(source) or len(route)!=len(source) or len(positive)!=len(source):raise ValueError('Branch input axes differ')
    result={m:{str(r):np.full(len(source),-1,np.int32) for r in resolutions} for m in ['A','C']}
    fine={m:np.zeros(len(source),bool) for m in ['A','C']};records=[]
    if fit is None:
        def fit(member,correct):
            embedding,_=pca_fit(x[member],seed,parameters['components'],parameters['pca_iterations'])
            if correct:embedding,_=harmonize(embedding,source[member],seed,parameters)
            return graph_grid(embedding,seed,parameters,resolutions)
    offsets={m:{str(r):0 for r in resolutions} for m in ['A','C']}
    def put(method,member,labels):
        for resolution,lab in labels.items():
            _,codes=np.unique(lab,return_inverse=True)
            result[method][resolution][member]=codes+offsets[method][resolution]
            offsets[method][resolution]+=len(np.unique(codes))
    for ds in np.unique(source):
        for parent in np.unique(route[source==ds]):
            member=np.flatnonzero((source==ds)&(route==parent)&positive)
            if not len(member):continue
            eligible=len(member)>=parameters['minimum_fit_group']
            labels=fit(member,False) if eligible else {str(r):np.zeros(len(member),np.int32) for r in resolutions}
            put('A',member,labels);fine['A'][member]=eligible
            record={'method':'A','source':str(ds),'route':str(parent),'n_cells':len(member),
                    'status':'source_lineage_fine_graph' if eligible else 'small_branch_coarse_only'}
            records.append(record)
            if checkpoint:checkpoint('source_lineage_branch_complete',**record)
    for parent in np.unique(route):
        if parent==UNRESOLVED:
            for ds in np.unique(source):
                member=np.flatnonzero((route==parent)&(source==ds)&positive)
                if not len(member):continue
                put('C',member,{str(r):result['A'][str(r)][member] for r in resolutions})
                fine['C'][member]=fine['A'][member]
                records.append({'method':'C','source':str(ds),'route':str(parent),'n_cells':len(member),
                                'status':'unknown_pool_reuses_source_local_A_not_forced_joint_alignment'})
            continue
        member=np.flatnonzero((route==parent)&positive)
        if not len(member):continue
        eligible=len(member)>=parameters['minimum_fit_group']
        labels=fit(member,True) if eligible else {str(r):np.zeros(len(member),np.int32) for r in resolutions}
        put('C',member,labels);fine['C'][member]=eligible
        record={'method':'C','route':str(parent),'n_cells':len(member),'n_sources':len(np.unique(source[member])),
                'status':'joint_lineage_fine_graph' if eligible else 'small_branch_coarse_only'}
        records.append(record)
        if checkpoint:checkpoint('joint_lineage_branch_complete',**record)
    for method,variants in result.items():
        for resolution,labels in variants.items():
            if not np.array_equal(labels>=0,positive):raise ValueError('Positive-RNA cells lost or zero-RNA cells modeled')
            for code in np.unique(labels[positive]):
                member=labels==code
                if len(np.unique(route[member]))!=1:raise ValueError('Candidate group crosses coarse routing branches')
                if (method=='A' or route[member][0]==UNRESOLVED) and len(np.unique(source[member]))!=1:
                    raise ValueError('Source-local group crosses sources')
    return result,fine,records


def main(root):
    started=time.monotonic();cfg=json.loads((root/'config/source_lineage_spec_v1.json').read_text())
    for binding in cfg['bindings'].values():
        if sha256_file(binding['path'])!=binding['sha256']:raise ValueError('Source-lineage frozen input changed')
    def checkpoint(stage,**extra):
        save(root/'checkpoints'/(str(time.time_ns())+'.json'),{'observed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'stage':stage,'wall_seconds':time.monotonic()-started,**extra})
    test=subprocess.run([sys.executable,'-B',str(root/'tests/test_source_lineage_discovery_v1.py')],capture_output=True,text=True,timeout=180)
    (root/'provenance/test_stdout.txt').write_text(test.stdout);(root/'provenance/test_stderr.txt').write_text(test.stderr)
    if test.returncode!=0 or 'Ran 4 tests' not in test.stderr or 'skipped' in test.stderr.lower():raise ValueError('Source-lineage branch counterexamples failed')
    save(root/'provenance/numerical_tests_v1.json',{'actual_exit_code':test.returncode,'n_tests':4,'no_skips':True})
    routing=Path(cfg['routing']);previous=Path(cfg['previous_granularity']);discovery=Path(cfg['discovery']);pilot=Path(cfg['pilot'])
    if json.loads((routing/'provenance/supervisor_terminal_v1.json').read_text())['actual_child_exit_code']!=0:raise ValueError('Routing stage incomplete')
    with np.load(routing/'outputs/provisional_native_routes_v1.npz',allow_pickle=False) as f:
        source,oldroute,positive,routecodes,reason=f['source'],f['inherited_route'],f['positive_native_RNA'],f['route_code'],f['reason_code']
    route=np.asarray([FAMILIES[k] if k>=0 else UNRESOLVED for k in routecodes])
    data=load_pilot(pilot,cfg['manifest'],checkpoint=checkpoint)
    if not np.array_equal(source,data['source']) or not np.array_equal(oldroute,data['route']) or not np.array_equal(positive,data['native_total']>0):
        raise ValueError('Pilot versus routing axes differ')
    selected=json.loads((discovery/'outputs/discovery_feature_selection_v1.json').read_text());features=selected['selected_feature_indices']
    if [data['genes'][i] for i in features]!=selected['selected_feature_names']:raise ValueError('Shared feature axis differs')
    x=data.pop('x')[:,features];parameters=json.loads((discovery/'config/model_spec_v1.json').read_text())
    variants,fine,records=fit_branches(x,source,route,positive,cfg['seed'],parameters,cfg['resolutions'],checkpoint=checkpoint)
    del x
    with np.load(previous/'outputs/granularity_partitions_v1.npz',allow_pickle=False) as f:
        variants['B']={str(r):f['B__'+str(r)] for r in cfg['resolutions']}
    if any(not np.array_equal(v>=0,positive) for v in variants['B'].values()):raise ValueError('Frozen B comparison mask differs')
    side=root/'outputs/source_lineage_partitions_v1.npz'
    arrays={m+'__'+r:a for m,v in variants.items() for r,a in v.items()}
    with side.open('xb') as f:np.savez_compressed(f,**arrays)
    with np.load(side,allow_pickle=False) as f:
        for key,v in arrays.items():np.testing.assert_array_equal(f[key],v)
    # Compatible first-seed entry point for a later native fine-group review.
    with (root/'outputs'/('candidate_partitions_'+str(cfg['seed'])+'.npz')).open('xb') as f:
        np.savez_compressed(f,**{m:v['1.0'] for m,v in variants.items()})
    with (root/'outputs/derived_axis_and_model_mask.npz').open('xb') as f:
        np.savez_compressed(f,source=source,inherited_route=oldroute,provisional_native_route=route,modeled=positive,
                            A_fine_discovery_eligible=fine['A'],C_fine_discovery_eligible=fine['C'],routing_reason_code=reason)
    save(root/'outputs/branch_fitting_summary_v1.json',records)
    checkpoint('all_branch_fits_and_partition_readback_complete',n_branch_records=len(records))
    reader=FullContextAtlas(cfg['manifest']);units=[]
    for entry in json.loads((pilot/'outputs/pilot_input_inventory_v1.json').read_text()):
        positions=np.load(pilot/entry['selector_file'],allow_pickle=False)
        units.extend(reader.frozen_expression_units(entry['dataset'])[positions].astype(str))
    units=np.asarray(units);del reader
    profiles=[];pairs=[];coverage=[]
    for context,parents in [('inherited',oldroute),('provisional_native',route)]:
        for ds in np.unique(source):
            for parent in np.unique(parents[source==ds]):
                member=(source==ds)&(parents==parent)&positive
                if member.sum()<cfg['minimum_stratum_cells']:continue
                u,depth=units[member],np.log1p(data['native_total'][member]);cache={};size={}
                for method,parts in variants.items():
                    for resolution,labels in parts.items():
                        key=method+'__'+resolution;lab=labels[member];size[key]=granularity(lab);cache[key]={}
                        row={'context':context,'source':str(ds),'parent_not_truth':str(parent),'method':method,'resolution':float(resolution),
                             'granularity':size[key],'n_expression_units_not_patients':len(np.unique(u)),
                             'fine_discovery_fraction':float(fine[method][member].mean()) if method in fine else 1.,
                             'gamma_RNA_candidate':gamma_group_retention(lab,data['gamma_RNA_candidate'][member]),'state_prediction':{}}
                        for state,v in data['state'].items():
                            result=unit_out_state(lab,u,depth,v[member],minimum_other_cells=cfg['minimum_other_unit_group_cells'],return_predictions=True)
                            cache[key][state]=result;row['state_prediction'][state]={k:v for k,v in result.items() if k not in ['prediction','baseline']}
                        profiles.append(row)
                for a,b in [('A','B'),('A','C'),('B','C')]:
                    nmatch=0
                    for ra in cfg['resolutions']:
                        for rb in cfg['resolutions']:
                            ka,kb=a+'__'+str(ra),b+'__'+str(rb);pa,pb=size[ka],size[kb]
                            ratio=max(pa['effective_group_count'],pb['effective_group_count'])/min(pa['effective_group_count'],pb['effective_group_count'])
                            if ratio>cfg['maximum_effective_group_ratio']:continue
                            nmatch+=1
                            for state,v in data['state'].items():
                                pairs.append({'context':context,'source':str(ds),'parent_not_truth':str(parent),
                                    'first_method':a,'second_method':b,'first_resolution':ra,'second_resolution':rb,'state':state,
                                    'effective_group_ratio':ratio,**paired_scores(v[member],u,cache[ka][state],cache[kb][state])})
                    coverage.append({'context':context,'source':str(ds),'parent_not_truth':str(parent),'first_method':a,'second_method':b,
                                     'n_matched_granularity_pairs':nmatch,'status':'overlap' if nmatch else 'no_overlap_in_frozen_grid'})
                checkpoint('source_parent_diagnostics_complete',context=context,source=str(ds),parent=str(parent))
    output=root/'outputs/source_lineage_comparison_profiles_v1.json.gz'
    with gzip.open(output,'xt') as f:json.dump({'profiles':profiles,'paired_same_query_results':pairs,'granularity_coverage':coverage,
        'state_gene_sets':data['state_genes'],'routes_are_proposals_not_truth':True,'clinical_outcomes_used':False,
        'discovery_includes_all_pilot_sources':True,'references_may_encode_RNA_readouts':True,'new_method_winner_selected':False,
        'A_is_source_by_lineage_fine_discovery':True,'A_cross_source_identity_coordination_pending':True,
        'C_known_routes_joint_unknown_pool_source_local':True,'expression_units_are_not_patients':True},f,allow_nan=False)
    size=sum(p.stat().st_size for p in root.rglob('*') if p.is_file())
    if size>cfg['output_cap_bytes']:raise ValueError('Source-lineage output cap exceeded')
    save(root/'provenance/terminal_v1.json',{'status':'SOURCE_LINEAGE_FINE_DISCOVERY_AND_PAIRED_DIAGNOSTICS_COMPLETE_NATIVE_FINE_REVIEW_PENDING',
        'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_tests_actual_no_skips':4,
        'n_profile_records':len(profiles),'n_pair_state_records':len(pairs),'partitions_sha256':sha256_file(side),
        'output_sha256':sha256_file(output),'wall_seconds':time.monotonic()-started,
        'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,'written_bytes_before_terminal':size,
        'native_counts_copied':False,'new_final_identity':False,'method_winner_selected':False})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)})
        raise
