"""Audit complete structure coverage and replay positive/negative/unestimable cases.

Selects numerical-audit cases by declared metric status, not biological ranking.
Local labels and heldout RNA programs are not independent identity truth.
"""
import argparse,datetime,json,time,resource
from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist
from annotation_representation_pilot_v1 import save
from annotation_graph_diagnostics_v1 import direct_unit_out_graph_state
from common_reference_structure_metrics_v1 import exact_neighbor_graph
from extract_full_source_local_marker_evidence_v1 import storage_gate
from full_context_atlas_v1 import sha256_file
from native_archive_collection_v1 import row_axis_sha256
from native_lineage_routing_v1 import FAMILIES


def scope_key(r):return (r['route'],r['source'],r['conditioning_guarded_type_not_truth'])
def profile_key(r):return (*scope_key(r),r['panel_scope'],r['variant'],r['state'])


def main(root):
    start=time.monotonic();cfg=json.loads((root/'config/structure_audit_spec_v1.json').read_text());producer=Path(cfg['producer']);sup=json.loads((producer/'provenance/supervisor_terminal_v1.json').read_text());terminal=json.loads((producer/'provenance/terminal_v1.json').read_text())
    if sup['actual_child_exit_code'] or sup['wall_timeout']:raise ValueError('Successful complete structure producer required')
    for b in cfg['bindings'].values():
        if sha256_file(b['path'])!=b['sha256']:raise ValueError('Structure audit frozen input changed')
    file=producer/'outputs/common_reference_native_structure_diagnostics_v1.json'
    if sha256_file(file)!=terminal['output_sha256']:raise ValueError('Structure output hash differs')
    v=json.loads(file.read_text());spec=json.loads((producer/'config/common_structure_spec_v1.json').read_text());models=Path(spec['producer']);mspec=json.loads((models/'config/common_reference_model_spec_v1.json').read_text())
    with np.load(models/'outputs/reference_union_axes_and_native_readouts_v1.npz',allow_pickle=False) as f:axes={k:f[k] for k in f.files}
    am=json.loads((models/'outputs/reference_union_axis_metadata_v1.json').read_text());sources=am['sources'];sampling=json.loads((producer/'outputs/private_fixed_diagnostic_sampling_v1.json').read_text());sample_map={scope_key(r):r for r in sampling};coverage={scope_key(r):r for r in v['coverage'] if 'source' in r};profiles={profile_key(r):r for r in v['state_profiles']}
    if len(profiles)!=len(v['state_profiles']) or len(sample_map)!=len(sampling) or len(coverage)!=sum('source' in r for r in v['coverage']):raise ValueError('Duplicate diagnostic keys')
    if len(profiles)!=terminal['n_profiles'] or len(v['paired_state_profiles'])!=terminal['n_paired_profiles'] or len(v['coverage'])!=terminal['n_coverage_records']:raise ValueError('Diagnostic record counts differ')
    sampled_total=0
    for k,r in sample_map.items():
        c=coverage[k];q=np.asarray(r['union_reference_rows'],np.int64)
        if len(q)!=len(set(q)) or len(q)<spec['minimum_stratum_cells'] or len(q)>spec['maximum_stratum_cells'] or np.any(q<0) or np.any(q>=len(axes['native_row'])):raise ValueError('Invalid diagnostic sampling rows')
        if row_axis_sha256(q)!=r['row_axis_sha256'] or len(q)!=c['n_sampled_reference_rows']:raise ValueError('Fixed sample hash/coverage differs')
        if not np.all(axes['source_code'][q]==sources.index(r['source'])) or len(np.unique(axes['expression_unit_code'][q]))!=r['n_expression_units_not_patients']:raise ValueError('Diagnostic source/unit axes differ')
        if r['route']=='global_overview':eligible=axes['global_reference_fit'][q]
        else:eligible=axes['lineage_reference_fit'][q]&(axes['route_code'][q]==FAMILIES.index(r['route']))
        if not eligible.all() or c['status']!='diagnostic_computed_not_biological_validation':raise ValueError('Diagnostic sampled outside model pool')
        sampled_total+=len(q)
    for k,r in profiles.items():
        q=sample_map[k[:3]]['union_reference_rows']
        if r['n_input_cells']!=len(q) or not 0<=r['n_evaluable']<=len(q) or r['n_sampled_reference_rows']!=len(q):raise ValueError('State diagnostic denominator differs')
        if r['state'] not in am['state_order'] or r['variant'].split('__')[0]!=r['route']:raise ValueError('State/model scope differs')
        if r['panel_scope']=='complete_frozen_panel_primary' and not coverage[k[:3]]['complete_state_panels'][r['state']] and r['n_evaluable']:raise ValueError('Incomplete panel treated as full measured state')
        if r['status']=='estimated' and (r['n_evaluable']<100 or r['depth_baseline_SSE']<1e-10):raise ValueError('Estimated state lacks denominator or variance')
    for r in v['paired_state_profiles']:
        a=profiles[(*scope_key(r),r['panel_scope'],r['first_variant'],r['state'])];b=profiles[(*scope_key(r),r['panel_scope'],r['second_variant'],r['state'])]
        if r['n_common_query_cells']>min(a['n_evaluable'],b['n_evaluable']):raise ValueError('Paired query denominator exceeds either arm')
        if r['status']=='estimated':np.testing.assert_allclose(r['second_minus_first'],r['score_second']-r['score_first'],rtol=1e-10,atol=1e-10)
    cases=[]
    for kind,predicate in [('positive',lambda r:r['status']=='estimated' and r['RNA_prediction_improvement_over_depth']>0),('negative',lambda r:r['status']=='estimated' and r['RNA_prediction_improvement_over_depth']<0),('not_estimable',lambda r:r['status']!='estimated')]:
        matches=sorted([r for r in v['state_profiles'] if r['panel_scope']=='complete_frozen_panel_primary' and predicate(r)],key=profile_key)
        if not matches:cases.append({'case':kind,'status':'no_case_available_in_fixed_profiles'});continue
        r=matches[0];q=np.asarray(sample_map[scope_key(r)]['union_reference_rows'],np.int64);identifier='__'.join(r['variant'].split('__')[:2]);record=next(x for x in json.loads((models/'outputs/common_reference_candidate_model_index_v1.json').read_text()) if x['model_id']==identifier)
        with np.load(models/'outputs'/record['coordinates'],allow_pickle=False) as f:
            rows=f['reference_union_row'];positions=np.searchsorted(rows,q);np.testing.assert_array_equal(rows[positions],q);z=f['raw_PCA_fit_embedding' if r['variant'].endswith('__raw_PCA') else 'reference_embedding'][positions]
        graph=exact_neighbor_graph(z);d=cdist(z.astype(float),z.astype(float));np.fill_diagonal(d,np.inf);k=min(30,len(q)-1);cutoff=np.partition(d,k-1,axis=1)[:,k-1]
        for i in range(len(q)):
            near=graph.indices[graph.indptr[i]:graph.indptr[i+1]]
            if i in near or len(near)!=k or np.any(d[i,near]>cutoff[i]+2e-4):raise ValueError('Nearest-neighbor graph violates independent distance cutoff')
        units=axes['expression_unit_code'][q];depth=np.log1p(axes['native_total_RNA'][q]);values=axes['state_readouts'][q,am['state_order'].index(r['state'])].astype(float)
        if not coverage[scope_key(r)]['complete_state_panels'][r['state']]:values[:]=np.nan
        pred,base=direct_unit_out_graph_state(graph,units,depth,values,spec['minimum_other_unit_neighbors'],spec['minimum_train_cells']);ok=np.isfinite(pred)
        if int(ok.sum())!=r['n_evaluable']:raise ValueError('Independent query availability differs')
        error=None
        if r['status']=='estimated':
            baseline=float(np.sum((values[ok]-base[ok])**2));sse=float(np.sum((values[ok]-pred[ok])**2));score=1-sse/baseline;error=abs(score-r['RNA_prediction_improvement_over_depth']);np.testing.assert_allclose(baseline,r['depth_baseline_SSE'],rtol=1e-8,atol=1e-8);np.testing.assert_allclose(sse,r['neighbor_plus_depth_SSE'],rtol=1e-8,atol=1e-8);np.testing.assert_allclose(score,r['RNA_prediction_improvement_over_depth'],rtol=1e-8,atol=1e-8)
        cases.append({'case':kind,'status':'independently_replayed','profile_key':list(profile_key(r)),'n_evaluable':int(ok.sum()),'score_absolute_difference':error,'independent_CSR_neighbor_distance_check':True})
    out=root/'outputs/common_reference_structure_audit_v1.json';save(out,{'all_coverage_and_profile_denominators_checked':True,'n_fixed_source_type_strata':len(sample_map),'n_sampled_exposures_not_independent_cells':sampled_total,'numeric_audit_cases':cases,'scientific_model_selection_or_identity_truth_claim':False})
    save(root/'provenance/terminal_v1.json',{'status':'COMMON_REFERENCE_STRUCTURE_DENOMINATORS_AND_DECLARED_NUMERIC_CASES_AUDITED','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'output_sha256':sha256_file(out),'producer_terminal_sha256':sha256_file(producer/'provenance/terminal_v1.json'),'n_state_profiles_checked':len(profiles),'n_paired_profiles_checked':len(v['paired_state_profiles']),'wall_seconds':time.monotonic()-start,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)});raise
