"""Compare frozen granularity matches on identical query cells.

No representation is refitted and no candidate is selected by RNA scores.
Per-cell readouts and predictions are RAM-only. Expression units are not
assumed to be patients; inherited routes are not annotation truth.
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

from annotation_granularity_metrics_v1 import granularity, unit_out_state, gamma_group_retention
from annotation_representation_pilot_v1 import ROUTES, save
from full_context_atlas_v1 import FullContextAtlas, sha256_file
from native_panel_aggregate_v1 import gather_csr_rows


def direct_unit_out(labels, units, depth, values, minimum_other_cells=20):
    """Independent centered OLS and direct group means, without group totals."""
    labels, units = np.asarray(labels), np.asarray(units)
    depth, values = np.asarray(depth, float), np.asarray(values, float)
    good = (labels >= 0) & np.isfinite(depth) & np.isfinite(values)
    prediction, baseline = np.full(len(labels), np.nan), np.full(len(labels), np.nan)
    for unit in np.unique(units[good]):
        train = good & (units != unit)
        query = good & (units == unit)
        if train.sum() < 30:continue
        x, y = depth[train], values[train]
        xc = x - x.mean();denom = float(xc @ xc)
        slope = float(xc @ (y - y.mean())) / denom if denom > 0 else 0.
        intercept = y.mean() - slope * x.mean()
        for label in np.unique(labels[query]):
            group = train & (labels == label)
            if group.sum() < minimum_other_cells:continue
            target = query & (labels == label)
            base = intercept + slope * depth[target]
            residual = values[group] - (intercept + slope * depth[group])
            baseline[target] = base;prediction[target] = base + residual.mean()
    return prediction, baseline


def paired_scores(values, units, a, b):
    values, units = np.asarray(values, float), np.asarray(units)
    common = (np.isfinite(values) & np.isfinite(a['prediction']) & np.isfinite(b['prediction'])
              & np.isfinite(a['baseline']) & np.isfinite(b['baseline']))
    out = {'status': 'insufficient_common_query_cells', 'n_input_cells': len(values),
           'n_common_query_cells': int(common.sum()), 'common_query_fraction': float(common.mean()),
           'n_common_expression_units_not_patients': len(np.unique(units[common])),
           'high_common_coverage_at_least90pct': bool(common.mean() >= .9)}
    if common.sum() < 100:return out
    np.testing.assert_allclose(a['baseline'][common], b['baseline'][common], atol=1e-9, rtol=1e-9,
                               err_msg='Compared partitions have different depth baseline training domains')
    baseline_sse = float(np.sum((values[common] - a['baseline'][common]) ** 2))
    out['depth_baseline_SSE'] = baseline_sse
    if baseline_sse < 1e-10:out['status'] = 'no_variable_depth_residual';return out
    score = []
    for item in [a, b]:
        score.append(1 - float(np.sum((values[common] - item['prediction'][common]) ** 2)) / baseline_sse)
    out.update(status='estimated', score_first=score[0], score_second=score[1],
               second_minus_first=score[1] - score[0], identical_query_cells_and_depth_baseline=True)
    return out


def compare_json_numeric(a, b, atol=1e-8):
    if isinstance(a, dict):
        if a.keys() != b.keys():raise ValueError('Diagnostic fields changed')
        for k in a:compare_json_numeric(a[k], b[k], atol)
    elif isinstance(a, (float, int)) and not isinstance(a, bool):
        np.testing.assert_allclose(a, b, atol=atol, rtol=1e-7)
    elif a != b:raise ValueError('Diagnostic categorical value changed')


def main(root):
    started = time.monotonic();cfg = json.loads((root/'config/paired_review_spec_v1.json').read_text())
    for binding in cfg['bindings'].values():
        if sha256_file(binding['path']) != binding['sha256']:raise ValueError('Paired-review frozen input changed')
    test = subprocess.run([sys.executable, '-B', str(root/'tests/test_paired_granularity_review_v1.py')],
                          capture_output=True, text=True, timeout=120)
    (root/'provenance/test_stdout.txt').write_text(test.stdout);(root/'provenance/test_stderr.txt').write_text(test.stderr)
    if test.returncode != 0 or 'Ran 4 tests' not in test.stderr or 'skipped' in test.stderr.lower():
        raise ValueError('Paired scoring counterexamples failed')
    gr = Path(cfg['granularity']);pilot = Path(cfg['pilot']);discovery = Path(cfg['discovery'])
    original_cfg = json.loads((gr/'config/granularity_spec_v1.json').read_text())
    for binding in original_cfg['bindings'].values():
        if sha256_file(binding['path']) != binding['sha256']:raise ValueError('Original frozen input changed')
    inventory = json.loads((pilot/'outputs/pilot_input_inventory_v1.json').read_text())
    audit_source = min(inventory, key=lambda r:(r['n_selected_cells'],r['dataset']))['dataset']
    plan = json.loads((pilot/'config/model_features_and_heldouts_v1.json').read_text())
    panels = {name:[g for g in genes if all(c[g] for c in plan['source_measurement_coverage'].values())]
              for name,genes in plan['heldout_RNA_readouts']['state'].items()}
    genes = sorted({g for p in plan['heldout_RNA_readouts'].values() for panel in p.values() for g in panel})
    with gzip.open(gr/'outputs/source_parent_granularity_profiles_v1.json.gz','rt') as f:original = json.load(f)
    if original['state_gene_sets'] != panels:raise ValueError('Frozen readout gene sets changed')
    profile_map = {(r['source'],r['inherited_parent_not_truth'],r['method'],str(r['resolution'])):r for r in original['profiles']}
    with np.load(gr/'outputs/granularity_partitions_v1.npz',allow_pickle=False) as f:parts = {k:f[k] for k in f.files}
    with np.load(discovery/'outputs/derived_axis_and_model_mask.npz',allow_pickle=False) as f:
        axis_source,axis_route,axis_modeled = f['source'], f['inherited_route'], f['modeled']
    if any(len(v) != len(axis_source) for v in parts.values()):raise ValueError('Partition axis differs')
    reader = FullContextAtlas(cfg['manifest']);offset = 0;pairs = [];audited_profiles = 0
    independent_predictions = 0;independent_max_error = 0.;native_state_values_checked = 0;native_max_error = 0.
    for item in inventory:
        ds = item['dataset'];selector = pilot/item['selector_file']
        if sha256_file(selector) != item['selector_sha256']:raise ValueError('Selector changed')
        positions = np.load(selector,allow_pickle=False);n = len(positions);sl = slice(offset,offset+n)
        if not np.all(axis_source[sl] == ds):raise ValueError('Source axis differs')
        route = np.asarray([ROUTES.get(str(v),'unresolved_context') for v in reader.frozen_labels(ds)['major'][positions]])
        if not np.array_equal(route,axis_route[sl]):raise ValueError('Inherited route axis differs')
        units = reader.frozen_expression_units(ds)[positions].astype(str)
        lookup = reader.feature_lookup(ds);available = np.array([len(lookup.get(g,set()))==1 for g in genes])
        cols = np.array([next(iter(lookup[g])) for j,g in enumerate(genes) if available[j]],np.int64)
        store = reader.store(ds);key = reader.legacy_config['counts_key'] if reader.is_legacy(ds) else reader.primary.partitions[ds]['counts_key']
        totals,readouts = [],[]
        for lo in range(0,n,2048):
            raw = gather_csr_rows(store[key],reader.rows(ds)[positions[lo:lo+2048]])
            if not np.isfinite(raw.data).all() or (raw.data<0).any() or not np.equal(raw.data,np.floor(raw.data)).all():
                raise ValueError('Native count domain changed')
            total = np.asarray(raw.sum(axis=1,dtype=np.float64)).ravel()
            y = np.full((raw.shape[0],len(genes)),np.nan,np.float32);y[:,available] = raw[:,cols].toarray()
            scale = np.zeros(len(total),np.float32);scale[total>0] = 10000/total[total>0]
            y = np.log1p(y*scale[:,None]);totals.append(total);readouts.append(y)
            if ds == audit_source:
                # Different arithmetic: explicit scalar row totals and direct CSC columns.
                independent_total = np.array([raw.data[raw.indptr[i]:raw.indptr[i+1]].sum(dtype=float) for i in range(raw.shape[0])])
                np.testing.assert_array_equal(total,independent_total)
                csc = raw.tocsc()
                for j in np.flatnonzero(available):
                    column = next(iter(lookup[genes[j]]));counts = np.asarray(csc[:,column].toarray()).ravel()
                    expected = np.log1p(np.divide(counts*10000.,independent_total,out=np.zeros(len(total)),where=independent_total>0))
                    error = float(np.max(np.abs(expected-y[:,j])))
                    native_max_error = max(native_max_error,error);native_state_values_checked += len(total)
                    np.testing.assert_allclose(expected,y[:,j],atol=2e-6,rtol=2e-6)
        totals,readouts = np.concatenate(totals),np.vstack(readouts)
        if not np.array_equal(totals>0,axis_modeled[sl]):raise ValueError('Modeled domain differs')
        state = {name:readouts[:,[genes.index(g) for g in panel]].mean(axis=1) if len(panel)>=3 else np.full(n,np.nan)
                 for name,panel in panels.items()}
        gg = np.column_stack([readouts[:,genes.index(g)] for g in ['TRDC','TRGC1','TRGC2','CD3D']])
        gamma = np.where(np.isfinite(gg).all(axis=1),(gg[:,0]>0)&((gg[:,1]>0)|(gg[:,2]>0))&(gg[:,3]>0),np.nan)
        for parent in np.unique(route):
            member = (route==parent)&(totals>0)
            if member.sum()<original_cfg['minimum_stratum_cells']:continue
            u,depth = units[member],np.log1p(totals[member]);cache = {};profiles = {}
            for variant,labels in parts.items():
                method,resolution = variant.split('__');lab = labels[sl][member]
                old = profile_map[(ds,str(parent),method,resolution)]
                if old['n_expression_units_not_patients'] != len(np.unique(u)):raise ValueError('Unit denominator differs')
                g = granularity(lab);compare_json_numeric(g,old['granularity'])
                compare_json_numeric(gamma_group_retention(lab,gamma[member]),old['gamma_RNA_candidate'])
                profiles[variant] = g;cache[variant] = {}
                for name,value in state.items():
                    v = value[member]
                    result = unit_out_state(lab,u,depth,v,minimum_other_cells=original_cfg['minimum_other_unit_group_cells'],return_predictions=True)
                    summary = {k:v for k,v in result.items() if k not in ['prediction','baseline']}
                    compare_json_numeric(summary,old['state_prediction'][name])
                    if ds == audit_source:
                        pred,base = direct_unit_out(lab,u,depth,v,original_cfg['minimum_other_unit_group_cells'])
                        for x,y in [(pred,result['prediction']),(base,result['baseline'])]:
                            np.testing.assert_array_equal(np.isfinite(x),np.isfinite(y))
                            ok = np.isfinite(x)
                            if ok.any():independent_max_error=max(independent_max_error,float(np.max(np.abs(x[ok]-y[ok]))))
                            np.testing.assert_allclose(x,y,atol=1e-8,rtol=1e-8,equal_nan=True)
                        independent_predictions += 1
                    cache[variant][name] = result
                audited_profiles += 1
            for a,b in [('A','B'),('A','C'),('B','C')]:
                rebuilt = []
                for ra in original_cfg['resolutions']:
                    for rb in original_cfg['resolutions']:
                        va,vb = a+'__'+str(ra),b+'__'+str(rb);pa,pb = profiles[va],profiles[vb]
                        if pa['status']!='estimated' or pb['status']!='estimated':continue
                        ratio=max(pa['effective_group_count'],pb['effective_group_count'])/min(pa['effective_group_count'],pb['effective_group_count'])
                        if ratio>original_cfg['maximum_effective_group_ratio']:continue
                        rebuilt.append({'resolution_'+a:ra,'resolution_'+b:rb,'effective_group_count_ratio':ratio})
                        for name,value in state.items():
                            r = paired_scores(value[member],u,cache[va][name],cache[vb][name])
                            pairs.append({'source':ds,'inherited_parent_not_truth':str(parent),'first_method':a,'second_method':b,
                                          'first_resolution':ra,'second_resolution':rb,'state':name,'effective_group_ratio':ratio,**r})
                oldmatch = next(r for r in original['granularity_overlap'] if r['source']==ds and r['inherited_parent_not_truth']==parent and r['methods']==[a,b])
                compare_json_numeric({'matches':rebuilt},{'matches':oldmatch['overlapping_grid_pairs']})
        offset += n
        save(root/'checkpoints'/(ds+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
                                           'source':ds,'profiles_replayed':audited_profiles,'paired_records':len(pairs)})
    if offset!=len(axis_source) or audited_profiles!=len(original['profiles']):raise ValueError('Incomplete profile replay')
    groups = collections.defaultdict(list)
    for r in pairs:groups[tuple(r[k] for k in ['source','inherited_parent_not_truth','first_method','second_method','state'])].append(r)
    summaries = []
    for key,rows in sorted(groups.items()):
        estimated = [r for r in rows if r['status']=='estimated'];high = [r for r in estimated if r['high_common_coverage_at_least90pct']]
        out = dict(zip(['source','inherited_parent_not_truth','first_method','second_method','state'],key))
        out.update(n_frozen_matched_pairs=len(rows),n_estimable_pairs=len(estimated),n_high_common_coverage_pairs=len(high),
                   resolution_pairs_are_not_independent_replicates=True)
        for label,subset in [('all_estimable',estimated),('high_common_coverage',high)]:
            d = [r['second_minus_first'] for r in subset]
            out[label] = {'n_pairs':len(d),'median_difference':float(np.median(d)) if d else None,
                          'minimum_difference':float(min(d)) if d else None,'maximum_difference':float(max(d)) if d else None,
                          'fraction_second_higher':float(np.mean(np.asarray(d)>0)) if d else None,
                          'minimum_common_query_fraction':float(min(r['common_query_fraction'] for r in subset)) if subset else None}
        summaries.append(out)
    payload={'pairs':pairs,'source_parent_summary':summaries,'no_method_selected':True,'clinical_outcomes_used':False,
             'discovery_transductive':True,'routes_inherited_not_truth':True,'independent_patients_not_assumed':True,
             'matching_used_only_frozen_granularity':True,'high_coverage_threshold':.9}
    output=root/'outputs/paired_granularity_review_v1.json.gz'
    with gzip.open(output,'xt') as f:json.dump(payload,f,allow_nan=False)
    size=sum(p.stat().st_size for p in root.rglob('*') if p.is_file())
    if size>cfg['output_cap_bytes']:raise ValueError('Paired review output cap exceeded')
    save(root/'provenance/terminal_v1.json',{'status':'PAIRED_QUERY_REVIEW_AND_SIZE_SELECTED_NUMERIC_AUDIT_COMPLETE',
         'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'output_sha256':sha256_file(output),
         'n_profiles_replayed':audited_profiles,'n_pairs_and_states':len(pairs),'n_tests_actual_no_skips':4,
         'independent_audit_source_selection':'smallest_pilot_source_before_outcome_inspection',
         'independent_audit_source':audit_source,'independent_prediction_cases':independent_predictions,
         'independent_prediction_max_absolute_error':independent_max_error,
         'native_readout_values_checked':native_state_values_checked,'native_readout_max_absolute_error':native_max_error,
         'wall_seconds':time.monotonic()-started,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
         'written_bytes_before_terminal':size,'new_authoritative_labels':False,'no_graphs_or_models_refitted':True})


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',required=True);root=Path(parser.parse_args().output)
    try:main(root)
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)})
        raise
