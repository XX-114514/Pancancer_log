"""Independent source-local descriptive marker arithmetic and native replay.

All-source array/marker/unit invariants plus one complete source count replay;
this verifies computation, not independent marker or cell-identity validation.
"""
import argparse
import datetime
import gc
import json
from pathlib import Path
import re
import resource
import subprocess
import sys
import time

import numpy as np
from annotation_representation_pilot_v1 import save
from audit_full_source_local_graphs_v1 import index_at, direct_native_rows
from full_context_atlas_v1 import sha256_file
from native_archive_collection_v1 import NativeArchiveCollection, row_axis_sha256, text_axis_sha256
from summarize_full_native_annotation_evidence_v1 import read_source


def verify_markers(meta, stat, units):
    n = stat['n']; genes = meta['native_measured_gene_names']; ng, p = len(n), len(genes)
    if p != len(set(genes)) or len(meta['groups']) != ng or np.any(n <= 0): raise ValueError('Native marker group/gene axis differs')
    for name in ['cp_sum', 'log_sum', 'detected']:
        a = stat[name]
        if a.shape != (ng, p) or not np.isfinite(a).all() or (a < 0).any(): raise ValueError('Native marker statistics invalid')
    if not np.equal(stat['detected'], np.floor(stat['detected'])).all() or np.any(stat['detected'] > n[:, None]):
        raise ValueError('Native marker detection denominator differs')
    if np.any(stat['log_sum'] > stat['cp_sum'] + 1e-5): raise ValueError('Invalid log-sum relation')
    family = stat['provisional_family_code']; index = dict(zip(genes, range(p)))
    pair_group, pair_unit, pair_n = (units[k] for k in ['pair_group_code', 'expression_unit_code', 'n_cells'])
    nu = len(meta['expression_unit_dictionary']); k = units['marker_gene_index'].shape[1]
    if np.any(pair_group < 0) or np.any(pair_group >= ng) or np.any(pair_unit < 0) or np.any(pair_unit >= nu) or np.any(pair_n <= 0):
        raise ValueError('Marker expression-unit pair axis differs')
    keys = pair_unit.astype(np.int64) * ng + pair_group
    if len(keys) != len(set(keys)) or np.any(np.diff(keys) <= 0): raise ValueError('Duplicate or unsorted marker-unit pair')
    if not np.equal(pair_n, np.floor(pair_n)).all(): raise ValueError('Noninteger unit cell count')
    np.testing.assert_array_equal(np.bincount(pair_group, weights=pair_n, minlength=ng), n)
    for name in ['cp_sum', 'detected']:
        a=units[name]
        if a.shape != (len(keys), k) or not np.isfinite(a).all() or (a < 0).any(): raise ValueError('Unit marker numeric values differ')
    if np.any(units['detected'] > pair_n[:, None]) or not np.equal(units['detected'], np.floor(units['detected'])).all():
        raise ValueError('Invalid unit marker detection count')
    own_by_group = np.zeros((ng,k)); det_by_group=np.zeros_like(own_by_group)
    np.add.at(own_by_group,pair_group,units['cp_sum']); np.add.at(det_by_group,pair_group,units['detected'])
    named_support = {(r['source_local_group_code'],r['gene']):r for r in meta['marker_unit_direction_summaries']}
    if len(named_support)!=len(meta['marker_unit_direction_summaries']): raise ValueError('Duplicate marker direction record')
    total_markers=0
    for group, info in enumerate(meta['groups']):
        if info['source_local_group_code']!=group or info['n_fit_cells']!=n[group] or info['provisional_family_code']!=family[group]:
            raise ValueError('Marker group metadata differs')
        same=family==family[group];same[group]=False
        count=int(n[same].sum());scope='other_source_local_groups_in_same_provisional_family'
        if count==0:
            same=np.arange(ng)!=group;count=int(n[same].sum());scope='other_source_local_groups_all_families_no_within_family_comparator'
        if info['n_comparator_cells']!=count or info['comparator_scope']!=scope:raise ValueError('Marker comparator differs')
        cp=stat['cp_sum'][group]/n[group];det=stat['detected'][group]/n[group]
        expected=[]
        if count:
            ocp=stat['cp_sum'][same].sum(axis=0)/count;odet=stat['detected'][same].sum(axis=0)/count
            effect=np.log2((cp+.1)/(ocp+.1));score=effect*np.maximum(det-odet,0)
            selected=[j for j,g in enumerate(genes) if det[j]>=.1 and det[j]>odet[j] and effect[j]>0 and not re.match(r'^(MT-|RPL|RPS|ENSG|IG[HKL][VDJC])',g)]
            expected=sorted(selected,key=lambda j:(-score[j],genes[j]))[:k]
        markers=info['own_group_markers']
        if [m['gene'] for m in markers]!=[genes[j] for j in expected]:raise ValueError('Own-group descriptive marker ranking differs')
        np.testing.assert_array_equal(units['marker_gene_index'][group,:len(markers)],expected)
        if np.any(units['marker_gene_index'][group,len(markers):]!=-1):raise ValueError('Missing marker slots fabricated')
        for slot,(j,m) in enumerate(zip(expected,markers)):
            for name,value in [('native_mean_CP10K',cp[j]),('native_detection_fraction',det[j]),('comparator_detection_fraction',odet[j]),('descriptive_log2_mean_ratio',effect[j])]:
                np.testing.assert_allclose(m[name],value,rtol=1e-12,atol=1e-12)
            if m['marker_evidence_type']!='own_group_expression_specificity_not_inferential_DEG':raise ValueError('Marker inferential status changed')
            np.testing.assert_allclose(own_by_group[group,slot],stat['cp_sum'][group,j],rtol=2e-10,atol=2e-7)
            np.testing.assert_array_equal(det_by_group[group,slot],stat['detected'][group,j])
            support=named_support[(group,m['gene'])];pairs=np.flatnonzero(pair_group==group)
            family_pairs=family[pair_group]==family[group]
            parent_n=np.bincount(pair_unit[family_pairs],weights=pair_n[family_pairs],minlength=nu)
            other_n=parent_n[pair_unit[pairs]]-pair_n[pairs];valid=other_n>=20
            valid &= pair_n[pairs]>=20
            if support['n_units_with_group']!=len(pairs) or support['n_comparable_expression_units']!=int(valid.sum()):raise ValueError('Unit support comparable denominator differs')
            if support['n_units_positive_mean_ratio']+support['n_units_negative_mean_ratio']>int(valid.sum()):raise ValueError('Unit support direction count exceeds domain')
            if not valid.any() and support['median_log2_mean_ratio'] is not None:raise ValueError('Noncomparable units assigned an effect')
            if not support['descriptive_not_replicate_aware_DEG'] or not support['expression_units_not_assumed_independent_patients']:raise ValueError('Unit evidence scope changed')
            total_markers+=1
        if np.any(units['cp_sum'][pair_group==group,len(markers):]) or np.any(units['detected'][pair_group==group,len(markers):]):raise ValueError('Missing marker slots have evidence')
    if len(named_support)!=total_markers:raise ValueError('Extra unit marker support record')
    return total_markers


def replay_complete_source(reader,ds,graph,reference,meta,stat,units):
    rows=reader.rows(ds);selected=np.flatnonzero(graph['fit_pool_candidate']);codes=graph['source_local_group_code'][selected]
    native=reader.store(ds);key=reader.legacy_config['counts_key'] if reader.is_legacy(ds) else reader.primary.partitions[ds]['counts_key']
    raw=direct_native_rows(native[key],rows[selected]);totals=np.asarray(raw.sum(axis=1,dtype=np.float64)).ravel()
    np.testing.assert_array_equal(totals,reference['native_total_RNA'][selected])
    x=raw[:,stat['native_columns']].astype(np.float64).toarray();det=x>0
    cp=x*(10000./totals[:,None]);log=np.log1p(cp)
    maximum=0.
    for group in range(len(stat['n'])):
        member=codes==group
        for name,value in [('cp_sum',cp[member].sum(axis=0)),('log_sum',log[member].sum(axis=0)),('detected',det[member].sum(axis=0))]:
            np.testing.assert_allclose(value,stat[name][group],rtol=2e-10,atol=2e-7)
            maximum=max(maximum,float(np.max(np.abs(value-stat[name][group]))))
    names,ucode=np.unique(reader.frozen_expression_units(ds)[selected].astype(str),return_inverse=True)
    if names.tolist()!=meta['expression_unit_dictionary']:raise ValueError('Direct unit dictionary differs')
    pg,pu=units['pair_group_code'],units['expression_unit_code'];expected_keys=ucode.astype(np.int64)*len(stat['n'])+codes
    keys,count=np.unique(expected_keys,return_counts=True)
    np.testing.assert_array_equal(keys,pu.astype(np.int64)*len(stat['n'])+pg);np.testing.assert_array_equal(count,units['n_cells'])
    for pair,(group,unit) in enumerate(zip(pg,pu)):
        member=(codes==group)&(ucode==unit);ix=units['marker_gene_index'][group];keep=ix>=0
        np.testing.assert_allclose(cp[member][:,ix[keep]].sum(axis=0),units['cp_sum'][pair,keep],rtol=2e-10,atol=2e-7)
        np.testing.assert_array_equal(det[member][:,ix[keep]].sum(axis=0),units['detected'][pair,keep])
    maximum_direction=0.
    for support in meta['marker_unit_direction_summaries']:
        group=support['source_local_group_code'];j=meta['native_measured_gene_names'].index(support['gene']);ratios=[]
        parent=stat['provisional_family_code'][codes]==stat['provisional_family_code'][group]
        for unit in range(len(names)):
            own=(codes==group)&(ucode==unit);other=parent&(codes!=group)&(ucode==unit)
            if own.sum()>=20 and other.sum()>=20:
                ratios.append(float(np.log2((cp[own,j].mean()+.1)/(cp[other,j].mean()+.1))))
        median=float(np.median(ratios)) if ratios else None
        if support['n_comparable_expression_units']!=len(ratios) or support['n_units_positive_mean_ratio']!=sum(r>0 for r in ratios) or support['n_units_negative_mean_ratio']!=sum(r<0 for r in ratios):raise ValueError('Direct unit direction replay differs')
        if median is None:
            if support['median_log2_mean_ratio'] is not None:raise ValueError('Direct missing unit effect differs')
        else:
            np.testing.assert_allclose(median,support['median_log2_mean_ratio'],rtol=2e-10,atol=2e-10)
            maximum_direction=max(maximum_direction,abs(median-support['median_log2_mean_ratio']))
    return {'source':ds,'n_fitted_native_rows_replayed':len(selected),'n_measured_genes_replayed':x.shape[1],
            'n_group_unit_pairs_replayed':len(pg),'maximum_native_stat_absolute_error':maximum,
            'maximum_unit_direction_absolute_error':maximum_direction,'not_independent_biological_validation':True}


def main(root):
    began=time.monotonic();spec=json.loads((root/'config/marker_audit_spec_v1.json').read_text());producer=Path(spec['producer'])
    config_path=producer/'config/native_marker_evidence_spec_v1.json'
    if sha256_file(config_path)!=spec['producer_config_sha256']:raise ValueError('Marker producer configuration changed')
    cfg=json.loads(config_path.read_text());request=json.loads((producer/'config/supervisor_request_v1.json').read_text())
    for path,h in request['script_SHA256'].items():
        if sha256_file(producer/path)!=h:raise ValueError('Marker producer deployed code changed')
    for binding in cfg['bindings'].values():
        if sha256_file(binding['path'])!=binding['sha256']:raise ValueError('Marker input binding changed')
    test=subprocess.run([sys.executable,'-B',str(root/'tests/test_native_marker_unit_support_v1.py')],capture_output=True,text=True,timeout=90)
    (root/'provenance/test_stderr.txt').write_text(test.stderr)
    if test.returncode or 'Ran 3 tests' not in test.stderr or 'skipped' in test.stderr.lower():raise ValueError('Marker numerical tests failed')
    save(root/'provenance/numerical_tests_v1.json',{'actual_exit_code':0,'n_tests':3,'no_skips':True})
    full=spec['scope']=='full_completed_archive'
    if full:
        sup=json.loads((producer/'provenance/supervisor_terminal_v1.json').read_text());terminal=json.loads((producer/'provenance/terminal_v1.json').read_text())
        ip=producer/'outputs/full_source_local_marker_evidence_index_v1.json'
        if sup['actual_child_exit_code'] or sup['wall_timeout'] or sha256_file(ip)!=terminal['index_sha256']:raise ValueError('Full marker producer incomplete')
        records=json.loads(ip.read_text())
    else:records=spec['completed_sources']
    reader=NativeArchiveCollection(cfg['archive_manifest']);reference=Path(cfg['reference']);ref=index_at(reference,'full_archive_reference_index_v1.json')
    if not records or len({r['source'] for r in records})!=len(records):raise ValueError('Empty/duplicate marker audit scope')
    if not {r['source'] for r in records}<=set(reader.archive_partitions):raise ValueError('Marker audit source scope differs')
    if full and {r['source'] for r in records}!=set(reader.archive_partitions):raise ValueError('Incomplete full marker audit')
    reports=[]
    for record in records:
        ds=record['source'];mp=producer/'outputs'/record['metadata']
        if sha256_file(mp)!=record['metadata_sha256']:raise ValueError('Marker metadata digest changed')
        meta=json.loads(mp.read_text());gp=Path(cfg['producer']);gm,graph=read_source(gp,meta['source_graph_record'])
        rm,rv=read_source(reference,ref[ds]);rows=reader.rows(ds)
        if meta['native_rows_sha256']!=row_axis_sha256(rows) or meta['native_ids_sha256']!=text_axis_sha256(reader.native_ids(ds)):raise ValueError('Marker native axis changed')
        np.testing.assert_array_equal(graph['native_row'],rows);np.testing.assert_array_equal(rv['native_row'],rows)
        if meta['n_archive_rows']!=len(rows) or meta['n_fitted_rows']!=int(graph['fit_pool_candidate'].sum()) or meta['n_unfitted_archive_rows_retained']!=int((~graph['fit_pool_candidate']).sum()):raise ValueError('Marker denominator changed')
        loaded=[]
        for filekey,digestkey in [('group_native_statistics','group_native_statistics_sha256'),('marker_unit_arrays','marker_unit_arrays_sha256')]:
            path=producer/'outputs'/meta[filekey]
            if sha256_file(path)!=meta[digestkey]:raise ValueError('Native statistic digest changed')
            with np.load(path,allow_pickle=False) as f:loaded.append({k:f[k] for k in f.files})
        stat,units=loaded
        genes=meta['native_measured_gene_names'];lookup=reader.feature_lookup(ds)
        if len(stat['native_columns'])!=len(genes) or any(lookup.get(g,set())!={int(c)} for g,c in zip(genes,stat['native_columns'])):raise ValueError('Marker native measured features changed')
        np.testing.assert_array_equal(stat['n'],np.bincount(graph['source_local_group_code'][graph['fit_pool_candidate']],minlength=len(stat['n'])))
        if meta['new_identity_assigned'] or meta['inferential_DEG_computed'] or meta['patient_independence_claimed'] or meta['parent_group_markers_borrowed']:raise ValueError('Marker evidence scope promoted')
        if text_axis_sha256(meta['expression_unit_dictionary'])!=meta['expression_unit_dictionary_sha256']:raise ValueError('Unit dictionary digest differs')
        count=verify_markers(meta,stat,units)
        reports.append({'source':ds,'n_archive_rows':len(rows),'n_fitted_rows':meta['n_fitted_rows'],'n_groups':len(stat['n']),
                        'n_markers_arithmetic_checked':count,'source_marker_metadata_sha256':record['metadata_sha256']})
        save(root/'checkpoints'/(ds+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':ds,'completed_sources':len(reports)})
        reader._fields.clear();del meta,stat,units,loaded,graph,rv;gc.collect()
    bounded=[r for r in records if r['n_fitted_rows']<=spec['maximum_full_replay_rows']]
    if not bounded:raise ValueError('No qualified bounded complete source for direct marker replay')
    chosen=min(bounded,key=lambda r:(r['n_fitted_rows'],r['source']));ds=chosen['source']
    meta=json.loads((producer/'outputs'/chosen['metadata']).read_text());gm,graph=read_source(Path(cfg['producer']),meta['source_graph_record']);rm,rv=read_source(reference,ref[ds])
    with np.load(producer/'outputs'/meta['group_native_statistics'],allow_pickle=False) as f:stat={k:f[k] for k in f.files}
    with np.load(producer/'outputs'/meta['marker_unit_arrays'],allow_pickle=False) as f:units={k:f[k] for k in f.files}
    replay=replay_complete_source(reader,ds,graph,rv,meta,stat,units)
    output=root/'outputs/native_local_marker_audit_v1.json';save(output,{'scope':spec['scope'],'sources':reports,'complete_source_native_replay':replay,
       'not_all_native_counts_recomputed':True,'new_final_identity':False})
    result={'status':'FULL_LOCAL_MARKER_AUDIT_COMPLETE' if full else 'COMPLETED_SOURCE_MARKER_PREFLIGHT_ONLY','scope':spec['scope'],
      'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_sources':len(reports),
      'n_groups_checked':sum(r['n_groups'] for r in reports),'n_markers_arithmetic_checked':sum(r['n_markers_arithmetic_checked'] for r in reports),
      'output_sha256':sha256_file(output),'wall_seconds':time.monotonic()-began,
      'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,'complete_source_native_replay':replay}
    if full:result['producer_terminal_sha256']=sha256_file(producer/'provenance/terminal_v1.json')
    save(root/'provenance/terminal_v1.json',result)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)});raise
