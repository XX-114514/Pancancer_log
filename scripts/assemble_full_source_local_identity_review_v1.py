"""Join local native markers, panel/reference/author evidence for identity review.

Proposals keep their evidence scope. Group assignments, inherited labels and
RNA state summaries are distinct; no final common identity is released here.
"""
import argparse
import collections
import datetime
import gc
import hashlib
import json
from pathlib import Path
import resource
import subprocess
import sys
import time

import numpy as np
from annotation_representation_pilot_v1 import save
from audit_full_source_local_graphs_v1 import index_at
from audit_native_author_label_evidence_v1 import author_context
from extract_full_source_local_marker_evidence_v1 import exact_process_live, storage_gate
from full_context_atlas_v1 import sha256_file
from native_archive_collection_v1 import NativeArchiveCollection, row_axis_sha256, text_axis_sha256
from native_annotation_review_v1 import classify_group
from native_lineage_routing_v1 import FAMILIES
from summarize_full_native_annotation_evidence_v1 import read_source


def category_counts(values, groups, n_groups):
    values=np.asarray(values).astype(str);groups=np.asarray(groups)
    if values.ndim!=1 or groups.shape!=values.shape or n_groups<1 or not np.issubdtype(groups.dtype,np.integer) or np.any(groups<0) or np.any(groups>=n_groups):
        raise ValueError('Category source/group axes differ')
    labels,inverse=np.unique(values,return_inverse=True)
    matrix=np.bincount(groups.astype(np.int64)*len(labels)+inverse,minlength=n_groups*len(labels)).reshape(n_groups,len(labels))
    return [{str(labels[j]):int(row[j]) for j in np.flatnonzero(row)} for row in matrix]


def grouped_sums(groups, values, n_groups):
    if groups.ndim!=1 or np.asarray(values).shape!=(len(groups),) or np.any(groups<0) or np.any(groups>=n_groups):raise ValueError('Local evidence aggregation axes differ')
    return np.bincount(groups,weights=values,minlength=n_groups)


def process_source(root,cfg,reader,ds,marker_record,refs,supp,modelmeta):
    began=time.monotonic();marker=Path(cfg['markers']);ref=Path(cfg['reference']);sup=Path(cfg['supplemental']);graphroot=Path(cfg['graphs'])
    mp=marker/'outputs'/marker_record['metadata']
    if sha256_file(mp)!=marker_record['metadata_sha256']:raise ValueError('Local marker snapshot changed')
    mm=json.loads(mp.read_text());gm,ga=read_source(graphroot,mm['source_graph_record']);rm,rv=read_source(ref,refs[ds]);sm,sv=read_source(sup,supp[ds])
    rows=reader.rows(ds)
    for side in [ga,rv,sv]:np.testing.assert_array_equal(side['native_row'],rows)
    if any(m['native_ids_sha256']!=text_axis_sha256(reader.native_ids(ds)) for m in [mm,gm,rm,sm]):raise ValueError('Identity review native axis differs')
    stats_path=marker/'outputs'/mm['group_native_statistics']
    if sha256_file(stats_path)!=mm['group_native_statistics_sha256']:raise ValueError('Native group stats changed')
    with np.load(stats_path,allow_pickle=False) as f:stat={k:f[k] for k in f.files}
    selected=np.flatnonzero(ga['fit_pool_candidate']);codes=ga['source_local_group_code'][selected];ng=gm['n_source_local_groups']
    n=np.bincount(codes,minlength=ng);np.testing.assert_array_equal(n,stat['n'])
    panel_counts={name:grouped_sums(codes,(rv['native_panel_bitmask'][selected]&(1<<j))!=0,ng).astype(np.int64)
                  for j,name in enumerate(rm['native_panel_bit_order'])}
    family_counts={name:grouped_sums(codes,(rv['native_family_bitmask'][selected]&(1<<j))!=0,ng).astype(np.int64)
                  for j,name in enumerate(rm['native_family_bit_order'])}
    supplemental={name:grouped_sums(codes,sv['cue_support'][selected,j],ng).astype(np.int64) for j,name in enumerate(sm['panel_order'])}
    inherited=reader.frozen_labels(ds)
    old_major=category_counts(inherited['major'][selected],codes,ng);old_fine=category_counts(inherited['fine'][selected],codes,ng)
    field=rm['author_field'];author=reader.field(ds,field).astype(str) if field else np.full(len(rows),'',dtype='<U1')
    authors=category_counts(author[selected],codes,ng)
    names,inv=np.unique(author[selected],return_inverse=True);contexts=np.asarray([author_context(v) for v in names])[inv]
    author_coarse=category_counts(contexts,codes,ng)
    references={}
    for model in rm['source_reference_models']:
        positive=rv[model+'__n_positive_labels'][selected];unique=rv[model+'__unique_code'][selected]
        scope=rv['immune_reference_scope'][selected] if model.startswith('immune_') else rv['organ_reference_scope'][selected]
        valid=unique>=0;labels=modelmeta[model]['labels']
        if np.any(valid&(positive!=1)) or np.any((positive>=0)&~scope):raise ValueError('Reference scope or unique status changed')
        labelcounts=category_counts(np.asarray(labels)[unique[valid]],codes[valid],ng) if valid.any() else [{} for _ in range(ng)]
        references[model]={'unique_label_counts':labelcounts,
                          'n_scored':grouped_sums(codes,positive>=0,ng),'n_in_scope':grouped_sums(codes,scope,ng),
                          'n_scored_multiple':grouped_sums(codes,positive>1,ng),'n_scored_no_positive':grouped_sums(codes,positive==0,ng)}
    state=[]
    for j,name in enumerate(rm['state_order']):
        values=rv['state_native_logCP10K'][selected,j];valid=np.isfinite(values)
        ns=np.bincount(codes[valid],minlength=ng);sums=np.bincount(codes[valid],weights=values[valid],minlength=ng)
        state.append((name,ns,sums))
    unit_values=reader.frozen_expression_units(ds)[selected].astype(str)
    unit_hist=category_counts(unit_values,codes,ng)
    genes=mm['native_measured_gene_names'];gene_index={g:i for i,g in enumerate(genes)}
    review_genes=rm['review_genes'];records=[]
    for group in range(ng):
        positions=selected[codes==group];members=rows[positions]
        stable=hashlib.sha256((cfg['local_label_version']+'\0'+ds+'\0'+row_axis_sha256(members)).encode()).hexdigest()
        panel={name:float(v[group]/n[group]) for name,v in panel_counts.items()}
        detection={gene:float(stat['detected'][group,gene_index[gene]]/n[group]) for gene in review_genes if gene in gene_index}
        proposal=classify_group(panel,detection,references['immune_low']['unique_label_counts'][group],int(n[group]))
        markerinfo=mm['groups'][group]
        if markerinfo['source_local_group_code']!=group:raise ValueError('Local marker group order differs')
        refs_for_group={name:{'unique_label_counts':v['unique_label_counts'][group],**{k:int(v[k][group]) for k in ['n_scored','n_in_scope','n_scored_multiple','n_scored_no_positive']},
                             'is_identity_truth':False,'reference_training_overlap_unresolved':True} for name,v in references.items()}
        supplemental_fraction={name:float(v[group]/n[group]) for name,v in supplemental.items()}
        supplemental_review=any(v>=.2 for v in supplemental_fraction.values())
        counts={key:int(value) for key,value in collections.Counter(author_coarse[group]).items()}
        display=markerinfo['own_group_markers'][0]['gene'] if markerinfo['own_group_markers'] else None
        # Provisional display combines evidence-supported type and a native own
        # marker. Source namespace and stable key stay in dedicated fields.
        label=proposal['local_celltype']
        proposal_display=(label+'_'+display+'_'+str(group+1)) if display else label+'_marker_unresolved_'+str(group+1)
        records.append({'source_local_group_code':group,'local_subcluster_id':ds+'::'+cfg['local_label_version']+'::'+str(group+1),
          'stable_local_group_sha256':stable,'native_members_row_sha256':row_axis_sha256(members),'n_fit_cells':int(n[group]),
          'n_expression_units_not_assumed_patients':len(unit_hist[group]),'provisional_routing_family':FAMILIES[stat['provisional_family_code'][group]] if stat['provisional_family_code'][group]>=0 else 'unresolved_context',
          'local_major_proposal':proposal['local_major'],'local_celltype_proposal':proposal['local_celltype'],
          'proposal_evidence_status':proposal['annotation_evidence_status'],'proposal_confidence':proposal['confidence'],
          'native_supported_families_in_legacy_rule':proposal['native_supported_families'],
          'proposal_uses_existing_pilot_rule_not_final_taxonomy':True,'supplemental_lineage_review_required':supplemental_review,
          'supplemental_native_multigene_fraction':supplemental_fraction,'native_panel_fraction':panel,
          'native_family_fraction':{name:float(v[group]/n[group]) for name,v in family_counts.items()},
          'native_review_gene_detection':detection,'reference_proposals':refs_for_group,'native_author_label_counts':authors[group],
          'author_coarse_context_counts':counts,'author_labels_are_not_identity_truth':True,
          'inherited_major_counts':old_major[group],'inherited_fine_counts':old_fine[group],
          'marker_origin':'own_source_local_group','marker_evidence_type':'own_group_expression_specificity_not_inferential_DEG',
          'display_marker_proposal':display,'display_subcluster_proposal':proposal_display,
          'native_marker_evidence_metadata_sha256':marker_record['metadata_sha256'],
          'descriptive_RNA_states':{name:{'n_measured':int(ns[group]),'mean_native_logCP10K':float(sums[group]/ns[group]) if ns[group] else None,
                                    'measurement_coverage':rm['state_measurement_coverage'][name],'functional_identity_inferred':False} for name,ns,sums in state},
          'malignancy_from_lineage_inferred':False,'malignancy_evidence_kept_in_parallel_registry':True,
          'new_authoritative_identity':False,'cross_source_common_identity_not_assigned':True})
    output=root/'outputs'/(ds+'_local_identity_review_v1.json')
    save(output,{'source':ds,'n_archive_rows':len(rows),'n_fit_rows':len(selected),'n_unfit_archive_rows_retained':len(rows)-len(selected),
      'n_local_groups':ng,'native_ids_sha256':text_axis_sha256(reader.native_ids(ds)),'groups':records,
      'source_marker_record':marker_record,'source_graph_record':mm['source_graph_record'],
      'shared_programs_and_source_taxonomy_require_review':True,'local_display_numbers_are_source_version_scoped':True,
      'RNA_state_and_annotated_identity_separate':True,'wall_seconds':time.monotonic()-began})
    return {'source':ds,'n_archive_rows':len(rows),'n_fit_rows':len(selected),'n_local_groups':ng,'output':output.name,'output_sha256':sha256_file(output)}


def main(root):
    began=time.monotonic();cfg=json.loads((root/'config/local_identity_review_spec_v1.json').read_text())
    for binding in cfg['bindings'].values():
        if sha256_file(binding['path'])!=binding['sha256']:raise ValueError('Identity review input binding changed')
    test_specs=[(['test_local_identity_evidence_join_v1.py'],3),(['test_native_author_label_evidence_v1.py'],5),
                (['test_native_annotation_review_v1.py','NativeReviewTests.test_unmeasured_marker_and_single_gene_do_not_support_lineage',
                  'NativeReviewTests.test_conflicting_lineages_and_minority_gamma_candidate_not_whole_group_identity'],2)]
    for arguments,count in test_specs:
        test=subprocess.run([sys.executable,'-B',str(root/'tests'/arguments[0])]+arguments[1:],capture_output=True,text=True,timeout=120)
        (root/'provenance'/(arguments[0]+'.stderr.txt')).write_text(test.stderr)
        if test.returncode or ('Ran %d tests'%count) not in test.stderr or 'skipped' in test.stderr.lower():raise ValueError('Local identity evidence counterexamples failed')
    save(root/'provenance/numerical_tests_v1.json',{'actual_exit_code':0,'n_tests':10,'no_skips':True})
    marker=Path(cfg['markers']);refs=index_at(Path(cfg['reference']),'full_archive_reference_index_v1.json')
    supp=index_at(Path(cfg['supplemental']),'supplemental_native_index_v1.json')
    reader=NativeArchiveCollection(cfg['archive_manifest']);expected=set(reader.archive_partitions)
    if set(refs)!=expected or set(supp)!=expected:raise ValueError('Identity review source scopes differ')
    modelmeta=json.loads((Path(cfg['reference'])/'outputs/reference_models_v1.json').read_text());results={}
    while set(results)!=expected:
        completed=set()
        for cp in (marker/'checkpoints').glob('*.json'):
            try:v=json.loads(cp.read_text())
            except json.JSONDecodeError:
                if time.time()-cp.stat().st_mtime<5:continue
                raise
            completed.add(v['source'])
        if not completed<=expected:raise ValueError('Unregistered local marker source')
        for ds in sorted(completed-set(results),key=lambda d:(reader.archive_partitions[d]['n_archive_rows'],d)):
            storage_gate(root,cfg);mp=marker/'outputs'/(ds+'_local_native_marker_evidence_v1.json');mm=json.loads(mp.read_text())
            record={'source':ds,'n_archive_rows':mm['n_archive_rows'],'n_fitted_rows':mm['n_fitted_rows'],'metadata':mp.name,'metadata_sha256':sha256_file(mp)}
            results[ds]=process_source(root,cfg,reader,ds,record,refs,supp,modelmeta)
            save(root/'checkpoints'/(ds+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':ds,'completed_sources':len(results)})
            reader._fields.clear();gc.collect()
        if set(results)!=expected:
            terminal=marker/'provenance/supervisor_terminal_v1.json'
            if terminal.exists():
                v=json.loads(terminal.read_text())
                if v['actual_child_exit_code'] or v['wall_timeout']:raise RuntimeError('Local marker producer failed; review snapshot preserved')
            elif not exact_process_live(marker):raise RuntimeError('Local marker producer process identity unavailable')
            if time.monotonic()-began>cfg['maximum_producer_wait_seconds']:raise TimeoutError('Bounded review producer wait exhausted')
            time.sleep(15)
    while not (marker/'provenance/supervisor_terminal_v1.json').exists():
        if time.monotonic()-began>cfg['maximum_producer_wait_seconds']:raise TimeoutError('Bounded final review producer wait exhausted')
        time.sleep(5)
    sup=json.loads((marker/'provenance/supervisor_terminal_v1.json').read_text());mt=json.loads((marker/'provenance/terminal_v1.json').read_text())
    ip=marker/'outputs/full_source_local_marker_evidence_index_v1.json'
    if sup['actual_child_exit_code'] or sup['wall_timeout'] or sha256_file(ip)!=mt['index_sha256']:raise ValueError('Full marker producer not accepted as complete')
    final={r['source']:r for r in json.loads(ip.read_text())}
    for ds,r in results.items():
        out=json.loads((root/'outputs'/r['output']).read_text());bound=out['source_marker_record']
        if bound['metadata_sha256']!=final[ds]['metadata_sha256']:raise ValueError('Final marker source differs from review snapshot')
    output=root/'outputs/full_local_identity_review_index_v1.json';save(output,[results[ds] for ds in sorted(expected)])
    save(root/'provenance/terminal_v1.json',{'status':'FULL_LOCAL_IDENTITY_REVIEW_JOIN_COMPLETE_NOT_FINAL_LABELS','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
      'n_sources':len(results),'n_archive_rows':sum(r['n_archive_rows'] for r in results.values()),'n_local_groups':sum(r['n_local_groups'] for r in results.values()),
      'index_sha256':sha256_file(output),'marker_producer_terminal_sha256':sha256_file(marker/'provenance/terminal_v1.json'),
      'wall_seconds':time.monotonic()-began,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
      'new_final_identity':False,'native_expression_reread':False})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)});raise
