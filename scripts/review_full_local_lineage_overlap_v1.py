"""All cross-lineage group pairs: same-cell support versus segregated support.

Descriptive, after inspecting broad-panel conflict coverage. No new cutoff,
identity reassignment, doublet diagnosis or statistical significance is used.
"""
import argparse
import collections
import datetime
import gc
import itertools
import json
from pathlib import Path
import resource
import subprocess
import sys
import time

import numpy as np
from annotation_representation_pilot_v1 import save
from audit_full_source_local_graphs_v1 import index_at
from full_context_atlas_v1 import sha256_file
from lineage_support_contingency_v1 import contingency
from summarize_full_native_annotation_evidence_v1 import read_source

PANELS={'T':['T'],'NK':['NK_related'],'B':['B'],'Plasma':['Plasma'],
 'Myeloid':['Myeloid','DC_related','Neutrophil_related'],'Mast':['Mast'],
 'Endothelial':['Endothelial','Lymphatic_related'],'Fibroblast':['Fibroblast'],'Mural':['Mural'],
 'Epithelial':['Epithelial'],'Glial':['Glial_Oligodendrocyte_related'],
 'Astrocyte':['Astrocyte_related'],'Neuron':['Neuron_related'],'Erythroid':['Erythroid_related']}


def main(root):
    began=time.monotonic();cfg=json.loads((root/'config/local_lineage_overlap_spec_v1.json').read_text())
    for binding in cfg['bindings'].values():
        if sha256_file(binding['path'])!=binding['sha256']:raise ValueError('Lineage overlap input binding changed')
    test=subprocess.run([sys.executable,'-B',str(root/'tests/test_lineage_support_contingency_v1.py')],capture_output=True,text=True,timeout=90)
    (root/'provenance/test_stderr.txt').write_text(test.stderr)
    if test.returncode or 'Ran 3 tests' not in test.stderr or 'skipped' in test.stderr.lower():raise ValueError('Lineage overlap counterexamples failed')
    save(root/'provenance/numerical_tests_v1.json',{'actual_exit_code':0,'n_tests':3,'no_skips':True})
    review=Path(cfg['review']);graphs=Path(cfg['graphs']);reference=Path(cfg['reference'])
    index=index_at(review,'full_local_identity_review_index_v1.json');gi=index_at(graphs,'full_native_source_local_graph_index_v1.json');ri=index_at(reference,'full_archive_reference_index_v1.json')
    if set(index)!=set(gi) or set(index)!=set(ri):raise ValueError('Lineage overlap source scope differs')
    reports=[]
    for ds in sorted(index):
        record=index[ds];path=review/'outputs'/record['output']
        if sha256_file(path)!=record['output_sha256']:raise ValueError('Review source changed')
        content=json.loads(path.read_text());gm,ga=read_source(graphs,gi[ds]);rm,rv=read_source(reference,ri[ds])
        np.testing.assert_array_equal(ga['native_row'],rv['native_row'])
        if content['native_ids_sha256']!=gm['native_ids_sha256'] or gm['native_ids_sha256']!=rm['native_ids_sha256']:raise ValueError('Lineage overlap native axes differ')
        order=rm['native_panel_bit_order'];support={}
        for name,panels in PANELS.items():
            if not set(panels)<=set(order):raise ValueError('Frozen lineage panel absent')
            mask=sum(1<<order.index(panel) for panel in panels)
            support[name]=(rv['native_panel_bitmask']&mask)!=0
        groups=[];aggregates={};unique_conflict_cells=0
        for g in content['groups']:
            if g['proposal_evidence_status']!='cross_lineage_marker_conflict_review':continue
            code=g['source_local_group_code'];member=(ga['source_local_group_code']==code)&ga['fit_pool_candidate'];n=int(member.sum())
            if n!=g['n_fit_cells']:raise ValueError('Lineage overlap group denominator differs')
            unique_conflict_cells+=n;pairs=[]
            for first,second in itertools.combinations(sorted(g['native_supported_families_in_legacy_rule']),2):
                stats=contingency(support[first][member],support[second][member]);key=first+' | '+second
                pairs.append({'first':first,'second':second,**stats})
                if key not in aggregates:aggregates[key]=collections.Counter()
                for name in ['n_cells','n_first_only','n_second_only','n_both','n_neither','n_first_supported','n_second_supported']:aggregates[key][name]+=stats[name]
                aggregates[key]['n_local_groups']+=1
            groups.append({'source_local_group_code':code,'stable_local_group_sha256':g['stable_local_group_sha256'],
              'n_fit_cells':n,'original_rule_supported_families':g['native_supported_families_in_legacy_rule'],'pairs':pairs,
              'new_identity_assigned':False})
        output=root/'outputs'/(ds+'_local_lineage_overlap_v1.json')
        save(output,{'source':ds,'review_output_sha256':record['output_sha256'],'n_local_groups':len(content['groups']),
          'n_conflict_groups':len(groups),'n_unique_cells_in_conflict_groups':unique_conflict_cells,'groups':groups,
          'all_supported_family_pairs_included':True,'not_replicate_level_statistical_inference':True,
          'pair_aggregate_denominators_repeated_across_different_pairs':True,'family_union_support_can_exceed_individual_panel_support':True,
          'no_threshold_or_route_changed':True})
        reports.append({'source':ds,'n_archive_rows':record['n_archive_rows'],'n_conflict_groups':len(groups),
          'n_unique_cells_in_conflict_groups':unique_conflict_cells,'pair_aggregates':{k:dict(v) for k,v in aggregates.items()},
          'output':output.name,'output_sha256':sha256_file(output)})
        save(root/'checkpoints'/(ds+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':ds,'completed_sources':len(reports)})
        if sum(p.stat().st_size for p in root.rglob('*') if p.is_file())>cfg['output_cap_bytes']:raise ValueError('Lineage overlap output cap exceeded')
        del ga,rv,content,support;gc.collect()
    output=root/'outputs/full_local_lineage_overlap_index_v1.json';save(output,reports)
    save(root/'provenance/terminal_v1.json',{'status':'FULL_LOCAL_LINEAGE_OVERLAP_DESCRIPTIVE_REVIEW_COMPLETE_NOT_RELABELING',
      'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_sources':len(reports),
      'n_conflict_groups':sum(r['n_conflict_groups'] for r in reports),'n_unique_cells_in_conflict_groups':sum(r['n_unique_cells_in_conflict_groups'] for r in reports),
      'index_sha256':sha256_file(output),'review_producer_terminal_sha256':sha256_file(review/'provenance/terminal_v1.json'),
      'wall_seconds':time.monotonic()-began,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
      'new_final_identity':False,'native_counts_reread':False,'new_statistical_significance_tests':False})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)});raise
