"""Frozen two-panel sensitivity; no cell removal or authoritative relabeling.

Globins alone and the S100A8/A9 pair are descriptive RNA evidence. Test whether
requiring an additional lineage anchor changes review flags. This is not an
ambient-RNA/doublet diagnosis and not a test of biological annotation accuracy.
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

from annotation_representation_pilot_v1 import save
from full_context_atlas_v1 import FullContextAtlas,sha256_file
from native_annotation_review_v1 import REVIEW_GENES,marker_calls,classify_group
from native_panel_aggregate_v1 import gather_csr_rows

EXTRA_GENES=['GYPA','GYPB','ALAS2','SLC4A1']
ERYTHROID_ANCHORS=['AHSP',*EXTRA_GENES]
NEUTROPHIL_ANCHORS=['CSF3R','FCGR3B']


def anchored_calls(detected,genes,measured,old):
    index={g:i for i,g in enumerate(genes)}
    def any_gene(names):
        cols=[index[g] for g in names if g in index and measured[index[g]]]
        return detected[:,cols].any(axis=1) if cols else np.zeros(len(detected),bool)
    # Restrict original support; additional genes cannot create new positive cells.
    return {'Erythroid_related':old['Erythroid_related'] & any_gene(ERYTHROID_ANCHORS),
            'Neutrophil_related':old['Neutrophil_related'] & any_gene(NEUTROPHIL_ANCHORS)}


def main(root,producer,audit):
    started=time.monotonic()
    if json.loads((audit/'provenance/supervisor_terminal_v1.json').read_text())['actual_child_exit_code']!=0:
        raise ValueError('Native review numeric audit incomplete')
    cfg=json.loads((producer/'config/review_spec_v1.json').read_text())
    for binding in cfg['bindings'].values():
        if sha256_file(binding['path'])!=binding['sha256']:raise ValueError('Frozen annotation input changed')
    test=subprocess.run([sys.executable,'-B',str(root/'tests/test_lineage_anchor_sensitivity_v1.py')],capture_output=True,text=True,timeout=120)
    (root/'provenance/test_stdout.txt').write_text(test.stdout);(root/'provenance/test_stderr.txt').write_text(test.stderr)
    if test.returncode!=0 or 'Ran 3 tests' not in test.stderr or 'skipped' in test.stderr.lower():
        raise ValueError('Lineage-anchor counterexamples failed')
    reader=FullContextAtlas(cfg['manifest']);pilot=Path(cfg['pilot']);discovery=Path(cfg['discovery'])
    inventory=json.loads((pilot/'outputs/pilot_input_inventory_v1.json').read_text())
    review_index=json.loads((producer/'outputs/source_annotation_review_index_v1.json').read_text())
    with np.load(discovery/'outputs'/('candidate_partitions_'+str(cfg['seed'])+'.npz'),allow_pickle=False) as f:parts={m:f[m] for m in ['A','B','C']}
    genes=sorted(set(REVIEW_GENES+EXTRA_GENES));offset=0;summaries=[]
    for inv,entry in zip(inventory,review_index):
        ds=inv['dataset'];selector=pilot/inv['selector_file']
        if ds!=entry['source'] or sha256_file(selector)!=inv['selector_sha256']:raise ValueError('Source or selector changed')
        path=producer/'outputs'/entry['output']
        if sha256_file(path)!=entry['output_sha256']:raise ValueError('Review changed')
        with gzip.open(path,'rt') as f:review=json.load(f)
        positions=np.load(selector,allow_pickle=False);n=len(positions);lookup=reader.feature_lookup(ds)
        measured=np.array([len(lookup.get(g,set()))==1 for g in genes]);columns=np.array([next(iter(lookup[g])) if measured[i] else 0 for i,g in enumerate(genes)])
        store=reader.store(ds);key=reader.legacy_config['counts_key'] if reader.is_legacy(ds) else reader.primary.partitions[ds]['counts_key']
        oldparts=collections.defaultdict(list);newparts=collections.defaultdict(list);positive=[]
        for lo in range(0,n,2048):
            raw=gather_csr_rows(store[key],reader.rows(ds)[positions[lo:lo+2048]])
            if not np.isfinite(raw.data).all() or (raw.data<0).any() or not np.equal(raw.data,np.floor(raw.data)).all():raise ValueError('Native count domain differs')
            x=raw[:,columns].toarray();x[:,~measured]=0;detected=x>0
            old=marker_calls(detected,genes,measured);new=anchored_calls(detected,genes,measured,old)
            for k,v in old.items():oldparts[k].append(v)
            for k,v in new.items():newparts[k].append(v)
            positive.append(np.asarray(raw.sum(axis=1)).ravel()>0)
        old={k:np.concatenate(v) for k,v in oldparts.items()};new={k:np.concatenate(v) for k,v in newparts.items()};positive=np.concatenate(positive)
        records=[];summary=collections.defaultdict(collections.Counter)
        transitions=collections.defaultdict(collections.Counter)
        for group in review['groups']:
            member=(parts[group['method']][offset:offset+n]==group['partition_code'])&positive;nrna=int(member.sum())
            if nrna!=group['n_native_RNA_evaluable']:raise ValueError('Group RNA denominator differs')
            fractions={k:float(v[member].mean()) if nrna else 0. for k,v in old.items()}
            for k,v in fractions.items():
                if not np.isclose(v,group['native_multi_marker_fraction'][k],atol=1e-12,rtol=0):raise ValueError('Original panel fraction differs')
            scenarios={}
            for scenario,keys in [('original',[]),('erythroid_anchor',['Erythroid_related']),('neutrophil_anchor',['Neutrophil_related']),('both_anchors',list(new))]:
                frac=dict(fractions)
                for k in keys:frac[k]=float(new[k][member].mean()) if nrna else 0.
                call=classify_group(frac,group['native_review_gene_detection'],group['reference_proposals']['immune_low']['unique_label_counts'],nrna)
                if group['partition_code']<0:call.update(local_major='unresolved',local_celltype='unresolved',confidence='unresolved',annotation_evidence_status='no_frozen_candidate_partition')
                if scenario=='original':
                    for k in ['local_major','local_celltype','annotation_evidence_status','native_supported_families']:
                        if call[k]!=group[k]:raise ValueError('Original group call not reproduced')
                scenarios[scenario]={k:call[k] for k in ['local_major','local_celltype','annotation_evidence_status','native_supported_families']}
                summary[(group['method'],scenario)][call['annotation_evidence_status']]+=group['n_pilot_cells']
                transitions[(group['method'],scenario)][group['local_major']+' -> '+call['local_major']]+=group['n_pilot_cells']
            records.append({'method':group['method'],'source_local_group_id':group['source_local_group_id'],'n_pilot_cells':group['n_pilot_cells'],
                            'scenarios':scenarios,'candidate_rule_changes_not_authoritative_labels':True})
        coverage={p:{'measured':[g for g in gg if measured[genes.index(g)]],'unmeasured_or_ambiguous':[g for g in gg if not measured[genes.index(g)]]}
                  for p,gg in [('erythroid_anchors',ERYTHROID_ANCHORS),('neutrophil_anchors',NEUTROPHIL_ANCHORS)]}
        source={'source':ds,'n_pilot_cells':n,'anchor_measurement':coverage,
                'cell_support':{k:{'original_positive':int(old[k].sum()),'anchor_supported':int(v.sum()),
                                   'original_without_anchor':int((old[k]&~v).sum())} for k,v in new.items()},
                'scenario_status_cell_counts':[{'method':m,'scenario':s,'status_counts':dict(v),'major_transitions':dict(transitions[(m,s)])}
                                               for (m,s),v in summary.items()],
                'all_original_group_calls_reproduced':True,'unresolved_reduction_is_not_accuracy':True,
                'not_ambient_RNA_or_doublet_diagnosis':True,'no_cells_removed':True,'new_labels_released':False}
        out=root/'outputs'/(ds+'_anchor_sensitivity_v1.json.gz')
        with gzip.open(out,'xt') as f:json.dump({'source_summary':source,'groups':records},f,allow_nan=False)
        summaries.append({**source,'private_output_sha256':sha256_file(out)})
        offset+=n;save(root/'checkpoints'/(ds+'.json'),{'source':ds,'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat()})
    if any(offset!=len(v) for v in parts.values()):raise ValueError('Incomplete source scope')
    output=root/'outputs/source_anchor_sensitivity_summary_v1.json';save(output,summaries)
    size=sum(p.stat().st_size for p in root.rglob('*') if p.is_file())
    if size>16*1024**2:raise ValueError('Anchor sensitivity output cap exceeded')
    save(root/'provenance/terminal_v1.json',{'status':'TWO_PANEL_ANCHOR_SENSITIVITY_COMPLETE_NOT_NEW_IDENTITY_OR_QC',
         'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'output_sha256':sha256_file(output),
         'n_sources':len(summaries),'n_pilot_cells':offset,'n_tests_actual_no_skips':3,
         'all_original_group_calls_reproduced':True,'wall_seconds':time.monotonic()-started,
         'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,'written_bytes_before_terminal':size})


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    for arg in ['output','producer','audit']:parser.add_argument('--'+arg,required=True)
    args=parser.parse_args();root=Path(args.output)
    try:main(root,Path(args.producer),Path(args.audit))
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)})
        raise
