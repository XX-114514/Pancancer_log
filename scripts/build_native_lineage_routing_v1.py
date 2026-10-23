"""Versioned pilot coarse routing; native counts are read-only and RAM-only."""
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
from native_annotation_review_v1 import REVIEW_GENES,marker_calls
from native_panel_aggregate_v1 import gather_csr_rows
from review_lineage_anchor_sensitivity_v1 import EXTRA_GENES,anchored_calls
from native_lineage_routing_v1 import FAMILIES,STATUS,UNRESOLVED,reference_family,local_family,native_family_support,route_cells


def main(root):
    started=time.monotonic();cfg=json.loads((root/'config/native_routing_spec_v1.json').read_text())
    for binding in cfg['bindings'].values():
        if sha256_file(binding['path'])!=binding['sha256']:raise ValueError('Frozen routing input changed')
    producer=Path(cfg['annotation']);anchor=Path(cfg['anchor']);pilot=Path(cfg['pilot']);discovery=Path(cfg['discovery'])
    for p in [producer,anchor]:
        if json.loads((p/'provenance/supervisor_terminal_v1.json').read_text())['actual_child_exit_code']!=0:raise ValueError('Dependent producer incomplete')
    test=subprocess.run([sys.executable,'-B',str(root/'tests/test_native_lineage_routing_v1.py')],capture_output=True,text=True,timeout=120)
    (root/'provenance/test_stdout.txt').write_text(test.stdout);(root/'provenance/test_stderr.txt').write_text(test.stderr)
    if test.returncode!=0 or 'Ran 6 tests' not in test.stderr or 'skipped' in test.stderr.lower():raise ValueError('Routing counterexamples failed')
    oldcfg=json.loads((producer/'config/review_spec_v1.json').read_text())
    for binding in oldcfg['bindings'].values():
        if sha256_file(binding['path'])!=binding['sha256']:raise ValueError('Frozen native review input changed')
    models=json.loads((producer/'outputs/reference_models_v1.json').read_text())
    maps={model:{label:reference_family(model,label) for label in meta['labels']} for model,meta in models.items()}
    save(root/'outputs/reference_family_dictionary_v1.json',{'maps':maps,'scope':'coarse_proposals_not_independent_truth',
          'ignored_reference_labels_not_negative_evidence':True,'normal_organ_reference_does_not_determine_malignancy':True})
    reader=FullContextAtlas(cfg['manifest']);inventory=json.loads((pilot/'outputs/pilot_input_inventory_v1.json').read_text())
    index=json.loads((producer/'outputs/source_annotation_review_index_v1.json').read_text())
    anchor_index={r['source']:r for r in json.loads((anchor/'outputs/source_anchor_sensitivity_summary_v1.json').read_text())}
    with np.load(discovery/'outputs/candidate_partitions_2026102301.npz',allow_pickle=False) as f:local_parts=f['A']
    with np.load(discovery/'outputs/derived_axis_and_model_mask.npz',allow_pickle=False) as f:source,inherited,modeled=f['source'],f['inherited_route'],f['modeled']
    genes=sorted(set(REVIEW_GENES+EXTRA_GENES));offset=0;summaries=[];side_parts=collections.defaultdict(list)
    for inv,entry in zip(inventory,index):
        ds=inv['dataset'];selector=pilot/inv['selector_file']
        if ds!=entry['source'] or sha256_file(selector)!=inv['selector_sha256']:raise ValueError('Source/selector changed')
        positions=np.load(selector,allow_pickle=False);n=len(positions);sl=slice(offset,offset+n)
        if not np.all(source[sl]==ds):raise ValueError('Source axis differs')
        path=producer/'outputs'/entry['output']
        if sha256_file(path)!=entry['output_sha256']:raise ValueError('Native proposal table changed')
        with gzip.open(path,'rt') as f:review=json.load(f)
        anchorpath=anchor/'outputs'/(ds+'_anchor_sensitivity_v1.json.gz')
        if sha256_file(anchorpath)!=anchor_index[ds]['private_output_sha256']:raise ValueError('Anchor group review changed')
        with gzip.open(anchorpath,'rt') as f:anchors=json.load(f)
        # Audit every source A context against its frozen proposed group axis.
        group_context={r['source_local_group_id']:r['scenarios']['both_anchors'] for r in anchors['groups'] if r['method']=='A'}
        local_code=np.full(n,-1,np.int16);covered=np.zeros(n,bool)
        for record in review['groups']:
            if record['method']!='A':continue
            member=local_parts[sl]==record['partition_code']
            if member.sum()!=record['n_pilot_cells'] or covered[member].any():raise ValueError('Local group context axis invalid')
            family=local_family(group_context[record['source_local_group_id']]['local_major'])
            if family:local_code[member]=FAMILIES.index(family)
            covered|=member
        if not covered.all():raise ValueError('Missing local context')
        lookup=reader.feature_lookup(ds);measured=np.array([len(lookup.get(g,set()))==1 for g in genes])
        cols=np.array([next(iter(lookup[g])) if measured[j] else 0 for j,g in enumerate(genes)])
        store=reader.store(ds);key=reader.legacy_config['counts_key'] if reader.is_legacy(ds) else reader.primary.partitions[ds]['counts_key']
        blocks=[];positive=[]
        for lo in range(0,n,2048):
            raw=gather_csr_rows(store[key],reader.rows(ds)[positions[lo:lo+2048]])
            if not np.isfinite(raw.data).all() or (raw.data<0).any() or not np.equal(raw.data,np.floor(raw.data)).all():raise ValueError('Native count domain differs')
            x=raw[:,cols].toarray();x[:,~measured]=0;calls=marker_calls(x>0,genes,measured)
            calls.update(anchored_calls(x>0,genes,measured,calls));family=native_family_support(calls)
            blocks.append(np.column_stack([family[k] for k in FAMILIES]));positive.append(np.asarray(raw.sum(axis=1)).ravel()>0)
        native=np.vstack(blocks);positive=np.concatenate(positive)
        if not np.array_equal(positive,modeled[sl]):raise ValueError('RNA eligibility changed')
        reference=np.zeros_like(native);path=producer/'outputs'/(ds+'_reference_proposals_v1.npz')
        if sha256_file(path)!=entry['reference_proposal_sha256']:raise ValueError('Reference sidecar changed')
        model_status={}
        with np.load(path,allow_pickle=False) as f:
            for model in ['immune_low','immune_high',oldcfg['source_organ_model'][ds]]:
                codes=f[model+'__unique_code'];npositive=f[model+'__n_positive_labels']
                if len(codes)!=n or not np.array_equal(codes>=0,npositive==1):raise ValueError('Reference unique-label axis invalid')
                families=np.array([-1 if maps[model][label] is None else FAMILIES.index(maps[model][label]) for label in models[model]['labels']])
                valid=codes>=0;mapped=np.full(n,-1,np.int16);mapped[valid]=families[codes[valid]]
                for j in range(len(FAMILIES)):reference[:,j]|=mapped==j
                model_status[model]={'unique_in_scope_family':int((mapped>=0).sum()),'unique_ignored_family':int((valid&(mapped<0)).sum()),
                                     'unassigned':int((npositive==0).sum()),'multiple':int((npositive>1).sum()),'outside_evaluated_scope':int((npositive<0).sum())}
        codes,reasons=route_cells(native,reference,local_code,positive)
        native_mask=(native.astype(np.int32)*(1<<np.arange(len(FAMILIES)))).sum(axis=1).astype(np.int32)
        reference_mask=(reference.astype(np.int32)*(1<<np.arange(len(FAMILIES)))).sum(axis=1).astype(np.int32)
        route=np.array([FAMILIES[k] if k>=0 else UNRESOLVED for k in codes])
        mismatch=(codes>=0)&(local_code>=0)&(codes!=local_code)
        summary={'source':ds,'n_pilot_cells':n,'coarse_route_counts':dict(collections.Counter(route)),
                 'reason_counts':{name:int((reasons==j).sum()) for j,name in enumerate(STATUS)},
                 'n_routed_cells_differing_from_group_context':int(mismatch.sum()),
                 'reference_scope_counts':model_status,'reference_models_are_not_independent_votes':True,
                 'author_labels_unchanged':True,'provisional_routing_not_final_identity':True,
                 'anchor_source_output_sha256':sha256_file(anchorpath),
                 'inherited_route_crosswalk':[{'inherited':str(a),'new_provisional':str(b),'n_cells':int(v)}
                      for (a,b),v in collections.Counter(zip(inherited[sl],route)).items()]}
        summaries.append(summary)
        for k,v in [('route_code',codes),('reason_code',reasons),('local_context_code',local_code),
                    ('native_family_bitmask',native_mask),('reference_family_bitmask',reference_mask)]:side_parts[k].append(v)
        offset+=n;save(root/'checkpoints'/(ds+'.json'),{'source':ds,'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat()})
    if offset!=len(source):raise ValueError('Incomplete routing scope')
    arrays={k:np.concatenate(v) for k,v in side_parts.items()};arrays.update(source=source,inherited_route=inherited,positive_native_RNA=modeled)
    side=root/'outputs/provisional_native_routes_v1.npz'
    with side.open('xb') as f:np.savez_compressed(f,**arrays)
    with np.load(side,allow_pickle=False) as f:
        for k,v in arrays.items():np.testing.assert_array_equal(f[k],v)
    save(root/'outputs/source_routing_summary_v1.json',summaries)
    save(root/'outputs/routing_dictionary_v1.json',{'family_codes':FAMILIES,'unresolved_code':-1,'reason_codes':STATUS,
          'native_only_is_low_confidence_provisional':True,'native_dropout_requires_reference_and_local_context':True,
          'malignancy_not_inferred':True,'new_final_authority':False})
    size=sum(p.stat().st_size for p in root.rglob('*') if p.is_file())
    if size>cfg['output_cap_bytes']:raise ValueError('Routing output cap exceeded')
    save(root/'provenance/terminal_v1.json',{'status':'SOURCE_NATIVE_COARSE_ROUTING_PROPOSALS_COMPLETE_NOT_FINAL_IDENTITY',
         'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_sources':len(summaries),'n_pilot_cells':offset,
         'n_tests_actual_no_skips':6,'output_sha256':sha256_file(side),'summary_sha256':sha256_file(root/'outputs/source_routing_summary_v1.json'),
         'wall_seconds':time.monotonic()-started,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
         'written_bytes_before_terminal':size,'new_final_identity':False,'clinical_outcomes_used':False,'native_counts_copied':False})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)})
        raise
