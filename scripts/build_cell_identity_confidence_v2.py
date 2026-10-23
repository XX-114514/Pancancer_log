"""Apply the same guarded identity policy to frozen A/B/C pilot partitions.

The source-private output is a categorical sidecar. Counts, graphs, original
labels and clinical outcomes are untouched. Marker references explicitly point
to parent discovery groups, pending markers for any newly guarded subset.
"""
import argparse
import collections
import datetime
import gzip
import hashlib
import json
from pathlib import Path
import resource
import subprocess
import sys
import time
import numpy as np

from annotation_representation_pilot_v1 import save
from full_context_atlas_v1 import FullContextAtlas, sha256_file
from native_panel_aggregate_v1 import gather_csr_rows
from native_annotation_review_v1 import REVIEW_GENES, marker_calls
from review_lineage_anchor_sensitivity_v1 import EXTRA_GENES, anchored_calls
from native_lineage_routing_v1 import FAMILIES, native_family_support, reference_family
from cell_identity_confidence_v2 import FINE_TYPES, FINE_PANEL, low_reference_identity, decide_cell


def frequencies(x):
    return {str(k): int(v) for k, v in sorted(collections.Counter(x).items())}


def bits(x):
    return (np.asarray(x, np.uint32) * (1 << np.arange(x.shape[1], dtype=np.uint32))).sum(axis=1).astype(np.uint32)


def categories(values):
    labels, codes = np.unique(np.asarray(values).astype(str), return_inverse=True)
    if len(labels) > 32767:raise ValueError('Categorical code capacity exceeded')
    return codes.astype(np.int16), labels.tolist()


def main(root):
    started = time.monotonic();cfg = json.loads((root/'config/cell_confidence_spec_v2.json').read_text())
    for binding in cfg['bindings'].values():
        if sha256_file(binding['path']) != binding['sha256']:raise ValueError('Confidence input changed')
    for stage in cfg['accepted_dependencies']:
        sup = json.loads((Path(stage)/'provenance/supervisor_terminal_v1.json').read_text())
        if sup['actual_child_exit_code'] != 0 or sup['wall_timeout']:raise ValueError('Dependency acceptance incomplete')
    t = subprocess.run([sys.executable,'-B',str(root/'tests/test_cell_identity_confidence_v2.py')],capture_output=True,text=True,timeout=120)
    (root/'provenance/test_stdout.txt').write_text(t.stdout);(root/'provenance/test_stderr.txt').write_text(t.stderr)
    if t.returncode != 0 or 'Ran 10 tests' not in t.stderr or 'skipped' in t.stderr.lower():raise ValueError('Cell identity guard counterexamples failed')
    save(root/'provenance/numerical_tests_v1.json',{'actual_exit_code':t.returncode,'n_tests':10,'no_skips':True})
    producer = Path(cfg['annotation']);pilot = Path(cfg['pilot']);discovery = Path(cfg['discovery']);routing = Path(cfg['routing'])
    nc = json.loads((producer/'config/review_spec_v1.json').read_text());models = json.loads((producer/'outputs/reference_models_v1.json').read_text())
    refmaps = {m:[reference_family(m,label) for label in x['labels']] for m,x in models.items()}
    lowmaps = [low_reference_identity(label) for label in models['immune_low']['labels']]
    inventory = json.loads((pilot/'outputs/pilot_input_inventory_v1.json').read_text())
    index = json.loads((producer/'outputs/source_annotation_review_index_v1.json').read_text())
    with np.load(discovery/'outputs/derived_axis_and_model_mask.npz',allow_pickle=False) as f:
        source,positive_all,route = f['source'],f['modeled'],f['provisional_native_route']
        fine_fit = {'A':f['A_fine_discovery_eligible'],'B':positive_all,'C':f['C_fine_discovery_eligible']}
    with np.load(discovery/'outputs/candidate_partitions_2026102301.npz',allow_pickle=False) as f:parts = {m:f[m] for m in ['A','B','C']}
    with np.load(routing/'outputs/provisional_native_routes_v1.npz',allow_pickle=False) as f:
        previous_native = f['native_family_bitmask'];previous_reference = f['reference_family_bitmask']
        if not np.array_equal(f['source'],source):raise ValueError('Routing source axis differs')
    reader = FullContextAtlas(cfg['manifest']);genes = sorted(set(REVIEW_GENES+EXTRA_GENES));offset = 0;outputs=[];summaries=[]
    for inv, entry in zip(inventory,index):
        ds = inv['dataset'];selector = pilot/inv['selector_file'];n = inv['n_selected_cells'];sl = slice(offset,offset+n)
        if ds != entry['source'] or sha256_file(selector) != inv['selector_sha256']:raise ValueError('Source selector changed')
        positions = np.load(selector,allow_pickle=False)
        if len(positions) != n or not np.all(source[sl] == ds):raise ValueError('Confidence source axis differs')
        rp = producer/'outputs'/entry['output'];sp = producer/'outputs'/(ds+'_reference_proposals_v1.npz')
        if sha256_file(rp) != entry['output_sha256'] or sha256_file(sp) != entry['reference_proposal_sha256']:raise ValueError('Accepted native proposals changed')
        with gzip.open(rp,'rt') as f:review = json.load(f)
        with np.load(sp,allow_pickle=False) as f:reference_side = {k:f[k] for k in f.files}
        lookup = reader.feature_lookup(ds);measured = np.array([len(lookup.get(g,set()))==1 for g in genes])
        cols = np.array([next(iter(lookup[g])) if measured[j] else 0 for j,g in enumerate(genes)])
        z=reader.store(ds);count_key=reader.legacy_config['counts_key'] if reader.is_legacy(ds) else reader.primary.partitions[ds]['counts_key']
        calls_parts=collections.defaultdict(list);positive=[];cd4=[]
        for lo in range(0,n,2048):
            raw=gather_csr_rows(z[count_key],reader.rows(ds)[positions[lo:lo+2048]])
            x=raw[:,cols].toarray();x[:,~measured]=0;det=x>0
            calls=marker_calls(det,genes,measured);calls.update(anchored_calls(det,genes,measured,calls))
            for key,values in calls.items():calls_parts[key].append(values)
            positive.append(np.asarray(raw.sum(axis=1,dtype=np.float64)).ravel()>0)
            cd4.append(det[:,genes.index('CD4')] & calls['T'])
        calls={k:np.concatenate(v) for k,v in calls_parts.items()};positive=np.concatenate(positive)
        np.testing.assert_array_equal(positive,positive_all[sl])
        native_dict=native_family_support(calls);native=np.column_stack([native_dict[k] for k in FAMILIES]);native_bits=bits(native)
        np.testing.assert_array_equal(native_bits,previous_native[sl])
        own=np.column_stack([np.concatenate(cd4) if t=='CD4_T' else calls[FINE_PANEL[t]] for t in FINE_TYPES])
        reference=np.zeros_like(native);low_codes=reference_side['immune_low__unique_code'];low_major=[];low_fine=[]
        for code in low_codes:
            major, fine = lowmaps[int(code)] if code >= 0 else (None,set())
            low_major.append(major);low_fine.append(fine)
        for model in ['immune_low','immune_high',nc['source_organ_model'][ds]]:
            codes=reference_side[model+'__unique_code'];npos=reference_side[model+'__n_positive_labels']
            if len(codes)!=n or not np.array_equal(codes>=0,npos==1):raise ValueError('Reference proposal axis changed')
            for j,family in enumerate(FAMILIES):
                accepted_codes=[k for k,v in enumerate(refmaps[model]) if v==family]
                reference[:,j]|=np.isin(codes,accepted_codes)
        native_sets=[frozenset(np.asarray(FAMILIES)[row]) for row in native]
        reference_sets=[frozenset(np.asarray(FAMILIES)[row]) for row in reference]
        own_sets=[frozenset(np.asarray(FINE_TYPES)[row]) for row in own]
        native_major=[]
        for i,families in enumerate(native_sets):
            family=next(iter(families)) if len(families)==1 else None
            if family=='T_NK':native_major.append('T' if calls['T'][i] else 'NK')
            elif family=='B_Plasma':native_major.append('Plasma' if calls['Plasma'][i] else 'B')
            else:native_major.append(family)
        # Full source group aggregates reproduce the already audited input panels.
        for g in review['groups']:
            member=(parts[g['method']][sl]==g['partition_code']) & positive
            if int(member.sum()) != g['n_native_RNA_evaluable']:raise ValueError('Guard input group denominator differs')
            for panel,old in g['native_multi_marker_fraction'].items():
                actual=float(calls[panel][member].mean()) if member.any() else 0.
                if abs(actual-old)>1e-12:raise ValueError('Guard native panel differs from accepted group review')
        out={'native_row':reader.rows(ds)[positions],'positive_native_RNA':positive,'native_family_bitmask':native_bits,
             'reference_family_bitmask':bits(reference),'own_fine_bitmask':bits(own),
             'reference_fine_bitmask':bits(np.array([[t in s for t in FINE_TYPES] for s in low_fine],bool))}
        dictionaries={};method_summaries={};group_records=[];crosswalk=[];method_labels={}
        for field,values in [('native_major',native_major),('low_reference_major',low_major),('discovery_route',route[sl])]:
            out[field+'_code'],dictionaries[field]=categories(['not_available' if v is None else v for v in values])
        old=reader.frozen_labels(ds);oldfine=old['fine'][positions].astype(str);oldmajor=old['major'][positions].astype(str)
        for method in ['A','B','C']:
            bycode={g['partition_code']:g for g in review['groups'] if g['method']==method}
            if set(bycode)!=set(np.unique(parts[method][sl])):raise ValueError('Source group table does not cover partition')
            records=[]
            for i,code in enumerate(parts[method][sl]):
                group=bycode[int(code)]
                records.append(decide_cell(native_sets[i],reference_sets[i],native_major[i],low_major[i],
                    group['local_major'],group['local_celltype'],own_sets[i],low_fine[i],positive[i],fine_fit[method][sl][i]))
            for field in ['consensus_major','consensus_celltype','coarse_evidence','fine_evidence','direct_fine_candidate']:
                out[method+'__'+field+'_code'],dictionaries[method+'__'+field]=categories([r[field] for r in records])
            for field in ['group_major_disagreement','group_type_disagreement','fine_discovery_fit_available']:
                out[method+'__'+field]=np.array([r[field] for r in records],bool)
            out[method+'__parent_discovery_group_code']=parts[method][sl]
            coarse=np.array([r['consensus_major'] for r in records]);fine=np.array([r['consensus_celltype'] for r in records]);method_labels[method]=fine
            # Across-method guard checks on the full actual source, not sampled cases.
            hard=(native.sum(axis=1)>1)|(reference.sum(axis=1)>1)|((native.sum(axis=1)==1)&(reference.sum(axis=1)==1)&(bits(native)!=bits(reference)))|~positive
            if (coarse[hard]!='unresolved').any():raise ValueError('Cell conflict overwritten by group')
            gamma=fine=='T_gamma_delta_RNA_candidate'
            if (gamma & ~calls['gamma_delta_T_related']).any():raise ValueError('Gamma identity lacks own RNA evidence')
            if not np.array_equal(out[method+'__fine_discovery_fit_available'],fine_fit[method][sl]):raise ValueError('Fine-fit exposure changed')
            method_summaries[method]={'n_cells':n,'major_counts':frequencies(coarse),'celltype_counts':frequencies(fine),
                'coarse_evidence_counts':frequencies(r['coarse_evidence'] for r in records),'fine_evidence_counts':frequencies(r['fine_evidence'] for r in records),
                'n_group_major_disagreements':int(out[method+'__group_major_disagreement'].sum()),'n_group_type_disagreements':int(out[method+'__group_type_disagreement'].sum()),
                'n_direct_fine_candidates_in_coarse_only_branches':sum(r['direct_fine_candidate']!='unresolved' and not r['fine_discovery_fit_available'] for r in records)}
            for (before,after),num in collections.Counter(zip(oldfine,fine)).items():crosswalk.append({'method':method,'old_fine':str(before),'guarded_celltype':str(after),'n_cells':num})
            for code,g in bycode.items():
                member=parts[method][sl]==code
                group_records.append({'method':method,'source_local_group_id':g['source_local_group_id'],'parent_group_sha256':g['group_sha256'],
                    'partition_code':code,'n_cells':int(member.sum()),'guarded_celltype_counts':frequencies(fine[member]),
                    'n_expression_units_not_patients':g['n_expression_units_not_patients'],'display_marker':g['display_marker'],'display_number':g['display_number'],
                    'marker_origin':'parent_discovery_group_not_recomputed_for_guarded_subset','marker_evidence_type':'parent_group_expression_specificity_not_DEG',
                    'original_own_group_marker_evidence_file':entry['output'],'subgroup_inference_not_new_clustering':True})
        for a,b in [('A','B'),('A','C'),('B','C')]:
            # Strong direct cell evidence is method-independent even in minority groups.
            def direct_mask(m):
                code=out[m+'__fine_evidence_code'];labels=np.array(dictionaries[m+'__fine_evidence'])
                return np.isin(labels[code],['own_native_and_unique_reference_supported','own_gamma_RNA_and_reference_candidate_requires_orthogonal_review'])
            da,db=direct_mask(a),direct_mask(b)
            np.testing.assert_array_equal(da,db);np.testing.assert_array_equal(method_labels[a][da],method_labels[b][db])
        index_field=z['obs'].attrs.get('_index','_index');ids=reader.field(ds,index_field)[positions].astype(str)
        if len(np.unique(ids))!=n:raise ValueError('Native cell keys not unique inside source')
        cell_key_sha=hashlib.sha256(('\0'.join(ds+'::'+v for v in ids)).encode()).hexdigest()
        side=root/'outputs'/(ds+'_cell_identity_confidence_v2.npz')
        with side.open('xb') as f:np.savez_compressed(f,**out)
        with np.load(side,allow_pickle=False) as f:
            for k,v in out.items():np.testing.assert_array_equal(f[k],v)
        payload={'source':ds,'n_pilot_cells':n,'cell_axis_sha256':cell_key_sha,'selector_sha256':inv['selector_sha256'],
                 'categorical_dictionaries':dictionaries,'native_family_bits':FAMILIES,'fine_type_bits':FINE_TYPES,'methods':method_summaries,
                 'parent_group_records':group_records,'old_fine_crosswalk':crosswalk,'native_author_field_provenance':sorted({g['author_field_origin'] for g in review['groups']}),
                 'label_version':cfg['label_version'],'candidate_labels_not_new_authority':True,'confidence_is_not_calibrated_probability':True,
                 'same_evidence_policy_for_all_methods':True,'new_reference_family_mask_changed_from_previous_routing_cells':int((bits(reference)!=previous_reference[sl]).sum()),
                 'RNA_state_and_malignancy_separate':True,'clinical_outcomes_used':False}
        desc=root/'outputs'/(ds+'_cell_identity_confidence_v2.json.gz')
        with gzip.open(desc,'xt') as f:json.dump(payload,f,allow_nan=False)
        outputs.append({'source':ds,'n_pilot_cells':n,'sidecar':side.name,'sidecar_sha256':sha256_file(side),'metadata':desc.name,'metadata_sha256':sha256_file(desc),'cell_axis_sha256':cell_key_sha})
        summaries.append({'source':ds,'n_pilot_cells':n,'methods':method_summaries,'reference_family_mask_changes':payload['new_reference_family_mask_changed_from_previous_routing_cells']})
        offset+=n;save(root/'checkpoints'/(ds+'.json'),{'source':ds,'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat()})
        if sum(p.stat().st_size for p in root.rglob('*') if p.is_file())>cfg['output_cap_bytes']:raise ValueError('Confidence sidecar cap exceeded')
    if offset!=len(source) or len(outputs)!=len(inventory):raise ValueError('Incomplete confidence source scope')
    save(root/'outputs/cell_identity_confidence_index_v2.json',outputs);save(root/'outputs/source_cell_identity_confidence_summary_v2.json',summaries)
    save(root/'provenance/terminal_v1.json',{'status':'COMMON_CELL_CONFIDENCE_V2_DROPOUT_CONFLICT_AND_CROSSWALK_COMPLETE_NOT_NEW_AUTHORITY',
        'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_pilot_cells':offset,'n_sources':len(outputs),'n_tests_actual_no_skips':10,
        'output_sha256':sha256_file(root/'outputs/cell_identity_confidence_index_v2.json'),'summary_sha256':sha256_file(root/'outputs/source_cell_identity_confidence_summary_v2.json'),
        'all_source_native_panel_group_aggregates_reproduced':True,'all_hard_cell_conflicts_retained':True,'all_direct_cell_fine_candidates_method_invariant':True,
        'all_sidecars_serialized_readback_identical':True,'wall_seconds':time.monotonic()-started,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
        'written_bytes_before_terminal':sum(p.stat().st_size for p in root.rglob('*') if p.is_file()),'new_authoritative_labels':False,'native_counts_copied':False})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)})
        raise
