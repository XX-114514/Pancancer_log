"""All-source review axes and smallest-source direct counts/reference replay.

The replay checks computation, not biological annotation truth. Scientific
arrays stay in RAM and source-private outputs; no patient or cell rows exported.
"""
import argparse
import collections
import datetime
import gzip
import json
from pathlib import Path
import resource
import time

import anndata as ad
import celltypist
from celltypist.models import Model
import numpy as np
import pandas as pd
import scipy.sparse as sp

from annotation_representation_pilot_v1 import save
from full_context_atlas_v1 import FullContextAtlas, sha256_file
from native_annotation_review_v1 import LINEAGE_PANELS, TYPE_PANELS, STATE_PANELS, REVIEW_GENES


def direct_rows(group, requested):
    """Direct CSR spans with an independent pointer/offset implementation."""
    unique, inverse = np.unique(requested, return_inverse=True)
    pieces = [];i = 0
    while i < len(unique):
        j = i + 1
        while j < len(unique) and unique[j] - unique[i] < 12000:j += 1
        lo, hi = int(unique[i]), int(unique[j - 1]) + 1
        ptr = np.asarray(group['indptr'][lo:hi + 1], np.int64)
        first, last = int(ptr[0]), int(ptr[-1])
        block = sp.csr_matrix((np.asarray(group['data'][first:last]), np.asarray(group['indices'][first:last]), ptr - first),
                              shape=(hi - lo, int(group.attrs['shape'][1])))
        pieces.append(block[unique[i:j] - lo]);i = j
    return sp.vstack(pieces, format='csr')[inverse]


def main(root, producer):
    start = time.monotonic()
    terminal = json.loads((producer / 'provenance/terminal_v1.json').read_text())
    supervision = json.loads((producer / 'provenance/supervisor_terminal_v1.json').read_text())
    if supervision['actual_child_exit_code'] != 0 or supervision['wall_timeout']:raise ValueError('Producer is not complete')
    cfg = json.loads((producer / 'config/review_spec_v1.json').read_text())
    for item in cfg['bindings'].values():
        if sha256_file(item['path']) != item['sha256']:raise ValueError('Frozen input changed')
    index_path = producer / 'outputs/source_annotation_review_index_v1.json'
    if sha256_file(index_path) != terminal['index_sha256']:raise ValueError('Producer index changed')
    index = json.loads(index_path.read_text());pilot = Path(cfg['pilot']);discovery = Path(cfg['discovery'])
    inventory = json.loads((pilot / 'outputs/pilot_input_inventory_v1.json').read_text())
    with np.load(discovery / 'outputs/derived_axis_and_model_mask.npz', allow_pickle=False) as f:
        sources = f['source'];modeled = f['modeled']
    with np.load(discovery / ('outputs/candidate_partitions_' + str(cfg['seed']) + '.npz'), allow_pickle=False) as f:
        partitions = {m: f[m] for m in ['A', 'B', 'C']}
    reader = FullContextAtlas(cfg['manifest'])
    meta = json.loads((producer / 'outputs/reference_models_v1.json').read_text())
    checked = [];offset = 0;selected_source = min(inventory, key=lambda x: (x['n_selected_cells'], x['dataset']))['dataset']
    selected = None
    for inv, row in zip(inventory, index):
        ds = inv['dataset'];n = inv['n_selected_cells'];sl = slice(offset, offset + n)
        if row['source'] != ds or not np.all(sources[sl] == ds):raise ValueError('Source axis mismatch')
        p = producer / 'outputs' / row['output']
        if sha256_file(p) != row['output_sha256']:raise ValueError('Source review changed')
        with gzip.open(p, 'rt') as f:review = json.load(f)
        proposal = producer / 'outputs' / (ds + '_reference_proposals_v1.npz')
        if sha256_file(proposal) != row['reference_proposal_sha256']:raise ValueError('Reference sidecar changed')
        with np.load(proposal, allow_pickle=False) as f:side = {k: f[k] for k in f.files}
        if any(len(x) != n for x in side.values()) or not np.array_equal(side['positive_RNA'], modeled[sl]):raise ValueError('Sidecar/model axis differs')
        refs = ['immune_low', 'immune_high', cfg['source_organ_model'][ds]]
        for name in refs:
            unique, best, pos = [side[name + '__' + k] for k in ['unique_code', 'best_code', 'n_positive_labels']]
            prob = side[name + '__maximum_probability']
            if ((unique < -1) | (unique >= len(meta[name]['labels']))).any():raise ValueError('Reference code bounds')
            if not np.array_equal(unique >= 0, pos == 1):raise ValueError('Unique reference assignment differs')
            if not np.array_equal(np.isfinite(prob), pos >= 0):raise ValueError('Reference missingness differs')
            if (prob[np.isfinite(prob)] < 0).any() or (prob[np.isfinite(prob)] > 1).any():raise ValueError('Probability bounds')
            if name.startswith('immune_') and (pos[~side['immune_reference_scope']] != -1).any():raise ValueError('Out of scope immune prediction')
        records = review['groups']
        if len({r['group_sha256'] for r in records}) != len(records):raise ValueError('Duplicate stable group identity')
        positions = np.load(pilot / inv['selector_file'], allow_pickle=False)
        labels = reader.frozen_labels(ds);fine = labels['fine'][positions].astype(str);major = labels['major'][positions].astype(str)
        for method in ['A','B','C']:
            group_records = [r for r in records if r['method'] == method]
            current = partitions[method][sl]
            if sorted(r['partition_code'] for r in group_records) != list(np.unique(current)):raise ValueError('Candidate group coverage differs')
            for r in group_records:
                member = current == r['partition_code'];nrna = int((member & side['positive_RNA']).sum())
                if r['n_pilot_cells'] != int(member.sum()) or r['n_native_RNA_evaluable'] != nrna:raise ValueError('Group denominator differs')
                if dict(collections.Counter(fine[member])) != r['inherited_fine_crosswalk_not_truth']:raise ValueError('Old fine crosswalk differs')
                if dict(collections.Counter(major[member])) != r['inherited_major_counts_not_truth']:raise ValueError('Old major crosswalk differs')
                if r['malignant_status'] != 'not_inferred_from_lineage' or not r['RNA_state_separate']:raise ValueError('Identity/state/malignancy conflated')
                for name in refs:
                    ref = r['reference_proposals'][name]
                    labels_here = [meta[name]['labels'][c] for c in side[name + '__unique_code'][member] if c >= 0]
                    if dict(collections.Counter(labels_here)) != ref['unique_label_counts']:raise ValueError('Group reference count differs')
                    if sum(ref['unique_label_counts'].values()) + ref['unassigned'] + ref['multiple_labels'] + ref['not_evaluable_or_out_of_scope'] != int(member.sum()):
                        raise ValueError('Reference denominator differs')
                for value in r['native_multi_marker_fraction'].values():
                    if not 0 <= value <= 1:raise ValueError('Marker fraction bounds')
                if r['marker_candidates'] and r['display_marker'] != r['marker_candidates'][0]['gene']:raise ValueError('Display marker differs')
                if r['marker_origin'] != 'own_source_pilot_group' or r['valid_replicate_aware_DEG_available_in_this_run']:raise ValueError('Marker origin or DEG claim differs')
        checked.append({'source': ds, 'n_pilot_cells': n, 'all_groups_crosswalk_denominators_reference_codes_checked': True})
        if ds == selected_source:selected = (inv, side, review, sl, positions)
        offset += n
    if len(index) != len(inventory) or offset != len(sources):raise ValueError('Full pilot source scope differs')
    inv, side, review, sl, positions = selected;ds = inv['dataset'];n = len(positions)
    z = reader.store(ds);key = reader.legacy_config['counts_key'] if reader.is_legacy(ds) else reader.primary.partitions[ds]['counts_key']
    pieces = []
    for lo in range(0, n, 2048):pieces.append(direct_rows(z[key], reader.rows(ds)[positions[lo:lo + 2048]]))
    raw = sp.vstack(pieces, format='csr');del pieces
    totals = np.asarray(raw.sum(axis=1, dtype=np.float64)).ravel();positive = totals > 0
    if not np.array_equal(positive, side['positive_RNA']):raise ValueError('Direct native denominator mask differs')
    lookup = reader.feature_lookup(ds)
    names = np.array(['NATIVE_NONMODEL_COLUMN_' + str(i) for i in range(raw.shape[1])], dtype=object)
    # Independently assemble physical source gene names from var + lookup.
    from anndata._io.specs import read_elem
    var = read_elem(z['var'])
    for i in range(len(var)):
        gene = str(var.iloc[i]['gene_symbol']) if 'gene_symbol' in var else str(var.index[i])
        if lookup.get(gene, set()) == {i}:names[i] = gene
    scale = np.zeros(n);scale[positive] = 10000 / totals[positive]
    cp = raw.astype(np.float64).multiply(scale[:, None]).tocsr();log = cp.copy();log.data = np.log1p(log.data)
    refs = ['immune_low','immune_high',cfg['source_organ_model'][ds]];max_probability_error = 0.;predicted_values = 0
    for name in refs:
        if sha256_file(cfg['models'][name]['path']) != cfg['models'][name]['sha256']:raise ValueError('Cached reference model changed')
        model = Model.load(cfg['models'][name]['path'])
        scope = side['immune_reference_scope'] if name.startswith('immune_') else positive
        if review['reference_coverage'][name]['feature_fraction'] < cfg['minimum_model_feature_fraction']:continue
        rows = np.flatnonzero(scope)
        for lo in range(0, len(rows), 2048):
            chosen = rows[lo:lo + 2048]
            native = ad.AnnData(log[chosen], var=pd.DataFrame(index=names))
            official = celltypist.annotate(native, model=model, majority_voting=False, mode='prob match', p_thres=.5)
            probabilities = official.probability_matrix.to_numpy();npositive = (probabilities > .5).sum(axis=1)
            best = official.decision_matrix.to_numpy().argmax(axis=1);unique = np.where(npositive == 1, best, -1)
            for field, expected in [('best_code',best),('unique_code',unique),('n_positive_labels',npositive)]:
                np.testing.assert_array_equal(side[name + '__' + field][chosen], expected)
            err = float(np.max(np.abs(side[name + '__maximum_probability'][chosen] - probabilities.max(axis=1))))
            max_probability_error = max(max_probability_error, err)
            np.testing.assert_allclose(side[name + '__maximum_probability'][chosen], probabilities.max(axis=1), atol=3e-6, rtol=3e-6)
            predicted_values += len(chosen)
    measured = {g: len(lookup.get(g,set())) == 1 for g in REVIEW_GENES}
    det = {g: (np.asarray(raw[:, next(iter(lookup[g]))].toarray()).ravel() > 0) if measured[g] else np.zeros(n,bool) for g in REVIEW_GENES}
    joint = {p: np.sum([det[g] for g in genes],axis=0) >= 2 for p,genes in {**LINEAGE_PANELS, **TYPE_PANELS}.items()}
    joint['B'] &= det['CD79A'] | det['MS4A1']
    joint['Treg_related'] = det['FOXP3'] & (det['IL2RA'] | det['CTLA4'])
    joint['gamma_delta_T_related'] = det['TRDC'] & (det['TRGC1'] | det['TRGC2']) & det['CD3D']
    joint['CD8_T_related'] = det['CD8A'] & det['CD8B'] & (det['CD3D'] | det['CD3E'])
    checked_marker_values = 0;max_marker_error = 0.
    for r in review['groups']:
        member = (partitions[r['method']][sl] == r['partition_code']) & positive;nrna = int(member.sum())
        for p, values in joint.items():
            expected = float(values[member].mean()) if nrna else 0.
            if abs(expected-r['native_multi_marker_fraction'][p]) > 1e-12:raise ValueError('Direct joint-marker replay differs')
            checked_marker_values += 1
        for gene, observed in r['native_review_gene_detection'].items():
            expected = float(det[gene][member].mean())
            if abs(expected-observed) > 1e-12:raise ValueError('Direct gene detection replay differs')
            checked_marker_values += 1
        others = [x for x in review['groups'] if x['method']==r['method'] and x['partition_code']!=r['partition_code']]
        if r['marker_comparator'].startswith('other_source_pilot_groups_in_same_native_supported_parent'):
            others = [x for x in others if x['local_major']==r['local_major']]
        contrast = np.isin(partitions[r['method']][sl], [x['partition_code'] for x in others]) & positive
        for marker in r['marker_candidates']:
            col = next(iter(lookup[marker['gene']]))
            avg = float(cp[member,col].mean());other = float(cp[contrast,col].mean())
            detection = float((raw[member,col]>0).mean());otherdet = float((raw[contrast,col]>0).mean())
            ratio = float(np.log2((avg+.1)/(other+.1)))
            values = [(avg,marker['native_mean_CP10K']), (detection,marker['native_detection_fraction']),
                      (otherdet,marker['comparator_detection_fraction']), (ratio,marker['descriptive_log2_mean_ratio'])]
            for expected,observed in values:
                max_marker_error = max(max_marker_error,abs(expected-observed))
                if not np.isclose(expected,observed,atol=1e-8,rtol=1e-8):raise ValueError('Own-group expression marker replay differs')
                checked_marker_values += 1
    output = {'status':'ALL_REVIEW_AXES_AND_SIZE_SELECTED_DIRECT_NATIVE_REFERENCE_REPLAY_PASS_NOT_BIOLOGICAL_TRUTH',
              'source_checks':checked,'representative_source':ds,
              'representative_selection':'smallest_frozen_pilot_source_not_annotation_result',
              'n_direct_native_cells':n,'n_official_reference_cell_model_predictions_replayed':predicted_values,
              'max_reference_probability_error':max_probability_error,'n_native_marker_values_replayed':checked_marker_values,
              'max_native_marker_error':max_marker_error,
              'all_source_reference_and_native_numeric_replay':False,'independent_biological_identity_validation':False,
              'marker_ranking_replayed':False,'state_panel_means_numeric_replay':False}
    save(root/'outputs/review_audit_v1.json',output)
    save(root/'provenance/terminal_v1.json',{'status':output['status'],
         'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
         'output_sha256':sha256_file(root/'outputs/review_audit_v1.json'),
         'producer_terminal_sha256':sha256_file(producer/'provenance/terminal_v1.json'),
         'wall_seconds':time.monotonic()-start,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);p.add_argument('--producer',required=True)
    a=p.parse_args();root=Path(a.output)
    try:main(root,Path(a.producer))
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)})
        raise
