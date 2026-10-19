"""Outcome-free source descriptions of fixed candidate capture versus RNA depth.

Reuse accepted count/depth metadata only. This does not classify new cells or
infer receptor identity. Individual arrays remain in server RAM; outputs are
source/type aggregates in task tmpfs. No model or clinical outcome is read.
"""
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def depth_bins(values, reference):
    edges = np.quantile(values[reference], [.25, .5, .75], method='linear')
    # Equal depths stay together; empty bins are retained, never jittered.
    return np.searchsorted(edges, values, side='right'), edges


def fixture():
    values = np.array([0, 1, 1, 1, 2, 3, 10, 10])
    bins, edges = depth_bins(values, np.ones(8, bool))
    assert len(edges) == 3 and np.all(bins[values == 1] == bins[1])
    constant, _ = depth_bins(np.ones(8), np.ones(8, bool))
    assert np.array_equal(np.bincount(constant, minlength=4), [0, 0, 0, 8])


def main():
    fixture()
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    out = Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert out.parent.parent == Path('/dev/shm') and not out.exists()
    source = root / 'NSCLC_full_native_gene_units_20261016_v1'
    hard = root / 'NSCLC_gamma_delta_RNA_candidates_20261019_v1'
    audit = root / 'NSCLC_gamma_delta_RNA_candidate_audit_20261019_v1/private_terminal.json'
    assert sha(audit) == 'e4d6c57c1727d33dd4d123178929125787706afa64964520c7c0ab1ff7429502'
    terminal = json.loads((source / 'private_terminal.json').read_text())
    hard_terminal = json.loads((hard / 'private_terminal.json').read_text())
    names = ['private_cell_unit_index.npy', 'private_full_RNA.npy',
             'private_detected_genes.npy', 'private_all_native_units_and_quality.tsv']
    hashes = {}
    for name in names:
        hashes[name] = sha(source / name)
        assert hashes[name] == terminal['outputs'][name]['sha256']
    sidecar = hard / 'private_candidate_identity_sidecar.npz'
    assert sha(sidecar) == hard_terminal['files'][sidecar.name]['sha256']
    out.mkdir(mode=0o700)
    save = lambda name, data: (out / name).write_text(json.dumps(data, indent=2) + '\n')
    save('plan_before_depth_descriptions.json', dict(
        created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
        input_SHA256=hashes, candidate_SHA256=sha(sidecar), source='GSE243013',
        axes=['native_total_RNA', 'native_detected_genes'],
        quantile_reference='All original author T/NK cells; no candidate/outcome selection',
        quantiles=[.25, .5, .75], quantile_method='linear', equal_edge_assignment='right',
        scopes=['all51_author_types', 'five_disjoint_author_roles'],
        eligible_author_record_minimum_cells_per_type_or_role_and_bin=20,
        zero_denominator='Missing proportion, never zero capture rate',
        rules='Four unchanged original rule bits plus their union',
        no_new_identity_classification=True, no_clinical_test=True,
        no_new_expression_or_sequence_read=True, aggregates_only=True))
    u = pd.read_csv(source / names[3], sep='\t', keep_default_na=False).sort_values('unit_index')
    assert u.unit_index.tolist() == list(range(11701))
    unit = np.load(source / names[0], allow_pickle=False)
    RNA = np.load(source / names[1], allow_pickle=False)
    detected = np.load(source / names[2], allow_pickle=False)
    h = dict(np.load(sidecar, allow_pickle=False))
    assert len(unit) == len(RNA) == len(detected) == 1254749
    assert np.all(np.diff(h['native_row']) > 0)
    bits = np.zeros(len(unit), np.uint8)
    bits[h['native_row']] = h['rule_bits']
    assert np.array_equal(unit[h['native_row']], h['unit_index'])
    assert not np.any(((bits & 4) > 0) & ((bits & 1) == 0))
    assert not np.any(((bits & 8) > 0) & ((bits & 2) == 0))
    labels = u.sub_cell_type.to_numpy(str)[unit]
    lymphoid = u.major_cell_type.eq('T/NK cell').to_numpy()[unit]
    roles = np.full(len(unit), 'Non_T_NK', dtype=object)
    roles[lymphoid] = 'Other_T_NK'
    roles[np.char.startswith(labels, 'NK_')] = 'Author_NK'
    roles[labels == 'T_gdT_TRDV1'] = 'Author_Vdelta1'
    roles[labels == 'T_gdT_TRDV2'] = 'Author_Vdelta2'
    record = pd.Categorical(u.sampleID).codes[unit]
    assert record.max() + 1 == 243
    rules = {f'original_rule_bit_{j}': (bits & (1 << j)) > 0 for j in range(4)}
    rules['fixed_RNA_candidate_union'] = bits != 0
    rows, bins_contract = [], []
    for axis, values in [('native_total_RNA', RNA), ('native_detected_genes', detected)]:
        assert np.isfinite(values).all() and (values >= 0).all()
        bins, edges = depth_bins(values, lymphoid)
        bins_contract.append(dict(source='GSE243013', axis=axis,
                                  reference_cells=int(lymphoid.sum()), boundaries=edges.tolist()))
        for scope, groups in [('all51_author_types', labels), ('five_disjoint_author_roles', roles)]:
            categories = sorted(set(groups))
            for group in categories:
                group_mask = groups == group
                # Four bins, including empty positions, form a complete partition.
                assert np.bincount(bins[group_mask], minlength=4).sum() == group_mask.sum()
                for depth_bin in range(4):
                    mask = group_mask & (bins == depth_bin)
                    n = int(mask.sum())
                    sizes = np.bincount(record[mask], minlength=243)
                    eligible = sizes >= 20
                    item = dict(source='GSE243013', scope=scope, author_group=group,
                                axis=axis, depth_bin=depth_bin + 1, cells=n,
                                author_records_with_cells=int((sizes > 0).sum()),
                                author_records_atleast20=int(eligible.sum()),
                                RNA_mean=float(RNA[mask].mean()) if n else None,
                                detected_genes_mean=float(detected[mask].mean()) if n else None,
                                patient_independence_not_established=True)
                    for name, flag in rules.items():
                        counts = np.bincount(record[mask & flag], minlength=243)
                        assert (counts <= sizes).all()
                        rates = counts[eligible] / sizes[eligible]
                        item[name + '_cells'] = int(counts.sum())
                        item[name + '_fraction'] = float(counts.sum() / n) if n else None
                        item[name + '_equal_record_median'] = float(np.median(rates)) if len(rates) else None
                    # Independent direct-mask counts check every aggregate numerator.
                    for name, flag in rules.items():
                        assert item[name + '_cells'] == np.count_nonzero(flag[mask])
                    rows.append(item)
    assert len(rows) == (51 + 5) * 2 * 4
    pd.DataFrame(rows).to_csv(out / 'private_all448_author_group_depth_descriptions.tsv', sep='\t', index=False)
    save('private_depth_bin_contract.json', bins_contract)
    # Source summaries are an exhaustive partition; no causal dropout inference.
    save('private_interpretation_boundary.json', dict(
        existing_candidate_flags_unchanged=True, author_labels_unchanged=True,
        counts_measure_rule_capture_not_receptor_identity=True,
        depth_association_not_proof_of_dropout_causality=True,
        NK_or_alpha_beta_RNA_not_automatic_gamma_exclusion=True,
        absent_author_labels_not_absent_gamma_cells=True,
        no_new_P_or_FDR=True, no_candidate_subset_chosen_using_clinical_effects=True))
    t = dict(status='COMPLETE_FIXED_CANDIDATE_DEPTH_DESCRIPTIONS_NO_IDENTITY_ADOPTION',
             completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
             aggregate_positions=len(rows), all_numerators_direct_mask_rechecked=True,
             private_individual_arrays_persisted=False,
             files={p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.iterdir()})
    assert sum(v['bytes'] for v in t['files'].values()) < 1024**2
    save('private_terminal.json', t)
    print(json.dumps(dict(status=t['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
