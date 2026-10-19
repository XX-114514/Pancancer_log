"""ARCHIVAL UNEXECUTED HARD-RULE EXTENSION; scientific gate not accepted.

Frozen gamma-delta RNA candidate sidecars on all qualified cached atlas arms.

Requires a separately reviewed calibration receipt before any expression read.
No original identity is modified; unmeasured genes or rules stay unavailable.
"""
import io
import json
import os
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
import zarr
from gamma_delta_RNA_rules_v1 import MARKERS, RULES, flags, digest_array, fixture
from NSCLC_external_lineage_inputs_v2 import load, cell_arms
from external_functional_coverage_inputs_v1 import read_inputs
from merge_public_preprocessed_zarr import csr_read
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def physical_columns(root, arm, dataset, axis):
    mapping = [axis['lookup'].get(g, set()) for g in MARKERS]
    indices = np.array([next(iter(v)) if len(v) == 1 else -1 for v in mapping], np.int64)
    known = indices >= 0
    assert len(np.unique(indices[known])) == known.sum()
    if arm == 'frozen_original_annotation':
        p = root.parents[2] / 'Project_v3/annotation_rebuild/20261010_uniform_expression_annotation_v1/by_gse' / dataset / 'genes.tsv'
        genes = pd.read_csv(p, sep='\t', keep_default_na=False)
        indices[known] = genes.core_column.to_numpy(np.int64)[indices[known]]
    return indices, known


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT']); guard = resource(root, 10 * 1024**2)
    review = root / 'NSCLC_gamma_delta_RNA_candidate_review_20261019_v1/private_candidate_extension_decision.json'
    decision = json.loads(review.read_text())
    assert decision['continue_fixed_candidate_only_extension'] is True
    audit = root / 'NSCLC_gamma_delta_RNA_candidate_audit_20261019_v1/private_terminal.json'
    assert sha(audit) == decision['calibration_audit_sha256']
    assert json.loads(audit.read_text())['status'] == 'ACCEPTED_COMPLETE_NATIVE_CANDIDATE_SIDECAR_AUTHOR_AGGREGATES_AND_INDEPENDENT_FIRST8192_RAW_CELLS'
    out = root / 'all61_gamma_delta_RNA_candidates_20261019_v1'; assert not out.exists(); out.mkdir(mode=0o700)
    save = lambda name, value: (out / name).write_text(json.dumps(value, indent=2) + '\n')
    data = load(); provider = read_inputs(root); axes = provider['source_axes']; assert len(axes) == 61
    lib = data['lib']; assert len(lib) == 1690
    save('plan_before_all61_native_candidate_RNA.json', dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), rule_module_sha256=sha(Path(__file__).with_name('gamma_delta_RNA_rules_v1.py')), calibration_review_sha256=sha(review), calibration_audit_sha256=sha(audit), input_metadata_hashes=data['hashes'], native_axis_input_hashes=provider['input_sha256'], markers=MARKERS, rules=RULES, resource_guard=guard, production_output_cap_bytes=8 * 1024**2, audit_output_cap_bytes=2 * 1024**2, scope='All61 existing atlas-source arms and1690 libraries. All canonical selected cells in every uniquely core5-measurable source, including frozen-parent discordance; no new raw/quarantined/reference entries. All unavailable source/rule positions retained as missing, never biological negative. No clinical outcome or fitted threshold used.', mapping='Original frozen native source column -> existing core column; immune accepted exact/stable alias unique physical column. No new aliases or outcome-based source selection.', identity='Exploratory RNA candidates only. Frozen T, NK/boundary/lymphoid and other parent groups reported separately. No formal TCR/state/identity adoption, original node/CM untouched.', raw_audit='Per source metadata-seeded64 cells plus first32 positive candidates, de-duplicated. Independent manual CSR row scan checks all nine counts and rule bits only on these bounded rows; full candidate counts/keys checked separately. Not full independent raw replay.', no_new_expression_download=True))
    save('private_rule_fixture.json', fixture())
    counts = np.zeros((len(lib), 4), np.int64)
    v1 = counts.copy(); v2 = counts.copy()
    measured = np.zeros((len(lib), 4), bool)
    parent_counts = np.zeros((len(lib), 4, 3), np.int64)
    source_rows = []; proofs = []; audit_rows = []
    saved_bytes = 0; source_number = 0
    for chunk in cell_arms(data):
        arm = chunk['arm']; original = arm == 'frozen_original_annotation'
        z = zarr.open_group(str(data['old_zarr'] if original else data['atlas']), mode='r')
        raw = z['X'] if original else z['layers/counts']
        if original:
            with np.load(data['old'] / 'private_canonical_row_selection.npz', allow_pickle=False) as selected: global_rows = selected['global_rows']
        else: global_rows = np.arange(chunk['cell_count'], dtype=np.int64)
        assert len(global_rows) == chunk['cell_count'] and len(np.unique(global_rows)) == len(global_rows)
        majors = np.asarray(chunk['major'])[chunk['node'].codes]
        parent = np.full(len(majors), 2, np.int8)
        parent[np.isin(majors, ['T cell', 'T'])] = 0
        parent[np.isin(majors, ['T/NK lymphoid', 'Lymphoid cell', 'T_NK_boundary', 'NK'])] = 1
        parts = []
        for dataset in sorted(str(v) for v in chunk['source'].categories):
            positions = np.flatnonzero(chunk['source'].codes == chunk['source'].categories.get_loc(dataset))
            library_positions = lib.index[lib.arm.eq(arm) & lib.dataset.eq(dataset)].to_numpy(int)
            assert len(positions) == int(lib.loc[library_positions, 'annotated_scope_cells'].sum())
            np.testing.assert_array_equal(np.bincount(chunk['library_index'][positions], minlength=len(lib))[library_positions], lib.loc[library_positions, 'annotated_scope_cells'])
            physical, known = physical_columns(root, arm, dataset, axes[(arm, dataset)])
            gates = np.array([known[:5].all(), known[:7].all(), known[:5].all(), known[:7].all()])
            measured[library_positions] = gates
            source_number += 1
            source = dict(source_index=source_number, arm=arm, dataset=dataset, canonical_cells=len(positions), libraries=len(library_positions), uniquely_measured_markers=int(known.sum()), physical_columns=physical.tolist(), core5_measurable=bool(gates[0]), supported7_measurable=bool(gates[1]), RNA_read=bool(gates[0]), original_identity_mutated=False)
            if not gates[0]:
                source['status'] = 'RNA_NOT_READ_CORE_GENE_UNMEASURABLE'; source_rows.append(source); continue
            rows = global_rows[positions]; ordering = np.argsort(rows); rows = rows[ordering]; local_positions = positions[ordering]
            assert np.array_equal(rows, np.unique(rows))
            x = np.zeros((len(rows), len(MARKERS)), np.int32)
            for block in np.unique(rows // 4096):
                lo = int(block) * 4096; hi = min(lo + 4096, int(raw.attrs['shape'][0]))
                take = np.flatnonzero((rows >= lo) & (rows < hi)); current = csr_read(raw, lo, hi)
                cell = current[rows[take] - lo]
                marker = cell[:, physical[known]].toarray()
                assert np.isfinite(marker).all() and (marker >= 0).all() and marker.max(initial=0) < 2**31 and np.array_equal(marker, np.floor(marker))
                x[np.ix_(take, np.flatnonzero(known))] = marker.astype(np.int32)
            bits = flags(x, known); positive = np.flatnonzero(bits)
            li = chunk['library_index'][local_positions]; parents = parent[local_positions]
            for j in range(4):
                found = (bits & (1 << j)) > 0
                counts[:, j] += np.bincount(li[found], minlength=len(lib))
                if known[7]: v1[:, j] += np.bincount(li[found & (x[:, 7] > 0)], minlength=len(lib))
                if known[8]: v2[:, j] += np.bincount(li[found & (x[:, 8] > 0)], minlength=len(lib))
                for p in range(3): parent_counts[:, j, p] += np.bincount(li[found & (parents == p)], minlength=len(lib))
            parts.append(dict(global_row=rows[positive].astype(np.int64), library_position=li[positive].astype(np.int16), frozen_node_code=chunk['node'].codes[local_positions[positive]].astype(np.int16), parent_group=parents[positive], rule_bits=bits[positive], TRDV1_detection=np.where(known[7], x[positive, 7] > 0, -1).astype(np.int8), TRDV2_detection=np.where(known[8], x[positive, 8] > 0, -1).astype(np.int8)))
            rng = np.random.default_rng(2026101900 + source_number)
            chosen = np.unique(np.concatenate([rng.choice(len(rows), min(64, len(rows)), replace=False), positive[:32]]))
            audit_rows.extend(dict(source_index=source_number, arm=arm, dataset=dataset, global_row=int(rows[i]), library_position=int(li[i])) for i in chosen)
            proofs.append(dict(source_index=source_number, arm=arm, dataset=dataset, measured_markers=known.tolist(), physical_columns=physical.tolist(), raw_audit_rows_sha256=digest_array(rows[chosen], '<i8'), raw_audit_count_sha256=digest_array(x[chosen], '<i4'), raw_audit_rule_bits_sha256=digest_array(bits[chosen], 'u1'), raw_audit_cells=len(chosen), all_candidate_global_rows_sha256=digest_array(rows[positive], '<i8'), all_candidate_rule_bits_sha256=digest_array(bits[positive], 'u1')))
            source['status'] = 'FIXED_RNA_CANDIDATES_MEASURED_NO_FORMAL_IDENTITY'; source_rows.append(source)
            save('private_progress.json', dict(stage='ALL61_EXISTING_RNA_CANDIDATE_MEASUREMENT', completed_source_entries=source_number, total_source_entries=61, observed_UTC=datetime.now(timezone.utc).isoformat()))
        if parts:
            arrays = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
            payload = io.BytesIO(); np.savez_compressed(payload, **arrays); content = payload.getvalue()
            saved_bytes += len(content); assert saved_bytes < 6 * 1024**2
            (out / ('private_' + arm + '_candidate_identity_sidecar.npz')).write_bytes(content)
        save('private_' + arm + '_frozen_node_categories.json', [str(v) for v in chunk['node'].categories])
    assert source_number == 61
    tables = []
    for j, rule in enumerate(RULES):
        frame = lib[['arm', 'dataset', 'library_id', 'annotated_scope_cells']].copy()
        frame['rule'] = rule; frame['RNA_rule_measurable'] = measured[:, j]
        for name, values in [('candidate_cells', counts[:, j]), ('frozen_T_parent_candidate_cells', parent_counts[:, j, 0]), ('frozen_boundary_or_NK_candidate_cells', parent_counts[:, j, 1]), ('frozen_other_parent_candidate_cells', parent_counts[:, j, 2])]:
            frame[name] = pd.array(values, dtype='Int64'); frame.loc[~measured[:, j], name] = pd.NA
        tables.append(frame)
    all_rows = pd.concat(tables, ignore_index=True)
    all_rows.to_csv(out / 'private_all1690_library_four_rule_counts.tsv.gz', sep='\t', index=False, compression='gzip')
    summary = []
    for (arm, dataset, rule), g in all_rows.groupby(['arm', 'dataset', 'rule'], sort=True):
        good = g.RNA_rule_measurable
        summary.append(dict(arm=arm, dataset=dataset, rule=rule, libraries=len(g), measurable_libraries=int(good.sum()), libraries_candidate_atleast20=int(g.loc[good, 'candidate_cells'].ge(20).sum()), candidate_cells=int(g.loc[good, 'candidate_cells'].sum()) if good.any() else None, frozen_T_parent_candidate_cells=int(g.loc[good, 'frozen_T_parent_candidate_cells'].sum()) if good.any() else None, frozen_boundary_or_NK_candidate_cells=int(g.loc[good, 'frozen_boundary_or_NK_candidate_cells'].sum()) if good.any() else None, frozen_other_parent_candidate_cells=int(g.loc[good, 'frozen_other_parent_candidate_cells'].sum()) if good.any() else None, clinical_or_formal_identity_adoption=False))
    pd.DataFrame(summary).to_csv(out / 'private_complete61_source_rule_summary.tsv', sep='\t', index=False)
    pd.DataFrame(source_rows).to_csv(out / 'private_complete61_native_measurement_scope.tsv', sep='\t', index=False)
    pd.DataFrame(audit_rows).to_csv(out / 'private_prespecified_source_raw_audit_rows.tsv.gz', sep='\t', index=False, compression='gzip')
    save('private_all_source_count_and_candidate_proofs.json', proofs)
    files = {p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()}
    assert sum(d['bytes'] for d in files.values()) < 8 * 1024**2
    receipt = dict(status='COMPLETE_ALL61_FIXED_GAMMA_DELTA_RNA_CANDIDATES_BOUNDED_RAW_AUDIT_PENDING', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), files=files, source_arms=61, library_arm_entries=1690, full_marker_arrays_persisted=False, clinical_or_formal_identity_adoption=False, scientific_outputs_public=False)
    save('private_terminal.json', receipt)
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
