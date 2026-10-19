"""Metadata-only full-source reference availability; no RNA or outcome fitting."""
import hashlib
import json
import os
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_external_lineage_inputs_v2 import load, cell_arms
from NSCLC_native_MFA_interpretation_inputs_v1 import sha

KINDS = ['EXPLICIT_GAMMA_DELTA_AUTHOR', 'DELTA_SUBTYPE_TEXT_ONLY',
         'GAMMA_MARKER_TEXT_ONLY', 'EXPLICIT_ALPHA_BETA_AUTHOR',
         'EXPLICIT_NK_AUTHOR', 'OTHER_T_AUTHOR', 'UNRESOLVED_OR_OTHER',
         'CONFLICTING_GAMMA_ALPHA_BETA_LABEL']


def label_kind(value):
    s = str(value).casefold()
    gamma = bool(re.search(r'γ\s*δ|gamma[\s_\-/]*delta|(?:^|[^a-z0-9])(?:gd[\s_-]*t|t[\s_-]*gd|gdt)(?:$|[^a-z0-9])', s))
    ab = bool(re.search(r'α\s*β|alpha[\s_\-/]*beta|(?:^|[^a-z0-9])(?:abt|tab)(?:$|[^a-z0-9])', s))
    if gamma and ab:
        return 7
    if gamma:
        return 0
    if re.search(r'v\s*δ\s*[12]|v[\s_-]*(?:delta|d)[\s_-]*[12](?:$|[^0-9])', s):
        return 1
    if re.search(r'(?:^|[^a-z0-9])tr(?:dc|gc[12]|dv[123]|gv9)(?:$|[^a-z0-9])', s):
        return 2
    if ab:
        return 3
    if re.search(r'(?:^|[^a-z0-9])nk(?:$|[^a-z0-9])|natural[\s_-]*killer', s):
        return 4
    if re.search(r'(?:^|[^a-z0-9])t(?:[\s_-]*cells?|[\s_-]*lymph)|(?:^|[^a-z0-9])(?:cd[48]|treg)(?:$|[^a-z0-9])', s):
        return 5
    return 6


def fixture():
    cases = {'T_gdT_TRDV1': 0, 'gdT': 0, 'gamma-delta T cells': 0,
             'TCRγδ T': 0, 'Gamma_delta_T': 0, 'Tgd': 0,
             'Vdelta1': 1, 'Vδ2 T': 1, 'TRDC+ T': 2, 'TRDV1': 2,
             'alpha-beta T': 3, 'αβ T': 3, 'NK': 4, 'T_NK_boundary': 4,
             'CD8 T': 5, 'T cells': 5, 'Fibroblast': 6, 'unknown': 6,
             'gdT/alpha-beta mixture': 7, 'CD74 myeloid': 6}
    for value, expected in cases.items():
        assert label_kind(value) == expected, 'SEMANTIC_FIXTURE_FAILED'
    return dict(status='PASS_EXPLICIT_IDENTITY_MARKER_ONLY_MIXED_AND_UNRESOLVED_LABELS', cases=len(cases))


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    out = Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert out.parent.parent == Path('/dev/shm') and not out.exists()
    assert shutil.disk_usage('/dev/shm').free > 512 * 1024**2
    out.mkdir(mode=0o700)
    save = lambda name, value: (out / name).write_text(json.dumps(value, indent=2) + '\n')
    save('plan_before_all_source_author_reference_counts.json', dict(
        created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
        all_source_arm_positions=61, all_library_arm_entries=1690, label_kinds=KINDS,
        rules_source_sha256=sha(__file__), calibration_thresholds_selected=False,
        scope='Exact inherited author strings on frozen canonical rows; original52 plus archived immune9. Missing explicit labels are not biological absence. No expression, clinical endpoint, new classifier or formal annotation write.',
        source_role='Source-arm entries and duplicated accessions do not establish independent studies or patients.',
        storage='Project data read-only, task-owned tmpfs, small private local copy; no new raw downloads.', output_cap_bytes=2 * 1024**2))
    save('private_semantic_fixture.json', fixture())
    data = load()
    lib = data['lib']
    counts = np.zeros((len(lib), len(KINDS)), np.int64)
    label_rows, source_rows, proof = [], [], []
    for chunk in cell_arms(data):
        authors = chunk['author']
        kinds = np.asarray([label_kind(v) for v in authors.categories], np.int8)
        code = np.full(len(authors), 6, np.int8)
        valid = authors.codes >= 0
        code[valid] = kinds[authors.codes[valid]]
        li = chunk['library_index']
        keys, totals = np.unique(li * len(KINDS) + code, return_counts=True)
        counts[keys // len(KINDS), keys % len(KINDS)] += totals
        # Separate dataframe groupby checks every nonzero source-library/kind count.
        alternate = pd.DataFrame({'library': li, 'kind': code}).groupby(['library', 'kind'], sort=True).size()
        assert len(alternate) == len(keys)
        for (i, k), total in alternate.items():
            assert counts[i, k] == total
        source = chunk['source']
        label_codes = authors.codes.astype(np.int64) + 1
        packed = source.codes.astype(np.int64) * (len(authors.categories) + 1) + label_codes
        skeys, stotals = np.unique(packed, return_counts=True)
        for key, total in zip(skeys, stotals):
            si, ai = divmod(int(key), len(authors.categories) + 1)
            label = str(authors.categories[ai - 1]) if ai else 'AUTHOR_LABEL_MISSING'
            label_rows.append(dict(arm=chunk['arm'], dataset=str(source.categories[si]), author_label=label,
                                   semantic_kind=KINDS[label_kind(label)], n_cells=int(total),
                                   author_label_is_TCR_truth=False))
        parent = np.asarray(chunk['major'])[chunk['node'].codes]
        for dataset in source.categories:
            at = source.codes == source.categories.get_loc(dataset)
            ids = np.flatnonzero(lib.arm.eq(chunk['arm']) & lib.dataset.eq(str(dataset)))
            expected = int(lib.loc[ids, 'annotated_scope_cells'].sum())
            assert int(at.sum()) == expected == int(counts[ids].sum())
            explicit = at & (code == 0)
            lymphoid = np.isin(parent, ['T cell', 'T', 'T/NK lymphoid', 'Lymphoid cell', 'T_NK_boundary', 'NK'])
            row = dict(arm=chunk['arm'], dataset=str(dataset), canonical_cells=expected,
                       libraries=len(ids), libraries_with_explicit_gamma=int((counts[ids, 0] > 0).sum()),
                       explicit_gamma_in_frozen_lymphoid=int((explicit & lymphoid).sum()),
                       explicit_gamma_outside_frozen_lymphoid=int((explicit & ~lymphoid).sum()),
                       author_labels_not_receptor_gold_standard=True)
            row.update({name: int(counts[ids, j].sum()) for j, name in enumerate(KINDS)})
            source_rows.append(row)
        proof.append(dict(arm=chunk['arm'], canonical_cells=chunk['cell_count'],
                          inherited_author_field=chunk['author_field'], author_values_sha256=chunk['author_values_sha256'],
                          all_nonzero_library_kind_counts_independent_groupby_checked=True))
        save('private_progress.json', dict(stage='ALL_SOURCE_ARM_METADATA_COUNT', completed_arms=len(proof), observed_UTC=datetime.now(timezone.utc).isoformat()))
    assert len(source_rows) == 61 and np.array_equal(counts.sum(1), lib.annotated_scope_cells.to_numpy())
    pd.DataFrame(source_rows).to_csv(out / 'private_complete61_source_author_reference_scope.tsv', sep='\t', index=False)
    pd.DataFrame(label_rows).to_csv(out / 'private_complete_source_author_label_semantics.tsv.gz', sep='\t', index=False, compression='gzip')
    # No library or patient identifiers persist: only source/context coverage.
    work = lib[['arm', 'dataset', 'cancer_label', 'material_class', 'phase_class', 'capture_scope', 'ex_vivo', 'reviewed_patient']].copy()
    for j, name in enumerate(KINDS):
        work[name] = counts[:, j]
    work['explicit_gamma_positive_libraries'] = counts[:, 0] > 0
    work['library_entries'] = 1
    fields = ['arm', 'dataset', 'cancer_label', 'material_class', 'phase_class', 'capture_scope', 'ex_vivo']
    context = work.groupby(fields, sort=True, dropna=False)[KINDS + ['library_entries', 'explicit_gamma_positive_libraries', 'reviewed_patient']].sum().reset_index()
    assert context.library_entries.sum() == 1690
    context.to_csv(out / 'private_complete_source_context_author_reference_scope.tsv.gz', sep='\t', index=False, compression='gzip')
    save('private_all_author_metadata_proofs.json', dict(input_hashes=data['hashes'], arms=proof, counts_total_verified=True,
         inference='Author string semantics only; no RNA or TCR/sequence/protein inference and no evaluated identity accuracy.'))
    t = dict(status='COMPLETE_ALL61_AUTHOR_REFERENCE_AVAILABILITY_NO_NEW_RNA_OR_CLINICAL_TEST',
             completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
             all_source_arm_positions=61, library_arm_entries=1690, all_library_kind_counts_alternate_groupby_checked=True,
             full_expression_or_cell_score_arrays_persisted=False, original_identity_or_inference_changed=False,
             scientific_outputs_public=False, files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()})
    assert sum(v['bytes'] for v in t['files'].values()) < 2 * 1024**2
    save('private_terminal.json', t)
    print(json.dumps(dict(status=t['status'], terminal_sha256=sha(out/'private_terminal.json'))))


if __name__ == '__main__':
    main()
