"""All reviewed pretreatment material, fixed interface; RAM-only RNA and inference."""
import hashlib
import json
import os
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
from scipy.sparse import csr_matrix
from NSCLC_external_lineage_inputs_v2 import load
from fixed_interface_rank_family_math_v1 import infer, close
from qualify_fixed_interface_all_cached_axes_v2 import stable, read_native_var


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def array_sha(a):
    a = np.ascontiguousarray(a)
    return hashlib.sha256(str(a.shape).encode() + str(a.dtype).encode() + a.tobytes()).hexdigest()


def main():
    config = json.load(sys.stdin)
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT']); out = Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert out.parent.parent == Path('/dev/shm') and not out.exists()
    out.mkdir(mode=0o700)
    def save(name, value):
        (out / name).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    data = load()
    aq = Path(os.environ['PANCANCER_AXIS_QUALIFICATION'])
    assert sha(aq / 'private_terminal.json') == config['accepted_axis_terminal_SHA256']
    at = json.loads((aq / 'private_terminal.json').read_text())
    for name in ['private_all61_current_source_axis_qualification.tsv', 'private_metadata_input_and_audit_proof.json']:
        assert sha(aq / name) == at['files'][name]['sha256']
    axes = pd.read_csv(aq / 'private_all61_current_source_axis_qualification.tsv', sep='\t')
    chosen = axes[axes.arm.eq('archived_immune_atlas') & axes.dataset.eq('GSE286827')]
    assert len(chosen) == 1 and bool(chosen.all_fixed_roles_unique.iloc[0])
    metadata = root / 'HNSCC_author_clinical_table_user_20261014_v1/private_key_and_response_audit_v1/all57_private_author_table_join.tsv'
    clinical = root / 'HNSCC_frozen_MCP_clinical_preparation_20261014_v1/private_patient_clinical_definition.tsv'
    assert sha(metadata) == '9942706870518c20c87a45fed44caff6bec22a1a512b6f194dc17d406df5cd22'
    assert sha(clinical) == 'ce3c924ead635e51d322669e9b3504b5abee38a30044394eb9c4f7778f949db6'
    author = pd.read_csv(metadata, sep='\t', dtype=str, keep_default_na=False).set_index('library_id')
    c = pd.read_csv(clinical, sep='\t', dtype=str, keep_default_na=False).set_index('patient_id')
    assert len(author) == 57 and author.index.is_unique and len(c) == 29 and c.index.is_unique
    save('private_plan_before_current_gene_expression_and_effects.json', dict(created_UTC=datetime.now(timezone.utc).isoformat(),
         source_sha256=sha(__file__), math_sha256=sha(Path(__file__).with_name('fixed_interface_rank_family_math_v1.py')),
         axis_qualification_terminal_sha256=sha(aq / 'private_terminal.json'), author_bridge_sha256=sha(metadata),
         clinical_definition_sha256=sha(clinical), fixed_private_interface=config, primary_questions=2,
         selection='Every current reviewed author pretreatment material; one per patient; both frozen T/Myeloid >=20, original continuous pathology/HPV/treatment complete. Not restricted to old26 MCP subset.',
         full_controls='Treatment, HPV, both meanlogRNA, meanloggenes, logcells, captured-cell fractions',
         finite_nulls='3 scenarios x256 x999, reject design if any nominal.05 rate>.10; no retries',
         inference='19999 HC3 restricted-null wild,2000 patient bootstrap treatment xHPV,>=1900 validCI,allLOO,BH2/BY2 and simultaneous.99 MCupperBH2',
         limits='Exploratory source complement, prior project use; broad current parents not original fine gamma endpoint, no new formal identity or causal/protein claim',
         arrays_RAM_only=True, original_followup_CSV_not_processed=True, output_cap_bytes=4 * 1024 ** 2))
    lib = data['lib'][data['lib'].arm.eq('archived_immune_atlas') & data['lib'].dataset.eq('GSE286827')].copy()
    assert lib.library_id.is_unique and set(lib.library_id) <= set(author.index)
    records = []; reasons = Counter()
    for row in lib.sort_values('library_id').itertuples(index=False):
        a = author.loc[row.library_id]; p = str(a.patient_registry)
        assert row.patient_block == 'GSE286827::author_patient::' + p
        assert row.timepoint == a.timepoint_registry
        cl = c.loc[p]
        source_y = float(pd.to_numeric(a['Tumor.regression'], errors='coerce'))
        clinical_y = float(pd.to_numeric(cl.regression, errors='coerce'))
        assert ((not np.isfinite(source_y) and not np.isfinite(clinical_y)) or
                (np.isfinite(source_y) and source_y == clinical_y))
        assert a.neo_type == cl.treatment and a['HPV/p16'] == cl.HPV
        if a.timepoint_registry != 'pre':
            reasons['NOT_ORIGINAL_PRETREATMENT'] += 1; continue
        if not row.reviewed_patient:
            reasons['PATIENT_NOT_REVIEWED'] += 1; continue
        y = clinical_y
        if not np.isfinite(y) or not 0 <= y <= 100:
            reasons['MISSING_OR_INVALID_CONTINUOUS_PATHOLOGY'] += 1; continue
        if cl.treatment not in ['D', 'D+T'] or cl.HPV not in ['negative', 'positive']:
            reasons['MISSING_DECLARED_CLINICAL_CONTROL'] += 1; continue
        records.append(dict(library_id=row.library_id, patient=p, regression=y,
                            treatment=float(cl.treatment == 'D+T'), HPV=float(cl.HPV == 'positive'),
                            binary_concordant=str(cl.author_binary_concordant).lower() == 'true',
                            total_cells=int(row.annotated_scope_cells)))
    baseline = pd.DataFrame(records)
    assert baseline.patient.is_unique
    z = zarr.open_group(str(data['atlas']), mode='r')
    catalog = data['catalogs']['archived_immune_atlas'].set_index('node_id')
    datasets = np.asarray(read_elem(z['obs/dataset']))
    source_rows = np.flatnonzero(datasets == 'GSE286827')
    libraries = np.asarray(read_elem(z['obs/library_id']))[source_rows].astype(str)
    nodes = np.asarray(read_elem(z['obs/node_id']))[source_rows].astype(str)
    assert set(nodes) <= set(catalog.index)
    actual = Counter(zip(libraries, nodes))
    prior = data['counts'][data['counts'].arm.eq('archived_immune_atlas') & data['counts'].dataset.eq('GSE286827')]
    expected = {(str(r.library_id), str(r.node_id)): int(r.n_cells) for r in prior.itertuples(index=False)}
    assert dict(actual) == expected
    major = catalog.loc[nodes, 'major'].to_numpy(str)
    frozen = catalog.loc[nodes, 'frozen_role'].to_numpy(int)
    pool = np.where(major == 'T', 0, np.where(major == 'Myeloid', 1, -1))
    for row in baseline.itertuples(index=False):
        take = libraries == row.library_id
        assert int(take.sum()) == row.total_cells
    support = {(str(row.library_id), kind): int(((libraries == row.library_id) & (pool == kind)).sum())
               for row in baseline.itertuples(index=False) for kind in [0, 1]}
    eligible = baseline.library_id.map(lambda key: min(support[(key, 0)], support[(key, 1)]) >= 20)
    reasons['FEWER_THAN20_EITHER_CURRENT_PARENT'] += int((~eligible).sum())
    frame = baseline[eligible].sort_values('patient').reset_index(drop=True)
    assert frame.patient.is_unique and len(frame) >= 20
    var = read_native_var(z['var'])
    order = np.asarray(read_elem(z['uns/source_order'])).astype(str).tolist()
    mask = np.asarray(read_elem(z['varm/measured_in_source']), bool)
    proof = json.loads((aq / 'private_metadata_input_and_audit_proof.json').read_text())
    assert hashlib.sha256(var.to_json(orient='split').encode()).hexdigest() == proof['full_atlas_gene_metadata_sha256']
    assert hashlib.sha256(mask.tobytes()).hexdigest() == proof['full_source_mask_sha256']
    measured = mask[:, order.index('GSE286827')]
    columns = []
    for symbol, identifier in zip(config['symbols'], config['stable_ids']):
        matches = []
        for i in np.flatnonzero(measured):
            v = var.iloc[i]
            if str(v.feature_identity_status) == 'ambiguous_source_feature' or str(v.canonical_gene_identity).startswith('AMBIGUOUS::'):
                continue
            values = {str(var.index[i]), str(v.gene_symbol), str(v.canonical_gene_identity)}
            if str(v.feature_identity_status) == 'stable_Ensembl_ID':
                values.update(str(v.source_gene_aliases).split('|'))
            if symbol in values or identifier in {stable(x) for x in values}:
                matches.append(int(i))
        assert len(matches) == 1
        columns.append(matches[0])
    assert len(set(columns)) == 3
    library_index = {r.library_id: i for i, r in enumerate(frame.itertuples(index=False))}
    group = np.array([2 * library_index[k] + p if k in library_index and p >= 0 else -1 for k, p in zip(libraries, pool)])
    take = group >= 0; rows = source_rows[take]; groups = group[take]; roles = frozen[take]
    order_idx = np.argsort(rows); rows = rows[order_idx]; groups = groups[order_idx]; roles = roles[order_idx]
    ng = 2 * len(frame); nf = len(var); counts = z['layers/counts']
    assert int(counts.attrs['shape'][1]) == nf
    assert int(counts.attrs['shape'][0]) == len(datasets)
    sums = np.zeros((ng, 9)); alternate_sums = np.zeros((ng, 9))
    totals = np.bincount(groups, minlength=ng)
    for i, row in enumerate(frame.itertuples(index=False)):
        for kind in [0, 1]:
            assert totals[2 * i + kind] == support[(row.library_id, kind)] >= 20
    manual = {}; fine = np.zeros((len(frame), 3))
    for i in range(len(frame)):
        t = groups == 2 * i; my = groups == 2 * i + 1
        fine[i] = [np.count_nonzero(t & (roles == 0)) / t.sum(),
                   np.count_nonzero(t & (roles == 1)) / t.sum(), np.count_nonzero(my & (roles == 2)) / my.sum()]
        for kind in [0, 1]:
            for physical in rows[groups == 2 * i + kind][:2]:
                a, b = np.asarray(counts['indptr'][int(physical):int(physical) + 2], int)
                ix = np.asarray(counts['indices'][a:b], int); values = np.asarray(counts['data'][a:b], float)
                assert len(np.unique(ix)) == len(ix)
                manual[int(physical)] = (float(values.sum()), int(np.count_nonzero(values)),
                                         np.array([values[ix == col].sum() for col in columns]))
    seen = 0; checked_manual = set()
    for number, chunk in enumerate(np.unique(rows // 1024)):
        begin = int(chunk) * 1024; end = min(begin + 1024, int(counts.attrs['shape'][0]))
        ptr = np.asarray(counts['indptr'][begin:end + 1], dtype=np.int64)
        data_values = np.asarray(counts['data'][ptr[0]:ptr[-1]], float)
        idx = np.asarray(counts['indices'][ptr[0]:ptr[-1]], int)
        raw = csr_matrix((data_values, idx, ptr - ptr[0]), shape=(end - begin, nf))
        lo, hi = np.searchsorted(rows, [begin, end]); selected = rows[lo:hi]
        B = raw[(selected - begin).astype(int), :].tocsr(); B.eliminate_zeros()
        assert np.isfinite(B.data).all() and (B.data > 0).all() and np.array_equal(B.data, np.rint(B.data))
        RNA = np.asarray(B.sum(1)).ravel(); detected = np.diff(B.indptr)
        assert (RNA > 0).all()
        vals = B[:, columns].toarray(); logvals = np.log1p(vals * 10000 / RNA[:, None])
        close(logvals, np.log(RNA[:, None] + vals * 10000) - np.log(RNA[:, None]), 1e-12)
        for k, physical in enumerate(selected):
            if int(physical) in manual:
                expected_rna, expected_detected, expected_values = manual[int(physical)]
                close(RNA[k], expected_rna, 0); close(vals[k], expected_values, 0)
                assert detected[k] == expected_detected
                checked_manual.add(int(physical))
        values = np.c_[np.log1p(RNA), np.log1p(detected), logvals, (vals > 0).astype(float),
                       ((vals[:, 0] > 0) & (vals[:, 1] > 0)).astype(float)]
        alt_values = np.c_[np.log(1 + RNA), np.log(1 + detected),
                           np.log(RNA[:, None] + vals * 10000) - np.log(RNA[:, None]),
                           (vals > 0).astype(float), ((vals[:, 0] > 0) & (vals[:, 1] > 0)).astype(float)]
        for j in range(9):
            sums[:, j] += np.bincount(groups[lo:hi], weights=values[:, j], minlength=ng)
        for group_code in np.unique(groups[lo:hi]):
            alternate_sums[group_code] += alt_values[groups[lo:hi] == group_code].sum(axis=0)
        seen += len(selected)
        if number % 20 == 0:
            save('private_progress.json', dict(stage='ALL_BASELINE_FIXED_PARENT_NATIVE_RNA',
                 completed_chunks=number + 1, total_chunks=len(np.unique(rows // 1024))))
    assert seen == len(rows) == int(totals.sum()) and checked_manual == set(manual)
    close(sums, alternate_sums, 1e-7)
    close(sums / totals[:, None], alternate_sums / totals[:, None], 1e-11)
    avg = sums / totals[:, None]
    T = avg[::2]; My = avg[1::2]
    sender = np.minimum(My[:, 2], My[:, 3]); receiver = T[:, 4]
    interface = np.log1p(sender * receiver)
    C = np.c_[T[:, :2], np.log1p(totals[::2]), totals[::2] / frame.total_cells.to_numpy(),
              My[:, :2], np.log1p(totals[1::2]), totals[1::2] / frame.total_cells.to_numpy(),
              frame.treatment.to_numpy(), frame.HPV.to_numpy()]
    names = np.array(['T_logRNA', 'T_loggenes', 'T_logcells', 'T_fraction', 'My_logRNA', 'My_loggenes',
                      'My_logcells', 'My_fraction', 'agent::treatment', 'agent::HPV'])
    questions = [('HNSCC_fixed_interface_20261020_01', interface, frame.regression.to_numpy(float)),
                 ('HNSCC_fixed_interface_20261020_02', sender, receiver)]
    save('private_source_only_native_capture_and_clinical_coverage.json', dict(source='GSE286827',
         all_current_materials=len(lib), original_author_patients=len(c), current_eligible_baseline_patients=len(frame),
         evaluated_baseline_materials=len(baseline), exclusion_reason_counts=dict(reasons),
         selected_parent_cells=int(totals.sum()), manually_reconstructed_full_CSR_rows=len(manual),
         unique_regression_values=int(frame.regression.nunique()),
         sender_joint_positive_quantiles=np.quantile(My[:, 8], [0, .25, .5, .75, 1]).tolist(),
         receiver_positive_quantiles=np.quantile(T[:, 7], [0, .25, .5, .75, 1]).tolist(),
         receptor_all_zero_patients=int((receiver == 0).sum()), sender_all_zero_patients=int((sender == 0).sum()),
         author_binary_concordant_patients=int(frame.binary_concordant.sum()),
         selected_rows_SHA256=array_sha(rows), selected_unit_order_SHA256=array_sha(frame.patient.to_numpy(str)),
         source_numeric_input_SHA256=array_sha(np.c_[C, interface, sender, receiver, frame.regression.to_numpy(float), fine]),
         maximum_independent_pool_mean_error=float(np.max(np.abs(sums / totals[:, None] - alternate_sums / totals[:, None]))),
         identities_frozen_not_gamma_specific=True, real_arrays_not_persisted=True))
    result = infer(questions, C, names, (2 * frame.treatment + frame.HPV).to_numpy(int), fine,
                   frame.binary_concordant.to_numpy(bool), save)
    save('private_complete2_HNSCC_source_associations.json', result)
    files = {p.name: dict(bytes=p.stat().st_size, sha256=sha(p)) for p in out.iterdir() if p.is_file()}
    assert sum(v['bytes'] for v in files.values()) < 4 * 1024 ** 2
    terminal = dict(status='COMPLETE_ALL_HNSCC_BASELINE_FIXED_INTERFACE_AND_TWO_QUESTIONS_NUMERIC_REPLAY_PENDING_REVIEW',
         completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), files=files,
         registered_questions=2, no_human_RNA_or_patient_input_or_random_arrays_persisted=True,
         new_external_clinical_table_or_RNA_download=False,
         scope='Shared native selection/metadata I-O, bounded manual CSR; all model/null/random values independently projected, not new raw-experiment or untouched clinical validation')
    save('private_terminal.json', terminal)
    print(json.dumps(dict(status=terminal['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
