"""All-universe functional readouts independently by two-group pooled t algebra."""
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import load, sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    data = load()
    root = data['root']
    src = root/'NSCLC_external_disjoint_functional_RNA_20261016_v1'
    h = root/'HNSCC_receiver_followup_job_control_20261016_v2/measure_NSCLC_external_disjoint_functional_RNA_v1'
    receipt = json.loads((h/'launch_receipt.json').read_text())
    deadline = time.monotonic()+3600
    while not (src/'private_terminal.json').exists():
        p = Path('/proc')/str(receipt['pid'])/'stat'
        assert p.exists(), 'Original functional producer missing without terminal'
        f = p.read_text().rsplit(')', 1)[1].split()
        assert f[0] != 'Z' and f[19] == str(receipt['proc_start_ticks'])
        assert time.monotonic() < deadline, 'Observation expired, do not infer termination'
        time.sleep(3)
    resource(root, 1048576)
    terminal = json.loads((src/'private_terminal.json').read_text())
    assert terminal['source_sha256'] == receipt['source_sha256']
    for name, item in terminal['files'].items():
        assert sha(src/name) == item['sha256']
    count = root/'NSCLC_native_MFA_external_frozen_role_count_inputs_20261016_v2'
    units = pd.read_csv(root/'NSCLC_external_fixed_gene_programme_assay_20261016_v1/private_all6760_programme_units_metadata_quality_TMM.tsv.gz', sep='\t')
    coverage = pd.read_csv(src/'private_all3050_source_function_gene_coverage.tsv', sep='\t').set_index(['arm', 'dataset', 'reference'])
    manifest = pd.read_csv(count/'private_all61_count_archive_manifest.tsv', sep='\t')
    with np.load(src/'private_all6760_by50_disjoint_functional_RNA.npz', allow_pickle=False) as z:
        names = z['functions'].astype(str).tolist()
        saved = {k: z[k] for k in ['ULM_score', 'target_mean', 'nontarget_mean', 'expressed_target_genes']}
        assert z['unit_index'].tolist() == list(range(6760))
    assert names == sorted(data['refs']) and len(coverage) == 3050
    expected = {k: np.full((6760, 50), np.nan) for k in saved}
    forbidden = set(int(i) for i in data['union'])
    max_error = {k: 0. for k in saved}
    checkpoint = root/'HNSCC_receiver_followup_job_control_20261016_v2'/Path(__file__).stem/'private_audit_checkpoint.json'
    for item in manifest.itertuples():
        with np.load(count/item.archive, allow_pickle=False) as z:
            obs, columns, rows = z['counts'], z['common_native_columns'], z['global_unit_indices']
        keep = np.array([int(i) not in forbidden for i in columns])
        genes = data['genes'][columns[keep]]
        assert not (set(int(i) for i in columns[keep]) & forbidden)
        valid = units.iloc[rows].RNA_measurement_available.to_numpy(bool)
        ids = rows[valid]
        denom = units.effective_common_library_size.to_numpy()[ids]
        x = np.log1p(obs[valid][:, keep].astype(float)/denom[:, None]*1e6)/np.log(2)
        g = x.shape[1]
        total = x.sum(axis=1)
        squares = np.einsum('ij,ij->i', x, x)
        for j, name in enumerate(names):
            target = np.array([i for i, gene in enumerate(genes) if gene in data['refs'][name]], int)
            k = len(target)
            good = k >= 20 and g-k >= 500
            cov = coverage.loc[(item.arm, item.dataset, name)]
            assert cov.remaining_measured_target_genes == k and cov.remaining_measured_nontarget_genes == g-k
            assert cov.excluded_measured_programme_genes == (~keep).sum()
            assert cov.qualified == good
            if not good or not len(ids):
                continue
            ts = x[:, target].sum(axis=1)
            bs = total-ts
            mt, mb = ts/k, bs/(g-k)
            within = squares-ts*ts/k-bs*bs/(g-k)
            with np.errstate(invalid='ignore', divide='ignore'):
                t = (mt-mb)/np.sqrt(within/(g-2)*(1/k+1/(g-k)))
            expected['ULM_score'][ids, j] = t
            expected['target_mean'][ids, j] = mt
            expected['nontarget_mean'][ids, j] = mb
            expected['expressed_target_genes'][ids, j] = np.count_nonzero(x[:, target], axis=1)
        checkpoint.write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(), completed_source_index=int(item.source_index), total_sources=61))+'\n')
    for key in saved:
        assert np.array_equal(np.isnan(expected[key]), np.isnan(saved[key])), key
        assert np.allclose(expected[key], saved[key], equal_nan=True, rtol=1e-9, atol=1e-9), key
        max_error[key] = float(np.nanmax(abs(expected[key]-saved[key])))
    proof = src/'private_independent_disjoint_functional_RNA_audit.json'
    assert not proof.exists()
    report = dict(status='PASS_ALL61_UNIVERSES3050_GATES6760_POOLS50_FUNCTIONS_DISJOINTNESS_T_ALGEBRA_AND_MISSINGNESS',
                  finished_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
                  producer_terminal_sha256=sha(src/'private_terminal.json'), maximum_errors=max_error,
                  upstream_count_archives_recomputed_not_new_raw_cell_replay=True,
                  coordination_or_clinical_effects_loaded=False, direct_communication_or_clinical_validation=False)
    proof.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(status=report['status'], finished_UTC=report['finished_UTC'], audit_sha256=sha(proof))))


if __name__ == '__main__':
    main()
