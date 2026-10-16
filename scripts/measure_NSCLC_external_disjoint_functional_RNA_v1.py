"""All-source functional RNA using genes excluded from every programme view."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import load, sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from disjoint_hallmark_ULM_kernel_v1 import score


def main():
    data = load()
    root = data['root']
    budget = resource(root, 32*1024**2)
    count = root/'NSCLC_native_MFA_external_frozen_role_count_inputs_20261016_v2'
    assay = root/'NSCLC_external_fixed_gene_programme_assay_20261016_v1'
    assert sha(count/'private_independent_native_count_input_audit.json') == '7bd028cc1c89ac7324ac25c65042b5b184cb377162a38f9544caf7bb24048c04'
    assert sha(assay/'private_independent_fixed_programme_assay_audit.json') == 'd745f978950b0882f5ff30f3c975cd6c1e626296b2a55e4ab4464b081d44d936'
    out = root/'NSCLC_external_disjoint_functional_RNA_20261016_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(), code_sha256=sha(__file__), budget=budget,
                all_scope='All61 sources/all6760 pool positions/all50 frozen Hallmark references, not selected by coordination or clinical effects.',
                assay='Binary-membership ULM of same-pool native pseudobulk log2(1+CPM), target and non-target means, pooled expressed-target counts. Reuse accepted same-pool TMM effective library sizes. Explicit external pseudobulk readout, not exact original mean-per-cell function values.',
                excluded='Union of all four frozen original feature sets excluded from target AND background; only physically measured unique source-common native genes retained. No new feature/gene-set learning or missing-gene zero fill.',
                coverage='>=20 remaining measured targets, >=500 measured non-target genes. All50 retained, absent/low-cell-count/RNA-unmeasurable pools kept missing; constant remaining RNA yields missing ULM, not pathway absence.',
                followup='Describe all50 functions on the exact same already eligible matched-library/patient pair cases, controlling the same own/partner full quality and category span. No new P or primary family. Gene disjointness is not independent-patient validation; whole-RNA denominator and shared sample biology can still couple readouts.',
                support='Pooled expressed targets are genes with count>0 in the pool, not per-cell positive percentages; low support is flagged descriptively, not used to select a favourable function.',
                reference_provenance_sha256=sha(data['refroot']/'reference_provenance.json'),
                excluded_union_index_sha256=__import__('hashlib').sha256(data['union'].astype('<i8').tobytes()).hexdigest(),
                ULM_kernel_sha256=sha(Path(__file__).with_name('disjoint_hallmark_ULM_kernel_v1.py')),
                coordination_or_clinical_effects_loaded=False, clinical_or_direct_LR_validation=False)
    (out/'plan_before_all_disjoint_functional_measurements.json').write_text(json.dumps(plan, indent=2)+'\n')
    axis = pd.read_csv(count/'private_complete_native_gene_axis.tsv', sep='\t')
    assert np.array_equal(axis.gene.to_numpy(str), data['genes'])
    units = pd.read_csv(assay/'private_all6760_programme_units_metadata_quality_TMM.tsv.gz', sep='\t')
    manifest = pd.read_csv(count/'private_all61_count_archive_manifest.tsv', sep='\t')
    names = sorted(data['refs'])
    output = {k: np.full((6760, 50), np.nan) for k in ['ULM_score', 'target_mean', 'nontarget_mean', 'expressed_target_genes']}
    coverage = []
    for item in manifest.itertuples():
        assert sha(count/item.archive) == item.sha256
        with np.load(count/item.archive, allow_pickle=False) as z:
            obs, columns, rows = z['counts'], z['common_native_columns'], z['global_unit_indices']
        keep = ~np.isin(columns, data['union'])
        columns = columns[keep]
        genes = data['genes'][columns]
        targets = []
        qualified = []
        for j, name in enumerate(names):
            members = np.flatnonzero(np.isin(genes, list(data['refs'][name])))
            good = len(members) >= 20 and len(genes)-len(members) >= 500
            coverage.append(dict(arm=item.arm, dataset=item.dataset, reference=name,
                                 frozen_reference_genes=len(data['refs'][name]),
                                 remaining_measured_target_genes=len(members),
                                 remaining_measured_nontarget_genes=len(genes)-len(members),
                                 excluded_measured_programme_genes=int((~keep).sum()), qualified=good))
            if good:
                qualified.append(j)
                targets.append(members)
        available = units.iloc[rows].RNA_measurement_available.to_numpy(bool)
        ids = rows[available]
        if len(ids) and targets:
            effective = units.effective_common_library_size.to_numpy()[ids]
            expression = np.log2(1+obs[available][:, keep].astype(float)*1e6/effective[:, None])
            result = score(expression, targets)
            for key in output:
                output[key][np.ix_(ids, np.asarray(qualified))] = result[key]
        (out/'private_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(), completed_source_index=int(item.source_index), total_sources=61))+'\n')
    np.savez_compressed(out/'private_all6760_by50_disjoint_functional_RNA.npz', **output,
                        functions=np.array(names), unit_index=np.arange(6760))
    pd.DataFrame(coverage).to_csv(out/'private_all3050_source_function_gene_coverage.tsv', sep='\t', index=False)
    terminal = dict(status='COMPLETE_ALL61_SOURCES6760_POOLS50_DISJOINT_FUNCTION_READOUTS_AUDIT_PENDING',
                    finished_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
                    files={p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()},
                    clinical_or_coordination_effects_loaded=False, new_P_or_primary_family=False)
    assert sum(v['bytes'] for v in terminal['files'].values()) < budget['output_cap_bytes']
    (out/'private_terminal.json').write_text(json.dumps(terminal, indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'], finished_UTC=terminal['finished_UTC'], terminal_sha256=sha(out/'private_terminal.json'))))


if __name__ == '__main__':
    main()
