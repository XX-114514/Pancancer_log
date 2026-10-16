"""Fixed public RNA diagnostics, blind to ongoing coordination effects."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource

YGENES = 'KDM5D DDX3Y EIF1AY UTY ZFY USP9Y TMSB4Y RPS4Y1'.split()
STRESS = 'FOS JUN JUNB JUND ATF3 EGR1 HSPA1A HSPA1B HSP90AA1 DNAJB1 DDIT3'.split()


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    budget = resource(root, 8*1024**2)
    count = root/'NSCLC_native_MFA_external_frozen_role_count_inputs_20261016_v2'
    assay = root/'NSCLC_external_fixed_gene_programme_assay_20261016_v1'
    assert sha(count/'private_independent_native_count_input_audit.json') == '7bd028cc1c89ac7324ac25c65042b5b184cb377162a38f9544caf7bb24048c04'
    assert sha(assay/'private_independent_fixed_programme_assay_audit.json') == 'd745f978950b0882f5ff30f3c975cd6c1e626296b2a55e4ab4464b081d44d936'
    out = root/'NSCLC_external_public_RNA_diagnostics_20261016_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(), code_sha256=sha(__file__), budget=budget,
                goal='All61 sources/6760 pools, fixed public RNA diagnostics before reading external coordination effects. No extra P, no alteration to primary family.',
                sex_linked_RNA_proxy='mean log2(1+CPM) of measured fixed Y-linked genes minus XIST log2(1+CPM); require XIST and >=6/8 Y genes physically measured. Not a clinical sex assignment.',
                Y_genes=YGENES, stress_genes=STRESS, stress='mean log2(1+CPM) of fixed IEG/stress panel; require >=9/11 genes physically measured. Zero expression is observed, unmeasured is missing.',
                normalization='Reuse accepted same-pool effective common-native TMM library size; no new reference/feature learning.',
                followup='On exactly the already eligible test cases, own-library matched and same patient averages, descriptively add both RNA sex proxies and/or both stress readouts to the full primary nuisance span. No alternate primary, no new P or case rescue; missing required readouts means unavailable. Shared biological states, handling and demographic explanations remain hypotheses, not contamination diagnoses.',
                coordination_results_loaded=False)
    (out/'plan_before_any_diagnostic_RNA_measurement.json').write_text(json.dumps(plan, indent=2)+'\n')
    axis = pd.read_csv(count/'private_complete_native_gene_axis.tsv', sep='\t')
    assert axis.native_column.tolist() == list(range(31831)) and axis.gene.is_unique
    lookup = {str(g): int(i) for i, g in enumerate(axis.gene)}
    assert {'XIST', *YGENES, *STRESS} <= set(lookup)
    units = pd.read_csv(assay/'private_all6760_programme_units_metadata_quality_TMM.tsv.gz', sep='\t')
    manifest = pd.read_csv(count/'private_all61_count_archive_manifest.tsv', sep='\t')
    output = np.full((6760, 2), np.nan)
    coverage = []
    for item in manifest.itertuples():
        assert sha(count/item.archive) == item.sha256
        with np.load(count/item.archive, allow_pickle=False) as z:
            obs, cols, rows = z['counts'], z['common_native_columns'], z['global_unit_indices']
        inverse = {int(g): i for i, g in enumerate(cols)}
        ygenes = [g for g in YGENES if lookup[g] in inverse]
        stress = [g for g in STRESS if lookup[g] in inverse]
        sex_gate = lookup['XIST'] in inverse and len(ygenes) >= 6
        stress_gate = len(stress) >= 9
        valid = units.iloc[rows].RNA_measurement_available.to_numpy(bool)
        ids = rows[valid]
        denom = units.effective_common_library_size.to_numpy()[ids]
        assert np.isfinite(denom).all() and (denom > 0).all()
        if sex_gate:
            yi = [inverse[lookup[g]] for g in ygenes]
            y = np.log2(1+obs[valid][:, yi]*1e6/denom[:, None]).mean(axis=1)
            x = np.log2(1+obs[valid, inverse[lookup['XIST']]]*1e6/denom)
            output[ids, 0] = y-x
        if stress_gate:
            si = [inverse[lookup[g]] for g in stress]
            output[ids, 1] = np.log2(1+obs[valid][:, si]*1e6/denom[:, None]).mean(axis=1)
        coverage.append(dict(arm=item.arm, dataset=item.dataset, fixed_Y_genes_measured=len(ygenes),
                             XIST_measured=lookup['XIST'] in inverse, fixed_stress_genes_measured=len(stress),
                             RNA_sex_proxy_gate=sex_gate, RNA_stress_gate=stress_gate))
    np.savez_compressed(out/'private_all6760_public_RNA_diagnostics.npz', values=output,
                        labels=np.array(['Y_minus_XIST_RNA_proxy', 'IEG_stress_RNA']), unit_index=np.arange(6760))
    pd.DataFrame(coverage).to_csv(out/'private_all61_public_diagnostic_gene_coverage.tsv', sep='\t', index=False)
    terminal = dict(status='COMPLETE_ALL61_PUBLIC_RNA_DIAGNOSTIC_MEASUREMENTS_NO_COORDINATION_OR_CLINICAL_EFFECTS_LOADED',
                    finished_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
                    files={p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()})
    assert sum(v['bytes'] for v in terminal['files'].values()) < budget['output_cap_bytes']
    (out/'private_terminal.json').write_text(json.dumps(terminal, indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'], finished_UTC=terminal['finished_UTC'], terminal_sha256=sha(out/'private_terminal.json'))))


if __name__ == '__main__':
    main()
