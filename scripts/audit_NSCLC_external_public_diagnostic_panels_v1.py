"""Recompute every public diagnostic readout from accepted count archives."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    resource(root, 1048576)
    src = root/'NSCLC_external_public_RNA_diagnostics_20261016_v1'
    terminal = json.loads((src/'private_terminal.json').read_text())
    for name, item in terminal['files'].items():
        assert sha(src/name) == item['sha256']
    plan = json.loads((src/'plan_before_any_diagnostic_RNA_measurement.json').read_text())
    count = root/'NSCLC_native_MFA_external_frozen_role_count_inputs_20261016_v2'
    assay = root/'NSCLC_external_fixed_gene_programme_assay_20261016_v1'
    genes = pd.read_csv(count/'private_complete_native_gene_axis.tsv', sep='\t').set_index('gene').native_column
    units = pd.read_csv(assay/'private_all6760_programme_units_metadata_quality_TMM.tsv.gz', sep='\t')
    manifest = pd.read_csv(count/'private_all61_count_archive_manifest.tsv', sep='\t')
    covered = pd.read_csv(src/'private_all61_public_diagnostic_gene_coverage.tsv', sep='\t').set_index(['arm', 'dataset'])
    with np.load(src/'private_all6760_public_RNA_diagnostics.npz', allow_pickle=False) as z:
        values = z['values']
        assert z['unit_index'].tolist() == list(range(6760))
        assert z['labels'].astype(str).tolist() == ['Y_minus_XIST_RNA_proxy', 'IEG_stress_RNA']
    expected = np.full((6760, 2), np.nan)
    for item in manifest.itertuples():
        with np.load(count/item.archive, allow_pickle=False) as z:
            obs, columns, ids = z['counts'], z['common_native_columns'], z['global_unit_indices']
        loc = {int(g): i for i, g in enumerate(columns)}
        ys = [loc[int(genes[g])] for g in plan['Y_genes'] if int(genes[g]) in loc]
        stress = [loc[int(genes[g])] for g in plan['stress_genes'] if int(genes[g]) in loc]
        cov = covered.loc[item.arm, item.dataset]
        xist = int(genes['XIST']) in loc
        assert cov.fixed_Y_genes_measured == len(ys) and cov.fixed_stress_genes_measured == len(stress)
        assert cov.XIST_measured == xist
        assert cov.RNA_sex_proxy_gate == (xist and len(ys) >= 6)
        assert cov.RNA_stress_gate == (len(stress) >= 9)
        for local, unit in enumerate(ids):
            if not units.RNA_measurement_available.iloc[unit]:
                continue
            denom = units.effective_common_library_size.iloc[unit]
            if xist and len(ys) >= 6:
                y = sum(float(np.log1p(float(obs[local, i])*1e6/denom)/np.log(2)) for i in ys)/len(ys)
                x = float(np.log1p(float(obs[local, loc[int(genes['XIST'])]])*1e6/denom)/np.log(2))
                expected[unit, 0] = y-x
            if len(stress) >= 9:
                expected[unit, 1] = sum(float(np.log1p(float(obs[local, i])*1e6/denom)/np.log(2)) for i in stress)/len(stress)
    assert np.array_equal(np.isnan(expected), np.isnan(values))
    assert np.allclose(expected, values, equal_nan=True, rtol=1e-12, atol=1e-12)
    proof = src/'private_independent_public_RNA_diagnostic_audit.json'
    assert not proof.exists()
    report = dict(status='PASS_ALL61_PANEL_COVERAGE_AND_ALL6760_POSITIONS_NUMERIC_AND_MISSINGNESS',
                  finished_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
                  producer_terminal_sha256=sha(src/'private_terminal.json'),
                  maximum_numeric_error=float(np.nanmax(abs(expected-values))),
                  no_coordination_or_clinical_effects_loaded=True)
    proof.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(status=report['status'], finished_UTC=report['finished_UTC'], audit_sha256=sha(proof))))


if __name__ == '__main__':
    main()
