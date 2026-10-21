"""Finite same-patient explanation of a previously reviewed wide-count result.

Eight descriptive points, no new significance tests. Count-level version
correspondence never becomes automatic per-cell identity assignment.
"""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import importlib.util
import json
import os

import numpy as np


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def import_source(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parent = Path(os.environ['PANCANCER_FINE_STRUCTURE_OUTPUT'])
    auditpath = Path(os.environ['PANCANCER_FINE_STRUCTURE_AUDIT'])
    acceptance = json.loads(auditpath.read_text())
    assert acceptance['parent_terminal_SHA256'] == sha(parent / 'private_terminal.json')
    producer_source = parent.parent / 'describe_CM02_fine_state_structure_v1.py'
    auditor_source = auditpath.parent.parent / 'audit_CM02_fine_state_structure_v1.py'
    assert sha(producer_source) == json.loads((parent / 'private_terminal.json').read_text())['source_SHA256']
    assert sha(auditor_source) == acceptance['source_SHA256']
    producer = import_source(producer_source, 'pinned_descriptive_producer')
    auditor = import_source(auditor_source, 'pinned_SVD_auditor')
    out = Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert out.parent.parent == Path('/dev/shm') and not out.exists()
    out.mkdir(mode=0o700)
    save = lambda name, value: (out / name).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    save('plan_before_same_patient_denominator_decomposition.json', dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_SHA256=sha(__file__), phases=['pre', 'post'], representations=['original_counts_legacy_fractions_and_covariates', 'current_counts_legacy_fractions_and_covariates', 'current_lineage_fractions_legacy_covariates', 'current_lineage_fractions_lineage_covariates'], shared_patient_contract='intersection of accepted source-fine reference patients and both original/current legacy parents >=20; never change samples between four representations', focal_pair='previous count-version correspondence only; not an identity or clinical discovery', parent_covariates='two parent log counts plus same treatment and HPV', all_eight_points_checked_by_previously_accepted_independent_SVD=True, no_P_bootstrap_or_new_pair_selection=True, no_patient_rows_saved=True))
    labels, tab, cur, contexts, unique, nodes, catalog = producer.load_inputs()
    focal = json.loads(os.environ['PANCANCER_FOCAL_STATE_PAIR'])
    assert focal[0] in labels.CD4[:4] and focal[1] in labels.CD8
    versions = json.loads((Path(os.environ['PANCANCER_PATIENT_SCOPE_OUTPUT']) / 'private_all_count_version_comparisons.json').read_text())
    assert any(r['phase'] == 'all' and r['count_scope'] == 'original_author_reference' and r['original_count_label'] == 'Teff' and r['fine_labels'] == [focal[1]] and r['all_libraries_exact'] for r in versions)
    old_other = labels.OLD[4:]
    legacy_correspondence = {k: k for k in labels.OLD}
    # This is a count correspondence verified in the full earlier version audit.
    legacy_correspondence['Teff'] = focal[1]
    fine_labels = labels.CD4 + labels.CD8
    output = []
    for phase in ['pre', 'post']:
        g, x, y, c, strata, denominators = producer.materialize_context(labels, tab, cur, contexts['source_fine_reference'], phase)
        original = g[labels.OLD].to_numpy(float)
        current_legacy = x[:, [fine_labels.index(legacy_correspondence[k]) for k in labels.OLD]]
        original_parents = np.column_stack([original[:, :4].sum(axis=1), original[:, 4:].sum(axis=1)])
        legacy_parents = np.column_stack([current_legacy[:, :4].sum(axis=1), current_legacy[:, 4:].sum(axis=1)])
        keep = (original_parents >= 20).all(axis=1) & (legacy_parents >= 20).all(axis=1)
        g, x, original, original_parents, legacy_parents, c = g.loc[keep], x[keep], original[keep], original_parents[keep], legacy_parents[keep], c[keep]
        nuisance = c[:, 3:]
        covariates = lambda parents: np.column_stack([np.ones(len(g)), np.log1p(parents), nuisance])
        original_pair = np.column_stack([original[:, labels.OLD.index(focal[0])], original[:, labels.OLD.index('Teff')]])
        current_pair = x[:, [fine_labels.index(k) for k in focal]]
        lineage_parents = np.column_stack([x[:, :7].sum(axis=1), x[:, 7:].sum(axis=1)])
        representations = [('original_counts_legacy_fractions_and_covariates', original_pair / original_parents, covariates(original_parents)), ('current_counts_legacy_fractions_and_covariates', current_pair / legacy_parents, covariates(legacy_parents)), ('current_lineage_fractions_legacy_covariates', current_pair / lineage_parents, covariates(legacy_parents)), ('current_lineage_fractions_lineage_covariates', current_pair / lineage_parents, covariates(lineage_parents))]
        for name, values, design in representations:
            matrix = np.zeros((len(g), 13))
            matrix[:, 0], matrix[:, 7] = values[:, 0], values[:, 1]
            r = producer.effect(matrix, design)[0, 0]
            independent = auditor.effect(matrix, design)[0, 0]
            assert np.isfinite(r) == np.isfinite(independent)
            if np.isfinite(r):
                assert abs(r-independent) < 1e-9
            output.append(dict(phase=phase, representation=name, focal_CD4_reference=focal[0], focal_CD8_reference=focal[1], shared_eligible_author_patients=len(g), partial_Spearman=float(r) if np.isfinite(r) else None, same_patient_set_across_all_representations=True, independent_SVD_checked=True, no_P_or_new_identity_claim=True))
    assert len(output) == 8
    save('private_all8_same_patient_representation_points.json', output)
    terminal = dict(status='COMPLETE_ALL8_SAME_PATIENT_REPRESENTATION_POINTS_WITH_INDEPENDENT_SVD', completed_UTC=datetime.now(timezone.utc).isoformat(), source_SHA256=sha(__file__), parent_acceptance_SHA256=sha(auditpath), no_individual_rows_exported=True, no_new_clinical_spatial_or_confirmatory_claim=True, files={p.name: dict(SHA256=sha(p), bytes=p.stat().st_size) for p in out.iterdir()})
    save('private_terminal.json', terminal)
    print(json.dumps(dict(status=terminal['status'], terminal_SHA256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
