"""Freeze outcome-free case/source partitions and source-native gene opportunities."""
import hashlib
import json
import os
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget = resource(root, 16 * 1024**2)
    ready = root / 'immune_CM_RNA_supported_multicell_views_preflight_20261018_v1'
    acceptance = ready / 'private_independent_count_readiness_acceptance_20261018_v2.json'
    assert sha(acceptance) == '3f1261b9083f1cefcf0e072daecdef7c3db5d1fbf7f44b4ac1ff237efb126c15'
    table = pd.read_csv(ready / 'private_all_original_libraries_multicell_view_counts.tsv', sep='\t', keep_default_na=False)
    assert len(table) == 236
    latest = root / 'NSCLC_native_MFA_external_lineage_qualification_20261016_v2'
    terminal = json.loads((latest / 'private_terminal.json').read_text())
    filename = 'private_all1690_library_metadata_and_label_diagnostics.tsv'
    assert sha(latest / filename) == terminal['files'][filename]['sha256']
    # Only identity and sampling/design fields; no outcome values used.
    fields = ['arm', 'dataset', 'library_id', 'patient_block', 'reviewed_patient', 'resampling_block', 'patient_evidence',
              'cancer_label', 'phase_class', 'capture_scope', 'ex_vivo', 'annotated_scope_cells', 'latest_sidecar']
    lib = pd.read_csv(latest / filename, sep='\t', usecols=fields, keep_default_na=False)
    lib = lib[lib.arm.eq('archived_immune_atlas')].copy()
    lib['sampleID'] = lib.dataset.astype(str) + '::' + lib.library_id.astype(str)
    assert lib.sampleID.is_unique
    joined = table.merge(lib.drop(columns=['dataset']), on='sampleID', how='left', validate='many_to_one', indicator=True)
    assert joined['_merge'].eq('both').all()
    joined['case_evidence'] = np.where(joined.reviewed_patient, joined.patient_evidence, 'unreviewed_library_case_fallback')
    joined['case_block'] = joined.dataset.astype(str) + '::' + joined.resampling_block.astype(str)
    assert (joined.loc[joined.reviewed_patient, 'patient_block'].astype(str).str.len() > 0).all()
    selected = joined[joined.threshold.eq(1) & joined.eligible].copy()
    counts = selected.groupby('dataset').size()
    qualifying_sources = set(counts[counts >= 10].index)
    discovery = selected[selected.reviewed_patient & selected.dataset.isin(qualifying_sources)].copy()
    assert set(discovery.dataset) == {'GSE286827', 'GSE294109'}
    phase_priority = {'pre': 0, 'untreated': 0, 'post': 1, 'on_treatment': 1, 'progression_or_recurrence': 2, 'unknown': 3}
    discovery['phase_priority'] = discovery.phase_class.map(phase_priority)
    assert discovery.phase_priority.notna().all()
    discovery['selection_hash'] = discovery.sampleID.map(lambda x: hashlib.sha256(('2026101802|' + str(x)).encode()).hexdigest())
    discovery = discovery.sort_values(['case_block', 'phase_priority', 'selection_hash']).drop_duplicates('case_block')
    discovery['partition'] = 'reviewed_case_discovery'
    held = selected[selected.dataset.eq('GSE305078')].copy();held['partition'] = 'TNBC_source_holdout_sample_identity_pending'
    small = selected[~selected.dataset.isin(qualifying_sources)].copy();small['partition'] = 'small_source_projection_only'
    chosen = pd.concat([discovery, held, small], ignore_index=True)
    assert chosen.sampleID.is_unique and len(discovery) >= 20 and len(held) >= 10
    assert discovery.case_block.is_unique and set(discovery.case_block).isdisjoint(held.case_block)
    chosen['fit_sample'] = ['MCP_SAMPLE_%04d' % i for i in range(len(chosen))]
    atlas = root.parents[1] / 'runs/20261013_immunotherapy_atlas_v1'
    fixed_path = root.parents[1] / 'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json'
    fixed = json.loads(fixed_path.read_text());group_genes = set(g for gg in fixed['panels'].values() for g in gg)
    z = zarr.open_group(str(atlas / 'outputs/immune_validation_authority_v1.zarr'), mode='r')
    order = list(np.asarray(read_elem(z['uns/source_order'])).astype(str))
    measured = np.asarray(read_elem(z['varm/measured_in_source']), bool)
    symbol = np.asarray(read_elem(z['var/gene_symbol'])).astype(str)
    canonical = np.asarray(read_elem(z['var/canonical_gene_identity'])).astype(str)
    identity = np.asarray(read_elem(z['var/feature_identity_status'])).astype(str)
    aliases = np.asarray(read_elem(z['var/source_gene_aliases'])).astype(str)
    source_names = sorted(chosen.dataset.unique())
    required = sorted(set(discovery.dataset) | set(held.dataset))
    name_maps = {};forbidden = set()
    for source in source_names:
        si = order.index(source);names = defaultdict(set)
        for col in np.flatnonzero(measured[:, si]):
            labels = {symbol[col], canonical[col]}
            if identity[col] == 'stable_Ensembl_ID': labels.update(aliases[col].split('|'))
            if labels & group_genes: forbidden.add(int(col))
            if canonical[col].startswith('AMBIGUOUS::') or identity[col] == 'ambiguous_source_feature': continue
            for name in labels:
                if re.fullmatch(r'[A-Za-z][A-Za-z0-9._-]*', name) and not name.startswith('ENSG') and name not in {'nan', 'None', '__SOURCE_VALUE_MISSING__'}:
                    names[name].add(int(col))
        name_maps[source] = {name: next(iter(cols)) for name, cols in names.items() if len(cols) == 1}
    common = set.intersection(*(set(name_maps[s]) for s in required))
    by_tuple = defaultdict(list)
    for gene in common:
        columns = tuple(name_maps[s][gene] for s in required)
        if set(columns) & forbidden or gene in group_genes or re.match(r'^(MT-|RPL|RPS)', gene): continue
        by_tuple[columns].append(gene)
    genes = []
    for columns, alternatives in by_tuple.items():
        genes.append(sorted(alternatives, key=lambda name: (-sum(symbol[col] == name for col in columns), name))[0])
    collisions = set()
    for s in required:
        repeated = pd.Series([name_maps[s][g] for g in genes]).duplicated(keep=False).to_numpy()
        collisions.update(g for g, dup in zip(genes, repeated) if dup)
    genes = sorted(set(genes) - collisions)
    assert len(genes) >= 2000 and not (set(genes) & group_genes)
    rows = []
    for gene in genes:
        row = dict(gene=gene)
        for s in source_names:
            col = name_maps[s].get(gene, -1)
            row[s] = -1 if col in forbidden else col
        rows.append(row)
    mapping = pd.DataFrame(rows)
    assert (mapping[required] >= 0).all().all()
    out = root / 'immune_RNA_MCP_discovery_scope_20261018_v1';assert not out.exists();out.mkdir(mode=0o700)
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), budget=budget,
        count_readiness_acceptance_sha256=sha(acceptance), current_library_metadata_sha256=sha(latest / filename),
        frozen_panel_sha256=sha(fixed_path), atlas_release_sha256=sha(atlas / 'manifests/immune_authority_release_public.json'),
        selection='Primary count threshold1, original eligible libraries. Reviewed discovery cases represented once: pre first among eligible, then post/on, progression, unknown, tie by fixed metadata-only hash. No response fields loaded.',
        discovery_sources=sorted(discovery.dataset.unique()), source_holdout='GSE305078 kept entirely outside all gene expression fitting; exact clinical patient alias pending. Its pre stage is supported by newer official sampling metadata, not old empty CM phase.',
        stage='Latest reviewed sampling phase overlay; old CM table unchanged. Discovery mixes known/unknown sampling context and cannot be called all pretreatment.',
        genes='Unique source-native symbol/canonical/stable-ID alias, shared across discovery+TNBC. Collapse identical physical tuples preferring deposited direct symbol; remove residual physical collisions. All34 grouping/clinical markers and ALL matching physical aliases excluded, plus mitochondrial/ribosomal genes. Optional RMC missing measurement stays-1.',
        candidate_features='Per view: training-only >=1% raw positive and positive logCP10K variance; at most5000 highest variance candidates,2000 PCA genes,20 PCs. Full native depth normalization. No held-out variance/clinical outcomes used.',
        sampling='Up to100 cells per chosen material/view, fixed hash per global row/view; original view count/proportions retained separately. R-only and effector-only disjoint, whole original Myeloid third view; both/neither held outside primary fitting.',
        prospective_model='Official isolated DIALOGUE runtime. Primary K2 with three fixed seeds20261018/20261019/20261020; K3 seed20261018 sensitivity. All4 fit statuses retained; stability/heldout assessment required, never clinical selection of K.',
        prospective_controls=['cell log fullRNA', 'cell log detected genes', 'discovery-source indicator', 'original full view fraction in all native material cells', 'log1p original full view cells'],
        caveat='Patient IDs not promoted by unique library strings. Source-heldout data already entered original CM definition and is not untouched independent CM validation. Old312/30/18 families remain closed. No new MCP or clinical fit here.')
    joined.drop(columns='_merge').to_csv(out / 'private_all236_current_case_phase_overlay.tsv', sep='\t', index=False)
    chosen.to_csv(out / 'private_metadata_only_case_source_partition.tsv', sep='\t', index=False)
    mapping.to_csv(out / 'private_marker_excluded_native_gene_map.tsv', sep='\t', index=False)
    (out / 'plan_before_expression_or_MCP_fitting.json').write_text(json.dumps(plan, indent=2) + '\n')
    summary = []
    for (s, partition), g in chosen.groupby(['dataset', 'partition'], sort=True):
        summary.append(dict(dataset=s, partition=partition, selected_materials=len(g), original_eligible_materials=int(counts[s]),
            reviewed_cases=int(g.loc[g.reviewed_patient, 'case_block'].nunique()), phase_counts=g.phase_class.value_counts().to_dict(),
            original_CM_phase_vs_latest_disagreements=int(g.timepoint.ne(g.phase_class).sum()), genes_complete=int((mapping[s] >= 0).sum())))
    (out / 'private_source_design_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    files = [p for p in out.iterdir() if p.is_file()]
    t = dict(status='FROZEN_REVIEWED_CASE_DISCOVERY_TNBC_SOURCE_HOLDOUT_AND_MARKER_EXCLUDED_GENE_UNIVERSE_NO_EXPRESSION_FIT',
        completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), common_candidate_genes=len(genes),
        files={p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in files}, no_MCP_or_clinical_fit=True)
    (out / 'private_terminal.json').write_text(json.dumps(t, indent=2) + '\n')
    print(json.dumps(dict(status=t['status'], completed_UTC=t['completed_UTC'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
