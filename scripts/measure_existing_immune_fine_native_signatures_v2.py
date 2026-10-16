"""V2: canonicalize sparse group products before unique-index accumulation.

All-gene means live only in RAM. Private small outputs retain every frozen
unit/slot, source-native availability, and distinct reasons for missingness.
"""
import hashlib
import json
import os
import subprocess
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
from scipy.sparse import csr_matrix


PANEL = ('CD3D CD3E TRAC CD4 CD8A CD8B IL7R CCR7 TCF7 LEF1 '
         'FOXP3 IL2RA CTLA4 IKZF2 TIGIT TNFRSF18 CCR8 GZMK GZMB PRF1 '
         'NKG7 GNLY CTSW KLRB1 KLRD1 FCGR3A ZNF683 CXCR6 CD69 '
         'PDCD1 TOX LAG3 HAVCR2 ENTPD1 NR4A1 IFNG CCL4 CCL5 MKI67 '
         'LYZ LST1 TYROBP FCER1G C1QA C1QB C1QC APOE CD68 CSF1R '
         'FCN1 S100A8 S100A9 VCAN MARCO FOLR2 SPP1 TREM2 CD14 '
         'CD74 HLA-DRA DNAJB1 HSPA1A PTPRC EPCAM KRT19 COL1A1 '
         'PECAM1 VWF MS4A1 CD79A ITGA5 ITGB1 LTA LTB LTBR '
         'CXCL12 CXCR4 TLR2 TLR4 LY96 HMGB1 SDC2 PTPRJ ADAM17').split()


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8388608), b''):
            h.update(block)
    return h.hexdigest()


def resource(root, cap=256 * 1024 ** 2):
    df = subprocess.check_output(['df', '-P', str(root.resolve())], text=True, timeout=20)
    filesystem = df.splitlines()[-1].split()[0]
    quota = subprocess.check_output(['quota', '-g', '-v', '-w', 'USER002'], text=True, timeout=20)
    rows = [s.split() for s in quota.splitlines() if s.split() and s.split()[0] == filesystem]
    assert len(rows) == 1
    remaining = (int(rows[0][3]) - int(rows[0][1].rstrip('*'))) * 1024
    assert remaining > 500000000000 + cap + 1000000000, 'Live quota reserve not satisfied'
    return remaining


def accumulate(x, unit, n_units, means, q, panel_cols, panel_out):
    """Sparse group projection, after canonicalizing only a working block."""
    x.sum_duplicates()
    x.eliminate_zeros()
    x.sort_indices()
    assert np.isfinite(x.data).all() and (x.data > 0).all()
    assert np.equal(x.data, np.floor(x.data)).all()
    totals = np.asarray(x.sum(axis=1)).ravel().astype(float)
    detected = np.diff(x.indptr).astype(float)
    assert (totals > 0).all()
    group = csr_matrix((np.ones(len(unit)), (unit, np.arange(len(unit)))), shape=(n_units, len(unit)))
    logged = x.astype(np.float64, copy=True)
    logged.data = np.log1p(logged.data * np.repeat(10000 / totals, np.diff(logged.indptr)))
    grouped_log = group @ logged
    grouped_log.sum_duplicates()
    grouped_log.sort_indices()
    assert grouped_log.has_canonical_format
    summed = grouped_log.tocoo()
    means[summed.row, summed.col] += summed.data
    q[:, 0] += np.bincount(unit, minlength=n_units)
    q[:, 1] += np.bincount(unit, weights=totals, minlength=n_units)
    q[:, 2] += np.bincount(unit, weights=np.log1p(totals), minlength=n_units)
    q[:, 3] += np.bincount(unit, weights=np.log1p(detected), minlength=n_units)
    for key, block in [('sum_raw', x[:, panel_cols]), ('sum_log', logged[:, panel_cols])]:
        z = (group @ block).tocoo()
        panel_out[key][z.row, z.col] += z.data
    for key, threshold in [('positive', 1), ('ge2', 2)]:
        block = x[:, panel_cols].copy()
        block.data = (block.data >= threshold).astype(np.int32)
        block.eliminate_zeros()
        z = (group @ block).tocoo()
        panel_out[key][z.row, z.col] += z.data.astype(np.int64)


def self_test():
    # Duplicate coordinates, explicit zero, unequal full-RNA denominators.
    x = csr_matrix((np.array([1, 2, 0, 1, 3, 2]), np.array([0, 0, 2, 1, 2, 1]),
                    np.array([0, 3, 5, 6])), shape=(3, 4))
    ids = np.array([0, 1, 0])
    m = np.zeros((2, 4)); q = np.zeros((2, 4))
    p = {k: np.zeros((2, 2), dtype=np.int64 if k in ['positive', 'ge2'] else float)
         for k in ['sum_raw', 'sum_log', 'positive', 'ge2']}
    accumulate(x.copy(), ids, 2, m, q, np.array([0, 1]), p)
    a = x.toarray(); expected = np.log1p(a * (10000 / a.sum(axis=1))[:, None])
    for i in range(2):
        assert np.allclose(m[i], expected[ids == i].sum(axis=0))
        assert np.array_equal(p['positive'][i], (a[ids == i][:, [0, 1]] >= 1).sum(axis=0))
        assert np.array_equal(p['ge2'][i], (a[ids == i][:, [0, 1]] >= 2).sum(axis=0))
    assert np.array_equal(q[:, 0], [2, 1])


def main():
    self_test()
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    scope = root / 'existing_immune_fine_native_measurement_scope_20261016_v1'
    external = root / 'existing_immune_cohort_candidate_extension_preflight_20261016_v2'
    audit = json.loads((external / 'private_independent_audit.json').read_text())
    assert audit['status'].startswith('PASS_ALL_ATLAS_')
    scope_terminal = json.loads((scope / 'private_terminal.json').read_text())
    for name, digest in scope_terminal['files'].items():
        assert sha(scope / name) == digest
    units = pd.read_csv(scope / 'private_all_frozen_T_myeloid_node_units.tsv', sep='\t',
                        dtype={'dataset': str, 'library_id': str, 'node_id': str})
    mapping = np.load(scope / 'private_original_cell_to_fine_node_unit.npy')
    nu = len(units)
    atlas_run = root.parents[2] / 'Project_v5/runs/20261013_immunotherapy_atlas_v1'
    z = zarr.open_group(str(atlas_run / 'outputs/immune_validation_authority_v1.zarr'), mode='r')
    counts = z['layers']['counts']
    n_cells, n_features = counts.attrs['shape']
    assert len(mapping) == n_cells and nu == 6957
    order = np.asarray(read_elem(z['uns']['source_order'])).astype(str).tolist()
    measured = np.asarray(read_elem(z['varm']['measured_in_source']), bool)
    symbol = np.asarray(read_elem(z['var']['gene_symbol'])).astype(str)
    canonical = np.asarray(read_elem(z['var']['canonical_gene_identity'])).astype(str)
    identity = np.asarray(read_elem(z['var']['feature_identity_status'])).astype(str)
    aliases = np.asarray(read_elem(z['var']['source_gene_aliases'])).astype(str)
    source_map = {}
    for si, source in enumerate(order):
        names = defaultdict(set)
        for col in np.flatnonzero(measured[:, si]):
            if canonical[col].startswith('AMBIGUOUS::') or identity[col] == 'ambiguous_source_feature':
                continue
            labels = {symbol[col], canonical[col]}
            if identity[col] == 'stable_Ensembl_ID':
                labels.update(aliases[col].split('|'))
            for gene in labels - {'', 'nan', 'None', '__SOURCE_VALUE_MISSING__'}:
                names[gene].add(int(col))
        source_map[source] = names
    params_file = root / 'NSCLC_context_frozen_parameters_20261016_v1/private_all53_frozen_gene_parameters.tsv'
    params = pd.read_csv(params_file, sep='\t')
    parents_file = root / 'NSCLC_automatic_ligand_response_review_20261016_v1/private_all_primary_BH_supported_rows.tsv'
    parents = pd.read_csv(parents_file, sep='\t').set_index('record')
    tf_scope_file = root / 'NSCLC_gene_disjoint_receiver_TF_preflight_20261016_v1/private_all_deduplicated_TF_slots.tsv'
    tf_scope = pd.read_csv(tf_scope_file, sep='\t')
    net_file = root / 'fullgene_receiver_target_code_v1/collectri_reference_v1/human_collectri_signed_complexes.tsv'
    net = pd.read_csv(net_file, sep='\t')
    targets_by_TF = {name: f for name, f in net.groupby('source', sort=False)}
    native_file = root / 'NSCLC_full_native_gene_units_20261016_v1/private_native_gene_axis.tsv'
    discovery_genes = pd.Index(pd.read_csv(native_file, sep='\t').gene)
    assert discovery_genes.is_unique and not net.duplicated(['source', 'target']).any()
    panel = sorted(set(PANEL) | set(tf_scope.TF) | set(parents.ligand) |
                   {g for rs in parents.qualified_receptors for rr in rs.split('|') for g in rr.split('_')})
    # One panel axis per physical column; unresolved aliases stay explicit below.
    panel_columns = sorted({next(iter(source_map[s][g])) for s in order for g in panel
                            if len(source_map[s][g]) == 1})
    panel_lookup = {col: i for i, col in enumerate(panel_columns)}
    out = root / 'existing_immune_fine_native_signatures_20261016_v2'
    assert not out.exists()
    remaining = resource(root)
    memory = int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines()
                      if line.startswith('MemAvailable:'))) * 1024
    assert memory > nu * n_features * 8 * 3 + 4 * 1024 ** 3
    out.mkdir(mode=0o700)
    # Metadata-only selection, fixed before inspecting expression or scores.
    priority = units.apply(lambda r: hashlib.sha256(('20261016|' + str(r.dataset) + '|' + str(r.library_id) + '|' + str(r.node_id)).encode()).hexdigest(), axis=1)
    audit_indices = set()
    for _, f in units.assign(priority=priority).groupby('dataset', sort=True):
        audit_indices.update(f.sort_values('priority').head(2).index)
    for _, f in units.assign(priority=priority).groupby('node_id', sort=True):
        audit_indices.add(f.sort_values('priority').index[0])
    audit_indices.update(units[units.n_cells.lt(20)].assign(priority=priority).sort_values('priority').head(2).index)
    for i in priority.sort_values().index:
        if len(audit_indices) >= 64:
            break
        audit_indices.add(i)
    assert len(audit_indices) == 64
    audit_units = units.loc[sorted(audit_indices)]
    audit_units.to_csv(out / 'private_prespecified64_raw_audit_units.tsv', sep='\t', index=False)
    plan = dict(created_UTC=now(), code_sha256=sha(__file__), scope_terminal_sha256=sha(scope / 'private_terminal.json'),
                external_audit_sha256=sha(external / 'private_independent_audit.json'), atlas_release_sha256=sha(atlas_run / 'manifests/immune_authority_release_public.json'),
                frozen_parameters_sha256=sha(params_file), frozen_parents_sha256=sha(parents_file), frozen_TF_scope_sha256=sha(tf_scope_file),
                signed_network_sha256=sha(net_file), discovery_gene_axis_sha256=sha(native_file),
                metadata_only64_raw_audit_unit_sha256=sha(out / 'private_prespecified64_raw_audit_units.tsv'),
                units=nu, selected_cells=int((mapping >= 0).sum()), original_nodes=29, source_count=len(order),
                normalization='Raw layers/counts, sum duplicate coordinates on working copy, remove stored zeros, full physical-native RNA denominator per cell, natural log1p(CP10K), then mean by ORIGINAL fine unit. X/ambient-corrected values never used.',
                scope='All frozen units including n<20; preserve every node. No endpoint promotion, pooling, clinical outcome, effect, P/BH or new download.',
                gene_matching='Unique physical source-native canonical/exact symbol or stable-ID pipe alias; ambiguous features and physical-column collisions unresolved. Forbidden programme/interface mappings exclude ALL matching physical columns.',
                programmes='Original per-gene centering/scaling and signed weights, sum only measurable unique genes without rescaling. Full sender/receiver and original receiver mask (all sender genes + ligand + all receptor components); original95% L1 and masked80%/20-gene gates.',
                regulators='Every original parent/TF slot retained. Same NSCLC disjoint target genes/signs and reference-native background gene names mapped to source physical columns; >=20 targets/500 background. Original discovery coverage failures remain unavailable. Record target/background retention explicitly; no expression-based filter. Native backgrounds can differ, therefore no automatic equality of score magnitudes or exact replication claim.',
                panel='Public lineage/identity markers, original literal interfaces and every original regulator transcript; source-unmeasured/unresolved distinct from observed zero.',
                memory='Full gene means only in RAM. Persist no complete gene-expression matrix; private score and limited panel arrays only.',
                output_cap_bytes=256 * 1024 ** 2, reserve_decimal_bytes=500000000000, starting_quota_remaining_bytes=remaining,
                prospective_independent_audit='Rebuild unit mapping from original obs, independent raw rows and all-gene means for fixed64 units across sources/nodes/low coverage, all programmes and all TF scores on audited units with official source algorithm; full matrix status/coverage and score invariants. Not an all-cell independent expression replay.',
                no_original_identity_or_CM_or_weight_changes=True)
    (out / 'plan_before_native_expression_measurement.json').write_text(json.dumps(plan, indent=2) + '\n')
    means = np.zeros((nu, n_features), dtype=np.float64)
    quality = np.zeros((nu, 4))
    panel_out = {name: np.zeros((nu, len(panel_columns)), dtype=np.int64 if name in ['positive', 'ge2'] else np.float64)
                 for name in ['sum_raw', 'sum_log', 'positive', 'ge2']}
    source_indices = np.asarray([order.index(s) for s in units.dataset], dtype=np.int8)
    block_size = 4096
    for start in range(0, n_cells, block_size):
        stop = min(n_cells, start + block_size)
        ptr = np.asarray(counts['indptr'][start:stop + 1], dtype=np.int64)
        x = csr_matrix((np.asarray(counts['data'][int(ptr[0]):int(ptr[-1])]),
                        np.asarray(counts['indices'][int(ptr[0]):int(ptr[-1])]), ptr - ptr[0]),
                       shape=(stop - start, n_features))
        selected = mapping[start:stop] >= 0
        if selected.any():
            ids = mapping[start:stop][selected]
            x = x[selected].copy()
            x.sum_duplicates(); x.eliminate_zeros(); x.sort_indices()
            assert measured[x.indices, np.repeat(source_indices[ids], np.diff(x.indptr))].all(), 'Counts outside source measurement mask'
            accumulate(x, ids, nu, means, quality, np.asarray(panel_columns), panel_out)
        if start // block_size % 32 == 0 or stop == n_cells:
            (out / 'private_checkpoint.json').write_text(json.dumps(dict(observed_UTC=now(), rows_read=stop, total_rows=n_cells,
                   selected_cells_accumulated=int(quality[:, 0].sum()), stage='raw_count_read_and_RAM_aggregation')) + '\n')
    assert np.array_equal(quality[:, 0], units.n_cells) and int(quality[:, 0].sum()) == plan['selected_cells']
    means /= quality[:, :1]
    assert np.isfinite(means).all()
    units['mean_full_RNA'] = quality[:, 1] / quality[:, 0]
    units['mean_logRNA'] = quality[:, 2] / quality[:, 0]
    units['mean_loggenes'] = quality[:, 3] / quality[:, 0]
    units['later_n20_eligible'] = units.n_cells.ge(20)
    units['within_original_major_fraction'] = units.n_cells / units.groupby(['dataset', 'library_id', 'major']).n_cells.transform('sum')
    units['exact_NSCLC_endpoint_mapping_validated'] = False
    definitions, programme_maps, programme_values, tf_maps, tf_definitions, tf_values, panel_maps = {}, [], [], [], {}, [], []
    for source in order:
        names = source_map[source]
        ui = units.index[units.dataset.eq(source)].to_numpy()
        lookup = lambda gene: next(iter(names[gene])) if len(names[gene]) == 1 else -1
        background_stats = {}
        for gene in panel:
            col = lookup(gene)
            panel_maps.append(dict(dataset=source, gene=gene, native_column=col, panel_column=panel_lookup.get(col, -1),
                                   source_native_unambiguous=col >= 0, physical_matches=len(names[gene])))
        for parent_id, parent in parents.iterrows():
            sender_type = 'typeB' if parent.receiver_type == 'typeA' else 'typeA'
            all_programme = params[params.context.eq(parent.context) & params.program.eq(parent.program)]
            sender_genes = set(all_programme.loc[all_programme.type.eq(sender_type), 'gene'])
            interfaces = {parent.ligand} | {g for rr in parent.qualified_receptors.split('|') for g in rr.split('_')}
            forbidden = set(all_programme.gene) | interfaces
            forbidden_columns = {col for gene in forbidden for col in names[gene]}
            for role, typ, masked in [('sender', sender_type, False), ('receiver', parent.receiver_type, False), ('receiver_masked', parent.receiver_type, True)]:
                pp = all_programme[all_programme.type.eq(typ)]
                original_l1 = float(pp.weight.abs().sum())
                if masked:
                    pp = pp[~pp.gene.isin(sender_genes | interfaces)]
                retained_l1 = float(pp.weight.abs().sum()) / original_l1
                cols = np.array([lookup(g) for g in pp.gene])
                collisions = Counter(cols[cols >= 0])
                keep = np.array([c >= 0 and collisions[c] == 1 for c in cols])
                measured_l1 = float(pp.weight.abs().to_numpy()[keep].sum() / pp.weight.abs().sum())
                status = 'MEASURED_FROZEN_RAW_GENE_SCORE_NOT_ENDPOINT_ACCEPTANCE'
                if measured_l1 < .95 or (masked and (retained_l1 < .8 or keep.sum() < 20)):
                    status = 'INSUFFICIENT_FROZEN_WEIGHT_OR_GENE_COVERAGE_NOT_NEGATIVE'
                key = (int(parent.context), typ, int(parent.program), tuple(pp.gene), tuple(cols[keep]), tuple(keep))
                if key not in definitions:
                    definitions[key] = len(programme_values)
                    values = np.full(nu, np.nan)
                    programme_values.append(values)
                j = definitions[key]
                if status.startswith('MEASURED_'):
                    active = pp.loc[keep]
                    programme_values[j][ui] = means[np.ix_(ui, cols[keep])] @ (active.weight / active.gene_sd).to_numpy() - float((active.weight * active.gene_mean / active.gene_sd).sum())
                programme_maps.append(dict(dataset=source, record=parent_id, role=role, context=parent.context, type=typ, program=parent.program,
                                           score_column=j, original_L1_retained_after_mask=retained_l1, source_measured_L1_fraction=measured_l1,
                                           measured_genes=int(keep.sum()), status=status, corrected_for_training_quality=False))
            bg_genes = [g for g in discovery_genes if g not in forbidden]
            bg_cols = [lookup(g) for g in bg_genes]
            cc = Counter(c for c in bg_cols if c >= 0)
            background = np.asarray(sorted({c for c in bg_cols if c >= 0 and cc[c] == 1 and c not in forbidden_columns}), dtype=int)
            background_key = hashlib.sha256(background.astype('<i8').tobytes()).hexdigest()
            bg_index = {col: i for i, col in enumerate(background)}
            if background_key not in background_stats and len(background) >= 500:
                expression = means[np.ix_(ui, background)]
                avg = expression.mean(axis=1)
                expression -= avg[:, None]
                ssx = np.einsum('ij,ij->i', expression, expression)
                background_stats[background_key] = (avg, ssx)
                del expression
            bg_gene_set = set(bg_genes)
            active_slots = tf_scope[tf_scope.record.eq(parent_id)]
            # Score all original slots; coverage failures remain explicit rows.
            for slot in active_slots.itertuples():
                f = targets_by_TF[slot.TF]
                f = f[f.target.isin(bg_gene_set)]
                cols = np.array([lookup(g) for g in f.target])
                ccs = Counter(cols[cols >= 0])
                keep = np.array([c in bg_index and ccs[c] == 1 for c in cols])
                status = 'MEASURED_DISJOINT_TARGET_PATTERN_NOT_PROTEIN_ACTIVITY'
                if slot.status != 'READY_FOR_INPUT_AUDIT_NO_ACTIVITY_COMPUTED':
                    status = 'ORIGINAL_DISCOVERY_TF_COVERAGE_FAILURE_RETAINED'
                elif keep.sum() < 20 or len(background) < 500:
                    status = 'INSUFFICIENT_SOURCE_TARGETS_OR_BACKGROUND_NOT_NEGATIVE'
                target_cols = cols[keep]
                weights = f.weight.to_numpy(float)[keep]
                key = (background_key, tuple(target_cols), tuple(weights))
                if key not in tf_definitions:
                    tf_definitions[key] = len(tf_values)
                    tf_values.append(np.full(nu, np.nan))
                j = tf_definitions[key]
                if status.startswith('MEASURED_') and np.isnan(tf_values[j][ui]).all():
                    avg, ssx = background_stats[background_key]
                    ssw = float(weights @ weights - weights.sum() ** 2 / len(background))
                    covariance = means[np.ix_(ui, target_cols)] @ weights - avg * weights.sum()
                    correlation = covariance / np.sqrt(ssx * ssw)
                    assert np.isfinite(correlation).all() and (np.abs(correlation) < 1).all()
                    tf_values[j][ui] = correlation * np.sqrt((len(background) - 2) / ((1 - correlation + 1e-16) * (1 + correlation + 1e-16)))
                tf_maps.append(dict(dataset=source, record=parent_id, TF=slot.TF, score_column=j, status=status,
                                    original_disjoint_targets=len(f), source_disjoint_targets=int(keep.sum()),
                                    target_retained_fraction=float(keep.sum() / len(f)) if len(f) else np.nan,
                                    original_background=len(bg_genes), source_disjoint_background=len(background),
                                    background_retained_fraction=len(background) / len(bg_genes), original_status=slot.status))
        (out / 'private_checkpoint.json').write_text(json.dumps(dict(observed_UTC=now(), stage='frozen_scores', completed_source=source,
                   programme_rows=len(programme_maps), TF_rows=len(tf_maps))) + '\n')
    resource(root)
    # No complete gene expression or learned target list is persisted here.
    units.to_csv(out / 'private_all_original_fine_units_and_quality.tsv', sep='\t', index=False)
    pd.DataFrame(programme_maps).to_csv(out / 'private_all_source_parent_programme_scope.tsv.gz', sep='\t', index=False, compression='gzip')
    pd.DataFrame(tf_maps).to_csv(out / 'private_all_source_parent_TF_scope.tsv.gz', sep='\t', index=False, compression='gzip')
    pd.DataFrame(panel_maps).to_csv(out / 'private_all_source_panel_gene_axis.tsv.gz', sep='\t', index=False, compression='gzip')
    np.savez_compressed(out / 'private_frozen_programme_and_TF_scores.npz', programme=np.column_stack(programme_values), TF=np.column_stack(tf_values))
    np.savez_compressed(out / 'private_limited_native_panel.npz', physical_columns=np.asarray(panel_columns),
                        mean_logCP10K=panel_out['sum_log'] / quality[:, :1], sum_raw_counts=panel_out['sum_raw'],
                        positive_cell_counts=panel_out['positive'], ge2_cell_counts=panel_out['ge2'])
    total = sum(p.stat().st_size for p in out.iterdir() if p.is_file())
    assert total < plan['output_cap_bytes'], 'Output budget exceeded; no authority changed'
    result = dict(status='COMPLETE_ALL_FROZEN_FINE_UNIT_NATIVE_MEASUREMENTS_INDEPENDENT_AUDIT_PENDING_NO_EFFECTS', finished_UTC=now(),
                  code_sha256=sha(__file__), units=nu, nodes=29, sources=len(order), selected_cells=int(quality[:, 0].sum()),
                  programme_scope_rows=len(programme_maps), TF_scope_rows=len(tf_maps), distinct_programme_score_columns=len(programme_values),
                  distinct_TF_score_columns=len(tf_values), private_bytes_before_terminal=total,
                  files={p.name: sha(p) for p in out.iterdir() if p.is_file()}, all_gene_matrix_persisted=False,
                  no_effects_P_BH_endpoint_promotion_or_new_download=True)
    (out / 'private_terminal.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
