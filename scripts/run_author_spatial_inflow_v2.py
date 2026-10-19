"""Fixed all-source spatial inflow follow-up; original fine endpoints stay frozen."""
import io
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from qualify_all14_author_spatial_lineages_v1 import T_LABELS, MPS_LABELS
from author_spatial_native_reader_v2 import load_native
from author_spatial_inflow_math_v1 import graph, bins, group_mean, describe, null, pvalue, seed, fixture
from NSCLC_immune_RNA_MCP_cohort_math_v1 import adjust

MODES = ['DEPTH_AREA', 'DEPTH_AREA_XY']
DRAWS = 1999


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root, 3 * 1024**2)
    qualification = root / 'all14_author_spatial_lineage_qualification_20261019_v1'
    assert sha(qualification / 'private_terminal.json') == 'e68ea3caf60fc05a09abd0072ed9155aa8e09b4d27a19b72fe82f010a5ec1c9a'
    qt = json.loads((qualification / 'private_terminal.json').read_text())
    for name, info in qt['files'].items():
        assert sha(qualification / name) == info['sha256']
    request = root / 'private_complementary_spatial_interface_scope_request_20261019_v1.json'
    assert sha(request) == 'f31d71f1efdc8450f346020733b0edea354a962eeec57bcc58328ba36bd890d0'
    targets = json.loads(request.read_text())['targets']
    assert len(targets) == 1
    target = targets[0]
    genes = list(dict.fromkeys(target['original_ligand_genes'] + target['original_receptor_genes']))
    lc = [genes.index(g) for g in target['original_ligand_genes']]
    rc = [genes.index(g) for g in target['original_receptor_genes']]
    inputs = pd.read_csv(qualification / 'private_complete14_author_lineage_and_assay_scope.tsv', sep='\t')
    selected = inputs[inputs.author_labels_available & inputs.original_resource_gene_proxy_measurable & inputs.selected_geometry_has_source_proof]
    assert len(inputs) == 14 and len(selected) == 3
    mf = root / 'immune_same_name_family_source_review_v1/spatial_native_input_repair_v2/analysis_input_manifest.tsv'
    assert sha(mf) == '6fda649c8b2bf4abe010ba8d53d50e495bae29ce8d30648675d9be9bf91dcf83'
    manifest = pd.read_csv(mf, sep='\t').set_index('sample_id')
    old = root / 'immune_same_name_family_source_review_v1/spatial_contact_native_preflight_v3'
    accepted = json.loads((old / 'independent_audit.json').read_text())
    assert sha(old / 'inventory.json') == accepted['producer_inventory_sha256']
    old_inventory = json.loads((old / 'inventory.json').read_text())
    originals = {v['sample_id']: v for v in old_inventory['inputs']}
    out = root / 'author_spatial_inflow_20261019_v2'
    assert not out.exists()
    out.mkdir(mode=0o700)
    save = lambda name, value: (out / name).write_text(json.dumps(value, indent=2) + '\n')
    save('plan_before_new_spatial_expression_or_scores.json', dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), module_sha256={n: sha(Path(__file__).with_name(n)) for n in ['author_spatial_inflow_math_v1.py', 'author_spatial_native_reader_v2.py', 'qualify_all14_author_spatial_lineages_v1.py']}, qualification_sha256=sha(qualification / 'private_terminal.json'), private_target_sha256=sha(request), all_source_inputs=14, all_source_cancer_strata=7, native_genes=genes, selected_by_metadata_not_gene_values=True, lineage_author_labels=dict(T=T_LABELS, MPS=MPS_LABELS), primary='Mean over author T cells of receptor times all-QC-neighbor-normalized ligand from author MPS. Primary log1p(1e4*max(count-FOV-negative-probe95,0)/full-native-panel-RNA).', geometry='Within each FOV all accepted QC cells, remove duplicate centroids, symmetric6NN capped at2*median nearest-neighbor distance. No self/crossFOV edges. Sensitivity12NN cap3. Absolute physical units unknown.', null='Independent gene-mark permutations within FOV and author lineage, RNA quartiles and area halves; second null adds lineage-specific medianX/medianY quadrants. Geometry and labels fixed, singleton strata fixed.1999draws each, two-sided add-one tails; jointP=max(two nullP). BH7 including failed P1.', aggregation='Equal weight eligible FOVs within each image, then equal image weight within GSE/cancer. Technical conditional association only, not independent patients or population inference.', sensitivity='Raw-logCP10K and wider graph: observed/analytic conditional-null means and four-term decomposition only; no new P or rescue of failed primary.', no_positive_expression_sample_filter=True, unavailable_fine_gamma_and_specified_assembly_preserved=True, full_expression_graph_or_cell_score_arrays_persisted=False, new_download=False, resource_guard=guard, output_cap_bytes=2 * 1024**2, reference='https://liana.readthedocs.io/en/stable/tutorials/notebooks/inflow_score.html', implementation_note='Explicit mathematical definition, not a claim that a different installed API was invoked.'))
    save('private_numeric_fixture.json', fixture())
    archive_checks = {}
    for source in selected.itertuples():
        archive = Path(originals[source.sample_id]['archive'])
        if str(archive) not in archive_checks:
            observed = sha(archive)
            assert observed == originals[source.sample_id]['archive_sha256']
            archive_checks[str(archive)] = observed
    fields, raw_proofs, image_results, excluded = [], [], [], []
    all_strata = sorted(set(zip(inputs.dataset_id, inputs.cancer)))
    for number, source in enumerate(selected.itertuples()):
        path = Path(manifest.loc[source.sample_id, 'analysis_path'])
        cache_rel = 'private_native_cell_caches/' + source.sample_id + '.npz'
        cache = old / cache_rel
        assert sha(path) == source.source_input_sha256 and sha(cache) == old_inventory['outputs'][cache_rel]
        data, proof = load_native(path, cache, Path(originals[source.sample_id]['archive']), genes)
        raw_proofs.append(dict(dataset_id=source.dataset_id, sample_id=source.sample_id, input_sha256=sha(path), **proof))
        labels, qc, fov = data['author_labels'].astype(str), data['qc'].astype(bool), data['fov'].astype(str)
        total, area, xy = data['native_total'].astype(float), data['area'].astype(float), data['coords'].astype(float)
        counts = data['counts'].astype(float)
        excess = np.maximum(counts - data['background_q95'][:, None], 0.)
        expressions = {name: np.log1p(np.divide(1e4 * value, total[:, None], out=np.zeros_like(value), where=total[:, None] > 0)) for name, value in [('BACKGROUND_EXCESS', excess), ('RAW', counts)]}
        null_sum = {mode: np.zeros(DRAWS) for mode in MODES}
        expected_sum = {mode: 0. for mode in MODES}
        image_fields = []
        for field in sorted(set(fov)):
            mask = (fov == field) & qc & np.isfinite(xy).all(1)
            ix = np.flatnonzero(mask)
            _, inv, multiplicity = np.unique(xy[ix], axis=0, return_inverse=True, return_counts=True)
            duplicate = multiplicity[inv] > 1
            ix = ix[~duplicate]
            T = np.isin(labels[ix], T_LABELS)
            M = np.isin(labels[ix], MPS_LABELS)
            if min(int(T.sum()), int(M.sum())) < 20:
                excluded.append(dict(sample_id=source.sample_id, FOV=field, reason='FEWER_THAN20_PER_AUTHOR_LINEAGE_AFTER_GEOMETRY_QC', T_cells=int(T.sum()), MPS_cells=int(M.sum())))
                continue
            ti, mi = np.flatnonzero(T), np.flatnonzero(M)
            original = expressions['BACKGROUND_EXCESS'][ix]
            L, R = original[M][:, lc].min(1), original[T][:, rc].min(1)
            g, geom = graph(xy[ix], 6, 2.)
            cross = g[ti][:, mi].tocsr()
            base = dict(dataset_id=source.dataset_id, sample_id=source.sample_id, cancer=source.cancer, FOV=field, all_graph_QC_cells=len(ix), T_cells=len(ti), MPS_cells=len(mi), duplicate_centroids_removed=int(duplicate.sum()), ligand_positive_MPS=int((L > 0).sum()), receptor_positive_T=int((R > 0).sum()), **geom)
            primary = describe(cross, L, R)
            image_fields.append(primary)
            for mode in MODES:
                lg = bins(total[ix][M], area[ix][M], xy[ix][M], mode)
                rg = bins(total[ix][T], area[ix][T], xy[ix][T], mode)
                values, expected = null(cross, L, R, lg, rg, DRAWS, seed(source.sample_id, field, mode))
                null_sum[mode] += values
                expected_sum[mode] += expected
                fields.append(dict(base, readout='BACKGROUND_EXCESS', scale='PRIMARY6NN_CAP2', null=mode, analytic_null_expected=expected, movable_MPS_cells=sum(len(v) for v in lg if len(v) > 1), movable_T_cells=sum(len(v) for v in rg if len(v) > 1), **primary))
            for readout, k, radius in [('RAW', 6, 2.), ('BACKGROUND_EXCESS', 12, 3.)]:
                if k == 6:
                    W = cross
                else:
                    whole, _ = graph(xy[ix], k, radius)
                    W = whole[ti][:, mi].tocsr()
                e = expressions[readout][ix]
                ll, rr = e[M][:, lc].min(1), e[T][:, rc].min(1)
                desc = describe(W, ll, rr)
                for mode in MODES:
                    lg = bins(total[ix][M], area[ix][M], xy[ix][M], mode)
                    rg = bins(total[ix][T], area[ix][T], xy[ix][T], mode)
                    expected = float(np.mean(group_mean(rr, rg) * (W @ group_mean(ll, lg))))
                    fields.append(dict(base, readout=readout, scale='PRIMARY6NN_CAP2' if k == 6 else 'WIDE12NN_CAP3', null=mode, analytic_null_expected=expected, movable_MPS_cells=np.nan, movable_T_cells=np.nan, **desc))
            save('private_progress.json', dict(stage='FIXED_FIELD_CONDITIONAL_INFLOW', completed_images=number, current_image_eligible_fields=len(image_fields), completed_field_readouts=len(fields), observed_UTC=datetime.now(timezone.utc).isoformat()))
        n = len(image_fields)
        if n:
            observed = float(np.mean([v['inflow'] for v in image_fields]))
            image_results.append(dict(dataset_id=source.dataset_id, sample_id=source.sample_id, cancer=source.cancer, eligible_FOVs=n, observed=observed, null_values={mode: null_sum[mode] / n for mode in MODES}, analytic_expectation={mode: expected_sum[mode] / n for mode in MODES}))
        save('private_progress.json', dict(stage='FIXED_IMAGE_COMPLETE', completed_images=number + 1, total_images=len(selected), observed_UTC=datetime.now(timezone.utc).isoformat()))
    rows = []
    streams = {}
    for index, (gse, cancer) in enumerate(all_strata):
        images = [v for v in image_results if (v['dataset_id'], v['cancer']) == (gse, cancer)]
        row = dict(dataset_id=gse, cancer=cancer, original_source_inputs=int(((inputs.dataset_id == gse) & (inputs.cancer == cancer)).sum()), measured_images=len(images), eligible_technical_FOVs=sum(v['eligible_FOVs'] for v in images), reviewed_independent_patients=0, status='PREREQUISITES_UNAVAILABLE_NOT_BIOLOGICAL_NEGATIVE', observed=np.nan, joint_P=1.)
        if images:
            row.update(status='CONDITIONAL_IMAGE_ASSOCIATION_ONLY', observed=float(np.mean([v['observed'] for v in images])))
            probabilities = []
            for mode in MODES:
                values = np.mean([v['null_values'][mode] for v in images], axis=0)
                expected = float(np.mean([v['analytic_expectation'][mode] for v in images]))
                probability = pvalue(row['observed'], values)
                if np.std(values) <= 1e-14:
                    probability = 1.
                    row['status'] = 'MEASURED_DEGENERATE_OR_UNINFORMATIVE_NULL'
                row[mode + '_expected'] = expected
                row[mode + '_null_SD'] = float(np.std(values, ddof=1))
                row[mode + '_observed_minus_expected'] = row['observed'] - expected
                row[mode + '_P'] = probability
                probabilities.append(probability)
                streams['stratum%02d_%s' % (index, mode)] = values
            row['joint_P'] = max(probabilities)
        rows.append(row)
    result = pd.DataFrame(rows)
    assert len(result) == 7
    result['BH7'] = adjust(result.joint_P.to_numpy())
    result['BY7'] = adjust(result.joint_P.to_numpy(), by=True)
    result.to_csv(out / 'private_complete7_conditional_source_inflow_results.tsv', sep='\t', index=False)
    pd.DataFrame(fields).to_csv(out / 'private_all_field_inflow_and_decompositions.tsv.gz', sep='\t', index=False, compression='gzip')
    pd.DataFrame(excluded).to_csv(out / 'private_all_field_exclusions.tsv.gz', sep='\t', index=False, compression='gzip')
    save('private_all_original_native_panel_and_metadata_proofs.json', raw_proofs)
    save('private_unique_source_archive_SHA256.json', archive_checks)
    buf = io.BytesIO()
    np.savez_compressed(buf, **streams)
    assert len(buf.getvalue()) < 512 * 1024
    (out / 'private_complete_source_aggregate_null_streams.npz').write_bytes(buf.getvalue())
    receipt = dict(status='COMPLETE_FIXED7_AUTHOR_SPATIAL_INFLOW_FULL_NATIVE_PARITY_AND_CONDITIONAL_NULLS_REVIEW_PENDING', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), qualification_sha256=sha(qualification / 'private_terminal.json'), resource_guard=guard, source_stratum_positions=7, selected_input_positions=3, conditional_draws_per_null=DRAWS, full_expression_graph_or_cell_scores_persisted=False, fine_gamma_or_full_specified_assembly_validated=False, original_annotations_or_old_inference_changed=False, scientific_outputs_public=False, files={p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()})
    assert sum(v['bytes'] for v in receipt['files'].values()) < 2 * 1024**2
    save('private_terminal.json', receipt)
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
