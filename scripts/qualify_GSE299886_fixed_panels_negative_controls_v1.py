"""All-source fixed-panel capture and original negative-probe QC; no effects.

Inputs and coordinates remain private; legacy matrices and labels are immutable.
"""
import gzip
import json
import os
import tarfile
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from diagnose_GSE299886_original_metadata_identity_v1 import strings
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from repair_GSE299886_original_coordinate_sidecars_v1 import dense_rows


def now():
    return datetime.now(timezone.utc).isoformat()


def verify(directory):
    t = json.loads((directory / 'private_terminal.json').read_text())
    for n, v in t['files'].items():
        assert sha(directory / n) == v['sha256']


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    budget = resource(root, 96 * 1024**2)
    coord = root / 'GSE299886_original_coordinate_sidecars_20261018_v1'
    verify(coord)
    reports = json.loads((coord / 'private_all_source_coordinate_reconstruction_receipt.json').read_text())
    assert len(reports) == 3
    headers = json.loads((root / 'GSE299886_original_metadata_identity_diagnostic_20261018_v1/private_original_source_member_headers.json').read_text())
    schemas = json.loads((root / 'CM02_existing_targeted_spatial_schema_20261018_v2/private_all_eligible_targeted_input_schema.json').read_text())
    state = Path(json.loads((root / 'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root'])
    ap = state / 'state_spatial_proxy_calibration_inputs_20261016_v1/fixed_all14_panel_gene_axes.json'
    axes = json.loads(ap.read_text())
    fp = root.parents[1] / 'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json'
    frozen = json.loads(fp.read_text())
    out = root / 'GSE299886_fixed_panels_negative_controls_20261018_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    plan = dict(created_UTC=now(), source_sha256=sha(__file__), budget=budget,
                coordinate_terminal_sha256=sha(coord / 'private_terminal.json'),
                frozen_panel_contract_sha256=sha(fp), axes_sha256=sha(ap),
                scope='All3 repaired assays, all7 unchanged panels, all omitted source negative probes. No response, survival, spatial mixing effects, P or label changes.',
                lineage_strata=['ALL', 'pan_T_any_gene', 'pan_T_at_least2_distinct_genes'],
                lineage_caveat='Technical RNA capture strata, not independently validated T identity or exact CM subtype.',
                controls='Require all omitted source features negative/control-like and exact original RNA ID order. Compare retained biological sum + all negative probes against deposited metadata total; retain every discrepancy.',
                geometry='Nearest-other-cell distance distributions only, no absolute-unit assumption or patient/core inference. FOV missing remains missing.')
    (out / 'plan_before_all_source_capture_and_negative_probe_QC.json').write_text(json.dumps(plan, indent=2) + '\n')
    capture = []
    receipt = []
    for number, report in enumerate(reports):
        sid = report['sample_id']
        src = next(x for x in headers if x['sample_id'] == sid)
        schema = next(x for x in schemas if x['sample_id'] == sid)
        axis = next(x for x in axes if x['sample_id'] == sid)
        matrix = Path(schema['analysis_file'])
        assert sha(matrix) == report['matrix_SHA256'] == axis['input_sha256']
        archive = Path(src['archive'])
        before = (archive.stat().st_size, archive.stat().st_mtime_ns)
        z = np.load(coord / report['output'], allow_pickle=False)
        ids = z['source_cell_id']
        n = len(ids)
        with h5py.File(matrix, 'r') as h:
            features = strings(h['var'][h['var'].attrs['_index']][...]).tolist()
            if not set(features) <= set(src['count_header']):
                from anndata._io.specs import read_elem
                features = strings(read_elem(h['var']['gene_symbol'])).tolist()
            extras = [x for x in src['count_header'] if x not in set(features) and x != 'cell']
            assert len(features) == 500 and len(extras) == 50
            assert all(any(s in x.lower() for s in ['blank', 'control', 'negative']) for x in extras)
            negative_total = np.zeros(n, np.int64)
            negative_detected = np.zeros(n, np.int16)
            processed = 0
            with tarfile.open(archive, 'r:') as tar:
                f = tar.extractfile(src['count_member'])
                stream = gzip.GzipFile(fileobj=f) if src['count_member'].endswith('.gz') else f
                for chunk in pd.read_csv(stream, usecols=['cell', *extras], dtype={'cell': str}, keep_default_na=False, chunksize=4096):
                    stop = processed + len(chunk)
                    assert np.array_equal(chunk['cell'].to_numpy(str), ids[processed:stop])
                    v = chunk[extras].to_numpy(float)
                    assert np.isfinite(v).all() and (v >= 0).all() and (v == np.floor(v)).all()
                    negative_total[processed:stop] = v.sum(axis=1).astype(np.int64)
                    negative_detected[processed:stop] = (v > 0).sum(axis=1)
                    processed = stop
            assert processed == n
            mapping = axis['symbol_to_feature_index']
            panels = frozen['panels']
            complete = {k: all(g in mapping for g in genes) for k, genes in panels.items()}
            assert complete['pan_T']
            detected = {k: np.zeros(n, np.int16) for k in panels if complete[k]}
            twice = {k: np.zeros(n, np.int16) for k in detected}
            total = np.zeros(n, np.int64)
            for start in range(0, n, 2048):
                stop = min(n, start + 2048)
                v = dense_rows(h['X'], start, stop, len(features))
                total[start:stop] = v.sum(axis=1, dtype=np.int64)
                for k in detected:
                    a = v[:, [mapping[g] for g in panels[k]]]
                    detected[k][start:stop] = (a > 0).sum(axis=1)
                    twice[k][start:stop] = (a >= 2).sum(axis=1)
            assert np.array_equal(total, z['deposited_panel_total'])
            strata = {'ALL': np.ones(n, bool), 'pan_T_any_gene': detected['pan_T'] >= 1,
                      'pan_T_at_least2_distinct_genes': detected['pan_T'] >= 2}
            for sk, mask in strata.items():
                for k, genes in panels.items():
                    capture.append(dict(sample_id=sid, stratum=sk, n=int(mask.sum()), panel=k,
                        required_genes=len(genes), measured_genes=sum(g in mapping for g in genes), complete=complete[k],
                        cells_any_gene=int((mask & (detected[k] >= 1)).sum()) if complete[k] else None,
                        cells_ge2_genes=int((mask & (detected[k] >= 2)).sum()) if complete[k] else None,
                        cells_ge2_genes_each_value_ge2=int((mask & (twice[k] >= 2)).sum()) if complete[k] else None,
                        median_biological_total=float(np.median(total[mask])) if mask.any() else None,
                        median_negative_total=float(np.median(negative_total[mask])) if mask.any() else None))
            known = np.isfinite(z['source_metadata_transcript_count'])
            discrepancy = z['source_metadata_transcript_count'] - total - negative_total
            nn = cKDTree(z['xy_source']).query(z['xy_source'], k=2, workers=1)[0][:, 1]
            np.savez_compressed(out / ('private_source_%02d_negative_QC.npz' % number),
                existing_row_index=np.arange(n), negative_total=negative_total, negative_detected=negative_detected,
                biological_total=total, metadata_minus_biological_minus_negative=discrepancy,
                existing_matrix_sha256=np.array(report['matrix_SHA256']))
            receipt.append(dict(sample_id=sid, rows=n, negative_features=len(extras),
                source_RNA_ID_order_exact=True, retained_biological_total_rechecked=True,
                metadata_total_known=int(known.sum()), metadata_total_mismatches_after_negative=int((known & (discrepancy != 0)).sum()),
                metadata_total_discrepancy_range=[float(discrepancy[known].min()), float(discrepancy[known].max())] if known.any() else None,
                cells_with_negative_probe=int((negative_total > 0).sum()), negative_transcripts=int(negative_total.sum()),
                biological_transcripts=int(total.sum()), negative_per_feature_to_biological_per_feature_ratio=float((negative_total.sum()/len(extras))/(total.sum()/len(features))) if total.sum() else None,
                nearest_other_cell_source_unit_quantiles=np.quantile(nn, [0, .01, .1, .5, .9, .99, 1]).tolist(),
                FOV_known_rows=int(z['FOV_known'].sum()), physical_unit='UNVERIFIED', patient_core='UNREVIEWED'))
        assert before == (archive.stat().st_size, archive.stat().st_mtime_ns)
        assert sha(matrix) == report['matrix_SHA256']
    pd.DataFrame(capture).to_csv(out / 'private_all63_fixed_panel_capture_strata.tsv', sep='\t', index=False)
    (out / 'private_all3_negative_probe_geometry_QC.json').write_text(json.dumps(receipt, indent=2) + '\n')
    files = [p for p in out.iterdir() if p.is_file()]
    assert sum(p.stat().st_size for p in files) < 96 * 1024**2
    terminal = dict(status='COMPLETE_ALL3_NEGATIVE_PROBE_AND_ALL_FIXED_PANEL_CAPTURE_NO_EFFECTS',
        completed_UTC=now(), source_sha256=sha(__file__), capture_positions=len(capture),
        files={p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in files},
        no_spatial_test_or_state_assignment=True)
    (out / 'private_terminal.json').write_text(json.dumps(terminal, indent=2) + '\n')
    print(json.dumps(dict(status=terminal['status'], completed_UTC=terminal['completed_UTC'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
