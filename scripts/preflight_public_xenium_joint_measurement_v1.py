"""RAM-only panel and cell-key preflight of all selected Xenium matrices/cells.

Reads processed public matrices and coordinate tables, no images/transcript or
TCR-specific analysis. Exports only source-level panels and technical summaries.
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
import gzip
import hashlib
import io
import json
import os
import re
import traceback
import urllib.parse
import urllib.request

import h5py
import numpy as np
import pandas as pd
import pyarrow
import pyarrow.parquet as pq

GENES = ['CD3D','CD3E','CD3G','TRAC','CD4','CD8A','CD8B','FOXP3','IL2RA','CTLA4','IKZF2','TNFRSF9','LAG3','CXCL13','PDCD1','TOX','HAVCR2','ENTPD1','TCF7','CCR7','IL7R','GZMB','PRF1','NKG7','GNLY','MKI67','TOP2A','STMN1','TYMS','UBE2C','PCNA','EPCAM','KRT5','KRT14','KRT17','KRT8','KRT18','KRT19','CLDN4','TP63','SOX2','LST1','LYZ','CD68','C1QA','MS4A1','CD79A','COL1A1','DCN','PECAM1','VWF']
CELL_FIELDS = ['cell_id','x_centroid','y_centroid','transcript_counts','control_probe_counts','control_codeword_counts','total_counts','cell_area','nucleus_area']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def decode(values):
    return [x.decode('utf-8') if isinstance(x, bytes) else str(x) for x in values]


def load_index():
    path = Path(os.environ['PANCANCER_TAR_INDEX'])
    assert sha(path) == os.environ['PANCANCER_TAR_INDEX_SHA256']
    index = json.loads(path.read_text())
    assert index['complete_two_zero_block_terminator']
    assert urllib.parse.urlparse(index['url']).hostname == 'ftp.ncbi.nlm.nih.gov'
    selected = {}
    for r in index['members']:
        if r['role'] in ['cell_feature_matrix', 'cell_coordinates_and_metrics']:
            key = re.match(r'GSM[0-9]+', r['name'])[0]
            assert r['role'] not in selected.setdefault(key, {})
            selected[key][r['role']] = r
    assert len(selected) == 18 and all(len(v) == 2 for v in selected.values())
    assert sum(r['bytes'] for v in selected.values() for r in v.values()) <= 160*1024**2
    return index, [selected[k] for k in sorted(selected)]


def read_member(index, member):
    start, n = member['payload_offset'], member['bytes']
    assert 0 < n <= 64*1024**2
    headers = {'Range': f'bytes={start}-{start+n-1}', 'Accept-Encoding': 'identity'}
    if index.get('source_ETag'):
        headers['If-Match'] = index['source_ETag']
    else:
        assert index.get('source_Last_Modified')
        headers['If-Unmodified-Since'] = index['source_Last_Modified']
    with urllib.request.urlopen(urllib.request.Request(index['url'], headers=headers), timeout=60) as response:
        assert response.status == 206
        assert response.headers.get('Content-Range') == f'bytes {start}-{start+n-1}/{index["archive_bytes"]}'
        assert response.headers.get('ETag') == index.get('source_ETag')
        assert response.headers.get('Last-Modified') == index.get('source_Last_Modified')
        raw = response.read(n+1)
    assert len(raw) == n
    return raw


def inspect_pair(index, pair, ordinal):
    matrix_bytes = read_member(index, pair['cell_feature_matrix'])
    coordinate_bytes = read_member(index, pair['cell_coordinates_and_metrics'])
    with h5py.File(io.BytesIO(matrix_bytes), 'r') as h:
        matrix = h['matrix']
        names = decode(matrix['features/name'][:])
        kinds = decode(matrix['features/feature_type'][:])
        cells = decode(matrix['barcodes'][:])
        shape = matrix['shape'][:].astype(int).tolist()
        assert len(names) == len(kinds) == shape[0] and len(cells) == shape[1]
        assert pd.Index(cells).is_unique and len(matrix['indptr']) == len(cells)+1
        panel = {g: sum(n == g and k == 'Gene Expression' for n, k in zip(names, kinds)) for g in GENES}
        assert all(n <= 1 for n in panel.values())
    raw = coordinate_bytes
    if raw[:2] == b'\x1f\x8b':
        with gzip.GzipFile(fileobj=io.BytesIO(raw)) as f:
            raw = f.read(512*1024**2+1)
        assert len(raw) <= 512*1024**2
    assert raw[:4] == raw[-4:] == b'PAR1'
    source = io.BytesIO(raw)
    schema = pq.read_schema(source)
    assert {'cell_id', 'x_centroid', 'y_centroid'}.issubset(schema.names)
    coordinates = pq.read_table(io.BytesIO(raw), columns=[c for c in CELL_FIELDS if c in schema.names]).to_pandas()
    cell_ids = coordinates.cell_id.astype(str)
    assert cell_ids.is_unique
    positions = pd.Index(cell_ids).get_indexer(cells)
    aligned = bool(len(coordinates) == len(cells) and (positions >= 0).all())
    finite = bool(np.isfinite(coordinates[['x_centroid', 'y_centroid']].to_numpy(float)).all())
    return dict(archive_sample_ordinal=ordinal, status='ALIGNED_PROCESSED_MATRIX_AND_COORDINATE_KEYS' if aligned and finite else 'UNRESOLVED_CELL_KEY_OR_COORDINATE', matrix_SHA256=hashlib.sha256(matrix_bytes).hexdigest(), compressed_coordinate_SHA256=hashlib.sha256(coordinate_bytes).hexdigest(), processed_bytes=len(matrix_bytes)+len(coordinate_bytes), matrix_features=shape[0], matrix_cells=shape[1], coordinate_cells=len(coordinates), exact_cell_key_alignment=aligned, finite_coordinates=finite, coordinate_schema={f.name: str(f.type) for f in schema}, whitelisted_gene_copies=panel, expression_values_not_analyzed=True, cell_rows_and_unlisted_feature_names_not_exported=True)


def main():
    out = Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert out.parent.parent == Path('/dev/shm') and not out.exists()
    out.mkdir(mode=0o700)
    save = lambda n, v: (out / n).write_text(json.dumps(v, indent=2) + '\n')
    save('plan_before_joint_measurement_preflight.json', dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_SHA256=sha(__file__), archive_index_SHA256=os.environ['PANCANCER_TAR_INDEX_SHA256'], selected_roles=['cell_feature_matrix', 'cell_coordinates_and_metrics'], all18_archive_samples_including_non_HNSCC_preserved_for_qualification=True, biological_exclusion_before_any_future_spatial_effect=True, known_gene_whitelist=GENES, coordinate_fields=CELL_FIELDS, no_images_full_transcripts_or_TCR_analysis=True, all_payloads_and_individual_arrays_RAM_only=True, no_new_spatial_or_clinical_effects=True, HTTP_range_bytes_cap=160*1024**2, runtime=dict(h5py=h5py.__version__, numpy=np.__version__, pandas=pd.__version__, pyarrow=pyarrow.__version__)))
    index, pairs = load_index()
    rows = []
    with ThreadPoolExecutor(max_workers=3) as pool:
        work = {pool.submit(inspect_pair, index, p, i): i for i, p in enumerate(pairs)}
        for future in as_completed(work):
            try:
                rows.append(future.result())
            except Exception as error:
                rows.append(dict(archive_sample_ordinal=work[future], status='FAILED_QUALIFICATION', error_class=type(error).__name__, source_lines=[f.lineno for f in traceback.extract_tb(error.__traceback__) if f.filename == __file__]))
            checkpoint = out / 'private_checkpoint.tmp'
            checkpoint.write_text(json.dumps(dict(completed_opportunities=len(rows), failed_opportunities=sum(r['status'] == 'FAILED_QUALIFICATION' for r in rows)), indent=2)+'\n')
            checkpoint.replace(out / 'private_checkpoint.json')
    rows.sort(key=lambda r: r['archive_sample_ordinal'])
    save('private_all18_joint_measurement_opportunities.json', rows)
    succeeded = [r for r in rows if r['status'] == 'ALIGNED_PROCESSED_MATRIX_AND_COORDINATE_KEYS']
    save('private_whitelist_gene_panel_coverage.json', [dict(gene=g, aligned_sample_files_with_one_copy=sum(r['whitelisted_gene_copies'][g] == 1 for r in succeeded), all_archive_samples=len(rows)) for g in GENES])
    terminal = dict(status='COMPLETE_ALL18_JOINT_MEASUREMENT_OPPORTUNITIES', completed_UTC=datetime.now(timezone.utc).isoformat(), source_SHA256=sha(__file__), all18_opportunities_preserved=len(rows) == 18, all18_aligned=len(succeeded) == 18, scientific_cohort_and_state_identity_not_yet_accepted=True, no_cell_rows_or_unlisted_probe_names_saved=True, no_new_spatial_or_clinical_effects=True, files={p.name: dict(SHA256=sha(p), bytes=p.stat().st_size) for p in out.iterdir()})
    save('private_terminal.json', terminal)
    print(json.dumps(dict(status=terminal['status'], all18_aligned=terminal['all18_aligned'])))


if __name__ == '__main__':
    main()
