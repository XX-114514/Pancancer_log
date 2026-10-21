"""Independent HDF5/Arrow key audit with bounded, resumable-size HTTP ranges.

All source payloads and cell arrays stay in RAM. Rechecks official CXCL13 aliases
and stable Ensembl ID, never replaces an unmeasured gene with a prediction.
"""
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from threading import Lock
import gzip
import hashlib
import importlib.util
import io
import json
import math
import os
import traceback
import urllib.error
import urllib.request

import h5py
import pyarrow.parquet as pq

CXCL13_ALIASES = {'CXCL13','ANGIE','ANGIE2','BCA-1','BCA1','BLC','BLR1L','SCYB13'}
CXCL13_ENSEMBL = 'ENSG00000156234'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parent = Path(os.environ['PANCANCER_JOINT_PREFLIGHT'])
    out = Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert out.parent.parent == Path('/dev/shm') and not out.exists()
    t = json.loads((parent / 'private_terminal.json').read_text())
    for n, v in t['files'].items():
        assert sha(parent / n) == v['SHA256']
    source = parent.parent / 'preflight_public_xenium_joint_measurement_v1.py'
    assert sha(source) == t['source_SHA256']
    spec = importlib.util.spec_from_file_location('pinned_original_preflight', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    index, pairs = module.load_index()
    previous = json.loads((parent / 'private_all18_joint_measurement_opportunities.json').read_text())
    previous = {r['archive_sample_ordinal']: r for r in previous}
    failures = {i for i, r in previous.items() if r['status'] != 'ALIGNED_PROCESSED_MATRIX_AND_COORDINATE_KEYS'}
    assert all(previous[i]['status'] == 'FAILED_QUALIFICATION' and previous[i]['error_class'] == 'AssertionError' and 69 in previous[i]['source_lines'] for i in failures)
    out.mkdir(mode=0o700)
    save = lambda n, v: (out / n).write_text(json.dumps(v, indent=2) + '\n')
    save('plan_before_joint_audit_and_transport_recovery.json', dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_SHA256=sha(__file__), parent_terminal_SHA256=sha(parent / 'private_terminal.json'), all18_opportunities_preserved=True, read_chunk_bytes=4*1024**2, maximum_attempts_per_exact_chunk=3, retry_only_network_or_short_range_read=True, source_total_and_modification_fixed=True, all_old_success_payload_SHAs_rechecked=True, independent_Counter_cell_keys_and_Python_finite_coordinates=True, original_failed_cases_replayed_from_same_complete_RAM_bytes=True, CXCL13_aliases=sorted(CXCL13_ALIASES), CXCL13_Ensembl=CXCL13_ENSEMBL, identity_source='https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000156234', no_expression_values_or_cell_rows_exported=True))
    lock, network = Lock(), dict(bytes=0, requests=0)
    def read(member):
        start, n = member['payload_offset'], member['bytes']
        assert 0 < n <= 64*1024**2
        chunks = []
        for begin in range(start, start+n, 4*1024**2):
            size = min(4*1024**2, start+n-begin)
            for attempt in range(3):
                headers = {'Range': f'bytes={begin}-{begin+size-1}', 'Accept-Encoding': 'identity'}
                headers.update({'If-Match': index['source_ETag']} if index.get('source_ETag') else {'If-Unmodified-Since': index['source_Last_Modified']})
                try:
                    with urllib.request.urlopen(urllib.request.Request(index['url'], headers=headers), timeout=45) as r:
                        assert r.status == 206 and r.headers.get('Content-Range') == f'bytes {begin}-{begin+size-1}/{index["archive_bytes"]}'
                        assert r.headers.get('ETag') == index.get('source_ETag') and r.headers.get('Last-Modified') == index.get('source_Last_Modified')
                        raw = r.read(size+1)
                    with lock:
                        network['requests'] += 1
                        network['bytes'] += len(raw)
                        assert network['bytes'] <= 320*1024**2 and network['requests'] <= 1024
                    if len(raw) != size:
                        raise EOFError('Short exact-range transport response')
                    chunks.append(raw)
                    break
                except (urllib.error.URLError, TimeoutError, EOFError):
                    if attempt == 2:
                        raise
        raw = b''.join(chunks)
        assert len(raw) == n
        return raw
    def inspect(ordinal, pair):
        matrix = read(pair['cell_feature_matrix'])
        coordinates = read(pair['cell_coordinates_and_metrics'])
        hashes = dict(matrix_SHA256=hashlib.sha256(matrix).hexdigest(), compressed_coordinate_SHA256=hashlib.sha256(coordinates).hexdigest())
        if ordinal not in failures:
            assert all(previous[ordinal][k] == v for k, v in hashes.items())
        text = lambda x: x.decode('utf-8') if isinstance(x, bytes) else str(x)
        with h5py.File(io.BytesIO(matrix), 'r') as h:
            m = h['matrix']
            names = [text(v) for v in m['features/name'][:]]
            ids = [text(v).split('.')[0] for v in m['features/id'][:]]
            types = [text(v) for v in m['features/feature_type'][:]]
            cells = [text(v) for v in m['barcodes'][:]]
            shape = [int(v) for v in m['shape'][:]]
            assert len(names) == len(ids) == len(types) == shape[0] and len(cells) == shape[1]
            counts = Counter(n for n, k in zip(names, types) if k == 'Gene Expression')
            panel = {g: counts[g] for g in module.GENES}
            assert all(v <= 1 for v in panel.values())
            alias = sum(k == 'Gene Expression' and n.strip().upper() in CXCL13_ALIASES for n, k in zip(names, types))
            ensembl = sum(k == 'Gene Expression' and i == CXCL13_ENSEMBL for i, k in zip(ids, types))
            prefix = sum(k == 'Gene Expression' and n.strip().upper().startswith('CXCL13') for n, k in zip(names, types))
        raw = coordinates
        if raw[:2] == b'\x1f\x8b':
            with gzip.GzipFile(fileobj=io.BytesIO(raw)) as f:
                raw = f.read(512*1024**2+1)
            assert len(raw) <= 512*1024**2
        assert raw[:4] == raw[-4:] == b'PAR1'
        schema = pq.read_schema(io.BytesIO(raw))
        table = pq.read_table(io.BytesIO(raw), columns=['cell_id','x_centroid','y_centroid'])
        cc = Counter(str(v) for v in table.column('cell_id').to_pylist())
        bc = Counter(cells)
        assert cc == bc and all(v == 1 for v in bc.values())
        finite = all(math.isfinite(float(v)) for field in ['x_centroid','y_centroid'] for v in table.column(field).to_pylist())
        assert finite
        if ordinal not in failures:
            assert previous[ordinal]['whitelisted_gene_copies'] == panel
            assert previous[ordinal]['matrix_features'] == shape[0] and previous[ordinal]['matrix_cells'] == shape[1] == previous[ordinal]['coordinate_cells']
            assert previous[ordinal]['coordinate_schema'] == {f.name: str(f.type) for f in schema}
        result = dict(archive_sample_ordinal=ordinal, status='ACCEPTED_EXACT_MATRIX_COORDINATE_CELL_KEY_MULTISET', **hashes, matrix_features=shape[0], matrix_cells=shape[1], coordinate_cells=len(cc), exact_cell_key_alignment=True, finite_coordinates=True, coordinate_schema={f.name: str(f.type) for f in schema}, feature_type_counts=dict(Counter(types)), whitelisted_gene_copies=panel, CXCL13_known_alias_matches=alias, CXCL13_Ensembl_id_matches=ensembl, CXCL13_prefix_name_candidates=prefix, all_original_success_checks_replayed=ordinal not in failures)
        if ordinal in failures:
            result['_RAM_recovery_payloads'] = (matrix, coordinates)
        return result
    rows = []
    with ThreadPoolExecutor(max_workers=3) as pool:
        futures = {pool.submit(inspect, i, p): i for i, p in enumerate(pairs)}
        for f in as_completed(futures):
            rows.append(f.result())
            temporary = out / 'private_checkpoint.tmp'
            temporary.write_text(json.dumps(dict(completed_opportunities=len(rows), network=dict(network)), indent=2)+'\n')
            temporary.replace(out / 'private_checkpoint.json')
    for row in rows:
        if row['archive_sample_ordinal'] in failures:
            matrix, coordinates = row.pop('_RAM_recovery_payloads')
            cache = dict(cell_feature_matrix=matrix, cell_coordinates_and_metrics=coordinates)
            original = module.read_member
            module.read_member = lambda ix, member: cache[member['role']]
            try:
                replay = module.inspect_pair(index, pairs[row['archive_sample_ordinal']], row['archive_sample_ordinal'])
            finally:
                module.read_member = original
            assert replay['status'] == 'ALIGNED_PROCESSED_MATRIX_AND_COORDINATE_KEYS'
            assert all(replay[k] == row[k] for k in ['matrix_SHA256','compressed_coordinate_SHA256','matrix_cells','coordinate_cells','whitelisted_gene_copies'])
            row['original_parser_replayed_on_complete_same_RAM_payloads'] = True
    rows.sort(key=lambda r: r['archive_sample_ordinal'])
    assert len(rows) == 18
    save('private_all18_independent_joint_measurement_acceptance.json', rows)
    save('private_whitelist_gene_panel_coverage.json', [dict(gene=g, measured_archive_files=sum(r['whitelisted_gene_copies'][g] == 1 for r in rows), all_archive_files=18) for g in module.GENES])
    save('private_terminal.json', dict(status='ACCEPTED_ALL18_JOINT_KEYS_WITH_TRANSPORT_RECOVERY_AND_GENE_IDENTITY_AUDIT', completed_UTC=datetime.now(timezone.utc).isoformat(), source_SHA256=sha(__file__), parent_terminal_SHA256=sha(parent / 'private_terminal.json'), all18_Counter_keys_and_finite_coordinates_replayed=True, all17_original_success_payload_SHAs_and_panel_schemas_match=True, original_failed_slots_replayed_with_complete_cached_RAM_inputs=True, CXCL13_Ensembl_alias_and_prefix_checks_complete=True, network=network, no_individual_arrays_unknown_probe_names_or_payloads_saved=True, no_new_spatial_clinical_or_cell_identity_claim=True, files={p.name: dict(SHA256=sha(p), bytes=p.stat().st_size) for p in out.iterdir()}))
    print(json.dumps(dict(status='ACCEPTED_ALL18_JOINT_KEYS_WITH_TRANSPORT_RECOVERY_AND_GENE_IDENTITY_AUDIT')))


if __name__ == '__main__':
    main()
