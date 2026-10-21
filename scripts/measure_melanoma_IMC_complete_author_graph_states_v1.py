"""Complete author graph in RAM, then fixed outcome-blind state organization.

No inferred physical Position units, author code execution, stored edge arrays,
individual exports, threshold tuning, clinical models or independent-cohort claim.
"""
from pathlib import Path
import csv
import hashlib
import importlib.util
import io
import json
import os
import re
import zipfile
import numpy as np
import pandas as pd


class DigestReader(io.RawIOBase):
    def __init__(self, stream):
        self.stream = stream
        self.digest = hashlib.sha256()
        self.count = 0

    def readable(self):return True

    def readinto(self, buffer):
        raw = self.stream.read(len(buffer))
        buffer[:len(raw)] = raw
        self.digest.update(raw)
        self.count += len(raw)
        return len(raw)


def load(path, expected, name):
    assert hashlib.sha256(path.read_bytes()).hexdigest() == expected
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def resolve_rows(rows, target, locations, lookup, aliases, registry):
    assert rows and all(len(row) == 4 for row in rows)
    images = []
    for row in rows:
        assert row[1] in registry or row[1] in aliases
        images.append(row[1] if row[1] in registry else aliases[row[1]])
    first = list(zip(images, [row[2] for row in rows]))
    second = list(zip(images, [row[3] for row in rows]))
    a = target.get_indexer(pd.MultiIndex.from_tuples(first))
    b = target.get_indexer(pd.MultiIndex.from_tuples(second))
    assert (a >= 0).all() and (b >= 0).all()
    # Entire edge stream independently resolved with ordinary Python keys.
    assert np.array_equal(a, [lookup[key] for key in first])
    assert np.array_equal(b, [lookup[key] for key in second])
    return np.asarray(images, dtype=object), locations[a], locations[b]


def selected_edges(reader, frame, target, locations, aliases, margin, role_labels):
    lookup = {key: index for index, key in enumerate(target)}
    registry = set(frame['sample_id'])
    types = frame['celltype20'].to_numpy()
    collected = {}
    total = selected = self_edges = 0
    batch = []
    def accept(rows):
        nonlocal total, selected, self_edges
        images, a, b = resolve_rows(rows, target, locations, lookup, aliases, registry)
        total += len(rows)
        assert total <= 20000000
        self_edges += int((a == b).sum())
        left_a = np.isin(types[a], role_labels[:2])
        left_b = np.isin(types[b], role_labels[:2])
        take = ((left_a & (types[b] == role_labels[2])) |
                (left_b & (types[a] == role_labels[2]))) & (a != b)
        take &= np.isin(images, list(margin))
        selected += int(take.sum())
        for roi in np.unique(images[take]):
            positions = take & (images == roi)
            collected.setdefault(roi, []).append(np.column_stack([a[positions], b[positions]]))
    for row in reader:
        batch.append(row)
        if len(batch) == 50000:
            accept(batch);batch = []
    if batch:accept(batch)
    inverse = np.empty(len(frame), dtype=np.int64)
    groups = frame.groupby('sample_id', sort=True).indices
    for indices in groups.values():inverse[indices] = np.arange(len(indices))
    result = {}
    for roi in margin:
        global_edges = np.concatenate(collected[roi]) if roi in collected else np.empty((0, 2), dtype=np.int64)
        global_edges = np.unique(np.sort(global_edges, axis=1), axis=0)
        if len(global_edges):
            assert (frame['sample_id'].to_numpy()[global_edges] == roi).all()
        result[roi] = inverse[global_edges]
    return result, dict(all_input_edges_resolved_twice=total,
        source_self_edges_excluded_from_state_pairs=self_edges,
        selected_raw_cross_role_edges=selected,
        selected_unique_undirected_cross_role_edges=sum(len(x) for x in result.values()),
        all_author_endpoints_and_ROIs_checked=True,
        duplicate_directions_counted_once=True)


def analyze(bundle, parameters):
    assert set(parameters) == {'regulatory_author_label', 'conventional_CD4_author_label', 'effector_author_label'}
    roles = [parameters[key] for key in ['regulatory_author_label', 'conventional_CD4_author_label', 'effector_author_label']]
    assert len(set(roles)) == 3
    folder = Path(__file__).parent
    root = folder.parent/'private_output'
    accepted = json.loads((root/'private_frozen_step_author_node_native_geometry_v1_source_summary.json').read_text())
    assert accepted['status'] == 'COMPLETE_AUTHOR_NODE_BIJECTION_AND_EDGE_PREFIX_NATIVE_COORDINATE_REVIEW'
    assert accepted['complete_node_member_CRC_checked'] and accepted['independent_CSV_Pandas_node_replay']
    node = load(folder/'qualify_melanoma_IMC_native_graph_key_bridge_local_v4.py',
        '427434209f230231011232b718e286eca77d6b779f72cdbd8d1ac6b182ff6eef', 'author_node_bind')
    transport = load(Path(os.environ['PANCANCER_IMC_BOUNDED_TRANSPORT']),
        '1536f71ebb4399ec3d1b3fe8239b02a90b892ca3a420c8f934fd183226bd8329', 'author_graph_transport')
    rf = transport.RangeFile(None, 'https://zenodo.org/api/records/6838169/files/Melanoma_code_data.zip/content',668227028,
        checkpoint=lambda value:(root/'private_author_graph_transport_progress.json').write_text(json.dumps(value,indent=2)+'\n'))
    def directory(archive):
        return hashlib.sha256(json.dumps([(e.filename,e.file_size,e.compress_size,e.CRC) for e in archive.infolist()],separators=(',',':')).encode()).hexdigest()
    with zipfile.ZipFile(rf) as archive:
        assert directory(archive) == '26339da50b39ffdf7beb641a67f50523c867f2d3a86b24f42175f3a16dfe8992'
        code = archive.read('Rscript/Spatial_analysis.Rmd')
        assert hashlib.sha256(code).hexdigest() == 'd8faddad268f2900a83f8cfa3f8352dd697128d408a528181fb997c5773dba56'
        members = {}
        for role in ['sc_data','n_dat']:
            line = next(line for line in code.decode('utf-8-sig').splitlines() if re.match(r'\s*'+role+r'\s*<-\s*read.csv',line))
            literal = re.search(r'''read.csv\(\s*(["'])([^"'\r\n]+)\1''',line).group(2)
            while literal.startswith('../') or literal.startswith('./'):literal = literal.split('/',1)[1]
            assert literal in archive.namelist()
            members[role] = literal
        assert archive.getinfo(members['sc_data']).file_size <= 40*1024**2
        blob = archive.read(members['sc_data'])
        assert hashlib.sha256(blob).hexdigest() == accepted['complete_node_member_SHA256']
        nodes = pd.read_csv(io.BytesIO(blob), dtype=str,keep_default_na=False).rename(columns={'Unnamed: 0':''})
        target, locations, _, _ = node.bind_nodes(bundle['frame'],nodes,bundle['raw_ROI_aliases'])
        del nodes,blob
        info = archive.getinfo(members['n_dat'])
        assert info.file_size <= 300*1024**2 and info.compress_size <= 48*1024**2
        # TextIOWrapper reads to EOF, enforcing complete ZIP CRC. No CSV
        # member, row key or graph array is written to a filesystem.
        with archive.open(info) as stream:
            digest_reader = DigestReader(stream)
            reader = csv.reader(io.TextIOWrapper(io.BufferedReader(digest_reader),encoding='utf-8-sig',newline=''))
            header = next(reader);assert len(header) == 4
            edges, qualification = selected_edges(reader,bundle['frame'],target,locations,
                bundle['raw_ROI_aliases'],bundle['ROI_to_source_patient'],roles)
            assert digest_reader.count == info.file_size
        qualification.update(complete_edge_member_CRC_checked=True,
            complete_edge_member_SHA256=digest_reader.digest.hexdigest(),
            complete_node_member_SHA256=accepted['complete_node_member_SHA256'],
            source_edge_member_uncompressed_bytes=info.file_size,
            actual_physical_distance_threshold_not_reconstructed=True,
            author_paper_describes_20um_centroid_graph=True)
    rf.clear_cache()
    with zipfile.ZipFile(rf) as replay:
        assert directory(replay) == '26339da50b39ffdf7beb641a67f50523c867f2d3a86b24f42175f3a16dfe8992'
    rf.close()
    qualification['fresh_wire_directory_replayed'] = True
    (root/'private_author_graph_input_qualification.json').write_text(json.dumps(qualification,indent=2)+'\n')
    core = load(folder/'analyze_melanoma_IMC_author_graph_states_v1.py', '57d4af556e450f5a267b7a2e10caf3a097250090e35a97669a62d2c104ee87be','author_graph_state_core')
    result = core.analyze(dict(bundle,author_neighbor_edges=edges),parameters)
    result['author_graph_input_qualification'] = qualification
    result['wrapper_SHA256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return result
