"""RAM-session module: full author node keys plus bounded graph-prefix audit.

Read only the explicitly referenced small processed node table and a bounded
edge prefix. A complete key AND label bijection must identify the main rows.
This is measurement provenance, not a spatial or clinical association test.
"""
from datetime import datetime, timezone
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


def sha(value):
    return hashlib.sha256(value).hexdigest()


def load(path, expected, name):
    assert sha(path.read_bytes()) == expected
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def canonical_images(values, aliases, registry):
    result = []
    for value in values:
        assert value in registry or value in aliases, 'unmapped_author_image_key'
        result.append(value if value in registry else aliases[value])
    return result


def bind_nodes(frame, nodes, aliases):
    """Accept only equivalent complete literal key/label bijections."""
    assert set(nodes) == {'', 'Imageid', 'type', 'cell_id'}
    assert len(nodes) == len(frame), 'different_node_universes'
    canonical = canonical_images(nodes['Imageid'], aliases, set(frame['sample_id']))
    target = pd.MultiIndex.from_arrays([canonical, nodes['cell_id']])
    assert target.is_unique, 'duplicate_author_node_keys'
    accepted = []
    mappings = []
    diagnostics = []
    for field in ['X', 'Unnamed: 0']:
        keys = pd.MultiIndex.from_arrays([frame['sample_id'], frame[field]])
        unique = keys.is_unique
        locations = keys.get_indexer(target) if unique else np.full(len(nodes), -1)
        complete = bool(unique and (locations >= 0).all())
        labels = bool(complete and np.array_equal(
            frame['celltype20'].to_numpy()[locations], nodes['type'].to_numpy()))
        diagnostics.append(dict(identity_field=field, keys_unique=bool(unique),
            complete_literal_key_bijection=complete, all_author_labels_equal=labels))
        if labels:
            accepted.append(field)
            mappings.append(locations)
    assert mappings, 'no_complete_key_and_label_bijection'
    assert all(np.array_equal(mappings[0], x) for x in mappings[1:]), 'ambiguous_row_identity'
    return target, mappings[0], accepted, diagnostics


def edge_prefix_summary(rows, target, locations, coords, aliases, registry):
    assert len(rows) > 0
    assert all(len(row) == 4 for row in rows)
    images = canonical_images([row[1] for row in rows], aliases, registry)
    first = pd.MultiIndex.from_arrays([images, [row[2] for row in rows]])
    second = pd.MultiIndex.from_arrays([images, [row[3] for row in rows]])
    a = target.get_indexer(first)
    b = target.get_indexer(second)
    assert (a >= 0).all() and (b >= 0).all(), 'edge_endpoint_not_in_complete_nodes'
    main_a, main_b = locations[a], locations[b]
    distances = np.sqrt(np.sum((coords[main_a]-coords[main_b])**2, axis=1))
    assert np.isfinite(distances).all()
    return dict(prefix_edges=len(rows), all_prefix_endpoints_have_unique_main_rows=True,
        self_edges=int((main_a == main_b).sum()),
        zero_native_distance_edges=int((distances == 0).sum()),
        native_distance_quantiles=np.quantile(distances, [.0, .25, .5, .75, 1.0]).tolist(),
        prefix_not_random_or_full_graph=True, physical_micron_scale_not_inferred=True)


def analyze(bundle, parameters):
    assert parameters == dict(maximum_edge_prefix_rows=1024,
        maximum_edge_prefix_bytes=1048576, maximum_node_member_bytes=40*1024**2)
    source = load(Path(os.environ['PANCANCER_IMC_TRANSPORT_HELPER']),
        '9fb200d98e321584fde8bca6cb92f72151ac33ec00d64e04d93e7e2b33e870e2', 'graph_source')
    transport = load(Path(os.environ['PANCANCER_IMC_BOUNDED_TRANSPORT']),
        '709f9b5ced875fea663557bad936f32d004e12d5daf291048b7de476064b0d68', 'graph_transport')
    record = json.loads(source.fetch(source.BASE+'api/records/6838169', 2*1024**2))
    item = next(x for x in record['files'] if x['key'] == 'Melanoma_code_data.zip')
    assert item['size'] == 668227028 and item['checksum'] == 'md5:da74b7228763b76dcd2d9830c1a10cd9'
    rf = transport.RangeFile(source, item['links']['self'], item['size'])
    frame = bundle['frame']
    with zipfile.ZipFile(rf) as archive:
        functions = archive.read('Rscript/functions.R')
        assert sha(functions) == 'bfb6f99c80947643655aa5d4d8ef310e5f0a1c24aac86326482cb32c123565ee'
        compact = re.sub(r'\s+', '', functions.decode('utf-8-sig')).replace('"', "'")
        assert "colnames(n_dat)[2:4]=c('ImageNumber','firstobjectnumber','secondobjectnumber')" in compact
        assert 'col4<-n_dat$firstobjectnumber' in compact and 'col7<-n_dat$secondobjectnumber' in compact
        code = archive.read('Rscript/Spatial_analysis.Rmd')
        assert sha(code) == 'd8faddad268f2900a83f8cfa3f8352dd697128d408a528181fb997c5773dba56'
        members = {}
        for role in ['sc_data', 'n_dat']:
            line = next(s for s in code.decode('utf-8-sig').splitlines()
                if re.match(r'\s*'+role+r'\s*<-\s*read.csv', s))
            match = re.search(r'''read.csv\(\s*(["'])([^"'\r\n]+)\1''', line)
            assert match is not None
            literal = match.group(2)
            while literal.startswith('../') or literal.startswith('./'):
                literal = literal.split('/', 1)[1]
            assert literal in archive.namelist()
            members[role] = literal
        info = archive.getinfo(members['sc_data'])
        assert info.file_size <= parameters['maximum_node_member_bytes']
        blob = archive.read(info)
        nodes = pd.read_csv(io.BytesIO(blob), dtype=str, keep_default_na=False)
        nodes = nodes.rename(columns={'Unnamed: 0': ''})
        scalar = list(csv.DictReader(io.StringIO(blob.decode('utf-8-sig'))))
        assert len(scalar) == len(nodes)
        assert list(scalar[0]) == list(nodes)
        assert all([r[c] for r in scalar] == nodes[c].tolist() for c in nodes)
        del scalar
        try:
            target, locations, fields, diagnostics = bind_nodes(frame, nodes, bundle['raw_ROI_aliases'])
        except AssertionError as error:
            allowed = {'different_node_universes', 'unmapped_author_image_key', 'duplicate_author_node_keys',
                'no_complete_key_and_label_bijection', 'ambiguous_row_identity', ''}
            assert str(error) in allowed
            rf.close()
            return dict(status='AUTHOR_NODE_KEY_OR_LABEL_ALIGNMENT_NOT_ACCEPTED',
                completed_UTC=datetime.now(timezone.utc).isoformat(),
                producer_SHA256=sha(Path(__file__).read_bytes()), source_nodes=len(nodes),
                main_cells=len(frame), complete_node_member_CRC_checked=True,
                complete_node_member_SHA256=sha(blob), independent_CSV_Pandas_node_replay=True,
                reason=str(error) or 'unexpected_node_schema', no_missing_records_filled_or_forced=True,
                no_individual_arrays_or_rows_written=True, no_spatial_or_clinical_effect=True)
        with archive.open(members['n_dat']) as stream:
            header = stream.readline(16385)
            assert len(header) <= 16384 and header.endswith(b'\n')
            assert len(next(csv.reader([header.decode('utf-8-sig')]))) == 4
            prefix = bytearray()
            for _ in range(parameters['maximum_edge_prefix_rows']):
                remaining = parameters['maximum_edge_prefix_bytes']-len(prefix)
                assert remaining > 0
                line = stream.readline(min(remaining+1, 16385))
                if not line:
                    break
                assert len(line) <= remaining and line.endswith(b'\n')
                prefix.extend(line)
        rows = list(csv.reader(io.StringIO(prefix.decode('utf-8-sig'))))
        geometry = edge_prefix_summary(rows, target, locations, bundle['native_coordinates'],
            bundle['raw_ROI_aliases'], set(frame['sample_id']))
    again = next(x for x in json.loads(source.fetch(source.BASE+'api/records/6838169', 2*1024**2))['files']
        if x['key'] == item['key'])
    assert all(item.get(k) == again.get(k) for k in ['id', 'size', 'checksum', 'links'])
    rf.close()
    return dict(status='COMPLETE_AUTHOR_NODE_BIJECTION_AND_EDGE_PREFIX_NATIVE_COORDINATE_REVIEW',
        completed_UTC=datetime.now(timezone.utc).isoformat(), producer_SHA256=sha(Path(__file__).read_bytes()),
        source_nodes=len(nodes), complete_node_member_CRC_checked=True,
        complete_node_member_SHA256=sha(blob), independent_CSV_Pandas_node_replay=True,
        main_identity_fields_with_equivalent_complete_maps=fields, identity_diagnostics=diagnostics,
        prefix_geometry=geometry, complete_edge_member_CRC_not_checked=True,
        no_full_author_graph_or_physical_unit_claim=True,
        no_new_state_spatial_or_clinical_effect=True, no_individual_arrays_or_rows_written=True)
