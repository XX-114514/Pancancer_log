"""Qualified public Xenium source context; raw objects remain in process RAM."""
from collections import Counter
from pathlib import Path
import gzip
import hashlib
import importlib.util
import io
import json
import os
import re
import tarfile
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET

import h5py
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.sparse import csc_matrix

REQUIRED = ['CD3D','CD3E','TRAC','CD4','CD8A','FOXP3','IL2RA','CTLA4','TNFRSF9','LAG3','PDCD1','TOX','HAVCR2','MKI67','TOP2A','UBE2C','TCF7','IL7R','MS4A1','CD79A','EPCAM','SOX2']
R_NAMES = ['TNFRSF9_positive','LAG3_positive','cycling_marker_positive']
E_NAMES = ['dysfunction_marker_positive','cycling_marker_positive','memory_marker_positive']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verified(path):
    t = json.loads((path / 'private_terminal.json').read_text())
    for name, v in t.get('files', {}).items():
        assert sha(path / name) == v['SHA256']
    return t


def read_member(index, member):
    start, n = member['payload_offset'], member['bytes']
    assert 0 < n <= 64*1024**2
    chunks = []
    for first in range(start, start+n, 4*1024**2):
        size = min(4*1024**2, start+n-first)
        for attempt in range(3):
            headers = {'Range': f'bytes={first}-{first+size-1}', 'Accept-Encoding':'identity'}
            headers.update({'If-Match':index['source_ETag']} if index.get('source_ETag') else {'If-Unmodified-Since':index['source_Last_Modified']})
            try:
                with urllib.request.urlopen(urllib.request.Request(index['url'], headers=headers), timeout=45) as response:
                    assert response.status == 206 and response.headers.get('Content-Range') == f'bytes {first}-{first+size-1}/{index["archive_bytes"]}'
                    assert response.headers.get('ETag') == index.get('source_ETag') and response.headers.get('Last-Modified') == index.get('source_Last_Modified')
                    b = response.read(size+1)
                if len(b) != size:
                    raise EOFError('Short exact source range')
                chunks.append(b)
                break
            except (urllib.error.URLError, TimeoutError, EOFError):
                if attempt == 2:
                    raise
    raw = b''.join(chunks)
    assert len(raw) == n
    return raw


def site_category(value):
    s = value.casefold()
    if 'lymph' in s or s.strip() in ['ln','lymph node']:
        return 'explicit_lymph_node_site'
    if 'neck' in s:
        return 'neck_site_unspecified'
    if 'oropharyn' in s or 'tonsil' in s or 'base of tongue' in s:
        return 'oropharyngeal_site'
    if 'oral' in s:
        return 'oral_cavity_site'
    if 'hypopharyn' in s:
        return 'hypopharyngeal_site'
    if 'laryn' in s:
        return 'laryngeal_site'
    return 'other_or_unresolved_site'


def load_context():
    parent = Path(os.environ['PANCANCER_JOINT_PREFLIGHT'])
    joint = Path(os.environ['PANCANCER_JOINT_AUDIT'])
    units = Path(os.environ['PANCANCER_SOURCE_UNIT_AUDIT'])
    old, accepted, source_units = verified(parent), verified(joint), verified(units)
    assert accepted['all18_Counter_keys_and_finite_coordinates_replayed']
    src = parent.parent / 'preflight_public_xenium_joint_measurement_v1.py'
    assert sha(src) == old['source_SHA256']
    spec = importlib.util.spec_from_file_location('pinned_panel_reader', src)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    index, pairs = module.load_index()
    proofs = json.loads((joint / 'private_all18_independent_joint_measurement_acceptance.json').read_text())
    assert len(proofs) == len(pairs) == 18
    url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE300nnn/GSE300147/miniml/GSE300147_family.xml.tgz'
    with urllib.request.urlopen(url, timeout=25) as r:
        b = r.read(2*1024**2+1)
    assert len(b) <= 2*1024**2 and hashlib.sha256(b).hexdigest() == source_units['MINiML_archive_SHA256']
    with tarfile.open(fileobj=io.BytesIO(b), mode='r:gz') as archive:
        document = None
        for entry in archive:
            if entry.isfile() and entry.name.endswith('.xml'):
                assert entry.size <= 16*1024**2
                document = archive.extractfile(entry).read(16*1024**2+1)
                break
    assert document is not None and len(document) <= 16*1024**2
    records = {}
    for element in ET.fromstring(document).findall('.//{*}Sample'):
        title = element.findtext('{*}Title') or ''
        m = re.fullmatch(r'Patient ([0-9]+), run ([0-9]+)', title); assert m
        source = (element.findtext('{*}Channel/{*}Source') or '').casefold()
        fields = {v.attrib.get('tag','').strip().casefold().replace(' ','_'): (v.text or '').strip() for v in element.findall('{*}Channel/{*}Characteristics') if v.attrib.get('tag','').strip().casefold().replace(' ','_') in ['recurrence','p16_status','site_of_tumor_resection']}
        recurrence = fields.get('recurrence','').casefold()
        recurrence = 'primary' if recurrence == 'primary' else ('recurrent' if recurrence in ['recurrent','recurrence','locoregional recurrence','local recurrence'] else 'other_or_unknown')
        diagnosis = 'ameloblastoma' if 'ameloblastoma' in source else ('HNSCC' if 'hnscc' in source or 'squamous' in source else 'other_or_unresolved')
        records[element.attrib['iid']] = dict(case=m[1], run=m[2], diagnosis=diagnosis, recurrence=recurrence, anatomical_site=site_category(fields.get('site_of_tumor_resection','')))
    for i, pair in enumerate(pairs):
        key = re.match(r'GSM[0-9]+',pair['cell_feature_matrix']['name'])[0]
        assert proofs[i]['archive_sample_ordinal'] == i
        pair['source_metadata'] = records[key]
    assert len(records) == len(pairs) == 18
    assert Counter(p['source_metadata']['diagnosis'] for p in pairs) == Counter(HNSCC=17, ameloblastoma=1)
    return dict(index=index, pairs=pairs, proofs=proofs, gene_whitelist=module.GENES, joint_audit_SHA256=sha(joint/'private_terminal.json'), unit_audit_SHA256=sha(units/'private_terminal.json'))


def load_pair(context, ordinal):
    pair, expected = context['pairs'][ordinal], context['proofs'][ordinal]
    mb = read_member(context['index'],pair['cell_feature_matrix'])
    cb = read_member(context['index'],pair['cell_coordinates_and_metrics'])
    assert hashlib.sha256(mb).hexdigest() == expected['matrix_SHA256'] and hashlib.sha256(cb).hexdigest() == expected['compressed_coordinate_SHA256']
    decode = lambda xs: [x.decode() if isinstance(x,bytes) else str(x) for x in xs]
    with h5py.File(io.BytesIO(mb),'r') as h:
        m=h['matrix']; names=decode(m['features/name'][:]); ids=decode(m['features/id'][:]); types=decode(m['features/feature_type'][:]); cells=decode(m['barcodes'][:]); shape=tuple(int(x) for x in m['shape'][:])
        ptr=np.asarray(m['indptr'][:],np.int64); indices=np.asarray(m['indices'][:],np.int64); raw_data=np.asarray(m['data'][:]); assert raw_data.dtype.kind in 'iuf' and np.isfinite(raw_data).all() and (raw_data>=0).all() and (raw_data==np.floor(raw_data)).all() and raw_data.max(initial=0)<2**63; data=raw_data.astype(np.int64)
    assert shape == (len(names),len(cells)) and len(ptr)==len(cells)+1 and ptr[0]==0 and ptr[-1]==len(data)==len(indices) and (np.diff(ptr)>=0).all()
    assert (data>=0).all() and (indices>=0).all() and (indices<shape[0]).all()
    canonical=np.array([t=='Gene Expression' and bool(re.fullmatch(r'ENSG[0-9]+(?:\.[0-9]+)?',i)) for i,t in zip(ids,types)])
    assert canonical.any()
    known=set(context['gene_whitelist']); canonical|=np.array([t=='Gene Expression' and n in known for n,t in zip(names,types)])
    selected=[]
    for g in REQUIRED:
        found=[i for i,(n,t) in enumerate(zip(names,types)) if n==g and t=='Gene Expression']; assert len(found)==1; selected.append(found[0])
    matrix=csc_matrix((data,indices,ptr),shape=shape,copy=True); matrix.sum_duplicates(); matrix.sort_indices()
    human=np.asarray(matrix[canonical].sum(axis=0)).ravel().astype(np.int64)
    values=matrix[selected].toarray().T.astype(np.int64)
    assert data.sum(dtype=np.float64)<2**53
    columns=np.repeat(np.arange(shape[1],dtype=np.int64),np.diff(ptr))
    human_entries=canonical[indices]
    direct_human=np.bincount(columns[human_entries],weights=data[human_entries],minlength=shape[1]).astype(np.int64)
    slots=np.full(shape[0],-1,np.int64);slots[selected]=np.arange(len(selected));wanted=slots[indices]>=0
    direct=np.bincount(columns[wanted]*len(selected)+slots[indices[wanted]],weights=data[wanted],minlength=shape[1]*len(selected)).reshape(shape[1],len(selected)).astype(np.int64)
    np.testing.assert_array_equal(human,direct_human);np.testing.assert_array_equal(values,direct)
    raw=cb
    if raw[:2]==b'\x1f\x8b':
        with gzip.GzipFile(fileobj=io.BytesIO(raw)) as f:raw=f.read(512*1024**2+1)
        assert len(raw)<=512*1024**2
    coordinates=pq.read_table(io.BytesIO(raw),columns=['cell_id','x_centroid','y_centroid','cell_area','nucleus_area']).to_pandas()
    assert coordinates.cell_id.astype(str).is_unique
    order=pd.Index(coordinates.cell_id.astype(str)).get_indexer(cells); assert (order>=0).all() and len(order)==len(coordinates)
    coordinates=coordinates.iloc[order].reset_index(drop=True).drop(columns='cell_id')
    assert np.isfinite(coordinates[['x_centroid','y_centroid']].to_numpy(float)).all()
    return dict(metadata=pair['source_metadata'], counts=values, human_gene_molecules=human, coordinates=coordinates, canonical_human_features=int(canonical.sum()), total_features=shape[0], matrix_cells=shape[1],raw_CSC_marker_columns_and_human_totals_independent_bincount_checked=True)


def observed_masks(data, minimum, parent_definition):
    x={g:data['counts'][:,i] for i,g in enumerate(REQUIRED)}
    area=data['coordinates'].cell_area.to_numpy(float)
    qc=(data['human_gene_molecules']>=minimum)&np.isfinite(area)&(area>0)
    t=(x['CD3D']>=1)|(x['CD3E']>=1)|(x['TRAC']>=1)
    regulatory=t&(x['FOXP3']>=1)&((x['IL2RA']>=1) if parent_definition=='dual_regulatory_evidence' else ((x['IL2RA']>=1)|(x['CTLA4']>=1)))&qc
    cd8=t&(x['CD8A']>=1)&qc; ambiguous=regulatory&cd8; regulatory&=~ambiguous; cd8&=~ambiguous
    cycling=(x['MKI67']>=1)&((x['TOP2A']>=1)|(x['UBE2C']>=1))
    rs=np.column_stack([x['TNFRSF9']>=1,x['LAG3']>=1,cycling])&regulatory[:,None]
    es=np.column_stack([(x['PDCD1']>=1)&((x['TOX']>=1)|(x['HAVCR2']>=1)),cycling,(x['TCF7']>=1)|(x['IL7R']>=1)])&cd8[:,None]
    b=(x['MS4A1']>=1)&(x['CD79A']>=1)&qc&~t; epithelial=(x['EPCAM']>=1)&qc&~t
    return dict(qc=qc,T_support=t&qc,regulatory=regulatory,CD8=cd8,ambiguous=ambiguous,R_states=rs,E_states=es,B_cycling=b&cycling,epithelial_cycling=epithelial&cycling)
