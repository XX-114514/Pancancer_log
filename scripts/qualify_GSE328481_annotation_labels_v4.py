"""Biological-label-qualified annotation vocabulary; labels are not identity truth.

Never read matrix, coordinates, observation identifiers or receptor fields.
Only labels passing the fixed biological-word / native-gene vocabulary gate may be emitted. Other strings remain RAM only.
"""
from collections import Counter
from pathlib import Path
from datetime import datetime,timezone
import hashlib
import json
import io
import numpy as np

FIELDS=['anno','anno2','classes','bin50 structural annotation','TLScluster']


def normalized(value):
    value=value.decode() if isinstance(value,bytes) else str(value)
    return ''.join(c for c in value.lower() if c.isalnum())


def vocabulary():
    groups={
        'T':['T cell','T cells','T_cell','T lymphocytes'],
        'T_NK':['T/NK cell','T/NK cells','T_NK'],
        'CD4_T':['CD4 T','CD4 T cell','CD4 T cells','CD4T'],
        'CD8_T':['CD8 T','CD8 T cell','CD8 T cells','CD8T'],
        'Treg':['Treg','Tregs','Regulatory T','Regulatory T cell','Regulatory T cells','CD4T_Treg'],
        'B':['B','B cell','B cells','B_cell'],
        'plasma':['plasma','plasma cell','plasma cells','plasmablast'],
        'NK':['NK','NK cell','NK cells','Natural killer cells'],
        'myeloid':['myeloid','myeloid cells','macrophage','macrophages','monocyte','monocytes'],
        'DC':['DC','DCs','Dendritic cell','Dendritic cells'],
        'endothelial':['Endothelial','Endothelial cell','Endothelial cells'],
        'fibroblast':['Fibroblast','Fibroblasts','Fibroblast cells'],
        'epithelial':['Epithelial','Epithelial cell','Epithelial cells'],
        'tumor_region':['Tumor','Tumour','Tumor region','Tumor cells','Cancer','Cancer cells'],
        'stroma_region':['stroma','stromal','immune-excluded stroma'],
        'iTLS_region':['iTLS','iTLSs','immature TLS','immature TLSs'],
        'TLS_region':['TLS','TLSs','TLS region','tertiary lymphoid structure'],
        'alveolar_region':['alveolar','alveolar tissue'],
        'cilia_region':['lung cilia','cilia'],
        'plasma_aggregate_region':['plasma cell aggregates','plasma cell aggregate'],
        'myeloid_enriched_region':['myeloid cell-enriched region','myeloid enriched region'],
        'normal_region':['normal','normal tissue'],
    }
    for gene in ['CCR8','FOXP3','IL2RA','MKI67','LAG3','TNFRSF9']:
        groups['Treg_'+gene]=[p+gene for p in ['Treg_','CD4T_Treg_','Regulatory T_']]
    for gene in ['MKI67','GZMK','IL7R','ZNF683','CXCL13','LAYN','ISG15']:
        groups['CD8_T_'+gene]=[p+gene for p in ['CD8T_','CD8_T_','CD8 T_']]
    out={}
    for key,values in groups.items():
        for value in values:
            n=normalized(value);assert n not in out or out[n]==key;out[n]=key
    return out


def summarize(h):
    import h5py
    obs=h['obs'];idx=obs.attrs.get('_index','_index');idx=idx.decode() if isinstance(idx,bytes) else str(idx)
    n=len(obs[idx]);aliases=vocabulary();fields={};budget=0
    var=h['var'];vi=var.attrs.get('_index','_index');vi=vi.decode() if isinstance(vi,bytes) else str(vi)
    assert isinstance(var[vi],h5py.Dataset) and len(var[vi])<=100000
    genes={str(v.decode() if isinstance(v,bytes) else v).strip().upper() for v in var[vi][:]}
    assert len(genes)==len(var[vi])
    gene_vocabulary_SHA256=hashlib.sha256(json.dumps(sorted(genes),separators=(',',':')).encode()).hexdigest()
    for key in FIELDS:
        if key not in obs:fields[key]=dict(status='FIELD_ABSENT');continue
        obj=obs[key];counts=Counter();missing=0
        if isinstance(obj,h5py.Group) and set(['categories','codes'])<=set(obj.keys()):
            categories=obj['categories']
            if not isinstance(categories,h5py.Dataset) or categories.dtype.kind not in ['O','S','U']:
                fields[key]=dict(status='NUMERIC_OR_NESTED_CATEGORIES_NOT_INTERPRETED');continue
            if len(obj['categories'])>256:fields[key]=dict(status='CATEGORY_COUNT_OVER_CAP');continue
            raw=list(obj['categories'][:]);assert all(len(str(v))<=512 for v in raw)
            assert len(obj['codes'])==n;rawcounts=np.zeros(len(raw),np.int64)
            for lo in range(0,n,4096):
                codes=np.asarray(obj['codes'][lo:lo+4096],int);assert ((codes>=-1)&(codes<len(raw))).all()
                missing+=int((codes<0).sum());rawcounts+=np.bincount(codes[codes>=0],minlength=len(raw))
            for label,count in zip(raw,rawcounts):
                label=label.decode() if isinstance(label,bytes) else str(label);counts[label]+=int(count)
        elif isinstance(obj,h5py.Dataset) and obj.ndim==1 and obj.dtype.kind in ['O','S','U']:
            assert len(obj)==n
            for lo in range(0,n,4096):
                values=obj[lo:lo+4096]
                budget+=sum(len(str(v).encode()) for v in values);assert budget<=8*1024**2
                for value in values:
                    value=value.decode() if isinstance(value,bytes) else str(value);counts[value]+=1
        else:
            fields[key]=dict(status='NUMERIC_OR_UNSUPPORTED_FIELD_NOT_INTERPRETED',shape=list(obj.shape) if isinstance(obj,h5py.Dataset) else None);continue
        recognized=Counter();qualified=Counter();unknown_categories=unknown_rows=0
        for label,count in counts.items():
            if normalized(label) in aliases:recognized[aliases[normalized(label)]]+=count
            else:unknown_categories+=1;unknown_rows+=count
            if safe_label(label,genes):qualified[label]+=count
        assert sum(recognized.values())+unknown_rows+missing==n
        fields[key]=dict(status='ALLOWLISTED_CATEGORY_COUNTS_ONLY',recognized_counts=dict(recognized),biological_vocabulary_qualified_label_counts=dict(qualified),
            unknown_category_count=unknown_categories,unknown_observation_count=unknown_rows,missing_observation_count=missing,
            unqualified_raw_values_never_emitted=True)
    return dict(source='GSE328481',scope='Single previously fixed smallest author cell-bin pilot, not full11 or patient validation',
        author_observations=n,fields=fields,only_annotation_aggregate_fields_read=True,
        no_expression_coordinates_identifiers_clonotypes_or_clinical_rows_read=True,
        category_name_match_not_identity_truth=True,gene_symbol_vocabulary_SHA256=gene_vocabulary_SHA256)


def synthetic():
    import h5py
    b=io.BytesIO()
    with h5py.File(b,'w') as h:
        v=h.create_group('var');v.create_dataset('_index',data=np.array([b'FOXP3',b'LAG3']))
        o=h.create_group('obs');o.attrs['_index']='_index';o.create_dataset('_index',data=np.array([b'id1',b'id2',b'id3']))
        g=o.create_group('anno');g.create_dataset('categories',data=np.array([b'Treg',b'unknown_private_string']))
        g.create_dataset('codes',data=np.array([0,1,-1],np.int8));o.create_dataset('TRAcdr3aa',data=np.array([b'never_read']*3))
        r=summarize(h);v=r['fields']['anno'];assert v['recognized_counts']=={'Treg':1} and v['unknown_observation_count']==1 and v['missing_observation_count']==1
        encoded=json.dumps(r);assert all(x not in encoded for x in ['unknown_private_string','never_read','id1'])
    b.close()
    b=io.BytesIO()
    with h5py.File(b,'w') as h:
        v=h.create_group('var');v.create_dataset('_index',data=np.array([b'FOXP3',b'LAG3']))
        o=h.create_group('obs');o.create_dataset('_index',shape=(3,),dtype='i8')
        a=o.create_group('anno');a.create_dataset('categories',data=np.array([b'Treg']))
        a.create_dataset('codes',data=np.array([0,0,-1],np.int8))
        t=o.create_group('TLScluster');c=t.create_group('categories');c.attrs['encoding-type']='nullable-integer'
        c.create_dataset('values',data=np.array([1]));c.create_dataset('mask',data=np.array([False]))
        t.create_dataset('codes',data=np.zeros(3,np.int8))
        a=o.create_group('classes');a.create_dataset('categories',data=np.array([1,2]))
        a.create_dataset('codes',data=np.array([0,0,1],np.int8))
        result=summarize(h)
        assert result['fields']['anno']['recognized_counts']=={'Treg':2}
        for key in ['TLScluster','classes']:
            assert result['fields'][key]['status']=='NUMERIC_OR_NESTED_CATEGORIES_NOT_INTERPRETED'
    b.close()


def read_source(root,plan):
    import h5py
    synthetic()
    folder=Path(root)/'GSE328481_public_pilot_stream_transfer_20261018_v2'
    terminal=json.loads((folder/'private_terminal.json').read_text());name=Path(terminal['chosen']['filename'])
    assert name.name==str(name) and terminal['download_sha256']==plan['expected_file_SHA256']
    path=folder/name;assert path.stat().st_size==plan['expected_file_bytes']
    digest=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(16*1024**2),b''):digest.update(block)
    assert digest.hexdigest()==plan['expected_file_SHA256']
    with h5py.File(path,'r') as h:result=summarize(h)
    result.update(completed_UTC=datetime.now(timezone.utc).isoformat(),existing_file_SHA256=digest.hexdigest(),
        existing_file_bytes=path.stat().st_size,synthetic_unknown_values_sequence_and_identifier_nonemission_checked=True,
        no_new_download_or_server_writes=True)
    return result
