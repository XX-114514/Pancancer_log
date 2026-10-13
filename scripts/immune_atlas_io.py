"""Native sparse source adapters shared by the immune atlas stages."""
import json
from pathlib import Path
import re
import h5py
import anndata as ad
from anndata._io.specs import read_elem
import numpy as np
import pandas as pd
from scipy import sparse
import zarr

from public_gene_identity_contract import build_contract, ID_PATTERN

_OPEN_CACHE={}


class Source:
    def __init__(self, record):
        self.record=record
        path=Path(record['path'])
        if str(path) not in _OPEN_CACHE:
            handle=h5py.File(path,'r') if path.suffix=='.h5ad' else zarr.open_group(str(path),mode='r')
            _OPEN_CACHE[str(path)]=(handle,read_elem(handle['obs']))
        self.handle,self.obs=_OPEN_CACHE[str(path)]
        self.matrix=ad.io.sparse_dataset(self.handle[record['counts_key']])
        self.rows=np.load(record['row_positions_ref'])
        self.var=read_elem(self.handle[record['var_key']])
        self.columns=np.arange(len(self.var))
        if record.get('new_parent_dataset'):
            order=list(read_elem(self.handle['uns'])['source_order'])
            measured=np.asarray(self.handle['varm/measured_in_source'])[:,order.index(record['new_parent_dataset'])]
            self.columns=np.flatnonzero(measured)
            self.var=self.var.iloc[self.columns].copy()
        self.symbols=None
        self.canonical_genes=None

    def close(self):
        if isinstance(self.handle,h5py.File):self.handle.close()

    def identity_view(self):
        var=self.var.copy()
        original=var.index.astype(str)
        symbols=pd.Series(original,index=np.arange(len(var)))
        for name in ['gene_symbol','feature_name','gene_name_original']:
            if name in var:
                names=var[name].astype(str).to_numpy()
                valid=~pd.Series(names).isin(['','nan','None','__SOURCE_VALUE_MISSING__']).to_numpy()
                symbols.loc[valid]=names[valid]
                break
        self.symbols=symbols.astype(str).to_numpy()
        if 'gene_id_original' not in var:
            var['gene_id_original']=[x if ID_PATTERN.fullmatch(x) else '' for x in original]
        # Preserve genuinely ambiguous parent keys instead of repairing them by name.
        for i,name in enumerate(original):
            if name.startswith('AMBIGUOUS::'):
                self.symbols[i]=name
        var.index=pd.Index(self.symbols)
        return var

    def raw(self, indices):
        native_rows=self.rows[np.asarray(indices)]
        raw=sparse.csr_matrix(self.matrix[native_rows,:])[:,self.columns].tocsr()
        if not np.isfinite(raw.data).all() or (raw.data<0).any() or (raw.data!=np.rint(raw.data)).any():
            raise ValueError('Native count contract failed: '+self.record['dataset'])
        return raw

    def analysis(self, indices):
        raw=self.raw(indices)
        total=np.asarray(raw.sum(axis=1)).ravel().astype(np.float64)
        corrected=raw.astype(np.float64)
        source_rows=self.rows[np.asarray(indices)]
        if self.record.get('ambient_removed_key'):
            removal=ad.io.sparse_dataset(self.handle[self.record['ambient_removed_key']])[source_rows,:]
            corrected=corrected-sparse.csr_matrix(removal)[:,self.columns]
            corrected.eliminate_zeros()
            if (corrected.data < -1e-7).any():raise ValueError('Correction exceeds source count')
            corrected.data=np.maximum(corrected.data,0)
            corrected.eliminate_zeros()
            total=np.asarray(corrected.sum(axis=1)).ravel()
        elif self.record.get('count_denominator_column'):
            native=pd.to_numeric(self.obs.iloc[source_rows][self.record['count_denominator_column']],errors='raise').to_numpy(dtype=float)
            if (native+1e-6<total).any():raise ValueError('Native normalization denominator below retained counts')
            total=native
        if not np.isfinite(total).all() or (total<=0).any():raise ValueError('Invalid analysis denominator')
        return raw,corrected,total


def projected(matrix, positions, width):
    coo=matrix.tocoo()
    result=sparse.coo_matrix((coo.data,(coo.row,np.asarray(positions)[coo.col])),shape=(matrix.shape[0],width)).tocsr()
    result.sum_duplicates()
    return result


def normalized(matrix, denominator):
    out=(sparse.diags(10000./np.asarray(denominator))@matrix).tocsr()
    out.data=np.log1p(out.data)
    return out.astype(np.float32)


def load_contract(path, roles=None):
    config=json.loads(Path(path).read_text())
    records=[r for r in config['records'] if roles is None or r['role'] in roles]
    sources=[Source(r) for r in records]
    native={s.record['dataset']:s.identity_view() for s in sources}
    mappings,var,identities,diagnostics=build_contract(native)
    replacements={key:key.removeprefix('SYMBOL::') for key in var.index
                  if key.startswith('SYMBOL::AMBIGUOUS::')}
    if replacements:
        var.index=pd.Index([replacements.get(key,key) for key in var.index],name=var.index.name)
        for source in sources:
            name=source.record['dataset']
            mappings[name]=[replacements.get(key,key) for key in mappings[name]]
            for item in identities[name]:
                if item['canonical_feature'] in replacements:
                    item['canonical_feature']=replacements[item['canonical_feature']]
                    item['identity_status']='preserved_parent_ambiguous_feature'
            for column,key in enumerate(mappings[name]):
                if key.startswith('AMBIGUOUS::'):
                    var.loc[key,'gene_symbol']=str(source.var.iloc[column].get('gene_symbol',source.symbols[column]))
                    var.loc[key,'feature_identity_status']='ambiguous_source_feature'
        var=var.sort_index()
        diagnostics['parent_ambiguous_features_preserved']=len(replacements)
    lookup={name:i for i,name in enumerate(var.index)}
    for source in sources:
        source.canonical_genes=mappings[source.record['dataset']]
        source.global_columns=np.array([lookup[name] for name in source.canonical_genes],dtype=np.int32)
    return config,sources,var,identities,diagnostics


def read_registry(config, sources):
    obs=pd.read_csv(config['all_rows_ref'],sep='\t',dtype=str,keep_default_na=False,index_col=0)
    parts=[]
    for source in sources:
        frame=obs.loc[obs.dataset.eq(source.record['dataset'])].copy()
        if len(frame)!=len(source.rows):raise ValueError('Source registry/row selector length differs')
        expected=source.record['dataset']+'::'+source.obs.index[source.rows].astype(str)
        if not np.array_equal(frame.index.astype(str),expected):raise ValueError('Source registry order differs')
        parts.append(frame)
    result=pd.concat(parts)
    if not result.index.is_unique:raise ValueError('Global cell IDs repeat')
    return result


def save_json_new(path,value):
    with Path(path).open('x') as stream:json.dump(value,stream,indent=2,allow_nan=False)
