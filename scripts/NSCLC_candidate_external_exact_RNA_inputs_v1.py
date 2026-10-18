"""Private frozen panel, source axes and row/role I/O for exact RNA extraction."""
import json
import os
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem

from NSCLC_external_lineage_inputs_v2 import load, cell_arms, VIEWS
from external_functional_coverage_inputs_v1 import read_inputs
from describe_NSCLC_native_interface_mechanistic_followup_v1 import read_scope, SCOPE_SHA
from NSCLC_native_MFA_interpretation_inputs_v1 import sha

OUT = 'NSCLC_candidate_external_exact_RNA_20261018_v1'


def inputs():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    route=root/'NSCLC_candidate_external_exact_RNA_preflight_20261018_v1'
    proof=route/'private_independent_exact_RNA_metadata_audit.json'
    assert sha(proof)=='638443f53df36f6741a86ebb5d98bacb6145e6c363f29f4a4f7553b82c219868'
    assert json.loads(proof.read_text())['producer_terminal_sha256']==sha(route/'private_terminal.json')
    manifest=json.loads((route/'private_terminal.json').read_text())
    for name,item in manifest['files'].items():assert sha(route/name)==item['sha256']
    scope,candidates=read_scope(dict(root=root))
    genes=np.array(sorted(set(scope['ligand_specs'])|{g for s in scope['ligand_specs'].values() for c in s['complexes'] for g in c}))
    complexes=sorted({tuple(c) for s in scope['ligand_specs'].values() for c in s['complexes']})
    data=load(); ref=read_inputs(root);lib=data['lib'];nu=len(lib)*4
    old=root/'NSCLC_native_MFA_external_frozen_role_count_inputs_20261016_v2'
    proof=old/'private_independent_native_count_input_audit.json'
    assert sha(proof)=='7bd028cc1c89ac7324ac25c65042b5b184cb377162a38f9544caf7bb24048c04'
    assert json.loads(proof.read_text())['producer_terminal_sha256']==sha(old/'private_terminal.json')
    terminal=json.loads((old/'private_terminal.json').read_text())
    for name in ['private_all6760_native_count_units_and_quality.tsv.gz','private_prespecified96_raw_count_audit_units.tsv.gz','plan_before_all_external_native_count_extraction.json']:
        assert sha(old/name)==terminal['files'][name]['sha256']
    units=pd.read_csv(old/'private_all6760_native_count_units_and_quality.tsv.gz',sep='\t',keep_default_na=False)
    selection=pd.read_csv(old/'private_prespecified96_raw_count_audit_units.tsv.gz',sep='\t',keep_default_na=False)
    assert units.unit_index.tolist()==list(range(nu)) and len(selection)==96 and selection.unit_index.is_unique
    assert units['view'].tolist()==list(np.tile(VIEWS,len(lib)))
    assert np.array_equal(units.library_id.astype(str),np.repeat(lib.library_id.astype(str).to_numpy(),4))
    oldplan=json.loads((old/'plan_before_all_external_native_count_extraction.json').read_text())
    assert sha(old/'private_prespecified96_raw_count_audit_units.tsv.gz')==oldplan['predetermined_raw_audit_units_sha256']
    assert len(selection[['arm','dataset']].drop_duplicates())==61
    routes=pd.read_csv(route/'private_all_requested_gene_mappings.tsv',sep='\t')
    infos=[];unit_source=np.full(nu,-1,int);measured=np.zeros((nu,len(genes)),bool)
    complex_measured=np.zeros((nu,len(complexes)),bool)
    for si,((arm,source),axis) in enumerate(sorted(ref['source_axes'].items())):
        positions=[next(iter(axis['lookup'][g])) if g in axis['lookup'] and len(axis['lookup'][g])==1 else -1 for g in genes]
        multiplicity=Counter(p for p in positions if p>=0)
        positions=np.array([p if p>=0 and multiplicity[p]==1 else -1 for p in positions],int)
        table=routes[routes.arm.eq(arm)&routes.dataset.eq(source)].set_index('gene')
        assert np.array_equal(positions,table.loc[genes,'source_column'].to_numpy(int))
        libraries=lib.index[lib.arm.eq(arm)&lib.dataset.eq(source)].to_numpy(int)
        rows=(libraries[:,None]*4+np.arange(4)).reshape(-1);unit_source[rows]=si
        measured[rows]=positions>=0
        for ci,rec in enumerate(complexes):complex_measured[rows,ci]=all(positions[list(genes).index(g)]>=0 for g in rec)
        infos.append(dict(arm=arm,dataset=source,source_columns=positions,global_units=rows))
    assert len(infos)==61 and (unit_source>=0).all()
    return dict(root=root,data=data,ref=ref,units=units,selection=selection,genes=genes,complexes=complexes,
                infos=infos,unit_source=unit_source,measured=measured,complex_measured=complex_measured,
                hashes=dict(scope=SCOPE_SHA,route_audit=sha(route/'private_independent_exact_RNA_metadata_audit.json'),
                    old_count_audit=sha(proof),units=sha(old/'private_all6760_native_count_units_and_quality.tsv.gz'),
                    prespecified96=sha(old/'private_prespecified96_raw_count_audit_units.tsv.gz'),
                    metadata=data['hashes'],gene_axes=ref['input_sha256']))


def raw_arms(state):
    data=state['data'];oldbase=state['root'].parents[2]/'Project_v3/annotation_rebuild/20261010_uniform_expression_annotation_v1'
    for arm in cell_arms(data):
        original=arm['arm']=='frozen_original_annotation'
        z=zarr.open_group(str(data['old_zarr'] if original else data['atlas']),mode='r')
        raw=z['X'] if original else z['layers']['counts'];nr,nf=map(int,raw.attrs['shape'])
        unit=np.full(nr,-1,np.int32);roles=arm['roles'][0];take=roles>=0
        rr=np.arange(len(roles))
        if original:
            with np.load(data['old']/'private_canonical_row_selection.npz',allow_pickle=False) as f:rr=f['global_rows']
        assert len(rr)==len(roles) and len(np.unique(rr))==len(rr)
        unit[rr[take]]=arm['library_index'][take]*4+roles[take]
        if not original:
            source_order=np.asarray(read_elem(z['uns']['source_order'])).astype(str).tolist()
            source_mask=np.asarray(read_elem(z['varm']['measured_in_source']),bool)
        infos={}
        for si,info in enumerate(state['infos']):
            if info['arm']!=arm['arm']:continue
            if original:
                path=oldbase/'by_gse'/info['dataset']/'genes.tsv'
                assert sha(path)==state['ref']['legacy_gene_axis_sha256'][info['dataset']]
                g=pd.read_csv(path,sep='\t',keep_default_na=False);columns=g.core_column.to_numpy(int)
                present=np.zeros(nf,bool);present[columns]=True
                physical=np.array([columns[c] if c>=0 else -1 for c in info['source_columns']],int)
            else:
                physical=info['source_columns'];present=source_mask[:,source_order.index(info['dataset'])]
            assert len(set(physical[physical>=0]))==int((physical>=0).sum())
            assert present[physical[physical>=0]].all()
            infos[si]=dict(physical=physical,present=present)
        yield arm['arm'],raw,unit,infos
