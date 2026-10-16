"""Independent bin-count aggregation of all frozen-role marker footprints.

Original sources: fresh reaggregation of accepted native unit caches, NOT a new
all-original-cell count extraction. Immune atlas: full fresh native count read.
"""
import json
from collections import Counter
from datetime import datetime,timezone

import numpy as np
import pandas as pd
import zarr
from scipy.sparse import load_npz,csr_matrix

from NSCLC_external_lineage_inputs_v2 import load,cell_arms,sha,matched_quota,VIEWS
from NSCLC_native_MFA_interpretation_inputs_v1 import PANELS
from external_functional_coverage_inputs_v1 import read_inputs


def main():
    data=load();root=data['root'];matched_quota(root);lib=data['lib'];n=len(lib)
    out=root/'NSCLC_native_MFA_external_role_marker_profiles_20261016_v2'
    terminal=out/'private_terminal.json';t=json.loads(terminal.read_text())
    assert t['status'].startswith('COMPLETE_ALL61_SOURCE6760_')
    for name,item in t['files'].items():assert sha(out/name)==item['sha256']
    plan=json.loads((out/'plan_before_any_external_role_marker_measurement.json').read_text())
    panels=['T_IDENTITY','CYTOTOXIC_TNK','TREG_IDENTITY','MYELOID_IDENTITY','B_IDENTITY','PLASMA_IDENTITY']
    markers=sorted(set('CD4 CD8A CD8B'.split())|{g for name in panels for g in PANELS[name]})
    assert markers==plan['markers'];ng=len(markers)
    with np.load(out/'private_fixed_marker_pool_sums.npz',allow_pickle=False) as z:
        saved={k:z[k] for k in z.files}
    assert saved['markers'].astype(str).tolist()==markers
    axes=read_inputs(root)['source_axes'];maps={};avail=np.zeros((n,ng),bool)
    for key,axis in axes.items():
        cols=[next(iter(axis['lookup'][g])) if g in axis['lookup'] and len(axis['lookup'][g])==1 else None for g in markers]
        count=Counter(c for c in cols if c is not None)
        cols=[c if c is not None and count[c]==1 else None for c in cols];maps[key]=cols
        ix=lib.index[lib.arm.eq(key[0])&lib.dataset.eq(key[1])]
        avail[ix]=np.array([c is not None for c in cols])
    assert np.array_equal(avail,saved['source_measured'])
    sizes=np.zeros(n*4,np.int64);logs=np.zeros((n*4,ng));positive=np.zeros((n*4,ng),np.int64)
    lookup={(r.arm,r.dataset,r.library_id):i for i,r in enumerate(lib.itertuples())}
    old=root.parents[2]/'Project_v3/annotation_rebuild/20261010_uniform_expression_annotation_v1'
    catalog=pd.read_csv(old/'all_node_number_catalog.tsv',sep='\t',keep_default_na=False).node_id.to_numpy(str)
    roles=data['catalogs']['frozen_original_annotation'].set_index('node_id').frozen_role
    cache_hashes={v['dataset']:v['manifest_sha256'] for v in t['cache_manifest_sha256']}
    assert len(cache_hashes)==52
    for source,expected in sorted(cache_hashes.items()):
        base=old/'by_gse'/source;assert sha(base/'manifest.json')==expected
        manifest=json.loads((base/'manifest.json').read_text())
        for f in ['units.tsv','genes.tsv','sum_log1p_CP10K.npz','positive_counts.npz']:
            assert sha(base/f)==manifest['output_sha256'][f]
        units=pd.read_csv(base/'units.tsv',sep='\t',keep_default_na=False)
        destination=np.full(len(units),-1,int)
        for r in units.itertuples():
            key=('frozen_original_annotation',source,r.sample_uid)
            node=catalog[r.node_index]
            if r.candidate_eligible and key in lookup and node in roles.index and roles.loc[node]>=0:
                destination[r.unit_row]=lookup[key]*4+int(roles.loc[node])
        take=np.flatnonzero(destination>=0)
        sizes+=np.bincount(destination[take],weights=units.n_cells.to_numpy()[take],minlength=n*4).astype(np.int64)
        columns=maps['frozen_original_annotation',source]
        for file,dest in [('sum_log1p_CP10K.npz',logs),('positive_counts.npz',positive)]:
            cache=load_npz(base/file).tocsr()
            for j,col in enumerate(columns):
                if col is None:continue
                value=cache[take,col].toarray().ravel()
                grouped=np.bincount(destination[take],weights=value,minlength=n*4)
                if file.startswith('positive'):
                    assert np.array_equal(grouped,grouped.astype(np.int64));grouped=grouped.astype(np.int64)
                dest[:,j]+=grouped
            del cache
    arm=next(a for a in cell_arms(data) if a['arm']=='archived_immune_atlas')
    z=zarr.open_group(str(data['atlas']),mode='r');raw=z['layers']['counts'];nr,nf=raw.attrs['shape']
    colmap=np.full((len(arm['source'].categories),nf),-1,np.int16)
    for si,source in enumerate(arm['source'].categories):
        for j,col in enumerate(maps['archived_immune_atlas',str(source)]):
            if col is not None:colmap[si,col]=j
    role=arm['roles'][0];selected_cells=0
    for start in range(0,nr,3072):
        stop=min(nr,start+3072);keep=role[start:stop]>=0
        if not keep.any():continue
        ptr=np.asarray(raw['indptr'][start:stop+1],np.int64);lo,hi=int(ptr[0]),int(ptr[-1])
        matrix=csr_matrix((np.asarray(raw['data'][lo:hi]),np.asarray(raw['indices'][lo:hi]),ptr-lo),shape=(stop-start,nf))[keep].copy()
        matrix.sum_duplicates();matrix.eliminate_zeros();matrix.sort_indices()
        ids=arm['library_index'][start:stop][keep]*4+role[start:stop][keep]
        source=arm['source'].codes[start:stop][keep]
        total=np.asarray(matrix.sum(axis=1)).ravel();assert (total>0).all()
        co=matrix.tocoo();m=colmap[source[co.row],co.col];valid=m>=0
        cells=co.row[valid];markers_index=m[valid]
        packed=ids[cells]*ng+markers_index
        values=np.log1p(co.data[valid].astype(float)*10000/total[cells])
        logs+=np.bincount(packed,weights=values,minlength=n*4*ng).reshape(n*4,ng)
        positive+=np.bincount(packed,minlength=n*4*ng).reshape(n*4,ng)
        sizes+=np.bincount(ids,minlength=n*4);selected_cells+=len(ids)
    assert np.array_equal(sizes,saved['n_cells']) and np.array_equal(positive,saved['positive_counts'])
    error=float(np.max(np.abs(logs-saved['sum_logCP10K'])))
    assert np.allclose(logs,saved['sum_logCP10K'],atol=1e-7,rtol=1e-10),error
    units=pd.read_csv(out/'private_all6760_frozen_role_units.tsv',sep='\t',keep_default_na=False)
    assert len(units)==n*4 and units.unit_index.tolist()==list(range(n*4))
    assert np.array_equal(units.n_cells,sizes)
    for i,r in enumerate(lib.itertuples()):
        u=units.iloc[i*4:i*4+4]
        assert u.arm.eq(r.arm).all() and u.dataset.eq(r.dataset).all() and u.library_id.eq(r.library_id).all()
        assert u['view'].tolist()==VIEWS
    table=pd.read_csv(out/'private_all_frozen_role_marker_measurements.tsv.gz',sep='\t',keep_default_na=False)
    assert len(table)==n*4*ng and not table.duplicated(['unit_index','gene']).any()
    for j,gene in enumerate(markers):
        rows=table[table.gene.eq(gene)].set_index('unit_index').loc[range(n*4)]
        measured=np.repeat(avail[:,j],4);nonzero=sizes>0;valid=measured&nonzero
        status=np.where(~measured,'GENE_UNMEASURED_OR_UNRESOLVED_NOT_NEGATIVE',np.where(nonzero,'MEASURED','CELL_POOL_ABSENT_NOT_NEGATIVE'))
        assert np.array_equal(status,rows.status.to_numpy())
        assert rows.loc[~valid,'mean_logCP10K'].eq('').all() and rows.loc[~valid,'positive_fraction'].eq('').all()
        assert np.allclose(rows.loc[valid,'mean_logCP10K'].astype(float),logs[valid,j]/sizes[valid],atol=1e-9,rtol=1e-9)
        assert np.allclose(rows.loc[valid,'positive_fraction'].astype(float),positive[valid,j]/sizes[valid],atol=1e-12)
    proof=out/'private_independent_all_role_marker_audit.json';assert not proof.exists()
    report=dict(status='PASS_ALL61_SOURCE6760_FROZEN_POOLS_ALL_FIXED_MARKER_VALUES_MISSINGNESS_AND_COUNTS',
        finished_UTC=datetime.now(timezone.utc).isoformat(),audit_source_sha256=sha(__file__),producer_terminal_sha256=sha(terminal),
        marker_value_rows_checked=len(table),maximum_sum_log_error=error,
        original_scope='ALL52 source accepted native unit caches independently reaggregated; not new all-original-cell raw replay',
        immune_scope='Fresh full atlas count stream, all frozen role-selected cells',immune_selected_cells=selected_cells,
        formal_identity_or_external_model_accepted=False,new_P=False)
    proof.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],finished_UTC=report['finished_UTC'],audit_sha256=sha(proof))))


if __name__=='__main__':main()
