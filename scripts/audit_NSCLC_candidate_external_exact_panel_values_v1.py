"""All pool invariants, all52 original caches, and fixed96 raw-pool replay.

The raw replay scope is96 pools, not all6760. No associations or clinical values.
"""
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, load_npz

from NSCLC_candidate_external_exact_RNA_inputs_v1 import inputs,raw_arms,OUT
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def compare(actual, expected, measured, exact=False):
    assert np.isnan(actual[~measured]).all()
    if exact:assert np.array_equal(actual[measured],expected[measured])
    else:assert np.allclose(actual[measured],expected[measured],atol=1e-6,rtol=1e-9)
    return float(np.max(abs(actual[measured]-expected[measured]))) if measured.any() else 0.


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);resource(root,1024**2)
    job=root/'HNSCC_receiver_followup_job_control_20261016_v2/extract_NSCLC_candidate_external_exact_RNA_v1'
    launch=json.loads((job/'launch_receipt.json').read_text());out=root/OUT
    assert launch['source_sha256']=='33d7d8cbdd0aa6cf009ce6b54cca0d0b4f04de4aed9af5c8763ebd77911c7a96'
    deadline=time.monotonic()+21600
    while not (out/'private_terminal.json').exists():
        p=Path('/proc')/str(launch['pid'])/'stat';assert p.exists(),'Bound producer missing without terminal'
        fields=p.read_text().rsplit(')',1)[1].split()
        assert fields[0]!='Z' and fields[19]==str(launch['proc_start_ticks'])
        assert time.monotonic()<deadline,'Observation expired, no restart or failure inferred'
        time.sleep(3)
    terminal=json.loads((out/'private_terminal.json').read_text())
    assert terminal['source_sha256']==launch['source_sha256']
    for name,item in terminal['files'].items():assert sha(out/name)==item['sha256']
    available=int(next(s for s in Path('/proc/meminfo').read_text().splitlines() if s.startswith('MemAvailable:')).split()[1])*1024
    assert available>3*1024**3
    state=inputs();plan=json.loads((out/'plan_before_any_external_exact_panel_expression.json').read_text())
    assert plan['frozen_inputs']==state['hashes']==terminal['frozen_inputs']
    assert terminal['helper_sha256']==sha(Path(__file__).with_name('NSCLC_candidate_external_exact_RNA_inputs_v1.py'))
    with np.load(out/'private_all_frozen_pool_exact_panel_measurements.npz',allow_pickle=False) as f:z={k:f[k] for k in f.files}
    units=pd.read_csv(out/'private_all6760_exact_panel_unit_axis.tsv.gz',sep='\t',keep_default_na=False)
    pd.testing.assert_frame_equal(units,state['units'])
    nu=len(units);ng=len(state['genes']);nrec=len(state['complexes'])
    assert np.array_equal(z['unit_index'],np.arange(nu)) and np.array_equal(z['genes'],state['genes'])
    assert z['complexes'].tolist()==['_'.join(v) for v in state['complexes']]
    assert np.array_equal(z['measured'],state['measured']) and np.array_equal(z['complex_measured'],state['complex_measured'])
    assert np.array_equal(z['n_cells'],units.n_cells) and np.array_equal(z['full_RNA_sum'],units.full_source_RNA_sum)
    assert np.array_equal(z['zero_RNA_cells'],units.zero_full_RNA_cells)
    for field,mask in [('sum_logCP10K',z['measured']),('positive_counts',z['measured']),
                       ('receptor_positive',z['complex_measured']),('receptor_positive2',z['complex_measured'])]:
        a=z[field];assert np.isnan(a[~mask]).all() and np.isfinite(a[mask]).all() and (a[mask]>=0).all()
        if field!='sum_logCP10K':assert np.equal(a[mask],np.floor(a[mask])).all()
    assert ((z['positive_counts']<=z['n_cells'][:,None])|~z['measured']).all()
    assert ((z['receptor_positive2']<=z['receptor_positive'])|~z['complex_measured']).all()
    ccols=[[list(state['genes']).index(g) for g in rec] for rec in state['complexes']]
    for j,cols in enumerate(ccols):
        take=z['complex_measured'][:,j]
        assert (z['receptor_positive'][take,j]<=z['positive_counts'][take][:,cols].min(axis=1)).all()
    for derived,numerator in [('mean_logCP10K','sum_logCP10K'),('positive_fraction','positive_counts'),('receptor_positive_fraction','receptor_positive')]:
        expected=np.divide(z[numerator],z['n_cells'][:,None],out=np.full_like(z[numerator],np.nan),where=z['n_cells'][:,None]>0)
        assert np.allclose(z[derived],expected,atol=1e-12,rtol=1e-12,equal_nan=True)
    # Independent exact aggregation from all52 accepted original full-gene caches.
    old=state['root'].parents[2]/'Project_v3/annotation_rebuild/20261010_uniform_expression_annotation_v1'
    original=pd.read_csv(old/'all_node_number_catalog.tsv',sep='\t',keep_default_na=False)
    catalog=state['data']['catalogs']['frozen_original_annotation'].set_index('node_id')
    lib=state['data']['lib'];lookup={(r.arm,r.dataset,r.library_id):i for i,r in enumerate(lib.itertuples())}
    cached_log=np.zeros((nu,ng));cached_positive=np.zeros((nu,ng),np.int64);cached_n=np.zeros(nu,np.int64)
    cache_proofs=[]
    for info in state['infos']:
        if info['arm']!='frozen_original_annotation':continue
        base=old/'by_gse'/info['dataset'];manifest=json.loads((base/'manifest.json').read_text())
        assert manifest['status']=='COMPLETE_EXACT_MEMBER_EXPRESSION_EVIDENCE'
        for name in ['units.tsv','genes.tsv','sum_log1p_CP10K.npz','positive_counts.npz']:
            assert sha(base/name)==manifest['output_sha256'][name]
        u=pd.read_csv(base/'units.tsv',sep='\t',keep_default_na=False)
        axis=pd.read_csv(base/'genes.tsv',sep='\t',keep_default_na=False).gene.to_numpy(str)
        assert u.unit_row.tolist()==list(range(len(u)))
        u['node_id']=original.node_id.to_numpy()[u.node_index.to_numpy(int)]
        libraries=set(lib.loc[lib.arm.eq(info['arm'])&lib.dataset.eq(info['dataset']),'library_id'])
        pick=u.candidate_eligible&u.sample_uid.isin(libraries)
        u=u[pick].copy();assert u.node_id.isin(catalog.index).all()
        inherited=state['data']['counts'];inherited=inherited[inherited.arm.eq(info['arm'])&inherited.dataset.eq(info['dataset'])].set_index(['library_id','node_id']).n_cells
        current=u.set_index(['sample_uid','node_id']).n_cells
        assert current.index.is_unique and set(current.index)==set(inherited.index)
        assert np.array_equal(current.to_numpy(),inherited.loc[current.index].to_numpy())
        u['role']=catalog.loc[u.node_id,'frozen_role'].to_numpy(int);u=u[u.role>=0]
        dest=np.array([lookup[(info['arm'],info['dataset'],r.sample_uid)]*4+r.role for r in u.itertuples()],int)
        np.add.at(cached_n,dest,u.n_cells.to_numpy(np.int64))
        positions=info['source_columns'];valid=np.flatnonzero(positions>=0)
        assert np.array_equal(axis[positions[valid]],state['genes'][valid])
        for name,result in [('sum_log1p_CP10K.npz',cached_log),('positive_counts.npz',cached_positive)]:
            a=load_npz(base/name).tocsr();assert a.shape==(len(pd.read_csv(base/'units.tsv',sep='\t',usecols=['unit_row'])),len(axis))
            values=a[u.unit_row.to_numpy(int)][:,positions[valid]].toarray()
            if name=='positive_counts.npz':
                assert np.array_equal(values,np.rint(values));values=values.astype(np.int64)
            np.add.at(result,(dest[:,None],valid[None,:]),values)
        cache_proofs.append(dict(dataset=info['dataset'],manifest_sha256=sha(base/'manifest.json')))
    assert len(cache_proofs)==52
    oi=units.arm.eq('frozen_original_annotation').to_numpy()
    assert np.array_equal(cached_n[oi],z['n_cells'][oi])
    cache_log_error=compare(z['sum_logCP10K'][oi],cached_log[oi],z['measured'][oi])
    compare(z['positive_counts'][oi],cached_positive[oi],z['measured'][oi],exact=True)
    cached_mean=np.divide(np.where(z['measured'][oi],cached_log[oi],np.nan),cached_n[oi,None],out=np.full_like(cached_log[oi],np.nan),where=cached_n[oi,None]>0)
    assert np.allclose(cached_mean,z['mean_logCP10K'][oi],atol=1e-9,rtol=1e-9,equal_nan=True)
    # Fresh all-panel/all-complex raw replay of the ORIGINAL fixed96 pools only.
    chosen=state['selection'].unit_index.to_numpy(int);inverse=np.full(nu,-1,int);inverse[chosen]=np.arange(96)
    audit_log=np.zeros((96,ng));audit_positive=np.zeros((96,ng),np.int64)
    audit_joint=np.zeros((96,nrec),np.int64);audit_joint2=np.zeros_like(audit_joint)
    audit_n=np.zeros(96,np.int64);audit_RNA=np.zeros(96,np.int64);audit_zero=np.zeros(96,np.int64)
    checkpoint=root/'HNSCC_receiver_followup_job_control_20261016_v2'/Path(__file__).stem/'private_audit_checkpoint.json'
    for arm,raw,unit,infos in raw_arms(state):
        destination=np.full(len(unit),-1,int);captured=unit>=0;destination[captured]=inverse[unit[captured]]
        nr,nf=map(int,raw.attrs['shape'])
        for start in range(0,nr,3072):
            stop=min(start+3072,nr);take=destination[start:stop]>=0
            if take.any():
                ptr=np.asarray(raw['indptr'][start:stop+1],np.int64);lo,hi=int(ptr[0]),int(ptr[-1])
                x=csr_matrix((np.asarray(raw['data'][lo:hi]),np.asarray(raw['indices'][lo:hi]),ptr-lo),shape=(stop-start,nf))[take].copy()
                x.sum_duplicates();x.eliminate_zeros();x.sort_indices()
                assert np.isfinite(x.data).all() and (x.data>0).all() and np.equal(x.data,np.floor(x.data)).all()
                ids=destination[start:stop][take];sources=state['unit_source'][chosen[ids]]
                for si in np.unique(sources):
                    keep=sources==si;m=x[keep];local=ids[keep];info=infos[int(si)]
                    assert info['present'][m.indices].all()
                    total=np.asarray(m.sum(axis=1)).ravel();ok=info['physical']>=0
                    panel=np.zeros((len(local),ng),np.int64);panel[:,ok]=m[:,info['physical'][ok]].toarray().astype(np.int64)
                    ratio=np.divide(panel*10000.,total[:,None],out=np.zeros_like(panel,dtype=float),where=total[:,None]>0)
                    logs=np.log1p(ratio)
                    for uid in np.unique(local):
                        rr=local==uid;a=panel[rr]
                        audit_n[uid]+=int(rr.sum());audit_RNA[uid]+=int(total[rr].sum());audit_zero[uid]+=int((total[rr]==0).sum())
                        audit_log[uid]+=logs[rr].sum(axis=0);audit_positive[uid]+=(a>0).sum(axis=0)
                        for j,cols in enumerate(ccols):
                            if not ok[cols].all():continue
                            minimum=a[:,cols].min(axis=1)
                            audit_joint[uid,j]+=int(np.count_nonzero(minimum>=1));audit_joint2[uid,j]+=int(np.count_nonzero(minimum>=2))
            if start//3072%64==0 or stop==nr:
                checkpoint.write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),stage='fixed96_exact_panel_raw_replay',arm=arm,rows_considered=stop,total_rows=nr,raw_cells_accumulated=int(audit_n.sum())))+'\n')
    assert np.array_equal(audit_n,z['n_cells'][chosen]) and np.array_equal(audit_RNA,z['full_RNA_sum'][chosen]) and np.array_equal(audit_zero,z['zero_RNA_cells'][chosen])
    raw_log_error=compare(z['sum_logCP10K'][chosen],audit_log,z['measured'][chosen])
    compare(z['positive_counts'][chosen],audit_positive,z['measured'][chosen],exact=True)
    compare(z['receptor_positive'][chosen],audit_joint,z['complex_measured'][chosen],exact=True)
    compare(z['receptor_positive2'][chosen],audit_joint2,z['complex_measured'][chosen],exact=True)
    raw_mean=np.divide(np.where(z['measured'][chosen],audit_log,np.nan),audit_n[:,None],out=np.full_like(audit_log,np.nan),where=audit_n[:,None]>0)
    assert np.allclose(raw_mean,z['mean_logCP10K'][chosen],atol=1e-9,rtol=1e-9,equal_nan=True)
    report=dict(status='PASS_ALL6760_AXES_INVARIANTS52_CACHE_PANELS_AND96_RAW_PANELS_COMPLETE_RECEPTORS',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),
        producer_terminal_sha256=sha(out/'private_terminal.json'),frozen_inputs=state['hashes'],
        original52_cache_manifest_proofs=cache_proofs,maximum_cached_sum_log_error=cache_log_error,
        maximum_fixed96_raw_sum_log_error=raw_log_error,raw_pools_independently_checked=96,
        raw_cells_in96_pools=int(audit_n.sum()),all6760_pools_independently_raw_replayed=False,
        limitations='Shared accepted row/role/source mapping I/O. Full pool axes/invariants and all52 original-panel caches checked; independent raw aggregation covers only predeclared96 pools. No effect, clinical or biological replication claim.')
    proof=out/'private_independent_exact_panel_audit.json';assert not proof.exists()
    proof.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(status=report['status'],completed_UTC=report['completed_UTC'],proof_sha256=sha(proof))))


if __name__=='__main__':main()
