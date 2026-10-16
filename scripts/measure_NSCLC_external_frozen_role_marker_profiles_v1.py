"""Fixed public identity-marker footprints for all source/frozen-pool units.

Existing original native caches are reused. Archived immune counts are streamed
once, retaining only the predefined marker panel and source-specific missingness.
No classifier, formal cell states, clinical effects or learned gene export.
"""
import json
from datetime import datetime,timezone

import numpy as np
import pandas as pd
import zarr
from scipy.sparse import csr_matrix,load_npz

from NSCLC_external_lineage_inputs_v2 import load,cell_arms,sha,matched_quota,VIEWS
from NSCLC_native_MFA_interpretation_inputs_v1 import PANELS
from external_functional_coverage_inputs_v1 import read_inputs

PANEL_NAMES=['T_IDENTITY','CYTOTOXIC_TNK','TREG_IDENTITY','MYELOID_IDENTITY','B_IDENTITY','PLASMA_IDENTITY']
MARKERS=sorted(set('CD4 CD8A CD8B'.split())|{g for name in PANEL_NAMES for g in PANELS[name]})


def main():
    data=load();root=data['root'];lib=data['lib'];quota=matched_quota(root)
    audit=root/'NSCLC_native_MFA_external_lineage_qualification_20261016_v2/private_independent_lineage_case_audit.json'
    assert sha(audit)=='b07bf0d9b682e8b6cc2ef72192a0ff110a722ba99bd7a931c5e7745c7640699a'
    out=root/'NSCLC_native_MFA_external_role_marker_profiles_20261016_v1'
    assert not out.exists();out.mkdir(mode=0o700)
    axes=read_inputs(root)['source_axes']
    old=root.parents[2]/'Project_v3/annotation_rebuild/20261010_uniform_expression_annotation_v1'
    original=pd.read_csv(old/'all_node_number_catalog.tsv',sep='\t',keep_default_na=False)
    catalogue=data['catalogs']['frozen_original_annotation'].set_index('node_id')
    selected=data['counts'];lookup={(r.arm,r.dataset,r.library_id):i for i,r in enumerate(lib.itertuples())}
    n=len(lib);ng=len(MARKERS)
    with np.load(audit.parent/'private_all1690_three_channel_role_counts.npz',allow_pickle=False) as f:
        expected=f['counts'][:,0,:].reshape(-1)
    nc=np.zeros(n*4,np.int64);logs=np.zeros((n*4,ng));positive=np.zeros((n*4,ng),np.int64)
    available=np.zeros((n,ng),bool)
    physical_maps={}
    for key,axis in axes.items():
        positions=[next(iter(axis['lookup'][g])) if g in axis['lookup'] and len(axis['lookup'][g])==1 else -1 for g in MARKERS]
        multiplicity=pd.Series([c for c in positions if c>=0]).value_counts().to_dict()
        positions=np.array([c if c>=0 and multiplicity[c]==1 else -1 for c in positions],int)
        physical_maps[key]=positions
        ix=lib.index[(lib.arm==key[0])&(lib.dataset==key[1])]
        available[ix]=positions>=0
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),quota=quota,
        upstream_audit_sha256=sha(audit),input_metadata_sha256=data['hashes'],
        scope='All61 source-arm entries/1690 libraries x4 frozen candidate roles, including zero/small units; no source or role selected by expression.',
        markers=MARKERS,panels={name:PANELS[name] for name in PANEL_NAMES},
        original='Reuse immutable source-native unit sum_log1p_CP10K and positive_counts caches; exact current eligible sample/node counts checked.',
        immune='Read raw layers/counts with original frozen metadata rows; canonicalize working sparse blocks only. Full physical RNA per-cell denominator, natural log1p(CP10K), sum by disjoint frozen role/library. Unique source-measured marker mappings only.',
        output='Marker means and positive-cell fractions; measured biological zeros distinct from unavailable mapping and absent cell pool. Every6760 unit retained. Descriptive identities, no P or automatic annotation/relabeling.',
        limits='Markers are prespecified lineage checks, not independent proof of subtype specificity; all original annotations and factor weights remain frozen.',
        prospective_audit='Exact ALL original cache-unit count/mapping checks plus fresh independent marker aggregation required before external effect models. Production completion alone is not biological validation.',
        external_factor_or_clinical_inference=False,new_download=False)
    (out/'plan_before_any_external_role_marker_measurement.json').write_text(json.dumps(plan,indent=2)+'\n')
    cache_proofs=[]
    for (arm,source),axis in sorted(axes.items()):
        if arm!='frozen_original_annotation':continue
        base=old/'by_gse'/source;manifest=json.loads((base/'manifest.json').read_text())
        assert manifest['status']=='COMPLETE_EXACT_MEMBER_EXPRESSION_EVIDENCE'
        files=['units.tsv','genes.tsv','sum_log1p_CP10K.npz','positive_counts.npz']
        for f in files:assert sha(base/f)==manifest['output_sha256'][f]
        units=pd.read_csv(base/'units.tsv',sep='\t',keep_default_na=False)
        genes=pd.read_csv(base/'genes.tsv',sep='\t').gene.astype(str).to_numpy()
        assert units.unit_row.tolist()==list(range(len(units)))
        assert all((c<0 or genes[c]==g) for c,g in zip(physical_maps[arm,source],MARKERS))
        assert pd.api.types.is_bool_dtype(units.candidate_eligible)
        units['node_id']=original.node_id.to_numpy()[units.node_index.to_numpy(int)]
        pick=units.candidate_eligible & units.sample_uid.isin(lib.loc[(lib.arm==arm)&(lib.dataset==source),'library_id'])
        u=units[pick].copy()
        assert u.node_id.isin(catalogue.index).all()
        inherited=selected[(selected.arm==arm)&(selected.dataset==source)].set_index(['library_id','node_id']).n_cells
        current=u.set_index(['sample_uid','node_id']).n_cells
        assert current.index.is_unique and set(current.index)==set(inherited.index)
        assert np.array_equal(current.to_numpy(),inherited.loc[current.index].to_numpy())
        u['role']=catalogue.loc[u.node_id,'frozen_role'].to_numpy(int);u=u[u.role>=0]
        rows=np.array([lookup[(arm,source,r.sample_uid)]*4+r.role for r in u.itertuples()],int)
        group=csr_matrix((np.ones(len(u),dtype=np.int64),(rows,u.unit_row.to_numpy(int))),shape=(n*4,len(units)))
        nc+=np.asarray(group@units.n_cells.to_numpy()).ravel().astype(np.int64)
        positions=physical_maps[arm,source];valid=np.flatnonzero(positions>=0)
        # Limit dense arrays to fixed public markers before group projection.
        for file,dest in [('sum_log1p_CP10K.npz',logs),('positive_counts.npz',positive)]:
            raw=load_npz(base/file).tocsr();assert raw.shape==(len(units),len(genes))
            values=(group@raw[:,positions[valid]]).toarray()
            if file.startswith('positive'):assert np.array_equal(values,values.astype(np.int64))
            dest[:,valid]+=values
            del raw,values
        cache_proofs.append(dict(dataset=source,manifest_sha256=sha(base/'manifest.json')))
        (out/'private_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),stage='original_cached_markers',completed_source=source,completed_sources=len(cache_proofs)))+'\n')
    assert len(cache_proofs)==52
    # Original cached counts must match frozen metadata before streaming atlas.
    oi=lib.index[lib.arm.eq('frozen_original_annotation')].to_numpy()
    assert np.array_equal(nc.reshape(n,4)[oi],expected.reshape(n,4)[oi])
    arm=next(a for a in cell_arms(data) if a['arm']=='archived_immune_atlas')
    z=zarr.open_group(str(data['atlas']),mode='r');raw=z['layers']['counts'];nr,nf=raw.attrs['shape'];assert nr==len(arm['node'])
    physical=sorted({int(c) for (a,s),cols in physical_maps.items() if a=='archived_immune_atlas' for c in cols if c>=0})
    pi={c:i for i,c in enumerate(physical)}
    panel_sum=np.zeros((n*4,len(physical)));panel_positive=np.zeros((n*4,len(physical)),np.int64)
    role=arm['roles'][0];libidx=arm['library_index'];unit=libidx*4+role
    for start in range(0,nr,4096):
        stop=min(nr,start+4096);take=role[start:stop]>=0
        if take.any():
            ptr=np.asarray(raw['indptr'][start:stop+1],np.int64);lo,hi=int(ptr[0]),int(ptr[-1])
            x=csr_matrix((np.asarray(raw['data'][lo:hi]),np.asarray(raw['indices'][lo:hi]),ptr-lo),shape=(stop-start,nf))[take].copy()
            x.sum_duplicates();x.eliminate_zeros();x.sort_indices()
            assert np.isfinite(x.data).all() and (x.data>0).all() and np.equal(x.data,np.floor(x.data)).all()
            total=np.asarray(x.sum(axis=1)).ravel();assert (total>0).all()
            ids=unit[start:stop][take]
            nc+=np.bincount(ids,minlength=n*4)
            group=csr_matrix((np.ones(len(ids)),(ids,np.arange(len(ids)))),shape=(n*4,len(ids)))
            block=x[:,physical].astype(float)
            binary=block.copy();binary.data[:]=1
            values=(group@binary).toarray();assert np.equal(values,np.floor(values)).all();panel_positive+=values.astype(np.int64)
            block.data=np.log1p(block.data*np.repeat(10000/total,np.diff(block.indptr)))
            panel_sum+=(group@block).toarray()
        if start//4096%32==0 or stop==nr:
            (out/'private_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),stage='atlas_raw_marker_stream',rows_read=stop,total_rows=nr))+'\n')
    for (a,source),positions in physical_maps.items():
        if a!='archived_immune_atlas':continue
        libs=lib.index[(lib.arm==a)&(lib.dataset==source)].to_numpy();rows=(libs[:,None]*4+np.arange(4)).reshape(-1)
        for j,c in enumerate(positions):
            if c>=0:logs[rows,j]=panel_sum[rows,pi[c]];positive[rows,j]=panel_positive[rows,pi[c]]
    assert np.array_equal(nc,expected)
    assert (positive<=nc[:,None]).all() and (positive>=0).all() and np.isfinite(logs).all() and (logs>=0).all()
    unit_rows=[];records=[]
    for i,r in enumerate(lib.itertuples()):
        for vi,view in enumerate(VIEWS):
            ix=i*4+vi;count=int(nc[ix]);unit_rows.append(dict(unit_index=ix,arm=r.arm,dataset=r.dataset,library_id=r.library_id,view=view,n_cells=count,n20=count>=20))
            for j,gene in enumerate(MARKERS):
                status='MEASURED' if available[i,j] and count else 'CELL_POOL_ABSENT_NOT_NEGATIVE' if available[i,j] else 'GENE_UNMEASURED_OR_UNRESOLVED_NOT_NEGATIVE'
                records.append(dict(unit_index=ix,gene=gene,status=status,mean_logCP10K=float(logs[ix,j]/count) if status=='MEASURED' else np.nan,
                    positive_fraction=float(positive[ix,j]/count) if status=='MEASURED' else np.nan))
    pd.DataFrame(unit_rows).to_csv(out/'private_all6760_frozen_role_units.tsv',sep='\t',index=False)
    pd.DataFrame(records).to_csv(out/'private_all_frozen_role_marker_measurements.tsv.gz',sep='\t',index=False,compression='gzip')
    np.savez_compressed(out/'private_fixed_marker_pool_sums.npz',n_cells=nc,sum_logCP10K=logs,positive_counts=positive,source_measured=available,markers=np.array(MARKERS))
    terminal=dict(status='COMPLETE_ALL61_SOURCE6760_FROZEN_POOL_FIXED_MARKER_PROFILES_AUDIT_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),
        source_sha256=sha(__file__),cache_manifest_sha256=cache_proofs,files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()},
        new_external_factor_scores_or_P=False,annotations_changed=False)
    assert sum(v['bytes'] for v in terminal['files'].values())<268435456
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'],finished_UTC=terminal['finished_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
