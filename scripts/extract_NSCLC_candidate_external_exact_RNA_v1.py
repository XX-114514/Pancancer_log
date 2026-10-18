"""Exact fixed-panel mean-log RNA and same-cell receptors for all existing pools."""
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

from NSCLC_candidate_external_exact_RNA_inputs_v1 import inputs,raw_arms,OUT
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    state=inputs();root=state['root'];budget=resource(root,32*1024**2)
    available=int(next(s for s in Path('/proc/meminfo').read_text().splitlines() if s.startswith('MemAvailable:')).split()[1])*1024
    assert available>3*1024**3
    out=root/OUT;assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),
        helper_sha256=sha(Path(__file__).with_name('NSCLC_candidate_external_exact_RNA_inputs_v1.py')),
        frozen_inputs=state['hashes'],budget=budget,memory_available_bytes=available,minimum_free_RAM_bytes=3*1024**3,
        workers=1,BLAS_threads=1,raw_chunk_rows=4096,
        scope='All original candidates and curated functional complexes across all61 sources1690 libraries6760 frozen pools, including empty/small/unknown-context units. No expression selection or outcome read.',
        normalization='Original physical full-source cell total; true cell-weighted mean(log1p(CP10K)), not log of summed pseudobulk. Preserve full denominator, literal source mapping and frozen row/role selection.',
        complexes='All components count>=1 in SAME cell main; all>=2 sensitivity. Every predefined complex retained.',
        missing='Unavailable genes/complexes are NaN with measured masks; measured biological zeros retained. Empty captured pools have zero cell counts and undefined means/fractions. Zero full-RNA cells retained and explicitly counted.',
        QC='Exact all6760 counts/fullRNA/zeroRNA vs accepted prior counts; original full-RNA/log-gene means also checked.',
        audit='All52 original unit mean-log/positive caches for the whole selected panel plus original metadata-selected96 raw pools, all panel genes and complete complexes. This is NOT all6760-pool independent raw replay.',
        no_associations_P_new_downloads_or_annotation_changes=True)
    (out/'plan_before_any_external_exact_panel_expression.json').write_text(json.dumps(plan,indent=2)+'\n')
    nu=len(state['units']);ng=len(state['genes']);nr=len(state['complexes'])
    sums=np.zeros((nu,ng));positive=np.zeros((nu,ng),np.int64)
    joint=np.zeros((nu,nr),np.int64);joint2=np.zeros_like(joint)
    nc=np.zeros(nu,np.int64);RNA=np.zeros(nu,np.int64);zero=np.zeros(nu,np.int64)
    logRNA=np.zeros(nu);loggenes=np.zeros(nu)
    ci=[np.array([list(state['genes']).index(g) for g in rec],int) for rec in state['complexes']]
    for arm,raw,unit,infos in raw_arms(state):
        nrow,ncol=map(int,raw.attrs['shape'])
        for start in range(0,nrow,4096):
            stop=min(start+4096,nrow);valid=unit[start:stop]>=0
            if valid.any():
                ptr=np.asarray(raw['indptr'][start:stop+1],np.int64);lo,hi=int(ptr[0]),int(ptr[-1])
                x=csr_matrix((np.asarray(raw['data'][lo:hi]),np.asarray(raw['indices'][lo:hi]),ptr-lo),shape=(stop-start,ncol))[valid].copy()
                x.sum_duplicates();x.eliminate_zeros();x.sort_indices()
                assert np.isfinite(x.data).all() and (x.data>0).all() and np.equal(x.data,np.floor(x.data)).all()
                ids=unit[start:stop][valid];sources=state['unit_source'][ids]
                for si in np.unique(sources):
                    take=sources==si;m=x[take];dest=ids[take];info=infos[int(si)]
                    assert info['present'][m.indices].all()
                    total=np.asarray(m.sum(axis=1)).ravel();assert (total>=0).all() and (total<2**53).all()
                    panel=np.zeros((len(dest),ng),np.int64);ok=info['physical']>=0
                    panel[:,ok]=m[:,info['physical'][ok]].toarray().astype(np.int64)
                    assert (panel>=0).all() and (panel<=total[:,None]).all()
                    scale=np.divide(10000.,total,out=np.zeros(len(total)),where=total>0)
                    transformed=np.log1p(panel*scale[:,None])
                    group=csr_matrix((np.ones(len(dest)),(dest,np.arange(len(dest)))),shape=(nu,len(dest)))
                    sums+=group@transformed
                    positive+=(group@(panel>0)).astype(np.int64)
                    nc+=np.bincount(dest,minlength=nu)
                    np.add.at(RNA,dest,total.astype(np.int64));np.add.at(zero,dest,total==0)
                    np.add.at(logRNA,dest,np.log1p(total));np.add.at(loggenes,dest,np.log1p(np.diff(m.indptr)))
                    for j,cols in enumerate(ci):
                        if not ok[cols].all():continue
                        hit=(panel[:,cols]>=1).all(axis=1);hit2=(panel[:,cols]>=2).all(axis=1)
                        joint[:,j]+=np.bincount(dest[hit],minlength=nu)
                        joint2[:,j]+=np.bincount(dest[hit2],minlength=nu)
            if start//4096%32==0 or stop==nrow:
                checkpoint=dict(observed_UTC=datetime.now(timezone.utc).isoformat(),stage='exact_panel_raw_stream',arm=arm,rows_read=stop,total_rows=nrow,selected_cells_accumulated=int(nc.sum()))
                (out/'private_checkpoint.json').write_text(json.dumps(checkpoint)+'\n')
        resource(root,32*1024**2)
    expected=state['units'];assert np.array_equal(nc,expected.n_cells)
    assert np.array_equal(RNA,expected.full_source_RNA_sum) and np.array_equal(zero,expected.zero_full_RNA_cells)
    for actual,name in [(logRNA,'mean_log_full_source_RNA'),(loggenes,'mean_log_full_source_detected_genes')]:
        means=np.divide(actual,nc,out=np.full(nu,np.nan),where=nc>0)
        source=np.array([float(v) if str(v) not in ['','nan'] else np.nan for v in expected[name]])
        assert np.allclose(means,source,atol=1e-10,rtol=1e-10,equal_nan=True)
    assert (positive<=nc[:,None]).all() and (joint2<=joint).all()
    for j,cols in enumerate(ci):assert (joint[:,j]<=positive[:,cols].min(axis=1)).all()
    masked_sums=np.where(state['measured'],sums,np.nan);masked_positive=np.where(state['measured'],positive,np.nan)
    masked_joint=np.where(state['complex_measured'],joint,np.nan);masked_joint2=np.where(state['complex_measured'],joint2,np.nan)
    mean=np.divide(masked_sums,nc[:,None],out=np.full_like(sums,np.nan),where=nc[:,None]>0)
    fraction=np.divide(masked_positive,nc[:,None],out=np.full_like(sums,np.nan),where=nc[:,None]>0)
    joint_fraction=np.divide(masked_joint,nc[:,None],out=np.full((nu,nr),np.nan),where=nc[:,None]>0)
    np.savez_compressed(out/'private_all_frozen_pool_exact_panel_measurements.npz',
        genes=state['genes'],complexes=np.array(['_'.join(v) for v in state['complexes']]),unit_index=np.arange(nu),
        n_cells=nc,full_RNA_sum=RNA,zero_RNA_cells=zero,measured=state['measured'],complex_measured=state['complex_measured'],
        sum_logCP10K=masked_sums,mean_logCP10K=mean,positive_counts=masked_positive,positive_fraction=fraction,
        receptor_positive=masked_joint,receptor_positive2=masked_joint2,receptor_positive_fraction=joint_fraction)
    expected.to_csv(out/'private_all6760_exact_panel_unit_axis.tsv.gz',sep='\t',index=False,compression='gzip')
    files=[p for p in out.iterdir() if p.is_file()]
    assert sum(p.stat().st_size for p in files)<32*1024**2
    terminal=dict(status='COMPLETE_ALL_FROZEN_POOLS_EXACT_PANEL_AND_SAME_CELL_RECEPTORS_AUDIT_PENDING',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),helper_sha256=plan['helper_sha256'],
        frozen_inputs=state['hashes'],files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in files},
        inference_clinical_or_biological_replication=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'],completed_UTC=terminal['completed_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
