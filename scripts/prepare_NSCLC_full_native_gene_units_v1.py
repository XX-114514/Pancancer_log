"""Full native gene summaries by original case/author state; no clinical effects.

Read the existing coordinate matrix once, retaining all native genes and units.
No new download, cell selection, annotation changes or expression imputation.
"""
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(16777216),b''):h.update(block)
    return h.hexdigest()


class TrackedReader:
    def __init__(self,stream):self.stream=stream;self.digest=hashlib.sha256();self.bytes=0
    def read(self,n=-1):
        block=self.stream.read(n);self.digest.update(block);self.bytes+=len(block);return block
    def tell(self):return self.bytes


def chunks(stream,size=16777216):
    tail=b''
    while True:
        block=stream.read(size)
        if not block:
            if tail.strip():yield np.fromstring(tail.decode('ascii'),sep=' ').reshape(-1,3)
            return
        block=tail+block;cut=block.rfind(b'\n')
        if cut<0:tail=block;continue
        tail=block[cut+1:]
        yield np.fromstring(block[:cut+1].decode('ascii'),sep=' ').reshape(-1,3)


def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    old=R/'native_unit_expression_v1';out=R/'NSCLC_full_native_gene_units_20261016_v1'
    assert not out.exists();out.mkdir(mode=0o700)
    save=lambda name,value:(out/name).write_text(json.dumps(value,indent=2)+'\n')
    state_path=R/'count_matrix_download_v1/download_state.json'
    download=json.loads(state_path.read_text());assert download['status']=='COMPLETE_DOWNLOAD_GZIP_INTEGRITY_AND_SHA256'
    matrix=Path(download['path']);assert matrix.stat().st_size==download['bytes']
    assert json.loads((old/'independent_expression_audit.json').read_text())['status'].startswith('PASS_')
    paths=dict(metadata=R/'GSE243013_NSCLC_immune_scRNA_metadata.csv.gz',barcodes=R/'GSE243013_barcodes.csv.gz',
        genes=R/'GSE243013_genes.csv.gz',units=old/'unit_registry.tsv',prior_audit=old/'independent_expression_audit.json',download=state_path)
    registry=pd.read_csv(paths['units'],sep='\t',keep_default_na=False).sort_values('unit_index')
    assert registry.unit_index.tolist()==list(range(len(registry)))
    assert not registry.duplicated(['sampleID','sub_cell_type']).any()
    genes=pd.read_csv(paths['genes'],keep_default_na=False).geneSymbol.str.strip().str.upper().to_numpy(str)
    assert len(set(genes))==len(genes)
    bytes_planned=len(registry)*len(genes)*28
    memory={s.split(':')[0]:s.split(':')[1].strip() for s in Path('/proc/meminfo').read_text().splitlines()}
    available=int(memory['MemAvailable'].split()[0])*1024;free=shutil.disk_usage(out).free
    quota=subprocess.run(['quota','-g','-s','USER002'],capture_output=True,text=True,timeout=20)
    assert quota.returncode==0
    q=re.search(r'^\s*([0-9.]+)([KMGT])\s+\S+\s+([0-9.]+)([KMGT])',quota.stdout,re.M)
    assert q,'Group quota not resolved'
    scale={'K':1024,'M':1024**2,'G':1024**3,'T':1024**4}
    group_remaining=float(q[3])*scale[q[4]]-float(q[1])*scale[q[2]]
    assert available>bytes_planned*1.8+3*1024**3
    assert free>500*1024**3+bytes_planned*1.5
    assert group_remaining>bytes_planned*1.5+20*1024**3
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),
        scope='All original native cells, case IDs and author states; all native physical genes; no outcome selection or old selected200 restriction.',
        no_new_download=True,original_author_labels_changed=False,clinical_effects_executed=False,
        source_dimensions=download['dimensions'],native_units=len(registry),native_features=len(genes),
        input_sha256={k:sha(p) for k,p in paths.items()},matrix_expected_sha256=download['sha256'],matrix_bytes=download['bytes'],
        matrix_gate='SHA256 tracked while streaming compressed bytes, gzip CRC at EOF, strict sorted unique coordinates, positive integer values and full native per-cell RNA total parity before accepting outputs.',
        statistics='Integer raw sum and positive cell count; mean per-cell CP10K and natural log1pCP10K, all native genes and units including small groups.',
        denominator='Original full native RNA per cell, verified against every all-gene row sum at EOF.',
        missingness='All physical native genes are measured. A true zero stays zero. Genes outside this axis are not added or filled.',
        quality='Native integer full RNA and detected-gene counts; mean(log1p(RNA)) and mean(log1p(genes)) per original unit.',
        all_old_overlap='All previously provided genes in every original unit, all three mean/detection outputs; no favorable subset.',
        output_layout='Separate .npy arrays for mmap access, one version only; no duplicate full-cell matrix or intermediate array copies persisted.',
        planned_array_bytes=bytes_planned,mem_available_bytes=available,filesystem_free_bytes=free,
        rounded_group_quota_remaining_bytes=group_remaining,filesystem_reserve_bytes=500*1024**3,
        independent_raw_audit='Separately recompute all native genes in a prespecified outcome-free unit sample, plus all old-cache overlaps and global conservation; production alone is not acceptance.')
    save('plan_before_full_native_extraction.json',plan)
    md=pd.read_csv(paths['metadata'],usecols=['cellID','sampleID','sub_cell_type','major_cell_type','total_counts','n_genes','n_genes_by_counts'],keep_default_na=False)
    barcodes=pd.read_csv(paths['barcodes'],keep_default_na=False).barcode.to_numpy(str)
    assert np.array_equal(barcodes,md.cellID.to_numpy(str))
    assert pd.MultiIndex.from_frame(md[['sampleID','cellID']]).is_unique
    total=md.total_counts.to_numpy(float);assert np.isfinite(total).all() and (total>0).all() and np.array_equal(total,np.floor(total))
    keys=pd.MultiIndex.from_frame(registry[['sampleID','sub_cell_type']]);observed=pd.MultiIndex.from_frame(md[['sampleID','sub_cell_type']])
    unit=keys.get_indexer(observed);assert (unit>=0).all()
    n=np.bincount(unit,minlength=len(registry));assert np.array_equal(n,registry.n_cells.to_numpy(int))
    shape=(len(registry),len(genes));raw_sum=np.zeros(shape,np.int64);positive=np.zeros(shape,np.int32)
    cp=np.zeros(shape,np.float64);logcp=np.zeros(shape,np.float64)
    cell_total=np.zeros(len(md),np.int64);cell_detected=np.zeros(len(md),np.int32)
    arrays=[raw_sum,positive,cp,logcp];flat=[a.reshape(-1) for a in arrays]
    started=time.monotonic();entries=0;last=-1;previous_checkpoint=-1
    with matrix.open('rb') as binary:
        tracked=TrackedReader(binary)
        with gzip.GzipFile(fileobj=tracked,mode='rb') as stream:
            assert stream.readline().startswith(b'%%MatrixMarket matrix coordinate real general')
            line=stream.readline()
            while line.startswith(b'%'):line=stream.readline()
            dims=tuple(map(int,line.split()));assert dims==(len(md),len(genes),2010550708)
            save('private_progress.json',dict(phase='STREAMING_ALL_NATIVE_GENES',entries=0,expected_entries=dims[2],observed_UTC=datetime.now(timezone.utc).isoformat()))
            for a in chunks(stream):
                assert np.isfinite(a).all() and np.array_equal(a,np.floor(a)) and (a[:,2]>0).all()
                row=a[:,0].astype(np.int64)-1;col=a[:,1].astype(np.int64)-1;value=a[:,2].astype(np.int64)
                assert row.min()>=0 and row.max()<len(md) and col.min()>=0 and col.max()<len(genes)
                coordinate=row*len(genes)+col
                assert coordinate[0]>last and (np.diff(coordinate)>0).all(),'Duplicate or unsorted source coordinates'
                last=int(coordinate[-1]);lo=int(row[0]);hi=int(row[-1])+1
                cell_total[lo:hi]+=np.bincount(row-lo,weights=value,minlength=hi-lo).astype(np.int64)
                cell_detected[lo:hi]+=np.bincount(row-lo,minlength=hi-lo).astype(np.int32)
                index=unit[row]*len(genes)+col
                normalized=value.astype(float)*10000/total[row]
                np.add.at(flat[0],index,value);np.add.at(flat[1],index,1)
                np.add.at(flat[2],index,normalized);np.add.at(flat[3],index,np.log1p(normalized))
                entries+=len(a);checkpoint=entries//25000000
                if checkpoint!=previous_checkpoint:
                    previous_checkpoint=checkpoint
                    save('private_progress.json',dict(phase='STREAMING_ALL_NATIVE_GENES',entries=entries,expected_entries=dims[2],
                        elapsed_seconds=time.monotonic()-started,observed_UTC=datetime.now(timezone.utc).isoformat()))
    assert entries==dims[2] and tracked.bytes==download['bytes'] and tracked.digest.hexdigest()==download['sha256']
    assert np.array_equal(cell_total,total.astype(np.int64)) and (cell_detected>0).all()
    assert np.array_equal(positive.sum(axis=1,dtype=np.int64),np.bincount(unit,weights=cell_detected,minlength=len(n)).astype(np.int64))
    assert np.array_equal(raw_sum.sum(axis=1),np.bincount(unit,weights=cell_total,minlength=len(n)).astype(np.int64))
    cp/=n[:,None];logcp/=n[:,None]
    assert np.allclose(cp.sum(axis=1),10000,atol=1e-6,rtol=0)
    assert np.isfinite(cp).all() and np.isfinite(logcp).all() and (cp>=0).all() and (logcp>=0).all()
    save('private_progress.json',dict(phase='ALL_OLD_CACHE_INTERSECTION_PARITY',entries=entries,observed_UTC=datetime.now(timezone.utc).isoformat()))
    gene_index={gene:i for i,gene in enumerate(genes)}
    with np.load(old/'native_unit_expression.npz',allow_pickle=False) as previous:
        flags=previous['provided'].astype(bool);old_genes=previous['genes'].astype(str)
        assert np.array_equal(previous['unit_index'],registry.unit_index) and np.array_equal(previous['n_cells'],n)
        old_columns=np.flatnonzero(flags);new_columns=np.array([gene_index[g] for g in old_genes[flags]])
        max_errors={}
        for key,new_array in [('mean_CP10K',cp),('mean_logCP10K',logcp),('detected_fraction',positive)]:
            expected=previous[key];error=0.
            for start in range(0,len(old_columns),256):
                oc=old_columns[start:start+256];nc=new_columns[start:start+256]
                current=new_array[:,nc].astype(float)
                if key=='detected_fraction':current/=n[:,None]
                error=max(error,float(np.max(np.abs(current-expected[:,oc]))))
                assert np.allclose(current,expected[:,oc],atol=1e-8 if key!='detected_fraction' else 1e-12,rtol=1e-9 if key!='detected_fraction' else 0)
            max_errors[key]=error;del expected
    registry['major_cell_type']=pd.Series(md.major_cell_type.to_numpy(),index=unit).groupby(level=0).agg(lambda x: x.iloc[0]).reindex(range(len(registry))).to_numpy()
    assert pd.DataFrame(dict(unit=unit,major=md.major_cell_type)).groupby('unit').major.nunique().eq(1).all()
    registry['mean_logRNA']=np.bincount(unit,weights=np.log1p(cell_total),minlength=len(n))/n
    registry['mean_loggenes']=np.bincount(unit,weights=np.log1p(cell_detected),minlength=len(n))/n
    registry['native_total_RNA']=raw_sum.sum(axis=1)
    registry['native_detected_gene_sum']=positive.sum(axis=1,dtype=np.int64)
    registry.to_csv(out/'private_all_native_units_and_quality.tsv',sep='\t',index=False)
    pd.DataFrame(dict(native_column=np.arange(len(genes)),gene=genes,source_measured=True)).to_csv(out/'private_native_gene_axis.tsv',sep='\t',index=False)
    outputs={}
    for filename,array in [('sum_raw_counts.npy',raw_sum),('positive_cell_counts.npy',positive),('mean_CP10K.npy',cp),('mean_logCP10K.npy',logcp),('private_cell_unit_index.npy',unit.astype(np.int32)),('private_full_RNA.npy',cell_total),('private_detected_genes.npy',cell_detected)]:
        save('private_progress.json',dict(phase='WRITING_FULL_NATIVE_ARRAY',file=filename,observed_UTC=datetime.now(timezone.utc).isoformat()))
        np.save(out/filename,array,allow_pickle=False);outputs[filename]=dict(sha256=sha(out/filename),shape=list(array.shape),dtype=str(array.dtype),bytes=(out/filename).stat().st_size)
    for filename in ['private_all_native_units_and_quality.tsv','private_native_gene_axis.tsv']:
        outputs[filename]=dict(sha256=sha(out/filename),bytes=(out/filename).stat().st_size)
    report=dict(status='ALL_NATIVE_NSCLC_FULL_GENE_UNIT_SUMMARIES_COMPLETE_INDEPENDENT_RAW_AUDIT_PENDING',
        finished_UTC=datetime.now(timezone.utc).isoformat(),native_cells=len(md),native_units=len(registry),native_genes=len(genes),
        native_samples=int(md.sampleID.nunique()),author_states=int(md.sub_cell_type.nunique()),nonzero_coordinates=entries,
        matrix_sha256=tracked.digest.hexdigest(),all_native_RNA_totals_exact=True,all_old_overlap_genes=len(old_columns),all_old_overlap_max_errors=max_errors,
        metadata_gene_count_mismatches={field:int(np.count_nonzero(cell_detected!=pd.to_numeric(md[field],errors='coerce').to_numpy())) for field in ['n_genes','n_genes_by_counts']},
        outputs=outputs,new_download=False,cell_labels_changed=False,clinical_tests_executed=False)
    save('private_terminal.json',report);print(json.dumps(report))


if __name__=='__main__':main()
