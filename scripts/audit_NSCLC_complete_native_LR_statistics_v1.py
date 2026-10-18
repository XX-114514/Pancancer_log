"""Full raw replay of every prospective LR gene in every original author unit.

Only comparison proofs are persisted; no new expression matrix is written.
"""
import gzip
import json
import os
import time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from prepare_NSCLC_full_native_gene_units_v1 import TrackedReader,chunks
from immune_RNA_MCP_new_source_io_v1 import wait


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,2*1024**2)
    out=root/'NSCLC_complete_native_LR_statistics_audit_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    accepted=wait(root,out,'audit_NSCLC_communication_tensor_inputs_v2',root/'NSCLC_communication_tensor_inputs_audit_20261019_v2/private_terminal.json',{'ACCEPTED_ALL243_BLIND_PARTITIONS_ALL51_TYPE_GATES_COMPLETE_RESOURCE_NATIVE_LR_AXIS_AND20_SLOTS'})
    src=root/'NSCLC_communication_tensor_inputs_20261019_v1';native=root/'NSCLC_full_native_gene_units_20261016_v1'
    scope=json.loads((src/'private_terminal.json').read_text());assert sha(src/'private_complete_primary_LR_gene_native_axis.tsv')==scope['files']['private_complete_primary_LR_gene_native_axis.tsv']['sha256']
    axis=pd.read_csv(src/'private_complete_primary_LR_gene_native_axis.tsv',sep='\t');cols=axis.native_column.to_numpy(int);assert len(cols)==len(set(cols))==1840
    u=pd.read_csv(native/'private_all_native_units_and_quality.tsv',sep='\t');assert u.unit_index.tolist()==list(range(11701))
    unit=np.load(native/'private_cell_unit_index.npy');total=np.load(native/'private_full_RNA.npy');detected=np.load(native/'private_detected_genes.npy');assert len(unit)==len(total)==1254749
    np.testing.assert_array_equal(np.bincount(unit,minlength=len(u)),u.n_cells)
    download=json.loads((root/'count_matrix_download_v1/download_state.json').read_text());assert download['status']=='COMPLETE_DOWNLOAD_GZIP_INTEGRITY_AND_SHA256';matrix=Path(download['path']);assert matrix.stat().st_size==download['bytes']
    save=lambda name,x:(out/name).write_text(json.dumps(x,indent=2)+'\n')
    memory={s.split(':')[0]:s.split(':')[1].strip() for s in Path('/proc/meminfo').read_text().splitlines()};assert int(memory['MemAvailable'].split()[0])*1024>4*1024**3
    save('plan_before_all_native_LR_replay.json',dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,input_audit=accepted,scope='Every requested1840native LRgene in all11701original units; all2010550708source coordinates and full1254749RNA/detected rowtotals. No unit sampling. No new expression array persisted.',source_expected_sha256=download['sha256'],source_bytes=download['bytes'],true_zero='Full zero-initialized measured gene axes, not only positive records.',arithmetic='Exact int64 raw sums and integer detection; independent bincount per compressed chunk for log1p CP10K sums; divide by exact all-cell unit denominator.'))
    raw=np.zeros((len(u),len(cols)),np.int64);positive=np.zeros_like(raw);logs=np.zeros_like(raw,dtype=float);rowtot=np.zeros(len(unit),np.int64);rowdet=np.zeros(len(unit),np.int64)
    reverse=np.full(31831,-1,int);reverse[cols]=np.arange(len(cols));entries=0;last=-1;lastcheckpoint=-1;start=time.monotonic()
    with matrix.open('rb') as binary:
        tracked=TrackedReader(binary)
        with gzip.GzipFile(fileobj=tracked,mode='rb') as stream:
            assert stream.readline().startswith(b'%%MatrixMarket matrix coordinate real general');line=stream.readline()
            while line.startswith(b'%'):line=stream.readline()
            dims=tuple(map(int,line.split()));assert dims==(1254749,31831,2010550708)
            for a in chunks(stream):
                assert np.isfinite(a).all() and np.equal(a,np.floor(a)).all() and (a[:,2]>0).all()
                rows=a[:,0].astype(np.int64)-1;gc=a[:,1].astype(np.int64)-1;values=a[:,2].astype(np.int64)
                assert rows.min()>=0 and rows.max()<len(unit) and gc.min()>=0 and gc.max()<31831
                index=rows*31831+gc;assert index[0]>last and np.all(index[1:]>index[:-1]);last=int(index[-1])
                lo=int(rows[0]);hi=int(rows[-1])+1
                rowtot[lo:hi]+=np.bincount(rows-lo,weights=values,minlength=hi-lo).astype(np.int64);rowdet[lo:hi]+=np.bincount(rows-lo,minlength=hi-lo)
                chosen=reverse[gc]>=0;r=rows[chosen];g=reverse[gc[chosen]];v=values[chosen];where=unit[r]*len(cols)+g
                # Reduce only occupied flattened positions, keeping working arrays bounded.
                occupied,inverse=np.unique(where,return_inverse=True)
                raw.flat[occupied]+=np.bincount(inverse,weights=v).astype(np.int64)
                positive.flat[occupied]+=np.bincount(inverse)
                logs.flat[occupied]+=np.bincount(inverse,weights=np.log1p(v.astype(float)*10000/total[r]))
                entries+=len(a);checkpoint=entries//50000000
                if checkpoint!=lastcheckpoint:
                    lastcheckpoint=checkpoint;save('private_progress.json',dict(phase='ALL_NATIVE_LR_RAW_REPLAY',entries=entries,expected_entries=dims[2],elapsed_seconds=time.monotonic()-start,observed_UTC=datetime.now(timezone.utc).isoformat()))
    assert entries==dims[2] and tracked.bytes==download['bytes'] and tracked.digest.hexdigest()==download['sha256']
    np.testing.assert_array_equal(rowtot,total);np.testing.assert_array_equal(rowdet,detected)
    logs/=u.n_cells.to_numpy()[:,None];errors={}
    for name,current,atol in [('sum_raw_counts',raw,0),('positive_cell_counts',positive,0),('mean_logCP10K',logs,1e-8)]:
        old=np.load(native/(name+'.npy'),mmap_mode='r');error=0.
        for begin in range(0,len(cols),64):
            end=min(begin+64,len(cols));a=current[:,begin:end];b=old[:,cols[begin:end]];error=max(error,float(np.max(np.abs(a-b))))
            if atol:np.testing.assert_allclose(a,b,atol=atol,rtol=1e-10)
            else:np.testing.assert_array_equal(a,b)
        errors[name]=error
    receipt=dict(status='ACCEPTED_ALL1840_NATIVE_LR_GENES_ALL11701_UNITS_FULL_RAW_COUNT_DETECTION_LOGCP10K_REPLAY',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),input_audit_terminal_sha256=sha(root/'NSCLC_communication_tensor_inputs_audit_20261019_v2/private_terminal.json'),producer_terminal_sha256=sha(src/'private_terminal.json'),matrix_sha256=tracked.digest.hexdigest(),source_entries=entries,native_LR_genes=len(cols),author_units=len(u),all_gene_unit_positions=len(cols)*len(u),max_absolute_errors=errors,full_RNA_and_detection_exact=True,resource_guard=guard,elapsed_seconds=time.monotonic()-start,new_matrix_persisted=False)
    save('private_terminal.json',receipt);print(json.dumps({'status':receipt['status'],'terminal_sha256':sha(out/'private_terminal.json')}))


if __name__=='__main__':main()
