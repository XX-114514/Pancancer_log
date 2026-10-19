"""Extract only frozen shared-panel rows; full-native tied ranks for two targets."""
import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from scipy.sparse import csr_matrix
from prepare_CM02_panel_bridge_reference_scope_v1 import budget


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
    return h.hexdigest()


def score(values, target, cap=1500):
    n=len(target);base=n*(n+1)/2
    order=np.sort(values[values>0])
    ranks=np.full(n,cap,dtype=np.float64)
    positive=target>0
    left=np.searchsorted(order,target[positive],side='left')
    right=np.searchsorted(order,target[positive],side='right')
    ranks[positive]=np.minimum(cap,1+len(order)-right+(right-left-1)/2)
    return 1-(ranks.sum()-base)/(n*cap-base)


def main(a):
    scope=a.root/'reference_scope_v1'
    contract=json.loads((scope/'private_frozen_reference_contract.json').read_text())
    assert sha(scope/'private_selected_reference_rows.tsv.gz')==contract['private_selected_rows_sha256']
    budget(a.root,contract['reserve']['estimated_additional_bytes'])
    rows=pd.read_csv(scope/'private_selected_reference_rows.tsv.gz',sep='\t',keep_default_na=False)
    original=rows.original_row.to_numpy(np.int64)
    assert np.all(np.diff(original)>0)
    z=zarr.open_group(str(a.reference),mode='r');raw=z['layers/counts'];indptr=raw['indptr'][:]
    nfeat=raw.attrs['shape'][1];genes=contract['shared_panel_genes'];core=contract['core14']
    out=a.root/'reference_native_inputs_v1';assert not out.exists();out.mkdir(mode=0o700)
    panel=np.lib.format.open_memmap(out/'private_panel_counts.npy',mode='w+',dtype=np.int32,shape=(len(rows),len(genes)))
    fulltot=np.zeros(len(rows),np.float64);fullgenes=np.zeros(len(rows),np.int32)
    corecounts=np.full((len(rows),len(core)),-1,np.int32)
    targets=np.full((len(rows),2),np.nan,np.float64)
    maps=contract['source_native_gene_columns']
    used=0
    for start in range(0,len(indptr)-1,5000):
        stop=min(len(indptr)-1,start+5000);lo=np.searchsorted(original,start);hi=np.searchsorted(original,stop)
        if lo==hi:continue
        low,high=int(indptr[start]),int(indptr[stop])
        block=csr_matrix((raw['data'][low:high],raw['indices'][low:high],indptr[start:stop+1]-low),shape=(stop-start,nfeat))
        block.sum_duplicates();block.eliminate_zeros()
        for source in rows.dataset.iloc[lo:hi].unique():
            chosen=np.flatnonzero(rows.dataset.iloc[lo:hi].to_numpy()==source)+lo
            piece=block[original[chosen]-start]
            assert np.isfinite(piece.data).all() and (piece.data>=0).all()
            assert np.equal(piece.data,np.rint(piece.data)).all()
            cols=[maps[source][g] for g in genes]
            dense=piece[:,cols].toarray();assert dense.max(initial=0)<2**31
            panel[chosen]=dense.astype(np.int32)
            fulltot[chosen]=np.asarray(piece.sum(1)).ravel()
            fullgenes[chosen]=np.diff(piece.indptr)
            for k,g in enumerate(core):
                col=maps[source][g]
                if col>=0:corecounts[chosen,k]=piece[:,col].toarray().ravel().astype(np.int32)
            complete=np.all(corecounts[chosen,:11]>=0,axis=1)
            for local,idx in enumerate(chosen):
                if rows.lineage.iloc[idx]=='T' and complete[local]:
                    values=piece.data[piece.indptr[local]:piece.indptr[local+1]]
                    targets[idx,0]=score(values,corecounts[idx,:5])
                    targets[idx,1]=score(values,corecounts[idx,5:11])
        used+=hi-lo
        if used//20000!=(used-(hi-lo))//20000:
            print(json.dumps(dict(phase='selected_native_panel_and_full_rank_targets',completed=used,total=len(rows))),flush=True)
    assert used==len(rows);panel.flush()
    np.savez_compressed(out/'private_full_native_target_arrays.npz',full_native_total=fulltot,
                        full_native_detected_genes=fullgenes,core14_native_counts=corecounts,
                        original_full_signature_UCell1500=targets)
    rows.to_csv(out/'private_selected_reference_metadata.tsv.gz',sep='\t',index=False,compression='gzip')
    terminal=dict(status='FROZEN_REFERENCE_SHARED_PANEL_NATIVE_COUNTS_AND_FULL_NATIVE_TARGETS_WRITTEN_AUDIT_PENDING',
                  finished_utc=datetime.now(timezone.utc).isoformat(),selected_cells=len(rows),
                  shared_panel_features=len(genes),source_count=len(contract['sources']),
                  source_contract_sha256=sha(scope/'private_frozen_reference_contract.json'),
                  producer_sha256=sha(Path(__file__)),native_array_bytes=(out/'private_panel_counts.npy').stat().st_size,
                  files={p.name:sha(p) for p in out.iterdir() if p.is_file()},
                  no_model_fit_no_spatial_effect=True)
    (out/'native_input_technical_public.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps({k:v for k,v in terminal.items() if k!='files'}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for field in ['root','reference']:p.add_argument('--'+field,type=Path,required=True)
    main(p.parse_args())
