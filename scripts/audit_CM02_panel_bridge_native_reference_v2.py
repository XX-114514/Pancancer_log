"""Independent physical-column, selected metadata, raw counts and rank replay."""
import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
from prepare_CM02_panel_bridge_reference_scope_v1 import budget


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
    return h.hexdigest()


def manual(values,queries):
    n=len(queries);cap=1500;base=n*(n+1)/2;total=0.
    for value in queries:
        if value<=0:rank=cap
        else:
            rank=min(cap,1+np.count_nonzero(values>value)+(np.count_nonzero(values==value)-1)/2)
        total+=rank
    return 1-(total-base)/(n*cap-base)


def main(a):
    budget(a.root,16*1024**2)
    contract=json.loads((a.root/'reference_scope_v1/private_frozen_reference_contract.json').read_text())
    producer=a.root/'reference_native_inputs_v2';proof=json.loads((producer/'native_input_technical_public.json').read_text())
    for name,digest in proof['files'].items():assert sha(producer/name)==digest
    meta=pd.read_csv(producer/'private_selected_reference_metadata.tsv.gz',sep='\t',keep_default_na=False)
    archive=np.load(producer/'private_full_native_target_arrays.npz')
    target={key:archive[key] for key in archive.files}
    panel=np.load(producer/'private_panel_counts.npy',mmap_mode='r')
    original=meta.original_row.to_numpy(np.int64)
    z=zarr.open_group(str(a.reference),mode='r');obs=read_elem(z['obs']);var=read_elem(z['var'])
    selected=obs.iloc[original]
    for field in ['dataset','patient_id','library_id','major','cancer_type','timepoint','assay','annotation_status']:
        assert selected[field].astype(str).tolist()==meta[field].astype(str).tolist()
    assert selected.patient_identity_verified.astype(str).str.lower().isin(['true','1']).all()
    assert selected.patient_id_source_status.astype(str).eq('exact_original_library_official_metadata_verified').all()
    assert ~selected.annotation_exclude_from_primary.astype(str).str.lower().isin(['true','1']).any()
    measured=np.asarray(read_elem(z['varm/measured_in_source']),bool)
    source_order=list(map(str,read_elem(z['uns/source_order'])))
    lookup={}
    requested=sorted(set(contract['shared_panel_genes']+contract['core14']+contract['LR_excluded_genes']))
    literal=pd.DataFrame({'label':list(map(str,var.index))+var.gene_symbol.astype(str).tolist(),
                          'col':list(range(len(var)))*2})
    extra=var.source_gene_aliases.astype(str).str.split('|').explode()
    aliases=pd.DataFrame({'label':extra.to_numpy(),'col':var.index.get_indexer(extra.index)})
    aliases=aliases[var.feature_identity_status.astype(str).to_numpy()[aliases.col]=='stable_Ensembl_ID']
    pairs=pd.concat([literal,aliases],ignore_index=True).drop_duplicates()
    groups=pairs.groupby('label',sort=False).col.apply(list).to_dict()
    for source in contract['sources']:
        mapping={}
        for gene in requested:
            # Expanded table, grouped literal membership, then source mask.
            hits=[int(j) for j in groups.get(gene,[]) if measured[j,source_order.index(source)]
                  and str(var.feature_identity_status.iloc[j])!='ambiguous_source_feature']
            mapping[gene]=hits[0] if len(hits)==1 else -1
            assert mapping[gene]==contract['source_native_gene_columns'][source][gene]
        # Build dense physical -> requested position maps for independent cell scans.
        idx=np.full(len(var),-1,np.int32)
        for j,g in enumerate(contract['shared_panel_genes']):idx[mapping[g]]=j
        lookup[source]=(mapping,idx)
    raw=z['layers/counts'];ip=raw['indptr'][:]
    all_values=0;rank_values=0;maxerr=0.;done=0
    for start in range(0,len(ip)-1,5000):
        stop=min(len(ip)-1,start+5000);lo=int(np.searchsorted(original,start));hi=int(np.searchsorted(original,stop))
        if lo==hi:continue
        offset=int(ip[start]);end=int(ip[stop]);values=raw['data'][offset:end];indices=raw['indices'][offset:end]
        for k in range(lo,hi):
            row=int(original[k]);v=values[int(ip[row])-offset:int(ip[row+1])-offset];cols=indices[int(ip[row])-offset:int(ip[row+1])-offset]
            positive=v>0;v=v[positive];cols=cols[positive]
            assert len(np.unique(cols))==len(cols), 'Original noncanonical CSR requires separately frozen canonicalization audit'
            mapping,positions=lookup[meta.dataset.iloc[k]]
            where=positions[cols];use=where>=0
            expect=np.zeros(panel.shape[1],np.int64);expect[where[use]]=v[use]
            assert np.array_equal(expect,panel[k]);all_values+=len(expect)
            assert float(v.sum())==target['full_native_total'][k]
            assert len(cols)==target['full_native_detected_genes'][k]
            core=np.array([v[cols==mapping[g]].sum() if mapping[g]>=0 else -1 for g in contract['core14']],np.int64)
            assert np.array_equal(core,target['core14_native_counts'][k])
            if meta.lineage.iloc[k]=='T' and np.all(core[:11]>=0):
                scores=np.array([manual(v,core[:5]),manual(v,core[5:11])])
                maxerr=max(maxerr,float(np.max(abs(scores-target['original_full_signature_UCell1500'][k]))));rank_values+=2
            else:assert np.isnan(target['original_full_signature_UCell1500'][k]).all()
            done+=1
        if done//30000!=(done-(hi-lo))//30000:
            print(json.dumps(dict(phase='independent_original_counts_targets',cells=done,total=len(meta))),flush=True)
    assert done==len(meta) and maxerr<1e-12
    out=a.root/'reference_native_audit_v2';assert not out.exists();out.mkdir(mode=0o700)
    result=dict(status='ALL_SELECTED_REFERENCE_RAW_COUNTS_PHYSICAL_AXES_METADATA_AND_NATIVE_UCELL_INDEPENDENTLY_ACCEPTED',
                finished_utc=datetime.now(timezone.utc).isoformat(),selected_cells=done,
                shared_panel_values_replayed=all_values,full_signature_scores_replayed=rank_values,
                maximum_UCell_difference=maxerr,producer_terminal_sha256=sha(producer/'native_input_technical_public.json'),
                auditor_sha256=sha(Path(__file__)),no_model_fit_or_spatial_effect=True)
    (out/'independent_native_audit_technical_public.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for field in ['root','reference']:p.add_argument('--'+field,type=Path,required=True)
    main(p.parse_args())
