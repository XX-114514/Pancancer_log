#!/usr/bin/env python3
"""Stream audited new-source shards into a compact AnnData Zarr.

Real cell rows and matrices remain private on the execution server. Original
counts are preserved; adopted correction is represented as a sparse removal
layer, avoiding another complete copy for cells whose counts were unchanged.
"""
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import sqlite3
import subprocess
import time

import anndata as ad
from anndata._io.specs import read_elem, write_elem
import numcodecs
import numpy as np
import pandas as pd
from scipy import sparse
import zarr

CAR='CXG_5b9b89d2_2e3e_4ad1_8965_666f64540fea'
NEURAL='CXG_843c1459_e852_429c_9339_17b5760e7716'
AUTHOR_CORRECTED={NEURAL,'CXG_67b6b9ac_b841_42ce_b2d4_ed2ec12327bd',
                  'CXG_34f5307e_7b4d_4a48_b68f_2ba844c6414b'}
RESERVE=500_000_000_000
COMPRESSOR=numcodecs.Blosc(cname='zstd',clevel=5,shuffle=numcodecs.Blosc.BITSHUFFLE)
MISSING='__SOURCE_VALUE_MISSING__'


def now():return datetime.now(timezone.utc).isoformat()


def save(path,value):
    with Path(path).open('x') as f:json.dump(value,f,indent=2,allow_nan=False)


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def csr_read(group,lo,hi):
    ptr=np.asarray(group['indptr'][lo:hi+1],dtype=np.int64)
    begin,end=int(ptr[0]),int(ptr[-1])
    return sparse.csr_matrix((np.asarray(group['data'][begin:end]),
        np.asarray(group['indices'][begin:end]),ptr-begin),shape=(hi-lo,int(group.attrs['shape'][1])))


def remap(matrix,positions,n_genes):
    result=sparse.csr_matrix((matrix.data.copy(),positions[matrix.indices],matrix.indptr.copy()),
                             shape=(matrix.shape[0],n_genes))
    result.sort_indices();return result


def adopted_counts(raw,estimate,obs,summary,dataset):
    adopt=bool(summary.get('ambient_validation_pass',False) and dataset not in AUTHOR_CORRECTED)
    loss=pd.to_numeric(obs.get('ambient_UMI_removed_fraction',pd.Series(0,index=obs.index))).to_numpy()
    selected=np.full(len(obs),adopt,dtype=bool)&(loss<=0.5)
    if selected.any():
        analysis=raw.astype(np.float64)+sparse.diags(selected.astype(np.float64))@(estimate.astype(np.float64)-raw)
    else:analysis=raw.astype(np.float64)
    removed=(raw.astype(np.float64)-analysis).tocsr();removed.eliminate_zeros()
    if removed.nnz and (not np.isfinite(removed.data).all() or removed.data.min() < -1e-4):
        raise ValueError('Invalid adopted removal values')
    if removed.nnz:removed.data=np.maximum(removed.data,0);removed.eliminate_zeros()
    analysis=(raw.astype(np.float64)-removed).tocsr()
    if analysis.nnz and analysis.data.min() < -1e-7:raise ValueError('Negative analysis counts')
    analysis.data=np.maximum(analysis.data,0);analysis.eliminate_zeros()
    total=np.asarray(analysis.sum(axis=1)).ravel()
    factors=np.divide(10000.,total,out=np.zeros(len(total)),where=total>0)
    normalized=(sparse.diags(factors)@analysis).tocsr()
    normalized.data=np.log1p(normalized.data).astype(np.float32)
    return removed,normalized,total,selected


def effective_capacity(path):
    physical=shutil.disk_usage(path).free
    result=subprocess.run(['quota','-g','-v','USER002'],capture_output=True,text=True,check=True)
    lines=result.stdout.splitlines()
    mount=subprocess.check_output(['df','--output=source',str(path)],text=True).splitlines()[-1].strip()
    for i,line in enumerate(lines):
        if mount not in line:continue
        suffix=line.split(mount,1)[1].strip()
        values=(suffix if suffix else lines[i+1].strip()).split()
        used=int(values[0].rstrip('*'));limit=int(values[2])
        quota=max(0,(limit-used)*1024) if limit else physical
        return min(physical,quota),dict(physical_free_bytes=physical,effective_free_bytes=min(physical,quota),
                                     quota_used_bytes=used*1024,quota_limit_bytes=limit*1024)
    raise RuntimeError('Cannot establish final-volume effective group quota')


class CSRAppender:
    def __init__(self,parent,key,dtype,n_genes):
        self.group=parent.create_group(key);self.group.attrs.update({'encoding-type':'csr_matrix','encoding-version':'0.1.0','shape':[0,n_genes]})
        self.data=self.group.create_dataset('data',shape=(0,),chunks=(1048576,),dtype=dtype,compressor=COMPRESSOR)
        self.indices=self.group.create_dataset('indices',shape=(0,),chunks=(1048576,),dtype='i4',compressor=COMPRESSOR)
        self.indptr=self.group.create_dataset('indptr',data=np.array([0],dtype='i8'),chunks=(65536,),compressor=COMPRESSOR)
        self.n_rows=0;self.n_genes=n_genes

    def append(self,matrix):
        matrix=matrix.tocsr()
        if matrix.shape[1]!=self.n_genes:raise ValueError('Feature count mismatch')
        if not np.isfinite(matrix.data).all():raise ValueError('Nonfinite matrix')
        if self.data.dtype.kind=='i':
            if (matrix.nnz and matrix.data.max()>np.iinfo(self.data.dtype).max) or not np.array_equal(matrix.data,np.rint(matrix.data)):
                raise ValueError('Integer count overflow/noninteger value')
        offset=len(self.data);self.data.append(matrix.data.astype(self.data.dtype,copy=False))
        self.indices.append(matrix.indices.astype(np.int32,copy=False))
        self.indptr.append(matrix.indptr[1:].astype(np.int64)+offset)
        self.n_rows+=matrix.shape[0];self.group.attrs['shape']=[self.n_rows,self.n_genes]


def discover(preprocess,audit_path):
    proof=json.loads(audit_path.read_text())
    if not proof['all_checks_pass'] or proof['n_datasets']!=9:raise ValueError('All nine sources require independent audit')
    records=[];gene_lists={};source_var={}
    for dataset in proof['datasets']:
        name=dataset['dataset'];directory=preprocess/'objects'/name
        source=json.loads((directory/'dataset_summary.json').read_text())
        if source['status']!='success' or source['pilot_only']:raise ValueError('Incomplete source')
        checked={r['library_token']:r for r in dataset['libraries']}
        for row in sorted(source['libraries'],key=lambda r:r['library_token']):
            token=row['library_token']
            if row['status']!='success' or not checked[token]['all_checks_pass']:raise ValueError('Incomplete shard')
            path=directory/token/'prepared.zarr';g=zarr.open_group(str(path),mode='r');var=read_elem(g['var'])
            genes=var.index.astype(str).tolist()
            if len(set(genes))!=len(genes):raise ValueError('Nonunique source features')
            if name in gene_lists and genes!=gene_lists[name]:raise ValueError('Source feature order drift')
            gene_lists[name]=genes
            source_var[name]=var[[c for c in ['source_feature_id','gene_id_original'] if c in var]].copy()
            records.append(dict(dataset=name,library_token=token,path=str(path),n_cells=row['n_output'],summary=row))
    return proof,records,gene_lists,source_var


def canonical_obs(frames):
    columns=sorted(set().union(*(set(f.columns) for f in frames)));kinds={};converted=[]
    for column in columns:
        present=[f[column] for f in frames if column in f]
        kinds[column]='boolean' if all(pd.api.types.is_bool_dtype(s.dtype) for s in present) else \
            'numeric' if all(pd.api.types.is_numeric_dtype(s.dtype) and not pd.api.types.is_bool_dtype(s.dtype) for s in present) else 'string'
        # Float canonicalization must not silently round large integer IDs.
        if kinds[column]=='numeric' and any(pd.api.types.is_integer_dtype(s.dtype) and
            ((s.dropna()>2**53).any() or (s.dropna() < -(2**53)).any()) for s in present):
            kinds[column]='string'
    for frame in frames:
        frame=frame.copy()
        for column in frame.columns:
            if kinds[column]=='numeric':frame[column]=pd.to_numeric(frame[column]).astype(np.float64)
            elif kinds[column]=='boolean':frame[column]=frame[column].astype('boolean')
            else:
                values=frame[column];converted_values=values.astype(str)
                if converted_values[~values.isna()].eq(MISSING).any():raise ValueError('Metadata missing sentinel collision')
                frame[column]=converted_values.where(~values.isna(),MISSING)
        converted.append(frame)
    obs=pd.concat(converted,axis=0,sort=False)
    if not obs.index.is_unique:raise ValueError('Global cell IDs are not unique')
    for column,kind in kinds.items():
        if kind=='string':obs[column]=pd.Categorical(obs[column].fillna(MISSING).astype(str))
        elif kind=='boolean':obs[column]=obs[column].astype('boolean')
    return obs,kinds


def fingerprints(records,gene_lists,collection,scratch):
    """Exact baseline and new-to-new matching; unresolvable aliases remain pending."""
    vocabulary=json.loads((collection/'manifests/fingerprint_gene_vocabulary.json').read_text())
    vocab={g:i for i,g in enumerate(vocabulary)}
    for g in sorted(set().union(*(set(v) for v in gene_lists.values()))-set(vocab)):vocab[g]=len(vocab)
    baseline=sqlite3.connect('file:'+str(collection/'manifests/count_fingerprints.sqlite')+'?mode=ro',uri=True)
    database=scratch/'combined_new_fingerprints.sqlite'
    if database.exists():raise ValueError('Existing fingerprint attempt preserved')
    new=sqlite3.connect(database);new.execute('CREATE TABLE seen(fp BLOB,barcode TEXT,cell_id TEXT,PRIMARY KEY(fp,barcode))')
    records_out=[];excluded=[]
    barcode_pattern=re.compile(r'([ACGT]{14,24})(?:-\d+)?')
    for record in records:
        g=zarr.open_group(record['path'],mode='r');obs=read_elem(g['obs'])
        positions=np.array([vocab[s] for s in gene_lists[record['dataset']]],dtype=np.int64)
        keep=np.ones(len(obs),dtype=bool)
        for lo in range(0,len(obs),1000):
            hi=min(lo+1000,len(obs));matrix=csr_read(g['layers']['counts'],lo,hi);matrix.eliminate_zeros()
            rows=[]
            for j in range(hi-lo):
                a,b=matrix.indptr[j:j+2];gi=positions[matrix.indices[a:b]];order=np.argsort(gi)
                pairs=np.column_stack((gi[order],matrix.data[a:b][order].astype('<i8'))).astype('<i8',copy=False)
                fingerprint=hashlib.blake2b(pairs.tobytes(),digest_size=24).digest()
                original=str(obs.original_cell_id.iloc[lo+j]);match=barcode_pattern.search(original)
                # Match the frozen baseline convention; do not invent a DNA barcode.
                barcode=match.group(1) if match else original
                rows.append((lo+j,fingerprint,barcode,str(obs.index[lo+j])))
            hashes=list({r[1] for r in rows if r[2]});found={}
            for at in range(0,len(hashes),400):
                hh=hashes[at:at+400]
                for cell,study,bc,fp in baseline.execute('SELECT cell_id,gse,barcode,fp FROM baseline WHERE fp IN ('+','.join('?' for _ in hh)+')',hh):
                    found[(fp,bc)]=(cell,study)
            for i,fp,barcode,cell in rows:
                if barcode is None:continue
                hit=found.get((fp,barcode));reason='exact_baseline_counts_and_barcode' if hit else None
                if not hit:
                    match=new.execute('SELECT cell_id FROM seen WHERE fp=? AND barcode=?',(fp,barcode)).fetchone()
                    if match:hit=(match[0],'new_collection');reason='exact_new_counts_and_barcode_review_quarantine'
                if hit:
                    keep[i]=False;excluded.append(dict(cell_id=cell,matched_cell_id=hit[0],matched_source=hit[1],reason=reason))
                else:new.execute('INSERT INTO seen VALUES(?,?,?)',(fp,barcode,cell))
            new.commit()
        record['keep_mask_ref']=str(scratch/(record['dataset']+'_'+record['library_token']+'_keep.npy'))
        np.save(record['keep_mask_ref'],keep);record['n_kept']=int(keep.sum());records_out.append(record)
        print(json.dumps(dict(utc=now(),phase='fingerprint_shard',dataset=record['dataset'],library=record['library_token'],n_kept=record['n_kept'],n_flagged=int((~keep).sum()))),flush=True)
    baseline.close();new.close()
    pd.DataFrame(excluded,columns=['cell_id','matched_cell_id','matched_source','reason']).to_csv(scratch/'duplicate_exclusions_private.tsv.gz',sep='\t',index=False)
    return records_out,len(excluded)


def calibration(records,gene_lists,global_genes):
    lookup={g:i for i,g in enumerate(global_genes)};totals={'counts':0,'X':0,'ambient_removed_counts':0};samples=[]
    for r in records:
        g=zarr.open_group(r['path'],mode='r');obs=read_elem(g['obs']);n=len(obs)
        mapping=np.array([lookup[x] for x in gene_lists[r['dataset']]],dtype=np.int32)
        lo=max(0,n//2-128);hi=min(n,lo+256)
        raw=csr_read(g['layers']['counts'],lo,hi);estimate=csr_read(g['layers']['counts_decontx'],lo,hi)
        removed,x,_,_=adopted_counts(raw,estimate,obs.iloc[lo:hi],r['summary'],r['dataset'])
        matrices={'counts':raw,'X':x,'ambient_removed_counts':removed};ratios={}
        for key,mat in matrices.items():
            mat=remap(mat,mapping,len(global_genes));dtype=np.int32 if key=='counts' else np.float32 if key=='X' else np.float64
            encoded=len(COMPRESSOR.encode(mat.data.astype(dtype)))+len(COMPRESSOR.encode(mat.indices.astype(np.int32)))
            ratio=encoded/max(1,raw.nnz);ratios[key]=ratio
            totals[key]+=int(ratio*int(g['layers']['counts']['data'].shape[0]))
        samples.append(dict(dataset=r['dataset'],library_token=r['library_token'],n_calibration_cells=hi-lo,bytes_per_raw_nonzero=ratios))
    matrix_bytes=sum(totals.values());planned=int(matrix_bytes*1.35)+12_000_000_000
    return dict(observed_utc=now(),matrix_estimate_bytes=matrix_bytes,matrix_estimate_by_layer=totals,
        headroom_multiplier=1.35,metadata_QC_and_chunk_overhead_allowance_bytes=12_000_000_000,
        expected_peak_bytes=planned,reserve_bytes=RESERVE,samples=samples)


def merge(preprocess,collection,final_root,scratch,code_sha,allow_subset=False):
    audit=preprocess/'manifests/independent_preparation_audit.json'
    proof,records,gene_lists,source_var=discover(preprocess,audit)
    if allow_subset:raise ValueError('Subset export cannot claim full production')
    sources=sorted(gene_lists);genes=sorted(set().union(*(set(v) for v in gene_lists.values())))
    budget=calibration(records,gene_lists,genes);free,capacity=effective_capacity(final_root)
    save(final_root/'manifests/final_storage_budget.json',dict(budget,capacity=capacity))
    if free<budget['expected_peak_bytes']+RESERVE:raise RuntimeError('Final compressed result plus 500GB reserve does not fit effective quota')
    records,n_excluded=fingerprints(records,gene_lists,collection,scratch)
    frames=[];source_coverage={};metadata_column_coverage={}
    for r in records:
        g=zarr.open_group(r['path'],mode='r');frame=read_elem(g['obs']);keep=np.load(r['keep_mask_ref'])
        frame=frame.iloc[np.flatnonzero(keep)].copy();frame['prepared_library_token']=r['library_token']
        frame['source_dataset']=r['dataset'];frame['source_final_AC']=False
        frame['source_alias_acceptance_status']='pending_primary_library_alias_review'
        frame['integration_branch']='manufactured_CAR_T_product' if r['dataset']==CAR else 'T_enriched_nuclei' if r['dataset']=='GSE292046' else 'tissue_or_comparator'
        frame['author_processing_or_neural_guard']=r['dataset'] in AUTHOR_CORRECTED
        frame['doublet_method_evaluable']=r['summary'].get('independent_doublet_status')=='computed'
        frame['downstream_release_eligible']=False
        metadata_column_coverage[r['dataset']]=sorted(set(metadata_column_coverage.get(r['dataset'],[]))|set(frame.columns))
        frames.append(frame)
    obs,kinds=canonical_obs(frames);del frames
    target=final_root/'outputs/public_preprocessed_union_v1.zarr'
    if target.exists():raise ValueError('Final attempt already exists')
    target.mkdir(mode=0o700)
    group=zarr.open_group(str(target),mode='w');group.attrs.update({'encoding-type':'anndata','encoding-version':'0.1.0'})
    layers=group.create_group('layers');layers.attrs.update({'encoding-type':'dict','encoding-version':'0.1.0'})
    arrays={'counts':CSRAppender(layers,'counts','i4',len(genes)),
            'ambient_removed_counts':CSRAppender(layers,'ambient_removed_counts','f8',len(genes)),
            'X':CSRAppender(group,'X','f4',len(genes))}
    lookup={g:i for i,g in enumerate(genes)};offset=0;private_index=[];sequence=0
    analysis_totals=np.zeros(len(obs),dtype=np.float64);adopted=np.zeros(len(obs),dtype=bool)
    for r in records:
        g=zarr.open_group(r['path'],mode='r');native_obs=read_elem(g['obs']);keep=np.load(r['keep_mask_ref'])
        mapping=np.array([lookup[x] for x in gene_lists[r['dataset']]],dtype=np.int32)
        start=offset
        for lo in range(0,len(native_obs),1000):
            hi=min(lo+1000,len(native_obs));mask=keep[lo:hi]
            if not mask.any():continue
            free,_=effective_capacity(final_root)
            if free<RESERVE+32_000_000:raise RuntimeError('Reached final-volume 500GB reserve during append')
            raw=csr_read(g['layers']['counts'],lo,hi)[mask].tocsr()
            estimate=csr_read(g['layers']['counts_decontx'],lo,hi)[mask].tocsr()
            removed,x,total,selected=adopted_counts(raw,estimate,native_obs.iloc[lo:hi].iloc[np.flatnonzero(mask)],r['summary'],r['dataset'])
            for key,mat in [('counts',raw),('ambient_removed_counts',removed),('X',x)]:arrays[key].append(remap(mat,mapping,len(genes)))
            analysis_totals[offset:offset+len(total)]=total;adopted[offset:offset+len(total)]=selected;offset+=len(total)
        private_index.append(dict(dataset=r['dataset'],library_token=r['library_token'],source_prepared_ref=r['path'],
            keep_mask_ref=r['keep_mask_ref'],row_start=start,row_end=offset,n_kept=offset-start,summary=r['summary']))
        save(final_root/'manifests'/('merge_progress_%04d.json'%sequence),dict(utc=now(),n_libraries=sequence+1,n_output_cells=offset,n_genes=len(genes)))
        sequence+=1
        print(json.dumps(dict(utc=now(),phase='merged_shard',dataset=r['dataset'],library=r['library_token'],n_output_cells=offset)),flush=True)
    if offset!=len(obs):raise ValueError('Cell order/row count mismatch')
    obs['analysis_total_counts']=analysis_totals;obs['analysis_correction_applied']=adopted
    obs['analysis_count_origin']=pd.Categorical(np.where(adopted,'validated_DecontX','published_source_counts'))
    var=pd.DataFrame(index=pd.Index(genes,name='gene_symbol'))
    measured=np.column_stack([var.index.isin(gene_lists[s]) for s in sources]).astype(bool)
    var['n_sources_measured']=measured.sum(axis=1)
    write_elem(group,'obs',obs);write_elem(group,'var',var)
    for name in ['obsm','obsp','varp']:write_elem(group,name,{})
    write_elem(group,'varm',{'measured_in_source':measured})
    uns=dict(collection_stage='all_retained_preprocessing_cells_source_acceptance_pending',
        source_order=sources,source_obs_column_coverage=metadata_column_coverage,metadata_canonical_kinds=kinds,
        string_missing_sentinel=MISSING,source_gene_coverage_semantics='unmeasured_zero_padding_is_not_biological_zero',
        counts_semantics='immutable_published_upstream_counts_after_audited_QC_doublet_selection',
        analysis_counts_recipe='layers[counts].astype(float64) - layers[ambient_removed_counts]',
        ambient_removed_counts_semantics='only_adopted_correction; diagnostic_rejected_estimates_not_applied',
        X_semantics='log1p(analysis_counts * 10000 / analysis_total_counts)',
        raw_counts_immutable=True,source_final_AC=False,batch_correction_performed=False,
        biological_annotation_performed=False,baseline_cells_unchanged=4676787,
        source_local_cluster_semantics='library_local_diagnostic_labels_not_global_cell_identity',
        exact_duplicate_rows_quarantined=n_excluded,code_git_sha=code_sha,
        preprocessing_audit_sha256=digest(audit),fullgene_join='outer_with_measured_in_source_mask')
    write_elem(group,'uns',uns)
    diagnostics=final_root/'diagnostics';diagnostics.mkdir(exist_ok=True)
    for r in records:
        folder=diagnostics/r['dataset']/r['library_token'];folder.mkdir(parents=True,exist_ok=True)
        for file in Path(r['path']).parent.iterdir():
            if file.is_file() and (file.suffix in ['.json','.png','.log','.tsv'] or file.name.endswith('.tsv.gz')):
                shutil.copyfile(file,folder/file.name)
    for source,frame in source_var.items():frame.to_csv(diagnostics/(source+'_source_features.tsv.gz'),sep='\t')
    shutil.copyfile(scratch/'duplicate_exclusions_private.tsv.gz',diagnostics/'duplicate_exclusions_private.tsv.gz')
    save(final_root/'manifests/combined_private_index.json',dict(gene_lists=gene_lists,records=private_index,canonical_kinds=kinds))
    result=dict(schema='combined_public_preprocessed_zarr_v1',completed_utc=now(),n_sources=len(sources),n_libraries=len(records),
        n_output_cells=len(obs),n_genes=len(genes),n_exact_duplicate_quarantined=n_excluded,
        n_cells_with_adopted_DecontX=int(adopted.sum()),source_final_AC=False,final_accepted_new_cells=0,
        batch_correction_performed=False,code_git_sha=code_sha,independent_final_audit='pending')
    save(final_root/'manifests/merge_terminal.json',result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--preprocess-root',required=True);p.add_argument('--collection',required=True)
    p.add_argument('--final-root',required=True);p.add_argument('--scratch-root',required=True);p.add_argument('--code-git-sha',required=True)
    a=p.parse_args();numcodecs.blosc.set_nthreads(2)
    print(json.dumps(merge(Path(a.preprocess_root),Path(a.collection),Path(a.final_root),Path(a.scratch_root),a.code_git_sha)),flush=True)
