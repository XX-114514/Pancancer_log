#!/usr/bin/env python3
"""Independently read every final row and compare to audited source shards."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re

import anndata as ad
from anndata._io.specs import read_elem
import numpy as np
import pandas as pd
from scipy import sparse
import zarr
from public_gene_identity_contract import build_contract


GUARDED={'CXG_843c1459_e852_429c_9339_17b5760e7716',
         'CXG_67b6b9ac_b841_42ce_b2d4_ed2ec12327bd',
         'CXG_34f5307e_7b4d_4a48_b68f_2ba844c6414b'}


def block(group,lo,hi):
    pointer=group['indptr'][lo:hi+1];begin,end=int(pointer[0]),int(pointer[-1])
    return sparse.csr_matrix((group['data'][begin:end],group['indices'][begin:end],pointer-begin),
                            shape=(hi-lo,int(group.attrs['shape'][1])))


def difference_ok(left,right,tolerance=0):
    diff=left-right;diff.eliminate_zeros()
    return bool(not diff.nnz or np.max(np.abs(diff.data))<=tolerance)


def project_expected(matrix,positions,width):
    # Independent COO assembly sums any native aliases of the same stable ID.
    triplets=matrix.tocoo()
    return sparse.coo_matrix((triplets.data,(triplets.row,positions[triplets.col])),
                              shape=(matrix.shape[0],width)).tocsr()


def audit(root):
    manifest=json.loads((root/'manifests/combined_private_index.json').read_text())
    terminal=json.loads((root/'manifests/merge_terminal.json').read_text())
    g=zarr.open_group(str(root/'outputs/public_preprocessed_union_v1.zarr'),mode='r')
    obs=read_elem(g['obs']);var=read_elem(g['var']);uns=read_elem(g['uns'])
    genes=var.index.astype(str);feature_pos={name:i for i,name in enumerate(genes)}
    source_order=list(uns['source_order']);measured=np.asarray(g['varm']['measured_in_source'])
    source_var={}
    for record in manifest['records']:
        if record['dataset'] not in source_var:
            source_var[record['dataset']]=read_elem(zarr.open_group(record['source_prepared_ref'],mode='r')['var'])
    reconstructed,expected_var,_,_=build_contract(source_var)
    global_checks=dict(global_ID_unique=obs.index.is_unique,global_features_unique=var.index.is_unique,
        exact_total_rows=len(obs)==terminal['n_output_cells'],exact_total_genes=len(var)==terminal['n_genes'],
        source_final_AC_stays_false=not bool(uns['source_final_AC']),
        no_batch_correction_claim=not bool(uns['batch_correction_performed']),
        no_new_annotation_claim=not bool(uns['biological_annotation_performed']),
        gene_measurement_mask_shape_exact=measured.shape==(len(var),len(source_order)))
    global_checks['stable_gene_identity_contract_exact']=bool(reconstructed==manifest['canonical_gene_lists'] and
        np.array_equal(var.index,expected_var.index) and
        all(np.array_equal(var[c].astype(str),expected_var[c].astype(str)) for c in expected_var.columns))
    for name in source_order:
        expected=genes.isin(manifest['canonical_gene_lists'][name])
        global_checks['source_gene_mask_'+name]=np.array_equal(expected,measured[:,source_order.index(name)])
    sources={};results=[];adopted_total=0;cells=0
    for record in manifest['records']:
        source=zarr.open_group(record['source_prepared_ref'],mode='r')
        native_obs=read_elem(source['obs']);native_var=read_elem(source['var'])
        keep=np.load(record['keep_mask_ref']);expected_obs=native_obs.iloc[np.flatnonzero(keep)]
        lo,hi=record['row_start'],record['row_end'];actual_obs=obs.iloc[lo:hi]
        mapping=np.array([feature_pos[str(x)] for x in manifest['canonical_gene_lists'][record['dataset']]],dtype=np.int32)
        checks=dict(ID_order_exact=np.array_equal(actual_obs.index.astype(str),expected_obs.index.astype(str)),
            row_count_exact=len(actual_obs)==int(keep.sum())==record['n_kept'],source_metadata_exact=True,
            raw_counts_exact=True,source_raw_count_totals_conserved=True,removal_nonnegative_finite=True,removal_not_above_original=True,
            adopted_counts_match_source=True,normalization_exact=True,normalization_total_exact=True,
            adopted_correction_flags_exact=True,unmeasured_features_never_nonzero=True,
            conservative_correction_choice=True)
        kinds=manifest['canonical_kinds'];missing=uns['string_missing_sentinel']
        for column in expected_obs.columns:
            if kinds[column]=='numeric':
                a=pd.to_numeric(expected_obs[column]).to_numpy(dtype=float)
                b=actual_obs[column].to_numpy(dtype=float)
                checks['source_metadata_exact'] &= bool(np.array_equal(a,b,equal_nan=True))
            elif kinds[column]=='boolean':
                left=expected_obs[column].astype('boolean');right=actual_obs[column].astype('boolean')
                checks['source_metadata_exact'] &= bool(np.array_equal(left.isna(),right.isna()) and
                    np.array_equal(left.to_numpy(dtype=bool,na_value=False),right.to_numpy(dtype=bool,na_value=False)))
            else:
                expected=expected_obs[column].astype(str).where(~expected_obs[column].isna(),missing)
                checks['source_metadata_exact'] &= bool(np.array_equal(expected,actual_obs[column].astype(str)))
        local_offset=0;summary=record['summary']
        model_adopt=bool(summary.get('ambient_validation_pass',False) and record['dataset'] not in GUARDED)
        if model_adopt:
            log=(Path(record['source_prepared_ref']).parent/'decontX_method.log').read_text()
            iterations=re.findall(r'Completed iteration:\s*(\d+)\s*\|\s*converge:\s*([0-9.eE+-]+)',log)
            checks['conservative_correction_choice'] &= bool(iterations and float(iterations[-1][1])<0.001)
        for start in range(0,len(native_obs),1000):
            end=min(start+1000,len(native_obs));mask=keep[start:end]
            if not mask.any():continue
            positions=np.flatnonzero(mask);n=len(positions);a=lo+local_offset;b=a+n
            raw=block(source['layers']['counts'],start,end)[positions].tocsr()
            estimate=block(source['layers']['counts_decontx'],start,end)[positions].tocsr()
            original=block(g['layers']['counts'],a,b);removed=block(g['layers']['ambient_removed_counts'],a,b)
            normalized=block(g['X'],a,b)
            expected_raw=project_expected(raw.astype(np.int64),mapping,len(var))
            checks['raw_counts_exact'] &= difference_ok(original,expected_raw)
            checks['source_raw_count_totals_conserved'] &= bool(np.array_equal(np.asarray(raw.sum(axis=1)),np.asarray(original.sum(axis=1))))
            expected_mask=np.zeros(len(var),dtype=bool);expected_mask[mapping]=True
            for matrix in [original,removed,normalized]:
                checks['unmeasured_features_never_nonzero'] &= bool(not matrix.nnz or expected_mask[matrix.indices].all())
            checks['removal_nonnegative_finite'] &= bool(np.isfinite(removed.data).all() and (not removed.nnz or removed.data.min()>=0))
            excess=removed-original
            checks['removal_not_above_original'] &= bool(not excess.nnz or excess.data.max()<=1e-7)
            slice_obs=native_obs.iloc[start:end].iloc[positions]
            losses=pd.to_numeric(slice_obs.get('ambient_UMI_removed_fraction',pd.Series(0,index=slice_obs.index))).to_numpy()
            flags=np.full(n,model_adopt,dtype=bool)&(losses<=0.5)
            expected=raw.astype(np.float64).tolil()
            if flags.any():expected[np.flatnonzero(flags)]=estimate[np.flatnonzero(flags)]
            expected=expected.tocsr()
            expected=project_expected(expected,mapping,len(var))
            adopted=original.astype(np.float64)-removed
            checks['adopted_counts_match_source'] &= difference_ok(adopted,expected,1e-4)
            totals=np.asarray(adopted.sum(axis=1)).ravel()
            checks['normalization_total_exact'] &= bool(np.allclose(totals,obs.analysis_total_counts.iloc[a:b],rtol=0,atol=1e-6))
            intended=(sparse.diags(np.divide(10000.,totals,out=np.zeros(n),where=totals>0))@adopted).tocsr()
            intended.data=np.log1p(intended.data)
            checks['normalization_exact'] &= difference_ok(normalized,intended,2e-5)
            checks['adopted_correction_flags_exact'] &= bool(np.array_equal(flags,obs.analysis_correction_applied.iloc[a:b].to_numpy()))
            adopted_total+=int(flags.sum());local_offset+=n
        cells+=len(actual_obs)
        result=dict(dataset=record['dataset'],library_token=record['library_token'],n_cells=len(actual_obs),
                    all_checks_pass=bool(all(checks.values())),checks={k:bool(v) for k,v in checks.items()})
        results.append(result);sources[record['dataset']]=sources.get(record['dataset'],0)+len(actual_obs)
        print(json.dumps(dict(phase='audited_combined_shard',dataset=record['dataset'],library=record['library_token'],n_cells_audited=cells,all_checks_pass=result['all_checks_pass'])),flush=True)
    global_checks['all_selected_rows_scanned']=cells==len(obs)
    global_checks['adopted_count_total_matches']=adopted_total==terminal['n_cells_with_adopted_DecontX']
    for key in ['X','counts','ambient_removed_counts']:
        matrix=g['X'] if key=='X' else g['layers'][key]
        global_checks['CSR_shape_and_terminal_pointer_'+key]=bool(list(matrix.attrs['shape'])==[len(obs),len(var)] and
            int(matrix['indptr'][-1])==len(matrix['data'])==len(matrix['indices']) and len(matrix['indptr'])==len(obs)+1)
        # Use the AnnData reader's sparse accessor, independently of our blocks.
        global_checks['AnnData_sparse_accessor_'+key]=ad.io.sparse_dataset(matrix).shape==(len(obs),len(var))
    return dict(schema='combined_public_zarr_independent_audit_v1',observed_utc=datetime.now(timezone.utc).isoformat(),
        n_sources=len(sources),n_libraries=len(results),n_cells=len(obs),n_genes=len(var),n_adopted_correction_cells=adopted_total,
        all_checks_pass=bool(all(global_checks.values()) and all(r['all_checks_pass'] for r in results)),
        global_checks={k:bool(v) for k,v in global_checks.items()},source_cell_counts=sources,libraries=results,
        source_final_AC=False,final_accepted_new_cells=0,batch_correction_performed=False,
        merge_manifest_sha256=hashlib.sha256((root/'manifests/merge_terminal.json').read_bytes()).hexdigest(),
        audit_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--final-root',required=True);p.add_argument('--report',required=True)
    args=p.parse_args();report=audit(Path(args.final_root))
    with open(args.report,'x') as f:json.dump(report,f,indent=2,allow_nan=False)
    print(json.dumps({k:v for k,v in report.items() if k not in ['libraries','global_checks']}),flush=True)
    raise SystemExit(0 if report['all_checks_pass'] else 1)
