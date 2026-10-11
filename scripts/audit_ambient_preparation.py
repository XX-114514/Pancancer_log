#!/usr/bin/env python3
"""Independently scan all prepared counts, normalization, identities and gates."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
import zarr
from anndata._io.specs import read_elem


def csr(group, ptr, lo, hi, genes):
    start,end=int(ptr[lo]),int(ptr[hi])
    return sparse.csr_matrix((np.asarray(group['data'][start:end]),np.asarray(group['indices'][start:end]),
                             ptr[lo:hi+1]-start),shape=(hi-lo,genes))


def audit(root, collection):
    results=[]; global_ids=set(); ready=[]; source_genes={}
    for directory in sorted(root.iterdir()):
        if not directory.is_dir() or not (directory/'dataset_summary.json').exists():
            continue
        summary=json.loads((directory/'dataset_summary.json').read_text())
        name=summary['dataset']; native_contract=directory/'native_QC_input_contract.json'
        native=json.loads(native_contract.read_text()) if native_contract.exists() else None
        upstream_path=Path(native['native_source_ref']) if native else collection/'objects'/('qc_'+name+'.zarr')
        upstream=zarr.open_group(str(upstream_path),mode='r')
        upobs=read_elem(upstream['obs']);upvar=read_elem(upstream['var'])
        source_genes[name]=set(upvar.index.astype(str))
        original_ids=upobs['original_cell_id'].astype(str).to_numpy()
        upstream_positions={c:i for i,c in enumerate(original_ids)}
        if len(upstream_positions)!=len(upobs):
            raise ValueError('Original source IDs not unique: '+name)
        upgroup=upstream['layers']['counts'] if 'counts' in upstream['layers'] else upstream['X']
        upmat=ad.io.sparse_dataset(upgroup)
        observed=0;dataset_checks=[]
        for libsummary in sorted(summary['libraries'],key=lambda x:x['library_token']):
            if libsummary['status']!='success':
                dataset_checks.append(dict(library_token=libsummary['library_token'],all_checks_pass=False,
                                           failed=['worker_failed'],n_output=0));continue
            lib=directory/libsummary['library_token'];z=zarr.open_group(str(lib/'prepared.zarr'),mode='r')
            obs=read_elem(z['obs']);var=read_elem(z['var']);ids=obs.index.astype(str).to_numpy()
            indices=np.array([upstream_positions[c] for c in obs.original_cell_id.astype(str)],dtype=np.int64)
            checks=dict(gene_order_exact=var.index.equals(upvar.index),
                ID_unique=len(set(ids))==len(ids) and not bool(set(ids)&global_ids),
                upstream_identity=True,raw_counts_exact=True,corrected_nonnegative_finite=True,
                corrected_no_added_counts=True,normalization_exact=True,QC_metadata_preserved=True,
                release_gate_conservative=True,complete_filter_mapping=True,QC_metrics_from_counts=True)
            score_ledger=pd.read_csv(lib/'doublet_scores_private.tsv.gz',sep='\t',keep_default_na=False)
            lib_names=obs.original_library_id.astype(str).unique()
            if len(lib_names)!=1:
                raise ValueError('Prepared shard mixes original libraries')
            native_ledger=pd.read_csv(directory/'native_QC_ledger_private.tsv.gz',sep='\t',keep_default_na=False) if native else None
            source_expected=upobs[upobs.original_library_id.astype(str).eq(lib_names[0])].original_cell_id.astype(str)
            if native:
                qc_ids=set(native_ledger.loc[native_ledger.filter_reason.eq(''),'original_cell_id'].astype(str))
                source_expected=source_expected[source_expected.isin(qc_ids)]
            extra=set(score_ledger.loc[score_ledger.additional_doublet_removed.astype(str).str.lower().eq('true'),'original_cell_id'].astype(str))
            expected_ids=set(source_expected)-extra
            checks['complete_filter_mapping']=set(obs.original_cell_id.astype(str))==expected_ids
            global_ids.update(ids)
            if not native:
                checks['upstream_identity']=np.array_equal(ids,upobs.index.to_numpy()[indices])
            else:
                checks['upstream_identity']=np.array_equal(ids,
                    np.array([name+'::'+str(i)+'::'+str(original_ids[i]) for i in indices]))
            # Study context is unchanged, except the separately audited patient overlay.
            for field in ['original_library_id','clinical_tier','response','response_source','treatment','sample_uid']:
                if field in upobs:
                    checks['QC_metadata_preserved'] &= np.array_equal(obs[field].astype(str).to_numpy(),
                                                           upobs.iloc[indices][field].astype(str).to_numpy())
            if name=='GSE286827':
                overlay=pd.read_csv(collection/'artifacts/QC_patient_mapping_repair_20261011_v2/GSE286827_QC_cell_metadata_overlay_v2.tsv.gz',sep='\t',keep_default_na=False).set_index('stable_cell_id')
                expected=upobs.iloc[indices].copy()
                for f in ['patient_uid','timepoint','patient_identity_status']:
                    expected[f]=expected[f].astype(str)
                    hit=expected.index.intersection(overlay.index);expected.loc[hit,f]=overlay.loc[hit,f]
                    checks['QC_metadata_preserved'] &= np.array_equal(obs[f].astype(str),expected[f].astype(str))
            elif 'patient_uid' in upobs:
                checks['QC_metadata_preserved'] &= np.array_equal(obs.patient_uid.astype(str),upobs.iloc[indices].patient_uid.astype(str))
            mats={key:z['layers'][key] for key in ['counts','counts_decontx']};mats['X']=z['X']
            pointers={key:np.asarray(v['indptr']) for key,v in mats.items()}
            selected=str(read_elem(z['uns']['preprocess_v1']['normalized_from']))
            for lo in range(0,len(obs),1000):
                hi=min(lo+1000,len(obs)); raw=csr(mats['counts'],pointers['counts'],lo,hi,len(var))
                reference=upmat[indices[lo:hi]].tocsr()
                difference=raw-reference;difference.eliminate_zeros()
                checks['raw_counts_exact'] &= difference.nnz==0
                raw_total=np.asarray(raw.sum(axis=1)).ravel();raw_genes=raw.getnnz(axis=1)
                checks['QC_metrics_from_counts'] &= bool(np.array_equal(raw_total,obs.total_counts.to_numpy()[lo:hi]) and
                                                          np.array_equal(raw_genes,obs.n_genes.to_numpy()[lo:hi]))
                corrected=csr(mats['counts_decontx'],pointers['counts_decontx'],lo,hi,len(var))
                checks['corrected_nonnegative_finite'] &= bool(np.isfinite(corrected.data).all() and np.all(corrected.data>=0))
                difference=corrected-raw
                checks['corrected_no_added_counts'] &= bool(not difference.nnz or np.max(difference.data)<=1e-4)
                actual=csr(mats['X'],pointers['X'],lo,hi,len(var))
                intended=(corrected if selected=='counts_decontx' else raw).astype(np.float32).copy()
                totals=np.asarray(intended.sum(axis=1)).ravel()
                intended=sparse.diags(np.divide(1e4,totals,out=np.zeros(len(totals)),where=totals>0))@intended
                intended=intended.tocsr();intended.data=np.log1p(intended.data)
                difference=actual-intended
                checks['normalization_exact'] &= bool(not difference.nnz or np.max(np.abs(difference.data))<2e-5)
            eligible=bool(libsummary['preintegration_release_eligible'])
            residual=obs.get('residual_doublet_flag',pd.Series(False,index=obs.index)).astype(bool)
            if eligible:
                method_status=libsummary.get('independent_doublet_status')
                if method_status is None and libsummary.get('n_pre_Scrublet',0)>=500:
                    method_status='computed' if np.isfinite(obs.doubletdetection_score).all() else 'unverified'
                checks['release_gate_conservative']=bool(libsummary.get('ambient_validation_pass',False) and
                    method_status=='computed' and not residual.any() and
                    not libsummary.get('doublet_overcap_review',False))
                log=(lib/'decontX_method.log').read_text()
                convergence=re.findall(r'Completed iteration:\s*(\d+)\s*\|\s*converge:\s*([0-9.eE+-]+)',log)
                checks['release_gate_conservative'] &= bool(convergence and float(convergence[-1][1])<0.001)
                if 'preintegration_cell_eligible' in obs:
                    correction_loss=obs.ambient_UMI_removed_fraction.to_numpy(dtype=float)
                    checks['release_gate_conservative'] &= np.array_equal(obs.preintegration_cell_eligible.astype(bool),correction_loss<=0.5)
            observed+=len(obs)
            all_pass=all(checks.values())
            dataset_checks.append(dict(library_token=libsummary['library_token'],all_checks_pass=bool(all_pass),
                checks={k:bool(v) for k,v in checks.items()},n_output=len(obs)))
            if eligible and all_pass:
                n_ready=int(obs.preintegration_cell_eligible.sum()) if 'preintegration_cell_eligible' in obs else len(obs)
                ready.append(dict(dataset=name,library_token=libsummary['library_token'],n_cells=n_ready,
                    counts_layer='counts',analysis_layer='counts_decontx',native_gene_coverage=len(var),
                    source_final_AC=False,integration_allowed=False,
                    readiness='technical_preparation_pass_source_alias_gate_pending'))
        results.append(dict(dataset=name,n_libraries=len(dataset_checks),n_output=observed,
            summary_count_matches=observed==summary['n_output'],
            all_checks_pass=all(c['all_checks_pass'] for c in dataset_checks),libraries=dataset_checks))
    intersections=[]
    tumor_names=[g for g in source_genes if not g.startswith('CXG_5b9b89d2')]
    if tumor_names:
        common=set.intersection(*(source_genes[g] for g in tumor_names))
        intersections.append(dict(branch='tumor_and_comparator_sources_candidate',n_sources=len(tumor_names),
                                  common_measured_genes=len(common),gene_missingness_must_be_preserved=True))
    return dict(schema='independent_ambient_preparation_audit_v1',observed_utc=datetime.now(timezone.utc).isoformat(),
        n_datasets=len(results),n_libraries=sum(r['n_libraries'] for r in results),
        n_output=sum(r['n_output'] for r in results),n_ready_technical_cells=sum(r['n_cells'] for r in ready),
        all_checks_pass=bool(results and all(r['all_checks_pass'] and r['summary_count_matches'] for r in results)),
        final_accepted_new_cells=0,cross_study_integrated=False,datasets=results,
        integration_candidate_index=ready,gene_join_diagnostics=intersections)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output-root',required=True);p.add_argument('--collection',required=True)
    p.add_argument('--report',required=True);args=p.parse_args()
    report=audit(Path(args.output_root),Path(args.collection))
    with open(args.report,'x') as f:json.dump(report,f,indent=2,allow_nan=False)
    print(json.dumps({k:v for k,v in report.items() if k not in ['datasets','integration_candidate_index']}),flush=True)
    raise SystemExit(0 if report['all_checks_pass'] else 1)
