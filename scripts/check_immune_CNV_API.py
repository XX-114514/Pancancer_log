#!/usr/bin/env python3
"""Bounded real-native-count API check; never promotes cohort malignancy labels."""
import argparse,json,logging,sys
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


def main(args):
    a=ad.read_zarr(args.prepared)
    groups=a.obs.ambient_diagnostic_group.astype(str)
    candidate=np.flatnonzero(groups.eq('epithelial').to_numpy())[:400]
    reference=np.flatnonzero(groups.isin(['T_NK','B_plasma','endothelial','fibroblast']).to_numpy())[:400]
    if len(candidate)<50 or len(reference)<50:raise ValueError('Insufficient technical pilot groups')
    positions=np.concatenate([candidate,reference]);a=a[positions].copy()
    a.X=sparse.csr_matrix(a.layers['counts'],dtype=np.int32)
    a.obs['CNV_role']=np.repeat(['candidate','reference'],[len(candidate),len(reference)])
    a.obs['label']=a.obs.CNV_role.astype(str);a.obs['batch']='single_original_library'
    args.output.mkdir(mode=0o700,exist_ok=False)
    annotation=pd.read_csv(args.gene_annotation,sep=r'\s+',header=None,names=['gene','chr','start','end'])
    annotation=annotation.loc[annotation.chr.str.fullmatch(r'chr(?:[1-9]|1[0-9]|2[0-2])')].drop_duplicates('gene')
    coord=args.output/'coordinate.tsv';annotation.to_csv(coord,sep='\t',header=False,index=False)
    sys.path.insert(0,str(args.core_dir))
    from malignant_inference import run_malignant_inference
    sc.settings.n_jobs=4
    logger=logging.getLogger('CNV_API_check');logger.addHandler(logging.StreamHandler());logger.setLevel(logging.INFO)
    summary=run_malignant_inference(a,label_col='label',category_col='CNV_role',
        table_dir=args.output/'tables',figure_dir=args.output/'figures',logger=logger,
        gene_annotation_path=coord,cfg=dict(candidate_category_values=['candidate'],reference_category_values=['reference'],
        candidate_label_terms=[],reference_label_terms=[],exclude_label_terms=[],min_candidate_cells=50,min_reference_cells=50,
        window_size=250,copykat={'enabled':False},random_seed=20261013),batch_col='batch',run_label='technical_API_check')
    if summary['status']!='success':raise RuntimeError('Core CNV API pilot did not complete')
    report=dict(real_native_count_cells=a.n_obs,technical_API_check_only=True,
                diagnostic_QC_labels_not_final_identity=True,formal_cohort_CNVevaluated=False,
                core_status=summary['status'],copykat_status=summary['copykat_status'],
                n_finite_input_cell_scores=int(np.isfinite(pd.to_numeric(a.obs.cnv_score)).sum()),
                raw_counts=True,same_original_library=True)
    with args.report.open('x') as f:json.dump(report,f,indent=2,allow_nan=False)
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['prepared','output','core-dir','gene-annotation','report']:p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
