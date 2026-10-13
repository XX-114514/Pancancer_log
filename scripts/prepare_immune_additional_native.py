#!/usr/bin/env python3
"""Create a private native-count input for the existing per-library QC worker."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import anndata as ad
import numpy as np
from scipy import sparse


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--input', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--dataset', required=True)
    parser.add_argument('--library-column', required=True)
    args=parser.parse_args()
    if args.output.exists():
        raise ValueError('Existing native attempt preserved')
    a=ad.read_h5ad(args.input)
    raw=sparse.csr_matrix(a.layers['counts'] if 'counts' in a.layers else a.X)
    if not a.obs_names.is_unique or not a.var_names.is_unique:
        raise ValueError('Native IDs/features must be unique')
    if (not np.isfinite(raw.data).all() or (raw.data<0).any() or
        (raw.data!=np.rint(raw.data)).any() or (raw.data>np.iinfo(np.int32).max).any()):
        raise ValueError('All native count entries must be nonnegative finite bounded integers')
    library=a.obs[args.library_column].astype(str)
    if library.isin(['', 'nan', 'None']).any():
        raise ValueError('Missing original library key')
    a.obs['original_cell_id']=a.obs_names.astype(str)
    a.obs['original_library_id']=library.to_numpy()
    a.obs['dataset_uid']=args.dataset
    a.obs['sample_uid']=args.dataset+'::'+library.to_numpy()
    a.obs['patient_uid']=''
    a.obs['clinical_tier']='study_context_only_patient_outcome_pending'
    a.layers.clear()
    a.layers['counts']=raw.astype(np.int32)
    # X is intentionally empty here; the existing worker uses layers/counts.
    a.X=sparse.csr_matrix(a.shape,dtype=np.int32)
    a.uns.clear()
    a.uns['native_source_semantics']='published_native_integer_counts; initial_sample_QC_pending'
    args.output.parent.mkdir(parents=True,exist_ok=True,mode=0o700)
    a.write_zarr(args.output)
    report=dict(observed_utc=datetime.now(timezone.utc).isoformat(),dataset=args.dataset,
                n_native_cells=a.n_obs,n_native_features=a.n_vars,n_libraries=library.nunique(),
                all_native_count_entries_validated=True,source_final_AC=False,
                counts_layer='counts',X_semantics='empty_input_placeholder_not_expression',
                prior_author_QC_not_assumed=True)
    with args.output.with_suffix('.native_input.json').open('x') as stream:json.dump(report,stream,indent=2)
    print(json.dumps(report),flush=True)


if __name__=='__main__':main()
