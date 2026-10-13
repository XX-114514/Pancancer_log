#!/usr/bin/env python3
"""Read-only representative UMAP QA; never displays patient or cell identifiers."""
import argparse,json
from pathlib import Path
import importlib.util
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem


def plot(args):
    g=zarr.open_group(str(args.object),mode='r');obs=read_elem(g['obs'])
    rng=np.random.default_rng(20261013)
    per=max(1,args.max_cells//obs.dataset.nunique())
    indices=np.sort(np.concatenate([rng.choice(rows,min(per,len(rows)),replace=False)
               for rows in obs.groupby('dataset',observed=True).indices.values()]))
    points=np.asarray(g['obsm/X_umap'])[indices];frame=obs.iloc[indices].copy()
    if not np.isfinite(points).all():raise ValueError('Invalid UMAP coordinates')
    figure,axes=plt.subplots(2,2,figsize=(15,11))
    columns=['dataset','cancer_type','major' if 'major' in frame else 'major_marker_proposal','capture_scope']
    titles=['Original study','Cancer context','Native multigene major proposal','Capture scope']
    for ax,column,title in zip(axes.ravel(),columns,titles):
        values=frame[column].astype(str).replace({'CXG_5b9b89d2_2e3e_4ad1_8965_666f64540fea':'CAR_T_product'})
        for j,name in enumerate(sorted(values.unique())):
            mask=values.eq(name).to_numpy()
            ax.scatter(points[mask,0],points[mask,1],s=.5,alpha=.5,color=plt.cm.tab20(j%20),rasterized=True,label=name)
        ax.set_title(title);ax.set_xlabel('UMAP1');ax.set_ylabel('UMAP2')
        ax.legend(fontsize=6,markerscale=5,loc='best',framealpha=.8)
    figure.suptitle(f'{args.stage}: {len(indices):,} representative cells / {len(obs):,} archived cells; proposals are not gold labels')
    figure.tight_layout();args.output.mkdir(mode=0o700,exist_ok=False)
    figure.savefig(args.output/'source_major_capture_UMAP.png',dpi=180);plt.close(figure)
    counts=pd.crosstab(obs.dataset.astype(str),obs['major' if 'major' in obs else 'major_marker_proposal'].astype(str))
    counts.to_csv(args.output/'source_major_aggregate.tsv',sep='\t')
    report=dict(stage=args.stage,n_archived_cells=len(obs),n_plot_cells=len(indices),
                sampling='seed20261013_balanced_source_display_only',all_cells_processed_not_claimed_by_plot=True,
                coordinates_finite=True,patient_or_cell_identifiers_displayed=False,major_labels='multigene_proposal_review_pending')
    with (args.output/'figure_scope_public.json').open('x') as f:json.dump(report,f,indent=2)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--object',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--max-cells',type=int,default=100000);p.add_argument('--stage',required=True)
    plot(p.parse_args())
