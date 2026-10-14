#!/usr/bin/env python3
"""Descriptive capture, assay and annotation-uncertainty supplements."""
import argparse
import textwrap
from pathlib import Path
import numpy as np
import pandas as pd
from plot_immune_publication_suite import Suite,MAJOR_COLORS,STUDY_COLORS,study


def main(a):
    s=Suite(a)
    capture=s.obs.capture_scope.astype(str).to_numpy();names=sorted(set(capture));assert len(names)<=9
    titles={'CAR_T_infusion_product':'CAR-T infusion product','CD45_enriched_immune':'CD45 enrichment',
      'ex_vivo_antigen_stimulated_lymphocytes':'Ex vivo antigen stimulation',
      'lymphocyte_enriched':'Lymphocyte enrichment','whole_tumor':'Whole tumour',
      'whole_tumor_patient_mapping_pending':'Whole tumour; patient pending',
      'serial_tumor_biopsies':'Serial tumour biopsies',
      'mixed_tumor_nuclei_NFAT_enrichment':'Mixed NFAT tumour nuclei'}
    fig,axes=s.figure(170,cols=3,rows=3)
    fig.subplots_adjust(left=.09,right=.94,bottom=.10,hspace=.68,wspace=.50)
    for ax,name in zip(axes.flat,names):
        ids=np.flatnonzero(capture==name);display=s.sample(ids,s.major[ids],60000)
        s.scatter(ax,display,s.major,MAJOR_COLORS,titles.get(name,name))
        ax.text(.02,.02,f'n={len(ids):,}; display={len(display):,}',transform=ax.transAxes,fontsize=5)
    if len(names)<9:
        ax=axes.flat[len(names)];counts=pd.Series(capture).value_counts().reindex(names)
        ax.barh(np.arange(len(names)),counts/1000,color='#777777')
        ax.set_yticks(np.arange(len(names)));ax.set_yticklabels([textwrap.fill(titles.get(x,x),18) for x in names],fontsize=5)
        ax.set_xlabel('Archived cells (thousands)');ax.set_title('Capture scope: full counts')
        for extra in list(axes.flat)[len(names)+1:]:extra.set_axis_off()
    pd.Series(capture).value_counts().to_csv(s.out/'tables/capture_scope_all_cell_counts.tsv',sep='\t')
    s.save(fig,'ED05_capture_scope_UMAP','Capture-specific atlas structure',
      'Every recorded capture scope in the unchanged atlas UMAP, using up to 60,000 major-balanced display cells per scope, fixed seed20261014. Full scope counts use all 2,670,057 cells. Enriched immune, stimulated lymphocyte, selected nuclei, serial biopsy and infusion-product captures are distinct observation contexts; these facets do not estimate comparable population prevalences.')
    assay=s.obs.assay.astype(str).to_numpy();assays=sorted(set(assay))
    colours=dict(zip(assays,STUDY_COLORS));assert len(assays)<=len(STUDY_COLORS)
    fig,axes=s.figure(165);fig.subplots_adjust(left=.13,right=.91,bottom=.22,wspace=.65,hspace=.65)
    s.scatter(axes[0,0],s.display,assay,colours,'Recorded assay / technology')
    from matplotlib.patches import Patch
    fig.legend(handles=[Patch(color=colours[x],label=x) for x in assays],ncol=3,
               loc='lower center',bbox_to_anchor=(.5,.02),frameon=False)
    table=pd.crosstab(s.obs.dataset,s.obs.assay).reindex(index=s.datasets,columns=assays,fill_value=0)
    frac=table.divide(table.sum(axis=1),axis=0)*100
    im=axes[0,1].imshow(frac,aspect='auto',cmap='cividis',vmin=0,vmax=100)
    axes[0,1].set_yticks(range(9));axes[0,1].set_yticklabels([study(x) for x in s.datasets])
    axes[0,1].set_xticks(range(len(assays)));axes[0,1].set_xticklabels(assays,rotation=45,ha='right')
    axes[0,1].set_title('Within-source technology');fig.colorbar(im,ax=axes[0,1],fraction=.04,pad=.02).set_label('Archived cells (%)')
    statuses=np.where(s.major=='Unresolved','Unresolved',np.where(s.major=='T_NK_boundary','T/NK boundary','Other major proposals'))
    pal={'Unresolved':'#D55E00','T/NK boundary':'#0072B2','Other major proposals':'#D6D6D6'}
    s.scatter(axes[1,0],s.display,statuses,pal,'Unresolved and boundary context')
    unresolved=pd.Series(s.major=='Unresolved').groupby(s.ds).mean().reindex(s.datasets)*100
    boundary=pd.Series(s.major=='T_NK_boundary').groupby(s.ds).mean().reindex(s.datasets)*100
    axes[1,1].barh(range(9),unresolved,color='#D55E00',label='Unresolved')
    axes[1,1].barh(range(9),boundary,left=unresolved,color='#0072B2',label='T/NK boundary')
    axes[1,1].set_yticks(range(9));axes[1,1].set_yticklabels([study(x) for x in s.datasets])
    axes[1,1].set_xlabel('Within-study archived cells (%)');axes[1,1].set_title('Full-cell annotation uncertainty')
    axes[1,1].legend(frameon=False)
    table.to_csv(s.out/'tables/assay_study_all_cell_counts.tsv',sep='\t')
    s.save(fig,'ED06_assay_and_uncertainty','Technology and annotation uncertainty',
      f'Recorded assays are retained exactly, including any source-specific or unknown labels. UMAP uses {len(s.display):,} study-balanced display points, fixed seed20261014. Assay proportions and unresolved/boundary fractions are full-cell counts. The 171,341 unresolved and 12,409 T/NK-boundary cells remain archived; their labels have not been forced or equated to malignant or normal truth.')
    s.finish()


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for k in ['atlas','run','stats','output']:p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--reuse-render',type=Path)
    main(p.parse_args())
