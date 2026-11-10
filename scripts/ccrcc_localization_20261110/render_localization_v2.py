"""Server-only raster figures from fixed identities, no row/coordinate exports."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm


def main(o):
    figdir=o/'figures_v2';figdir.mkdir(exist_ok=False)
    x=np.load(o/'raster_PRIVATE.npz');step=int(x['step']);qc=x['QC'];origin=x['origin']
    extent=(0,qc.shape[1]*step/1000,qc.shape[0]*step/1000,0)
    plt.rcParams.update({'font.size':9,'axes.titlesize':10})
    def density_plot(panels,name,title):
        fig,axes=plt.subplots(2,3,figsize=(13,8.4),layout='constrained')
        for ax,(key,label) in zip(axes.ravel(),panels):
            a=np.where(qc>=10,np.log1p(x[key]),np.nan)
            im=ax.imshow(a,extent=extent,origin='upper',cmap='magma',vmin=0,interpolation='nearest')
            ax.set_title(label);ax.set_xlabel('Tissue x (mm)');ax.set_ylabel('Tissue y (mm)')
            fig.colorbar(im,ax=ax,label='log(1 + cells / 100 x 100 um bin)',shrink=.8)
        fig.suptitle(title+'\nOriginal identities; 100 um bins; bins with <10 total QC cells hidden')
        fig.savefig(figdir/name,dpi=170);plt.close(fig)
    density_plot([(k,label) for k,label in [('Blood_EC','Blood endothelial'),('NK','NK'),('B','B'),('Mural','Mural'),('T','T'),('Lymphatic','Lymphatic')]],'01_existing_identity_density.png','ccRCC: tissue composition on unchanged identities')
    density_plot([('Mural_CCL19','CCL19+ Mural'),('B_CCR7','CCR7+ B'),('T_CCR7','CCR7+ T'),('Mural_CXCL12','CXCL12+ Mural'),('B_CXCR4','CXCR4+ B'),('T_CXCR4','CXCR4+ T')],'02_sender_receiver_marker_density.png','ccRCC: measured ligand / receptor location (RNA >=1 UMI)')
    cm=pd.read_csv(o/'CM_regional_PRIVATE.tsv.gz',sep='\t');lr=pd.read_csv(o/'LR_regional_PRIVATE.tsv.gz',sep='\t')
    target=[('EC_NK_edges','EC - NK edges'),('EC_B_edges','EC - B edges'),('same_EC_B_and_NK','One EC with B and NK')]
    panels=[]
    for rel,title in target:
        a=cm[(cm.relation==rel)&(cm.radius_um==50)&(cm.tile_um==250)];panels.append((a,title))
    for sender,receiver,ligand,receptor in [('Blood_EC','T','CCL19','CCR7'),('Mural','B','CCL19','CCR7'),('Mural','T','CCL19','CCR7'),('Mural','B','CXCL12','CXCR4'),('Mural','T','CXCL12','CXCR4')]:
        a=lr[(lr.sender==sender)&(lr.receiver==receiver)&(lr.ligand==ligand)&(lr.radius_um==50)&(lr.tile_um==250)]
        panels.append((a,f'{ligand}+ {sender} -> {receptor}+ {receiver}'))
    # Last panel states the interpretation limit, keeping all five prior candidates visible.
    fig,axes=plt.subplots(3,3,figsize=(13,11.5),layout='constrained')
    allblocks=pd.concat([cm[['block_x','block_y']],lr[['block_x','block_y']]]).drop_duplicates()
    bmin=allblocks.min().to_numpy(int);bmax=allblocks.max().to_numpy(int)
    shape=bmax-bmin+1;ex=((bmin[0]*500-origin[0])/1000,((bmax[0]+1)*500-origin[0])/1000,((bmax[1]+1)*500-origin[1])/1000,(bmin[1]*500-origin[1])/1000)
    for ax,(a,title) in zip(axes.ravel(),panels):
        a=a[a.opportunities>=20];grid=np.full(shape[::-1],np.nan)
        grid[a.block_y.to_numpy(int)-bmin[1],a.block_x.to_numpy(int)-bmin[0]]=100*a.excess.to_numpy()/a.opportunities.to_numpy()
        vmax=max(float(np.nanmax(np.abs(grid))) if np.isfinite(grid).any() else 0,1e-3)
        im=ax.imshow(grid,extent=ex,origin='upper',cmap='RdBu_r',norm=TwoSlopeNorm(0,-vmax,vmax),interpolation='nearest')
        ax.set_xlim(extent[0],extent[1]);ax.set_ylim(extent[2],extent[3]);ax.set_title(title);ax.set_xlabel('Tissue x (mm)');ax.set_ylabel('Tissue y (mm)');fig.colorbar(im,ax=ax,label='Observed - exact expected (pp)',shrink=.7)
    axes.ravel()[-1].axis('off');axes.ravel()[-1].text(.04,.95,'Descriptive localization only\n\n50 um graph cutoff; 250 um shuffle grid\n500 um display blocks, >=20 opportunities\n\nRed: excess over original local null\nBlue: lower than original local null\n\nEach panel has its own colour scale.\nBlocks are not independent patients.\nNo new regional P values.\nNot a reconstructed vessel or histology.',va='top',transform=axes.ravel()[-1].transAxes)
    fig.suptitle('ccRCC: where the pre-existing spatial effects arise\nExact expectations of the unchanged null; identities and graph fixed')
    fig.savefig(figdir/'03_exact_null_regional_effects.png',dpi=170);plt.close(fig)
    manifest={p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in figdir.glob('*.png')}
    (figdir/'manifest.json').write_text(json.dumps({'raster_bin_um':100,'effect_block_um':500,'all_per_cell_data_private':True,'all_figures_share_same_coordinate_origin_and_extent':True,'files':manifest},indent=2)+'\n');print(json.dumps(manifest),flush=True)


if __name__=='__main__':main(Path(sys.argv[1]))
