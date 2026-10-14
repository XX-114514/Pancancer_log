#!/usr/bin/env python3
"""Read-only manuscript atlas figures; all-cell summaries and explicit display samples."""
import argparse, gc, hashlib, html, json, os, platform, shutil, traceback
from datetime import datetime, timezone
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
import networkx as nx
import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem

SEED=20261014
MAJORS=['T','NK','B','Plasma','Myeloid','Neutrophil','Mast','Epithelial','Fibroblast',
        'Endothelial','Pericyte_SMC','Erythroid','Glial','Neuron','Oligodendrocyte',
        'Melanocytic','T_NK_boundary','Unresolved']
COLORS=['#0072B2','#56B4E9','#AA4499','#CC79A7','#009E73','#88CCAA','#44AA99',
        '#D55E00','#E69F00','#332288','#999933','#882255','#117733','#88CCEE',
        '#661100','#CC6677','#777777','#C5C5C5']
MAJOR_COLORS=dict(zip(MAJORS,COLORS))
STUDY_COLORS=['#0072B2','#E69F00','#009E73','#CC79A7','#56B4E9','#D55E00','#332288','#999933','#666666']
STATUS_ORDER=['CNV_supported_malignant','CNV_uncertain','CNV_negative_candidate',
              'lineage_reference_not_tumor_candidate','not_evaluated_skipped','not_evaluated']
STATUS_COLORS=dict(zip(STATUS_ORDER,['#D55E00','#E69F00','#0072B2','#56B4E9','#D0D0D0','#777777']))
PANELS={'T':['CD3D','CD3E','TRAC'],'NK':['GNLY','KLRD1','NKG7'],
 'B':['MS4A1','CD79A','CD74'],'Plasma':['MZB1','JCHAIN','SDC1'],
 'Myeloid':['LYZ','LST1','CSF1R'],'Neutrophil':['CSF3R','FCGR3B','S100A8'],
 'Mast':['KIT','TPSAB1','CPA3'],'Epithelial':['EPCAM','KRT8','KRT19'],
 'Fibroblast':['COL1A1','DCN','LUM'],'Endothelial':['PECAM1','VWF','KDR'],
 'Pericyte_SMC':['RGS5','ACTA2','TAGLN'],'Erythroid':['HBB','HBA1','GYPA'],
 'Glial':['GFAP','AQP4','SLC1A3'],'Neuron':['SNAP25','SYT1','RBFOX3'],
 'Oligodendrocyte':['MBP','PLP1','MOG'],'Melanocytic':['PMEL','MLANA','MITF'],
 'T_NK_boundary':['CD3D','TRAC','GNLY'],'Unresolved':['EPCAM','LYZ','CD3D']}
DETAIL_PANELS={
 'T':['CD3D','CD3E','TRAC','CD4','CD8A','CD8B','TCF7','CCR7','IL7R','LEF1',
      'FOXP3','IL2RA','CTLA4','GZMK','GZMB','PRF1','PDCD1','TOX','HAVCR2','LAG3',
      'IFIT1','ISG15','MKI67','TOP2A'],
 'B':['MS4A1','CD79A','CD79B','CD37','CD74','HLA-DRA','MZB1','JCHAIN','CD3D','TRAC','LST1','GNLY'],
 'Myeloid':['LYZ','LST1','CSF1R','CD14','FCGR3A','S100A8','S100A9','C1QA','C1QB',
            'C1QC','APOE','TREM2','SPP1','CD68','CD163','IL1B','HLA-DRA','FCER1A',
            'CD1C','CLEC9A','CST3','TCF4','LILRA4'],
 'NK':['GNLY','NKG7','KLRD1','FCGR3A','NCAM1','XCL1','XCL2','KLRB1','CD3D','TRAC']}

def now():return datetime.now(timezone.utc).isoformat()
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
def yes(v):return v.astype(str).str.lower().isin(['true','1'])
def study(v):return 'CAR-T product' if str(v).startswith('CXG_') else str(v)
def label(v):return str(v).replace('_',' ')
def short_node(v):
    bits=str(v).rsplit('_',2)
    return bits[-1]+' '+bits[-2] if len(bits)==3 else str(v)
def event(name,**data):print(json.dumps(dict(utc=now(),phase=name,**data)),flush=True)


class Suite:
    def __init__(self,a):
        self.a=a;self.out=a.output;self.out.mkdir(mode=0o700,exist_ok=False)
        self.figdir=self.out/'figures';self.figdir.mkdir();(self.out/'tables').mkdir()
        self.entries=[];self.excluded_genes=set()
        plt.rcParams.update({'font.family':'DejaVu Sans','font.size':6.2,'axes.labelsize':6.2,
          'axes.titlesize':7,'xtick.labelsize':5.6,'ytick.labelsize':5.6,'legend.fontsize':5.7,
          'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.linewidth':.55,
          'xtick.major.width':.5,'ytick.major.width':.5,'lines.linewidth':.7,
          'savefig.facecolor':'white','figure.facecolor':'white','axes.spines.top':False,
          'axes.spines.right':False,'axes.grid':False})
        self.g=zarr.open_group(str(a.atlas),mode='r')
        self.obs=read_elem(self.g['obs']);self.var=read_elem(self.g['var'])
        self.xy=np.asarray(self.g['obsm/X_umap']);assert np.isfinite(self.xy).all()
        assert len(self.obs)==2670057 and len(self.var)==49563 and self.obs.index.is_unique
        self.catalog=pd.read_csv(a.run/'tables/subcluster_catalog_refined_private.tsv',sep='\t',keep_default_na=False)
        self.stats_nodes=pd.read_csv(a.stats/'nodes_private.tsv',sep='\t',keep_default_na=False)
        genes=pd.read_csv(a.stats/'genes_private.tsv',sep='\t',index_col=0,keep_default_na=False)
        assert genes.index.equals(self.var.index)
        assert self.stats_nodes.node_id.tolist()==self.catalog.node_id.tolist()
        assert self.stats_nodes.subcluster.tolist()==self.catalog.subcluster.tolist()
        assert set(self.catalog.subcluster)==set(self.obs.subcluster.astype(str))
        assert self.obs.groupby('subcluster',observed=True).size().reindex(self.catalog.subcluster).tolist()==self.catalog.n_cells.tolist()
        self.arr={k:np.load(a.stats/n,mmap_mode='r') for k,n in
          [('sum','sum_CP10K.npy'),('positive','positive_counts.npy'),('measured','measured_counts.npy')]}
        assert all(v.shape==(270,49563) for v in self.arr.values())
        symbols=self.var.gene_symbol.astype(str)
        good=~symbols.duplicated(keep=False)&self.var.feature_identity_status.ne('ambiguous_source_feature')
        self.gene_index={str(symbols.iloc[i]):i for i in np.flatnonzero(good)}
        self.datasets=sorted(self.obs.dataset.astype(str).unique())
        self.ds_colors=dict(zip(self.datasets,STUDY_COLORS))
        self.major=self.obs.major.astype(str).to_numpy();self.ds=self.obs.dataset.astype(str).to_numpy()
        self.node=self.obs.subcluster.astype(str).to_numpy()
        self.bounds=[self.xy[:,0].min(),self.xy[:,0].max(),self.xy[:,1].min(),self.xy[:,1].max()]
        self.display=self.sample(np.arange(len(self.obs)),self.ds,200000)
        self.node_colors={n:plt.get_cmap('tab20')(i%20) if i<20 else plt.get_cmap('tab20b')((i-20)%20)
                          for i,n in enumerate(self.catalog.subcluster)}
        self.cm_index=json.loads((a.run/'manifests/CM_recovery_v3_public.json').read_text())
        assert self.cm_index['all_ready_scopes_completed']
        self.audit_marker_statistics()

    def sample(self,indices,groups,cap,seed=SEED):
        rng=np.random.default_rng(seed);frame=pd.DataFrame({'row':indices,'group':groups})
        each=max(1,cap//frame.group.nunique())
        picked=np.concatenate([rng.choice(x.row.to_numpy(),min(each,len(x)),replace=False)
            for _,x in frame.groupby('group',sort=True)])
        return picked[rng.permutation(len(picked))]

    def audit_marker_statistics(self):
        errors=[]
        for i,r in self.catalog.iterrows():
            if r.marker not in self.gene_index:continue
            j=self.gene_index[r.marker];den=float(self.arr['measured'][i,j])
            if den==0:raise ValueError('Named marker is unmeasured')
            mean=float(self.arr['sum'][i,j])/den;pct=float(self.arr['positive'][i,j])/den
            if abs(mean-float(r.marker_mean_CP10K))>1e-6 or abs(pct-float(r.marker_positive_fraction))>1e-9:
                raise ValueError('Marker sufficient statistics differ from final catalog')
            errors.append(abs(mean-float(r.marker_mean_CP10K)))
        self.marker_checks=dict(n_checked_named_markers=len(errors),maximum_mean_error=max(errors))

    def figure(self,h=150,cols=2,rows=2):
        fig,axes=plt.subplots(rows,cols,figsize=(183/25.4,h/25.4),squeeze=False)
        fig.subplots_adjust(left=.10,right=.96,bottom=.16,top=.93,wspace=.48,hspace=.56)
        for i,ax in enumerate(axes.flat):
            ax.text(-.12,1.035,chr(97+i),transform=ax.transAxes,fontsize=8,fontweight='bold')
        return fig,axes

    def save(self,fig,key,title,caption,main=False):
        width,height=fig.get_size_inches()*25.4
        assert abs(width-183)<1e-5 and height<=170.001
        fig.canvas.draw()
        texts=[t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text()]
        assert all(5<=t.get_fontsize()<=8.1 for t in texts)
        renderer=fig.canvas.get_renderer();bounds=fig.bbox
        not_drawn_ticks=set()
        for ax in fig.axes:
            for axis,limits in [(ax.xaxis,ax.get_xlim()),(ax.yaxis,ax.get_ylim())]:
                for tick in axis.get_major_ticks()+axis.get_minor_ticks():
                    if not min(limits)<=tick.get_loc()<=max(limits):
                        not_drawn_ticks.update([tick.label1,tick.label2])
        clipped=[]
        for t in texts:
            if t in not_drawn_ticks:continue
            b=t.get_window_extent(renderer)
            if b.x0<bounds.x0-1 or b.y0<bounds.y0-1 or b.x1>bounds.x1+1 or b.y1>bounds.y1+1:
                clipped.append(t.get_text())
        if clipped:raise ValueError('Figure text exceeds canvas: '+repr(clipped))
        for extension in ['pdf','svg','png']:
            path=self.figdir/(key+'.'+extension)
            assert not path.exists()
            fig.savefig(path,dpi=600)
        self.entries.append(dict(id=key,title=title,caption=caption,main=main,width_mm=183,
          height_mm=float(height),PNG_dpi=600,PDF_editable_text=True,SVG_editable_text=True))
        plt.close(fig);gc.collect();event('figure_saved',id=key)

    def scatter(self,ax,ids,values,palette,title,global_extent=True,legend=False):
        rng=np.random.default_rng(SEED+7);ids=np.asarray(ids);ids=ids[rng.permutation(len(ids))]
        v=np.asarray(values)[ids]
        colors=np.asarray([matplotlib.colors.to_rgba(palette[str(x)]) for x in v],dtype=np.float32)
        ax.scatter(self.xy[ids,0],self.xy[ids,1],c=colors,s=.12 if len(ids)>500000 else .4,
          alpha=.85,linewidths=0,rasterized=True)
        ax.set_title(title);ax.set_xlabel('UMAP 1 (a.u.)');ax.set_ylabel('UMAP 2 (a.u.)')
        ax.xaxis.set_major_locator(MaxNLocator(3));ax.yaxis.set_major_locator(MaxNLocator(3))
        if global_extent:ax.set_xlim(self.bounds[:2]);ax.set_ylim(self.bounds[2:])
        if legend:
            handles=[Line2D([],[],marker='o',linestyle='',color=palette[x],markersize=3,label=label(x)) for x in sorted(set(v))]
            ax.legend(handles=handles,loc='upper left',bbox_to_anchor=(1.02,1),frameon=False)

    def expression(self,node_rows,requested,pool=False):
        genes=[x for x in dict.fromkeys(requested) if x in self.gene_index]
        self.excluded_genes.update(set(requested)-set(genes))
        cols=[self.gene_index[x] for x in genes]
        sums=np.asarray(self.arr['sum'][:,cols]);positive=np.asarray(self.arr['positive'][:,cols]);measured=np.asarray(self.arr['measured'][:,cols])
        if pool:
            sums=np.vstack([sums[i].sum(axis=0) for i in node_rows]);positive=np.vstack([positive[i].sum(axis=0) for i in node_rows]);measured=np.vstack([measured[i].sum(axis=0) for i in node_rows])
        else:sums=sums[node_rows];positive=positive[node_rows];measured=measured[node_rows]
        mean=np.divide(sums,measured,out=np.full(sums.shape,np.nan),where=measured>0)
        pct=np.divide(positive,measured,out=np.full(sums.shape,np.nan),where=measured>0)
        assert np.nanmin(pct)>=0 and np.nanmax(pct)<=1
        colour=np.log1p(mean);top=np.nanmax(colour,axis=0)
        colour=np.divide(colour,top,out=np.zeros_like(colour),where=top>0)
        colour[measured==0]=np.nan
        return genes,mean,pct,colour

    def dots(self,fig,ax,rows,genes,rowlabels,pool=False):
        genes,mean,pct,colour=self.expression(rows,genes,pool)
        yy,xx=np.indices(mean.shape)
        width_pt=fig.get_figwidth()*72*ax.get_position().width*.94
        scale=min(6.7,width_pt/max(len(genes),1)*.82)**2
        dots=ax.scatter(xx.ravel(),yy.ravel(),s=(pct*scale).ravel(),c=colour.ravel(),cmap='cividis',vmin=0,vmax=1,linewidths=0)
        ax.set_xticks(np.arange(len(genes)));ax.set_xticklabels(genes,rotation=90,fontsize=5.2)
        ax.set_yticks(np.arange(len(rowlabels)));ax.set_yticklabels(rowlabels,fontsize=5.5)
        ax.set_xlim(-.8,len(genes)-.2);ax.set_ylim(len(rowlabels)-.4,-.8)
        ax.set_xlabel('Native marker gene');ax.set_ylabel('')
        cbar=fig.colorbar(dots,ax=ax,fraction=.025,pad=.015);cbar.set_label('Relative log mean CP10K')
        for spine in ax.spines.values():spine.set_visible(False)
        handles=[ax.scatter([],[],s=p*scale,color='#777777',label=str(int(p*100))+'%') for p in [.25,.5,.75,1]]
        position=ax.get_position()
        fig.legend(handles=handles,title='Detected / measured cells',ncol=4,frameon=False,
          loc='lower center',bbox_to_anchor=(position.x0+position.width/2,.075 if pool else .025),
          fontsize=5.5,title_fontsize=5.5,columnspacing=.7,handletextpad=.25)
        return mean,pct,genes

    def overview(self):
        fig,axes=self.figure(170);fig.subplots_adjust(bottom=.24,hspace=.65)
        self.scatter(axes[0,0],np.arange(len(self.obs)),self.major,MAJOR_COLORS,'Major proposals; all 2,670,057 cells')
        self.scatter(axes[0,1],self.display,self.ds,self.ds_colors,'Studies; balanced display')
        counts=self.obs.dataset.astype(str).value_counts().reindex(self.datasets)
        axes[1,0].barh(np.arange(9),counts/1000,color=[self.ds_colors[x] for x in self.datasets]);axes[1,0].set_yticks(np.arange(9));axes[1,0].set_yticklabels([study(x) for x in self.datasets]);axes[1,0].set_xlabel('Archived cells (thousands)');axes[1,0].set_title('Cohort size; full archive')
        eligible=yes(self.obs.clinical_validation_eligible).groupby(self.obs.dataset,observed=True).sum().reindex(self.datasets)
        axes[1,1].barh(np.arange(9),eligible/counts*100,color='#0072B2');axes[1,1].set_yticks(np.arange(9));axes[1,1].set_yticklabels([study(x) for x in self.datasets]);axes[1,1].set_xlim(0,100);axes[1,1].set_xlabel('Native endpoint-eligible cells (%)');axes[1,1].set_title('Source-specific endpoint coverage')
        mh=[Patch(color=MAJOR_COLORS[x],label=label(x)) for x in MAJORS]
        fig.legend(handles=mh,ncol=6,loc='lower center',bbox_to_anchor=(.5,.045),frameon=False,fontsize=5.5)
        self.save(fig,'Fig01_atlas_and_cohorts','An immune-treatment-associated atlas',
          f'a, All {len(self.obs):,} archived cells at unchanged atlas UMAP coordinates, shuffled drawing order. b, {len(self.display):,} study-balanced display cells, seed {SEED}; counts in c and eligibility in d use every cell. Study colours match cohort bars. Endpoint eligibility is source-specific and does not count independent patients. Major labels include proposals, unresolved and boundary categories. UMAP fit used 100,000 landmarks and transformed the remaining 2,570,057 cells.',True)

    def composition_and_markers(self):
        ct=pd.crosstab(self.obs.dataset,self.obs.major).reindex(index=self.datasets,columns=MAJORS,fill_value=0)
        ct.to_csv(self.out/'tables/study_major_all_cell_counts.tsv',sep='\t')
        fig,axes=self.figure(165,cols=1,rows=2);fig.subplots_adjust(left=.12,right=.89,bottom=.25,hspace=.65)
        left=np.zeros(9)
        for major in MAJORS:
            values=ct[major].to_numpy()/ct.sum(axis=1).to_numpy()*100
            axes[0,0].barh(np.arange(9),values,left=left,color=MAJOR_COLORS[major],height=.72);left+=values
        axes[0,0].set_yticks(np.arange(9));axes[0,0].set_yticklabels([study(x) for x in self.datasets]);axes[0,0].set_xlim(0,100);axes[0,0].set_xlabel('Within-study archived cells (%)');axes[0,0].set_title('Composition reflects capture scope')
        node_rows=[np.flatnonzero(self.catalog.major.eq(m)) for m in MAJORS]
        mean,pct,genes=self.dots(fig,axes[1,0],node_rows,sum([PANELS[m] for m in MAJORS],[]),[label(m) for m in MAJORS],True)
        axes[1,0].set_title('Multigene evidence from all-cell native statistics')
        pd.DataFrame(mean,index=MAJORS,columns=genes).to_csv(self.out/'tables/major_marker_mean_CP10K.tsv',sep='\t')
        pd.DataFrame(pct,index=MAJORS,columns=genes).to_csv(self.out/'tables/major_marker_detection_fraction.tsv',sep='\t')
        fig.legend(handles=[Patch(color=MAJOR_COLORS[m],label=label(m)) for m in MAJORS],ncol=6,
                   frameon=False,loc='lower center',bbox_to_anchor=(.5,.005),fontsize=5.2)
        self.save(fig,'Fig02_composition_and_identity_evidence','Composition and multigene identity evidence',
          'All archived cells contribute to composition and marker statistics. Dot area is the positive fraction among cells in sources that measured that gene; colour is log1p(mean native chosen-count CP10K), scaled to each gene maximum across displayed rows. Unmeasured genes remain missing. Capture-biased composition is descriptive. Proposals require multigene review; naming markers alone are not definitive identities.',True)

    def lineage(self,major,key,main=False):
        rows=np.flatnonzero(self.catalog.major.eq(major));cat=self.catalog.iloc[rows]
        ids=np.flatnonzero(self.major==major);cap=200000 if major=='T' else 40000
        sample=self.sample(ids,self.node[ids],cap)
        colours={n:self.node_colors[n] for n in cat.subcluster}
        fig,axes=self.figure(150,cols=2,rows=1);fig.subplots_adjust(left=.08,right=.88,bottom=.22,wspace=.7)
        self.scatter(axes[0,0],sample,self.node,colours,f'{label(major)} descriptive nodes',False)
        centres=[];numbers=[]
        for row in cat.itertuples():
            group=sample[self.node[sample]==row.subcluster]
            if len(group):
                xy=np.median(self.xy[group],axis=0)
                centres.append(xy);numbers.append(str(row.number))
        ax=axes[0,0];centres=np.asarray(centres)
        ax.get_xlim();ax.get_ylim();fig.canvas.draw()
        origin=ax.transAxes.inverted().transform(ax.transData.transform(centres))
        positions=origin.copy();rng=np.random.default_rng(SEED)
        for step in range(250):
            force=np.zeros_like(positions)
            for i in range(len(positions)):
                for j in range(i):
                    delta=(positions[i]-positions[j])/np.array([.085,.045]);distance=np.linalg.norm(delta)
                    if distance<1:
                        direction=delta/max(distance,1e-8) if distance>1e-8 else rng.normal(size=2)
                        push=direction*(1-distance)*.018
                        force[i]+=push;force[j]-=push
            positions=np.clip(positions+force+(origin-positions)*.001,.035,.965)
        for number,xy,position in zip(numbers,centres,positions):
            ax.annotate(number,xy=xy,xycoords='data',xytext=position,textcoords='axes fraction',fontsize=5,
                ha='center',va='center',bbox=dict(facecolor='white',edgecolor='none',alpha=.8,pad=.2),
                arrowprops=dict(arrowstyle='-',color='#888888',lw=.4),annotation_clip=True)
        genes=list(dict.fromkeys(DETAIL_PANELS.get(major,PANELS[major])+cat.marker.tolist()))
        self.dots(fig,axes[0,1],rows,genes,[str(r.number)+' '+r.marker for r in cat.itertuples()])
        axes[0,1].set_title('Own-node marker evidence; all cells')
        self.save(fig,key,f'{label(major)}: descriptive subclusters',
          f'{len(ids):,} archived {major} cells and {len(cat)} descriptive nodes; {len(sample):,} node-balanced display points, seed {SEED}, in the original atlas UMAP (no new lineage embedding). Numerical labels identify catalog numbers; labels are moved with leader lines for legibility while data coordinates stay unchanged. Complete labels and full-cell expression evidence are in the dotplot. Colour is relative per-gene log mean CP10K; area is detection/measured fraction. Fine candidates remain unpromoted; stress/sex markers are descriptors.',main)

    def states(self):
        cols=[x for x in self.obs if x.startswith('state_score_')];assert len(cols)==5
        ids=np.flatnonzero(self.major=='T');sample=self.sample(ids,self.ds[ids],200000)
        fig,axes=self.figure(160,cols=3,rows=2);fig.subplots_adjust(left=.09,right=.86,bottom=.17,wspace=.75,hspace=.65)
        for ax,key in zip(axes.flat,cols):
            values=self.obs[key].to_numpy();vmax=float(np.nanquantile(values[ids],.99))
            points=ax.scatter(self.xy[sample,0],self.xy[sample,1],c=values[sample],s=.3,
              cmap='cividis',vmin=0,vmax=max(vmax,1e-6),linewidths=0,rasterized=True)
            ax.set_title(label(key.removeprefix('state_score_')));ax.set_xlabel('UMAP 1 (a.u.)');ax.set_ylabel('UMAP 2 (a.u.)')
            ax.xaxis.set_major_locator(MaxNLocator(3));ax.yaxis.set_major_locator(MaxNLocator(3))
            cbar=fig.colorbar(points,ax=ax,fraction=.04,pad=.02);cbar.set_label('Score (a.u.)')
        t=self.obs.iloc[ids];means=t.groupby('subcluster',observed=True)[cols].mean().reindex(self.catalog.loc[self.catalog.major.eq('T'),'subcluster'])
        z=(means-means.mean())/means.std(ddof=0).replace(0,np.nan)
        ax=axes[1,2];im=ax.imshow(z,cmap='PuOr',vmin=-2,vmax=2,aspect='auto')
        ax.set_xticks(range(5));ax.set_xticklabels([label(x.removeprefix('state_score_')) for x in cols],rotation=90)
        ax.set_yticks(range(len(means)));ax.set_yticklabels([short_node(x) for x in means.index]);ax.set_title('Mean state score by T node')
        fig.colorbar(im,ax=ax,fraction=.04,pad=.02).set_label('Across-node z-score')
        means.to_csv(self.out/'tables/T_node_full_cell_state_means.tsv',sep='\t')
        self.save(fig,'Fig04_T_continuous_states','Continuous T-cell state axes',
          f'Five saved native panel-score axes, {len(sample):,} study-balanced T display cells. Each score colour range runs from zero to its all-T-cell 99th percentile; saturation is visual only. The heatmap uses all {len(ids):,} T cells and descriptive node means standardized separately within each axis, clipped visually at ±2. Scores are state axes, not mutually exclusive cell types or treatment-response tests.',True)

    def CNV(self):
        fig,axes=self.figure(160);fig.subplots_adjust(bottom=.22,hspace=.7)
        status=self.obs.malignant_status.astype(str).to_numpy()
        self.scatter(axes[0,0],np.arange(len(self.obs)),status,STATUS_COLORS,'CNV states; all archived cells')
        counts=pd.crosstab(self.obs.dataset,self.obs.malignant_status).reindex(index=self.datasets,columns=STATUS_ORDER,fill_value=0)
        left=np.zeros(9)
        for key in STATUS_ORDER:
            v=counts[key].to_numpy()/counts.sum(axis=1).to_numpy()*100
            axes[0,1].barh(np.arange(9),v,left=left,color=STATUS_COLORS[key]);left+=v
        axes[0,1].set_yticks(range(9));axes[0,1].set_yticklabels([study(x) for x in self.datasets]);axes[0,1].set_xlim(0,100);axes[0,1].set_xlabel('Archived cells (%)');axes[0,1].set_title('Evidence and non-evaluable capture')
        libs=self.obs[['dataset','sampleID','CNV_method_status']].drop_duplicates()
        assert len(libs)==400
        methods=pd.crosstab(libs.dataset,libs.CNV_method_status).reindex(index=self.datasets,columns=['success','skipped'],fill_value=0)
        x=np.arange(9);axes[1,0].barh(x,methods.success,color='#0072B2',label='Success');axes[1,0].barh(x,methods.skipped,left=methods.success,color='#B8B8B8',label='Skipped')
        axes[1,0].set_yticks(x);axes[1,0].set_yticklabels([study(s) for s in self.datasets]);axes[1,0].set_xlabel('Original libraries (n)');axes[1,0].set_title('149 success; 251 skipped');axes[1,0].legend(frameon=False)
        references=self.obs.loc[self.obs.tumor_reference_cell,['sampleID','cnv_score']].groupby('sampleID',observed=True).cnv_score.median()
        candidates=self.obs.loc[self.obs.cnv_evaluable,['sampleID','dataset','cnv_score']].groupby(['sampleID','dataset'],observed=True).cnv_score.median().reset_index()
        candidates['reference']=candidates.sampleID.map(references)
        valid=np.isfinite(candidates.reference)&candidates.reference.gt(0)
        candidates=candidates.loc[valid].copy();candidates['burden']=np.log2((candidates.cnv_score+1e-8)/(candidates.reference+1e-8))
        ds=sorted(candidates.dataset.astype(str).unique());rng=np.random.default_rng(SEED)
        ax=axes[1,1]
        for i,s in enumerate(ds):
            y=candidates.loc[candidates.dataset.astype(str).eq(s),'burden'].to_numpy()
            ax.scatter(i+rng.uniform(-.12,.12,len(y)),y,s=5,color=self.ds_colors[s],alpha=.6,rasterized=True)
            ax.plot([i-.23,i+.23],[np.median(y)]*2,color='black',lw=1)
        ax.axhline(0,color='#999999',lw=.5);ax.set_xticks(range(len(ds)));ax.set_xticklabels([study(s) for s in ds],rotation=45,ha='right');ax.set_ylabel('log2(candidate/reference median burden)');ax.set_title('Per-library burden; descriptive')
        fig.legend(handles=[Patch(color=STATUS_COLORS[s],label=label(s).replace('CNV ','').replace('not evaluated skipped','Skipped')) for s in STATUS_ORDER],ncol=3,loc='lower center',bbox_to_anchor=(.5,.015),frameon=False,fontsize=5.5)
        counts.to_csv(self.out/'tables/study_CNV_status_all_cells.tsv',sep='\t');methods.to_csv(self.out/'tables/study_CNV_library_methods.tsv',sep='\t')
        self.save(fig,'Fig05_CNV_evidence','Malignancy candidates and evaluation scope',
          'All-cell atlas CNV states and source-specific cell counts; original-library methods, not cell counts, define panel c. Panel d shows each evaluable library median candidate burden divided by its own reference median (log2), with library median across points as black line; no new threshold or hypothesis test. 52,657 single-method supported malignant candidates, 57,958 uncertain cells. Skips and negative candidates do not establish normal identity.',True)

    def CM_detail(self,result,key,main=False):
        scope=result['scope'];folder=self.a.run/'CM'/scope/result['discovery_attempt']
        W=pd.read_csv(folder/'W_private.tsv',sep='\t',index_col=0)
        H=pd.read_csv(folder/'H_patient_block_private.tsv',sep='\t',index_col=0)
        assert list(W.columns)==list(H.index) and np.isfinite(W.to_numpy()).all()
        factor_mass=W.sum(axis=0);wn=W.divide(factor_mass,axis=1)
        chosen=set()
        for c in wn:chosen.update(wn[c].nlargest(3 if len(wn.columns)>4 else 8).index)
        rows=sorted(chosen,key=lambda n:(str(self.catalog.set_index('subcluster').loc[n,'major']),str(n)))
        fig=plt.figure(figsize=(183/25.4,168/25.4))
        grid=fig.add_gridspec(3,2,left=.16,right=.87,bottom=.12,top=.90,wspace=.75,hspace=.95)
        axes=np.array([[fig.add_subplot(grid[:2,0]),fig.add_subplot(grid[0,1])],
                       [fig.add_subplot(grid[2,0]),fig.add_subplot(grid[1:,1])]])
        for i,ax in enumerate(axes.flat):ax.text(-.12,1.08,chr(97+i),transform=ax.transAxes,fontsize=8,fontweight='bold')
        ax=axes[0,0];im=ax.imshow(wn.loc[rows],aspect='auto',cmap='cividis',vmin=0,vmax=max(float(wn.loc[rows].to_numpy().max()),1e-8))
        ax.set_yticks(range(len(rows)));ax.set_yticklabels([short_node(x) for x in rows],fontsize=5)
        ax.set_xticks(range(len(W.columns)));ax.set_xticklabels(W.columns,rotation=90);ax.set_title('Selected node loadings');fig.colorbar(im,ax=ax,fraction=.04,pad=.02).set_label('W / factor mass')
        metadata=pd.read_csv(self.a.run/'CM'/scope/'inputs/sample_metadata_private.tsv',sep='\t',keep_default_na=False)
        mapping=metadata[['independent_unit','dataset']].drop_duplicates();assert not mapping.independent_unit.duplicated().any()
        source=mapping.set_index('independent_unit').dataset.reindex(H.columns);assert source.notna().all()
        contribution=H.mul(factor_mass,axis=0);contribution=contribution.divide(contribution.sum(axis=0),axis=1)
        med=contribution.T.groupby(source.to_numpy()).median().T
        ax=axes[0,1];im=ax.imshow(med,aspect='auto',cmap='cividis',vmin=0,vmax=1)
        ax.set_xticks(range(len(med.columns)));ax.set_xticklabels([study(x) for x in med.columns],rotation=45,ha='right');ax.set_yticks(range(len(H)));ax.set_yticklabels(H.index);ax.set_title('Source median factor contribution');fig.colorbar(im,ax=ax,fraction=.04,pad=.02).set_label('Relative contribution')
        profile=pd.read_csv(folder/'rank_survey.tsv',sep='\t');axes[1,0].plot(profile['rank'],profile.cophenetic,'o-',ms=3,color='#0072B2');axes[1,0].axvline(result['selected_rank'],color='#D55E00',lw=.8);axes[1,0].set_xlabel('Surveyed NMF rank');axes[1,0].set_ylabel('Consensus cophenetic');axes[1,0].set_title('Frozen rank survey');axes[1,0].xaxis.set_major_locator(MaxNLocator(integer=True,nbins=5))
        try:edges=pd.read_csv(folder/'global_edge.tsv',sep='\t')
        except pd.errors.EmptyDataError:edges=pd.DataFrame(columns=['subCluster1','subCluster2','correlation'])
        try:nodes=pd.read_csv(folder/'global_node.tsv',sep='\t')
        except pd.errors.EmptyDataError:nodes=pd.DataFrame()
        ax=axes[1,1];ax.set_title('Supported composition covariance')
        if len(edges):
            graph=nx.from_pandas_edgelist(edges,'subCluster1','subCluster2',edge_attr='correlation')
            pos=nx.spring_layout(graph,seed=SEED,weight='correlation')
            nx.draw_networkx_edges(graph,pos,ax=ax,alpha=.4,width=[float(x['correlation'])*1.5 for _,_,x in graph.edges(data=True)],edge_color='#777777')
            lookup=self.catalog.set_index('subcluster');nodecolors=[MAJOR_COLORS[str(lookup.loc[n,'major'])] for n in graph]
            nx.draw_networkx_nodes(graph,pos,ax=ax,node_size=12,node_color=nodecolors,linewidths=0)
            names=[n for n in rows if n in graph];ax.get_xlim();ax.get_ylim();fig.canvas.draw()
            anchors=np.array([pos[n] for n in names])
            origin=ax.transAxes.inverted().transform(ax.transData.transform(anchors));positions=origin.copy()
            rng=np.random.default_rng(SEED)
            for step in range(250):
                force=np.zeros_like(positions)
                for i in range(len(positions)):
                    for j in range(i):
                        delta=(positions[i]-positions[j])/np.array([.085,.045]);distance=np.linalg.norm(delta)
                        if distance<1:
                            direction=delta/max(distance,1e-8) if distance>1e-8 else rng.normal(size=2)
                            push=direction*(1-distance)*.018;force[i]+=push;force[j]-=push
                positions=np.clip(positions+force+(origin-positions)*.001,.035,.965)
            for name,anchor,position in zip(names,anchors,positions):
                ax.annotate(str(lookup.loc[name,'number']),xy=anchor,xycoords='data',xytext=position,
                    textcoords='axes fraction',fontsize=5,ha='center',va='center',
                    bbox=dict(facecolor='white',edgecolor='none',alpha=.8,pad=.15),
                    arrowprops=dict(arrowstyle='-',color='#888888',lw=.4))
            majors=sorted({str(lookup.loc[n,'major']) for n in graph})
            fig.legend(handles=[Patch(color=MAJOR_COLORS[m],label=label(m)) for m in majors],
                frameon=False,ncol=6,loc='lower center',bbox_to_anchor=(.5,.015),fontsize=5.2)
            ax.text(.02,.02,f'{len(graph)} nodes; {len(edges)} edges',transform=ax.transAxes,fontsize=5.5)
        else:ax.text(.5,.5,'No supported network\nNMF factors retained',ha='center',va='center',transform=ax.transAxes,fontsize=7)
        ax.set_axis_off()
        fig.suptitle(scope.replace('_',' ')+f" | K={result['selected_rank']}; supported member modules={result['n_modules_with_network_members']}",fontsize=7,y=.995)
        self.save(fig,key,f'CM scope: {scope}',
          f"Actual saved {scope} nsNMF, {result['n_nodes']} nodes, rank {result['selected_rank']}, {result['n_modules_with_network_members']} factors with supported network members. Rows are the union of top loading nodes per factor (three per factor for K>4, eight otherwise), at most 36 readable rows; all rows remain in saved W. W is scaled by factor mass; source median H contributions use matching inverse scaling, then within-block fractions. Blocks include explicitly exploratory unverified libraries where applicable. Numeric network labels refer to catalog node numbers; marker names appear in loading rows or the full catalog. Network edges are saved source-adjusted composition correlations after original FDR/specificity/support gates; layout is illustrative, not spatial or ligand-receptor evidence. Rank-boundary/LOSO-unestimable results remain in source reports. No new clinical tests.",main)

    def scope_overview(self):
        results=self.cm_index['results'];fig,axes=self.figure(155,cols=2,rows=1);fig.subplots_adjust(left=.28,right=.96,bottom=.15,wspace=.45)
        names=[r['scope'].replace('_',' ') for r in results];y=np.arange(len(names))
        rank=np.array([r['selected_rank'] for r in results]);supported=np.array([r['n_modules_with_network_members'] for r in results])
        axes[0,0].barh(y,rank,color='#DDDDDD',label='NMF rank');axes[0,0].barh(y,supported,color='#0072B2',height=.5,label='With network members');axes[0,0].set_yticks(y);axes[0,0].set_yticklabels(names,fontsize=5.2);axes[0,0].set_xlabel('Factors (n)');axes[0,0].set_title('All 15 completed scopes');axes[0,0].legend(frameon=False,loc='lower right')
        axes[0,1].barh(y,[r['n_network_edges'] for r in results],color=['#0072B2' if x else '#BBBBBB' for x in supported]);axes[0,1].set_yticks(y);axes[0,1].set_yticklabels([]);axes[0,1].set_xlabel('Supported global edges (n)');axes[0,1].set_title('Positive and unsupported scopes')
        self.save(fig,'ED01_CM_scope_coverage','Complete CM discovery scope',
          'All 15 completed composition scopes: 10 have supported network-member modules, five do not; another six prepared scopes were insufficient. Selected NMF rank and supported-member module count are distinct. Scopes overlap and their ranks must not be summed as independent modules. No unsupported result is hidden.')

    def technical(self):
        fig,axes=self.figure(160);fig.subplots_adjust(bottom=.20)
        ids=self.display
        for row,name in enumerate(['X_pca','X_pca_harmony']):
            coords=np.asarray(self.g['obsm'][name][:,:2]);values=self.ds[ids]
            colours=np.asarray([matplotlib.colors.to_rgba(self.ds_colors[x]) for x in values])
            axes[row,0].scatter(coords[ids,0],coords[ids,1],c=colours,s=.3,rasterized=True,linewidths=0)
            colours=np.asarray([matplotlib.colors.to_rgba(MAJOR_COLORS[x]) for x in self.major[ids]])
            axes[row,1].scatter(coords[ids,0],coords[ids,1],c=colours,s=.3,rasterized=True,linewidths=0)
            for ax in axes[row]:ax.set_xlabel('Component 1 (a.u.)');ax.set_ylabel('Component 2 (a.u.)')
            axes[row,0].set_title(('PCA' if row==0 else 'Harmony')+' by study');axes[row,1].set_title(('PCA' if row==0 else 'Harmony')+' by major proposal')
        fig.legend(handles=[Patch(color=self.ds_colors[x],label=study(x)) for x in self.datasets],ncol=5,loc='lower center',frameon=False)
        self.save(fig,'ED02_PCA_Harmony','Batch and biological structure diagnostics',
          f'Identical {len(ids):,} study-balanced display cells before and after Harmony in the first two component coordinates, not a before/after UMAP. All RNA counts remain unchanged. The original 10,000-cell neighbor diagnostic reported entropy 0.4127→0.8077 and marker-proposal purity 0.8223→0.8136; these are descriptive indicators, not a ground-truth annotation or proof of absent overcorrection.')
        fig,axes=self.figure(130,cols=2,rows=1);fig.subplots_adjust(left=.12,right=.96,bottom=.28,wspace=.55)
        lib=self.obs.groupby(['dataset','sampleID'],observed=True).native_analysis_denominator.median()
        rng=np.random.default_rng(SEED)
        for i,s in enumerate(self.datasets):
            vals=np.log10(lib.loc[s].to_numpy());axes[0,0].scatter(i+rng.uniform(-.15,.15,len(vals)),vals,s=5,color=self.ds_colors[s],alpha=.65)
            axes[0,0].plot([i-.25,i+.25],[np.median(vals)]*2,color='black',lw=1)
        axes[0,0].set_xticks(range(9));axes[0,0].set_xticklabels([study(x) for x in self.datasets],rotation=90);axes[0,0].set_ylabel('log10(median chosen RNA counts / library)');axes[0,0].set_title('All 400 original libraries')
        correction=yes(self.obs.analysis_correction_applied).groupby(self.obs.dataset,observed=True).mean().reindex(self.datasets)*100
        excluded=self.obs.annotation_exclude_from_primary.groupby(self.obs.dataset,observed=True).mean().reindex(self.datasets)*100
        x=np.arange(9);axes[0,1].bar(x-.18,correction,width=.36,color='#0072B2',label='Correction adopted');axes[0,1].bar(x+.18,excluded,width=.36,color='#E69F00',label='Primary annotation exclusion');axes[0,1].set_xticks(x);axes[0,1].set_xticklabels([study(x) for x in self.datasets],rotation=90);axes[0,1].set_ylabel('Archived cells (%)');axes[0,1].set_title('Available correction / exclusion provenance');axes[0,1].legend(frameon=False)
        self.save(fig,'ED03_depth_and_provenance','Available depth and correction provenance',
          'Left: per-original-library medians of the complete native chosen-count denominator; points are libraries and black segments are study medians. Right: full-cell correction-adopted and primary-fit exclusion fractions. This is available provenance, not a replacement for source QC; original mitochondrial/nGene fields are not universally available and have not been fabricated. Prior-QC sources without new ambient assessment retain their documented status.')

    def facets(self):
        fig,axes=self.figure(170,cols=3,rows=3);fig.subplots_adjust(left=.07,right=.97,bottom=.10,hspace=.65,wspace=.45)
        for ax,s in zip(axes.flat,self.datasets):
            ids=np.flatnonzero(self.ds==s);sample=self.sample(ids,self.major[ids],60000)
            self.scatter(ax,sample,self.major,MAJOR_COLORS,study(s))
            ax.text(.02,.02,f'n={len(ids):,}; display={len(sample):,}',transform=ax.transAxes,fontsize=5)
        self.save(fig,'ED04_all_cohort_UMAP','Every cohort in the shared atlas coordinates',
          f'All nine sources, each with unchanged common UMAP limits. Up to 60,000 major-balanced display points per source, seed {SEED}, preserving rare observed classes; this sampling is not used for composition proportions. Titles map to GEO IDs and the CAR-T CXG source in the cohort registry.')

    def finish(self):
        self.entries.sort(key=lambda e:(not e['main'],e['id']))
        gallery=['<!doctype html><html><head><meta charset="utf-8"><title>Immune atlas manuscript figures</title><style>body{font:16px system-ui;max-width:1100px;margin:32px auto;padding:12px;color:#17212b}img{width:100%;border:1px solid #ddd}article{margin:32px 0;padding:16px;background:#fafafa}a{margin-right:14px}p{line-height:1.55}</style></head><body><h1>Immune atlas: manuscript figure suite</h1><p>Full archive: 2,670,057 cells; nine sources; 400 original libraries. Major/fine identities and single-method CNV remain proposals. New CMs describe composition covariance. See captions for sampling, measurement and scope.</p>']
        for e in self.entries:
            gallery.append('<article><h2>'+html.escape(e['id']+' · '+e['title'])+'</h2><p>'+html.escape(e['caption'])+'</p><p>'+''.join('<a href="figures/'+e['id']+'.'+x+'">'+x.upper()+'</a>' for x in ['pdf','svg','png'])+'</p><img loading="lazy" src="figures/'+e['id']+'.png"></article>')
        gallery.append('</body></html>');(self.out/'index.html').write_text('\n'.join(gallery))
        captions=['# Manuscript figure legends','']
        for e in self.entries:captions.extend(['## '+e['id']+' — '+e['title'],'',e['caption'],''])
        (self.out/'FIGURE_LEGENDS.md').write_text('\n'.join(captions))
        files=[]
        from PIL import Image
        for p in sorted(self.figdir.iterdir()):
            item=dict(name=p.name,size_bytes=p.stat().st_size,sha256=sha(p))
            if p.suffix=='.png':
                with Image.open(p) as im:
                    item['pixel_size']=list(im.size);item['DPI']=list(im.info.get('dpi',[]));assert im.width>=4200
            if p.suffix=='.svg':assert '<text' in p.read_text()
            files.append(item)
        report=dict(observed_utc=now(),status='FIGURES_RENDERED_NUMERICAL_AND_FILE_GATES_PASSED_VISUAL_REVIEW_PENDING',
          n_atlas_cells=len(self.obs),n_sources=9,n_original_libraries=400,n_nodes=270,
          n_main_figures=sum(e['main'] for e in self.entries),n_figure_sets=len(self.entries),
          figures=self.entries,files=files,full_native_statistics_marker_audit=self.marker_checks,
          excluded_missing_or_ambiguous_requested_markers=sorted(self.excluded_genes),
          input_release_sha256=sha(self.a.run/'manifests/immune_authority_release_public.json'),
          producer_sha256=sha(Path(__file__)),seed=SEED,patient_or_cell_IDs_rendered=False,
          original_archive_read_only=True,no_new_clinical_hypothesis_tests=True,
          main_UMAP_all_cells=True,detail_display_sampling_explicit=True,
          vector_text_and_600dpi_raster_verified=True,visual_review_complete=False,
          python_version=platform.python_version(),matplotlib_version=matplotlib.__version__)
        with (self.out/'figure_manifest_public.json').open('x') as f:json.dump(report,f,indent=2)
        event('suite_rendered',n_figures=len(self.entries))


def main(a):
    suite=Suite(a)
    preserved=['Fig01_atlas_and_cohorts','Fig05_CNV_evidence','ED01_CM_scope_coverage',
               'ED02_PCA_Harmony','ED03_depth_and_provenance','ED04_all_cohort_UMAP']
    if a.reuse_render:
        old=json.loads((a.reuse_render/'figure_manifest_public.json').read_text())
        for e in old['figures']:
            if e['id'] not in preserved:continue
            suite.entries.append(e)
            for ext in ['pdf','svg','png']:
                src=a.reuse_render/'figures'/(e['id']+'.'+ext);dst=suite.figdir/src.name
                shutil.copy2(src,dst);assert sha(src)==sha(dst)
        assert len(suite.entries)==6
        for name in ['study_CNV_status_all_cells.tsv','study_CNV_library_methods.tsv']:
            src=a.reuse_render/'tables'/name;dst=suite.out/'tables'/name
            shutil.copy2(src,dst);assert sha(src)==sha(dst)
        event('reused_exact_unaffected_figures',n_figures=6)
    else:suite.overview();suite.CNV();suite.scope_overview();suite.technical();suite.facets()
    suite.composition_and_markers();suite.lineage('T','Fig03_T_descriptive_subclusters',True)
    suite.states()
    whole=next(r for r in suite.cm_index['results'] if r['scope']=='whole_TME')
    suite.CM_detail(whole,'Fig06_whole_TME_covariance',True)
    for i,m in enumerate(MAJORS):suite.lineage(m,f'ED_lineage_{i+1:02d}_{m}')
    for i,r in enumerate(suite.cm_index['results']):suite.CM_detail(r,f'ED_CM_{i+1:02d}_{r["scope"]}')
    suite.finish()


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for k in ['atlas','run','stats','output']:p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--reuse-render',type=Path)
    a=p.parse_args()
    try:main(a)
    except Exception:
        event('failed',traceback=traceback.format_exc());raise
