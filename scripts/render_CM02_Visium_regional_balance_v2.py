"""V2 separates measured mid-scores from grey unscorable backgrounds."""
import argparse,json,html
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource

def main(a):
 budget=resource(a.root,64*1024**2)
 proof=json.loads((a.audit/'independent_audit_public.json').read_text());assert proof['n_slots']==76 and proof['producer_terminal_sha256']==sha(a.producer/'terminal_public.json')
 assert not a.output.exists();a.output.mkdir(mode=0o700)
 d=pd.read_csv(a.producer/'private_all76_regional_descriptors.tsv',sep='\t');summaries=[]
 for (model,contrast),s in d.groupby(['model','contrast']):
  v=s.partial_Spearman_rho.dropna().to_numpy();loo=[np.median(np.delete(v,i)) for i in range(len(v))] if len(v)>1 else [np.nan]
  summaries.append(dict(model=model,contrast=contrast,total_libraries=len(s),estimable_libraries=len(v),positive_libraries=int((v>0).sum()),negative_libraries=int((v<0).sum()),median=float(np.median(v)) if len(v) else np.nan,q25=float(np.quantile(v,.25)) if len(v) else np.nan,q75=float(np.quantile(v,.75)) if len(v) else np.nan,leave_one_library_median_min=float(min(loo)),leave_one_library_median_max=float(max(loo)),no_P_or_patient_CI=True))
 pd.DataFrame(summaries).to_csv(a.output/'private_all4_complete_regional_summary.tsv',sep='\t',index=False)
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':6,'axes.titlesize':7,'axes.labelsize':6,'xtick.labelsize':5,'ytick.labelsize':5,'pdf.fonttype':42,'svg.fonttype':'none'})
 figs=[]
 def save(fig,name):
  for ext in ['pdf','svg','png']:fig.savefig(a.output/(name+'.'+ext),dpi=400)
  plt.close(fig);figs.append(name)
 fig,axes=plt.subplots(1,4,figsize=(183/25.4,140/25.4));fig.subplots_adjust(left=.13,right=.99,bottom=.18,top=.85,wspace=.3)
 codes=sorted(d.library_code.unique())
 for j,(contrast,model) in enumerate((x,y) for x in ['same_spot_R_C','neighbor_R_C_given_own_R'] for y in ['quality_T_IFN_cycle_regions_xy','plus_incomplete_legacy_composition']):
  ax=axes[j];part=d[(d.contrast==contrast)&d.model.eq(model)].set_index('library_code');ax.axvline(0,color='#999',lw=.8)
  for i,c in enumerate(codes):
   r=part.loc[c]
   if np.isfinite(r.partial_Spearman_rho):ax.plot(r.partial_Spearman_rho,i,'o',ms=2.5,color='#0072B2' if j%2==0 else '#D55E00')
   else:ax.text(0,i,'NA',va='center',ha='center',color='#999',fontsize=5)
  ax.set_xlim(-1,1);ax.set_ylim(len(codes)-.5,-.5);ax.set_yticks(range(len(codes)));ax.set_yticklabels([f'Library {c+1}' for c in codes] if j==0 else ['']*len(codes));ax.set_xlabel('Partial Spearman rho');ax.set_title(('Same region' if j<2 else 'Neighbor R / own C')+'\n'+('Fixed controls' if j%2==0 else '+ partial composition'))
 fig.suptitle('Fixed CM02 functional panels: complete regional descriptors',fontsize=8)
 fig.text(.07,.06,'All19 ccRCC Visium library units, all76 fixed slots. Points are descriptive, with no P or patient confidence intervals.\nT/IFN/cycle, RNA quality, author TLS/tumor regions and xy trends are fixed controls; neighboring R also adjusts own R.\nMixed regions and incomplete legacy composition do not identify exact cells, treatment response or direct communication.',fontsize=5.5)
 save(fig,'All76_regional_conditional_RNA_descriptors')
 V=Path(json.loads((a.root/'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root']);z=np.load(V/'GSE273952_spatial_neighbourhood_preflight_20261016_v1/private_spatial_graph_and_partial_composition.npz',allow_pickle=False);s=np.load(a.producer/'private_fixed34_counts_and7_scores.npz',allow_pickle=False)['scores'];degree=np.diff(z['graph_indptr'])
 for page,start in enumerate(range(0,len(codes),5)):
  group=codes[start:start+5];fig,axes=plt.subplots(len(group),3,figsize=(183/25.4,230/25.4),squeeze=False);fig.subplots_adjust(left=.07,right=.92,top=.94,bottom=.09,wspace=.15,hspace=.45)
  for row,lib in enumerate(group):
   allrows=np.where(z['library_codes']==lib)[0];ix=np.where((z['library_codes']==lib)&~z['folded_excluded']&(degree>0))[0];x=z['pixel_x'][allrows];y=z['pixel_y'][allrows];limits=(float(np.min(x)),float(np.max(x)),float(np.min(y)),float(np.max(y)))
   for col,(values,title) in enumerate([(s[:,0],'R fixed panel'),(s[:,1],'C fixed panel'),(np.minimum(s[:,0],s[:,1]),'min(R,C)')]):
    ax=axes[row,col];ax.scatter(z['pixel_x'][allrows],z['pixel_y'][allrows],c='#dddddd',s=.3,linewidths=0,rasterized=True);im=ax.scatter(z['pixel_x'][ix],z['pixel_y'][ix],c=values[ix],s=.8,cmap='viridis',vmin=0,vmax=1,linewidths=0,rasterized=True)
    ax.set_xlim(limits[0],limits[1]);ax.set_ylim(limits[3],limits[2]);ax.set_aspect('equal');ax.set_xticks([]);ax.set_yticks([]);ax.set_title(f'Library {lib+1} / '+title)
  cax=fig.add_axes([.94,.25,.012,.5]);fig.colorbar(im,cax=cax).set_label('Mean fixed-gene percentile within library')
  fig.suptitle('Complete ccRCC tissue regions: fixed functional RNA panels',fontsize=8)
  fig.text(.06,.025,'Every spot displayed; grey includes folded/disconnected or unscorable positions. Same pixel limits within each library.\nRelative within-library gene percentiles, not absolute molecule levels or CM W/S/H. Mixed spots do not show cell contact.',fontsize=5.5)
  save(fig,f'All19_libraries_fixed_panel_maps_page{page+1}')
 page=['<!doctype html><meta charset="utf-8"><title>CM02 mixed-region RNA descriptors</title><style>body{font:16px/1.7 system-ui;max-width:1100px;margin:25px auto;padding:20px}img{width:100%}</style><h1>CM02：全部区域功能描述与19文库图</h1><p>本图册保留全部76槽与19文库。区域RNA、部分相关和原细胞CM分别解释，没有P或患者置信区间。</p>']
 for name in figs:page+=['<h2>'+html.escape(name)+'</h2><a href="'+name+'.pdf">PDF</a> · <a href="'+name+'.svg">SVG</a><img src="'+name+'.png">']
 (a.output/'index.html').write_text('\n'.join(page));assert sum(p.stat().st_size for p in a.output.iterdir() if p.is_file())<64*1024**2
 receipt=dict(status='ALL5_REGIONAL_FIGURES_ALL19_LIBRARIES_RENDERED_VISUAL_REVIEW_PENDING',finished_utc=datetime.now(timezone.utc).isoformat(),n_figures=len(figs),n_libraries_shown=len(codes),n_model_slots=76,no_new_P=True,audit_sha256=sha(a.audit/'independent_audit_public.json'),source_sha256=sha(__file__));(a.output/'render_terminal_public.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))

if __name__=='__main__':
 p=argparse.ArgumentParser()
 for name in ['root','producer','audit','output']:p.add_argument('--'+name,type=Path,required=True)
 main(p.parse_args())
