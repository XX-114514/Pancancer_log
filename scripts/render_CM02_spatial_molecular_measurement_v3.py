"""V3 overlays primary binding specificity review without changing any counts."""
import argparse,json,html
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from NSCLC_native_MFA_interpretation_inputs_v1 import sha

def main(a):
 proof=json.loads((a.audit/'independent_audit_public.json').read_text());assert proof['n_slots']==83232
 assert sha(a.producer/'terminal_public.json')==proof['producer_terminal_sha256']
 assert not a.output.exists();a.output.mkdir(mode=0o700)
 d=pd.read_csv(a.producer/'private_all83232_spatial_molecular_opportunities.tsv.gz',sep='\t')
 plan=json.loads((a.producer/'private_bound_plan_before_molecular_counts.json').read_text())
 panels=json.loads((a.root.parents[1]/'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json').read_text())['panels']
 defining={k:set(panels['pan_T']+panels['regulatory_associated']+panels['cytotoxic']) for k in ['R','C']};defining['A']=set(panels['pan_T']+['HLA-DRA','CD74'])
 d['reference_supported_qualified']=d.measurement_qualified & d.n_reference_supported.gt(0) & ~d.curated_reference_rejected
 review={('TGFB1','TGFBR1_TGFBR2'):'canonical_complete_receptor_supported_not_activity',('IFNG','IFNGR1_IFNGR2'):'canonical_complete_receptor_supported_not_activity',('CCL5','CCR5'):'canonical_functional_pair_supported_not_activity',('CCL5','CCR4'):'specificity_disputed_do_not_promote',('CD80','CTLA4'):'binding_supported_not_current_cell_source',('CD86','CTLA4'):'binding_supported_not_current_cell_source'}
 d['primary_binding_review']=[review.get((l,r),'not_individually_reviewed_this_task') for l,r in zip(d.ligand,d.receptor)]
 d['group_defining_gene_overlap']=[';'.join(sorted(set(l.split('_')+r.split('_'))&(defining[di[0]]|defining[di[-1]]))) for l,r,di in zip(d.ligand,d.receptor,d.direction)]
 d.to_csv(a.output/'private_all83232_with_group_gene_overlap.tsv.gz',sep='\t',index=False)
 # Retain the entire resource table; qualification is a measurement gate, no rank.
 q=d[d.reference_supported_qualified];q.to_csv(a.output/'private_all_reference_supported_qualified_positions.tsv.gz',sep='\t',index=False)
 summary=d.groupby(['source_index','direction'])[['subunits_measured','measurement_qualified','reference_supported_qualified']].sum().reset_index()
 summary['all_resource_rows']=4624;summary.to_csv(a.output/'private_complete_measurement_funnel.tsv',sep='\t',index=False)
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':6,'axes.titlesize':7,'axes.labelsize':6,'xtick.labelsize':5,'ytick.labelsize':5,'pdf.fonttype':42,'svg.fonttype':'none'})
 figures=[]
 def save(fig,name):
  for ext in ['pdf','svg','png']:fig.savefig(a.output/(name+'.'+ext),dpi=400)
  plt.close(fig);figures.append(name)
 ordered=[(s,k) for s in range(3) for k in ['R_to_C','C_to_R','R_to_A','A_to_R','C_to_A','A_to_C']]
 x=summary.set_index(['source_index','direction']).reindex(ordered)
 values=x[['all_resource_rows','subunits_measured','measurement_qualified','reference_supported_qualified']].to_numpy(int)
 fig,ax=plt.subplots(figsize=(183/25.4,145/25.4));fig.subplots_adjust(left=.23,right=.94,top=.87,bottom=.19)
 im=ax.imshow(np.log1p(values),cmap='cividis',aspect='auto');ax.set_xticks(range(4));ax.set_xticklabels(['All resource','All subunits\nmeasured','Expression/group\ngates passed','+legacy record\nsupport']);ax.set_yticks(range(18));ax.set_yticklabels([f'Assay{s+1} / {k.replace("_to_"," > ")}' for s,k in ordered])
 for i in range(18):
  for j in range(4):ax.text(j,i,str(values[i,j]),ha='center',va='center',fontsize=6,color='white' if values[i,j]<100 else 'black')
 ax.set_title('All 83,232 spatial molecular measurement positions');fig.colorbar(im,ax=ax,fraction=.025,pad=.025).set_label('log1p(count)')
 fig.text(.07,.06,'R/C: mutually exclusive fixed functional RNA capture in T-capture cells.\nA: non-T-capture HLA-DRA+CD74 RNA group; unavailable when either defining gene is unmeasured.\nQualification and nearby RNA opportunities are not enrichment, receptor assembly, signaling or independent patients.',fontsize=5.5)
 save(fig,'All_resource_measurement_funnel')
 coord=a.root/'GSE299886_original_coordinate_sidecars_20261018_v1';reports=json.loads((coord/'private_all_source_coordinate_reconstruction_receipt.json').read_text())
 old=a.root/'CM02_MERFISH_conditional_spatial_mixing_20261018_v1'
 fig,axes=plt.subplots(2,3,figsize=(183/25.4,150/25.4));fig.subplots_adjust(left=.07,right=.99,top=.91,bottom=.17,wspace=.25,hspace=.28)
 for i,rec in enumerate(reports):
  xy=np.load(coord/rec['output'],allow_pickle=False)['xy_source'];z=np.load(old/('private_source_%02d_graph_state_and_null.npz'%i),allow_pickle=False);unit=float(z['unit']);xy=(xy-xy.min(0))/unit
  ti=z['existing_T_row_index'];state=z['state_threshold1'];name=['Lung TMA1','Lung TMA2','Mesothelioma TMA2'][i]
  for row in range(2):
   ax=axes[row,i]
   if row==0:ax.scatter(xy[:,0],xy[:,1],s=.12,color='#d5d9df',alpha=.3,rasterized=True,linewidths=0)
   else:ax.scatter(xy[ti,0],xy[ti,1],s=.6,color='#c9ced4',rasterized=True,linewidths=0)
   for value,color,label in [(1,'#0072B2','R only'),(2,'#D55E00','C only'),(3,'#CC79A7','Both')]:
    index=ti[state==value];ax.scatter(xy[index,0],xy[index,1],s=1.2,color=color,rasterized=True,linewidths=0)
   ax.set_xlim(-.03*xy[:,0].max(),1.03*xy[:,0].max());ax.set_ylim(-.03*xy[:,1].max(),1.03*xy[:,1].max());ax.set_aspect('equal');ax.set_title(name+('\nAll-cell context' if row==0 else '\nFull T-capture positions'));ax.set_xlabel('Source x / median NN');ax.set_ylabel('Source y / median NN')
 legend=[Line2D([],[],marker='o',ls='',color=c,markersize=3,label=l) for c,l in [('#0072B2','R only'),('#D55E00','C only'),('#CC79A7','Both'),('#c9ced4','Other T capture')]]
 fig.legend(handles=legend,loc='lower center',bbox_to_anchor=(.5,.085),ncol=4,frameon=False)
 fig.text(.07,.035,'Every accepted source coordinate shown, no selected tissue region or sampled cells. Relative units; not verified micrometres.\nTechnical TMA inputs can contain overlapping people. Capture patterns are not independently validated Treg/CD8 identities.',fontsize=5.5)
 save(fig,'All3_complete_tissue_R_C_maps')
 # All qualifying molecular records remain in the delivered full table. These
 # named canonical receptor diagnostics are explicit, with no statistical claims.
 diagnostics=[('TGFB1','TGFBR1_TGFBR2'),('IFNG','IFNGR1_IFNGR2'),('CCL5','CCR4'),('CCL5','CCR5'),('CD80','CTLA4'),('CD86','CTLA4')]
 col=[(s,k) for s in range(3) for k in ['R_to_C','C_to_R']];v=np.full((len(diagnostics),len(col)),np.nan)
 diagnostic_rows=[]
 for i,(l,r) in enumerate(diagnostics):
  for j,(source,di) in enumerate(col):
   rec=d[(d.source_index==source)&d.direction.eq(di)&d.ligand.eq(l)&d.receptor.eq(r)]
   if len(rec):
    row=rec.iloc[0];diagnostic_rows.append(row.to_dict())
    if row.reference_supported_qualified:v[i,j]=row.ligand_to_complete_receptor_edges
 pd.DataFrame(diagnostic_rows).to_csv(a.output/'private_named_receptor_diagnostics_all_positions.tsv',sep='\t',index=False)
 fig,ax=plt.subplots(figsize=(183/25.4,88/25.4));fig.subplots_adjust(left=.26,right=.94,top=.85,bottom=.31)
 cmap=plt.get_cmap('cividis').copy();cmap.set_bad('#eeeeee');im=ax.imshow(np.log1p(v),cmap=cmap,aspect='auto')
 ax.set_yticks(range(len(diagnostics)));ax.set_yticklabels([l+' > '+r.replace('_','+')+(' [disputed]' if (l,r)==('CCL5','CCR4') else '') for l,r in diagnostics]);ax.set_xticks(range(6));ax.set_xticklabels([f'Assay{s+1}\n{k.replace("_to_"," > ")}' for s,k in col]);ax.set_title('Curated-entry RNA diagnostics: nearby positive pairs')
 for i in range(len(diagnostics)):
  for j in range(6):ax.text(j,i,'NA' if np.isnan(v[i,j]) else str(int(v[i,j])),ha='center',va='center',fontsize=6,color='white' if np.isfinite(v[i,j]) and np.log1p(v[i,j])<np.nanmax(np.log1p(v))/2 else 'black')
 fig.colorbar(im,ax=ax,fraction=.025,pad=.025).set_label('log1p(pair count)')
 fig.text(.06,.08,'NA retains missing measurement or failed group/positive-cell gates. No P, activity or treatment-effect claim.\nThis named display explains receptor completeness; all 4,624 entries and 83,232 positions are retained in the full table.\nCCL5-CCR4 specificity disputed: PMID9169480/9430724. CTLA4 overlaps group-defining genes.\nPair counts reflect abundance, capture and distance; they are not effect sizes, ligand activity or binding validation.',fontsize=5.5)
 save(fig,'Named_complete_receptor_diagnostics')
 page=['<!doctype html><meta charset="utf-8"><title>CM02 spatial molecular measurement</title><style>body{font:16px/1.7 system-ui;max-width:1100px;margin:25px auto;padding:20px}img{width:100%}</style><h1>CM02：空间分子测量与完整来源图</h1><p>图中RNA捕获、测量资格和邻近计数均不等于直接通讯、受体活化或临床验证。原18槽邻近检验无校正后稳定富集，保持原结论。</p>']
 for name in figures:page+=['<h2>'+html.escape(name)+'</h2><a href="'+name+'.pdf">PDF</a> · <a href="'+name+'.svg">SVG</a><img src="'+name+'.png">']
 (a.output/'index.html').write_text('\n'.join(page))
 receipt=dict(status='ALL3_AGGREGATE_FIGURES_RENDERED_VISUAL_REVIEW_PENDING',finished_utc=datetime.now(timezone.utc).isoformat(),n_figures=3,producer_terminal_sha256=sha(a.producer/'terminal_public.json'),independent_audit_sha256=sha(a.audit/'independent_audit_public.json'),no_new_inference=True,group_defining_gene_overlap_marked_all_positions=True,source_sha256=sha(__file__))
 (a.output/'render_terminal_public.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))

if __name__=='__main__':
 p=argparse.ArgumentParser()
 for name in ['root','producer','audit','output']:p.add_argument('--'+name,type=Path,required=True)
 main(p.parse_args())
