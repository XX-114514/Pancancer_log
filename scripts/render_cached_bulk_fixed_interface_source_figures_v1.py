"""Source-only complete coverage and model summaries; no individual data."""
from pathlib import Path
import os,json,hashlib
from datetime import datetime,timezone
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm
from matplotlib.patches import Patch

def main():
 root=Path(os.environ['PANCANCER_RAM_OUTPUT']);out=root/'private_source_figures_v1';assert not out.exists();out.mkdir(mode=0o700)
 a=json.loads((root/'private_complete84_source_clinical_results.json').read_text());d=json.loads((root/'private_complete39_source_RNA_coordination_descriptions.json').read_text());assert len(a)==84 and len(d)==39
 fits=[r for r in a if r['P'] is not None];assert len(fits)>0
 plt.rcParams.update({'font.size':9,'svg.fonttype':'none'})
 fig,ax=plt.subplots(figsize=(9,max(5,len(fits)*.3+1.8)),layout='constrained')
 for i,r in enumerate(fits):
  y=len(fits)-1-i
  ax.errorbar(r['effect'],y,xerr=np.array([[r['effect']-r['CI_low']],[r['CI_high']-r['effect']]]),fmt='o',markersize=4,color='#246b8e',capsize=2,linewidth=1.1)
 labels=[r['cohort']+' / '+r['endpoint']+(' (OR)' if r['kind']=='logistic' else ' (HR)')+'   n='+str(r['n_units']) for r in fits]
 ax.set_yticks(range(len(fits)),labels[::-1]);ax.axvline(1,color='0.3',lw=.8,ls='--');ax.set_xscale('log')
 ax.set_xlim(min(r['CI_low'] for r in fits)*.8,max(r['CI_high'] for r in fits)*1.25);ax.set_xlabel('Effect per within-model SD; model-based 95% CI')
 ax.set_title('Whole-tissue fixed-interface RNA proxy: all estimable clinical models\nFull clinical family BH/BY correction retains all 84 planned positions',fontsize=11)
 ax.grid(axis='x',color='0.9');ax.spines[['top','right']].set_visible(False)
 fig.savefig(out/'private_all_estimable_clinical_forest.png',dpi=180);fig.savefig(out/'private_all_estimable_clinical_forest.svg');plt.close(fig)
 sources=[r['cohort'] for r in d];keys={r['cohort']:i for i,r in enumerate(d)}
 categories={'ESTIMABLE':1,'FIXED_N_OR_EVENTS_PER_PREDICTOR_GATE':2,'REQUIRED_BACKGROUND_GENE_MISSING':3,'INTERFACE_NOT_IN_CURRENT_SERVER_PANEL':4,'REQUIRED_STAGE_COVERAGE_LT70_PERCENT':5,'CONSTANT_REQUIRED_PREDICTOR_NO_TERMS_DROPPED':6,'NO_DECLARED_PRIMARY_TUMOUR_UNITS':7}
 cols=['response','OS','progression'];M=np.zeros((39,3),int);planned=set()
 for r in a:
  j=0 if r['endpoint']=='response' else 1 if r['endpoint']=='OS' else 2
  i=keys[r['cohort']];assert (i,j) not in planned;planned.add((i,j))
  M[i,j]=1 if r['P'] is not None else categories[r['reason']]
 assert len(planned)==84 and int((M==1).sum())==len(fits)
 colors=['#fafafa','#2166ac','#c7c7c7','#fdae61','#d73027','#dfc27d','#8c6bb1','#525252']
 fig,(ax,bx)=plt.subplots(1,2,figsize=(11,15),gridspec_kw={'width_ratios':[1,2]},layout='constrained')
 ax.imshow(M,cmap=ListedColormap(colors),norm=BoundaryNorm(np.arange(-.5,8.5),8),aspect='auto',interpolation='none')
 ax.set_yticks(range(39),sources,fontsize=8);ax.set_xticks(range(3),['Response','OS','PFS / PFI'],rotation=25);ax.set_title('All 84 clinical positions')
 for i,j in planned:ax.text(j,i,'+' if M[i,j]==1 else 'x',ha='center',va='center',color='white' if M[i,j] in [1,4,6,7] else 'black',fontsize=8)
 bx.axvline(0,color='0.5',lw=.8,ls='--');bx.set_xlim(-.5,.5);bx.set_ylim(38.5,-.5);bx.set_yticks(range(39),['']*39);bx.grid(axis='x',color='0.9')
 for i,r in enumerate(d):
  if r.get('partial_rank_effect') is not None:
   bx.plot(r['partial_rank_effect'],i,'o',color='#2166ac' if i<6 else '#4d4d4d',markersize=4)
   bx.text(.51,i,'n='+str(r['n_units']),va='center',fontsize=7,clip_on=False)
  else:bx.text(-.48,i,'Not estimable',fontsize=7,va='center',color='0.5')
 bx.set_xlabel('Partial-rank point estimate; no P value or CI');bx.set_title('All 39 source RNA-coordination descriptions')
 legend=['Not planned','Estimable','N / event gate','Required background absent','Interface outside cached panel','Stage coverage','Constant required predictor','No primary-tumour units']
 fig.legend(handles=[Patch(facecolor=c,label=l) for c,l in zip(colors,legend)],loc='outside lower center',ncol=2,fontsize=8)
 fig.suptitle('Fixed whole-tissue proxy: coverage and descriptive RNA coordination\nFive fixed RNA backgrounds; source-specific associations, not cell communication validation',fontsize=11)
 fig.savefig(out/'private_complete_family_coverage_and_coordination.png',dpi=160);fig.savefig(out/'private_complete_family_coverage_and_coordination.svg');plt.close(fig)
 files={p.name:dict(bytes=p.stat().st_size,SHA256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in out.iterdir()}
 assert sum(v['bytes'] for v in files.values())<3*1024**2
 terminal=dict(status='SOURCE_ONLY_FIGURES_RENDERED_ACTUAL_VISUAL_REVIEW_PENDING',completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),complete_clinical_positions=84,complete_source_descriptions=39,all_estimable_forest_rows=len(fits),files=files,no_patient_rows_or_gene_weights=True)
 (out/'private_figure_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(terminal))
if __name__=='__main__':main()
