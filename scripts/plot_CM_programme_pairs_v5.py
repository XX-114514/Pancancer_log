"""Plot every BH record and all observed paired-context coverage reasons."""
from pathlib import Path
import os,json,hashlib,re
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);S=V/'original_CM_programme_pair_associations_20261015_v2';iv=json.loads((S/'inventory.json').read_text());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
for name,h in iv['output_sha256'].items():assert sha(S/name)==h
d=pd.read_csv(S/'all_programme_pair_associations.tsv',sep='\t');s=d[d.BH_q<.05].sort_values(['gse_id','context','programme_A','programme_B']).reset_index(drop=True);assert len(s)==36;c=pd.read_csv(S/'all_context_coverage_reasons.tsv',sep='\t');O=V/'CM_programme_pair_figures_20261016_v4';assert not O.exists();O.mkdir()
plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
colors={'GSE178341':'#2878B5','GSE161529':'#339966','GSE270680':'#CC7935'}
replacements={'B_cell':'B','T_cell':'T','Macrophage':'Mac','Endothelial_cell':'Endo','Fibroblast':'Fib','Mural_cell':'Mural'}
def short(x):
 x=x.removeprefix('ACN_')
 for a,b in replacements.items():x=x.replace(a,b)
 return x
fig,ax=plt.subplots(figsize=(13,12));y=np.arange(len(s))[::-1]
for i,r in s.iterrows():
 ax.hlines(y[i],r.CI_low,r.CI_high,color=colors[r.gse_id],lw=1.2,alpha=.8);ax.scatter(r.adjusted_partial_rank,y[i],s=23,color=colors[r.gse_id],zorder=3)
labels=[f'{short(r.programme_A)} / {short(r.programme_B)} | {r.gse_id} c{r.context} n={r.analysis_blocks}' for r in s.itertuples()];ax.set_yticks(y,labels,fontsize=7);ax.axvline(0,color='#777777',lw=.8,ls='--');ax.set_xlim(-1.03,1.03);ax.set_ylim(-1,len(s));ax.set_xlabel('Adjusted partial rank correlation (95% block-bootstrap interval)');ax.grid(axis='x',alpha=.16);ax.set_title('Within-identity programme associations: all 36 exploratory BH records',pad=60)
handles=[Line2D([0],[0],color=color,marker='o',lw=1,label=f'{sid} ({name})') for sid,color,name in [('GSE178341',colors['GSE178341'],'colorectal'),('GSE161529',colors['GSE161529'],'breast'),('GSE270680',colors['GSE270680'],'gastric')]];ax.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,1.04),ncol=3,fontsize=7,frameon=False)
fig.text(.015,.014,'11,716 tests; 20,000 block permutations; one BH family. Cross-K records overlap. c = context; n = analysis blocks.\nBreast contexts use sample fallback blocks. These are programme associations, with no inferred ligand-receptor direction.',fontsize=8);fig.subplots_adjust(left=.46,right=.98,top=.91,bottom=.065)
for ext in ['pdf','svg','png']:fig.savefig(O/('all36_programme_associations_forest.'+ext),dpi=180)
plt.close(fig)
counts=c.status.value_counts();order=['INSUFFICIENT_BLOCKS','CONSTANT_PRESPECIFIED_COVARIATE','ESTIMABLE'];names=['Fewer than 20 analysis blocks','Constant prespecified covariate','Estimable full design'];vals=[int(counts.get(k,0)) for k in order];fig,ax=plt.subplots(figsize=(9.5,3.8));bars=ax.barh(np.arange(3),vals,color=['#aeb7c0','#d6a75a','#2878b5']);ax.invert_yaxis();ax.set_yticks(np.arange(3),names);ax.set_xscale('log');ax.set_xlim(.8,15000);ax.set_xlabel('Observed paired source-stratum contexts (log scale)');ax.set_title('Why most fixed CM pair contexts cannot be tested')
for b,v in zip(bars,vals):ax.text(v*1.08,b.get_y()+b.get_height()/2,f'{v:,}',va='center')
fig.text(.02,.02,'4,331 contexts with at least one paired library. Contexts with no paired capture are outside this table.\n60 fixed identity pairs; 8 pairs / 3 GSEs have at least one estimable context.',fontsize=8);fig.subplots_adjust(left=.36,right=.95,top=.85,bottom=.22)
for ext in ['pdf','svg','png']:fig.savefig(O/('context_coverage_reasons.'+ext),dpi=180)
plt.close(fig);out=dict(status='PLOTS_COMPLETE_ALL36_AND_ALL4331_CONTEXTS',source_inventory_sha256=sha(S/'inventory.json'),producer_sha256=sha(__file__),files={p.name:sha(p) for p in O.iterdir() if p.is_file()});(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
