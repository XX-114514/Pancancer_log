"""Plot all prespecified target terms from reviewed aggregate statistics only."""
from pathlib import Path
import os,json,hashlib
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'svg.fonttype':'none'})
r=Path(os.environ['BREAST_AGGREGATE_ROOT']);t=pd.read_csv(r/'all24_prespecified_target_terms.tsv',sep='\t');out=r/'figures';out.mkdir(exist_ok=True)
keys=[('primary21',x) for x in ['CM03','CM04','CM07','CM09','CM11','Vascular_B_NK_equal_mean','Vascular_B_NK_min_component']]+[('secondary_incremental3','Vascular_B_NK_min_component')]
labels=['CM03 program *','CM04 program','CM07 program','CM09 program *','CM11 program','B + EC + NK mean','B / EC / NK minimum','Minimum, adjusted for components']
fig,axes=plt.subplots(1,3,figsize=(13.8,5.8),sharey=True)
for col,(scope,title) in enumerate([('pembro','Pembrolizumab + chemotherapy (n=69)'),('control','Chemotherapy control (n=179)'),('interaction','Treatment interaction (n=248)')]):
 ax=axes[col];ax.axvline(1,color='#666',lw=.9,ls='--');ax.axhline(.5,color='#bdbdbd',lw=.8)
 for i,(family,score) in enumerate(keys):
  z=t[(t.family==family)&(t.score==score)&(t.target==scope)].iloc[0];y=7-i
  if z.status!='estimable':ax.text(1,y,'Unavailable',ha='center',va='center',color='#888',fontsize=8);continue
  ax.errorbar(z.OR,y,xerr=[[z.OR-z.CI_low],[z.CI_high-z.OR]],fmt='o',ms=4,color='#176d79' if family=='primary21' else '#aa5d22',capsize=2,lw=1.1)
  ax.text(1,y-.27,f'q={z.BH:.3f}',ha='center',va='center',fontsize=8,color='#555',bbox=dict(facecolor='white',edgecolor='none',pad=.4))
 ax.set_xscale('log');ax.set_xlim(.07,6);ax.set_xticks([.1,.25,.5,1,2,4]);ax.set_xticklabels(['0.1','0.25','0.5','1','2','4']);ax.set_title(title,fontsize=10,pad=13);ax.set_xlabel('Odds ratio per pooled-cohort SD (95% Wald CI)' if scope!='interaction' else 'Ratio of odds ratios (95% Wald CI)',fontsize=9);ax.grid(axis='x',alpha=.14);ax.spines[['top','right']].set_visible(False)
axes[0].set_yticks(np.arange(8));axes[0].set_yticklabels(labels[::-1],fontsize=9);axes[0].set_ylim(-.6,7.6)
fig.suptitle('I-SPY2 breast cancer: complete prespecified RNA-proxy analysis',fontsize=14,y=.985)
fig.text(.025,.035,'* Full CM03 / CM09 programs are unmeasurable on this array. Pan-T control in this separate analysis uses measured CD3D/CD3E.\nPrimary 21 and secondary 3 target slots have separate BH families, including unavailable slots. No target passes BH q<0.05.\nRNA proxies do not measure cell fractions or spatial organization. Control concurrency is not certified from GEO fields alone.',fontsize=8.5,va='bottom')
fig.subplots_adjust(left=.255,right=.99,bottom=.24,top=.85,wspace=.14)
for ext in ['png','svg']:
 p=out/('ISPY2_all24_target_terms.'+ext);fig.savefig(p,dpi=200)
 if ext=='svg':p.write_text('\n'.join(line.rstrip() for line in p.read_text().splitlines())+'\n')
plt.close(fig)
(r/'figure_manifest.json').write_text(json.dumps({'scope':'all24 slots from aggregate tables; no individual data','files':[{'name':str(p.relative_to(r)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size} for p in sorted(out.iterdir())]},indent=2)+'\n')
