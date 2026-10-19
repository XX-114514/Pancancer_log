"""Complete source-level internal loss diagnostics, no individual predictions."""
from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

def main():
 root=Path(os.environ['PANCANCER_RAM_OUTPUT']);out=root/'private_source_figures_v2';assert not out.exists();out.mkdir(mode=0o700)
 cases=json.loads((root/'private_all12_internal_CV_summary.json').read_text());rows=json.loads((root/'private_all24_performance_contrast_positions.json').read_text());assert len(cases)==12 and len(rows)==24
 fig,axes=plt.subplots(1,2,figsize=(14.2,9.6),sharey=True);pairs=[('joint_panT','panT_only'),('joint_arms_panT','arms_panT')]
 colors={'CM02_proxy':'#236c94','CM04_proxy':'#cf7a20'};valid=[r for r in rows if 'pooled_paired_logloss_improvement' in r]
 lo=min(r['conditional_CI_low'] for r in valid);hi=max(r['conditional_CI_high'] for r in valid);margin=(hi-lo)*.09
 labels=[c['cohort']+' / '+c['module_proxy'].replace('_proxy','')+'  (n='+str(c['n_units'])+')' for c in cases]
 for ax,pair,title in zip(axes,pairs,['Module proxy + pan-T versus pan-T','Balance + arms + pan-T versus arms + pan-T']):
  for i,c in enumerate(cases):
   r=next(r for r in rows if r['cohort']==c['cohort'] and r['module_proxy']==c['module_proxy'] and r['candidate']==pair[0])
   if 'pooled_paired_logloss_improvement' not in r:
    ax.text(lo,i,'Not estimable: N / minority gate',va='center',fontsize=8,color='0.5');continue
   x=r['pooled_paired_logloss_improvement'];a=r['conditional_CI_low'];b=r['conditional_CI_high']
   ax.errorbar(x,i,xerr=np.array([[x-a],[b-x]]),fmt='o',color=colors[c['module_proxy']],markersize=5,linewidth=1.3,capsize=3)
  ax.axvline(0,color='0.35',ls='--',lw=.8);ax.set_xlim(lo-margin,hi+margin);ax.set_title(title,fontsize=11)
  ax.grid(axis='x',color='0.9');ax.set_xlabel('Paired OOF logloss reduction; positive = better');ax.spines[['top','right']].set_visible(False)
  for k in [1.5,3.5,5.5,7.5,9.5]:ax.axhline(k,color='0.94',lw=.7)
 axes[0].set_yticks(range(12),labels,fontsize=9);axes[0].set_ylim(11.6,-.7)
 fig.subplots_adjust(left=.22,right=.985,bottom=.19,top=.855,wspace=.13)
 fig.suptitle('2026-10-14 immune atlas CM proxies (v1, T_whole_TME_primary_tumor)\nInternal prediction diagnostics with train-only gene CDFs: all 12 contexts / 24 contrasts',fontsize=12)
 fig.legend(handles=[Line2D([0],[0],marker='o',color=v,label=k.replace('_proxy',''),linestyle='none') for k,v in colors.items()],loc='upper center',bbox_to_anchor=(.5,.93),ncol=2,frameon=False)
 fig.text(.025,.035,'n = response-known units. 5-fold x 10 repeats; repeats are not independent patients.\nIntervals: pointwise paired bootstrap of fixed mean OOF predictions, conditional on fitted CV models. No new P/FDR or external validation.',fontsize=9)
 fig.savefig(out/'private_all12_context24_contrast_internal_loss.png',dpi=220);fig.savefig(out/'private_all12_context24_contrast_internal_loss.svg');plt.close(fig)
 files={p.name:dict(bytes=p.stat().st_size,SHA256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in out.iterdir()};assert sum(v['bytes'] for v in files.values())<2*1024**2
 result=dict(status='SOURCE_ONLY_COMPLETE_CV_FIGURE_RENDERED_VISUAL_REVIEW_PENDING',completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),registered_contexts=12,registered_contrasts=24,files=files,no_individual_predictions_gene_weights_or_matrices=True)
 (out/'private_figure_terminal.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
