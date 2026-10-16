"""Show every BH LR record, including missing intervals, without mechanism labels."""
from pathlib import Path
import os,json,hashlib
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);D=V/'state_after_DIALOGUE_LR_response_20261016_v1';p=D/'all_estimable_exploratory_LR_receiver_associations.tsv';sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();iv=json.loads((D/'inventory.json').read_text());assert sha(p)==iv['files'][p.name];t=pd.read_csv(p,sep='\t');t=t[t.BH_q<.05].sort_values(['context','source','target','MCP','ligand','receptor']).reset_index(drop=True);assert len(t)==14;O=V/'state_DLG_LR_review_figures_20261016_v1';assert not O.exists();O.mkdir();names=[];plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'svg.fonttype':'none'})
short=lambda s:s.replace('_cell','').replace('Blood_vascular_endothelial','Endothelial').replace('Regulatory_T','Treg').replace('Macrophage','Mac').replace('Epithelial','Epi')
fig,(ax,bx)=plt.subplots(1,2,figsize=(16,9),sharey=True,gridspec_kw={'width_ratios':[3,1]});y=np.arange(len(t))
for k,r in t.iterrows():
 ok=np.isfinite(r.CI_low) and np.isfinite(r.CI_high);c='#167d8d' if ok else '#956f49';ax.plot([r.CI_low,r.CI_high],[k,k],color=c,lw=2) if ok else None;ax.scatter(r.adjusted_masked_rho,k,edgecolors=c,marker='o' if ok else 'D',s=40,facecolors=c if ok else 'none');bx.barh(k,r.blocks,color=c,alpha=.65);bx.text(r.blocks+.6,k,str(r.blocks),va='center',fontsize=8);names.append(f'C{r.context} {r.MCP} | {short(r.source)} > {short(r.target)} | {r.ligand} - {r.receptor}')
ax.set_yticks(y,names);ax.invert_yaxis();ax.axvline(0,color='.5',lw=.8);ax.set_xlim(-1.03,1.03);ax.set_xlabel('Adjusted partial-rank association with LR-gene-masked receiver MCP');ax.grid(axis='x',alpha=.2);bx.set_xlabel('Observed analysis blocks');bx.set_xlim(0,t.blocks.max()*1.2);bx.tick_params(axis='y',left=False,labelleft=False)
fig.suptitle('October 2026 CM (v20261010): all 14 BH LR records',fontsize=15);fig.text(.02,.025,'Teal: qualifying 95% block-bootstrap interval. Hollow brown: interval unavailable (<1900/2000 valid resamples).\nSame-cohort exploratory associations; analysis blocks include reviewed patients and sample fallbacks. Not independent axes or causal validation.',fontsize=9);fig.subplots_adjust(left=.46,right=.98,top=.91,bottom=.14,wspace=.1)
for ext in ['png','pdf','svg']:fig.savefig(O/('all14_BH_LR_forest_and_coverage.'+ext),dpi=160)
plt.close(fig);t.to_csv(O/'all14_plotted_records.tsv',sep='\t',index=False);out=dict(status='ALL14_BH_RECORDS_PLOTTED_WITH_MISSING_INTERVALS_EXPLICIT',producer_sha256=sha(__file__),primary_inventory_sha256=sha(D/'inventory.json'),records=14,qualifying_intervals=int(t.CI_low.notna().sum()),files={x.name:sha(x) for x in O.iterdir() if x.is_file()});(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
