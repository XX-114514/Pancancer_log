from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,time
import numpy as np,pandas as pd
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);D=V/'GSE273952_spatial_RNA_descriptors_20261016_v1';A=V/'GSE273952_spatial_RNA_descriptors_audit_20261016_v1';sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
for t in range(480):
 if (A/'inventory.json').exists():break
 e=V/(A.name+'.execution.json')
 if e.exists() and json.loads(e.read_text()).get('exit_code') not in [None,0]:raise RuntimeError('Descriptor audit failed')
 time.sleep(30)
else:raise TimeoutError('Audit wait')
assert json.loads((A/'inventory.json').read_text())['status'].startswith('PASS_ALL');m=pd.read_csv(D/'public_all2451_programme_spatial_descriptors.tsv',sep='\t');lr=pd.read_csv(D/'public_all9_neighbour_NAMPT_target_RNA_descriptors.tsv',sep='\t');recs=json.loads((V/'all_lineage_state_candidate_dictionary_20261014_v4/frozen_programme_dictionary.json').read_text())['programmes'];ids=[r['programme_id'] for r in recs];w=m[m.programme_id.isin(ids)].pivot(index='programme_id',columns='model',values='median');O=V/'GSE273952_spatial_RNA_figures_20261016_v1';assert not O.exists();O.mkdir(mode=0o700);primary='quality_author_regions';sensitivity=['plus_incomplete_legacy6_composition','plus_spatial_polynomial_gradient'];plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
fig,axs=plt.subplots(1,2,figsize=(10.8,5.2))
for ax,key,title in zip(axs,sensitivity,['Add incomplete legacy6 composition','Add spatial polynomial gradient']):
 t=w[[primary,key]].dropna();ax.scatter(t[primary],t[key],s=13,alpha=.5,color='#176b87',edgecolors='none');lim=[min(-.05,t.min().min()-.03),max(.1,t.max().max()+.03)] if len(t) else [-.05,.1];ax.plot(lim,lim,'--',lw=1,color='grey');ax.set(xlim=lim,ylim=lim,xlabel='Primary median residual Moran I',ylabel='Sensitivity median residual Moran I',title=title+'\npaired measurable programmes: '+str(len(t)))
fig.suptitle('GSE273952: frozen programme spatial RNA structure',fontsize=13);fig.text(.5,.025,'Mixed spots; cohort medians across available libraries. Descriptors, no significance or cell-intrinsic-state claim.',ha='center',fontsize=9);fig.tight_layout(rect=(0,.05,1,.92))
for ext in ['png','pdf','svg']:fig.savefig(O/('all_frozen_programme_spatial_sensitivity.'+ext),dpi=180)
plt.close(fig);fig,ax=plt.subplots(figsize=(10.5,7));models=[primary]+sensitivity;colors=['#176b87','#d16a25','#6c4e9b'];labels={primary:'Quality + author regions',sensitivity[0]:'+ incomplete legacy6 fractions',sensitivity[1]:'+ spatial gradient'};order=['RNA_DDIT3_TARGETS_SIGNED','RNA_DDIT3_TARGETS_POS','RNA_DDIT3_TARGETS_NEG'];ticks=[];ticklabels=[]
for i,endpoint in enumerate(order):
 for k,model in enumerate(models):
  row=lr[(lr.endpoint==endpoint)&(lr.model==model)].iloc[0];yy=(2-i)*4+(2-k);ticks.append(yy);ticklabels.append(endpoint.replace('RNA_DDIT3_TARGETS_','DDIT3 ')+': '+labels[model]);label=f"{int(row.libraries_estimable)}/19; +{int(row.positive_libraries)}/-{int(row.negative_libraries)}"
  if np.isfinite(row['median']):ax.plot([row.q25,row.q75],[yy,yy],lw=3,color=colors[k]);ax.scatter([row['median']],[yy],s=38,color=colors[k],zorder=3)
  ax.text(1.015,yy,label,transform=ax.get_yaxis_transform(),va='center',fontsize=9)
ax.axvline(0,color='grey',ls='--',lw=1);ax.set(yticks=ticks,yticklabels=ticklabels,xlim=(-1,1),ylim=(-1,12),xlabel='Partial Spearman rho: median and library IQR (not confidence interval)',title='Neighbour NAMPT RNA versus fixed DDIT3 target RNA');fig.text(.5,.02,'19 specimen/library units; patient bridge unreviewed. Same-spot receptor/own NAMPT adjusted. No direct binding or cell-pair claim.',ha='center',fontsize=8.5);fig.tight_layout(rect=(0,.06,.91,.96))
for ext in ['png','pdf','svg']:fig.savefig(O/('fixed_NAMPT_target_RNA_effect_spread.'+ext),dpi=180)
plt.close(fig);out=dict(status='TWO_FROZEN_SPATIAL_DESCRIPTOR_FIGURES_RENDERED_PENDING_PIXEL_REVIEW',observed_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),descriptor_audit_sha256=sha(A/'inventory.json'),pvalues_or_CI_claim=False,files={p.name:sha(p) for p in O.iterdir() if p.is_file()});(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
