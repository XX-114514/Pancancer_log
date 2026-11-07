"""Standalone figures from independently accepted aggregate results only."""
import os,json,hashlib
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
r=Path(os.environ['V4_REVIEW_ROOT']);a=r/'completed_bulk_v1';out=r/'figures';out.mkdir(exist_ok=True)
proof=json.loads((a/'independent_measurement_audit.json').read_text());assert proof['status'].startswith('PASS_');c=pd.read_csv(a/'all2100_clinical_slots.tsv',sep='\t',float_precision='round_trip');cms=['CM%02d'%i for i in range(1,15)];cancers=json.loads((r/'configuration/plan_before_effects.json').read_text())['bulk']['TCGA_cancers'];cohorts=['IMvigor210','Riaz2017','Hugo2016','Jung2019','Rose2021','SU2C_MARK'];entries=[]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
def save(fig,name,caption):
 for ext in ['png','svg']:fig.savefig(out/(name+'.'+ext),dpi=230)
 plt.close(fig);entries.append({'name':name,'caption':caption})
def matrix(frame,rows,ax,title,unit):
 x=frame.pivot(index='cohort',columns='axis',values='log_effect').reindex(index=rows,columns=cms)/np.log(2);q=frame.pivot(index='cohort',columns='axis',values='BH_q').reindex(index=rows,columns=cms);cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#e3e5e7');im=ax.imshow(x.to_numpy(),aspect='auto',cmap=cmap,vmin=-2,vmax=2)
 for ii,jj in zip(*np.where(q.to_numpy()<.05)):ax.plot(jj,ii,'o',color='black',ms=3)
 ax.add_patch(Rectangle((1.5,-.5),1,len(rows),fill=False,lw=1.2,edgecolor='#303d52'));ax.set_xticks(range(14),cms,rotation=65,ha='right');ax.set_yticks(range(len(rows)),rows);ax.set_title(title);cb=plt.colorbar(im,ax=ax,fraction=.035,pad=.025,extend='both');cb.set_label(unit+' per within-cohort SD; clipped +/-2');return im
fig,axes=plt.subplots(1,2,figsize=(14,12));fig.subplots_adjust(left=.08,right=.95,top=.91,bottom=.14,wspace=.28)
for ax,ep in zip(axes,['PFI','OS']):matrix(c[c.family.eq('TCGA_survival')&c.endpoint.eq(ep)&c.model.eq('clinical_basic')],cancers,ax,'TCGA '+ep+': clinical basic model','log2 HR')
fig.suptitle('2026-11-06 V4: all 14 frozen CMT gene-program proxies',fontsize=15);fig.text(.08,.045,'Black dot: BH q<0.05 across all 1,848 TCGA slots. Outline: prespecified CM03. Grey: not estimable.\nAge, variable sex and pan-T adjustment; stage sensitivity remains in complete tables.\nThese programs are not exact cellular-module abundance. Single-cancer prognosis is not ICI benefit.',fontsize=10)
save(fig,'all14_TCGA_PFI_OS','All33 cancers and14 factors. Basic Cox model; BH across complete1848-slot family; unavailable retained. Displayed log2 HR per population SD, colors clipped; no cutpoint. Full stage sensitivity is in tables.')
fig,axes=plt.subplots(1,3,figsize=(16,5.5));fig.subplots_adjust(left=.09,right=.96,top=.80,bottom=.29,wspace=.39)
for ax,ep in zip(axes,['response','OS','PFS']):matrix(c[~c.family.eq('TCGA_survival')&c.endpoint.eq(ep)],cohorts,ax,ep,'log2 OR' if ep=='response' else 'log2 HR')
for ax in axes[1:]:ax.set_yticklabels(['']*len(cohorts))
fig.suptitle('All 14 CMT programs in all 6 ICI cohorts: no complete-family BH q<0.05',fontsize=13);fig.text(.09,.045,'All 84 response and 168 survival slots retained; respective full-family BH. Grey is unavailable, not a null effect.\nModels adjust for pan-T. Response OR>1 and survival HR>1 have different clinical meanings. No optimized thresholds.',fontsize=9)
save(fig,'all14_six_ICI_endpoints','Complete6 cohorts and14 factors in response/OS/PFS. No BH-positive associations in either complete family; grey denotes unavailable. Gene-program proxy, not validated prediction or treatment-specific benefit.')
fig,axes=plt.subplots(1,3,figsize=(13,5));fig.subplots_adjust(left=.14,right=.96,top=.83,bottom=.26,wspace=.36)
for ax,ep in zip(axes,['response','OS','PFS']):
 ax.axvline(1,color='#929292',lw=.8);ax.set_xscale('log');ax.set_xlim(.125,8);ax.set_ylim(5.7,-.7);ax.set_yticks(range(6),cohorts if ep=='response' else ['']*6);ax.set_title(ep+(' (OR)' if ep=='response' else ' (HR)'));ax.set_xlabel('Effect per cohort SD (95% Wald CI)')
 for yy,co in enumerate(cohorts):
  row=c[c.axis.eq('CM03')&c.cohort.eq(co)&c.endpoint.eq(ep)]
  assert len(row)==1;z=row.iloc[0]
  if z.status=='exploratory_model_numeric_audit_pass':
   lo,hi,mu=z.CI_low,z.CI_high,z.effect;ax.plot([max(lo,.125),min(hi,8)],[yy,yy],color='#173d63',lw=1.4);ax.plot(np.clip(mu,.125,8),yy,'o',mfc='#173d63' if z.BH_q<.05 else 'white',mec='#173d63',ms=6)
   if lo<.125:ax.plot(.125,yy,'<',color='#173d63',clip_on=False)
   if hi>8:ax.plot(8,yy,'>',color='#173d63',clip_on=False)
  else:ax.text(1,yy,'NA',va='center',ha='center',color='#888888',bbox={'facecolor':'white','edgecolor':'none'})
fig.suptitle('Prespecified CM03-associated 10-gene program: ICI associations',fontsize=14);fig.text(.14,.055,'No estimable CM03 ICI endpoint reaches complete-family BH q<0.05; all 95% Wald intervals cross 1.\nNA: fixed measurement, sample-size or endpoint gate. This is not exact B-endothelial-NK module activity.',fontsize=9)
save(fig,'CM03_six_ICI_forest','All6 cohorts and3 prespecified ICI endpoints for CM03 program. Pan-T adjusted; displayed Wald CI and complete-family BH, unavailable retained. No ICI support for this frozen program; not evidence that every possible module assay has no clinical association.')
manifest={'status':'AGGREGATE_FIGURES_RENDERED_PENDING_VISUAL_QA','UTC':datetime.now(timezone.utc).isoformat(),'figures':entries,'files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file()},'no_patient_ids_or_rows':True};(r/'figure_manifest_v1.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps({'n_figures':len(entries),'output':str(out)}))
