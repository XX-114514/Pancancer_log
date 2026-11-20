#!/usr/bin/env python3
"""Render a multi-release evidence story from externally reviewed aggregates.

No embeddings, clinical models, module fits or inferential tests are refitted.
All numerical research inputs are supplied outside the public code repository.
"""
import argparse,hashlib,json,shutil
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd

COL=['#567F98','#147D92','#8F849D','#B8664D','#C89362']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,
 'axes.titleweight':'bold','axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,
 'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,
 'svg.fonttype':'none','figure.facecolor':'white','savefig.facecolor':'white'})

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def panel(ax,letter,title):
 ax.set_title(title,loc='left',pad=12);ax.text(-.09,1.08,letter,transform=ax.transAxes,fontweight='bold',fontsize=13,va='top')
def blank(ax):
 ax.set_xticks([]);ax.set_yticks([])
 for s in ax.spines.values():s.set_visible(False)
def box(ax,x,y,w,h,title,body,color='#EDF3F7',size=8):
 ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.012',fc=color,ec='#8798A7',lw=.8))
 ax.text(x+w/2,y+h-.05,title,ha='center',va='top',fontweight='bold',fontsize=9)
 ax.text(x+w/2,y+h/2-.02,body,ha='center',va='center',fontsize=size,linespacing=1.45)

def render(root,prior,output,configuration=None):
 output.mkdir(exist_ok=False,parents=True);sources=output/'source_tables';sources.mkdir()
 ev=root/'evidence';j=json.loads((ev/'integrated_aggregate_evidence_v2.json').read_text())
 cfg=json.loads((configuration or ev/'figure_configuration_PRIVATE.json').read_text());records=[]
 def save(fig,name):
  for ext in ['png','pdf','svg']:fig.savefig(output/(name+'.'+ext),dpi=330,bbox_inches='tight',pad_inches=.14)
  plt.close(fig);records.append({'id':name,'files':{e:sha(output/(name+'.'+e)) for e in ['png','pdf','svg']}});print('RENDERED',name,flush=True)
 def table(name,x):x.to_csv(sources/(name+'.tsv'),sep='\t',index=False)
 # 1. Distinct resources, with accession overlap rather than summed cells.
 fig=plt.figure(figsize=(9.1,9.0));gs=fig.add_gridspec(2,2,height_ratios=[1.1,1],hspace=.55,wspace=.5)
 ax=fig.add_subplot(gs[0,:]);blank(ax);ax.set_xlim(0,1);ax.set_ylim(0,1);panel(ax,'a','Linked resources with distinct measurement scopes')
 for i,row in enumerate(cfg['resources']):
  box(ax,.015+i*.333,.56,.29,.31,row['name'],row['display'],size=8)
 for i,(title,body) in enumerate(cfg['branches']):
  box(ax,.015+i*.333,.08,.29,.30,title,body,color='#F7F7F7',size=7.6)
  ax.annotate('',xy=(.16+i*.333,.4),xytext=(.16+i*.333,.55),arrowprops={'arrowstyle':'->','color':'#8798A7'})
 ov=pd.read_csv(ev/'source_accession_overlap.tsv',sep='\t');patterns=ov.groupby(['Legacy_TME','Immune_v1','R4']).size().sort_values(ascending=False)
 ax=fig.add_subplot(gs[1,0]);labs=[' + '.join(n for n,flag in zip(['Legacy','Immune','R4'],k) if flag) for k in patterns.index]
 ax.barh(range(len(patterns)),patterns.values,color=COL[0]);ax.set_yticks(range(len(patterns)),labs);ax.invert_yaxis();ax.set_xlabel('Accession / collection labels');panel(ax,'b','Source overlap is explicit')
 for i,n in enumerate(patterns.values):ax.text(n+.5,i,str(n),va='center',fontsize=8)
 ax.set_xlim(0,max(patterns.values)*1.2)
 ax=fig.add_subplot(gs[1,1]);blank(ax);panel(ax,'c','Separate questions for each evidence layer')
 roles=[('Composition','Which identities covary?'),('RNA programs','Which states vary together?'),('Spatial organization','Which states are locally enriched?'),('Clinical readouts','Which fixed scores associate with outcomes?')]
 for i,(title,body) in enumerate(roles):
  ax.text(0,.92-i*.23,title,fontweight='bold',fontsize=10);ax.text(0,.84-i*.23,body,fontsize=8.6)
 fig.text(.09,.01,'Shared accessions do not establish exact cell or patient overlap. Cell totals are never added across releases.',fontsize=8)
 save(fig,'Figure1_multilayer_resources');table('source_overlap',ov)
 # 3. Legacy broad-lineage composition, source-specific evidence and robustness.
 fig=plt.figure(figsize=(9.2,8.6));gs=fig.add_gridspec(2,2,hspace=.62,wspace=.58)
 members=pd.read_csv(ev/'legacy_public/fig4_cm_top10_membership.tsv',sep='\t');members=members[members.top10.astype(str).str.lower().eq('true')]
 ct=pd.crosstab(members.cm,members.majorCluster);ct=ct.loc[:,ct.sum()>0]
 ax=fig.add_subplot(gs[0,0]);im=ax.imshow(ct,aspect='auto',cmap='Blues',vmin=0,vmax=max(1,ct.to_numpy().max()));ax.set_yticks(range(len(ct)),['L-'+x for x in ct.index]);ax.set_xticks(range(len(ct.columns)),ct.columns,rotation=70,ha='right',fontsize=7);fig.colorbar(im,ax=ax,fraction=.025,pad=.025,label='Top-ten node count');panel(ax,'a','Legacy TME modules span multiple lineages')
 cross=pd.read_csv(ev/'legacy172_composition_cross_source.tsv',sep='\t');counts=[len(cross),int((cross.ci_low>0).sum()),int((cross.ci_high<0).sum())]
 ax=fig.add_subplot(gs[0,1]);ax.barh(range(3),counts,color=[COL[0],COL[1],COL[3]]);ax.set_yticks(range(3),['All eligible groups','Positive interval','Negative interval']);ax.invert_yaxis();ax.set_xlabel('Pair x cancer x material groups');panel(ax,'b','Within-cancer, multiple-source composition')
 for i,n in enumerate(counts):ax.text(n+.4,i,str(n),va='center',fontsize=8)
 ax.set_xlim(0,max(counts)*1.2);ax.text(0,-.35,'Legacy 172-node release; intervals are descriptive,\nnot multiplicity-controlled replication claims.',transform=ax.transAxes,fontsize=7)
 bal=pd.DataFrame(j['balance']['component_comparisons']);bal=bal[bal['mode']=='equal_source_weight'];ax=fig.add_subplot(gs[1,0]);modules=sorted(bal.reference_module.unique())
 for i,m in enumerate(modules):
  q=bal[bal.reference_module==m];ax.scatter(np.full(len(q),i)+np.linspace(-.08,.08,len(q)),q.effective_A_cosine,color=COL[1],s=23);ax.plot([i-.18,i+.18],[q.effective_A_cosine.median()]*2,color='#222222',lw=1)
 ax.set_xticks(range(len(modules)),['I-T '+x for x in modules],rotation=35,ha='right');ax.set_ylabel('Matched effective-basis cosine');ax.set_ylim(.7,1.02);panel(ax,'c','Immune T modules under equal source weight');ax.text(0,-.32,'Four seeds share the same discovery data.\nThe smallest source receives greater relative weight.',transform=ax.transAxes,fontsize=7)
 ax=fig.add_subplot(gs[1,1]);m=j['R4']['composition']['masked_cosine_medians'];xs=list(m);ax.scatter(range(len(xs)),[m[x] for x in xs],color=COL[2],s=24);ax.set_xticks(range(len(xs)),[x.replace('CM','') for x in xs]);ax.set_ylim(0,1.04);ax.set_ylabel('Matched-cosine median');ax.set_xlabel('R4 CM');panel(ax,'d','R4 response to unavailable parent populations');ax.text(0,-.32,'This denominator-sensitivity design differs from\nthe source-weight sensitivity in panel c.',transform=ax.transAxes,fontsize=7)
 save(fig,'Figure3_composition_landscapes');table('immune_source_balance_all32',pd.DataFrame(j['balance']['component_comparisons']));table('legacy172_cross_source_composition',cross);table('legacy_top10_parent_counts',ct.rename_axis('cm').reset_index())
 # 4. The stable immune core and its clinical proxy, with complete context.
 clinical=pd.read_csv(ev/'immune_clinical_all_model_slots.tsv',sep='\t');fig=plt.figure(figsize=(9.4,8.4));gs=fig.add_gridspec(2,2,hspace=.62,wspace=.7)
 ax=fig.add_subplot(gs[0,0]);blank(ax);panel(ax,'a','A frozen regulatory / cytotoxic T-cell theme')
 ax.text(0,.93,'I-T CM02: seven retained core members',fontweight='bold',fontsize=9)
 for i,t in enumerate(cfg['immune_core']):ax.text(.03,.80-i*.09,t,fontsize=8.4,color=COL[1])
 ax.text(0,.04,'Clinical proxy: min(regulatory score, cytotoxic score)\nwith pan-T adjustment; not a W/H projection.',fontsize=8)
 ax=fig.add_subplot(gs[0,1]);panel(ax,'b','Fixed-proxy associations across ICI contexts')
 selected=[]
 for cohort,family,endpoint in cfg['clinical_focus']:
  q=clinical[(clinical.cohort==cohort)&(clinical.family==family)&(clinical.endpoint==endpoint)&(clinical.axis=='pCM02_regulatory_cytotoxic')];assert len(q)==1;selected.append(q.iloc[0])
 for i,r in enumerate(selected):
  ax.plot([r.CI_low,r.CI_high],[i,i],color=COL[1],lw=2);ax.scatter(r.effect,i,color=COL[1],s=27);ax.text(.03,i-.18,f'n={int(r.n)}, events={int(r.events)}, BH={r.BH_q:.3g}',transform=ax.get_yaxis_transform(),fontsize=7)
 ax.axvline(1,color='#AAB3BB',ls='--');ax.set_xscale('log');ax.set_yticks(range(len(selected)),['IMvigor response OR','SU2C OS HR']);ax.invert_yaxis();ax.set_ylim(1.6,-.6);ax.set_xlim(.2,10);ax.set_xlabel('Effect per proxy SD (95% interval)')
 ax=fig.add_subplot(gs[1,0]);cv=pd.DataFrame(j['CV_contrasts']);cv=cv[(cv.module_proxy=='CM02_proxy')&cv.candidate.isin(['joint_panT','joint_arms_panT'])].dropna(subset=['pooled_paired_logloss_improvement']);assert len(cv)==4
 for i,r in cv.reset_index(drop=True).iterrows():ax.plot([r.conditional_CI_low,r.conditional_CI_high],[i,i],color=COL[2],lw=1.6);ax.scatter(r.pooled_paired_logloss_improvement,i,color=COL[2],s=23)
 labels=[r.cohort.replace('_MARK','')+' | '+('vs pan-T' if r.candidate=='joint_panT' else 'vs both arms') for r in cv.itertuples()];ax.set_yticks(range(4),labels,fontsize=7);ax.invert_yaxis();ax.axvline(0,color='#AAB3BB',ls='--');ax.set_xlabel('Fixed-OOF log-loss improvement (conditional 95% interval)');panel(ax,'c','Association does not establish predictive increment')
 ax=fig.add_subplot(gs[1,1]);family=clinical.groupby('family').agg(planned=('status','size'),estimable=('p',lambda z:z.notna().sum()),BH_positive=('BH_q',lambda z:(z<.05).sum())).reset_index();family=family.set_index('family').reindex(cfg['clinical_family_order']).reset_index()
 y=np.arange(len(family));ax.barh(y,family.estimable,color=COL[0],label='Estimable');ax.barh(y,family.planned-family.estimable,left=family.estimable,color='#D9DDE1',label='Unavailable');ax.set_yticks(y,[s.replace('_',' ') for s in family.family],fontsize=7);ax.invert_yaxis();ax.set_xlabel('Planned model slots');ax.legend(frameon=False,fontsize=7,loc='lower right');panel(ax,'d','All frozen clinical families remain visible')
 fig.text(.07,-.005,'Each clinical family retains its original correction. OOF intervals condition on saved predictions; all four CM02 intervals include zero.',fontsize=7.5)
 save(fig,'Figure4_immune_core_and_clinical_proxy');table('immune_clinical_all312',clinical);table('immune_clinical_family_coverage',family);table('immune_predictive_increment_all24',pd.DataFrame(j['CV_contrasts']))
 # 5. Explicitly distinct spatial sources, preserving absent measurements.
 fig=plt.figure(figsize=(9.2,8));gs=fig.add_gridspec(2,2,hspace=.62,wspace=.7);ns=pd.DataFrame(j['NSCLC_LN_primary']);td=pd.DataFrame(j['melanoma_TDLN_primary']);ends=cfg['T_spatial_endpoints'];names=cfg['T_spatial_names']
 ax=fig.add_subplot(gs[0,0]);blank(ax);panel(ax,'a','A specific tissue question for T-cell states')
 ax.text(0,.92,'NSCLC regional lymph nodes',fontweight='bold');ax.text(0,.80,'Curative-intent resection | relative IMC states',fontsize=8)
 box(ax,.03,.45,.88,.27,'CD4-like anchor','FOXP3 / CTLA4 relative-high candidate',size=8)
 for i,name in enumerate(names):ax.text(.05,.36-i*.11,'CD8-like: '+name,fontsize=8)
 ax.text(0,-.06,'No CD3 in this panel. Relative state labels\ndo not certify Treg identity or functional exhaustion.',fontsize=7.6)
 ax=fig.add_subplot(gs[0,1]);panel(ax,'b','Within-NSCLC, MKI67-stratified association')
 for i,e in enumerate(ends):
  r=ns[(ns.endpoint==e)&(ns.estimand=='combined_MH')].iloc[0];ax.plot(np.exp([r.bootstrap_lower,r.bootstrap_upper]),[i,i],color=COL[1],lw=2);ax.scatter(np.exp(r.median_effect),i,color=COL[1],s=28)
 ax.axvline(1,color='#AAB3BB',ls='--');ax.set_yticks(range(3),names);ax.invert_yaxis();ax.set_xlabel('Patient-summary spatial OR (95% interval)');ax.text(0,-.27,cfg['T_spatial_coverage_note'],transform=ax.transAxes,fontsize=7.5)
 ax=fig.add_subplot(gs[1,0]);panel(ax,'c','Local conditional-null comparisons')
 for i,e in enumerate(ends):
  for k,(df,col,off) in enumerate([(ns,COL[1],-.13),(td,COL[2],.13)]):
   r=df[(df.endpoint==e)&(df.estimand=='local_composition_degree_and_cycle_conditional_null')].iloc[0]
   if pd.notna(r.median_effect):ax.plot(np.exp([r.bootstrap_lower,r.bootstrap_upper]),[i+off]*2,color=col,lw=1.6);ax.scatter(np.exp(r.median_effect),i+off,color=col,s=22,label=['NSCLC LN','Melanoma TDLN'][k] if i==0 else None)
   else:ax.text(.86,i+off,'Not measurable',color='#777777',fontsize=7)
 ax.axvline(1,color='#AAB3BB',ls='--');ax.set_yticks(range(3),names);ax.invert_yaxis();ax.set_xlabel('Null-centered enrichment ratio (95% interval)');ax.legend(frameon=False,fontsize=7,loc='upper left');ax.set_ylim(2.6,-.8)
 ax=fig.add_subplot(gs[1,1]);panel(ax,'d','Complete spatial families and qualification')
 keys=list(zip(ns.endpoint,ns.estimand));state=np.zeros((len(keys),2))
 for i,(e,t) in enumerate(keys):
  for c,df in enumerate([ns,td]):
   r=df[(df.endpoint==e)&(df.estimand==t)].iloc[0];state[i,c]=0 if pd.isna(r.BH) else 3 if r.BY<.05 else 2 if r.BH<.05 else 1
 from matplotlib.colors import ListedColormap
 ax.imshow(state,aspect='auto',vmin=0,vmax=3,cmap=ListedColormap(['#D9DDE1','#F0F2F4','#9EBECD','#147D92']));ax.set_xticks([0,1],['NSCLC','Melanoma']);ax.set_yticks(range(len(keys)),[e.replace('E0_','').replace('E1_CD8_','').replace('E2_','')+' | '+{'noncycling':'MKI67-low','cycling':'MKI67-high','combined_MH':'MH','cycling_minus_noncycling':'difference','local_composition_degree_and_cycle_conditional_null':'local null'}[t] for e,t in keys],fontsize=5.7)
 ax.text(0,-.18,'Dark: BY < .05; blue: BH only; light: tested;\ngray: structurally unavailable. No cross-source difference test.',transform=ax.transAxes,fontsize=7)
 save(fig,'Figure5_T_state_spatial_context');table('T_state_NSCLC_all15',ns);table('T_state_melanoma_all15',td)
 # 8. Historical LR recurrence and adjusted communication are separate estimands.
 fig=plt.figure(figsize=(9.5,8.6));gs=fig.add_gridspec(2,2,hspace=.65,wspace=.85)
 old=pd.read_csv(ev/'legacy_public/fig5_focus_axes_replication.tsv',sep='\t');old=old[old.axis.isin(cfg['legacy_LR_axes'])].copy();assert len(old)==len(cfg['legacy_LR_axes'])
 ax=fig.add_subplot(gs[0,0]);panel(ax,'a','Legacy within-CM LR detection coverage');y=np.arange(len(old));ax.barh(y,old.n_samples_standard_0_05 if 'n_samples_standard_0_05' in old else old['n_samples_standard_0.05'],color=COL[0]);ax.set_yticks(y,[x.replace('→','->').split('  ')[0] for x in old.label],fontsize=8);ax.invert_yaxis();ax.set_xlabel('Samples meeting historical dual-rank threshold')
 for i,r in old.reset_index(drop=True).iterrows():ax.text(float(r['n_samples_standard_0.05'])+.5,i,f"{int(r['n_samples_standard_0.05'])} / {int(r.n_eligible_samples)}",va='center',fontsize=7)
 ax.set_xlim(0,float(old['n_samples_standard_0.05'].max())*1.4);ax.text(0,-.28,'Coverage is conditional on the frozen within-CM filter.\nIt is not an abundance-adjusted association or mechanism.',transform=ax.transAxes,fontsize=7)
 ax=fig.add_subplot(gs[0,1]);blank(ax);panel(ax,'b','Legacy 172-node sample-level follow-up')
 for i,(title,body) in enumerate(cfg['legacy172_communication_notes']):ax.text(0,.92-i*.23,title,fontweight='bold',fontsize=10);ax.text(0,.83-i*.23,body,fontsize=8)
 overlap=pd.read_csv(prior/'figures/v4/source_tables/Figure6_CPDB_shared_events.tsv',sep='\t');q=overlap[overlap.core_module=='CM12'].copy();ax=fig.add_subplot(gs[1,0]);y=np.arange(len(q));ax.barh(y,q.focal_shared_with_other,color=COL[0],label='Shared');ax.barh(y,q.focal_only_significant,left=q.focal_shared_with_other,color=COL[4],label='Focal-only');ax.set_yticks(y,[r.source_identity.replace('Endo_','').replace('Str_','').replace('Mye_','').replace('_',' ')+' -> '+r.target_identity.replace('Endo_','').replace('Str_','').replace('Mye_','').replace('_',' ') for r in q.itertuples()],fontsize=6);ax.invert_yaxis();ax.set_xlabel('Joint-positive CPDB event combinations');ax.legend(frameon=False,fontsize=7);panel(ax,'c','R4 CM12 detections are mostly shared')
 pat=pd.read_csv(prior/'figures/v4/source_tables/Figure6_patient_LR_coverage.tsv',sep='\t');q=pat.groupby('source',as_index=False).tested_variable_scores.sum();ax=fig.add_subplot(gs[1,1]);ax.barh(range(len(q)),q.tested_variable_scores,color=COL[2]);ax.set_yticks(range(len(q)),q.source,fontsize=8);ax.invert_yaxis();ax.set_xlabel('Patient-qualified variable LR hypotheses');panel(ax,'d','R4 patient-qualified differential tests');ax.text(0,-.28,'No test passed the complete joint family.\nThe focal EC-PC pair qualified only in the renal source.',transform=ax.transAxes,fontsize=7)
 save(fig,'Figure8_communication_across_designs');table('legacy_fixed_LR_detection',old);table('R4_shared_focal_CPDB',overlap);table('R4_patient_LR_complete_coverage',pat)
 # Additional complete source-context forest; no selected cancer or sign.
 allpairs=pd.read_csv(ev/'legacy172_composition_by_stratum.tsv',sep='\t')[['pair_id','source','target']].drop_duplicates();allpairs['pair_number']=range(1,len(allpairs)+1);pairnum=allpairs.set_index('pair_id').pair_number.to_dict();table('legacy172_pair_dictionary',allpairs)
 fig,ax=plt.subplots(figsize=(9.4,10.3));q=cross.sort_values(['cancerType','pair_id','sample_class']).reset_index(drop=True)
 for i,r in q.iterrows():ax.plot([r.ci_low,r.ci_high],[i,i],color=COL[1] if r.ci_low>0 else COL[3] if r.ci_high<0 else '#9FA9B2',lw=1);ax.scatter(r.rho,i,color=COL[1] if r.ci_low>0 else COL[3] if r.ci_high<0 else '#9FA9B2',s=11)
 ax.axvline(0,color='#BBBBBB',ls='--');ax.set_yticks(range(len(q)),[str(r.cancerType)+' | P'+str(pairnum[r.pair_id]).zfill(2)+' | '+str(r.n_GSE)+' GSE | '+{'primary_or_unspecified_tumor':'Tumor','metastatic_tumor':'Metastatic','normal_or_adjacent':'Normal'}.get(r.sample_class,r.sample_class) for r in q.itertuples()],fontsize=6);ax.invert_yaxis();ax.set_xlabel('Cross-source composition correlation (descriptive interval)');ax.set_xlim(-1,1);ax.set_title('Legacy 172-node release | all within-cancer source groups',loc='left');fig.text(.13,.015,'All eligible groups retained. Intervals are not multiplicity-controlled discovery claims; shared pairs and cohorts remain dependent.',fontsize=7)
 save(fig,'FigureS14_legacy_source_contexts')
 (output/'figure_manifest_PRIVATE.json').write_text(json.dumps({'status':'RENDERED_PENDING_VISUAL_REVIEW','figures':records,'code_sha256':sha(Path(__file__)),'input_sha256':j['input_sha256'],'new_inferential_tests':False},indent=2)+'\n')

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--prior',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--configuration',type=Path);a=p.parse_args();render(a.root,a.prior,a.output,a.configuration)
