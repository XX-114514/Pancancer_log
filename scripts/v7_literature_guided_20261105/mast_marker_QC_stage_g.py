"""Marker-to-target annotated plots and sample-aware Mast QC review.
Creates review flags and a published-QC sensitivity pool, never final exclusions.
"""
import os,json,hashlib,math
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd,scipy.sparse as sp,anndata as ad,zarr,h5py
from anndata.io import read_elem
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
ROOT=Path(os.environ.get('PROJECT_ROOT',str(globals().get('PROJECT_ROOT',''))));assert ROOT.is_absolute()
r=ROOT/'annotation_rebuild/20261105_v7_literature_guided_reannotation_v1';o=r/'stage_g_mast_marker_QC_v2';o.mkdir(exist_ok=False)
started=datetime.now(timezone.utc).isoformat();spec=json.loads((r/'config/mast_marker_panels_v2.json').read_text())
a=ad.read_h5ad(r/'stage_b_mast_fullgene/mast_v7_fullgene_raw.h5ad');X=a.X.tocsr();md=pd.read_csv(r/'stage_c2_mast_initial_diagnostic/mast_initial_metadata_PRIVATE.tsv.gz',sep='\t',index_col='cell_id');assert md.index.equals(a.obs_names)
cl=md['RNA_snn_res.0.4'].astype(int).to_numpy();obs=a.obs.copy();n=a.n_obs;assert n==45036
panels=spec['panels'];refs=spec['panMC_reference_panels'];allgenes=list(dict.fromkeys(g for panel in panels+refs for g in panel['genes']));assert set(allgenes)<=set(a.var_names)
j=a.var_names.get_indexer(allgenes);raw=X[:,j].toarray();tot=np.asarray(X.sum(axis=1)).ravel();ng=np.asarray((X>0).sum(axis=1)).ravel();cp=raw*(10000/np.maximum(tot,1)[:,None]);where={g:i for i,g in enumerate(allgenes)}
gses=sorted(obs.gse_id.astype(str).unique());av=[]
for gs in gses:
 p=ROOT/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects/harmonized_per_gse'/(gs+'_symbol_counts.h5ad')
 with h5py.File(p,'r') as f:genes=pd.Index(read_elem(f['var']).index.astype(str))
 av.append(a.var_names.isin(genes))
av=np.asarray(av).T;gidx=pd.Index(gses).get_indexer(obs.gse_id.astype(str));measured=av[j][:,gidx].T
H=sp.csr_matrix((np.ones(n,np.int32),(cl,np.arange(n))),shape=(27,n));S=sp.csr_matrix((np.ones(n,np.int32),(cl,gidx)),shape=(27,len(gses)))
den=np.asarray(S@av[j].T.astype(np.int32));nn=np.asarray(H@(raw>0).astype(np.int32));ss=np.asarray(H@cp);groupn=np.bincount(cl,minlength=27)
expression=[]
for ci in range(27):
 for gj,g in enumerate(allgenes):
  d=int(den[ci,gj]);expression.append({'cluster':ci,'gene':g,'n_cells':int(groupn[ci]),'n_measured':d,'coverage_fraction':float(d/groupn[ci]),'n_expressing':int(nn[ci,gj]),'pct_expressing_measured':float(nn[ci,gj]/d) if d else np.nan,'mean_cp10k_measured':float(ss[ci,gj]/d) if d else np.nan})
e=pd.DataFrame(expression);e.to_csv(o/'mast27_expanded_measured_markers.tsv',sep='\t',index=False,na_rep='NA')
# Actual original-sample/method calibration, using all cells in each sample, not only Mast.
core=zarr.open_group(str(ROOT/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects/pancancer_core_fullgene_v1.zarr'),mode='r');pos=obs.global_position.to_numpy().astype(int)
assert pd.Index(read_elem(core['obs/_index']))[pos].equals(obs.index)
cal=pd.DataFrame({key:read_elem(core['obs/'+key]) for key in ['gse_id','sampleID','doublet_method','doublet_score']});cal['doublet_score']=pd.to_numeric(cal.doublet_score,errors='coerce')
for key in ['gse_id','sampleID']:assert np.array_equal(cal[key].astype(str).to_numpy()[pos],obs[key].astype(str).to_numpy())
grp=cal.groupby(['gse_id','sampleID','doublet_method'],dropna=False,observed=True)['doublet_score'];rank=grp.rank(pct=True,method='average');size=grp.transform('count')
obs['doublet_within_sample_method_percentile']=rank.to_numpy()[pos];obs['doublet_calibration_n']=size.to_numpy()[pos];obs['doublet_calibration_evaluable']=(obs.doublet_calibration_n>=100)&obs.doublet_within_sample_method_percentile.notna()
obs['raw_library_total']=tot;obs['raw_library_nFeature']=ng
program={}
for panel in panels:
 ids=[where[g] for g in panel['genes']];name=panel['id'];program[name]={}
 for k in [1,2,3]:
  hit=(raw[:,ids]>=k).sum(axis=1);obs[name+'_n_genes_ge'+str(k)+'UMI']=hit;program[name][k]=hit>=2
 obs[name+'_all_genes_measured']=measured[:,ids].all(axis=1)
 obs[name+'_marker_UMI']=raw[:,ids].sum(axis=1)
parent=program['mast_identity'][1];strong_parent=program['mast_identity'][2]
competitors=['T_competitor','myeloid_competitor','stromal_competitor','epithelial_competitor']
anycomp=np.zeros(n,bool)
for c in competitors:anycomp|=program[c][2]
obs['strong_Mast_and_competitor']=strong_parent&anycomp
for pct in [.9,.95,.99]:obs['joint_cross_program_doublet_q'+str(int(pct*100))]=strong_parent&anycomp&obs.doublet_calibration_evaluable.to_numpy()&(obs.doublet_within_sample_method_percentile.to_numpy()>=pct)
joint=obs['joint_cross_program_doublet_q95'].to_numpy();ambiguous=strong_parent&anycomp&~joint
obs['review_route_v2']=np.select([~parent,joint,ambiguous],['Mast_parent_low_marker_review','joint_doublet_candidate_review','cross_program_without_joint_doublet_review'],default='Mast_parent_supported_subtype_open')
# This is the exact published numeric QC sensitivity, not a universal final-QC gate.
mtgenes=['MT-ATP6','MT-ATP8','MT-CO1','MT-CO2','MT-CO3','MT-CYB','MT-ND1','MT-ND2','MT-ND3','MT-ND4','MT-ND4L','MT-ND5','MT-ND6'];mtidx=a.var_names.get_indexer(mtgenes);assert (mtidx>=0).all()
mtcov=av[mtidx].sum(axis=0)[gidx];mt=pd.to_numeric(obs.pct_counts_mt,errors='coerce').to_numpy();qknown=np.isfinite(mt)&(mtcov>=10)
paperqc=qknown&(tot>=600)&(tot<=120000)&(ng>=400)&(ng<=8000)&(mt<=10)
obs['panMC_published_QC_evaluable']=qknown;obs['panMC_published_QC_pass']=paperqc
obs['Mast_source_robust_refit_candidate']=parent&paperqc&~joint
obs['mast_review_cluster_res0.4']=cl
# Non-exclusive positive marker-associated candidate states. No CMA1-negative call.
obs['MCTC_marker_candidate_raw1']=parent&program['MCTC'][1];obs['MCTC_marker_candidate_raw2']=strong_parent&program['MCTC'][2]
obs['MCT_marker_candidate_raw1']=parent&program['MCT'][1];obs['MCT_marker_candidate_raw2']=strong_parent&program['MCT'][2]
obs['cycling_marker_candidate']=parent&program['cycling'][1]
obs.to_pickle(o/'mast_cell_QC_marker_review_PRIVATE.pkl.gz')
# Source robustness is an explicit modeling sensitivity; full RNA/markers stay intact.
coverage=av.mean(axis=1);features=pd.DataFrame({'gene_symbol':a.var_names,'n_Mast_sources_measured':av.sum(axis=1),'fraction_Mast_sources_measured':coverage,'source_robust_model_candidate':(coverage>=.8)&~a.var_names.str.match('^HB[AB]')})
features.to_csv(o/'source_feature_coverage_PRIVATE.tsv',sep='\t',index=False)
features.loc[features.source_robust_model_candidate,['gene_symbol']].to_csv(o/'source_robust_model_genes.tsv',sep='\t',index=False)
rows=[]
for ci in range(27):
 t=obs.iloc[np.flatnonzero(cl==ci)];z={'cluster':ci,'n_cells':len(t),'n_sources':t.gse_id.nunique(),'n_original_sample_units':t.groupby(['gse_id','sampleID'],observed=True).ngroups}
 for c in ['doublet_calibration_evaluable','strong_Mast_and_competitor','joint_cross_program_doublet_q90','joint_cross_program_doublet_q95','joint_cross_program_doublet_q99','panMC_published_QC_evaluable','panMC_published_QC_pass','Mast_source_robust_refit_candidate','MCTC_marker_candidate_raw1','MCTC_marker_candidate_raw2','MCT_marker_candidate_raw1','MCT_marker_candidate_raw2','cycling_marker_candidate']:z['n_'+c]=int(t[c].sum())
 for c in competitors:z['n_strong_Mast_AND_'+c]=int((strong_parent&(program[c][2]))[cl==ci].sum())
 for c in ['doublet_within_sample_method_percentile','decontx_contamination','pct_counts_mt']:z[c+'_median']=float(pd.to_numeric(t[c],errors='coerce').median())
 rows.append(z)
summary=pd.DataFrame(rows);summary.to_csv(o/'mast27_sample_aware_QC_review.tsv',sep='\t',index=False,na_rep='NA')
counts=obs.review_route_v2.value_counts().rename_axis('review_route_v2').reset_index(name='n_cells');assert counts.n_cells.sum()==n;counts.to_csv(o/'mast_review_routes_v2.tsv',sep='\t',index=False)
# Study/sample identity stays private; only aggregate candidate reproduction is public.
rep=[]
for flag in ['MCTC_marker_candidate_raw1','MCTC_marker_candidate_raw2','MCT_marker_candidate_raw1','MCT_marker_candidate_raw2','cycling_marker_candidate','joint_cross_program_doublet_q95']:
 t=obs[obs[flag]];rep.append({'candidate':flag,'n_cells':len(t),'n_sources':t.gse_id.nunique(),'n_original_sample_units':t.groupby(['gse_id','sampleID'],observed=True).ngroups})
pd.DataFrame(rep).to_csv(o/'positive_program_source_replication.tsv',sep='\t',index=False)
review_labels={ci:'Mast (subtype open)' for ci in range(27)}
review_labels.update({7:'Mast + epithelial program [review]',10:'Mast + T program [review]',13:'Mast + stromal/epithelial [review]',15:'Mast + myeloid program [review]',16:'Cycling-enriched Mast [state]',17:'Mast; other program / QC [review]',25:'Mast + keratin program [review]',26:'Mast + prostate/epithelial [review]'})
pd.DataFrame([{'cluster':k,'display_annotation':v,'annotation_status':'working_review_not_final'} for k,v in review_labels.items()]).to_csv(o/'mast27_display_annotation_v2.tsv',sep='\t',index=False)
maprows=[]
for panel in panels+refs:
 for gene in panel['genes']:maprows.append({'panel_id':panel['id'],'target_annotation':panel['target'].replace('\n',' '),'level':panel.get('level','published_reference_state'),'gene':gene,'source':panel.get('source',spec['panMC_source']),'transfer_status':'review association only; not a transferred label'})
pd.DataFrame(maprows).to_csv(o/'marker_to_type_state_map.tsv',sep='\t',index=False)
def plot(groups,filename,title):
 slots=[];bounds=[];x=0
 for panel in groups:
  lo=x
  for gene in panel['genes']:slots.append((x,gene));x+=1
  bounds.append((lo,x-1,panel));x+=.8
 fig=plt.figure(figsize=(max(22,x*.43+7),13));ax=fig.add_axes([.23,.14,.65,.67]);palette={'parent':'#c6dbef','candidate_subtype':'#dadaeb','state_axis':'#fee8c8','state':'#fee8c8','state_program':'#fee8c8','competing_parent':'#fdd0c0','published_reference_state':'#d9ead3'}
 gene_stats={g:(np.log1p(e.loc[e.gene.eq(g),'mean_cp10k_measured']).mean(),np.log1p(e.loc[e.gene.eq(g),'mean_cp10k_measured']).std(ddof=0)) for _,g in slots}
 xs=[];ys=[];sizes=[];colors=[];missing=[]
 for xx,g in slots:
  mu,sd=gene_stats[g]
  for z in e[e.gene.eq(g)].to_dict('records'):
   if z['coverage_fraction']<.8 or not np.isfinite(z['pct_expressing_measured']):missing.append((xx,z['cluster']));continue
   xs.append(xx);ys.append(z['cluster']);sizes.append(140*z['pct_expressing_measured']);colors.append(np.clip((np.log1p(z['mean_cp10k_measured'])-mu)/(sd or 1),-2,2))
 sc=ax.scatter(xs,ys,s=sizes,c=colors,cmap='RdBu_r',vmin=-2,vmax=2,edgecolors='none')
 if missing:mx,my=zip(*missing);ax.scatter(mx,my,s=14,c='#999999',marker='x')
 ax.set_xticks([p[0] for p in slots],[p[1] for p in slots],rotation=65,ha='right',fontsize=10)
 ax.set_yticks(range(27),[f'V7 C{k:02d} | {review_labels[k]}' for k in range(27)],fontsize=10);ax.set_ylim(26.8,-.8);ax.set_xlim(-.7,x-.9)
 for txt,ci in zip(ax.get_yticklabels(),range(27)):
  if '[review]' in review_labels[ci]:txt.set_color('#a23b32')
 ax.set_ylabel('Working identity / state and review status',fontsize=12);ax.grid(axis='y',alpha=.13);ax.set_axisbelow(True)
 top=ax.secondary_xaxis('top');top.set_xticks([(lo+hi)/2 for lo,hi,p in bounds],[p['target'] for lo,hi,p in bounds],fontsize=10);top.tick_params(length=0,pad=16)
 for lo,hi,panel in bounds:
  color=palette[panel.get('level','published_reference_state')];ax.add_patch(Rectangle((lo-.42,-.075),hi-lo+.84,.022,transform=ax.get_xaxis_transform(),color=color,clip_on=False))
  ax.plot([lo-.42,hi+.42],[1.008,1.008],transform=ax.get_xaxis_transform(),color='#666666',lw=1,clip_on=False)
 cbax=fig.add_axes([.90,.29,.012,.36]);cb=fig.colorbar(sc,cax=cbax);cb.set_label('Gene-scaled log1p(mean CP10K)',fontsize=11)
 handles=[ax.scatter([],[],s=140*f,color='#777777',label=str(int(f*100))+'% detected') for f in [.1,.5,1]];ax.legend(handles=handles,loc='upper left',bbox_to_anchor=(1.025,1.05),frameon=False,fontsize=10)
 fig.text(.23,.955,title,fontsize=19,weight='bold');fig.text(.23,.91,'Column groups name the corresponding cell type, candidate subtype or state; rows show working annotations.',fontsize=11)
 fig.text(.23,.028,'Dot size = detection among source-measured cells; x = source coverage <80%. Blue is relative expression, not absence.\nCandidate subtype/state panels are literature associations, not accepted labels. Ref C0-C9 are author state IDs, not the V7 cluster IDs.',fontsize=10)
 fig.savefig(o/(filename+'.png'),dpi=170);fig.savefig(o/(filename+'.pdf'));plt.close(fig)
plot(panels,'Mast_marker_targets_v2','V7 Mast | identity, candidate subtypes, states and competing lineages')
plot([panels[0]]+refs+[next(p for p in panels if p['id']=='T_competitor')],'Mast_panMC_reference_compass_v2','V7 Mast | comparison with the published pan-MC state marker panels')
manifest={'status':'completed_marker_targets_and_sample_aware_review','UTC_start':started,'UTC_end':datetime.now(timezone.utc).isoformat(),'n_cells':n,'n_sources':len(gses),'n_marker_genes':len(allgenes),'n_marker_cluster_records':len(e),'n_original_sample_method_groups_all_core':int(grp.ngroups),'n_evaluable_doublet_calibration':int(obs.doublet_calibration_evaluable.sum()),'n_joint_q95_doublet_candidates':int(joint.sum()),'n_Mast_refit_candidates_published_QC':int(obs.Mast_source_robust_refit_candidate.sum()),'n_source_robust_model_genes':int(features.source_robust_model_candidate.sum()),'review_route_counts':dict(zip(counts.review_route_v2,counts.n_cells.astype(int))),'new_final_annotation_changes':0,'limits':['Reference panels are not transferred labels.','Doublet percentiles are within original sample and method; joint flags remain candidates.','Published numeric QC is a sensitivity pool, not universal final QC.','No corrected full-gene decontX expression was assumed.','Source coverage filtering is for model features only; full library RNA and marker axes retained.']}
manifest['outputs']=[{'name':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in o.iterdir() if p.is_file()]
(o/'complete.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(manifest))
