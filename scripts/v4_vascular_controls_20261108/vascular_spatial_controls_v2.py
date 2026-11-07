"""Prospective vascular controls: patient summaries, regional ES, conditional LR edges.
Input cells must already have native coordinates, QC and reviewed broad identity.
All cell/region/FOv axes and patient-level rows remain private.
"""
from pathlib import Path
from itertools import product
import json,math
import numpy as np,pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import rankdata
AXES=[('CCL19','CCR7'),('CCL21','CCR7'),('CXCL12','CXCR4')]
SENDERS=['Mural','Fibroblast','Blood_EC'];RECEIVERS=['B','T','NK'];RADII=[25,50,100,200]
def fdr(p,by=False):
 p=np.asarray(p,float);o=np.argsort(p,kind='stable');rr=np.arange(1,len(p)+1);q=p[o]*len(p)/rr
 if by:q*=np.sum(1/rr)
 result=np.empty(len(p));result[o]=np.minimum(1,np.minimum.accumulate(q[::-1])[::-1]);return result
def paired_p_greater(d):
 d=np.asarray(d,float);d=d[np.isfinite(d)&(np.abs(d)>1e-12)]
 if not len(d):return 1.
 if len(d)>20:raise ValueError('Exact patient sign enumeration is bounded to20; use a separately reviewed method above20')
 ranks=rankdata(np.abs(d),method='average');observed=ranks[d>0].sum();values=np.array([0.])
 for v in ranks:values=np.r_[values,values+v]
 return float(np.mean(values>=observed-1e-12))
def neighborhood(c,genes,counts):
 rows=[];dist=np.full(len(c),np.nan);kind=np.full(len(c),-1,int);blood=c.cell_type.eq('Blood_EC').to_numpy();ack=counts[:,genes.index('ACKR1')] if 'ACKR1' in genes else None
 if ack is None:return pd.DataFrame([{'radius_um':r,'receiver':t,'family':'neighborhood_primary2' if r==50 and t in ['B','NK'] else 'neighborhood_sensitivity10','status':'unavailable_ACKR1_in_input_matrix','log_odds_enrichment':np.nan} for r in RADII for t in RECEIVERS])
 for field in c.fov.unique():
  ii=np.flatnonzero(c.fov.eq(field).to_numpy());ec=ii[blood[ii]]
  if not len(ec):continue
  d,j=cKDTree(c.loc[ec,['x_um','y_um']]).query(c.loc[ii,['x_um','y_um']]);dist[ii]=d;kind[ii]=(ack[ec[j]]>=1).astype(int)
 nonendo=~c['is_endothelial_reference'].to_numpy(bool) if 'is_endothelial_reference' in c else ~c.cell_type.isin(['Blood_EC','Lymphatic_EC']).to_numpy();lymph=c.cell_type.isin(RECEIVERS).to_numpy();n_target=int((blood&(ack>=1)).sum());n_other=int((blood&(ack<1)).sum())
 for radius in RADII:
  g1=nonendo&(kind==1)&(dist<=radius);g0=nonendo&(kind==0)&(dist<=radius);n1=int(g1.sum());n0=int(g0.sum());l1=int((g1&lymph).sum());l0=int((g0&lymph).sum())
  for typ in RECEIVERS:
   mask=c.cell_type.eq(typ).to_numpy();a=int((g1&mask).sum());b=int((g0&mask).sum());eligible=n_target>=20 and n_other>=20 and n1>=100 and n0>=100 and mask.sum()>=20
   log_or=math.log((a+.5)/(n1-a+.5))-math.log((b+.5)/(n0-b+.5)) if n1 and n0 else np.nan
   rows.append({'radius_um':radius,'receiver':typ,'family':'neighborhood_primary2' if radius==50 and typ in ['B','NK'] else 'neighborhood_sensitivity10','status':'eligible_descriptive_effect' if eligible else 'unavailable_count_or_identity_gate','n_target_EC':n_target,'n_other_EC':n_other,'n_target_nonEC':n1,'n_other_nonEC':n0,'target_receiver':a,'other_receiver':b,'target_fraction':a/n1 if n1 else np.nan,'other_fraction':b/n0 if n0 else np.nan,'log_odds_enrichment':log_or if eligible else np.nan,'target_total_lymphocytes':l1,'other_total_lymphocytes':l0,'target_conditional_lymph_fraction':a/l1 if l1 else np.nan,'other_conditional_lymph_fraction':b/l0 if l0 else np.nan,'distance_definition':'closest endothelial centroid within sameFOV; target=ACKR1 detected, other=undetected; not contact or vessel-wall distance'})
 return pd.DataFrame(rows)
def regional_ES(c,genes,counts):
 keys=pd.MultiIndex.from_arrays([c.fov,np.floor(c.x_um/100).astype(int),np.floor(c.y_um/100).astype(int)]);bin_id,levels=pd.factorize(keys);n=len(levels);nc=np.bincount(bin_id,minlength=n);rna=np.bincount(bin_id,weights=c.native_total,minlength=n)
 xy=np.column_stack([np.bincount(bin_id,weights=c[a],minlength=n)/nc for a in ['x_um','y_um']]);field=np.array([v[0] for v in levels]);blood=c.cell_type.eq('Blood_EC').to_numpy();near=np.zeros(n,bool)
 for f in np.unique(field):
  jj=np.flatnonzero(field==f);ec=np.flatnonzero(blood&c.fov.eq(f).to_numpy())
  if len(ec):near[jj]=cKDTree(c.loc[ec,['x_um','y_um']]).query(xy[jj])[0]<=100
 typecount={t:np.bincount(bin_id,weights=c.cell_type.eq(t).to_numpy(),minlength=n) for t in SENDERS+RECEIVERS};strata=pd.factorize(pd.MultiIndex.from_arrays([np.floor(np.log2(nc)),np.floor(np.log2(np.maximum(rna,1)))]))[0]
 gsum={g:np.bincount(bin_id,weights=counts[:,j],minlength=n) for j,g in enumerate(genes) if g in set(sum(([a,b] for a,b in AXES),[]))};rows=[]
 for sender,receiver,(lig,rec) in product(SENDERS,RECEIVERS,AXES):
  row={'sender':sender,'receiver':receiver,'ligand':lig,'receptor':rec,'status':'unavailable','bin_um':100,'vascular_band_um':100,'gene_detection_UMI':1,'two_independent_scRNA_CellChat_support':'pending_not_satisfied_by_two_algorithms_on_same_cells'}
  if lig not in gsum or rec not in gsum:row['reason']='LR_unavailable_in_input_matrix';rows.append(row);continue
  aa=typecount[sender];bb=typecount[receiver];groups=[near&(nc>=10)&(aa>=2)&(bb>=2),near&(nc>=10)&(aa>=2)&(bb==0),near&(nc>=10)&(aa==0)&(bb>=2)];positive=(gsum[lig]>=1)&(gsum[rec]>=1)
  for name,g in zip(['AB','A_only','B_only'],groups):row[name+'_bins']=int(g.sum());row[name+'_positive_bins']=int((g&positive).sum());row[name+'_ES_raw']=float(positive[g].mean()) if g.any() else np.nan
  common=[]
  for st in np.unique(strata):
   ix=[g&(strata==st) for g in groups];ns=[int(v.sum()) for v in ix]
   if min(ns)>=1:common.append((min(ns),[float(positive[v].mean()) for v in ix]))
  matched=[sum(int((g&(strata==st)).sum()) for st in np.unique(strata) if all((h&(strata==st)).any() for h in groups)) for g in groups]
  row['minimum_common_depth_matched_bins']=min(matched)
  if min(g.sum() for g in groups)<5 or min(matched)<5 or not common:row['reason']='insufficient_three_group_or_common_depth_support';rows.append(row);continue
  ww=np.array([v[0] for v in common]);es=np.average(np.array([v[1] for v in common]),axis=0,weights=ww)
  row.update(status='eligible_descriptive_effect',reason='',AB_ES_standardized=float(es[0]),A_only_ES_standardized=float(es[1]),B_only_ES_standardized=float(es[2]),AB_minus_A=float(es[0]-es[1]),AB_minus_B=float(es[0]-es[2]),standardization='common log2(totalUMI) x log2(cellcount) bins weighted by minimum support across3 groups; residual confounding remains')
  rows.append(row)
 t=pd.DataFrame(rows)
 for col in ['AB_minus_A','AB_minus_B']:
  if col not in t:t[col]=np.nan
 return t
def edge_expression_null(c,genes,counts,nboot=999,seed=2026110801):
 rng=np.random.default_rng(seed);rows=[];geneix={g:i for i,g in enumerate(genes)};ct=c.cell_type.to_numpy();xy=c[['x_um','y_um']].to_numpy();strata=pd.factorize(pd.MultiIndex.from_arrays([c.fov,np.floor(c.x_um/250).astype(int),np.floor(c.y_um/250).astype(int),np.floor(np.log2(np.maximum(c.native_total,1)))]))[0]
 for sender,receiver in product(SENDERS,RECEIVERS):
  si=np.flatnonzero(ct==sender);ri=np.flatnonzero(ct==receiver);ee=[];ff=[]
  for field in c.fov.unique():
   ss=np.flatnonzero(c.iloc[si].fov.eq(field).to_numpy());rr=np.flatnonzero(c.iloc[ri].fov.eq(field).to_numpy())
   if not len(ss) or not len(rr):continue
   lists=cKDTree(xy[ri[rr]]).query_ball_point(xy[si[ss]],50);nn=np.array([len(v) for v in lists]);
   if nn.sum():ee.append(np.repeat(ss,nn));ff.append(rr[np.concatenate([np.asarray(v,int) for v in lists if len(v)])])
  left=np.concatenate(ee) if ee else np.array([],int);right=np.concatenate(ff) if ff else np.array([],int)
  sg=[np.flatnonzero(strata[si]==st) for st in np.unique(strata[si])];rg=[np.flatnonzero(strata[ri]==st) for st in np.unique(strata[ri])]
  for ligand,receptor in AXES:
   row={'sender':sender,'receiver':receiver,'ligand':ligand,'receptor':receptor,'status':'unavailable','n_edges_50um':len(left),'nboot':nboot,'null':'expression labels permuted within each identity, sameFOV,250um tile and log2 nativeUMI stratum; geometry/types fixed','two_independent_scRNA_CellChat_support':'pending'}
   if ligand not in geneix or receptor not in geneix:row['reason']='LR_unavailable_in_input_matrix';rows.append(row);continue
   if len(si)<10 or len(ri)<10 or not len(left):row['reason']='identity_or_edge_support';rows.append(row);continue
   lp=counts[si,geneix[ligand]]>=1;rp=counts[ri,geneix[receptor]]>=1;observed=int(np.count_nonzero(lp[left]&rp[right]));null=np.zeros(nboot,int)
   variable_s=[ix for ix in sg if len(ix)>1 and lp[ix].any() and not lp[ix].all()];variable_r=[ix for ix in rg if len(ix)>1 and rp[ix].any() and not rp[ix].all()]
   if lp.any() and rp.any():
    for b in range(nboot):
     ll=lp.copy();rr=rp.copy()
     for ix in variable_s:ll[ix]=rng.permutation(lp[ix])
     for ix in variable_r:rr[ix]=rng.permutation(rp[ix])
     null[b]=np.count_nonzero(ll[left]&rr[right])
   # Even when observed=0, simulate unless a partner is absent; null expectation may be positive.
   p=(1+int((null>=observed).sum()))/(nboot+1)
   row.update(status='eligible_conditional_spatial_null',reason='',observed_LR_positive_edges=observed,observed_edge_fraction=observed/len(left),null_mean_edge_fraction=float(null.mean()/len(left)),excess_edge_fraction=float((observed-null.mean())/len(left)),P_MC_add1=p,positive_sender_cells=int(lp.sum()),positive_receiver_cells=int(rp.sum()),null_is_identically_zero=not(lp.any() and rp.any()))
   rows.append(row)
 t=pd.DataFrame(rows)
 if 'excess_edge_fraction' not in t:t['excess_edge_fraction']=np.nan
 if 'P_MC_add1' not in t:t['P_MC_add1']=np.nan
 p=t.get('P_MC_add1',pd.Series(np.nan,index=t.index));q=fdr(p.fillna(1));t['BH_within_case27']=np.where(p.notna(),q,np.nan);return t

def patient_summary(case_dirs,donor_keys,out):
 """Requires an explicit unique patient bridge; never treats FOVs/cells aspatients."""
 out=Path(out);out.mkdir(parents=True,exist_ok=False);assert len(case_dirs)==len(donor_keys) and len(set(donor_keys))==len(donor_keys)
 assert len({str(Path(d).resolve()) for d in case_dirs})==len(case_dirs), 'Repeated case directory cannot create extra patients'
 source_hashes=[json.loads((Path(d)/'complete.json').read_text()).get('native_source_file_sha256') for d in case_dirs]
 assert all(source_hashes) and len(set(source_hashes))==len(source_hashes), 'Native source bindings missing or duplicated'
 nr=[];rr=[];ee=[]
 for k,folder in zip(donor_keys,case_dirs):
  folder=Path(folder)
  for file,bag in [('neighborhood.tsv',nr),('regional_ES.tsv',rr),('conditional_LR_edges.tsv',ee)]:
   z=pd.read_csv(folder/file,sep='\t');z['private_donor_key']=k;bag.append(z)
 n=pd.concat(nr);r=pd.concat(rr);e=pd.concat(ee);results=[]
 for (radius,receiver),t in n.groupby(['radius_um','receiver']):
  d=t.log_odds_enrichment.dropna().to_numpy();p=paired_p_greater(d) if len(d)>=5 else np.nan
  results.append({'family':'neighborhood_primary2' if radius==50 and receiver in ['B','NK'] else 'neighborhood_sensitivity10','hypothesis':f'{receiver}_{radius}um','n_patients':len(d),'effect_median_log_odds':float(np.median(d)) if len(d) else np.nan,'n_positive':int((d>0).sum()),'P':p})
 for (sender,receiver,ligand,receptor),t in r.groupby(['sender','receiver','ligand','receptor']):
  use=t[t.status.eq('eligible_descriptive_effect')];d1=use.get('AB_minus_A',pd.Series(dtype=float)).to_numpy();d2=use.get('AB_minus_B',pd.Series(dtype=float)).to_numpy();p1=paired_p_greater(d1) if len(d1)>=5 else np.nan;p2=paired_p_greater(d2) if len(d2)>=5 else np.nan
  results.append({'family':'regional_LR27','hypothesis':f'{sender}>{receiver}:{ligand}_{receptor}','n_patients':len(d1),'AB_minus_A_median':np.median(d1) if len(d1) else np.nan,'AB_minus_B_median':np.median(d2) if len(d2) else np.nan,'P_vs_A':p1,'P_vs_B':p2,'P':max(p1,p2) if np.isfinite(p1+p2) else np.nan})
 for (sender,receiver,ligand,receptor),t in e.groupby(['sender','receiver','ligand','receptor']):
  d=t.loc[t.status.eq('eligible_conditional_spatial_null'),'excess_edge_fraction'].dropna().to_numpy();results.append({'family':'conditional_edge_LR27','hypothesis':f'{sender}>{receiver}:{ligand}_{receptor}','n_patients':len(d),'effect_median_excess_fraction':np.median(d) if len(d) else np.nan,'P':paired_p_greater(d) if len(d)>=5 else np.nan})
 t=pd.DataFrame(results);expected={'neighborhood_primary2':2,'neighborhood_sensitivity10':10,'regional_LR27':27,'conditional_edge_LR27':27}
 for family,total in expected.items():
  ix=t.index[t.family.eq(family)];assert len(ix)==total,(family,len(ix));p=t.loc[ix,'P'];t.loc[ix,'BH']=np.where(p.notna(),fdr(p.fillna(1)),np.nan);t.loc[ix,'BY']=np.where(p.notna(),fdr(p.fillna(1),True),np.nan)
 t.to_csv(out/'all66_patient_hypotheses.tsv',sep='\t',index=False);return t
