"""Independent sparse tied-rank all-cell proxy/full/outside scores and all descriptors."""
from pathlib import Path
from datetime import datetime,timezone
from collections import defaultdict
import os,json,hashlib,time,gc
import numpy as np,pandas as pd,anndata as ad
from scipy.stats import rankdata
from scipy.sparse import csr_matrix

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def deltarank(X,cap):
 X=X.tocsr();data=np.zeros(X.nnz,float)
 for i in range(X.shape[0]):
  lo,hi=X.indptr[i:i+2];v=X.data[lo:hi];ix=np.flatnonzero(v>0);data[lo+ix]=np.minimum(cap,rankdata(-v[ix].astype(float),method='average'))-cap
 r=csr_matrix((data,X.indices.copy(),X.indptr.copy()),shape=X.shape);r.eliminate_zeros();return r
def manual(R,genes,sigs,cap):
 lookup={g:i for i,g in enumerate(genes)};row=[];col=[];ns=[]
 for j,s in enumerate(sigs):
  ns.append(len(s))
  for g in s:row.append(lookup[g]);col.append(j)
 S=csr_matrix((np.ones(len(row)),(row,col)),shape=(len(genes),len(sigs)));n=np.array(ns);base=n*(n+1)/2;de=(R@S).toarray();y=1-(n*cap+de-base)/(n*cap-base);y[:,n==0]=np.nan;return y,S

def vcorr(x,y):
 x=np.asarray(x,float);y=np.asarray(y,float)
 if y.ndim==1:y=np.repeat(y[:,None],x.shape[1],axis=1)
 a=rankdata(x,axis=0,nan_policy='propagate');b=rankdata(y,axis=0,nan_policy='propagate');a-=a.mean(0);b-=b.mean(0);den=np.linalg.norm(a,axis=0)*np.linalg.norm(b,axis=0);return np.divide((a*b).sum(0),den,out=np.full(x.shape[1],np.nan),where=den>1e-12)
V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);D=V/'state_spatial_partial_panel_calibration_20261016_v1'
for tick in range(480):
 if (D/'inventory.json').exists():break
 e=V/'state_spatial_partial_panel_calibration_20261016_v1.execution.json'
 if e.exists() and json.loads(e.read_text()).get('exit_code') is not None:raise RuntimeError('Calibration terminated without result inventory')
 time.sleep(30)
else:raise TimeoutError('Four-hour calibration dependency wait')
iv=json.loads((D/'inventory.json').read_text());I=V/'state_spatial_proxy_native_cells_20261016_v1';ii=json.loads((I/'inventory.json').read_text());assert sha(I/'inventory.json')==iv['native_input_inventory_sha256'];J=V/'state_spatial_proxy_calibration_inputs_20261016_v1';panels=json.loads((J/'fixed_all14_panel_gene_axes.json').read_text());dictionary=json.loads((V/'all_lineage_state_candidate_dictionary_20261014_v4/frozen_programme_dictionary.json').read_text());recs={r['programme_id']:r for r in dictionary['programmes']};pub=pd.read_csv(D/'public_all_source_programme_axis_summary.tsv',sep='\t').set_index(['gse_id','major','programme_id','axis_id']);files={x['relative_file']:x for x in iv['private_outputs']}
for name,h in iv['public_outputs'].items():assert sha(D/name)==h
O=V/'state_spatial_partial_panel_independent_audit_20261016_v1';assert not O.exists();O.mkdir(mode=0o700);seen=0;scorevalues=0;referencevalues=0;descriptor_records=0;max_score_error=0;max_descriptor_error=0;checks=[];checked_public=set();assert pub.index.is_unique
for info in ii['private_native_files']:
 p=I/info['relative_file'];assert sha(p)==info['sha256'];a=ad.read_h5ad(p);ids=list(a.uns['programme_ids']);genes=list(a.var_names);native=set(genes);ref=a.obsm['accepted_full_native_UCell1500'];Rfull=deltarank(a.X,1500);fsigs=[]
 for k in ids:
  original=recs[k]['primary_genes'];present=[g for g in original if g in native];fsigs.append(present if len(present)>=10 and len(present)/len(original)>=.8 else [])
 full,_=manual(Rfull,genes,fsigs,1500);assert np.allclose(full,ref,atol=1e-6,equal_nan=True);referencevalues+=int(np.isfinite(ref).sum());max_score_error=max(max_score_error,float(np.nanmax(abs(full-ref))) if np.isfinite(ref).any() else 0);axgroups=defaultdict(list)
 for pan in panels:axgroups[pan['gene_axis_sha256']].append(pan)
 for ax,pp in sorted(axgroups.items()):
  pan=pp[0];features=defaultdict(set)
  for gene,j in pan['symbol_to_feature_index'].items():
   if gene in native:features[j].add(gene)
  resolved=[]
  for j,names in sorted(features.items()):
   preferred=pan['preferred_symbol_by_feature'][j]
   if preferred in names:resolved.append(preferred)
   elif len(names)==1:resolved.append(next(iter(names)))
  names=sorted(g for g in resolved if resolved.count(g)==1);cap=min(1500,len(names));sigs=[[g for g in recs[k]['primary_genes'] if g in set(names)] for k in ids];panel=a[:,names];R=deltarank(panel.X,cap);proxy,S=manual(R,names,sigs,cap);outside=[]
  for k,s in zip(ids,sigs):
   remaining=[g for g in recs[k]['primary_genes'] if g not in set(names)];present=[g for g in remaining if g in native];outside.append(present if len(s) and len(present)>=10 and len(present)/max(1,len(remaining))>=.8 else [])
  out,_=manual(Rfull,genes,outside,1500);stem=Path(info['relative_file']).stem;rel=f'private_cell_scores/{stem}_{ax[:10]}.npz';p=D/rel;assert sha(p)==files[rel]['sha256'];z=np.load(p,allow_pickle=False);assert np.array_equal(z['global_rows'],a.obs.core_position) and np.array_equal(z['unit_index'],a.obs.unit_index) and list(z['programme_ids'])==ids;assert np.allclose(proxy,z['proxy'],atol=1e-6,equal_nan=True) and np.allclose(out,z['outside_reference'],atol=1e-6,equal_nan=True);er1=float(np.nanmax(abs(proxy-z['proxy']))) if np.isfinite(proxy).any() else 0;er2=float(np.nanmax(abs(out-z['outside_reference']))) if np.isfinite(out).any() else 0;max_score_error=max(max_score_error,er1,er2);scorevalues+=int(np.isfinite(proxy).sum()+np.isfinite(out).sum());support=((panel.X>0).astype(int)@S).toarray();fraction=(panel.X@S).toarray()/a.obs.full_RNA_total.to_numpy()[:,None];assert np.array_equal(support,z['partial_detected_genes']) and np.allclose(fraction,z['partial_full_RNA_fraction'],atol=1e-7)
  rel=f'private_identity_metrics/{stem}_{ax[:10]}.tsv.gz';p=D/rel;assert sha(p)==files[rel]['sha256'];tab=pd.read_csv(p,sep='\t');tab.node_id=tab.node_id.astype(str);indexed=tab.set_index(['node_id','stratum','programme_id']);assert len(indexed)==len(tab) and indexed.index.is_unique;groups=a.obs.groupby(['node_id','stratum'],sort=True,observed=True).indices;expected_status={}
  for (node,stratum),ix in groups.items():
   n=len(ix);blocks=a.obs.effective_resampling_block.astype(str).to_numpy()[ix];uniq=np.unique(blocks);bix=[np.flatnonzero(blocks==b) for b in uniq];nr=len(uniq);pr=z['proxy'][ix].astype(float);rr=ref[ix].astype(float);oo=z['outside_reference'][ix].astype(float);RNA=a.obs.full_RNA_total.to_numpy()[ix];panelRNA=np.asarray(panel.X[ix].sum(1)).ravel();det=a.obs.full_detected_genes.to_numpy()[ix];br=np.array([rr[j].mean(0) for j in bix]);bp=np.array([pr[j].mean(0) for j in bix]);bo=np.array([oo[j].mean(0) for j in bix]);calc={'cell_full_reference_rho':vcorr(rr,pr),'cell_outside_panel_reference_rho':vcorr(oo,pr),'block_mean_full_reference_rho':vcorr(br,bp),'block_mean_outside_panel_reference_rho':vcorr(bo,bp),'proxy_full_RNA_rho':vcorr(pr,RNA),'proxy_panel_RNA_rho':vcorr(pr,panelRNA),'proxy_detected_genes_rho':vcorr(pr,det)}
   extra=dict(proxy_IQR=np.quantile(pr,.75,axis=0)-np.quantile(pr,.25,axis=0),full_reference_IQR=np.quantile(rr,.75,axis=0)-np.quantile(rr,.25,axis=0),partial_support_atleast1_fraction=(support[ix]>=1).mean(0),partial_support_atleast5_fraction=(support[ix]>=5).mean(0),partial_signature_full_RNA_fraction=np.median(fraction[ix],axis=0))
   for j,k in enumerate(ids):
    r=indexed.loc[str(node),stratum,k];nativeok=np.isfinite(rr[:,j]).all();has=len(sigs[j]);status='FULL_NATIVE_REFERENCE_UNMEASURABLE' if not nativeok else 'NO_SIGNATURE_GENE_ON_NATIVE_PANEL' if not has else 'PARTIAL_PANEL_SCORE_UNDEFINED' if not np.isfinite(pr[:,j]).all() else 'FEWER_THAN20CELLS_OR5BLOCKS_TECHNICAL_COVERAGE' if n<20 or nr<5 else 'CONSTANT_PROXY_OR_REFERENCE' if np.ptp(pr[:,j])<=1e-12 or np.ptp(rr[:,j])<=1e-12 else 'DESCRIPTIVE_IDENTITY_INTERNAL_CALIBRATION';assert r.status==status and r.cells==n and r.analysis_blocks==nr and r.partial_signature_genes==has and not r.biological_validation and not r.full_signature_spatial_qualified;usable=status=='DESCRIPTIVE_IDENTITY_INTERNAL_CALIBRATION'
    for key,vals in calc.items():
     value=vals[j] if usable else np.nan;assert np.isclose(value,r[key],atol=1e-8,equal_nan=True),(key,value,r[key])
     if np.isfinite(value):max_descriptor_error=max(max_descriptor_error,abs(value-r[key]))
    for metric,vals in extra.items():assert np.isclose(vals[j],r[metric],atol=1e-7,equal_nan=True)
    assert r.native_panel_features==len(names) and r.proxy_maxrank==cap
    descriptor_records+=1
  for k,g in tab.groupby('programme_id',sort=True):
   pubkey=(info['gse_id'],info['major'],k,ax);checked_public.add(pubkey);row=pub.loc[pubkey];use=g[g.status.eq('DESCRIPTIVE_IDENTITY_INTERNAL_CALIBRATION')];assert row.identity_strata==len(g) and row.calibratable_identity_strata==len(use)
   for metric in ['cell_full_reference_rho','cell_outside_panel_reference_rho','block_mean_full_reference_rho','block_mean_outside_panel_reference_rho','proxy_full_RNA_rho','proxy_panel_RNA_rho','proxy_detected_genes_rho']:
    v=use[metric].dropna();assert row[metric+'_estimable_strata']==len(v);assert np.isclose(row[metric+'_median'],v.median() if len(v) else np.nan,atol=1e-8,equal_nan=True)
   for status,count in g.status.value_counts().items():assert row['status_'+status]==count
  checks.append(dict(gse_id=info['gse_id'],major=info['major'],axis_id=ax,cells=a.n_obs,manual_proxy_max_error=er1,manual_outside_reference_max_error=er2,identity_descriptor_records=len(tab)));del panel,R,proxy,out,z;gc.collect()
 seen+=a.n_obs;(O/'progress.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),cells_checked=seen,manual_partial_and_outside_score_values=scorevalues,manual_full_reference_values=referencevalues,identity_descriptor_records=descriptor_records))+'\n');del a,Rfull,full;gc.collect()
assert seen==ii['cells']==91447 and checked_public==set(pub.index);pd.DataFrame(checks).to_csv(O/'public_all_source_axis_independent_checks.tsv',sep='\t',index=False);result=dict(status='PASS_ALL91447_CELL_SPARSE_MANUAL_FULL_PARTIAL_OUTSIDE_UCELL_AND_ALL_IDENTITY_RHOS_SCOPE_SOURCE_SUMMARIES',observed_UTC=datetime.now(timezone.utc).isoformat(),auditor_sha256=sha(__file__),calibration_inventory_sha256=sha(D/'inventory.json'),cells_checked=seen,manual_partial_and_outside_score_values=scorevalues,manual_full_reference_values=referencevalues,identity_descriptor_records=descriptor_records,max_manual_score_error=max_score_error,max_descriptor_error=max_descriptor_error,public_output_sha256=sha(O/'public_all_source_axis_independent_checks.tsv'),formal_states_or_biological_validation=False);(O/'inventory.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
