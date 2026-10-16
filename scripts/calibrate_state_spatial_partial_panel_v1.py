"""All frozen programmes x all panels: descriptive identity-internal proxy calibration."""
from pathlib import Path
from datetime import datetime,timezone
from collections import defaultdict
import os,json,hashlib,time,gc,shutil,importlib.metadata
import numpy as np,pandas as pd,anndata as ad
from scipy.stats import rankdata
import pyucell as uc

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def rho(a,b):
 if len(a)<3 or not np.isfinite(a).all() or not np.isfinite(b).all():return np.nan
 a=rankdata(a);b=rankdata(b);a-=a.mean();b-=b.mean();d=np.linalg.norm(a)*np.linalg.norm(b)
 return float(a@b/d) if d>1e-12 else np.nan
def score(a,sigs,cap):
 if not sigs:return {}
 uc.compute_ucell_scores(a,signatures=sigs,max_rank=cap,ties_method='average',missing_genes='skip',chunk_size=250,n_jobs=1,device='cpu',suffix='_proxy')
 return {k:a.obs[k+'_proxy'].to_numpy().copy() for k in sigs}
def resolved_panel(panel,native):
 aliases=defaultdict(list)
 for name,j in panel['symbol_to_feature_index'].items():
  if name in native:aliases[j].append(name)
 take=[]
 for j,names in sorted(aliases.items()):
  preferred=panel['preferred_symbol_by_feature'][j]
  if preferred in names:take.append(preferred)
  elif len(names)==1:take.append(names[0])
 repeated={g for g in take if take.count(g)>1}
 return sorted(g for g in take if g not in repeated)
V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);planp=Path('/tmp/state_spatial_partial_panel_calibration_plan_v1.json');plan=json.loads(planp.read_text());assert sha(__file__)==plan['producer_sha256'];assert importlib.metadata.version('pyucell')=='0.8.0';I=V/'state_spatial_proxy_native_cells_20261016_v1'
for tick in range(480):
 if (I/'inventory.json').exists():break
 e=V/'state_spatial_proxy_native_cells_20261016_v1.execution.json'
 if e.exists() and json.loads(e.read_text()).get('exit_code') is not None:raise RuntimeError('Native extraction terminal without inventory')
 time.sleep(30)
else:raise TimeoutError('Four-hour extraction dependency wait')
iv=json.loads((I/'inventory.json').read_text());J=V/'state_spatial_proxy_calibration_inputs_20261016_v1';ji=json.loads((J/'inventory.json').read_text());panp=J/'fixed_all14_panel_gene_axes.json';assert sha(panp)==ji['files'][panp.name];panels=json.loads(panp.read_text());assert len(panels)==14;D=V/'all_lineage_state_candidate_dictionary_20261014_v4/frozen_programme_dictionary.json';assert sha(D)==plan['dictionary_sha256'];recs={r['programme_id']:r for r in json.loads(D.read_text())['programmes']};O=V/'state_spatial_partial_panel_calibration_20261016_v1';assert not O.exists();O.mkdir(mode=0o700);(O/'private_cell_scores').mkdir(mode=0o700);(O/'private_identity_metrics').mkdir(mode=0o700);summary=[];files=[];maps=[];nc=0;checks=[]
for info in iv['private_native_files']:
 assert shutil.disk_usage(V).free>500*2**30;p=I/info['relative_file'];assert sha(p)==info['sha256'];a=ad.read_h5ad(p);assert a.n_obs==info['cells'];ids=list(a.uns['programme_ids']);assert all(info['major'] in recs[k]['allowed_majors'] for k in ids);native=set(a.var_names);reference=a.obsm['accepted_full_native_UCell1500'];axes={}
 for pan in panels:
  axis=pan['gene_axis_sha256'];resolved=resolved_panel(pan,native)
  if axis in axes:assert axes[axis]['genes']==resolved
  else:axes[axis]=dict(genes=resolved,samples=[])
  axes[axis]['samples'].append(pan['sample_id']);maps.append(dict(gse_id=info['gse_id'],major=info['major'],sample_id=pan['sample_id'],axis_id=axis,native_resolved_panel_features=len(resolved),original_features=pan['features'],panel_missing_not_zero_filled=True))
 assert len(axes)==4;outside_signatures={};outside_keys={};partial={}
 for ax,axinfo in sorted(axes.items()):
  bg=set(axinfo['genes']);partial[ax]={}
  for k in ids:
   original=recs[k]['primary_genes'];inside=[g for g in original if g in bg];remaining=[g for g in original if g not in bg];outside=[g for g in remaining if g in native];partial[ax][k]=inside
   if len(inside) and len(outside)>=10 and len(outside)/max(1,len(remaining))>=.8:
    key='outside_'+ax[:10]+'_'+k;outside_signatures[key]=outside;outside_keys[ax,k]=key
 outs=score(a,outside_signatures,1500);groupcodes,groupnames=pd.factorize(pd.MultiIndex.from_frame(a.obs[['node_id','stratum']]),sort=True);blockcodes,blocknames=pd.factorize(a.obs.effective_resampling_block.astype(str),sort=True);tot=a.obs.full_RNA_total.to_numpy(float);det=a.obs.full_detected_genes.to_numpy(float)
 for ax,axinfo in sorted(axes.items()):
  names=axinfo['genes'];cap=min(1500,len(names));proxy=np.full(reference.shape,np.nan,np.float32);outarr=np.full(reference.shape,np.nan,np.float32);support=np.zeros(reference.shape,np.uint8);frac=np.zeros(reference.shape,np.float32);sig={k:partial[ax][k] for k in ids if len(partial[ax][k]) and cap>len(partial[ax][k])};panel=a[:,names].copy();scores=score(panel,sig,cap) if cap>=2 else {};panelRNA=np.asarray(panel.X.sum(1)).ravel()
  for j,k in enumerate(ids):
   if k in scores:proxy[:,j]=scores[k]
   if (ax,k) in outside_keys:outarr[:,j]=outs[outside_keys[ax,k]]
   gs=partial[ax][k]
   if gs:x=a.X[:,a.var_names.get_indexer(gs)];support[:,j]=np.asarray((x>0).sum(1)).ravel();frac[:,j]=np.asarray(x.sum(1)).ravel()/tot
  # Independent manual tied-rank checks for three actual cells on each resolved panel.
  chosen=np.unique(np.linspace(0,a.n_obs-1,min(3,a.n_obs),dtype=int));err=0;num=0
  for cell in chosen:
   row=panel.X[cell].toarray().ravel().astype(float);ranks=np.full(len(row),cap,float);nz=np.flatnonzero(row>0);ranks[nz]=np.minimum(cap,rankdata(-row[nz],method='average'));lookup={g:i for i,g in enumerate(names)}
   for j,k in enumerate(ids):
    gs=partial[ax][k]
    if k not in scores:continue
    n=len(gs);base=n*(n+1)/2;manual=1-(sum(float(ranks[lookup[g]]) for g in gs)-base)/(n*cap-base);err=max(err,abs(manual-float(proxy[cell,j])));num+=1
  assert err<1e-6;checks.append(dict(gse_id=info['gse_id'],major=info['major'],axis_id=ax,manual_scores=num,max_manual_error=err));records=[]
  for code,gn in enumerate(groupnames):
   ii=np.flatnonzero(groupcodes==code);blocks=blockcodes[ii];unique=np.unique(blocks);n=len(ii);nb=len(unique);blockix=[np.flatnonzero(blocks==v) for v in unique]
   for j,k in enumerate(ids):
    ref=reference[ii,j];pr=proxy[ii,j];out=outarr[ii,j];npresent=len(partial[ax][k]);native_ok=np.isfinite(ref).all();okay=n>=plan['min_cells_per_identity_stratum'] and nb>=plan['min_blocks_per_identity_stratum']
    if not native_ok:status='FULL_NATIVE_REFERENCE_UNMEASURABLE'
    elif not npresent:status='NO_SIGNATURE_GENE_ON_NATIVE_PANEL'
    elif not np.isfinite(pr).all():status='PARTIAL_PANEL_SCORE_UNDEFINED'
    elif not okay:status='FEWER_THAN20CELLS_OR5BLOCKS_TECHNICAL_COVERAGE'
    elif np.ptp(pr)<=1e-12 or np.ptp(ref)<=1e-12:status='CONSTANT_PROXY_OR_REFERENCE'
    else:status='DESCRIPTIVE_IDENTITY_INTERNAL_CALIBRATION'
    usable=status=='DESCRIPTIVE_IDENTITY_INTERNAL_CALIBRATION';block_ref=np.array([ref[x].mean() for x in blockix]);block_proxy=np.array([pr[x].mean() for x in blockix]);block_out=np.array([out[x].mean() for x in blockix]);rows=dict(gse_id=info['gse_id'],major=info['major'],node_id=str(gn[0]),stratum=gn[1],programme_id=k,axis_id=ax,cells=n,analysis_blocks=nb,partial_signature_genes=npresent,original_signature_genes=len(recs[k]['primary_genes']),native_panel_features=len(names),proxy_maxrank=cap,status=status,cell_full_reference_rho=rho(ref,pr) if usable else np.nan,cell_outside_panel_reference_rho=rho(out,pr) if usable else np.nan,block_mean_full_reference_rho=rho(block_ref,block_proxy) if usable else np.nan,block_mean_outside_panel_reference_rho=rho(block_out,block_proxy) if usable else np.nan,proxy_full_RNA_rho=rho(pr,tot[ii]) if usable else np.nan,proxy_panel_RNA_rho=rho(pr,panelRNA[ii]) if usable else np.nan,proxy_detected_genes_rho=rho(pr,det[ii]) if usable else np.nan,proxy_IQR=float(np.quantile(pr,.75)-np.quantile(pr,.25)) if np.isfinite(pr).all() else np.nan,full_reference_IQR=float(np.quantile(ref,.75)-np.quantile(ref,.25)) if native_ok else np.nan,partial_support_atleast1_fraction=float((support[ii,j]>=1).mean()),partial_support_atleast5_fraction=float((support[ii,j]>=5).mean()),partial_signature_full_RNA_fraction=float(np.median(frac[ii,j])),full_signature_spatial_qualified=False,biological_validation=False);records.append(rows)
  tab=pd.DataFrame(records);rel=f'private_identity_metrics/{Path(info["relative_file"]).stem}_{ax[:10]}.tsv.gz';tab.to_csv(O/rel,sep='\t',index=False,compression='gzip');files.append(dict(relative_file=rel,sha256=sha(O/rel),records=len(tab)));rel=f'private_cell_scores/{Path(info["relative_file"]).stem}_{ax[:10]}.npz';np.savez_compressed(O/rel,global_rows=a.obs.core_position.to_numpy(),unit_index=a.obs.unit_index.to_numpy(),programme_ids=np.array(ids),proxy=proxy,outside_reference=outarr,partial_detected_genes=support,partial_full_RNA_fraction=frac);files.append(dict(relative_file=rel,sha256=sha(O/rel),cells=a.n_obs))
  for k,g in tab.groupby('programme_id',sort=True):
   use=g[g.status.eq('DESCRIPTIVE_IDENTITY_INTERNAL_CALIBRATION')];item=dict(gse_id=info['gse_id'],major=info['major'],programme_id=k,axis_id=ax,identity_strata=len(g),calibratable_identity_strata=len(use),partial_signature_genes=int(g.partial_signature_genes.iloc[0]),native_panel_features=len(names),proxy_maxrank=cap,full_signature_spatial_qualified=False,biological_validation=False)
   for col in ['cell_full_reference_rho','cell_outside_panel_reference_rho','block_mean_full_reference_rho','block_mean_outside_panel_reference_rho','proxy_full_RNA_rho','proxy_panel_RNA_rho','proxy_detected_genes_rho']:
    vals=use[col].dropna();item[col+'_median']=float(vals.median()) if len(vals) else np.nan;item[col+'_estimable_strata']=len(vals)
   for status,count in g.status.value_counts().items():item['status_'+status]=int(count)
   summary.append(item)
  del panel,proxy,outarr,support,frac,scores;gc.collect()
 nc+=a.n_obs;(O/'progress.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),cells_calibrated=nc,source_majors_completed=len({(x['gse_id'],x['major']) for x in summary}),summary_records=len(summary)))+'\n');del a,outs;gc.collect()
assert nc==iv['cells'];pd.DataFrame(summary).to_csv(O/'public_all_source_programme_axis_summary.tsv',sep='\t',index=False);pd.DataFrame(maps).to_csv(O/'public_all14_native_axis_resolution.tsv',sep='\t',index=False);pd.DataFrame(checks).to_csv(O/'public_actual_cell_manual_proxy_checks.tsv',sep='\t',index=False);out=dict(status='ALL_FROZEN_NONDISCOVERY_CELLS_ALL814_PROGRAMMES_ALL14_PANELS_DESCRIPTIVE_PROXY_CALIBRATION_NOT_SPATIAL_VALIDATION',observed_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(planp),native_input_inventory_sha256=sha(I/'inventory.json'),selection_inventory_sha256=sha(J/'inventory.json'),cells=nc,programmes=814,panels=14,unique_symbol_axes=4,public_summary_records=len(summary),private_outputs=files,public_outputs={p.name:sha(p) for p in O.glob('public_*.tsv')},manual_proxy_scores=sum(x['manual_scores'] for x in checks),max_manual_proxy_error=max(x['max_manual_error'] for x in checks),formal_states_or_spatial_validation=False);(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ['private_outputs','public_outputs']}))
