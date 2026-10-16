"""Float64 block means only; preserve every other TSV field and all score bytes."""
from pathlib import Path
from datetime import datetime,timezone
from collections import defaultdict
import os,json,hashlib,shutil,gzip
import numpy as np,pandas as pd,anndata as ad
from scipy.stats import rankdata
V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);D=V/'state_spatial_partial_panel_calibration_20261016_v1';I=V/'state_spatial_proxy_native_cells_20261016_v1';sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();di=json.loads((D/'inventory.json').read_text());ii=json.loads((I/'inventory.json').read_text());O=V/'state_spatial_partial_panel_calibration_20261016_v2';assert not O.exists();O.mkdir(mode=0o700);(O/'private_cell_scores').mkdir(mode=0o700);(O/'private_identity_metrics').mkdir(mode=0o700);source_files={x['relative_file']:x for x in di['private_outputs']};pubname='public_all_source_programme_axis_summary.tsv';pubheader=(D/pubname).read_text().splitlines()[0].split('\t');pubrows=[x.split('\t') for x in (D/pubname).read_text().splitlines()[1:]];assert all(len(x)==len(pubheader) for x in pubrows);pidx={x:i for i,x in enumerate(pubheader)};pkeys=['gse_id','major','programme_id','axis_id'];pindex={tuple(r[pidx[k]] for k in pkeys):i for i,r in enumerate(pubrows)};assert len(pindex)==len(pubrows);newpub=[x.copy() for x in pubrows];files=[];changes=0;public_changes=0;checked=0
metrics=['block_mean_full_reference_rho','block_mean_outside_panel_reference_rho']
def vcorr(x,y):
 a=rankdata(x,axis=0,nan_policy='propagate');b=rankdata(y,axis=0,nan_policy='propagate');a-=a.mean(0);b-=b.mean(0);den=np.linalg.norm(a,axis=0)*np.linalg.norm(b,axis=0);return np.divide((a*b).sum(0),den,out=np.full(a.shape[1],np.nan),where=den>1e-12)
def fmt(x):return str(float(x)) if np.isfinite(x) else ''
for info in ii['private_native_files']:
 p=I/info['relative_file'];assert sha(p)==info['sha256'];a=ad.read_h5ad(p);ids=list(a.uns['programme_ids']);stem=p.stem;groups=a.obs.groupby(['node_id','stratum'],sort=True,observed=True).indices;ref=a.obsm['accepted_full_native_UCell1500'].astype(float)
 for cp in sorted((D/'private_cell_scores').glob(stem+'_*.npz')):
  rel=str(cp.relative_to(D));assert sha(cp)==source_files[rel]['sha256'];shutil.copyfile(cp,O/rel);assert sha(O/rel)==sha(cp);files.append({**source_files[rel],'sha256':sha(O/rel)});z=np.load(cp,allow_pickle=False);proxy=z['proxy'].astype(float);outside=z['outside_reference'].astype(float);ax=None;mp=D/'private_identity_metrics'/(cp.stem+'.tsv.gz');rel=str(mp.relative_to(D));assert sha(mp)==source_files[rel]['sha256']
  with gzip.open(mp,'rt') as f:lines=f.read().splitlines()
  header=lines[0].split('\t');ix={k:j for j,k in enumerate(header)};rows=[s.split('\t') for s in lines[1:]];assert all(len(s)==len(header) for s in rows);lookup={(s[ix['node_id']],s[ix['stratum']],s[ix['programme_id']]):j for j,s in enumerate(rows)};assert len(lookup)==len(rows);new=[s.copy() for s in rows];vals=defaultdict(lambda:defaultdict(list));ax=rows[0][ix['axis_id']]
  for (node,stratum),take in groups.items():
   b=a.obs.effective_resampling_block.astype(str).to_numpy()[take];idx=[take[np.flatnonzero(b==k)] for k in np.unique(b)];pr=np.array([proxy[j].mean(0,dtype=np.float64) for j in idx]);rr=np.array([ref[j].mean(0,dtype=np.float64) for j in idx]);oo=np.array([outside[j].mean(0,dtype=np.float64) for j in idx]);calc=[vcorr(pr,rr),vcorr(pr,oo)]
   for j,k in enumerate(ids):
    ri=lookup[str(node),stratum,k];eligible=rows[ri][ix['status']]=='DESCRIPTIVE_IDENTITY_INTERNAL_CALIBRATION'
    for metric,v in zip(metrics,calc):
     value=float(v[j]) if eligible else np.nan;old=rows[ri][ix[metric]];new[ri][ix[metric]]=fmt(value)
     if old!=fmt(value):changes+=1
     if eligible and np.isfinite(value):vals[k][metric].append(value)
  allowed={ix[x] for x in metrics};assert all(all(old[j]==r[j] for j in range(len(header)) if j not in allowed) for old,r in zip(rows,new))
  with gzip.open(O/rel,'wt') as f:f.write('\t'.join(header)+'\n');f.writelines('\t'.join(r)+'\n' for r in new)
  files.append({**source_files[rel],'sha256':sha(O/rel)})
  for k in ids:
   ri=pindex[info['gse_id'],info['major'],k,ax]
   for metric in metrics:
    vs=vals[k][metric];newpub[ri][pidx[metric+'_median']]=fmt(np.median(vs)) if vs else '';newpub[ri][pidx[metric+'_estimable_strata']]=str(len(vs));public_changes+=int(newpub[ri][pidx[metric+'_median']]!=pubrows[ri][pidx[metric+'_median']])
  checked+=len(rows)
 (O/'progress.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),identity_records_checked=checked,block_metric_fields_rewritten=changes))+'\n')
allowed={pidx[x+s] for x in metrics for s in ['_median','_estimable_strata']};assert all(all(old[j]==r[j] for j in range(len(pubheader)) if j not in allowed) for old,r in zip(pubrows,newpub));(O/pubname).write_text('\t'.join(pubheader)+'\n'+''.join('\t'.join(r)+'\n' for r in newpub))
for name,h in di['public_outputs'].items():
 if name!=pubname:assert sha(D/name)==h;shutil.copyfile(D/name,O/name)
out={**di,'status':'ALL91447_PARTIAL_PANEL_SCORES_UNCHANGED_BLOCK_MEANS_FLOAT64_REPAIR_PENDING_INDEPENDENT_AUDIT','observed_UTC':datetime.now(timezone.utc).isoformat(),'producer_sha256':sha(__file__),'parent_inventory_sha256':sha(D/'inventory.json'),'parent_directory':D.name,'all_cell_score_files_byte_identical':True,'all_other_TSV_fields_exactly_preserved':True,'block_mean_precision':'float64 accumulation of unchanged float32 cell scores','identity_records_checked':checked,'block_metric_fields_rewritten':changes,'source_summary_block_median_fields_rewritten':public_changes,'private_outputs':files,'public_outputs':{name:sha(O/name) for name in di['public_outputs']}};(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ['private_outputs','public_outputs']}))
