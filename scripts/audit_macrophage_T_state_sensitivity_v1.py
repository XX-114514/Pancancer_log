"""Independent family/aggregation audit and balanced direct full-X CR2 fixtures.
Not an independent biological validation or complete replay of every coefficient.
"""
import os,json,hashlib
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd
from scipy.stats import t

def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def bh(p):
 p=np.asarray(p,float);order=np.argsort(p,kind='stable');q=np.empty(len(p));q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return q

def explicit(Z,x,Y,blocks):
 X=np.column_stack((Z,x));B=np.linalg.pinv(X,rcond=1e-10);coef=B@Y;c=B[-1];H=X@B;E=Y-X@coef;vectors=[];scores=[]
 for block in np.unique(blocks):
  ix=np.flatnonzero(blocks==block);e,V=np.linalg.eigh(np.eye(len(ix))-H[np.ix_(ix,ix)])
  if e.min()<=1e-8:return None
  A=(V/np.sqrt(e))@V.T;v=np.zeros(len(x));v[ix]=A@c[ix];vectors.append((np.eye(len(x))-H)@v);scores.append(c[ix]@A@E[ix])
 P=np.column_stack(vectors);C=P.T@P;df=np.trace(C)**2/np.sum(C*C);beta=coef[-1];se=np.sqrt(np.sum(np.array(scores)**2,axis=0));p=2*t.sf(abs(beta/se),df);crit=t.ppf(.975,df)
 return dict(beta=beta,se=se,df=df,p_exploratory=p,lower=beta-crit*se,upper=beta+crit*se)

def main():
 root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])/'immune_same_name_family_source_review_v1/macrophage_T_state_sensitivity_v1';plan=json.loads(Path(os.environ['PANCANCER_MAC_T_STATE_PLAN']).read_text());assert digest(__file__)==plan['auditor_sha256'];inv=json.loads((root/'inventory.json').read_text())
 for f,h in inv['outputs'].items():assert digest(root/f)==h
 read=lambda n:pd.read_csv(root/n,sep='\t',keep_default_na=True)
 models=read('all27744_pooled_LR_models.tsv.gz');comp=read('all_pooled_TF_components.tsv.gz');src=read('all_source_TF_effects.tsv.gz');loso=read('all_LOSO_TF_effects.tsv.gz');keys=['direction','variant','LR_index'];assert len(models)==27744 and not models.duplicated(keys).any();assert not comp.duplicated(keys+['TF']).any();assert not src.duplicated(keys+['TF','gse_id']).any();assert not loso.duplicated(keys+['TF','omitted_gse']).any()
 assert set(models.direction)==set(plan['direction_contexts']);assert (comp.n_GSEs>=2).all();assert comp.max_source_information_share.between(0,1).all()
 grouped={k:g for k,g in comp.groupby(keys)}
 for row in models.itertuples(index=False):
  k=(row.direction,row.variant,row.LR_index);g=grouped.get(k)
  assert row.eligible_TFs==(0 if g is None else len(g))
  if g is not None:
   valid=g.dropna(subset=['p_exploratory'])
   if len(valid):
    lead=valid.sort_values('p_exploratory',kind='stable').iloc[0];assert np.isclose(row.joint_p_exploratory,min(1,len(g)*lead.p_exploratory),rtol=1e-10);assert row.lead_TF==lead.TF
   else:assert pd.isna(row.joint_p_exploratory)
  else:assert pd.isna(row.joint_p_exploratory)
 for variant,g in models.groupby('variant'):
  assert len(g)==9248
  for _,d in g.groupby('direction'):assert set(d.LR_index)==set(range(4624))
  q=bh(g.joint_p_exploratory.fillna(1));valid=g.status.eq('pooled_exploratory_estimated');np.testing.assert_allclose(g.loc[valid,'joint_BH_exploratory'],q[valid],atol=1e-12);np.testing.assert_allclose(g.loc[valid,'joint_BY_exploratory'],np.minimum(1,q[valid]*np.sum(1/np.arange(1,9249))),atol=1e-12)
 cindex=comp.set_index(keys+['TF']);sindex=src.set_index(keys+['TF','gse_id']);lindex=loso.set_index(keys+['TF','omitted_gse']);fixture_counts={};checks=0;maximum=0.
 for file in sorted((root/'private_pooled_input_cases').glob('*.npz')):
  a=np.load(file);direction=str(a['direction']);variant=str(a['variant']);li=int(a['LR_index']);key=(direction,variant,li);fixture_counts[(direction,variant)]=fixture_counts.get((direction,variant),0)+1
  def compare(effect,table,suffix=()):
   nonlocal checks,maximum
   if effect is None:return
   for j,tf in enumerate(a['TFs']):
    row=table.loc[key+(str(tf),)+suffix]
    for metric in ['beta','se','df','lower','upper','p_exploratory']:
     if metric not in row.index:continue
     expected=effect[metric] if metric=='df' else effect[metric][j];actual=float(row[metric]);np.testing.assert_allclose(actual,expected,rtol=2e-7,atol=2e-9,equal_nan=True);maximum=max(maximum,abs(actual-expected))
    checks+=1
  compare(explicit(a['Z'],a['x'],a['Y'],a['blocks']),cindex)
  for gse in np.unique(a['gse']):
   m=a['gse']==gse;compare(explicit(a['Z'][m],a['x'][m],a['Y'][m],a['blocks'][m]),sindex,(str(gse),));compare(explicit(a['Z'][~m],a['x'][~m],a['Y'][~m],a['blocks'][~m]),lindex,(str(gse),))
 assert sum(fixture_counts.values())==plan['audit_input_cases'] and len(fixture_counts)==6 and all(v==4 for v in fixture_counts.values())
 result=dict(status='PASS_FAMILY_AND_AGGREGATION_WITH_BALANCED_NUMERIC_FIXTURES',finished_UTC=datetime.now(timezone.utc).isoformat(),inventory_sha256=digest(root/'inventory.json'),models=len(models),TF_components=len(comp),numeric_fixture_groups=24,numeric_component_checks=checks,max_absolute_difference=maximum,all_coefficients_independently_replayed=False,source_and_LOSO_numeric_checks_included=True,auditor_sha256=digest(__file__))
 (root/'independent_audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':main()
