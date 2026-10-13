"""Independent full-family BH/coverage checks and direct full-X numeric fixtures."""
from pathlib import Path
import os,json,hashlib,importlib.util
from datetime import datetime,timezone
import numpy as np,pandas as pd

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def bh(p):
 p=np.asarray(p,float);ix=np.argsort(p,kind='stable');q=np.empty(len(p));q[ix]=np.minimum(1,np.minimum.accumulate((p[ix]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return q

def main():
 root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])/'immune_same_name_family_source_review_v1/mac_T_receiver_programs_v1';plan=json.loads(Path(os.environ['PANCANCER_RECEIVER_PROGRAM_PLAN']).read_text());assert sha(__file__)==plan['auditor_sha256'];ref=Path(os.environ['PANCANCER_EXPLICIT_POOL_AUDITOR']);assert sha(ref)==plan['explicit_auditor_sha256'];spec=importlib.util.spec_from_file_location('ref',ref);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);inv=json.loads((root/'inventory.json').read_text())
 for f,h in inv['outputs'].items():assert sha(root/f)==h
 df=pd.read_csv(root/'all147968_receiver_program_models.tsv.gz',sep='\t');keys=['direction','variant','LR_index','program'];assert len(df)==147968 and not df.duplicated(keys).any()
 for (direction,variant,program),g in df.groupby(keys[:2]+['program']):assert len(g)==4624 and set(g.LR_index)==set(range(4624))
 for v,g in df.groupby('variant'):
  assert len(g)==36992;q=bh(g.p_exploratory.fillna(1));ok=g.status.eq('inference_estimable');np.testing.assert_allclose(g.loc[ok,'BH_exploratory'],q[ok],atol=1e-12);np.testing.assert_allclose(g.loc[ok,'BY_exploratory'],np.minimum(1,q[ok]*sum(1/np.arange(1,36993))),atol=1e-12);assert g.loc[ok,['beta','se','lower','upper','p_exploratory']].notna().all().all()
 src=pd.read_csv(root/'all_source_effects.tsv.gz',sep='\t');leave=pd.read_csv(root/'all_LOSO_effects.tsv.gz',sep='\t');assert not src.duplicated(keys+['gse_id']).any() and not leave.duplicated(keys+['omitted_gse']).any();lookup=df.set_index(keys);sx=src.set_index(keys+['gse_id']);lx=leave.set_index(keys+['omitted_gse']);count=0;maximum=0.;balanced={};extras=0
 for file in sorted((root/'private_numeric_fixtures').glob('*.npz')):
  a=np.load(file);key=(str(a['direction']),str(a['variant']),int(a['LR_index']));effect=module.explicit(a['Z'],a['x'],a['Y'],a['blocks']);assert effect is not None
  if file.name.startswith('CASE'):balanced[key[:2]]=balanced.get(key[:2],0)+1
  else:extras+=1
  def compare(e,table,suffix=()):
   nonlocal count,maximum
   if e is None:return
   for j,name in enumerate(a['programs']):
    row=table.loc[key+(str(name),)+suffix]
    for metric in ['beta','se','lower','upper','df','p_exploratory']:
     if metric not in row.index or pd.isna(row[metric]):continue
     expected=e[metric] if metric=='df' else e[metric][j];np.testing.assert_allclose(row[metric],expected,rtol=2e-7,atol=2e-9,equal_nan=True);maximum=max(maximum,abs(row[metric]-expected))
    count+=1
  compare(effect,lookup)
  for sid in np.unique(a['gse']):
   m=a['gse']==sid;compare(module.explicit(a['Z'][m],a['x'][m],a['Y'][m],a['blocks'][m]),sx,(str(sid),));compare(module.explicit(a['Z'][~m],a['x'][~m],a['Y'][~m],a['blocks'][~m]),lx,(str(sid),))
 assert len(balanced)==8 and all(v==plan['fixture_groups_per_direction_variant'] for v in balanced.values())
 j=dict(status='PASS_FULL_FAMILY_AND_BALANCED_PLUS_FIXED_INTERFACE_NUMERIC_FIXTURES',finished_UTC=datetime.now(timezone.utc).isoformat(),inventory_sha256=sha(root/'inventory.json'),models=len(df),balanced_fixture_groups=sum(balanced.values()),extra_fixed_interface_groups=extras,numeric_component_checks=count,max_absolute_difference=maximum,all_coefficients_independently_replayed=False,independent_biological_validation=False,auditor_sha256=sha(__file__));(root/'independent_audit.json').write_text(json.dumps(j,indent=2)+'\n');print(json.dumps(j),flush=True)
if __name__=='__main__':main()
