"""Resolve current exact memberships to already accepted native gene moments.

Profiles are referenced rather than copied to a new expression archive. Source
gene axes and whole native RNA denominators retain their original meaning.
"""
import hashlib,json
from pathlib import Path
import numpy as np

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

class AcceptedNativeProfiles:
 def __init__(self,working_acceptance,complete,complete_audit,changed,changed_audit):
  self.complete=Path(complete);self.changed=Path(changed);self.cache=None
  for producer,auditor in [(self.complete,Path(complete_audit)),(self.changed,Path(changed_audit))]:
   for stage in [producer,auditor]:
    if load(stage/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']!=0:raise ValueError('Successful profile production and audit required')
   if load(auditor/'provenance/terminal_v1.json')['producer_terminal_SHA256']!=sha(producer/'provenance/terminal_v1.json'):raise ValueError('Profile audit binding differs')
  accepted=load(working_acceptance)
  if sha(accepted['producer_manifest']['path'])!=accepted['producer_manifest']['SHA256']:raise ValueError('Working manifest changed')
  self.working=load(accepted['producer_manifest']['path'])
  for k in ['index','group_catalog']:
   if sha(self.working[k]['path'])!=self.working[k]['SHA256']:raise ValueError('Working profile axis changed')
  self.groups=load(self.working['group_catalog']['path'])['groups']
  ep=self.complete/'outputs/TME_full_native_evidence_index_v1.json';sp=self.changed/'outputs/TCR_changed_membership_native_index_v1.json'
  if sha(ep)!=load(self.complete/'provenance/terminal_v1.json')['index_SHA256'] or sha(sp)!=load(self.changed/'provenance/terminal_v1.json')['index_SHA256']:raise ValueError('Native profile index changed')
  self.ei={r['source']:r for r in load(ep)};self.si={r['source']:r for r in load(sp)}
 def source(self,source):
  if self.cache and self.cache['source']==source:return self.cache
  self.cache=None;entry=self.ei[source];mp=self.complete/'outputs'/entry['metadata']
  if sha(mp)!=entry['metadata_SHA256']:raise ValueError('Complete native metadata changed')
  e=load(mp);p=e['group_statistics']
  if sha(p['path'])!=p['SHA256']:raise ValueError('Complete native statistics changed')
  with np.load(p['path'],allow_pickle=False) as f:es={k:f[k] for k in ['n','cp_sum','log_sum','detected','native_columns']}
  genes=e['native_measured_gene_names'];cols=np.asarray(e['native_columns']);np.testing.assert_array_equal(es['native_columns'],cols)
  available={g['stable_group_sha256']:(es,g['source_annotation_group_code'],{'statistics_path':p['path'],'statistics_SHA256':p['SHA256'],'profile_row':g['source_annotation_group_code'],'metadata_path':str(mp),'metadata_SHA256':sha(mp),'stable_membership_SHA256':g['stable_group_sha256']}) for g in e['groups']}
  if source in self.si:
   sr=self.si[source];sm=self.changed/'outputs'/sr['metadata']
   if sha(sm)!=sr['metadata_SHA256']:raise ValueError('Changed native metadata changed')
   s=load(sm);p=s['statistics']
   if sha(p['path'])!=p['SHA256']:raise ValueError('Changed native statistics changed')
   with np.load(p['path'],allow_pickle=False) as f:ss={k:f[k] for k in ['n','cp_sum','log_sum','detected','native_columns']}
   assert s['native_measured_genes']==genes;np.testing.assert_array_equal(ss['native_columns'],cols)
   for g in s['groups']:
    key=g['stable_membership_SHA256'];j=g['measurement_subset_code']
    if key in available:
     old,k,_=available[key]
     for field in ['n','detected']:np.testing.assert_array_equal(old[field][k],ss[field][j])
     for field in ['cp_sum','log_sum']:np.testing.assert_allclose(old[field][k],ss[field][j],rtol=2e-10,atol=1e-7)
    available[key]=(ss,j,{'statistics_path':p['path'],'statistics_SHA256':p['SHA256'],'profile_row':j,'metadata_path':str(sm),'metadata_SHA256':sha(sm),'stable_membership_SHA256':key})
  groups=sorted([g for g in self.groups if g['source']==source],key=lambda g:g['source_annotation_group_code']);assert [g['source_annotation_group_code'] for g in groups]==list(range(len(groups)))
  stat={'n':np.zeros(len(groups),np.int64),**{k:np.zeros((len(groups),len(genes)),np.float64) for k in ['cp_sum','log_sum','detected']}};references=[]
  for j,g in enumerate(groups):
   if g['stable_group_sha256'] not in available:raise ValueError('Current exact group has no accepted own native profile')
   a,k,reference=available[g['stable_group_sha256']]
   if int(a['n'][k])!=g['n_scope_group_rows']:raise ValueError('Own native profile membership count differs')
   for field in stat:stat[field][j]=a[field][k]
   references.append(reference)
  tn=np.asarray([bool(g['major_distribution']) and set(g['major_distribution']).issubset({'T','NK','T_NK'}) for g in groups],bool)
  self.cache={'source':source,'genes':genes,'native_columns':cols,'groups':groups,'stat':stat,'profile_references':references,'T_NK_group_mask':tn,'native_metadata_path':str(mp),'native_metadata_SHA256':sha(mp)}
  return self.cache
