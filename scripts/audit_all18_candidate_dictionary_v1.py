"""Independent signature reconstruction and scope acceptance for measurement only."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib
import numpy as np,pandas as pd

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def frame(p):
 with np.load(p,allow_pickle=True) as z:return pd.DataFrame(z['data'],index=z['index'],columns=z['columns'])
def selected(v,genes):return [str(genes[j]) for j in sorted(range(len(v)),key=lambda j:(-float(v[j]),j)) if v[j]>0][:50]
R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);D=V/'all_lineage_state_candidate_dictionary_20261014_v4';dp=D/'frozen_programme_dictionary.json';d=json.loads(dp.read_text());iv=json.loads((D/'inventory.json').read_text());assert iv['dictionary_sha256']==sha(dp) and d['schema']=='all18_candidate_state_dictionary_v4';assert d['formal_states'] is False and len(d['input_scopes'])==24
plan=json.loads(Path('/tmp/all_lineage_state_candidate_dictionary_plan_v4.json').read_text());assert sha(Path('/tmp/all_lineage_state_candidate_dictionary_plan_v4.json'))==d['plan_sha256'];accepted=R/'all_lineage_state_native_input_audit_20261014_v2/accepted_input_inventory.json';assert sha(accepted)==d['accepted_input_inventory_sha256'];majors=sorted(r['major'] for r in json.loads(accepted.read_text())['files']);lookup={r['programme_id']:r for r in d['programmes']};assert len(lookup)==len(d['programmes']);expected=set();checks=[]
for i,m in enumerate(majors):
 diag=R/'all_lineage_cNMF_diagnostics_20261014_v1/major_00' if i==0 else V/'all_lineage_cNMF_diagnostics_20261014_v2'/f'major_{i:02}';tab=pd.read_csv(diag/'all45_candidates.tsv',sep='\t');assert len(tab)==45;base=V/'all_lineage_cNMF_state_contrast_20261014_v2'/f'major_{i:02}'/'cnmf/cnmf_tmp';ci=json.loads((diag/'inventory.json').read_text())
 for k in [10,15,20]:
  cp=base/f'cnmf.gene_spectra_score.k_{k}.dt_0_5.df.npz';tp=base/f'cnmf.gene_spectra_tpm.k_{k}.dt_0_5.df.npz'
  for p in [cp,tp]:assert sha(p)==ci['input_output_sha256'][str(p.relative_to(base.parent.parent))]
  c=frame(cp);t=frame(tp)
  for rr in tab[tab.K.eq(k)].itertuples():
   expected.add(rr.candidate_id);r=lookup[rr.candidate_id];assert r['allowed_majors']==[m] and r['method']=='cNMF';assert r['primary_genes']==selected(c.loc[rr.component].to_numpy(),c.columns);assert r['TPM_sensitivity_genes']==selected(t.loc[rr.component].to_numpy(),t.columns);assert r['algorithm_stable_candidate']==bool(rr.repeat_and_crossK_candidate);assert r['native_training_supported_cells']==rr.informative_gene_support_cells;checks.append(dict(programme_id=rr.candidate_id,method='cNMF',primary_genes=len(r['primary_genes']),TPM_genes=len(r['TPM_sensitivity_genes']),signature_reconstructed=True))
sp=json.loads(Path('/tmp/all_lineage_spectra_diagnostic_plan_v2.json').read_text())
for block,mm in sp['blocks'].items():
 diag=R/'all_lineage_Spectra_diagnostics_20261014_v1'/block if block=='lymphoid_T' else V/'all_lineage_Spectra_diagnostics_20261014_v2'/block;tab=pd.read_csv(diag/'all_candidate_diagnostics.tsv',sep='\t',keep_default_na=False);si=json.loads((diag/'inventory.json').read_text());root=V/'all_lineage_Spectra_state_discovery_20261014_v2';weights=[]
 for seed in sp['seeds']:
  p=root/block/f'seed{seed}'/'factors_and_discovery_scores.npz';assert sha(p)==si['input_seed_sha256'][str(p.relative_to(root))]
  with np.load(p,allow_pickle=False) as z:weights.append(z['factors'].copy());genes=z['genes'].copy()
 for rr in tab.itertuples():
  if not rr.stable_context_candidate:continue
  expected.add(rr.candidate_id);r=lookup[rr.candidate_id];v=[np.asarray(weights[j][ix],float) for j,ix in enumerate([rr.anchor_factor,rr.seed2_factor,rr.seed3_factor])];v=[x/np.sqrt(np.sum(x*x)) for x in v];middle=np.sort(np.stack(v),axis=0)[1];sig=[str(genes[j]) for j in sorted(range(len(genes)),key=lambda j:(-float(middle[j]),j))[:50]];assert r['primary_genes']==sig and r['allowed_majors']==(mm if rr.scope=='block_global' else [rr.scope]);assert r['stability_main_candidate']==bool(rr.stable_shared_or_lineage_candidate);checks.append(dict(programme_id=rr.candidate_id,method='Spectra',primary_genes=len(sig),TPM_genes=0,signature_reconstructed=True))
assert expected==set(lookup) and sum(r['method']=='cNMF' for r in lookup.values())==810
owners={'epithelial_RNA':{'Epithelial cell'},'stromal_RNA':{'Fibroblast','Mesenchymal cell'},'mural_RNA':{'Mural cell','Mesenchymal cell'},'neural_RNA':{'Glial cell'}}
for r in lookup.values():
 assert 0<len(r['primary_genes'])<=50 and len(set(r['primary_genes']))==len(r['primary_genes']);assert not r['formal_state_labels'] and not r['cross_method_agreement_assumed'];gs=set(r['primary_genes'])
 for m in r['allowed_majors']:
  foreign=[k for k in owners if len(gs.intersection(plan['diagnostic_panels'][k]))>=3 and m not in owners[k]];own=[k for k in owners if len(gs.intersection(plan['diagnostic_panels'][k]))>=3 and m in owners[k]];assert foreign==r['lineage_RNA_review_by_major'][m] and own==r['own_lineage_identity_panels_by_major'][m]
 assert r['foreign_RNA_sensitive']==any(r['lineage_RNA_review_by_major'].values())
for name,h in iv['public_output_sha256'].items():assert sha(D/name)==h
O=V/'all18_candidate_dictionary_acceptance_20261014_v1';assert not O.exists();O.mkdir();pd.DataFrame(checks).to_csv(O/'all_signature_reconstruction_checks.tsv',sep='\t',index=False);out=dict(status='ACCEPTED_ALL18_DICTIONARY_FOR_NATIVE_PROJECTION_NOT_FORMAL_STATES',observed_UTC=datetime.now(timezone.utc).isoformat(),audit_code_sha256=sha(__file__),dictionary_sha256=sha(dp),programmes=len(lookup),cNMF_programmes=810,all_signatures_independently_reconstructed=True,all18_original_major_scope=True,lineage_context_flags_checked=True,biological_state_review_complete=False,formal_states=False,checks_sha256=sha(O/'all_signature_reconstruction_checks.tsv'));(O/'accepted_dictionary_inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
