"""Audit control construction from independent cell/cache inputs, no producer import."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,gc
import numpy as np,pandas as pd


def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()


def main():
 start=datetime.now(timezone.utc).isoformat();R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);P=R/'immune_same_name_family_source_review_v1';O=P/'fullgene_receiver_target_association_v1';F=P/'fullgene_receiver_target_readiness_v1';C=P/'source_reviewed_family_LR_cache_v1';assert not (O/'independent_controls_audit_v1.json').exists()
 repair=json.loads((R/'fullgene_receiver_target_code_v1/target_association_audit_repair_v2.json').read_text());assert sha(__file__)==repair['controls_auditor_sha256'];assert sha(O/'inventory.json')==repair['producer_inventory_sha256']
 scopes=pd.read_csv(P/'patient_block_provenance_overlay_v1/private_all1290_patient_provenance_overlay.tsv',sep='\t',keep_default_na=False).set_index('sampleID')
 models=pd.read_csv(O/'all1860_target_association_models.tsv',sep='\t',keep_default_na=False);details={j['test_id']:j for j in json.loads((O/'score_provenance.json').read_text())}
 panels={}
 scoretable=pd.read_csv(C/'all_sample_family_gene_programme_scores.tsv.gz',sep='\t')
 for r in scoretable[['program','panel']].drop_duplicates().itertuples():panels.setdefault(r.program,set()).update(r.panel.split(';'))
 mix=pd.read_csv(C/'all_sample_family_member_mixture.tsv.gz',sep='\t');mix=mix[mix['mode'].eq('marker_supported_main') & mix.included.eq(True)]
 with np.load(C/'source_reviewed_native_family_expression_cache.npz',allow_pickle=False) as c:
  mode=list(c['modes']).index('marker_supported_main');families=list(c['families']);samples=pd.Index(c['samples']);genes=pd.Index(c['genes']);fraction=c['within_major_fractions'][mode];mean_cp=c['mean_CP10K'][mode];mean_log=c['mean_log1p_CP10K'][mode];gene_measured=c['gene_measured']
 checks=0;control_values=0
 for rec in json.loads((P/'supported_LR_native_cell_panel_v1/inventory.json').read_text())['source_records']:
  sid=rec['gse_id'];assert sha(rec['selected_cache'])==rec['cache_sha256'];z=dict(np.load(rec['selected_cache'],allow_pickle=False));cell_fam=z['families'][z['family_code']];rg={g:i for i,g in enumerate(z['genes'])};ss=pd.Index(z['samples']);scope=scopes.loc[ss];cg=samples.get_indexer(ss);assert (cg>=0).all()
  for test in models[models.gse_id.eq(sid)].test_id.unique():
   file=O/'private_designs'/(test+'.npz')
   if not file.exists():continue
   saved=dict(np.load(file,allow_pickle=False));row=models[models.test_id.eq(test)].iloc[0];item=details[test];ids=ss.get_indexer(saved['sample_keys']);assert (ids>=0).all();mf=cell_fam==row.sender_family;mr=z['subCluster']==row.receiver_original_node
   source_counts=np.bincount(z['sample_code'][mf],minlength=len(ss));receiver_counts=np.bincount(z['sample_code'][mr],minlength=len(ss));fields=[]
   for k in ids:
    sender=np.flatnonzero(mf & (z['sample_code']==k));receiver=np.flatnonzero(mr & (z['sample_code']==k));L=np.mean(z['raw_counts'][sender][:,[rg[g] for g in row.ligand.split('_')]]*(10000/z['full_RNA_total'][sender])[:,None],axis=0);V=np.mean(z['raw_counts'][receiver][:,[rg[g] for g in row.receptor.split('_')]]*(10000/z['full_RNA_total'][receiver])[:,None],axis=0)
    fields.append([np.log1p(L).mean(),np.log1p(V).mean(),np.log1p(len(sender)),np.log1p(len(receiver)),np.log1p(z['full_RNA_total'][sender].mean()),np.log1p(z['full_RNA_total'][receiver].mean()),np.log1p(z['native_detected_genes'][sender].mean()),np.log1p(z['native_detected_genes'][receiver].mean())])
   a0=np.asarray(fields);assert np.allclose(a0,saved['raw_a_M0_expression_depth'],rtol=1e-8,atol=1e-9)
   sender_f=families.index(row.sender_family);receiver_f=families.index(row.receiver_family)
   # Reconstruct major totals from the independently audited family fraction.
   with np.load(F/'private_source_unit_arrays'/(sid+'.npz'),allow_pickle=False) as b:
    bn=b['n_cells'].reshape(len(b['nodes']),len(ss));fam_total=bn[b['node_family']==row.receiver_family].sum(axis=0)
   parent_major=fam_total[ids]/fraction[receiver_f,cg[ids]];rfrac=receiver_counts[ids]/parent_major
   sub=mix[mix.gse_id.eq(sid) & mix.family.eq(row.sender_family)];context_keys=set(scope.index[scope.stratum.eq(row.stratum) & scope.analysis_eligible]);present=[]
   for node,g in sub[sub.sampleID.isin(context_keys)].groupby('subCluster'):
    if g.n_cells_native.sum()>0:present.append(node)
   present.sort(key=lambda n:int(n.rsplit('_',1)[1]));reference=present[-1];assert reference==item['sender_mixture_reference'];others=[n for n in present if n!=reference]
   columns=[]
   for node in others:
    values=sub[sub.subCluster.eq(node)].set_index('sampleID').n_cells_native.reindex(saved['sample_keys'],fill_value=0).to_numpy()/source_counts[ids];columns.append(values)
   mix_values=np.column_stack(columns) if columns else np.empty((len(ids),0));a1=np.column_stack([a0,fraction[sender_f,cg[ids]],rfrac,mix_values]);assert np.allclose(a1,saved['raw_a_M1_abundance_identity'],rtol=1e-8,atol=1e-9)
   exclude=set(row.ligand.split('_')+row.receptor.split('_')+['CCL5','IL7R','MS4A1','FCER1G']);applicable=list(item['sender_programme_retained_genes'])
   expected_applicable=set(panels) if row.sender_family.startswith('T_cell_') else {'IFN_response8','cycling8','antigen_presentation8'};assert set(applicable)==expected_applicable
   programs=[]
   for label in applicable:
    remaining=sorted(panels[label]-exclude);measured=[g for g in remaining if g in genes and gene_measured[cg[ids[0]],genes.get_loc(g)]];assert set(measured)==set(item['sender_programme_retained_genes'][label].split(';'))
    good=len(measured)>=4 and len(measured)>=.6*len(remaining)
    programs.append(mean_log[sender_f,cg[ids]][:,genes.get_indexer(measured)].mean(axis=1) if good else np.full(len(ids),np.nan))
   a2=np.column_stack([a1,np.column_stack(programs)]);assert np.allclose(a2,saved['raw_a_M2_sender_shared_programmes'],rtol=1e-8,atol=1e-9,equal_nan=True)
   # Names and ordering are verified separately from numeric values.
   assert list(saved['names_M1_abundance_identity'][-len(others):])==['sender_mix_'+n for n in others] if others else True
   assert list(saved['names_M2_sender_shared_programmes'][-len(applicable):])==['sender_program_'+n for n in applicable]
   checks+=1;control_values+=a0.size+a1.size+a2.size
  del z;gc.collect()
 assert checks==325
 audit=dict(status='PASS_ALL325_NATIVE_CELL_LR_COUNTS_RNA_DEPTH_MAJOR_DENOMINATORS_SENDER_MIXTURE_AND_SHARED_PROGRAMME_CONTROLS',started_UTC=start,finished_UTC=datetime.now(timezone.utc).isoformat(),auditor_sha256=sha(__file__),producer_inventory_sha256=sha(O/'inventory.json'),cases_checked=checks,raw_control_values_checked=control_values,method='Native cell cache L/R/RNA/nGenes; accepted independent major fractions and member TSV; old1933gene family programme means with exact exclusions; no producer import',independent_biological_validation=False,scientific_release=False,axis_validation_complete=False)
 (O/'independent_controls_audit_v1.json').write_text(json.dumps(audit,indent=2)+'\n');print(json.dumps(audit),flush=True)


if __name__=='__main__':main()
