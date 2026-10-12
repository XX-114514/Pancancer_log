"""Targeted independent full-X refits of all three full-gene-only state-control BH records.
Same-data numerical verification, never an independent biological validation.
"""
import os,json,importlib.util,hashlib
from pathlib import Path
from datetime import datetime,timezone
from collections import defaultdict
import numpy as np,pandas as pd
from scipy.linalg import block_diag
from scipy.sparse import load_npz

def standardize(v):
 v=np.asarray(v,float);sd=np.std(v,axis=0);return (v-v.mean(axis=0))/np.where(sd>1e-12,sd,1)

def main():
 R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);P=R/'immune_same_name_family_source_review_v1';N=P/'major_pair_native_assembly_v1';T=P/'major_receiver_TF_preflight_v3';D=P/'major_LR_TF_design_preflight_v1';S=P/'macrophage_T_state_sensitivity_v1';plan=json.loads(Path(os.environ['PANCANCER_MAC_T_POOL_PLAN']).read_text());audit=Path(os.environ['PANCANCER_EXPLICIT_POOL_AUDITOR']);spec=importlib.util.spec_from_file_location('explicit',audit);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);table=pd.read_csv(S/'all27744_pooled_LR_models.tsv.gz',sep='\t');selected=table[(table.variant=='full_gene_sensitivity')&(table.joint_BH_exploratory<.05)];assert len(selected)==3;component=pd.read_csv(S/'all_pooled_TF_components.tsv.gz',sep='\t').set_index(['direction','variant','LR_index','TF']);scope=pd.read_csv(D/'all_context_model_readiness.tsv',sep='\t').set_index('context_id');counts=pd.read_csv(AR/'20261010_expression_annotation_downstream_v1/inputs/sample_major_subcluster_counts.tsv',sep='\t').groupby('sampleID').n_cells.sum();cache={}
 for cid in plan['contexts']:
  row=scope.loc[cid];a=np.load(D/'private_context_arrays'/(cid+'.npz'));b=np.load(T/'private_context_arrays'/(cid+'.npz'));g=np.load(N/'private_group_arrays'/(row.group_id+'.npz'));m=g['corrected_material_for_next_inference']==row.material;ns=g['source_n'][m];nt=g['target_n'][m];samples=a['samples'];assert np.array_equal(samples,b['samples']) and np.array_equal(samples,g['samples'][m]);base=np.column_stack((ns/counts.loc[samples],nt/counts.loc[samples],np.log1p(ns),np.log1p(nt),g['source_mean_detected_nonLR_genes'][m],g['target_mean_detected_nonLR_genes'][m]));cache[cid]=(row,a,b,g,m,base)
 state={};sc={}
 for cid,(row,a,b,g,m,base) in cache.items():
  sid=row.gse_id
  if sid not in sc:
   units=pd.read_csv(N/'private_source_major_arrays'/(sid+'_units.tsv'),sep='\t');genes=np.load(N/'private_source_major_arrays'/(sid+'_axes.npz'))['genes'];scores=np.load(T/'private_source_TF_arrays'/(sid+'.npz'));lookup={(r.sampleID,r.majorCluster):r.unit_row for r in units.itertuples()};sc[sid]=(lookup,{str(g):i for i,g in enumerate(genes)},load_npz(N/'private_source_major_arrays'/(sid+'_meanlog.npz')),scores['primary_scores'][:,np.flatnonzero(scores['TFs']=='HSF1')[0]])
  lookup,genes,mean,h= sc[sid];si=[lookup[(v,row.source)] for v in a['samples']];ri=[lookup[(v,row.target)] for v in a['samples']];state[cid]=(h[si],genes,mean[ri])
 checks=0;maximum=0.;records=[]
 for r in selected.itertuples():
  for variant in ['abundance_complexity','subtype_mixture','full_gene_sensitivity']:
   eligibility=defaultdict(list)
   for cid in plan['direction_contexts'][r.direction]:
    row,a,b,g,m,base=cache[cid]
    if a['mixture_design_status' if variant=='subtype_mixture' else 'primary_design_status'][r.LR_index]!=0:continue
    for tf in a['candidate_TF_index'][a['candidate_LR_index']==r.LR_index]:eligibility[int(tf)].append(cid)
   patterns=defaultdict(list)
   for tf,ids in eligibility.items():
    if len({cache[c][0].gse_id for c in ids})>=2:patterns[tuple(ids)].append(tf)
   for ids,TFs in patterns.items():
    ZZ=[];xx=[];yy=[];bb=[]
    for cid in ids:
     row,a,b,g,m,base=cache[cid];hsf,genes,mean=state[cid];own=mean[:,[genes[gene] for gene in r.ligand.split('_')]].toarray().min(axis=1);Z=np.column_stack((np.ones(len(a['samples'])),standardize(np.column_stack((base,g['receptor_meanlog'][m,r.LR_index],hsf,own)))))
     if variant=='subtype_mixture':Z=np.column_stack((Z,standardize(a['source_target_mixture_PCs'])))
     ZZ.append(Z);xx.append(standardize(g['ligand_meanlog'][m,r.LR_index]));yy.append(standardize(b['sensitivity_TF_scores' if variant=='full_gene_sensitivity' else 'primary_TF_scores'][:,TFs]));bb.append(b['blocks'])
    effect=module.explicit(block_diag(*ZZ),np.concatenate(xx),np.vstack(yy),np.concatenate(bb))
    if effect is None:continue
    for j,tf in enumerate(TFs):
     name=str(cache[ids[0]][2]['TFs'][tf]);actual=component.loc[(r.direction,variant,r.LR_index,name)]
     for metric in ['beta','se','df','lower','upper','p_exploratory']:
      if metric=='p_exploratory' and pd.isna(actual[metric]):continue
      expected=effect[metric] if metric=='df' else effect[metric][j];np.testing.assert_allclose(actual[metric],expected,rtol=2e-7,atol=2e-9,equal_nan=True);maximum=max(maximum,abs(actual[metric]-expected))
     checks+=1
   records.append(dict(direction=r.direction,ligand=r.ligand,receptor=r.receptor,variant=variant,patterns=len(patterns)))
 result=dict(status='PASS_ALL3_FULLGENE_ONLY_STATE_RECORDS_ALL_TFS_THREE_VARIANTS',finished_UTC=datetime.now(timezone.utc).isoformat(),component_refits=checks,maximum_absolute_difference=maximum,records=records,numeric_only_not_independent_biological_validation=True,code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest());(S/'all3_sensitivity_positive_fullX_numeric_audit_v1.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='records'}))
if __name__=='__main__':main()
