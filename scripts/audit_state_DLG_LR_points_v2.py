"""Independent all-point/LR-mask/BH audit plus prespecified resampling diagnostics."""
from pathlib import Path
from datetime import datetime,timezone
from collections import Counter
import os,json,hashlib
import numpy as np,pandas as pd
from scipy.stats import rankdata
from scipy.linalg import lstsq
V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();P=V/'state_after_DIALOGUE_LR_response_20261016_v1';iv=json.loads((P/'inventory.json').read_text());assert iv['estimable_tests']==4012
for f,h in iv['files'].items():assert sha(P/f)==h
T=pd.read_csv(P/'all_estimable_exploratory_LR_receiver_associations.tsv',sep='\t');K=['context','source','target','MCP','ligand','receptor'];assert not T.duplicated(K).any();N=len(T);order=np.argsort(T.permutation_p.to_numpy(),kind='stable');q=np.empty(N);q[order]=np.minimum(1,np.minimum.accumulate((T.permutation_p.iloc[order].to_numpy()*N/np.arange(1,N+1))[::-1])[::-1]);assert np.allclose(q,T.BH_q,atol=1e-12);assert not T.ligand.eq('C1QA').mul(T.receptor.eq('CD93')).any();assert not (T.ligand.eq('INHBA')&T.receptor.eq('ACTR2')).any();assert (T.valid_bootstraps<=2000).all();assert T.loc[T.valid_bootstraps<1900,['CI_low','CI_high']].isna().all().all();assert T.loc[T.valid_bootstraps>=1900,['CI_low','CI_high']].notna().all().all()
S=V/'state_after_DIALOGUE_LR_scope_detail_20261016_v2';si=json.loads((S/'inventory.json').read_text());assert si['parent_inventory_sha256']==sha(P/'inventory.json') and si['all_effects_p_CI_BH_byte_identical'];assert sha(S/'all_estimable_exploratory_LR_receiver_associations.tsv')==iv['files']['all_estimable_exploratory_LR_receiver_associations.tsv'];scope=pd.read_csv(S/'all258944_hypothesis_scope_counts.tsv',sep='\t');assert scope.records.sum()==258944 and scope.loc[scope.status.eq('ESTIMABLE_EXPLORATORY_FULL_BLOCK_ASSOCIATION'),'records'].sum()==4012
G=V/'state_after_DIALOGUE_gene_components_20261016_v1';gi=json.loads((G/'inventory.json').read_text());assert sha(G/'inventory.json')==iv['component_inventory_sha256'];gene=G/'private_all_block_gene_components.tsv.gz';assert sha(gene)==gi['files'][gene.name];comp=pd.read_csv(gene,sep='\t');E=R/'immune_same_name_family_source_review_v1/CM_edge_LR_extraction_v2';ei=json.loads((E/'inventory.json').read_text());groups=pd.read_csv(E/'all68_directed_group_summary.tsv',sep='\t');Q=V/'all_source_native_UCell_unit_sensitivity_20261014_v1';qi=json.loads((Q/'inventory.json').read_text());qc={}
for info in qi['private_outputs']:
 p=Q/info['relative_file'];assert sha(p)==info['sha256'];z=np.load(p,allow_pickle=False)
 for k,r,d in zip(z['unit_index'],z['mean_full_RNA'],z['mean_full_detected_genes']):qc[int(k)]=(float(r),float(d))
def rho(x,y,C=None):
 a=rankdata(x);b=rankdata(y);a-=a.mean();b-=b.mean();na=np.linalg.norm(a);nb=np.linalg.norm(b)
 if C is not None:a-=C@lstsq(C,a,lapack_driver='gelsd')[0];b-=C@lstsq(C,b,lapack_driver='gelsd')[0]
 if np.linalg.norm(a)<=1e-12*max(na,1) or np.linalg.norm(b)<=1e-12*max(nb,1):return np.nan
 return float(a@b/(np.linalg.norm(a)*np.linalg.norm(b)))
O=V/'state_DLG_LR_point_audit_20261016_v2';assert not O.exists();O.mkdir(mode=0o700);errs=[];helderrs=[];diag=[];checked=0;rng=np.random.default_rng(2026101606)
for ctx,ct in T.groupby('context',sort=True):
 root=V/'state_after_DIALOGUE_inputs_20261016_v3'/f'context_{ctx:04d}';ci=json.loads((root/'private_input_inventory.json').read_text());parts=pd.read_csv(root/'private_block_library_partition.tsv',sep='\t');cache={}
 for (source,target,mcp),st in ct.groupby(['source','target','MCP'],sort=True):
  rev=source==ci['endpoints'][1];typ='typeA' if rev else 'typeB';g=groups[groups.stratum.eq(ci['stratum'])&groups.source.eq(source)&groups.target.eq(target)];assert len(g)==1;p=E/'private_group_arrays'/(g.iloc[0].group_id+'.npz');assert sha(p)==ei['outputs'][str(p.relative_to(E))];z=np.load(p,allow_pickle=False);ii=pd.Index(z['samples'].astype(str)).get_indexer(parts.sample_uid.astype(str));assert (ii>=0).all() and list(z['patient_blocks'][ii])==list(parts.effective_resampling_block);qs=np.array([qc[int(k)] for k in parts['unit_index_B' if rev else 'unit_index_A']]);qt=np.array([qc[int(k)] for k in parts['unit_index_A' if rev else 'unit_index_B']]);cv=np.column_stack([np.log1p(z['source_n'][ii]),np.log1p(z['target_n'][ii]),z['source_within_major'][ii],z['target_within_major'][ii],np.log1p(qs[:,0]),np.log1p(qt[:,0]),np.log1p(qs[:,1]),np.log1p(qt[:,1])]);lrindex={(str(l),str(r)):i for i,(l,r) in enumerate(zip(z['ligands'],z['receptors']))};cc=comp[(comp.context==ctx)&comp.type.eq(typ)&comp.program.eq(mcp)];w=cc.groupby('gene').weight.first();B=cc.pivot(index='sample',columns='gene',values='contribution').reindex(parts.effective_resampling_block);assert B.notna().all().all();weights=w.reindex(B.columns).to_numpy();full=B.to_numpy()@weights
  for row in st.itertuples():
   j=lrindex[(row.ligand,row.receptor)];values=z['scores'][ii,j];mask=np.isfinite(values[:,0]);assert mask.sum()==row.blocks and mask.sum()/len(parts)>=.8 and z['native_measured'][j];remove=set(row.ligand.split('_')+row.receptor.split('_'));keep=~B.columns.isin(remove);assert keep.sum()==row.retained_programme_genes>=20;ret=np.abs(weights[keep]).sum()/np.abs(weights).sum();assert abs(ret-row.retained_L1_fraction)<1e-12 and ret>=.8;assert ';'.join(sorted(remove&set(B.columns)))==(row.masked_genes if isinstance(row.masked_genes,str) else '');y=(B.to_numpy()[:,keep]@weights[keep])[mask];x=np.log1p(values[mask,0]);rawcov=cv[mask];sd=rawcov.std(0);assert (sd>1e-12).all();D=np.column_stack([np.ones(len(x)),(rawcov-rawcov.mean(0))/sd]);assert np.linalg.matrix_rank(D)==9 and len(x)-10>=10
   calc=[rho(x,y),rho(x,y,D),rho(x,full[mask],D),rho(np.log1p(values[mask,1]),y,D),rho(-values[mask,2],y,D),rho(-values[mask,3],y,D)];saved=[row.raw_masked_rho,row.adjusted_masked_rho,row.adjusted_unmasked_rho,row.lr_means_adjusted_rho,row.negative_specificity_rank_adjusted_rho,row.negative_magnitude_rank_adjusted_rho];assert np.allclose(calc,saved,atol=1e-8,equal_nan=True);errs.extend([abs(a-b) for a,b in zip(calc,saved) if np.isfinite(a) and np.isfinite(b)]);tr=parts.partition.eq('discovery').to_numpy()[mask];te=parts.partition.eq('holdout').to_numpy()[mask];assert tr.sum()==row.discovery_observed_blocks and te.sum()==row.internal_holdout_observed_blocks
   if tr.sum()>=20 and te.sum()>=10 and np.linalg.matrix_rank(D[tr])==9:
    rx=x[te]-D[te]@lstsq(D[tr],x[tr],lapack_driver='gelsd')[0];ry=y[te]-D[te]@lstsq(D[tr],y[tr],lapack_driver='gelsd')[0];hh=rho(rx,ry);assert np.isclose(hh,row.internal_holdout_rho,atol=1e-8,equal_nan=True)
    if np.isfinite(hh):helderrs.append(abs(hh-row.internal_holdout_rho))
   else:assert pd.isna(row.internal_holdout_rho)
   # Three prespecified actual cases: robust positive, robust negative, fragile interval.
   select=(ctx==716 and source=='Blood_vascular_endothelial_cell_RAMP2_10' and row.ligand=='LGALS1' and row.receptor=='PTPRC' and mcp=='MCP1') or (ctx==1060 and source=='Macrophage_IFI6_18' and row.ligand=='MIF' and row.receptor=='CD74_CXCR4' and mcp=='MCP2') or (ctx==1611 and source=='Epithelial_cell_KRT18_52' and row.ligand=='ANXA2' and row.receptor=='TLR2' and mcp=='MCP2')
   if select:
    bs=[];bad=0;constant=Counter();names=['source_log_cells','target_log_cells','source_within_major','target_within_major','source_log_RNA','target_log_RNA','source_log_detected','target_log_detected']
    for k in range(2000):
     ix=rng.integers(0,len(x),len(x))
     if np.linalg.matrix_rank(D[ix])<9:
      bad+=1
      for nm in np.array(names)[rawcov[ix].std(0)<=1e-12]:constant[str(nm)]+=1
      continue
     a=rho(x[ix],y[ix],D[ix])
     if np.isfinite(a):bs.append(a)
    fullcounts=[]
    for k,nm in enumerate(names):
     vals,counts=np.unique(rawcov[:,k],return_counts=True);fullcounts.append(dict(covariate=nm,unique_values=len(vals),dominant_level_blocks=int(counts.max()),other_blocks=int(len(x)-counts.max())))
    diag.append(dict(context=int(ctx),source=source,target=target,MCP=mcp,ligand=row.ligand,receptor=row.receptor,original_valid_bootstraps=int(row.valid_bootstraps),independent_seed=2026101606,independent_valid_bootstraps=len(bs),independent_rank_deficient_bootstraps=bad,independent_bootstrap_CI=np.quantile(bs,[.025,.975]).tolist() if len(bs)>=1900 else None,independent_percentiles_descriptive_only=np.quantile(bs,[.025,.975]).tolist(),constant_covariate_counts_in_failed_resamples=dict(constant),covariate_level_coverage=fullcounts,not_replacement_of_primary_bootstrap_or_P=True))
   checked+=1
 (O/'progress.json').write_text(json.dumps(dict(context=int(ctx),point_records_checked=checked,observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n')
assert checked==4012 and len(diag)==3;pd.DataFrame([dict(context=d['context'],ligand=d['ligand'],receptor=d['receptor'],original_valid_bootstraps=d['original_valid_bootstraps'],independent_valid_bootstraps=d['independent_valid_bootstraps'],independent_rank_deficient_bootstraps=d['independent_rank_deficient_bootstraps']) for d in diag]).to_csv(O/'three_actual_resampling_cases.tsv',sep='\t',index=False);(O/'resampling_design_diagnostics.json').write_text(json.dumps(diag,indent=2)+'\n');b=T[T.BH_q<.05];out=dict(status='PASS_ALL4012_POINT_MASK_SENSITIVITY_HELDOUT_BH_AND_CI_GATES_WITH_THREE_INDEPENDENT_RESAMPLING_DIAGNOSTICS',observed_UTC=datetime.now(timezone.utc).isoformat(),auditor_sha256=sha(__file__),producer_inventory_sha256=sha(P/'inventory.json'),scope_detail_inventory_sha256=sha(S/'inventory.json'),point_records_checked=checked,max_point_error=max(errs),max_heldout_error=max(helderrs),BH_records=len(b),BH_with_CI_excluding_zero=int(((b.CI_low*b.CI_high)>0).sum()),primary_permutations_independently_replayed_all=False,all_primary_bootstrap_endpoint_values_independently_replayed=False,independent_resampling_cases=3,files={p.name:sha(p) for p in O.iterdir() if p.is_file()},formal_states=False,causal_validation=False);(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
