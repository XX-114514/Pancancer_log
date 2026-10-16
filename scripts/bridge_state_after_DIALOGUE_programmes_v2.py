"""Descriptive bridge for every fixed MCP slot; never select on heldout/LR/clinical p."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,time
import numpy as np,pandas as pd
from scipy.stats import rankdata,spearmanr
V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);O=V/'state_after_DIALOGUE_programme_bridge_20261016_v2';assert not O.exists();O.mkdir(mode=0o700)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
plan_path=Path('/tmp/state_after_DIALOGUE_bridge_plan_v2.json');plan=json.loads(plan_path.read_text());assert sha(__file__)==plan['producer_sha256'];P=V/'state_after_DIALOGUE_frozen_projection_20261016_v1';start=time.monotonic()
while not (P/'inventory.json').exists():
 e=json.loads((V/'state_after_DIALOGUE_projection_chain_20261016_v1.execution.json').read_text());p=Path('/proc',str(e['child_pid']),'cmdline');assert p.exists() and e['producer'].encode() in p.read_bytes(),e;assert time.monotonic()-start<32*3600;time.sleep(20)
iv=json.loads((P/'inventory.json').read_text())
for f,h in iv['files'].items():assert sha(P/f)==h,f
D=V/'all_lineage_state_candidate_dictionary_20261014_v4/frozen_programme_dictionary.json';assert sha(D)==plan['dictionary_sha256'];dictionary=json.loads(D.read_text())['programmes'];assert len(dictionary)==814
C=V/'original_CM_programme_pair_associations_20261015_v2/all_context_coverage_reasons.tsv';assert sha(C)==plan['context_coverage_sha256'];contexts=pd.read_csv(C,sep='\t',keep_default_na=False);contexts=contexts[contexts.status.eq('ESTIMABLE')];assert len(contexts)==14
A=V/'all_source_native_UCell_audit_aggregation_20261014_v2';assert sha(A/'inventory.json')==plan['ucell_audit_sha256'];ai=json.loads((A/'inventory.json').read_text());F=R/'all_CM_state_native_input_preflight_20261014_v1';fi=json.loads((F/'inventory.json').read_text());U=F/'private_all_frozen_sample_node_units.tsv.gz';assert sha(U)==fi['private_unit_table_sha256'];units=pd.read_csv(U,sep='\t',keep_default_na=False).set_index('unit_index',drop=False);major_map=units.groupby('subCluster').major.agg(lambda s:sorted(set(s)));assert major_map.map(len).eq(1).all();major_map=major_map.map(lambda x:x[0]);aggregates={}
for info in ai['private_aggregate_files']:
 p=A/info['relative_file'];assert sha(p)==info['sha256'];z=dict(np.load(p,allow_pickle=False));aggregates[info['major']]=(z,{int(u):i for i,u in enumerate(z['unit_index'])},{str(v):i for i,v in enumerate(z['programme_ids'])})
slotpath=P/'public_slot_scope.tsv';scope=pd.read_csv(slotpath,sep='\t');gp=P/'public_fixed_gene_weights.tsv';genes=pd.read_csv(gp,sep='\t') if gp.exists() else pd.DataFrame();sp=P/'private_block_mean_scores.tsv';scores=pd.read_csv(sp,sep='\t') if sp.exists() else pd.DataFrame();slotrows=[];overlap=[];correlations=[];maxrhoerror=0.;checked=0
for ctx in contexts.itertuples():
 root=V/'state_after_DIALOGUE_inputs_20261016_v3'/f'context_{ctx.context:04d}';ci=json.loads((root/'private_input_inventory.json').read_text());partpath=root/'private_block_library_partition.tsv';assert sha(partpath)==ci['files'][partpath.name];parts=pd.read_csv(partpath,sep='\t',keep_default_na=False)
 for side in ['A','B']:
  typ='type'+side;endpoint=getattr(ctx,'endpoint_'+side);major=major_map[endpoint];ag,unitix,pix=aggregates[major];refs=[p for p in dictionary if major in p['allowed_majors']]
  for mcp in ['MCP1','MCP2']:
   q=scope[(scope.context==ctx.context)&scope.type.eq(typ)&scope.program.eq(mcp)];assert len(q)<=1;ok=len(q)==1 and q.iloc[0].status=='FROZEN_RECONSTRUCTION_PASS';slotrows.append(dict(context=ctx.context,gse_id=ctx.gse_id,CMs=ctx.CMs,endpoint=endpoint,major=major,type=typ,MCP=mcp,reconstruction_pass=ok,reference_programmes=len(refs),formal_state=False))
   if not ok:continue
   gg=genes[(genes.context==ctx.context)&genes.type.eq(typ)&genes.program.eq(mcp)];assert len(gg)>0 and not gg.gene.duplicated().any();up=set(gg.loc[gg.weight>0,'gene']);down=set(gg.loc[gg.weight<0,'gene']);allg=up|down
   for ref in refs:
    rid=ref['programme_id'];rg=set(ref['primary_genes']);common=allg&rg;upcommon=up&rg;downcommon=down&rg;assert len(common)==len(upcommon)+len(downcommon);weights=dict(zip(gg.gene,gg.weight));den=float(np.abs(gg.weight).sum());overlap.append(dict(context=ctx.context,gse_id=ctx.gse_id,CMs=ctx.CMs,type=typ,MCP=mcp,endpoint=endpoint,reference_programme=rid,reference_method=ref['method'],reference_stable=ref['algorithm_stable_candidate'],reference_descriptors=';'.join(ref['descriptive_programmes']),reference_basic_expression=ref['basic_expression_review_flag'],reference_foreign_RNA_sensitive=ref['foreign_RNA_sensitive'],MCP_genes=len(allg),reference_genes=len(rg),up_overlap=len(upcommon),down_overlap=len(downcommon),jaccard=len(common)/len(allg|rg),absolute_weight_fraction_in_reference=sum(abs(weights[g]) for g in common)/den,overlap_genes=';'.join(sorted(common))))
    col=pix[rid]
    for partition in ['discovery','holdout']:
     pp=parts[parts.partition.eq(partition)];v=scores[(scores.context==ctx.context)&scores.type.eq(typ)&scores.program.eq(mcp)&scores.partition.eq(partition)][['sample','score']];j=pp.merge(v,left_on='effective_resampling_block',right_on='sample',validate='one_to_one');assert len(j)==len(pp)==len(v);rows=np.array([unitix[int(k)] for k in j['unit_index_'+side]],int);uv=ag['programme_mean'][rows,col];valid=(ag['valid_score_cells'][rows,col]>=20)&np.isfinite(uv)&np.isfinite(j.score.to_numpy());x=j.score.to_numpy()[valid];y=uv[valid];n=len(x);rho=np.nan;reason='FEWER_THAN_10_PAIRED_BLOCKS'
     if n>=10:
      rx=rankdata(x);ry=rankdata(y);rx-=rx.mean();ry-=ry.mean();norm=np.linalg.norm(rx)*np.linalg.norm(ry)
      if norm>0:
       rho=float(rx@ry/norm);check=float(spearmanr(x,y).statistic);err=abs(check-rho);assert err<1e-12;maxrhoerror=max(maxrhoerror,err);checked+=1;reason='DESCRIPTIVE_SHARED_RNA_NOT_INDEPENDENT_EVIDENCE'
      else:reason='CONSTANT_SCORE'
     correlations.append(dict(context=ctx.context,gse_id=ctx.gse_id,CMs=ctx.CMs,type=typ,MCP=mcp,endpoint=endpoint,reference_programme=rid,partition=partition,blocks_before_measurement=len(pp),paired_measured_blocks=n,Spearman_rho=rho,status=reason,median_fraction_cells_expressing_at_least5_signature_genes=float(np.median(ag['informative_gene_support_cells'][rows[valid],col]/ag['valid_score_cells'][rows[valid],col])) if n else np.nan,median_signature_RNA_fraction=float(np.median(ag['programme_RNA_fraction_sum'][rows[valid],col]/ag['valid_score_cells'][rows[valid],col])) if n else np.nan,reference_score_IQR=float(np.quantile(y,.75)-np.quantile(y,.25)) if n else np.nan,MCP_score_IQR=float(np.quantile(x,.75)-np.quantile(x,.25)) if n else np.nan,reference_stable=ref['algorithm_stable_candidate'],reference_native_training_support_fraction=ref['native_training_supported_fraction'],reference_basic_expression=ref['basic_expression_review_flag'],reference_foreign_RNA_sensitive=ref['foreign_RNA_sensitive']))
assert len(slotrows)==56
pd.DataFrame(slotrows).to_csv(O/'all56_MCP_slot_CM_scope.tsv',sep='\t',index=False);pd.DataFrame(overlap).to_csv(O/'all_MCP_reference_gene_overlap.tsv',sep='\t',index=False);pd.DataFrame(correlations).to_csv(O/'all_MCP_UCell_block_correspondence.tsv',sep='\t',index=False);files={p.name:sha(p) for p in O.iterdir() if p.is_file()};out=dict(status='ALL56_FIXED_MCP_SLOTS_BRIDGED_DESCRIPTIVE_NOT_NEW_STATE_OR_LR_VALIDATION',observed_UTC=datetime.now(timezone.utc).isoformat(),plan_sha256=sha(plan_path),projection_inventory_sha256=sha(P/'inventory.json'),dictionary_sha256=sha(D),scope=56,reference_pair_records=len(overlap),score_correspondence_records=len(correlations),independent_rank_checks=checked,max_rank_error=maxrhoerror,files=files,no_programme_selection_by_P=True,formal_states=False,independent_validation=False);(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
