from pathlib import Path
import os,json
import pandas as pd
V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);P=V/'state_after_DIALOGUE_frozen_projection_20261016_v1';B=V/'state_after_DIALOGUE_programme_bridge_20261016_v1';T=V/'state_after_DIALOGUE_internal_holdout_20261016_v1';t=pd.read_csv(T/'public_internal_holdout_results.tsv',sep='\t');scope=pd.read_csv(B/'all56_MCP_slot_CM_scope.tsv',sep='\t');g=pd.read_csv(P/'public_fixed_gene_weights.tsv',sep='\t');o=pd.read_csv(B/'all_MCP_reference_gene_overlap.tsv',sep='\t');c=pd.read_csv(B/'all_MCP_UCell_block_correspondence.tsv',sep='\t');out=[]
for r in t[t.BH_q<.05].itertuples():
 item=dict(context=int(r.context),MCP=r.program,heldout_blocks=int(r.blocks),rho=float(r.holdout_rho),CI=[float(r.ci_low),float(r.ci_high)],BH_q=float(r.BH_q),endpoints=[])
 for typ in ['typeA','typeB']:
  ss=scope[(scope.context==r.context)&scope.type.eq(typ)&scope.MCP.eq(r.program)].iloc[0];gg=g[(g.context==r.context)&g.type.eq(typ)&g.program.eq(r.program)];ov=o[(o.context==r.context)&o.type.eq(typ)&o.MCP.eq(r.program)].sort_values(['jaccard','absolute_weight_fraction_in_reference'],ascending=False).head(3);co=c[(c.context==r.context)&c.type.eq(typ)&c.MCP.eq(r.program)&c.partition.eq('holdout')].copy();co['abs_rho']=co.Spearman_rho.abs();co=co.sort_values('abs_rho',ascending=False).head(3)
  item['endpoints'].append(dict(type=typ,endpoint=ss.endpoint,CMs=ss.CMs,major=ss.major,positive_genes=gg[gg.weight>0].sort_values('weight',ascending=False).head(12)[['gene','weight']].to_dict('records'),negative_genes=gg[gg.weight<0].sort_values('weight').head(12)[['gene','weight']].to_dict('records'),top_reference_overlap=ov[['reference_programme','reference_stable','reference_descriptors','reference_basic_expression','reference_foreign_RNA_sensitive','jaccard','up_overlap','down_overlap','overlap_genes']].to_dict('records'),top_heldout_UCell_correspondence=co[['reference_programme','paired_measured_blocks','Spearman_rho']].to_dict('records')))
 out.append(item)
print(json.dumps(dict(purpose='Descriptive review only; all fixed programmes remain in LR scope regardless of p',records=out),ensure_ascii=False))
