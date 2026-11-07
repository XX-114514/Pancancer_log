"""Canonical IL17E/IL25 alias review only; original DE, enrichment and network kept."""
import json
from datetime import datetime, timezone
import numpy as np,pandas as pd,anndata as ad

r=PROJECT_ROOT/'runs/20261107_V4_CoVarNet_communication_v1';out=r/'audit/producer_alias_review_v2';out.mkdir(exist_ok=False)
source='https://www.ncbi.nlm.nih.gov/gene/64806'
a=ad.read_h5ad(r/'inputs/full_common_PRIVATE.h5ad');response=pd.read_csv(r/'results/cytokine/response_full_family.tsv',sep='\t');mem=pd.read_csv(r/'audit/network_members.tsv',sep='\t')
mapping=pd.read_csv(r/'results/cytokine/cytokine_producer_symbol_mapping.tsv',sep='\t').fillna('')
row=mapping.cytokine_symbol_mouse.eq('IL17E');assert row.sum()==1
mapping.loc[row,'human_subunits']='IL25';mapping.loc[row,'status']='canonical_IL17E_alias_to_IL25_verified_NCBI_HGNC'
mapping['canonical_subunits_measured_in_common_input']=mapping.human_subunits.map(lambda s:bool(s) and all(x in a.var_names for x in s.split(';')))
mapping.to_csv(out/'producer_symbol_mapping_reviewed_v2.tsv',sep='\t',index=False)
network=pd.read_csv(r/'results/cytokine/production_response_network.tsv',sep='\t');added=[]
for branch,qkey in [('author_unsigned_annotation_universe','q_within_CM_member'),('directional_measured_universe','q_full_family')]:
    valid=response[response.branch.eq(branch)&response[qkey].lt(.05)]
    for e in valid[valid.cytokine_symbol_mouse.eq('IL17E')].itertuples():
        if 'IL25' not in a.var_names:continue
        responders=set(valid.loc[valid.cm.eq(e.cm)&valid.cytokine.eq(e.cytokine),'node'])
        for prod in mem.loc[mem.cm.eq(e.cm),'subCluster']:
            if prod in responders:continue
            mask=a.obs.subCluster.eq(prod)&a.obs.CMT.eq(e.cm)
            if not mask.any():continue
            mean=float(a[mask,'IL25'].X.mean())
            if mean>.1:added.append({'branch':branch,'cm':e.cm,'producer':prod,'cytokine':e.cytokine,'subunit':'IL25','receiver':e.node,'direction':e.direction,'mean_expression':mean,'response_q':getattr(e,qkey),'all_required_subunits_mean_gt_0_1':True,'meaning':'RNA production-response association; canonical IL17E/IL25 alias, not protein or intervention'})
reviewed=pd.concat([network[network.subunit.ne('IL17E')],pd.DataFrame(added)],ignore_index=True)
reviewed.to_csv(out/'production_response_network_reviewed_v2.tsv',sep='\t',index=False)
assert reviewed.mean_expression.gt(.1).all() and reviewed.response_q.lt(.05).all()
result={'status':'PASS_canonical_IL17E_IL25_alias_review_no_DE_or_enrichment_refit',
    'UTC':datetime.now(timezone.utc).isoformat(),'primary_NCBI_HGNC_source':source,
    'canonical_IL25_measured_in_all13_source_common_input':'IL25' in a.var_names,
    'n_original_subunit_edges':len(network),'n_alias_added_edges':len(added),'n_reviewed_subunit_edges':len(reviewed),
    'original_outputs_preserved':True,'n_DE_or_enrichment_models_rerun':0,
    'direction_semantics':'up/down refer to the gene-signature direction in cytokine perturbation; down-signature concordance is not proof of decreased cytokine activity',
    'causal_or_protein_claim':False}
(out/'complete.json').open('x').write(json.dumps(result,indent=2));print(json.dumps(result))
