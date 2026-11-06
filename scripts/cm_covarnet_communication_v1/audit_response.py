from pathlib import Path
import pandas as pd,numpy as np,json,anndata as ad
from scipy.stats import ttest_ind_from_stats,hypergeom
r=Path(__file__).resolve().parents[1];o=r/'results/cytokine';a=pd.read_csv(o/'response_full_family.tsv',sep='\t');tested=a[a.status=='tested'];expect=hypergeom.sf(tested.overlap-1,tested.universe,tested.signature_size,tested.DE_size);err=float(np.max(np.abs(expect-tested.p)));assert err<1e-12
assert ((tested.q_full_family>=0)&(tested.q_full_family<=1)).all();edges=pd.read_csv(o/'production_response_network.tsv',sep='\t');assert (edges.mean_expression>.1).all() and (edges.response_q<.05).all();assert not (edges.producer==edges.receiver).any()
# Preserve full tables; publish an explicit NA-masked review view for untested slots.
a.loc[a.status!='tested',['q_within_identity','q_full_family']]=np.nan;a.to_csv(o/'response_full_family_reviewed.tsv',sep='\t',index=False)
z=np.load(r/'private/CMT_identity_stats.npz');ix=pd.read_csv(r/'private/CMT_identity_stats_index.tsv',sep='\t');b=ad.read_h5ad(r/'inputs/full_LR.h5ad',backed='r');checks=[]
for _,row in ix[(ix.CMT=='CM01')&(ix.identity.isin(['Bcell_Naive_TCL1A','Macrophage_SPP1_TREM2_SPP1']))].iterrows():
 mask=(b.obs.CMT==row.CMT)&(b.obs.identity==row.identity);x=b[mask].X;gi=pd.Index(z['genes']).get_indexer(b.var_names);num=z['n'][row.name];assert mask.sum()==num;mean=np.asarray(x.mean(axis=0)).ravel();saved=z['sums'][row.name,gi]/num;delta=float(np.max(np.abs(mean-saved)));assert delta<.001;checks.append({'identity':row.identity,'n_cells':int(num),'max_mean_error':delta})
res={'status':'response_inference_readback_pass','hypergeometric_max_error':err,'n_tests':len(tested),'n_unqualified':int((a.status!='tested').sum()),'response_full_family_q_lt_05':a[a.q_full_family<.05].groupby('branch').size().to_dict(),'n_network_edges':len(edges),'actual_expression_probes':checks,'patient_inference':False,'causal_inference':False};(r/'audit/response_acceptance.json').write_text(json.dumps(res,indent=2));print(res)
