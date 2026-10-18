#!/usr/bin/env python3
"""Read-only anonymous evidence for scope-qualified immune CM names; no outcome fitting."""
import argparse, hashlib, json
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import pandas as pd

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main(a):
    catpath=a.run/'tables/subcluster_catalog_refined_private.tsv'
    cat=pd.read_csv(catpath,sep='\t',keep_default_na=False).set_index('subcluster')
    indexpath=a.run/'manifests/CM_recovery_v3_public.json'
    index=json.loads(indexpath.read_text());items=[];files=[]
    for r in index['results']:
        d=a.run/'CM'/r['scope']/r['discovery_attempt']
        W=pd.read_csv(d/'W_private.tsv',sep='\t',index_col=0)
        H=pd.read_csv(d/'H_patient_block_private.tsv',sep='\t',index_col=0)
        members=pd.read_csv(d/'each_node.tsv',sep='\t')
        edges=pd.read_csv(d/'global_edge.tsv',sep='\t')
        meta=pd.read_csv(a.run/'CM'/r['scope']/'inputs/sample_metadata_private.tsv',sep='\t',keep_default_na=False)
        assert list(W.columns)==list(H.index) and W.index.is_unique and H.columns.is_unique
        assert set(W.index)<=set(cat.index) and np.isfinite(W.to_numpy()).all() and (W.to_numpy()>=0).all()
        assert len(W.columns)==r['selected_rank']
        mass=W.sum();assert (mass>0).all();wn=W.divide(mass,axis=1)
        h=H.mul(mass,axis=0);h=h.divide(h.sum(axis=0),axis=1)
        assert np.allclose(wn.sum(),1) and np.allclose(h.sum(),1)
        mapping=meta[['independent_unit','dataset']].drop_duplicates()
        assert not mapping.independent_unit.duplicated().any()
        source=mapping.set_index('independent_unit').dataset.reindex(H.columns);assert source.notna().all()
        sm=h.T.groupby(source.to_numpy()).median()
        source_n=source.value_counts().to_dict()
        for c in W.columns:
            nodes=members.loc[members.cm.eq(c),'subCluster'].tolist()
            assert len(nodes)==len(set(nodes)) and set(nodes)<=set(W.index)
            internal=edges.loc[edges.subCluster1.isin(nodes)&edges.subCluster2.isin(nodes)] if len(edges) else edges
            material=set(wn.index[wn[c]>=.01])
            material_internal=internal.loc[internal.subCluster1.isin(material)&internal.subCluster2.isin(material)] if len(internal) else internal
            majors=wn[c].groupby(cat.major.reindex(W.index)).sum().sort_values(ascending=False)
            top=[]
            for n,w in wn[c].nlargest(8).items():
                row=cat.loc[n]
                top.append(dict(node=n,weight=float(w),major=str(row.major),marker=str(row.marker),
                  subtype_candidate=str(row.subtype_candidate),identity_promoted=str(row.identity_promoted),
                  network_member=n in nodes,n_cells=int(row.n_cells)))
            core=[]
            for n in sorted(nodes,key=lambda x:-float(wn.loc[x,c])):
                row=cat.loc[n];core.append(dict(node=n,weight=float(wn.loc[n,c]),major=str(row.major),marker=str(row.marker),subtype_candidate=str(row.subtype_candidate)))
            items.append(dict(scope=r['scope'],factor=c,version='2026-10-14 immune atlas v1',
              discovery_attempt=r['discovery_attempt'],n_sources=r['n_sources'],n_blocks=r['n_composition_blocks'],
              n_verified_patient_blocks=r['n_verified_patient_blocks'],n_unverified_library_units=r['n_unverified_library_units'],
              rank_at_survey_boundary=r['rank_at_survey_boundary'],LOSO_refit_count=r['LOSO_refit_count'],
              LOSO_unestimable_count=r['LOSO_unestimable_count'],
              n_network_members=len(nodes),n_internal_edges=len(internal),
              diagnostic_weight_cutoff=.01,n_nodes_weight_ge_1pct=len(material),
              maximum_node_weight=float(wn[c].max()),effective_loading_node_number=float(1/(wn[c]**2).sum()),
              n_internal_edges_both_endpoints_weight_ge_1pct=len(material_internal),
              supported_member_weight=float(wn.loc[nodes,c].sum()),
              major_weight={str(k):float(v) for k,v in majors.items()},top_loading_nodes=top,network_members=core,
              source_median_relative_contribution={str(k):float(v) for k,v in sm[c].items()},
              source_block_counts={str(k):int(v) for k,v in source_n.items()}))
        for name in ['W_private.tsv','H_patient_block_private.tsv','each_node.tsv','global_edge.tsv']:
            p=d/name;files.append(dict(logical_path=f"CM/{r['scope']}/{r['discovery_attempt']}/{name}",sha256=sha(p)))
    assert len(items)==79 and sum(x['n_network_members']>0 for x in items)==62
    result=dict(observed_utc=datetime.now(timezone.utc).isoformat(),input_version='20261013_immunotherapy_atlas_v1',
      namespace='2026-10-14 immune atlas CM v1; scope required; distinct from v20261010 nine CMs',
      no_outcome_tests=True,no_patient_or_cell_identifiers=True,original_models_read_only=True,
      n_factor_instances=len(items),n_instances_with_network_members=sum(x['n_network_members']>0 for x in items),
      n_instances_with_internal_edges=sum(x['n_internal_edges']>0 for x in items),
      overlapping_scopes_not_independent_modules=True,catalog_sha256=sha(catpath),CM_index_sha256=sha(indexpath),
      producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest() if '__file__' in globals() and Path(__file__).is_file() else None,
      input_files=files,modules=items)
    text=json.dumps(result,indent=2,ensure_ascii=False)
    if a.output=='-':print(text)
    else:
        with Path(a.output).open('x') as f:f.write(text+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--output',default='-')
    main(p.parse_args())
