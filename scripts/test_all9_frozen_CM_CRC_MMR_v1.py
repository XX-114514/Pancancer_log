"""All nine frozen composition CMs: source-cohort MMR, capture gated.
Original H is a descriptive activity, not independent module recovery.
All real results/patient keys remain on server.
"""
from datetime import datetime,timezone
import hashlib,json,os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from cm_immune_recovery_diagnostic_v1 import COMPARTMENT
from test_CRC_DIALOGUE_MMR_cohort_v1 import rank_effect,hc3


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    downstream=Path(os.environ['PANCANCER_FROZEN_DOWNSTREAM'])
    out=root/'CRC_all9_frozen_CM_MMR_20261014_v1'
    assert not out.exists();out.mkdir()
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),
      scope='original_unsorted_primary_CRC_MAJ158',modules=['CM%02d'%i for i in range(1,10)],
      core_definition='actual frozen graph edge endpoint union; original W relative core weights',
      capture_gate='at least half core nodes>=20cells, at least min(3,coremajors),min(2,corecompartments),coreweight>=0.6',
      patient_selection='first lexicographic capture-eligible library per reviewed source patient; no MMR/activity selection',
      outcome='original H continuous score; no module refit or formal clinical validation',
      rank_min_each_MMR=5,rank_min_patients=20,permutations=20000,bootstrap=2000,
      adjusted_model='MMRd+right_site+pathological_T_numeric+log1p_all_cell_count',
      adjusted_min_patients=20,residual_df_min=10,full_rank_required=True,
      families='9slots each primary_rank and exploratory_HC3, failed slots p1 internally; missing published as NA',
      limits=['posthoc same-source ecological association','MMR not therapy response','no source sex/age covariates',
              'pathological T stage is not complete clinical stage','capture selection can differ across MMR groups'],
      private_outputs=True)
    (out/'plan_before_results.json').write_text(json.dumps(plan,indent=2)+'\n')
    parent=root/'CM_interface_cancer_association_v1/source_unit_review_v1'
    side=pd.read_csv(parent/'public_library_patient_tissue_processing_sidecar.tsv',sep='\t',keep_default_na=False)
    from audit_CRC_DIALOGUE_MMR_metadata_v2 import parse_soft
    author=parse_soft(parent/'GSE178341_family.soft.gz')
    native=np.load(root/'immune_same_name_family_source_review_v1/major_pair_native_assembly_v1/private_group_arrays/MAJ158.npz')
    scope=pd.DataFrame({'sampleID':native['samples'],'analysis_block':native['corrected_blocks_for_next_inference'],
                        'reviewed':native['corrected_identity_reviewed_for_next_inference']})
    scope=scope[scope.reviewed]
    scope=scope.merge(side[side.mapping_status.eq('exact_public_library_match')],on='sampleID',validate='one_to_one')
    assert len(scope)==int(native['corrected_identity_reviewed_for_next_inference'].sum())
    assert scope.analysis_block.eq(scope.reviewed_analysis_block).all()
    scope=scope.merge(author[['gsm','mmrstatus','tissuesitesimple','tumorstage']],on='gsm',validate='many_to_one',suffixes=('_accepted','_author'))
    assert scope.mmrstatus_accepted.eq(scope.mmrstatus_author).all()
    scope['MMRd']=scope.mmrstatus_author.map({'MMRd':1.,'MMRp':0.})
    scope['right_site']=scope.tissuesitesimple.map({'right':1.,'left':0.})
    scope['T_stage']=pd.to_numeric(scope.tumorstage.str.lower().str.extract(r'^pt([1-4])(?:[ab])?$')[0],errors='coerce')
    paths=dict(counts=downstream/'inputs/sample_major_subcluster_counts.tsv',
               catalog=downstream/'downstream_review_v2/number_to_annotation_dictionary.tsv',
               W=downstream/'covarnet/tables/nmf_W_subcluster_by_cm.tsv',
               H=downstream/'covarnet/tables/nmf_H_cm_by_sample.tsv',
               edges=downstream/'downstream_review_v2/corrected_CM_edges_source_and_marker_evidence.tsv')
    counts=pd.read_csv(paths['counts'],sep='\t');counts=counts[counts.sampleID.isin(scope.sampleID)]
    catalog=pd.read_csv(paths['catalog'],sep='\t').set_index('subCluster')
    W=pd.read_csv(paths['W'],sep='\t',index_col=0);H=pd.read_csv(paths['H'],sep='\t',index_col=0)
    assert W.shape==(172,9) and H.shape==(9,1290) and set(W.columns)==set(plan['modules'])
    assert set(H.index)==set(W.columns) and set(W.index)==set(catalog.index)
    edges=pd.read_csv(paths['edges'],sep='\t');assert len(edges)==120
    node=counts.pivot_table(index='sampleID',columns='subCluster',values='n_cells',aggfunc='sum',fill_value=0).reindex(index=scope.sampleID,columns=W.index,fill_value=0)
    assert node.index.is_unique and not node.isna().any().any()
    total=counts.groupby('sampleID').n_cells.sum();scope['log_total_cells']=np.log1p(scope.sampleID.map(total))
    rows=[];coverage=[]
    for cm in plan['modules']:
        edge=edges[edges.cm.eq(cm)];core=sorted(set(edge.subCluster1)|set(edge.subCluster2));assert core
        weights=W.loc[core,cm];weights=weights/weights.sum()
        majors=catalog.loc[core,'majorCluster'];compartments=majors.map(COMPARTMENT)
        if compartments.isna().any():
            coverage.append(dict(cm=cm,core_nodes=len(core),status='unmapped_compartment_not_imputed'))
            for method in ['primary_rank','exploratory_adjusted_HC3']:
                rows.append(dict(cm=cm,method=method,status='unmapped_compartment_not_imputed',n_patients=0,
                                 n_MMRd=0,n_MMRp=0,effect=np.nan,p=np.nan,lower=np.nan,upper=np.nan,invalid_bootstrap=0))
            continue
        present=node[core]>=20;weighted=present.to_numpy()@weights.to_numpy()
        nnode=present.sum(1).to_numpy()
        nmajor=np.array([majors[present.iloc[i]].nunique() for i in range(len(node))])
        ncomp=np.array([compartments[present.iloc[i]].nunique() for i in range(len(node))])
        capture=(nnode>=int(np.ceil(len(core)/2)))&(weighted>=.6)&(nmajor>=min(3,majors.nunique()))&(ncomp>=min(2,compartments.nunique()))
        candidates=scope.loc[capture].sort_values(['analysis_block','sampleID']).drop_duplicates('analysis_block').copy()
        candidates['activity']=H.loc[cm,candidates.sampleID].to_numpy()
        assert candidates.activity.notna().all()
        candidates.to_csv(out/(cm+'_private_capture_selected_patients.tsv'),sep='\t',index=False)
        coverage.append(dict(cm=cm,core_nodes=len(core),core_majors=majors.nunique(),core_compartments=compartments.nunique(),
                             available_libraries=len(scope),captured_libraries=int(capture.sum()),selected_patients=len(candidates)))
        for method in ['primary_rank','exploratory_adjusted_HC3']:
            cols=['activity','MMRd']+(['right_site','T_stage','log_total_cells'] if method.endswith('HC3') else [])
            use=candidates.dropna(subset=cols);g=use.MMRd.to_numpy().astype(bool);y=use.activity.to_numpy();n=len(y)
            row=dict(cm=cm,method=method,status='insufficient_capture_or_MMR_groups',n_patients=n,n_MMRd=int(g.sum()),
                     n_MMRp=int((~g).sum()),effect=np.nan,p=np.nan,lower=np.nan,upper=np.nan,invalid_bootstrap=0)
            if n>=20 and min(int(g.sum()),int((~g).sum()))>=5:
                rng=np.random.default_rng(20261014+int(hashlib.sha256((cm+method).encode()).hexdigest()[:8],16));boots=[]
                if method=='primary_rank':
                    effect=rank_effect(y,g);ranks=rankdata(y);n1=int(g.sum());n0=n-n1;extreme=0
                    for _ in range(20000):
                        u=ranks[rng.choice(n,n1,replace=False)].sum()-n1*(n1+1)/2
                        extreme+=abs(2*u/(n1*n0)-1)>=abs(effect)-1e-12
                    row.update(status='estimated_source_cohort_exploration',effect=effect,p=(extreme+1)/20001)
                    for _ in range(2000):
                        ix=np.r_[rng.choice(np.flatnonzero(g),n1,replace=True),rng.choice(np.flatnonzero(~g),n0,replace=True)]
                        boots.append(rank_effect(y[ix],g[ix]))
                else:
                    X=np.column_stack((np.ones(n),g.astype(float),use[['right_site','T_stage','log_total_cells']].to_numpy()))
                    fit,status=hc3(y,X);row['status']=status
                    if fit is not None:
                        row.update(effect=fit['beta'],p=fit['HC3_p'])
                        for _ in range(2000):
                            ix=rng.integers(0,n,n);b,_=hc3(y[ix],X[ix])
                            if b is None:row['invalid_bootstrap']+=1
                            else:boots.append(b['beta'])
                if len(boots)>=1800:row['lower'],row['upper']=map(float,np.quantile(boots,[.025,.975]))
            rows.append(row)
    result=pd.DataFrame(rows);assert len(result)==18
    for _,ix in result.groupby('method').groups.items():
        assert len(ix)==9;p=result.loc[ix,'p'].fillna(1).to_numpy();order=np.argsort(p);q=np.empty(9)
        q[order]=np.minimum(1,np.minimum.accumulate((p[order]*9/np.arange(1,10))[::-1])[::-1])
        result.loc[ix,'BH_nine_slot_family']=np.where(result.loc[ix,'p'].notna(),q,np.nan)
    result.to_csv(out/'all18_private_module_MMR_results.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(coverage).to_csv(out/'all9_private_capture_coverage.tsv',sep='\t',index=False)
    (out/'private_terminal.json').write_text(json.dumps(dict(status='ALL9_FROZEN_CM_CRC_MMR_COMPLETE_PRIVATE_NOT_INDEPENDENT_VALIDATION',
      finished_UTC=datetime.now(timezone.utc).isoformat(),inputs_sha256={k:sha(p) for k,p in paths.items()},
      no_new_module_fit=True,not_verified_interaction=True),indent=2)+'\n')
    print('ALL9_FROZEN_CM_CRC_MMR_COMPLETE_RESULTS_REMAIN_SERVER')


if __name__=='__main__':
    main()
