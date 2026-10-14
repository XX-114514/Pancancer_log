"""Frozen exploratory MMR cohort test; all scientific outputs stay on server.

Separate discovery / holdout / combined; original CM discovery used this cohort.
MMR is molecular subtype, not treatment response. No automatic covariate removal.
"""
import hashlib
import itertools
import json
import math
import os
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata._io.specs import read_elem
from scipy.stats import rankdata, t


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def rank_effect(y, exposed):
    n1, n0 = int(exposed.sum()), int((~exposed).sum())
    ranks = rankdata(y)
    u = ranks[exposed].sum() - n1 * (n1 + 1) / 2
    return float(2 * u / (n0 * n1) - 1)


def hc3(y, X):
    if len(y) < 20 or len(y) - X.shape[1] < 10:
        return None, 'insufficient_patients_or_residual_df'
    if np.linalg.matrix_rank(X) != X.shape[1] or np.linalg.cond(X) > 1e8:
        return None, 'unidentifiable_design_no_covariate_dropped'
    bread = np.linalg.inv(X.T @ X)
    coef = bread @ X.T @ y
    residual = y - X @ coef
    leverage = np.sum((X @ bread) * X, axis=1)
    if np.max(leverage) >= 1 - 1e-8:
        return None, 'saturated_leverage'
    cov = bread @ (X.T @ (((residual / (1 - leverage)) ** 2)[:, None] * X)) @ bread
    se = np.sqrt(cov[1, 1])
    if not np.isfinite(se) or se <= 0:
        return None, 'undefined_HC3'
    return dict(beta=float(coef[1]), HC3_p=float(2*t.sf(abs(coef[1]/se), len(y)-X.shape[1]))), 'estimated'


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    origin = root / 'CRC_fine_pair_DIALOGUE_20261014_v1'
    dest = origin / 'private_MMR_cohort_test_v1'
    assert not dest.exists()
    dest.mkdir()
    # Write the contract before reading outcomes or performing tests.
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(), producer_sha256=sha(__file__),
                seed=20261014, programs=['MCP1','MCP2'], endpoints=['Tnode9','Macnode18'],
                partitions=['discovery','holdout','combined_exploratory'],
                primary='two_sided_patient_label_permutation_rank_biserial_MMRd_minus_MMRp',
                rank_minimum_each_group=5, permutations=20000, bootstrap=2000,
                adjustment='MMRd+right_site+numeric_pathological_T_stage+both_original_fine_major_fractions',
                adjustment_minimum_patients=20, residual_df_minimum=10,
                no_automatic_covariate_removal=True, sex_missing_not_inferred=True,
                BH_family='all four program_endpoint slots per partition per method; unsupported p=1 internally',
                interpretation='cohort molecular subtype association; neither independent project validation nor therapy prediction',
                clinical_results_private=True)
    (dest / 'plan_before_results.json').write_text(json.dumps(plan, indent=2)+'\n')
    metadata_path = origin/'private_MMR_metadata_audit_v2/all38_exact_private_author_metadata.tsv'
    scores_path = origin/'private_frozen_projection_v1/all_patient_mean_scores_private.tsv'
    md = pd.read_csv(metadata_path, sep='\t', keep_default_na=False)
    scores = pd.read_csv(scores_path, sep='\t')
    src = root/'state_singlecell_inputs_20261014_v1/private_native_sources/GSE178341.h5ad'
    with h5py.File(src) as f:
        obs = read_elem(f['obs'])
    node = obs.groupby(['sample_uid','subCluster'], observed=True).size()
    major = obs.groupby(['sample_uid','major'], observed=True).size()
    labels = [('CD8_T_cell_NKG7_9','T cell'), ('Macrophage_IFI6_18','Macrophage')]
    for i, (label, parent) in enumerate(labels):
        md['fraction'+str(i)] = [node.loc[s,label]/major.loc[s,parent] for s in md.sample_uid]
    md['MMRd'] = md.mmrstatus_author.map({'MMRd':1.,'MMRp':0.})
    md['right_site'] = md.tissuesitesimple.map({'right':1.,'left':0.})
    md['T_stage'] = pd.to_numeric(md.tumorstage.str.lower().str.extract(r'^pt([1-4])(?:[ab])?$')[0], errors='coerce')
    rows = []
    for partition in plan['partitions']:
        scope = md if partition == 'combined_exploratory' else md[md.partition.eq(partition)]
        for program, endpoint in itertools.product(plan['programs'], plan['endpoints']):
            gs = scores[scores.program.eq(program)&scores.cell_type.eq(endpoint)]
            data = scope.merge(gs[['sample','score']], left_on='analysis_block', right_on='sample', how='left', validate='one_to_one')
            assert not data.analysis_block.duplicated().any()
            for method in ['rank_permutation','adjusted_HC3_exploratory']:
                row = dict(partition=partition, program=program, endpoint=endpoint, method=method,
                           status='no_gene_supported_score', n_scope=len(data), n_complete=0, n_MMRd=0, n_MMRp=0,
                           effect=np.nan, lower=np.nan, upper=np.nan, p=np.nan, invalid_bootstrap=0)
                if data.score.notna().any():
                    cols = ['score','MMRd'] + (['right_site','T_stage','fraction0','fraction1'] if method.startswith('adjusted') else [])
                    use = data.dropna(subset=cols)
                    n = len(use); g = use.MMRd.to_numpy().astype(bool); y = use.score.to_numpy()
                    row.update(n_complete=n,n_MMRd=int(g.sum()),n_MMRp=int((~g).sum()),n_missing=len(data)-n)
                    row['status'] = 'insufficient_MMR_group_patients'
                    if min(row['n_MMRd'],row['n_MMRp']) >= 5:
                        seed = 20261014 + int(hashlib.sha256((partition+program+endpoint+method).encode()).hexdigest()[:8],16)
                        rng = np.random.default_rng(seed)
                        boots = []
                        if method == 'rank_permutation':
                            effect = rank_effect(y,g); ranks = rankdata(y); n1 = int(g.sum()); n0=n-n1
                            combinations = math.comb(n,n1)
                            extreme=0
                            if combinations <= 20000:
                                for selected in itertools.combinations(range(n),n1):
                                    e=2*(ranks[list(selected)].sum()-n1*(n1+1)/2)/(n0*n1)-1
                                    extreme += abs(e) >= abs(effect)-1e-12
                                p=extreme/combinations; row['permutation_mode']='exact_all_assignments'
                            else:
                                for _ in range(20000):
                                    e=2*(ranks[rng.choice(n,n1,replace=False)].sum()-n1*(n1+1)/2)/(n0*n1)-1
                                    extreme += abs(e) >= abs(effect)-1e-12
                                p=(extreme+1)/20001; row['permutation_mode']='Monte_Carlo_20000'
                            for _ in range(2000):
                                ix=np.r_[rng.choice(np.flatnonzero(g),n1,replace=True),rng.choice(np.flatnonzero(~g),n0,replace=True)]
                                boots.append(rank_effect(y[ix],g[ix]))
                            row.update(status='estimated',effect=effect,p=p)
                        else:
                            covs=use[['right_site','T_stage','fraction0','fraction1']].to_numpy()
                            covs[:,1:]=(covs[:,1:]-covs[:,1:].mean(0))/np.where(covs[:,1:].std(0)>0,covs[:,1:].std(0),1)
                            X=np.column_stack((np.ones(n),g.astype(float),covs))
                            fit, status=hc3(y,X); row['status']=status
                            if fit is not None:
                                row.update(effect=fit['beta'],p=fit['HC3_p'])
                                for _ in range(2000):
                                    ix=rng.integers(0,n,n); b,_=hc3(y[ix],X[ix])
                                    if b is None:row['invalid_bootstrap']+=1
                                    else:boots.append(b['beta'])
                        if len(boots)>=1800:
                            row['lower'],row['upper']=map(float,np.quantile(boots,[.025,.975]))
                        elif boots:row['interval_status']='too_many_unidentifiable_bootstraps'
                rows.append(row)
    result=pd.DataFrame(rows)
    result['BH_four_slot_family']=np.nan
    for _,ix in result.groupby(['partition','method']).groups.items():
        assert len(ix)==4
        p=result.loc[ix,'p'].fillna(1).to_numpy(); order=np.argsort(p)
        q=np.empty(4);q[order]=np.minimum(1,np.minimum.accumulate((p[order]*4/np.arange(1,5))[::-1])[::-1])
        result.loc[ix,'BH_four_slot_family']=np.where(result.loc[ix,'p'].notna(),q,np.nan)
    result.to_csv(dest/'all24_private_results.tsv',sep='\t',index=False,na_rep='NA')
    inventory=dict(status='MMR_COHORT_EXPLORATION_COMPLETE_RESULTS_PRIVATE_NOT_VERIFIED_INTERACTION',
                   finished_UTC=datetime.now(timezone.utc).isoformat(), producer_sha256=sha(__file__),
                   inputs_sha256={'metadata':sha(metadata_path),'frozen_scores':sha(scores_path)},
                   no_new_MCP_or_LR_fitting=True, no_patient_export=True, publication_pending=True,
                   limitations=['Sex/age not available in cached source','MMR is not ICI response',
                                'Combined/discovery are exploratory','Original CM discovery already used this source',
                                'Rank permutation assumes exchangeable patients; adjusted HC3 is separate exploratory model'])
    (dest/'private_inventory.json').write_text(json.dumps(inventory,indent=2)+'\n')
    print('MMR_COHORT_EXPLORATION_COMPLETE_RESULTS_REMAIN_SERVER')


if __name__=='__main__':
    main()
