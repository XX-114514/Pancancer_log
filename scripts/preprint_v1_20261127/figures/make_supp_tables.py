"""Supplementary Tables S1 (modules) and S2 (source registry) for preprint v1.
Inputs: public Pancancer_log inventories (S1a, S2) and private rerun outputs (S1b, S1c; supplement only, not for Git)."""
import pandas as pd, numpy as np, sys, os
LOG, W0, OUT = sys.argv[1], sys.argv[2], sys.argv[3]; os.makedirs(OUT, exist_ok=True)
inv = LOG + '/inventories'
# ---- S1a: R4/K15 module members (public aggregate) ----
m = pd.read_csv(inv + '/analysis/2026-11-19/WEHI_R4_CM_validation/results/CM_validation_20261118/R4_member_comparison_eligibility.tsv', sep='\t')
m.insert(2, 'lineage', m.node.str.split('_').str[0])
m.rename(columns={'cm': 'R4_K15_module', 'node': 'identity'}).to_csv(OUT + '/Table_S1a_R4_K15_module_members.tsv', sep='\t', index=False)
summ = m.groupby('cm').agg(n_members=('node', 'size'), sources_with_high_units=('high_sources', 'first'),
                           largest_single_source_share=('largest_high_source_fraction', 'first'),
                           members=('node', lambda s: '; '.join(s)))
summ['Tph_Treg_Tex_all_members'] = m.groupby('cm').node.apply(lambda s: {'CD4T_Tfh-like_CXCL13', 'CD4T_Treg_FOXP3', 'CD8T_Tex_LAG3'} <= set(s))
summ.index.name = 'R4_K15_module'; summ.to_csv(OUT + '/Table_S1a_R4_K15_module_summary.tsv', sep='\t')
# ---- S1b/c: curated CoVarNet rerun (supplement only) ----
rows = []
for d in ['A_original_units', 'B_baseline', 'R_replication']:
    p = f'{W0}/covarnet_rerun/{d}/discovery_v1'
    W = pd.read_csv(p + '/W_private.tsv', sep='\t', index_col=0); W = W / W.sum(0)
    for cm in W: rows += [dict(design=d, module=cm, state=s, W_share=round(v, 4)) for s, v in W[cm].items()]
pd.DataFrame(rows).to_csv(OUT + '/Table_S1b_curated_CoVarNet_W_shares_SUPPLEMENT_ONLY.tsv', sep='\t', index=False)
e = []
for d in ['A_original_units', 'B_baseline', 'R_replication']:
    c = pd.read_csv(f'{W0}/covarnet_rerun/{d}/discovery_v1/source_adjusted_pair_correlations_private.tsv', sep='\t')
    c = c[['subCluster1', 'subCluster2', 'correlation', 'pval_fdr', 'spe', 'n_positive_sources']]; c.insert(0, 'design', d); e.append(c.round(4))
pd.concat(e).to_csv(OUT + '/Table_S1c_curated_source_adjusted_edges_SUPPLEMENT_ONLY.tsv', sep='\t', index=False)
# ---- S2: source registry (public aggregate) ----
dd = pd.read_csv(inv + '/tme_v1/dataset_detail_by_source.tsv', sep='\t')
sr = pd.read_csv(inv + '/source_registry/source_registry_v1.tsv', sep='\t')
keep_d = ['source', 'analysis_main', 'technology', 'units', 'patients', 'cells_integration_v2', 'malignant', 'uncertain', 'non_malignant',
          'tme_v1_status', 'cells_tme_v1', 'capture_scope', 'treatment_status', 'sampling_site', 'fresh_frozen', 'suspension', 'patient_id_rule']
keep_s = ['source', 'pmid', 'journal', 'pub_year', 'first_author', 'paper_title', 'chemistry_from_protocols', 'geo_samples_total', 'units_used']
S2 = dd[keep_d].merge(sr[keep_s], on='source', how='left')
ex = pd.read_csv(inv + '/tme_v1/scope_exclusions_v1.tsv', sep='\t')
S2['scope_exclusion'] = S2.source.map(ex[ex.level == 'source'].set_index('source').reason)
dec = pd.read_csv(inv + '/tme_v1/source_decisions_v1.tsv', sep='\t')
S2['capture_decision'] = S2.source.map(dec.groupby('source').decision.agg('; '.join))
S2.to_csv(OUT + '/Table_S2_source_registry.tsv', sep='\t', index=False)
print('S1a modules', len(summ), '| S1b rows', len(rows), '| S2 sources', len(S2), '| cells v2', int(S2.cells_integration_v2.sum()), '| tme', int(S2.cells_tme_v1.fillna(0).sum()))
print(S2.tme_v1_status.value_counts().to_dict())
