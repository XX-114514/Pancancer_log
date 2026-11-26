import pandas as pd, numpy as np, json, os
R = os.environ['TFINE_ROOT'] + '/covarnet_rerun'
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 40)
Ws = {}
for s in ['A_original_units', 'B_baseline', 'R_replication']:
    o = f'{R}/{s}/discovery_v1'
    if not os.path.exists(o + '/CM_public.json'):
        print(s, 'not finished'); continue
    rep = json.load(open(o + '/CM_public.json'))
    print(f'\n===== {s}: K={rep["selected_rank"]} blocks={rep["n_composition_blocks"]} sources={rep["n_sources"]} edges={rep["n_network_edges"]}')
    print(pd.read_csv(o + '/rank_survey.tsv', sep='\t').round(3).T.to_string())
    W = pd.read_csv(o + '/W_private.tsv', sep='\t', index_col=0); Ws[s] = W
    share = W / W.sum(0)
    print('W share per CM (columns sum to 1):'); print(share.round(2).to_string())
    for f in ['CM_members.tsv', 'CM_members_topn3.tsv', 'CM_members_K4_topn3.tsv']:
        p = f'{o}/{f}'
        if os.path.exists(p):
            m = pd.read_csv(p, sep='\t')
            print(f, {cm: [f'{a}({w:.2g})' for a, w in zip(g.subCluster, g.weight)] for cm, g in m.groupby('cm')})
    if os.path.exists(o + '/W_K4_private.tsv'):
        W4 = pd.read_csv(o + '/W_K4_private.tsv', sep='\t', index_col=0)
        print('K4 W share:'); print((W4 / W4.sum(0)).round(2).to_string())
    e = pd.read_csv(o + '/source_adjusted_pair_correlations_private.tsv', sep='\t')
    key = e[((e.subCluster1 == 'Treg') | (e.subCluster2 == 'Treg')) | ((e.subCluster1 == 'Tex') | (e.subCluster2 == 'Tex'))]
    print('source-adjusted pair correlations involving Treg/Tex:')
    print(key[['subCluster1', 'subCluster2', 'correlation', 'pval_fdr', 'spe', 'n_positive_sources']].sort_values('correlation', ascending=False).round(3).head(12).to_string(index=False))
    if os.path.exists(o + '/LOSO_stability.tsv'):
        print('LOSO best cosine:'); print(pd.read_csv(o + '/LOSO_stability.tsv', sep='\t').pivot(index='cm', columns='heldout_source', values='best_cosine').round(3).to_string())

def cos(A, B):
    A = A / np.linalg.norm(A, axis=0); B = B / np.linalg.norm(B, axis=0); return A.T @ B
for a, b in [('A_original_units', 'R_replication'), ('B_baseline', 'R_replication'), ('A_original_units', 'B_baseline')]:
    if a in Ws and b in Ws:
        A, B = Ws[a], Ws[b].loc[Ws[a].index]
        print(f'\ncosine {a} (rows) vs {b} (cols)'); print(pd.DataFrame(cos(A.values, B.values), index=A.columns, columns=B.columns).round(2).to_string())
