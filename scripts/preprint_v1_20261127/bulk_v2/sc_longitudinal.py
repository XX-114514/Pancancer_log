"""PREREG_bulk_v2 §5: single-cell pre vs post (curated states)."""
import os
import pandas as pd, numpy as np
from scipy.stats import wilcoxon, spearmanr
W = os.environ['TFINE_ROOT']
st = pd.read_csv(W + '/cells_phaseB_states.tsv.gz', sep='\t', index_col=0)
ob = pd.read_csv(W + '/cells_phaseA.tsv.gz', sep='\t', index_col=0, usecols=['row', 'dataset', 'sampleID', 'patient_id', 'tissue', 'timepoint', 'treatment', 'source_response'])
d = st.join(ob)
def fractions(x, key):
    c4 = x[x.lineage == 'CD4'].groupby(key).state.value_counts().unstack(fill_value=0)
    c8 = x[x.lineage == 'CD8'].groupby(key).state.value_counts().unstack(fill_value=0)
    ok = c4.index.intersection(c8.index); ok = [u for u in ok if c4.loc[u].sum() >= 50 and c8.loc[u].sum() >= 50]
    f = pd.DataFrame({'Treg_CD4': c4.loc[ok, 'Treg'] / c4.loc[ok].sum(1), 'Tph_CD4': c4.loc[ok].get('CD4_CXCL13_Tph', 0) / c4.loc[ok].sum(1),
                      'Tex_CD8': c8.loc[ok, 'Tex'] / c8.loc[ok].sum(1)})
    return f
def run(name, pre, post, extra=None):
    a, b = fractions(pre, 'patient_id'), fractions(post, 'patient_id'); pts = a.index.intersection(b.index)
    out = [f'== {name}: patients pre={len(a)} post={len(b)} paired={len(pts)}']
    for k in ['Treg_CD4', 'Tph_CD4', 'Tex_CD8']:
        dl = b.loc[pts, k] - a.loc[pts, k]
        out.append(f'  {k}: median pre {a.loc[pts,k].median():.3f} post {b.loc[pts,k].median():.3f} delta {dl.median():+.3f} Wilcoxon p={wilcoxon(dl).pvalue if len(pts)>=5 else np.nan:.4f}')
    rng = np.random.default_rng(20261126)
    for e1, e2 in [('Treg_CD4', 'Tex_CD8'), ('Tph_CD4', 'Tex_CD8'), ('Treg_CD4', 'Tph_CD4')]:
        r0 = spearmanr(a.loc[pts, e1], a.loc[pts, e2])[0]; r1 = spearmanr(b.loc[pts, e1], b.loc[pts, e2])[0]
        bs = []
        for _ in range(2000):
            s = rng.choice(pts, len(pts), replace=True)
            x0 = spearmanr(a.loc[s, e1], a.loc[s, e2])[0]; x1 = spearmanr(b.loc[s, e1], b.loc[s, e2])[0]
            if np.isfinite(x0) and np.isfinite(x1): bs.append(x1 - x0)
        lo, hi = np.percentile(bs, [2.5, 97.5])
        out.append(f'  edge {e1}~{e2}: pre rho {r0:.2f}  post rho {r1:.2f}  diff {r1-r0:+.2f} (95% boot {lo:+.2f} to {hi:+.2f}, valid {len(bs)})')
    if extra is not None:
        for grp, ids in extra.items():
            ids = [p for p in pts if p in ids]
            if len(ids) >= 4:
                s = '  ' + grp + f' (n={len(ids)}): ' + '; '.join(f'{k} delta {(b.loc[ids,k]-a.loc[ids,k]).median():+.3f}' for k in ['Treg_CD4', 'Tph_CD4', 'Tex_CD8'])
                out.append(s)
    print('\n'.join(out)); return a.join(b, lsuffix='_pre', rsuffix='_post', how='inner')
h = d[d.dataset == 'GSE286827']
T1 = run('GSE286827 HNSCC neoadjuvant ICI (pre vs post)', h[h.timepoint == 'pre'], h[h.timepoint == 'post'])
t = d[(d.dataset == 'GSE169246') & (d.tissue != 'blood')]
print('GSE169246 treatment values:', t.drop_duplicates('patient_id').treatment.value_counts().to_dict(), '| response:', t.drop_duplicates('patient_id').source_response.value_counts().to_dict())
pinfo = t.drop_duplicates('patient_id').set_index('patient_id')
extra = {f'arm={a}': set(pinfo.index[pinfo.treatment == a]) for a in pinfo.treatment.dropna().unique()}
extra.update({f'resp={r}': set(pinfo.index[pinfo.source_response == r]) for r in pinfo.source_response.dropna().unique()})
T2 = run('GSE169246 TNBC tumour (Pre vs Post)', t[t.timepoint == 'Pre'], t[t.timepoint == 'Post'], extra)
for a in pinfo.treatment.dropna().unique():
    ids = pinfo.index[pinfo.treatment == a]; tt = t[t.patient_id.isin(ids)]
    try: run(f'GSE169246 arm {a}', tt[tt.timepoint == 'Pre'], tt[tt.timepoint == 'Post'])
    except Exception as e: print('  arm', a, 'not estimable:', e)
pd.concat({'GSE286827': T1, 'GSE169246': T2}).to_csv(W + '/bulk_v2/sc_pre_post_patient_fractions.tsv', sep='\t')
