"""Phase B: assign curated T states (cluster review 2026-11-26) and run pre-specified tests (PREREG_phaseB.md)."""
import os
import anndata as ad, numpy as np, pandas as pd
from scipy.stats import spearmanr, norm
W = os.environ['TFINE_ROOT']
A = ad.read_h5ad(W + '/T_cells_phaseA.h5ad', backed='r')
obs = A.obs.copy()
gi = {g: list(A.var_names).index(g) for g in ['CD4', 'CD40LG', 'CD8A', 'CD8B', 'FOXP3', 'IL2RA']}
C = A.layers['counts'] if 'counts' in A.layers else A.X
E = {}
for g, j in gi.items():
    col = C[:, j]
    E[g] = np.asarray(col.toarray() if hasattr(col, 'toarray') else col).ravel()
cd8 = (E['CD8A'] > 0) | (E['CD8B'] > 0)
cd4 = ((E['CD4'] > 0) | (E['CD40LG'] > 0)) & ~cd8
treg_cell = ((E['FOXP3'] > 0) | (E['IL2RA'] > 0)) & ~cd8
cl = obs.leiden_fine.astype(str).values
fixed = {'13': ('CD4', 'Treg'), '17': ('CD4', 'CD4_CXCL13_Tph'), '4': ('CD4', 'CD4_memory'), '16': ('CD4', 'CD4_memory'),
         '15': ('CD4', 'CD4_memory'), '7': ('CD4', 'CD4_memory'), '5': ('CD4', 'CD4_memory'),
         '1': ('CD8', 'Tex'), '8': ('CD8', 'CD8_GZMK_Tem'), '3': ('CD8', 'CD8_GZMK_Tem'), '2': ('CD8', 'CD8_Tem_mixed'),
         '0': ('CD8', 'CD8_TRM'), '10': ('CD8', 'CD8_Temra'),
         '9': ('innate', 'gdT'), '14': ('innate', 'MAIT'), '11': ('innate', 'MAIT_Th17like')}
lin = np.full(len(cl), 'unknown', dtype=object); st = np.full(len(cl), 'ambiguous', dtype=object)
for k, (l, s) in fixed.items():
    m = cl == k; lin[m] = l; st[m] = s
# mixed clusters: cell-level lineage
for k, s4, s8 in [('12', 'CD4_naive', 'CD8_naive'), ('6', 'CD4_ISG', 'CD8_ISG'), ('18', 'CD4_ambiguous', 'CD8_ambiguous'),
                  ('20', 'CD4_ambiguous', 'CD8_ambiguous')]:
    m = cl == k
    lin[m & cd8] = 'CD8'; st[m & cd8] = s8
    lin[m & cd4] = 'CD4'; st[m & cd4] = s4
# cycling cluster 19: FOXP3/IL2RA+ non-CD8 -> Treg (cycling); CD8 with Tex score>0 -> Tex (cycling); other CD8 -> CD8_cycling_other
m = cl == '19'
tex = obs.sc_Tex.values > 0
lin[m & treg_cell] = 'CD4'; st[m & treg_cell] = 'Treg'
lin[m & cd4 & ~treg_cell] = 'CD4'; st[m & cd4 & ~treg_cell] = 'CD4_cycling_other'
lin[m & cd8 & tex] = 'CD8'; st[m & cd8 & tex] = 'Tex'
lin[m & cd8 & ~tex] = 'CD8'; st[m & cd8 & ~tex] = 'CD8_cycling_other'
obs['lineage'] = lin; obs['state'] = st
obs['cycling'] = obs.sc_Cycling.values > 0.5
obs['stress_flag'] = obs.sc_Stress.values > 1.0
obs[['leiden_fine', 'lineage', 'state', 'cycling', 'stress_flag', 'ct_Low', 'ct_Low_mv', 'subcluster']].to_csv(W + '/cells_phaseB_states.tsv.gz', sep='\t')
print(pd.crosstab(obs.state, obs.lineage).to_string())

# ---- analysis strata ----
ds, tis, tp = obs.dataset.astype(str), obs.tissue.astype(str), obs.timepoint.astype(str)
obs['pat'] = np.where(obs.patient_id.astype(str).isin(['', 'nan', 'None']), obs.sampleID.astype(str), obs.patient_id.astype(str))
strata = {
 'D_GSE286827_pre': (ds == 'GSE286827') & (tp == 'pre'),
 'D_GSE294109_primary': (ds == 'GSE294109') & tis.str.startswith('Primary'),
 'D_GSE305078_unit': (ds == 'GSE305078'),
 'D_GSE256326': (ds == 'GSE256326'),
 'R_GSE169246_tumour_pre': (ds == 'GSE169246') & (tp == 'Pre') & (tis != 'blood'),
 'R_GSE176021_tumour': (ds == 'GSE176021') & (tis == 'tumor'),
 'C_GSE169246_blood_pre': (ds == 'GSE169246') & (tp == 'Pre') & (tis == 'blood'),
 'C_GSE176021_blood': (ds == 'GSE176021') & (tis == 'blood'),
 'C_GSE176021_normal_lung': (ds == 'GSE176021') & (tis == 'adjacent normal lung'),
 'C_GSE176021_LN': (ds == 'GSE176021') & (tis == 'LN'),
 'X_GSE286827_post': (ds == 'GSE286827') & (tp == 'post'),
}
unit_key = lambda name: 'sampleID' if name == 'D_GSE305078_unit' else 'pat'
cd8_states = ['Tex', 'CD8_GZMK_Tem', 'CD8_Tem_mixed', 'CD8_TRM', 'CD8_Temra', 'CD8_naive']
rows = []
for name, m in strata.items():
    o = obs[m.values]
    if len(o) == 0: continue
    u = o[unit_key(name)].astype(str)
    c4 = o[o.lineage == 'CD4'].groupby(u[o.lineage == 'CD4']).state.value_counts().unstack(fill_value=0)
    c8 = o[o.lineage == 'CD8'].groupby(u[o.lineage == 'CD8']).state.value_counts().unstack(fill_value=0)
    ok = c4.index.intersection(c8.index)
    ok = [x for x in ok if c4.loc[x].sum() >= 50 and c8.loc[x].sum() >= 50]
    if len(ok) < 5:
        rows.append(dict(stratum=name, n=len(ok))); continue
    f4 = c4.loc[ok].div(c4.loc[ok].sum(1), axis=0); f8 = c8.loc[ok].div(c8.loc[ok].sum(1), axis=0)
    for s in cd8_states + ['CD4_CXCL13_Tph']:
        if s not in f8 and s not in f4: continue
    treg = f4.get('Treg', pd.Series(0, index=ok))
    tph = f4.get('CD4_CXCL13_Tph', pd.Series(0, index=ok))
    res = dict(stratum=name, n=len(ok), median_Treg_CD4=round(treg.median(), 3), median_Tex_CD8=round(f8.get('Tex', pd.Series(0, index=ok)).median(), 3))
    for s in cd8_states:
        r, p = spearmanr(treg, f8.get(s, pd.Series(0, index=ok)))
        res['Treg~' + s] = round(r, 2); res['p_Treg~' + s] = p
    r, p = spearmanr(tph, f8.get('Tex', pd.Series(0, index=ok))); res['Tph~Tex'] = round(r, 2); res['p_Tph~Tex'] = p
    r, p = spearmanr(treg, tph); res['Treg~Tph'] = round(r, 2)
    # partial Spearman: rank residuals on stress fraction and log T count
    sf = o.groupby(u).stress_flag.mean().loc[ok]; nT = np.log(o.groupby(u).size().loc[ok])
    Xc = np.column_stack([np.ones(len(ok)), pd.Series(sf).rank(), pd.Series(nT).rank()])
    res_ = lambda y: pd.Series(y).rank().values - Xc @ np.linalg.lstsq(Xc, pd.Series(y).rank().values, rcond=None)[0]
    from scipy.stats import pearsonr
    pr = pearsonr(res_(treg.values), res_(f8.get('Tex', pd.Series(0, index=ok)).values))
    res['partial_Treg~Tex'] = round(pr[0], 2); res['p_partial'] = pr[1]
    # old-node comparison (within-T fraction of old CTLA4 vs GZMB nodes)
    old = o.groupby(u).subcluster.value_counts(normalize=True).unstack(fill_value=0).loc[ok]
    if 'T_CTLA4_235' in old and 'T_GZMB_231' in old:
        res['old_CTLA4~GZMB'] = round(spearmanr(old['T_CTLA4_235'], old['T_GZMB_231'])[0], 2)
    rows.append(res)
    pd.concat([f4.add_prefix('CD4frac_'), f8.add_prefix('CD8frac_')], axis=1).to_csv(W + f'/fractions_{name}.tsv', sep='\t')
R = pd.DataFrame(rows)
R.to_csv(W + '/phaseB_results.tsv', sep='\t', index=False)
pd.set_option('display.width', 400); pd.set_option('display.max_columns', 60)
print(R.to_string())
# pooled replication (Fisher z)
rep = R[R.stratum.str.startswith('R_')]
z = np.arctanh(rep['Treg~Tex'].clip(-0.999, 0.999)); w = rep.n - 3
zp = (z * w).sum() / w.sum(); se = 1 / np.sqrt(w.sum())
print('pooled replication rho', round(np.tanh(zp), 3), 'CI', round(np.tanh(zp - 1.96 * se), 3), round(np.tanh(zp + 1.96 * se), 3))
