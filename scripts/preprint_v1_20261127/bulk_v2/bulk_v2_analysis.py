"""Bulk readout v2 (PREREG_bulk_v2.md + Amendment A)."""
import os
import json, numpy as np, pandas as pd, statsmodels.api as sm, warnings
from statsmodels.duration.hazard_regression import PHReg
from scipy.stats import wilcoxon, chi2, norm
warnings.filterwarnings('ignore')
B = os.environ['TFINE_ROOT'] + '/bulk_v2'; DD = B + '/data'
SU = os.environ['PROJECT_V5_ROOT'] + '/runs/20261018_priority_immune_CM_bulk_v1'
P = json.load(open(B + '/frozen_arm_panels.json'))
PAN = {'Treg': P['Treg'], 'Tex': P['Tex'], 'panT': ['CD3D', 'CD3E', 'TRAC'], 'CYT': ['GZMA', 'PRF1'],
       'GEP': ['CCL5','CD27','CD274','CD276','CD8A','CMKLR1','CXCL9','CXCR6','HLA-DQA1','HLA-DRB1','HLA-E','IDO1','LAG3','NKG7','PDCD1LG2','PSMB10','STAT1','TIGIT'],
       'v1_reg': ['FOXP3','IL2RA','CTLA4','IKZF2','TNFRSF18'], 'v1_cyt': ['GZMB','PRF1','CTSW','NKG7','CD8A','CD8B']}

def score(E):
    R = E.rank(axis=1, pct=True, method='average'); S = pd.DataFrame(index=E.columns)
    for k, g in PAN.items():
        ok = [x for x in g if x in E.index and E.loc[x].notna().all() and E.loc[x].nunique() > 1]
        S[k] = R.loc[ok].mean() if len(ok) == len(g) else np.nan
    S['R1'] = np.minimum(S.Treg, S.Tex); S['v1_pCM02'] = np.minimum(S.v1_reg, S.v1_cyt)
    return S

def ev(s): return s.astype(str).str.upper().str[0].map({'1': 1, '0': 0, 'D': 1, 'L': 0, 'T': 1, 'F': 0})

def load_cbio(s, keep_pre=True):
    E = pd.read_csv(f'{DD}/{s}_expr.tsv', sep='\t', index_col=0)
    E.columns = E.columns.astype(str)
    S = pd.read_csv(f'{DD}/{s}_clin_SAMPLE.tsv', sep='\t', dtype={'id': str, 'patient': str}).set_index('id')
    Pt = pd.read_csv(f'{DD}/{s}_clin_PATIENT.tsv', sep='\t', dtype={'id': str, 'patient': str}).set_index('patient')
    S = S.loc[S.index.intersection(E.columns)]
    return E, S, Pt

def patient_table(E, S, Pt, phase='Pre'):
    S = S[S.SAMPLE_TREATMENT.astype(str) == phase]
    sc = score(E)  # percentile within cohort across all samples of both phases
    sc = sc.loc[S.index].groupby(S.patient).mean()
    c = pd.DataFrame(index=sc.index)
    resp = S.groupby('patient').RESPONSE.first() if 'RESPONSE' in S else Pt.RESPONSE
    c['response'] = resp.reindex(c.index).map({'Complete Response': 1, 'Partial Response': 1, 'Stable Disease': 0, 'Progressive Disease': 0})
    for ep in ['OS', 'PFS']:
        if f'{ep}_MONTHS' in Pt:
            c[ep] = pd.to_numeric(Pt[f'{ep}_MONTHS'], errors='coerce').reindex(c.index)
            c[ep + '_event'] = ev(Pt[f'{ep}_STATUS']).reindex(c.index)
    for k in ['PRIOR_ICI_RX', 'ICI_RX']:
        src = S.groupby('patient')[k].first() if k in S else (Pt[k] if k in Pt else None)
        if src is not None: c[k] = src.reindex(c.index)
    return sc.join(c)

def load_su2c():
    g = pd.read_csv(DD + '/SU2C-MARK_Harmonized_rnaseqc_tpm_v1.gct', sep='\t', skiprows=2, index_col=1)
    g = g.drop(columns=['Name']); genes = set(sum(PAN.values(), []))
    E = g[g.index.isin(genes)]; E = E.groupby(level=0).mean()
    c = pd.read_excel(SU + '/SU2C_MARK_panel_ready_v1/Table_S1_Clinical_Annotations.xlsx', header=2)
    sid = 'Harmonized_SU2C_RNA_Tumor_Sample_ID_v2'
    cc = c[c.RNA_All.eq(1) & c['Pre-treatment_RNA_Sample_QC'].eq('Keep')].set_index(sid)
    E = E[cc.index.intersection(E.columns)]; sc = score(E); cc = cc.loc[E.columns]
    out = pd.DataFrame(index=E.columns)
    out['response'] = cc.Harmonized_Confirmed_BOR.map({'CR': 1, 'PR': 1, 'SD': 0, 'PD': 0})
    for ep in ['OS', 'PFS']:
        out[ep] = pd.to_numeric(cc['Harmonized_' + ep + '_Days'], errors='coerce')
        e = pd.to_numeric(cc['Harmonized_' + ep + '_Event'], errors='coerce'); out[ep + '_event'] = e.where(e.isin([0, 1]))
    return sc.join(out)

def z(x): return (x - x.mean()) / x.std()

def fit(d, ep, pred, adj):
    cols = [pred] + adj
    if ep == 'response':
        x = d[cols + ['response']].dropna()
        if len(x) < 20 or x.response.sum() < 5 or (1 - x.response).sum() < 5: return dict(status='unavailable', n=len(x))
        X = sm.add_constant(x[cols].apply(z)); m = sm.Logit(x.response, X).fit(disp=0)
        b, se = m.params[pred], m.bse[pred]
        return dict(status='ok', n=len(x), events=int(x.response.sum()), effect='OR', est=np.exp(b), lo=np.exp(b - 1.96 * se), hi=np.exp(b + 1.96 * se), p=m.pvalues[pred], beta=b, se=se, llf=m.llf)
    x = d[cols + [ep, ep + '_event']].dropna(); x = x[x[ep] > 0]
    if len(x) < 20 or x[ep + '_event'].sum() < 10: return dict(status='unavailable', n=len(x))
    m = PHReg(x[ep].values, x[cols].apply(z).values, status=x[ep + '_event'].values, ties='efron').fit()
    b, se = m.params[0], m.bse[0]
    return dict(status='ok', n=len(x), events=int(x[ep + '_event'].sum()), effect='HR', est=np.exp(b), lo=np.exp(b - 1.96 * se), hi=np.exp(b + 1.96 * se), p=m.pvalues[0], beta=b, se=se, llf=m.llf)

cohorts = {}
for s, name, role in [('mel_iatlas_liu_2019', 'Liu2019', 'held_out'), ('mel_iatlas_gide_2019', 'Gide2019', 'held_out'),
                      ('blca_iatlas_imvigor210_2017', 'IMvigor210', 'reanalysis'), ('mel_iatlas_riaz_nivolumab_2017', 'Riaz2017', 'reanalysis'),
                      ('mel_iatlas_hugo_ucla_2016', 'Hugo2016', 'reanalysis')]:
    E, S, Pt = load_cbio(s); cohorts[name] = (role, patient_table(E, S, Pt), (E, S, Pt))
cohorts['SU2C_MARK'] = ('reanalysis', load_su2c(), None)

rows = []
for name, (role, d, _) in cohorts.items():
    for ep in ['response', 'OS', 'PFS']:
        for readout, pred, adj in [('R1_Treg_Tex', 'R1', ['panT']), ('R2_Tph_Tex', None, ['panT']),
                                   ('arm_Treg', 'Treg', ['panT']), ('arm_Tex', 'Tex', ['panT']), ('v1_pCM02', 'v1_pCM02', ['panT']),
                                   ('GEP', 'GEP', []), ('CYT', 'CYT', []), ('R1_adj_GEP', 'R1', ['panT', 'GEP'])]:
            if pred is None: r = dict(status='unavailable_no_bulk_specific_Tph_genes')
            elif ep not in d and ep != 'response': r = dict(status='unavailable_endpoint_absent')
            else: r = fit(d, ep, pred, adj)
            rows.append(dict(cohort=name, role=role, endpoint=ep, readout=readout, **r))
        if name == 'Liu2019' and ep == 'response':
            dd = d.copy(); dd['priorCTLA4'] = dd.PRIOR_ICI_RX.notna().astype(float)
            rows.append(dict(cohort=name, role=role, endpoint=ep, readout='R1_sens_priorCTLA4', **fit(dd, ep, 'R1', ['panT', 'priorCTLA4'])))
        if name == 'Gide2019' and ep == 'response':
            dd = d.copy(); dd['combo'] = dd.ICI_RX.astype(str).str.contains('Ipi').astype(float)
            rows.append(dict(cohort=name, role=role, endpoint=ep, readout='R1_sens_combo', **fit(dd, ep, 'R1', ['panT', 'combo'])))
R = pd.DataFrame(rows)
fam = R.readout.isin(['R1_Treg_Tex', 'R2_Tph_Tex'])
pv = R.loc[fam, 'p'].fillna(1).values; o = np.argsort(pv); n = len(pv)
q = np.minimum.accumulate((pv[o] * n / np.arange(1, n + 1))[::-1])[::-1]; qq = np.empty(n); qq[o] = np.clip(q, 0, 1)
R.loc[fam, 'BH_q_family36'] = qq
R.to_csv(B + '/bulk_v2_all_models.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 300)
show = R[R.status == 'ok'][['cohort', 'role', 'endpoint', 'readout', 'n', 'events', 'effect', 'est', 'lo', 'hi', 'p', 'BH_q_family36']]
print(show.round(3).to_string(index=False))
# H1 pooled (held-out, response)
h = R[(R.role == 'held_out') & (R.endpoint == 'response') & (R.readout == 'R1_Treg_Tex') & (R.status == 'ok')]
w = 1 / h.se ** 2; b = (h.beta * w).sum() / w.sum(); se = 1 / np.sqrt(w.sum()); Q = ((h.beta - b) ** 2 * w).sum()
print(f'\nH1 pooled held-out OR {np.exp(b):.2f} ({np.exp(b-1.96*se):.2f}-{np.exp(b+1.96*se):.2f}) p={2*norm.sf(abs(b/se)):.4f} Cochran Q={Q:.2f} p_het={chi2.sf(Q, len(h)-1):.3f}')
for rd in ['arm_Treg', 'arm_Tex', 'GEP', 'v1_pCM02']:
    h = R[(R.role == 'held_out') & (R.endpoint == 'response') & (R.readout == rd) & (R.status == 'ok')]
    w = 1 / h.se ** 2; b = (h.beta * w).sum() / w.sum(); se = 1 / np.sqrt(w.sum())
    print(f'  descriptive pooled held-out {rd}: OR {np.exp(b):.2f} ({np.exp(b-1.96*se):.2f}-{np.exp(b+1.96*se):.2f})')
# LRT R1 over panT+GEP (response)
for name, (role, d, _) in cohorts.items():
    x = d[['R1', 'panT', 'GEP', 'response']].dropna()
    if len(x) < 20: continue
    m0 = sm.Logit(x.response, sm.add_constant(x[['panT', 'GEP']].apply(z))).fit(disp=0)
    m1 = sm.Logit(x.response, sm.add_constant(x[['R1', 'panT', 'GEP']].apply(z))).fit(disp=0)
    print(f'  LRT R1 | panT+GEP {name}: chi2={2*(m1.llf-m0.llf):.2f} p={chi2.sf(2*(m1.llf-m0.llf),1):.4f}')
# Amendment A2: bulk pre vs on (Riaz, Gide)
out = []
for name in ['Riaz2017', 'Gide2019']:
    E, S, Pt = cohorts[name][2]; sc = score(E)
    pre = S[S.SAMPLE_TREATMENT == 'Pre']; on = S[S.SAMPLE_TREATMENT == 'On']
    a = sc.loc[pre.index].groupby(pre.patient).mean(); b2 = sc.loc[on.index].groupby(on.patient).mean()
    pts = a.index.intersection(b2.index); resp = cohorts[name][1].response.reindex(pts)
    if resp.isna().any():
        rr = S.groupby('patient').RESPONSE.first().reindex(pts).map({'Complete Response': 1, 'Partial Response': 1, 'Stable Disease': 0, 'Progressive Disease': 0}); resp = resp.fillna(rr)
    for k in ['Treg', 'Tex', 'R1', 'panT', 'GEP']:
        dlt = b2.loc[pts, k] - a.loc[pts, k]
        for grp, m in [('all', resp.notna() | resp.isna()), ('responder', resp == 1), ('non_responder', resp == 0)]:
            x = dlt[m.values]
            if len(x) >= 5:
                out.append(dict(cohort=name, score=k, group=grp, n=len(x), median_pre=a.loc[x.index, k].median(), median_on=b2.loc[x.index, k].median(), median_delta=x.median(), p=wilcoxon(x).pvalue))
L = pd.DataFrame(out); L.to_csv(B + '/bulk_pre_on_change.tsv', sep='\t', index=False)
print('\nBulk pre->on change (exploratory):'); print(L.round(3).to_string(index=False))
