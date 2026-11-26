"""Preprint v1 figures 1, 3, 4, 6 and S-bulk from frozen tables (no model fitting except the toy simulation in Fig 3a)."""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt, json, os
from scipy.stats import spearmanr
from sklearn.decomposition import NMF
mpl.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.titlesize': 7.5, 'axes.titleweight': 'bold', 'axes.titlelocation': 'left', 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
W0 = os.environ['TFINE_ROOT']; F = W0 + '/drafts/figs_v1'; I = F + '/inputs'
C = dict(treg='#2a6f97', tex='#c8553d', tph='#6a994e', grey='#9a9a9a', dark='#333333', light='#d9d9d9', acc='#e09f3e')
def save(fig, name):
    for ext in ['png', 'pdf', 'svg']: fig.savefig(f'{F}/{name}.{ext}', dpi=300 if ext == 'png' else None, bbox_inches='tight')
    plt.close(fig)
def lab(ax, s): ax.annotate(s, xy=(0, 1), xycoords='axes fraction', xytext=(-16, 5), textcoords='offset points', fontsize=10, fontweight='bold', va='bottom', ha='right')

# ---------------- Figure 1: release lineage ----------------
fig, ax = plt.subplots(figsize=(8.2, 3.6)); ax.axis('off'); ax.set_xlim(0, 10.2); ax.set_ylim(0, 6)
boxes = [
 (0.2, 4.2, 'Discovery collection (0906)', '54 GEO series · 43 cancer labels\n1,322 samples · 4,676,787 cells\nper-sample QC, scIB comparison,\nV7 tiered annotation, per-sample CNV'),
 (0.2, 1.4, 'Integration v2', '72 sources · 28 cancers · 1,014 patients\n4,720,931 cells · decontX ambient correction\nlineage-sorted captures excluded'),
 (3.6, 1.4, 'TME v1', 'malignancy calls (scATOMIC + inferCNV)\nnon-malignant cells retained\n3,347,659 cells · Harmony (source, unit)'),
 (7.0, 1.4, 'R4 resolved view', '1,974,309 cells · 44 identities\n70 sources · confidence tiers A/B/B2/C1\nCoVarNet K=15'),
 (3.6, 4.2, 'Immune atlas v1', '9 collections · 400 libraries\n2,670,057 cells · capture scopes kept separate'),
 (7.0, 4.2, 'Curated T states', '276,787 T cells (6 sources)\nstress/cycle genes excluded from identity\nCellTypist cross-check · 12 CD4/CD8 states')]
for x, y, t, s in boxes:
    ax.add_patch(mpl.patches.FancyBboxPatch((x, y), 2.9, 1.5, boxstyle='round,pad=0.03,rounding_size=0.08', fc='#f4f4f4', ec=C['dark'], lw=0.7))
    ax.text(x + 0.1, y + 1.33, t, fontweight='bold', va='top', fontsize=7.5); ax.text(x + 0.1, y + 1.0, s, va='top', fontsize=5.6, linespacing=1.35)
for (x1, y1, x2, y2) in [(3.1, 2.15, 3.6, 2.15), (6.5, 2.15, 7.0, 2.15), (6.5, 4.95, 7.0, 4.95)]:
    ax.annotate('', (x2, y2), (x1, y1), arrowprops=dict(arrowstyle='->', lw=0.8, color=C['dark']))
ax.annotate('', (1.65, 2.9), (1.65, 4.2), arrowprops=dict(arrowstyle='->', lw=0.8, color=C['dark'], ls='--'))
ax.text(1.75, 3.45, 're-ingested and re-QCed with new sources;\n36 accessions shared with R4; not additive', fontsize=5.6, color=C['dark'])
ax.text(0.2, 0.75, 'Cell totals are reported per release and are not additive. Module identifiers (CM numbers) are valid only within their release.', fontsize=6, style='italic')
ax.text(0.2, 0.35, 'Discovery → audit (Fig. 3) → R4 module catalogue (Fig. 4) → vascular positive control (Fig. 5) → Tex-centred T-cell module (Fig. 6)', fontsize=6)
save(fig, 'Fig1_release_lineage')

# ---------------- Figure 3: audit ----------------
fig = plt.figure(figsize=(7.2, 5.2)); gs = fig.add_gridspec(2, 6, hspace=0.9, wspace=1.6)
# 3a toy simulation (minmax amplifies rare nodes)
rng = np.random.default_rng(0); shares = {'raw': [], 'minmax': []}
for rep in range(10):
    H = rng.gamma(2, 1, (3, 200)); Wt = np.zeros((12, 3))
    for k in range(3): Wt[3 * k:3 * k + 3, k] = rng.uniform(0.5, 1.5, 3)
    X = Wt @ H + rng.gamma(1, 0.05, (12, 200)); X[9:] = rng.gamma(0.3, 0.02, (3, 200)) * (rng.random((3, 200)) < 0.15)
    for mode in ['raw', 'minmax']:
        Y = X if mode == 'raw' else (X - X.min(1, keepdims=True)) / (X.max(1, keepdims=True) - X.min(1, keepdims=True))
        Wf = NMF(3, beta_loss='kullback-leibler', solver='mu', init='nndsvda', max_iter=2000, random_state=rep).fit_transform(Y)
        Wf = Wf / Wf.sum(0); shares[mode].append(Wf[9:].sum() / Wf.sum())
ax = fig.add_subplot(gs[0, 0:2]); lab(ax, 'a')
ax.bar([0, 1], [np.mean(shares['raw']) * 100, np.mean(shares['minmax']) * 100], color=[C['light'], C['acc']], width=0.6)
for i, k in enumerate(['raw', 'minmax']): ax.scatter(np.full(10, i) + rng.uniform(-0.12, 0.12, 10), np.array(shares[k]) * 100, s=5, color=C['dark'], zorder=3)
ax.set_xticks([0, 1], ['raw frequency', 'per-node min–max']); ax.set_ylabel('rare-node W share (%)'); ax.set_title('Scaling amplifies rare nodes\n(simulation, 10 replicates)')
json.dump({k: [float(x) for x in v] for k, v in shares.items()}, open(I + '/sim_rare_shares.json', 'w'))
# 3b rank surveys
ax = fig.add_subplot(gs[0, 2:4]); lab(ax, 'b')
rs = {'original I-T (15 marker nodes)': (I + '/orig_rank_survey.tsv', C['grey']), 'curated, original units': (W0 + '/covarnet_rerun/A_original_units/discovery_v1/rank_survey.tsv', C['tex']),
      'curated, baseline only': (W0 + '/covarnet_rerun/B_baseline/discovery_v1/rank_survey.tsv', C['treg']), 'curated, independent cohorts': (W0 + '/covarnet_rerun/R_replication/discovery_v1/rank_survey.tsv', C['tph'])}
for name, (p, col) in rs.items():
    r = pd.read_csv(p, sep='\t'); ax.plot(r['rank'], r.cophenetic, '-o', ms=2.5, lw=0.9, color=col, label=name)
    d = np.where(np.diff(r.cophenetic) < -0.02)[0]; k = r['rank'].iloc[d[0]] if len(d) else r['rank'].max()
    ax.scatter([k], [r.cophenetic[r['rank'] == k].iloc[0]], s=40, facecolor='none', edgecolor=col, lw=1.1, zorder=4)
ax.set_xlabel('rank K'); ax.set_ylabel('cophenetic correlation'); ax.set_title('“First drop” rank rule\n(circle = selected K)'); ax.legend(frameon=False, fontsize=4.8, loc='upper center', bbox_to_anchor=(0.5, -0.28), ncol=1)
# 3c original CM02: W vs network membership
ax = fig.add_subplot(gs[0, 4:6]); lab(ax, 'c')
Wo = pd.read_csv(I + '/orig_W.tsv', sep='\t', index_col=0); mem = pd.read_csv(I + '/orig_members.tsv', sep='\t')
w = (Wo['CM02'] / Wo['CM02'].sum()).sort_values(ascending=True); m2 = set(mem[mem.cm == 'CM02'].subCluster)
ax.barh(range(len(w)), w.values, color=[C['acc'] if n in m2 else C['light'] for n in w.index], height=0.7)
ax.set_yticks(range(len(w)), [n[2:] if n.startswith('T_') else n for n in w.index], fontsize=5.5); ax.set_xlabel('W share in original I-T CM02')
ax.set_title('Network members ≠ loaded nodes\n(orange = 7 network members)')
# 3d annotation: original node -> curated state
st = pd.read_csv(W0 + '/cells_phaseB_states.tsv.gz', sep='\t', index_col=0); ob = pd.read_csv(W0 + '/cells_phaseA.tsv.gz', sep='\t', index_col=0, usecols=['row', 'cohort'])
st = st.join(ob); st = st[st.cohort == 'discovery']
grp = {'Treg': 'Treg', 'CD4_CXCL13_Tph': 'CXCL13+ CD4 (Tph)', 'Tex': 'Tex', 'CD4_memory': 'CD4 memory/naive', 'CD4_naive': 'CD4 memory/naive'}
st['g'] = st.state.map(grp).fillna(pd.Series(np.where(st.lineage == 'CD8', 'other CD8', 'other/ambiguous'), index=st.index))
nodes = ['T_CTLA4_235', 'T_ATOX1_238', 'T_GZMB_231', 'T_TYMS_233', 'T_DNAJB1_236', 'T_CCL4_234']
tab = pd.crosstab(st.subcluster, st.g, normalize='index').reindex(nodes)
order = ['Treg', 'CXCL13+ CD4 (Tph)', 'Tex', 'CD4 memory/naive', 'other CD8', 'other/ambiguous']; tab = tab.reindex(columns=order).fillna(0)
cols = [C['treg'], C['tph'], C['tex'], '#a3c4dc', '#e8b4a8', C['light']]
ax = fig.add_subplot(gs[1, 0:3]); lab(ax, 'd'); left = np.zeros(len(tab))
for c, col in zip(order, cols):
    ax.barh(range(len(tab)), tab[c].values, left=left, color=col, label=c, height=0.7); left += tab[c].values
ax.set_yticks(range(len(tab)), [(n[2:] if n.startswith('T_') else n) + ' node' for n in tab.index]); ax.invert_yaxis(); ax.set_xlabel('fraction of node cells by curated state')
ax.set_title('Original marker nodes mix curated states'); ax.legend(frameon=False, fontsize=5, ncol=3, loc='upper center', bbox_to_anchor=(0.5, -0.28))
tab.to_csv(I + '/node_to_state.tsv', sep='\t')
# 3e LOSO
ax = fig.add_subplot(gs[1, 3:6]); lab(ax, 'e')
lo = pd.read_csv(I + '/orig_LOSO.tsv', sep='\t'); lo = lo[lo.cm == 'CM02'].set_index('heldout_source').best_cosine
la = pd.read_csv(W0 + '/covarnet_rerun/A_original_units/discovery_v1/LOSO_stability.tsv', sep='\t'); la = la[la.cm == 'CM01'].set_index('heldout_source').best_cosine
lb = pd.read_csv(W0 + '/covarnet_rerun/B_baseline/discovery_v1/LOSO_stability.tsv', sep='\t'); lb = lb[lb.cm == 'CM01'].set_index('heldout_source').best_cosine
src = list(lo.index); x = np.arange(len(src))
for k, (s, col, lab_) in enumerate([(lo, C['grey'], 'original CM02 (marker nodes)'), (la, C['tex'], 'curated module, original units'), (lb, C['treg'], 'curated module, baseline only')]):
    ax.bar(x + (k - 1) * 0.26, s.reindex(src).values, width=0.25, color=col, label=lab_)
ax.set_xticks(x, src, rotation=20); ax.set_ylim(0.5, 1.01); ax.set_ylabel('best-matched cosine'); ax.set_title('Leave-one-source-out refits\n(x = removed source)'); ax.legend(frameon=False, fontsize=5, loc='upper center', bbox_to_anchor=(0.5, -0.3), ncol=1)
save(fig, 'Fig3_module_discovery_audit')

# ---------------- Figure 4: R4 15-module catalogue ----------------
m = pd.read_csv(I + '/R4_member_comparison_eligibility.tsv', sep='\t')
def lin(n):
    for p, l in [('CD4T', 'CD4 T'), ('CD8T', 'CD8 T'), ('ILC', 'NK/ILC'), ('Mye', 'Myeloid'), ('B_', 'B/plasma'), ('Plasma', 'B/plasma'), ('Endo', 'Endothelial'), ('Str', 'Stromal')]:
        if n.startswith(p): return l
    return 'other'
m['lin'] = m.node.map(lin)
cat = pd.crosstab(m.cm, m.lin); lins = ['CD4 T', 'CD8 T', 'NK/ILC', 'Myeloid', 'B/plasma', 'Endothelial', 'Stromal']
cat = cat.reindex(columns=lins).fillna(0)
info = m.groupby('cm').agg(high_sources=('high_sources', 'first'), largest=('largest_high_source_fraction', 'first'),
                           triad=('node', lambda s: all(x in set(s) for x in ['CD4T_Tfh-like_CXCL13', 'CD4T_Treg_FOXP3', 'CD8T_Tex_LAG3'])),
                           nonT=('node', lambda s: ', '.join(x.split('_', 1)[1] for x in s if not x.startswith(('CD4T', 'CD8T', 'ILC')))))
fig = plt.figure(figsize=(7.2, 4.6)); gs = fig.add_gridspec(1, 3, width_ratios=[2.2, 0.9, 2.6], wspace=0.08)
ax = fig.add_subplot(gs[0]); lab(ax, 'a'); left = np.zeros(15)
pal = [C['tph'], C['tex'], '#b07aa1', C['acc'], '#4e79a7', '#d37295', '#8c6d31']
for l, col in zip(lins, pal): ax.barh(range(15), cat[l].values, left=left, color=col, label=l, height=0.7); left += cat[l].values
ax.set_yticks(range(15), cat.index); ax.invert_yaxis(); ax.set_xlabel('network members (count)'); ax.legend(frameon=False, fontsize=5, ncol=4, loc='upper left', bbox_to_anchor=(0, -0.12)); ax.set_title('Member lineages')
ax2 = fig.add_subplot(gs[1], sharey=ax); ax2.barh(range(15), info.high_sources.values, color=C['grey'], height=0.7); ax2.set_xlabel('sources with\nhigh-activity units'); ax2.tick_params(labelleft=False); ax2.set_title('Coverage')
for i, (h, l) in enumerate(zip(info.high_sources, info.largest)): ax2.text(h + 0.5, i, f'{l:.2f}', va='center', fontsize=5)
ax3 = fig.add_subplot(gs[2], sharey=ax); ax3.axis('off'); ax3.set_title('Non-lymphoid members (● = Tph, Treg and Tex all members)', loc='left')
for i, (cm, r) in enumerate(info.iterrows()):
    ax3.text(0.0, i, ('● ' if r.triad else '○ ') + (r.nonT if r.nonT else '—'), fontsize=5.5, va='center', transform=ax3.get_yaxis_transform())
fig.text(0.13, -0.08, 'Numbers next to coverage bars: largest single-source fraction among high-activity units. Membership uses the original CoVarNet top-N network filter (see Fig. 3c);\nTreg and Tex are members of all 15 modules, so membership alone does not localise them to one module.', fontsize=5.3)
save(fig, 'Fig4_R4_module_catalogue'); info.join(cat).to_csv(I + '/R4_catalogue.tsv', sep='\t')

# ---------------- Figure 6: Tex-centred module ----------------
P = pd.read_csv(W0 + '/phaseB_results.tsv', sep='\t')
lbl = {'D_GSE286827_pre': 'HNSCC GSE286827 pre (26 pts)', 'D_GSE294109_primary': 'ccRCC GSE294109 (13 pts)', 'D_GSE305078_unit': 'TNBC GSE305078 (32 units)',
       'R_GSE169246_tumour_pre': 'TNBC GSE169246 pre (13 pts) †', 'R_GSE176021_tumour': 'NSCLC GSE176021 (15 pts) †',
       'C_GSE176021_normal_lung': 'Adjacent normal lung (12 pts)', 'C_GSE169246_blood_pre': 'TNBC blood pre (19 pts)', 'X_GSE286827_post': 'HNSCC GSE286827 post (28 pts)'}
fig = plt.figure(figsize=(7.4, 6.2)); gs = fig.add_gridspec(2, 3, hspace=0.6, wspace=0.95)
ax = fig.add_subplot(gs[0, 0:2]); lab(ax, 'a'); keys = list(lbl)
for i, k in enumerate(keys):
    r = P[P.stratum == k].iloc[0]; n = r.n; z = np.arctanh(r['Treg~Tex']); se = 1 / np.sqrt(n - 3)
    lo_, hi_ = np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)
    col = C['tex'] if k.startswith(('D_', 'R_')) else C['grey']
    ax.plot([lo_, hi_], [i, i], color=col, lw=1); ax.scatter([r['Treg~Tex']], [i], color=col, s=14, zorder=3, marker='D' if k.startswith('R_') else 'o')
    ax.scatter([r['Tph~Tex']], [i], color=C['tph'], s=10, marker='|', zorder=3)
ax.axvline(0, color=C['dark'], lw=0.5, ls='--'); ax.set_yticks(range(len(keys)), [lbl[k] for k in keys]); ax.invert_yaxis()
ax.set_xlabel('Spearman ρ, Treg/CD4 vs Tex/CD8 (Fisher-z 95% CI); green tick = Tph/CD4 vs Tex/CD8'); ax.set_title('Patient-level coupling with lineage-internal denominators')
ax.text(0.01, 0.02, '† not used for module discovery; pooled ρ 0.53 (0.17–0.76)', transform=ax.transAxes, fontsize=5.3)
ax = fig.add_subplot(gs[0, 2]); lab(ax, 'b')
states = ['Treg', 'CD4_CXCL13_Tph', 'Tex', 'CD8_TRM', 'CD4_ISG', 'CD8_ISG', 'CD8_GZMK_Tem', 'CD8_Tem_mixed', 'CD8_Temra', 'CD4_memory', 'CD4_naive', 'CD8_naive']
Wm = {}
for d, cm, nm in [('A_original_units', 'CM01', 'A'), ('B_baseline', 'CM01', 'B'), ('R_replication', 'CM02', 'R')]:
    w = pd.read_csv(f'{W0}/covarnet_rerun/{d}/discovery_v1/W_private.tsv', sep='\t', index_col=0)[cm]; Wm[nm] = (w / w.sum()).reindex(states)
Wm = pd.DataFrame(Wm); im = ax.imshow(Wm.values, cmap='Reds', vmin=0, vmax=0.45, aspect='auto')
ax.set_yticks(range(len(states)), [s.replace('CD4_CXCL13_Tph', 'Tph (CXCL13+ CD4)').replace('_', ' ') for s in states], fontsize=5.5); ax.set_xticks([0, 1, 2], ['A', 'B', 'R'], fontsize=6)
for i in range(len(states)):
    for j in range(3): ax.text(j, i, f'{Wm.values[i, j]:.2f}', ha='center', va='center', fontsize=4.8, color='white' if Wm.values[i, j] > 0.3 else 'black')
ax.set_title('Module W share\nA original units · B baseline · R independent', fontsize=6.5); fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02)
ax = fig.add_subplot(gs[1, 0]); lab(ax, 'c'); ax.axis('off'); ax.set_xlim(-1.4, 1.4); ax.set_ylim(-1.3, 1.5)
pos = {'Tex': (0, 1.0), 'Treg': (-1, -0.6), 'Tph': (1, -0.6)}
for a, b, txt, ls in [('Tex', 'Treg', 'r 0.38–0.43 (A/B)\n0.56 (R)', '-'), ('Tex', 'Tph', 'r 0.37–0.43 (A/B)\n0.62 (R)', '-'), ('Treg', 'Tph', '\n\nr −0.03 to 0.16, no edge', ':')]:
    (x1, y1), (x2, y2) = pos[a], pos[b]; ax.plot([x1, x2], [y1, y2], ls=ls, color=C['dark'], lw=1.2 if ls == '-' else 0.8)
    ax.text((x1 + x2) / 2 + (0.1 if a == 'Tex' and b == 'Tph' else -0.1 if a == 'Tex' else 0), (y1 + y2) / 2 + (-0.25 if b == 'Tph' and a == 'Treg' else 0), txt, fontsize=5.3, ha='left' if (a, b) == ('Tex', 'Tph') else 'right' if (a, b) == ('Tex', 'Treg') else 'center')
for k, (x, y) in pos.items():
    ax.scatter([x], [y], s=520, color={'Tex': C['tex'], 'Treg': C['treg'], 'Tph': C['tph']}[k], zorder=3); ax.text(x, y, k, ha='center', va='center', color='white', fontweight='bold', fontsize=7, zorder=4)
ax.set_title('Source-adjusted edges (curated)', fontsize=7.5)
ax = fig.add_subplot(gs[1, 1]); lab(ax, 'd'); ax.axis('off'); ax.set_xlim(-1.4, 1.4); ax.set_ylim(-1.3, 1.5)
pos = {'Tex\n(LAG3)': (0, 1.0), 'Treg\n(FOXP3)': (-1, -0.6), 'Tfh-like\n(CXCL13)': (1, -0.6), 'NK\n56bright': (0, -0.1)}
for a, b in [('Tex\n(LAG3)', 'Treg\n(FOXP3)'), ('Tex\n(LAG3)', 'Tfh-like\n(CXCL13)'), ('NK\n56bright', 'Tex\n(LAG3)'), ('NK\n56bright', 'Treg\n(FOXP3)'), ('NK\n56bright', 'Tfh-like\n(CXCL13)')]:
    (x1, y1), (x2, y2) = pos[a], pos[b]; ax.plot([x1, x2], [y1, y2], color=C['dark'], lw=1)
(x1, y1), (x2, y2) = pos['Treg\n(FOXP3)'], pos['Tfh-like\n(CXCL13)']; ax.plot([x1, x2], [y1, y2], ls=':', color=C['grey'], lw=0.8); ax.text(0, -0.85, 'no edge', fontsize=5.3, ha='center', color=C['grey'])
for k, (x, y) in pos.items():
    col = C['tex'] if 'Tex' in k else C['treg'] if 'Treg' in k else C['tph'] if 'Tfh' in k else '#b07aa1'
    ax.scatter([x], [y], s=520, color=col, zorder=3); ax.text(x, y, k, ha='center', va='center', color='white', fontweight='bold', fontsize=5, zorder=4)
ax.set_title('R4 network (separate atlas)', fontsize=7.5)
ax = fig.add_subplot(gs[1, 2]); lab(ax, 'e')
sc = pd.read_csv(W0 + '/bulk_v2/sc_pre_post_patient_fractions.tsv', sep='\t', index_col=[0, 1])
h = sc.loc['GSE286827']
for j, (k, col) in enumerate([('Treg_CD4', C['treg']), ('Tph_CD4', C['tph']), ('Tex_CD8', C['tex'])]):
    for _, r in h.iterrows(): ax.plot([j * 2.4, j * 2.4 + 1], [r[k + '_pre'], r[k + '_post']], color=col, lw=0.4, alpha=0.5)
    ax.scatter([j * 2.4, j * 2.4 + 1], [h[k + '_pre'].median(), h[k + '_post'].median()], color=col, s=18, zorder=3, edgecolor='black', lw=0.4)
ax.set_xticks([0, 1, 2.4, 3.4, 4.8, 5.8], ['pre', 'post'] * 3, fontsize=5.5); ax.set_ylabel('fraction within lineage')
for j, (t, p) in enumerate([('Treg/CD4', 'P=0.62'), ('Tph/CD4', 'P=0.022'), ('Tex/CD8', 'P=0.012')]): ax.text(j * 2.4 + 0.5, 1.02, f'{t}\n{p}', ha='center', fontsize=5.3, transform=ax.get_xaxis_transform())
ax.set_title('HNSCC after ICB (26 paired pts)\nexploratory', fontsize=7.5, pad=18)
save(fig, 'Fig6_Tex_centred_module')

# ---------------- Supplementary: bulk v2 ----------------
B = pd.read_csv(W0 + '/bulk_v2/bulk_v2_all_models.tsv', sep='\t'); B = B[(B.status == 'ok') & (B.endpoint == 'response')]
fig, ax = plt.subplots(figsize=(6.0, 3.4)); rows = []
for co in ['Liu2019', 'Gide2019', 'IMvigor210', 'SU2C_MARK', 'Riaz2017', 'Hugo2016']:
    for rd, col in [('R1_Treg_Tex', C['tex']), ('arm_Treg', C['treg']), ('GEP', C['grey']), ('v1_pCM02', C['acc'])]:
        r = B[(B.cohort == co) & (B.readout == rd)]
        if len(r): rows.append((co, rd, col, *r[['est', 'lo', 'hi']].iloc[0]))
for i, (co, rd, col, e, l, h_) in enumerate(rows):
    ax.plot([l, h_], [i, i], color=col, lw=0.9); ax.scatter([e], [i], color=col, s=10, zorder=3)
ax.set_xscale('log'); ax.axvline(1, color=C['dark'], lw=0.5, ls='--'); ax.set_yticks(range(len(rows)), [f'{co} · {rd}' for co, _, rd, *_ in [(r[0], None, r[1]) for r in rows]], fontsize=5)
ax.invert_yaxis(); ax.set_xlabel('OR per SD for CR/PR (95% CI); R1, Treg arm and v1 adjusted for pan-T'); ax.set_title('Bulk readouts: Liu and Gide are held-out; others previously analysed')
save(fig, 'FigS_bulk_readout_v2')
print('figures written', sorted(os.listdir(F)))
