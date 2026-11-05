"""Join V7 CoVarNet CM membership with per-section blisa hotspots/CCI.
Per section x CM: mean CCI over ordered within-CM identity pairs vs a lineage-matched random-member null."""
import sys, json, re
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import binomtest

run = Path(sys.argv[1]); branch = sys.argv[2]; seed = sys.argv[3]; bin_um = int(sys.argv[4]); membership = sys.argv[5]
MIN_CELLS, NPERM, ALPHA = 20, 999, 0.05
rng = np.random.default_rng(1337)
cmdir = run / f'outputs/covarnet_{branch}/seed_{seed}'
PREFIX = re.compile(r'^(Malignant|Malignancy_uncertain|Malignancy_not_evaluable)::')
if membership == 'network':
    mem = pd.read_csv(cmdir / 'author_compatibility_each_node_v1.tsv', sep='\t')
else:  # W top10
    W = pd.read_csv(cmdir / 'private_W_v1.tsv', sep='\t', index_col=0)
    mem = pd.concat([pd.DataFrame({'cm': c, 'subCluster': W[c].nlargest(10).index}) for c in W.columns])
mem['identity'] = mem.subCluster.str.replace(PREFIX, '', regex=True)
mem['prefix'] = mem.subCluster.str.extract(r'^(.*)::')[0].fillna('none')
ref_major = pd.read_csv(run / 'xenium/label_major.tsv', sep='\t', index_col=0).major

def bh_count(sig_p, m):
    """Exact BH rejections at ALPHA from the p<=0.05 subset (larger p cannot be rejected)."""
    p = np.sort(np.concatenate(sig_p)) if len(sig_p) else np.array([])
    if not len(p): return 0, 0.0
    k = np.nonzero(p <= ALPHA * np.arange(1, len(p) + 1) / m)[0]
    return (k.max() + 1 if len(k) else 0), (p[k.max()] if len(k) else 0.0)

rows, lr_rows = [], []
for d in sorted((run / 'xenium').iterdir()):
    b = d / f'blisa_bin{bin_um}'
    if not (b / 'summary.tsv').exists(): continue
    summ = pd.read_csv(b / 'summary.tsv', sep='\t').iloc[0]
    cells = pd.read_csv(d / 'cells.tsv.gz', sep='\t', usecols=['group', 'identifiable'])
    n_by = cells.group.value_counts(); ok_groups = set(n_by.index[(n_by >= MIN_CELLS)]) - {g for g in n_by.index if g.startswith('Unresolved_')}
    lr = pd.read_csv(b / 'LR_results.tsv.gz', sep='\t')
    m = max(1, len(lr) * int(summ.n_bins))
    sigp = [np.array([float(x) for x in str(s).split(';') if x not in ('', 'nan')]) for s in lr.sig_pval.fillna('')]
    k, thr = bh_count(sigp, m)
    lr['n_BH_hotspot_bins'] = [int((p <= thr).sum()) if k else 0 for p in sigp]
    lr_rows.append(lr[['interaction_name', 'ligand', 'receptor', 'annotation', 'ccc_mode', 'sig_numbers', 'n_BH_hotspot_bins']].assign(section=d.name, n_bins=int(summ.n_bins), m_tests=m))
    cci = pd.read_csv(b / 'CCI_scores.tsv.gz', sep='\t')
    lrcols = [c for c in cci.columns if c not in ('Sender', 'Receiver')]
    cci['score_all'] = cci[lrcols].sum(1)
    bh_lr = set(lr.interaction_name[lr.n_BH_hotspot_bins > 0])
    cci['score_BH'] = cci[[c for c in lrcols if c in bh_lr]].sum(1) if bh_lr else 0.0
    S = cci.set_index(['Sender', 'Receiver'])
    groups = sorted(ok_groups & set(cci.Sender) & set(cci.Receiver))
    gmajor = pd.Series({g: ref_major.get(g, g.replace('Other_', '')) for g in groups})
    def pair_mean(members, col):
        pr = [(s, r) for s in members for r in members if s != r and (s, r) in S.index]
        return S.loc[pr, col].mean() if pr else np.nan
    for cm, g in mem.groupby('cm'):
        present = sorted(set(g.identity) & set(groups))
        rec = dict(section=d.name, bin_um=bin_um, cm=cm, n_members=len(set(g.identity)), n_present_identifiable=len(present),
                   members_present=';'.join(present), n_BH_LR=len(bh_lr), n_LR=len(lr))
        if len(present) >= 2:
            for col in ('score_all', 'score_BH'):
                obs = pair_mean(present, col); null = []
                for _ in range(NPERM):
                    pick = []
                    for mbr in present:  # lineage-matched: same V7 major, not already picked
                        pool = [x for x in groups if gmajor[x] == gmajor[mbr] and x not in pick] or [x for x in groups if x not in pick]
                        pick.append(rng.choice(pool))
                    null.append(pair_mean(pick, col))
                null = np.array(null, float)
                rec.update({f'{col}_obs': obs, f'{col}_null_mean': np.nanmean(null), f'{col}_null_sd': np.nanstd(null),
                            f'{col}_z': (obs - np.nanmean(null)) / np.nanstd(null) if np.nanstd(null) > 0 else np.nan,
                            f'{col}_p_upper': (1 + np.sum(null >= obs)) / (1 + np.sum(np.isfinite(null)))})
            top = S.loc[[(s, r) for s in present for r in present if s != r and (s, r) in S.index], lrcols].sum().nlargest(5)
            rec['top_LR_within_CM'] = ';'.join(f'{i}:{v:.2f}' for i, v in top.items() if v > 0)
        rows.append(rec)
out = run / f'results/cm_blisa_{branch}_seed{seed}_{membership}_bin{bin_um}'; out.mkdir(parents=True, exist_ok=True)
res = pd.DataFrame(rows); res.to_csv(out / 'section_by_CM.tsv', sep='\t', index=False)
pd.concat(lr_rows).to_csv(out / 'LR_hotspots_by_section.tsv.gz', sep='\t', index=False)
# patient units: S18 carries patient ids; S07 reps are one patient; other sections are their own unit
units = json.load(open(run / 'config/xenium_patient_units.json'))
res['unit'] = res.section.map(units)
summ = []
for cm, g in res.dropna(subset=['score_all_z']).groupby('cm'):
    u = g.groupby('unit').score_all_z.mean(); ub = g.groupby('unit').score_BH_z.mean().dropna()
    summ.append(dict(cm=cm, n_units=len(u), units_z_pos=int((u > 0).sum()), sign_p=binomtest(int((u > 0).sum()), len(u), 0.5, alternative='greater').pvalue if len(u) else np.nan,
                     median_z=u.median(), n_units_BH=len(ub), units_BH_z_pos=int((ub > 0).sum()), median_z_BH=ub.median() if len(ub) else np.nan))
pd.DataFrame(summ).to_csv(out / 'CM_summary_across_units.tsv', sep='\t', index=False)
print(pd.DataFrame(summ).to_string())
