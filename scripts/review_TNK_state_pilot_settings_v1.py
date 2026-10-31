#!/usr/bin/env python3
"""Cross-setting review table for the data-driven T/NK state pilot (trial + sensitivity runs).

Read-only aid for human review; it names nothing, freezes nothing and starts no CM.
For every kept meta-program (MP) of every setting it records the tier, the reliability and
confounding diagnostics already computed by the pilot, and membership of fixed marker-gene
families in its top genes. MPs are then linked across settings by top-gene Jaccard: each
reference-setting MP is matched to the best MP of every other setting (match only if
Jaccard >= --min-jaccard, fixed at 0.3 before reading results), so stability across settings
can be read row by row.

Usage: review_TNK_state_pilot_settings_v1.py OUT_DIR REF_NAME=RUN_DIR NAME=RUN_DIR [...]
Writes OUT_DIR/{all_settings_meta_programs.tsv, cross_setting_matches.tsv, review_summary.json}.
"""
import argparse, json, os
import pandas as pd

# Marker families are descriptive flags for review, not state calls.
FAMILIES = {
    'interferon_ISG': {'ISG15', 'IFI6', 'IFI44L', 'IFIT1', 'IFIT2', 'IFIT3', 'IFITM1', 'MX1', 'MX2',
                       'OAS1', 'OAS2', 'OAS3', 'RSAD2', 'STAT1', 'IRF7', 'XAF1', 'LY6E', 'IFI27'},
    'heat_shock_stress': {'HSPA1A', 'HSPA1B', 'HSPA6', 'HSPA8', 'HSPH1', 'HSP90AA1', 'HSP90AB1',
                          'DNAJB1', 'DNAJA1', 'HSPD1', 'HSPE1', 'BAG3', 'FOS', 'FOSB', 'JUN',
                          'JUNB', 'ATF3', 'EGR1', 'IER2', 'DUSP1'},
    'cell_cycle': {'MKI67', 'TOP2A', 'CDK1', 'UBE2C', 'BIRC5', 'NUSAP1', 'TYMS', 'RRM2', 'ASPM',
                   'DLGAP5', 'CENPF', 'TUBA1B', 'STMN1', 'HMGB2'},
    'ribosomal_translation': None,   # RPL*/RPS*/EEF1*/EIF*
    'mitochondrial': None,           # MT-*
}


def family_hits(genes):
    out = {}
    for fam, members in FAMILIES.items():
        if fam == 'ribosomal_translation':
            hit = [g for g in genes if g.startswith(('RPL', 'RPS', 'EEF1', 'EIF'))]
        elif fam == 'mitochondrial':
            hit = [g for g in genes if g.startswith('MT-')]
        else:
            hit = [g for g in genes if g in members]
        out[fam] = len(hit)
    return out


def load_setting(name, run):
    t = pd.read_csv(os.path.join(run, 'outputs', 'triage.tsv'), sep='\t')
    summ = json.load(open(os.path.join(run, 'outputs', 'summary_v1.json')))
    t['setting'] = name
    t['genes'] = t['top_genes'].str.split(',')
    fams = pd.DataFrame([family_hits(g) for g in t['genes']])
    return pd.concat([t.reset_index(drop=True), fams.add_prefix('n_top_')], axis=1), summ


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out'); ap.add_argument('runs', nargs='+')
    ap.add_argument('--min-jaccard', type=float, default=0.3)
    a = ap.parse_args()
    runs = [r.split('=', 1) for r in a.runs]
    tabs, summ = {}, {}
    for name, run in runs:
        tabs[name], summ[name] = load_setting(name, run)
    allmp = pd.concat(tabs.values(), ignore_index=True)
    os.makedirs(a.out, exist_ok=True)
    allmp.drop(columns='genes').to_csv(os.path.join(a.out, 'all_settings_meta_programs.tsv'), sep='\t', index=False)
    ref = runs[0][0]
    rows = []
    for _, r in tabs[ref].iterrows():
        row = {'ref_mp': r['mp'], 'ref_tier': r['tier'], 'ref_top_genes': r['top_genes']}
        stable = 0
        for name, _ in runs[1:]:
            best, bj = None, 0.0
            for _, s in tabs[name].iterrows():
                j = len(set(r['genes']) & set(s['genes'])) / len(set(r['genes']) | set(s['genes']))
                if j > bj:
                    best, bj = s, j
            if best is not None and bj >= a.min_jaccard:
                row[f'{name}_mp'] = best['mp']; row[f'{name}_tier'] = best['tier']; row[f'{name}_jaccard'] = round(bj, 3)
                stable += 1
            else:
                row[f'{name}_mp'] = ''; row[f'{name}_tier'] = 'no_match'; row[f'{name}_jaccard'] = round(bj, 3)
        row['n_other_settings_matched'] = stable
        rows.append(row)
    m = pd.DataFrame(rows)
    m.to_csv(os.path.join(a.out, 'cross_setting_matches.tsv'), sep='\t', index=False)
    rep = {'reference_setting': ref, 'settings': [n for n, _ in runs], 'min_jaccard': a.min_jaccard,
           'n_meta_programs_per_setting': {n: int(len(t)) for n, t in tabs.items()},
           'tiers_per_setting': {n: summ[n].get('tiers') for n in summ},
           'review_aid_only_no_naming_no_freezing_no_CM': True,
           'marker_families_are_descriptive_flags': sorted(FAMILIES)}
    json.dump(rep, open(os.path.join(a.out, 'review_summary.json'), 'w'), indent=2)
    print(m[['ref_mp', 'ref_tier', 'n_other_settings_matched']].to_string(index=False))


if __name__ == '__main__':
    main()
