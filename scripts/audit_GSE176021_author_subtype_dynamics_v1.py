"""Independently rebuild all subtype fractions, contributions and covariance terms."""
import json
import os
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def close(a, b):
    assert np.allclose(a, b, atol=1e-10, rtol=0, equal_nan=True)


def covariance(a, b):
    return float(np.dot(a-a.mean(), b-b.mean())/(len(a)-1))


def spearman(a, b):
    if min(np.ptp(a), np.ptp(b)) == 0:
        return np.nan
    x, y = rankdata(a), rankdata(b)
    return float(np.dot(x-x.mean(), y-y.mean())/np.sqrt(np.sum((x-x.mean())**2)*np.sum((y-y.mean())**2)))


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root, 2*1024**2)
    src = root/'GSE176021_author_subtype_dynamics_decomposition_20261019_v1'
    assert sha(src/'private_terminal.json') == '3861070378c76028cc7e3aad19adedd36b6d64f20ba4f3e4420896a6113df40a'
    for folder, digest in [
        ('GSE176021_author_lineage_native_RNA_audit_20261019_v1', '6950e1c15fda41b93962cfd763348898802d8875f1eef6719ec395b25392730d'),
        ('GSE176021_paired_tissue_input_audit_20261019_v1', '1998d626c57e0165ab4b66573170d16a904c600e442823b86dfb37dd39ea5f89'),
        ('GSE176021_exact_original_author_annotation_audit_20261019_v2', 'd8d21a1c27465d9b8248857c0028e5106509e43426cbdfef9f44ccb1d71c4169')]:
        assert sha(root/folder/'private_terminal.json') == digest
    term = json.loads((src/'private_terminal.json').read_text())
    for name, info in term['files'].items():
        assert sha(src/name) == info['sha256']
    out = root/'GSE176021_author_subtype_dynamics_audit_20261019_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    scope = root/'GSE176021_author_lineage_RNA_scope_20261019_v2'
    inp = root/'GSE176021_paired_tissue_inputs_20261019_v1'
    cases = np.load(inp/'private_AP01_inputs.npz')['cases'].tolist()
    pairs = pd.read_csv(scope/'private_all_scope_material_pair_positions.tsv', sep='\t')
    pairs = pairs[pairs.scope.eq('all_author_lineage') & pairs.case_block.isin(cases) & pairs.both_views_eligible & pairs.tissue.isin(['tumor', 'adjacent normal lung'])]
    assert len(pairs) == 60 and not pairs.sampleID.duplicated().any()
    material = {r.sampleID:(r.case_block, r.tissue) for r in pairs.itertuples(index=False)}
    nmat = Counter(material.values())
    pools = pd.read_csv(scope/'private_all_eligible_pools.tsv', sep='\t')
    pools = pools[pools.scope.eq('all_author_lineage') & pools.sampleID.isin(material)]
    assert len(pools) == 120
    selected = np.load(root/'GSE176021_author_lineage_native_RNA_20261019_v1/private_selected_cell_quality_scores.npz')
    membership = np.load(scope/'private_all_scope_cell_memberships.npz')
    side = root/'GSE176021_exact_original_author_annotation_sidecar_20261019_v2'
    bridge = np.load(side/'private_all_atlas_cells_exact_author_bridge.npz')
    labels = json.loads((side/'private_author_label_lookup_and_coverage.json').read_text())
    scales = pd.read_csv(inp/'private_original_discovery_case_scales.tsv', sep='\t').set_index(['view', 'program']).discovery_case_sd
    totals_in = pd.read_csv(inp/'private_all_scope_case_tissue_values.tsv', sep='\t')
    totals_in = totals_in[totals_in.scope.eq('all_author_lineage')].set_index(['case_block', 'tissue', 'view', 'program', 'mode']).value
    fraction = defaultdict(float)
    score = defaultdict(lambda:np.zeros(4))
    full_fraction = defaultdict(float)
    selected_n = Counter()
    seen = []
    for p in pools.itertuples(index=False):
        cc, tissue = material[p.sampleID]
        sel = np.flatnonzero(selected['pool_index'] == p.pool_index)
        whole = np.flatnonzero(membership['pool_index'] == p.pool_index)
        assert len(sel) == p.selected_cells and len(whole) == p.cells
        selected_n[(cc, tissue, p.view)] += len(sel)
        seen.extend(selected['global_row'][sel].tolist())
        for is_selected, rows in [(True, selected['global_row'][sel]), (False, membership['global_row'][whole])]:
            pos = np.searchsorted(bridge['global_atlas_row'], rows)
            assert np.array_equal(bridge['global_atlas_row'][pos], rows)
            c3 = bridge['CD3_label_code'][pos]
            c8 = bridge['CD8_label_code'][pos]
            assert np.all(c3 >= 0)
            refs = {'CD3':np.asarray(labels['CD3']['labels'], object)[c3]}
            refs['CD8_refined'] = refs['CD3'] if p.view == 'R_FOXP3_capture' else np.asarray([labels['CD8']['labels'][x] if x >= 0 else 'UNAVAILABLE_IN_REFINED_CD8_REFERENCE' for x in c8], object)
            weight = 1/(nmat[(cc, tissue)]*len(rows))
            for ref, names in refs.items():
                for name in sorted(set(names)):
                    mask = names == name
                    key = (cc, tissue, p.view, ref, name)
                    if is_selected:
                        fraction[key] += int(mask.sum())*weight
                        score[key] += selected['scores'][sel[mask]].sum(axis=0)*weight
                    else:
                        full_fraction[key] += int(mask.sum())*weight
    assert len(seen) == 11939 and len(set(seen)) == len(seen)
    expected_parts = {}
    expected_totals = {}
    for ref in ['CD3', 'CD8_refined']:
        for pr in ['MCP1', 'MCP2']:
            for st, variant in enumerate(['full', 'fixed_stress_complement']):
                col = 2*(int(pr[-1])-1)+st
                for cc in cases:
                    for view in ['R_FOXP3_capture', 'E_effector_capture']:
                        def values(tissue):
                            return {k[-1]:(fraction[k], score[k][col]/scales.loc[(view, pr)]) for k in fraction if k[:4] == (cc, tissue, view, ref)}
                        t, n = values('tumor'), values('adjacent normal lung')
                        mt, mn = sum(v[1] for v in t.values()), sum(v[1] for v in n.values())
                        center = (mt+mn)/2
                        mode = 'primary' if st == 0 else 'fixed15_stress_complement'
                        close(mt, totals_in.loc[(cc, 'tumor', view, pr, mode)]/scales.loc[(view, pr)])
                        close(mn, totals_in.loc[(cc, 'adjacent normal lung', view, pr, mode)]/scales.loc[(view, pr)])
                        contributions = []
                        shared, one, misleading = 0, 0, 0
                        for label in sorted(set(t)|set(n)):
                            pt, wt = t.get(label, (0., 0.))
                            pn, wn = n.get(label, (0., 0.))
                            ut, un = wt/pt if pt else np.nan, wn/pn if pn else np.nan
                            if pt and pn:
                                c, w, u = (pt-pn)*((ut+un)/2-center), (pt+pn)*(ut-un)/2, 0.
                                status = 'SHARED_SELECTED_SUBTYPE'
                                shared += 1
                            else:
                                c, w, u = np.nan, np.nan, wt-wn-(pt-pn)*center
                                status = 'ONLY_TUMOUR_SELECTED' if pt else 'ONLY_NORMAL_SELECTED'
                                one += 1
                            ft = full_fraction.get((cc, 'tumor', view, ref, label), 0.)
                            fn = full_fraction.get((cc, 'adjacent normal lung', view, ref, label), 0.)
                            misleading += int(not (pt and pn) and ft > 0 and fn > 0)
                            key = (ref, pr, variant, cc, view, label)
                            expected_parts[key] = ([pt, pn, ut, un, center, c, w, u, ft, fn], status)
                            contributions.append([c, w, u])
                        sums = np.nansum(contributions, axis=0)
                        close(sums.sum(), mt-mn)
                        expected_totals[(ref, pr, variant, cc, view)] = [*sums, mt-mn, sums.sum()-mt+mn, shared, one, misleading]
    parts = pd.read_csv(src/'private_all_subtype_case_components.tsv', sep='\t')
    pkeys = ['reference', 'program', 'variant', 'case_block', 'view', 'author_subtype']
    pfields = ['p_tumour_selected', 'p_normal_selected', 'mean_tumour', 'mean_normal', 'common_center', 'composition_shared', 'within_shared', 'unshared_observed', 'p_tumour_all_cells', 'p_normal_all_cells']
    assert len(parts) == len(expected_parts) == 928 and not parts.duplicated(pkeys).any()
    for _, row in parts.iterrows():
        vals, status = expected_parts[tuple(row[pkeys])]
        close(row[pfields].to_numpy(float), vals)
        assert row.status == status
    totals = pd.read_csv(src/'private_all160_case_view_decomposition_totals.tsv', sep='\t')
    tkeys = pkeys[:-1]
    tfields = ['composition_shared', 'within_shared', 'unshared_observed', 'total_delta', 'closure_error', 'shared_subtypes', 'one_sided_selected_subtypes', 'one_sided_but_seen_both_all_cell_subtypes']
    assert len(totals) == len(expected_totals) == 160 and not totals.duplicated(tkeys).any()
    for _, row in totals.iterrows():
        close(row[tfields].to_numpy(float), expected_totals[tuple(row[tkeys])])
    sampling = pd.read_csv(src/'private_all80_selected_vs_all_cell_reference_coverage.tsv', sep='\t')
    assert len(sampling) == 80 and not sampling.duplicated(['reference', 'case_block', 'tissue', 'view']).any()
    for r in sampling.itertuples(index=False):
        prefix = (r.case_block, r.tissue, r.view, r.reference)
        labels_here = {k[-1] for k in full_fraction if k[:4] == prefix} | {k[-1] for k in fraction if k[:4] == prefix}
        f = np.array([full_fraction.get((*prefix, l), 0.) for l in labels_here])
        s = np.array([fraction.get((*prefix, l), 0.) for l in labels_here])
        close([f.sum(), s.sum()], [1., 1.])
        close([r.fraction_TV, r.all_cell_fraction_not_seen_in_selected], [.5*np.abs(f-s).sum(), f[s == 0].sum()])
        assert r.selected_cells == selected_n[prefix[:3]]
    descriptors = pd.read_csv(src/'private_all8_source_covariance_decompositions_no_P.tsv', sep='\t')
    assert len(descriptors) == 8 and not descriptors.duplicated(['reference', 'program', 'variant']).any()
    for _, row in descriptors.iterrows():
        base = (row.reference, row.program, row.variant)
        rd = np.array([expected_totals[(*base, cc, 'R_FOXP3_capture')][3] for cc in cases])
        ed = np.array([expected_totals[(*base, cc, 'E_effector_capture')][3] for cc in cases])
        close([row.Pearson_total_covariance, row.Spearman_total], [covariance(rd, ed), spearman(rd, ed)])
        for j, name in enumerate(['composition_shared', 'within_shared', 'unshared_observed']):
            a = np.array([expected_totals[(*base, cc, 'E_effector_capture')][j] for cc in cases])
            loo = [covariance(np.delete(rd, z), np.delete(a, z)) for z in range(10)]
            close([row[name+'_covariance_with_R'], row[name+'_Spearman_with_R'], row[name+'_mean_E_delta'], row[name+'_LOO_cov_min'], row[name+'_LOO_cov_max']], [covariance(rd, a), spearman(rd, a), a.mean(), min(loo), max(loo)])
        assert row.case_units == 10
    receipt = {'status':'ACCEPTED_ALL160_TOTALS_928_SUBTYPE_CONTRIBUTIONS_80_COVERAGE_AND8_COVARIANCE_DESCRIPTORS', 'completed_UTC':datetime.now(timezone.utc).isoformat(), 'source_sha256':sha(__file__), 'producer_terminal_sha256':sha(src/'private_terminal.json'), 'resource_guard':guard, 'all_positions_recomputed':[160,928,80,8], 'no_new_P_or_gene_fit':True, 'private_individual_rows_not_exported':True}
    (out/'private_terminal.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({'status':receipt['status'], 'terminal_sha256':sha(out/'private_terminal.json')}))


if __name__ == '__main__':
    main()
