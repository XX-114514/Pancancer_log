"""Replay every source-level natural-library metric; no individual RNA refit."""
from datetime import datetime, timezone
from pathlib import Path
from itertools import product
import hashlib
import json
import math
import os


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def close(a, b, atol=2e-8):
    assert a is not None and b is not None and math.isfinite(a) and math.isfinite(b)
    assert math.isclose(a, b, rel_tol=2e-6, abs_tol=atol)


def centered(n, first, second):
    v = second-first*first/n
    assert v >= -1e-10
    return max(0., v)


def replay(r):
    n = r['libraries']; assert n > 0
    a, b = r['sum_truth'], r['sum_estimate']
    va = centered(n, a, r['sum_truth2'])
    vb = centered(n, b, r['sum_estimate2'])
    cov = r['sum_product']-a*b/n
    close(r['true_mean'], a/n); close(r['estimate_mean'], b/n)
    close(r['MAE'], r['sum_absolute_error']/n)
    close(r['RMSE'], math.sqrt(r['sum_squared_error']/n))
    close(r['bias'], (b-a)/n)
    assert r['RMSE']+1e-12 >= r['MAE'] >= abs(r['bias'])-1e-12
    assert 0 <= r['true_mean'] <= 1+1e-12 and r['estimate_mean'] >= 0
    # Raw NNLS coefficients can exceed one. No fraction upper bound is imposed.
    if r['coefficient_view'] == 'normalized_supported_reference_sum':
        assert r['estimate_mean'] <= 1+1e-12
    if r['Pearson'] is not None:
        assert va > 0 and vb > 0
        close(r['Pearson'], cov/math.sqrt(va*vb))
        ra, rb = r['sum_rank_truth'], r['sum_rank_estimate']
        vra = centered(n, ra, r['sum_rank_truth2'])
        vrb = centered(n, rb, r['sum_rank_estimate2'])
        close(r['Spearman'], (r['sum_rank_product']-ra*rb/n)/math.sqrt(vra*vrb))
    else:
        assert r['Spearman'] is None and min(va, vb) <= 1.01e-14
    if r['calibration_slope'] is not None:
        close(r['calibration_slope'], cov/va)
    else:
        assert va <= 1.01e-14
    assert 0 <= r['positive_RNA_ge001_mixtures'] <= n
    assert 0 <= r['designed_zero_mixtures'] <= n
    assert r['positive_RNA_ge001_mixtures']+r['designed_zero_mixtures'] <= n
    for num, fields in [('positive_RNA_ge001_mixtures', ['positive_RNA_ge001_MAE']),
                        ('designed_zero_mixtures', ['designed_zero_leakage_gt005', 'designed_zero_max_contribution'])]:
        for k in fields:
            assert (r[k] is None) == (r[num] == 0)
            if r[k] is not None: assert math.isfinite(r[k]) and r[k] >= 0


def main():
    root = Path(os.environ['PANCANCER_PARENT_OUTPUT'])
    out = Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert out.parent == Path('/dev/shm') and not out.exists()
    terminal = json.loads((root/'private_terminal.json').read_text())
    assert terminal['status'] == 'COMPLETE_NATURAL_LIBRARY_REFERENCE_DIAGNOSTIC_SOURCE_REVIEW_PENDING'
    for n, info in terminal['files'].items():
        assert Path(n).name == n and sha(root/n) == info['SHA256'] and (root/n).stat().st_size == info['bytes']
    read = lambda n: json.loads((root/n).read_text())
    plan = read('private_plan_before_native_expression.json')['plan']
    sources = plan['prior_reference_plan']['source_scope']
    domains = plan['source_domains']; architectures = plan['reference_architectures']
    scales = plan['scalings']; views = plan['coefficient_views']
    keys = ['source','domain','architecture','scaling','coefficient_view']
    key = lambda r: tuple(r[k] for k in keys)
    expected = set(product(sources, domains, architectures, scales, views))
    coverage = read('private_all_context_coverage.json')
    assert len(coverage) == len(expected) == 64
    by_context = {key(r): r for r in coverage}; assert set(by_context) == expected
    coarse = plan['prior_reference_plan']['T_classes']+['NK','Other_supported_reference_pool','Boundary_supported_reference_pool']
    basis = {'separate_original_major': plan['supported_reference_classes'], 'merge_supported_nonT_background': coarse}
    targets = {a: b+list(plan['marginal_targets'])+['paired_R_only_E_minimum'] for a,b in basis.items()}
    rows = read('private_all_target_recovery_metrics.json'); assert len(rows) == 1376
    indexed = {(key(r), r['target']): r for r in rows}; assert len(indexed) == len(rows)
    assert set(indexed) == {(k,t) for k in expected for t in targets[k[2]]}
    failures = read('private_all_context_failures.json')
    for k,c in by_context.items():
        n = c['evaluated_libraries']; assert 0 <= c['conservative_blocks'] <= n <= c['expected_libraries']
        missing = sum(r['libraries'] for r in failures if tuple(r[x] for x in keys[:-1]) == k[:-1])
        assert missing+n == c['expected_libraries']
        for t in targets[k[2]]:
            r = indexed[k,t]; assert r['libraries'] == n
            if n: replay(r)
            else:
                assert r['status'] == c['status'] == 'UNESTIMABLE_ALL_LIBRARIES'
                assert all(r[x] is None for x in ['MAE','RMSE','bias','Pearson','Spearman'])
    paired = read('private_paired_background_comparisons.json')
    seen = set()
    for r in paired:
        k = tuple(r[x] for x in ['source','domain','scaling','coefficient_view','target'])
        assert k not in seen; seen.add(k)
        a = (k[0],k[1],'separate_original_major',k[2],k[3])
        b = (k[0],k[1],'merge_supported_nonT_background',k[2],k[3])
        assert r['target'] in set(targets[a[2]]) & set(targets[b[2]])
        assert 0 < r['matched_libraries'] <= min(by_context[a]['evaluated_libraries'],by_context[b]['evaluated_libraries'])
        close(r['fine_minus_merged_MAE'], r['fine_MAE']-r['merged_MAE'])
        assert r['inferential_P_or_FDR'] is False
        if r['matched_libraries'] == by_context[a]['evaluated_libraries'] == by_context[b]['evaluated_libraries']:
            close(r['fine_MAE'], indexed[a,r['target']]['MAE']); close(r['merged_MAE'], indexed[b,r['target']]['MAE'])
    source_rows = read('private_all_source_uncovered_RNA.json')
    assert len(source_rows) == 4 and {r['source'] for r in source_rows} == set(sources)
    meta = read('private_metadata_and_selection_contract.json')
    assert sum(r['material_cells'] for r in source_rows) == meta['material_cells']
    assert sum(r['material_libraries'] for r in source_rows) == meta['material_libraries']
    assert sum(r['conservative_blocks'] for r in source_rows) == meta['conservative_claimed_identity_components']
    for r in source_rows:
        assert 0 <= r['unsupported_RNA_fraction_mean'] <= r['unsupported_RNA_fraction_max'] <= 1
        assert 0 <= r['common_gene_RNA_fraction_mean'] <= 1
        assert r['all24_truth_mass_conserved'] and r['patient_independence_not_newly_verified']
        assert all(c['expected_libraries'] == r['material_libraries'] for c in coverage if c['source'] == r['source'])
    proof = read('private_full_numeric_proofs.json'); fixture = read('private_numerical_fixture.json')
    assert proof['model_fits'] == sum(c['evaluated_libraries'] for c in coverage)//2
    assert proof['KKT_max'] < 1e-7 and proof['direct_max_coefficient_error'] < 1e-6
    assert proof['direct_checks'] == sum(c['evaluated_libraries']>0 for c in coverage)//2
    assert all(proof[k] for k in ['all_material_T_cache_counts_depth_detected_exact','full_gene_aggregation_sum_proofs','all24_native_RNA_mass_conserved','all_original_selection_hashes_exact'])
    assert all(fixture[k] for k in ['flat_reference_indices_beyond_int8_are_exact','exact_known_coefficients_recovered','unnormalized_mass_distinct_from_normalized_share','same_source_claimed_identity_linked','cross_source_and_missing_tokens_not_linked'])
    result = dict(status='ACCEPTED_ALL_NATURAL_LIBRARY_SOURCE_AGGREGATES', completed_UTC=datetime.now(timezone.utc).isoformat(), source_SHA256=sha(__file__), parent_terminal_SHA256=sha(root/'private_terminal.json'), all_context_positions=64, all_target_positions=1376, all_source_metric_moments_replayed=True, matched_architecture_algebra_and_available_full_case_means_checked=True, complete_source_denominators_and_failure_accounting=True, no_individual_rows_or_second_raw_RNA_refit=True, raw_coefficients_not_assumed_fractions=True)
    out.mkdir(mode=0o700); (out/'private_independent_acceptance.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__': main()
