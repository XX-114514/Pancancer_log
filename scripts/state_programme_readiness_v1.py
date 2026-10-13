"""Reference coverage/overlap/complexity diagnostics on accepted major-unit caches.
Not Spectra training, UCell, a cell-state annotation or a communication test.
"""
from pathlib import Path
from datetime import datetime, timezone
import os, json, hashlib, math
import numpy as np
import pandas as pd
from scipy.sparse import load_npz
from scipy.stats import spearmanr


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8388608), b''):
            h.update(block)
    return h.hexdigest()


def main():
    started = datetime.now(timezone.utc).isoformat()
    R = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    O = R / 'cell_state_programme_readiness_20261013_v1'
    N = R / 'immune_same_name_family_source_review_v1/major_pair_native_assembly_v1'
    plan = json.loads(Path(os.environ['PANCANCER_STATE_PLAN']).read_text())
    reference = Path(os.environ['PANCANCER_STATE_REFERENCE'])
    assert not O.exists() and sha(__file__) == plan['producer_sha256']
    assert sha(reference) == plan['reference_sha256']
    inv = json.loads((N / 'inventory.json').read_text())
    audit = json.loads((N / 'independent_audit.json').read_text())
    assert sha(N / 'inventory.json') == plan['native_inventory_sha256'] == audit['producer_inventory_sha256']
    assert audit['status'] == 'PASS_ALL276_MAJOR_DIRECTIONS_1276224_CASES_AND_NATIVE_FINE_UNIT_AGGREGATION'
    resource = Path(os.environ['PANCANCER_FROZEN_LR_RESOURCE'])
    assert sha(resource) == plan['LR_sha256']
    lr = pd.read_csv(resource, sep='\t')
    assert len(lr) == 4624
    lrgenes = {g for c in ['ligand', 'receptor'] for v in lr[c] for g in v.split('_')}
    raw = json.loads(reference.read_text())
    programs = [(scope, name, sorted(set(gs))) for scope, sets in raw.items() for name, gs in sets.items()]
    assert len(programs) == plan['reference_programs'] == 181
    assert all(gs for _, _, gs in programs)
    O.mkdir()
    rows, catalogue, overlap, fixtures = [], [], [], []
    for scope, name, gs in programs:
        catalogue.append(dict(reference_scope=scope, program=name, genes=len(gs),
                              global_LR_overlap=len(set(gs) & lrgenes),
                              gene_sha256=hashlib.sha256(';'.join(gs).encode()).hexdigest(),
                              is_formal_cell_state=False))
    for i, (scope, name, gs) in enumerate(programs):
        for scope2, name2, gs2 in programs[i+1:]:
            inter = len(set(gs) & set(gs2))
            overlap.append(dict(program1=scope+'::'+name, program2=scope2+'::'+name2,
                                shared_genes=inter, jaccard=inter/len(set(gs) | set(gs2)),
                                overlap_min_fraction=inter/min(len(gs), len(gs2))))
    pd.DataFrame(catalogue).to_csv(O / 'all181_reference_program_catalogue.tsv', sep='\t', index=False)
    pd.DataFrame(overlap).to_csv(O / 'all16290_reference_overlap.tsv.gz', sep='\t', index=False)
    total_units = 0
    for sid in plan['sources']:
        paths = {k: N/'private_source_major_arrays'/(sid+suffix) for k, suffix in
                 [('units','_units.tsv'), ('axes','_axes.npz'), ('mean','_meanlog.npz'), ('pct','_pct_positive.npz')]}
        for p in paths.values():
            assert sha(p) == inv['outputs'][str(p.relative_to(N))]
        unit = pd.read_csv(paths['units'], sep='\t')
        axes = np.load(paths['axes'], allow_pickle=False)
        genes = list(map(str, axes['genes'])); gi = {g: j for j, g in enumerate(genes)}
        assert len(gi) == len(genes) and np.array_equal(unit.unit_row, np.arange(len(unit)))
        mean = load_npz(paths['mean']).tocsr(); pct = load_npz(paths['pct']).tocsr()
        depth = axes['mean_detected_nonLR_genes']
        assert mean.shape == pct.shape == (len(unit), len(genes)) and len(depth) == len(unit)
        assert (unit.n_cells >= 20).all()
        total_units += len(unit)
        for major, group in unit.groupby('majorCluster', sort=True):
            ii = group.unit_row.to_numpy()
            for scope, name, gs in programs:
                matched = [g for g in gs if g in gi]; cols = [gi[g] for g in matched]
                nonlr = [g for g in matched if g not in lrgenes]
                score = np.asarray(mean[ii][:,cols].mean(axis=1)).ravel() if cols else None
                detected = (pct[ii][:,cols].toarray() >= .01).sum(axis=1) if cols else None
                required = max(1, min(5, math.ceil(.3*len(matched))))
                fraction = float((detected >= required).mean()) if cols else None
                rho = float(spearmanr(score, depth[ii]).statistic) if cols and len(ii)>=10 and np.ptp(score)>1e-12 and np.ptp(depth[ii])>1e-12 else None
                if rho is not None and not np.isfinite(rho): rho = None
                row = dict(gse_id=sid, majorCluster=major, reference_scope=scope, program=name,
                           libraries=len(ii), cells=int(group.n_cells.sum()), reference_genes=len(gs),
                           matched_symbol_genes=len(matched), symbol_match_fraction=len(matched)/len(gs),
                           matched_genes_after_global_LR_exclusion=len(nonlr),
                           fraction_genes_lost_by_global_LR_exclusion=(len(matched)-len(nonlr))/len(matched) if matched else None,
                           detection_required_genes=required, fraction_libraries_meeting_detection=fraction,
                           descriptive_coverage_gate=bool(len(matched)>=5 and len(matched)/len(gs)>=.7 and fraction is not None and fraction>=.8),
                           mean_log_program_score=float(score.mean()) if cols else None,
                           score_sd=float(score.std()) if cols else None,
                           descriptive_complexity_spearman=rho,
                           applicability='GLOBAL_REFERENCE' if scope=='global' else 'FINE_IDENTITY_APPLICABILITY_REQUIRES_REVIEW',
                           state_assignment=False, patient_inference=False)
                rows.append(row)
                # Independent direct scalar loop for one fixed program in every source-major.
                if (scope, name) == ('MDC','Mac_IFNG_response'):
                    expected = [sum(float(mean[r, gi[g]]) for g in matched)/len(matched) for r in ii] if matched else None
                    assert expected is None or np.allclose(expected, score, rtol=1e-12, atol=1e-12)
                    if matched:
                        expected_detect = [sum(float(pct[r,gi[g]])>=.01 for g in matched) for r in ii]
                        assert np.array_equal(expected_detect, detected)
                    fixtures.append(dict(gse_id=sid,majorCluster=major,libraries=len(ii),
                                         program=name,status='PASS_DIRECT_SCALAR_SCORE_AND_DETECTION'))
        print(json.dumps(dict(gse_id=sid,units=len(unit),major_groups=unit.majorCluster.nunique())), flush=True)
    assert total_units == 2028
    table = pd.DataFrame(rows)
    assert len(table) == plan['expected_source_major_groups']*181
    assert not table[['gse_id','majorCluster','reference_scope','program']].duplicated().any()
    table.to_csv(O/'all_source_major_program_readiness.tsv.gz',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(fixtures).to_csv(O/'direct_scalar_diagnostic_checks.tsv',sep='\t',index=False)
    summary = dict(schema='cell_state_programme_readiness_v1',status='COMPLETE_REFERENCE_COVERAGE_DIAGNOSTICS_ONLY',
                   started_UTC=started,finished_UTC=datetime.now(timezone.utc).isoformat(),
                   producer_sha256=sha(__file__),reference_sha256=sha(reference),
                   native_inventory_sha256=sha(N/'inventory.json'),programs=181,sources=len(plan['sources']),
                   original_major_units=total_units,source_major_groups=plan['expected_source_major_groups'],
                   rows=len(table),descriptive_gate_rows=int(table.descriptive_coverage_gate.sum()),
                   LR_exclusion_loses_at_least_half_rows=int(table.fraction_genes_lost_by_global_LR_exclusion.ge(.5).sum()),
                   complexity_rho_abs_ge_point5_rows=int(table.descriptive_complexity_spearman.abs().ge(.5).sum()),
                   complexity_rho_estimable_rows=int(table.descriptive_complexity_spearman.notna().sum()),
                   scalar_checks=len(fixtures),matrix_factorization_executed=False,
                   cell_level_UCell_executed=False,formal_state_labels=False,new_LR_tests=False,
                   scientific_release=False,axis_validation_complete=False,
                   limits=['Only existing selected major/library aggregates; not all cells or subtype-conditional state inference',
                           'Reference matching and arbitrary diagnostic gates do not prove state presence or absence',
                           'Unmatched symbols may reflect aliases or upstream feature filtering, never silently filled zero',
                           'Scores are descriptive mean log1pCP10K, not Spectra/UCell or activation thresholds',
                           'Complexity correlation may be biological or technical; no p-values or patient independence claimed',
                           'Reference fine immune labels do not promote broad T/NK into CD4/CD8/Treg identities'],
                   outputs={f.name:sha(f) for f in O.iterdir() if f.is_file()})
    (O/'inventory.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps(summary),flush=True)


if __name__=='__main__':main()
