"""Full source-result conservation plus independent standard metric fixtures."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os

import numpy as np
import pandas as pd
from sklearn.metrics import balanced_accuracy_score, f1_score


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    source = Path(os.environ['PANCANCER_RAM_OUTPUT'])
    out = Path(os.environ['PANCANCER_RAM_AUDIT_OUTPUT'])
    producer = Path(os.environ['PANCANCER_REFERENCE_PRODUCER'])
    assert out.parent == Path('/dev/shm') and not out.exists()
    terminal = json.loads((source / 'private_terminal.json').read_text())
    assert sha(producer) == terminal['source_SHA256']
    for name, v in terminal['files'].items(): assert sha(source / name) == v['SHA256']
    contexts = json.loads((source / 'private_all8_reference_contexts.json').read_text())
    scores = pd.DataFrame(json.loads((source / 'private_all_source_class_prediction_descriptions.json').read_text()))
    counts = pd.DataFrame(json.loads((source / 'private_all_source_candidate_confusion_counts.json').read_text()))
    numeric = pd.DataFrame(json.loads((source / 'private_all_reference_numeric_proofs.json').read_text()))
    assert len(contexts) == 8 and len({(x['task'], x['heldout_source']) for x in contexts}) == 8
    assert not scores.duplicated(['task','heldout_source','model','candidate']).any()
    assert not counts.duplicated(['task','heldout_source','model','observed_candidate','predicted_candidate']).any()
    reviewed = 0
    for x in contexts:
        assert sum(v['train_cells'] for v in x['class_support']) == x['training_cells']
        assert sum(v['heldout_cells'] for v in x['class_support']) == x['heldout_cells']
        if x['status'] != 'COMPLETE_FIXED_SOURCE_HELDOUT_DIAGNOSTIC': continue
        complete = all(v['heldout_cells'] >= 20 and v['heldout_material_units'] >= 2 for v in x['class_support'])
        assert complete == x['all_heldout_classes_support_full_macro']
        assert all(v['train_cells'] >= 20 and v['train_sources'] >= 2 and v['train_material_units'] >= 3 for v in x['class_support'])
        for model in ['QC_only','RNA_PCs_plus_QC']:
            g = counts[counts.task.eq(x['task']) & counts.heldout_source.eq(x['heldout_source']) & counts.model.eq(model)]
            s = scores[scores.task.eq(x['task']) & scores.heldout_source.eq(x['heldout_source']) & scores.model.eq(model)].set_index('candidate')
            assert int(g.cells.sum()) == x['heldout_cells']
            for support in x['class_support']:
                name = support['candidate']; part = g[g.observed_candidate.eq(name)]
                assert int(part.cells.sum()) == support['heldout_cells'] == int(s.loc[name,'heldout_cells'])
                assert int(s.loc[name,'heldout_material_units']) == support['heldout_material_units']
                hit = int(part.loc[part.predicted_candidate.eq(name), 'cells'].sum())
                if support['heldout_cells']: assert abs(hit / support['heldout_cells'] - float(s.loc[name,'recall'])) < 1e-12
            metric = x['metrics'][model]
            assert metric['balanced_material_class_logloss'] >= 0 and metric['sampled_cell_logloss'] >= 0
            if complete: assert 0 <= metric['macroF1'] <= 1 and 0 <= metric['balanced_accuracy'] <= 1
            else: assert metric['macroF1'] is None and metric['balanced_accuracy'] is None
            n = numeric[numeric.task.eq(x['task']) & numeric.heldout_source.eq(x['heldout_source']) & numeric.model.eq(model)]
            assert len(n) == 1 and n.iloc[0].training_scaled_KKT_error < 1e-5 and n.iloc[0].heldout_source_excluded_from_training
            assert n.iloc[0].feature_indices_SHA256 == x['training_feature_indices_SHA256']
            reviewed += 1
        difference = x['metrics']['QC_only']['balanced_material_class_logloss'] - x['metrics']['RNA_PCs_plus_QC']['balanced_material_class_logloss']
        assert abs(difference - x['paired_class_material_balanced_logloss_improvement']) < 1e-12
        boot = x['conditional_library_bootstrap']
        assert boot['declared_draws'] == boot['valid_draws'] + boot['failed_missing_class_draws'] == 2000
        if 'CI_low' in boot: assert boot['valid_draws'] >= 1800 and np.isfinite([boot['CI_low'],boot['CI_high']]).all() and boot['CI_low'] <= boot['CI_high']
    for task in sorted({x['task'] for x in contexts}):
        rows = [x for x in contexts if x['task'] == task]
        total = sum(x['heldout_cells'] for x in rows)
        assert all(x['training_cells'] + x['heldout_cells'] == total for x in rows)
    spec = importlib.util.spec_from_file_location('accepted_reference_producer', producer)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    rng = np.random.default_rng(2026102007)
    for case, k in enumerate([3, 4, 7]):
        classes = ['class_' + str(i) for i in range(k)]
        truth = np.tile(np.array(classes), 100)
        p = rng.uniform(.01, 1, (len(truth), k))
        if case == 1: p[:,0] += 2
        p /= p.sum(axis=1, keepdims=True)
        weights = rng.lognormal(size=len(truth))
        metric, _, _ = module.metrics(truth, p, weights, classes)
        pred = np.array(classes)[p.argmax(axis=1)]
        np.testing.assert_allclose(metric['macroF1'], f1_score(truth,pred,average='macro',labels=classes,sample_weight=weights,zero_division=0),atol=1e-12)
        np.testing.assert_allclose(metric['balanced_accuracy'],balanced_accuracy_score(truth,pred,sample_weight=weights),atol=1e-12)
    out.mkdir(mode=0o700)
    result = dict(status='ACCEPTED_ALL_SOURCE_SUPPORT_COUNTS_RECALL_LOSS_DIFFERENCES_NUMERIC_PROOFS_AND_STANDARD_METRIC_FIXTURES', completed_UTC=datetime.now(timezone.utc).isoformat(), source_SHA256=sha(__file__), parent_terminal_SHA256=sha(source/'private_terminal.json'), registered_contexts=8, model_aggregate_checks=reviewed, source_confusion_records=len(counts), source_class_metric_records=len(scores), independent_weighted_macro_metric_fixture_families=3, real_PCA_or_full_RNA_refit=False, per_cell_rows_or_gene_weights_exported=False)
    (out/'private_source_aggregate_acceptance.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__': main()
