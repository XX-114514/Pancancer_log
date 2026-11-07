"""Independent saved-reference readback; training completion is not identity proof."""
import json, hashlib
from datetime import datetime, timezone
import numpy as np
import pandas as pd

r = PROJECT_ROOT / 'runs/20261107_v4_Visium_author_validation_v2'
p = r / 'reference'
receipt = json.loads((p / 'complete.json').read_text())
terminal = json.loads((r / 'launch/fit_reference_retry2_terminal.json').read_text())
assert terminal['exit_code'] == 0
assert receipt['epochs_recorded'] == receipt['epochs_requested'] == 250
assert receipt['posterior_samples'] == 1000
sig = pd.read_csv(p / 'signatures_PRIVATE.tsv', sep='\t', index_col=0)
expected = pd.read_csv(p / 'node_reference_coverage.tsv', sep='\t')
assert sig.shape == (receipt['n_genes_after_recommended_filter'], 50)
assert set(sig.columns) == set(expected.node)
assert sig.index.is_unique and sig.columns.is_unique
assert np.isfinite(sig.to_numpy()).all() and (sig.to_numpy() >= 0).all()
assert (sig.sum(axis=0) > 0).all()
digest = hashlib.sha256((p / 'signatures_PRIVATE.tsv').read_bytes()).hexdigest()
assert digest == receipt['signatures_sha256']
loss = np.asarray(np.load(p / 'ELBO_history_PRIVATE.npy', allow_pickle=True), dtype=float)
assert len(loss) == 250 and np.isfinite(loss).all()
assert float(loss[-1]) == receipt['ELBO_final']
assert set(pd.read_csv(p / 'model_genes.tsv', sep='\t').gene) == set(sig.index)
separation = pd.read_csv(p / 'all_panel_label_heldout_separability.tsv', sep='\t')
for panel, group in separation.groupby('gene_panel_key'):
    assert len(group) == 50 and set(group.node) == set(expected.node)
result = {
    'status': 'PASS_independent_fitted_reference_50_labels_250_epoch_finite_signatures_saved_loss_and_panel_scope',
    'UTC': datetime.now(timezone.utc).isoformat(),
    'n_reference_labels': 50, 'n_model_genes': len(sig),
    'epochs_readback': len(loss), 'posterior_samples_requested_by_executed_code': 1000,
    'ELBO_initial': float(loss[0]), 'ELBO_final': float(loss[-1]),
    'ELBO_last20_mean_relative_to_previous20': float(np.mean(loss[-20:]) / np.mean(loss[-40:-20]) - 1),
    'n_native_gene_panels': separation.gene_panel_key.nunique(),
    'signatures_SHA256': digest,
    'posterior_mean_signature_finite_nonnegative': True,
    'formal_convergence_or_independent_subtype_accuracy_claim': False,
}
(r / 'audit/fitted_reference_audit.json').open('x').write(json.dumps(result, indent=2))
print(json.dumps(result))
