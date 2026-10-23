"""Restrict candidate identities to an accepted native measurement universe.

Only gene metadata are used. Training expression still determines feature
ranking; no missing feature is created, imputed or converted into a zero.
"""
import hashlib
import json
from pathlib import Path


def restrict_candidates(candidates, heldout, binding, minimum=3000):
    if len(candidates) != len(set(candidates)):
        raise ValueError('Duplicate input feature identities')
    if binding is None:
        return list(candidates)
    path = Path(binding['path']); raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != binding['sha256']:
        raise ValueError('Frozen common measurement contract changed')
    contract = json.loads(raw)
    common = contract['common_nonheldout_feature_candidates']
    if len(common) != len(set(common)) or set(common) & set(heldout):
        raise ValueError('Common feature universe contains duplicated or heldout identities')
    if not contract['no_gene_absence_filled_as_zero'] or contract['counts_read'] or contract['feature_selection_fitted']:
        raise ValueError('Expected measurement metadata, not fitted query-expression features')
    available = set(common)
    selected = [g for g in candidates if g in available]
    if len(selected) < minimum:
        raise ValueError('Insufficient common candidate genes for frozen model')
    return selected
