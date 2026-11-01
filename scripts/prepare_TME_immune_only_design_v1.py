"""Restrict a frozen qualified-TME design to immune compartments, without refitting labels.

This is an analysis-domain filter, not a new QC or a malignancy classifier.
Parent-normalized frequencies and missingness are handled by the unchanged
accepted panel preparation and independent recount programs downstream.
"""
import argparse, datetime, hashlib, json, shutil, sys
from pathlib import Path
import numpy as np
import pandas as pd


def load(p): return json.loads(Path(p).read_text())
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def bound(b):
    p = Path(b['path'])
    assert sha(p) == b['SHA256']
    return p

def save(p, x):
    with Path(p).open('x') as f: json.dump(x, f, indent=2, allow_nan=False)


def immune_scope_mask(common, parent_codes, parents, allowed, known_malignant=None):
    common = np.asarray(common, bool)
    parent_codes = np.asarray(parent_codes)
    assert parent_codes.shape == common.shape
    assert len(allowed) == len(set(allowed)) and set(allowed) <= {x['name'] for x in parents}
    allowed_codes = [x['code'] for x in parents if x['name'] in allowed]
    mask = common & np.isin(parent_codes, allowed_codes)
    if known_malignant is not None:
        malignant = np.asarray(known_malignant)
        assert malignant.shape == common.shape and malignant.dtype.kind == 'b'
        mask &= ~malignant
    return mask


def main(R):
    cfg = load(R/'config/spec_v1.json')
    sys.path.insert(0, str(Path(cfg['joined'])/'scripts'))
    from extract_full_source_local_marker_evidence_v1 import storage_gate
    storage_gate(R, cfg)
    old = Path(cfg['original_preparation'])
    design = Path(cfg['original_design'])
    assert load(old/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code'] == 0
    assert load(design/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code'] == 0
    spec = load(old/'config/spec_v1.json')
    common_path = old/'outputs/private_common_assessed_cell_mask_v1.npy'
    common = np.load(common_path)
    original = np.load(design/'outputs/private_technical_material_candidate_mask_v1.npy')
    assert common.dtype.kind == 'b' and np.all(~common | original)
    policy = cfg['scope_policy']
    assert hashlib.sha256(json.dumps(policy, sort_keys=True).encode()).hexdigest() == cfg['scope_policy_SHA256']
    assert policy['unit_eligibility_inherited_no_enriched_material_reentry']
    assert set(policy['allowed_parents']) == {'T_NK', 'B_Plasma', 'Myeloid', 'Mast'}
    memberships = [load(bound(b)) for b in spec['schemes']]
    reference = next(x for x in memberships if x['scheme_id'] == 'CT_identity_cluster')
    parents = reference['parents']
    codes = np.load(bound(reference['parent_codes']))
    malignancy = np.load(bound(cfg['known_malignant_cell_mask'])) if cfg.get('known_malignant_cell_mask') else None
    mask = immune_scope_mask(common, codes, parents, policy['allowed_parents'], malignancy)
    for m in memberships:
        assert m['parents'] == parents
        axis = np.load(bound(m['parent_codes']))
        np.testing.assert_array_equal(axis[common], codes[common])
        assert np.isin(axis[mask], [x['code'] for x in parents if x['name'] in policy['allowed_parents']]).all()
    # The fixed downstream interface calls this a technical mask; this derived
    # file explicitly combines inherited technical eligibility and immune scope.
    out = R/'outputs/private_technical_material_candidate_mask_v1.npy'
    np.save(out, mask)
    src = design/'outputs/private_analysis_unit_context_and_coverage_v1.tsv'
    dst = R/'outputs'/src.name
    shutil.copyfile(src, dst)
    assert sha(src) == sha(dst)
    meta = pd.read_csv(dst, sep='\t')
    assert 'cancer_type' in meta
    assert set(meta.cancer_type.dropna().astype(str)) <= set(policy['reviewed_recorded_cancer_labels'])
    deferred = set(policy['deferred_hematologic_cancer_labels'])
    assert deferred <= set(meta.cancer_type.dropna().astype(str)), 'Do not invent unobserved cancer labels'
    for p in (design/'outputs').glob('private_*_unit_mask_v1.npy'):
        a = np.load(p)
        assert a.dtype.kind == 'b' and a.shape == (len(meta),)
        # Unit exclusion uses reviewed recorded cancer labels, never cell proportions.
        a = a & ~meta.cancer_type.astype(str).isin(deferred).to_numpy()
        np.save(R/'outputs'/p.name, a)
    counts = {x['name']: int(np.count_nonzero(mask & (codes == x['code']))) for x in parents}
    assert sum(counts.values()) == int(mask.sum())
    assert all(n == 0 for key, n in counts.items() if key not in policy['allowed_parents'])
    report = dict(schema='TME_immune_only_design_v1',scope_policy=cfg['scope_policy'],
        original_common_mask={'path':str(common_path),'SHA256':sha(common_path)},
        derived_analysis_mask={'path':str(out),'SHA256':sha(out)},
        original_material_metadata_SHA256=sha(src),n_common_cells=int(common.sum()),
        n_immune_assessed_cells=int(mask.sum()),parent_cell_counts=counts,
        excluded_nonimmune_cells=int(np.count_nonzero(common & ~np.isin(codes,[x['code'] for x in parents if x['name'] in policy['allowed_parents']]))),
        explicit_malignancy_mask_available=malignancy is not None,
        no_malignancy_inference_from_immune_identity_or_cell_proportions=True,
        deferred_hematologic_cancer_labels=sorted(deferred),
        no_annotation_QC_or_authoritative_matrix_overwrite=True,
        legacy_technical_mask_filename_is_derived_analysis_domain_not_new_QC=True)
    p = R/'outputs/accepted_immune_only_design_v1.json'
    save(p, report)
    storage_gate(R,cfg)
    save(R/'provenance/terminal_v1.json',dict(status='IMMUNE_ONLY_ASSESSED_DOMAIN_PREPARED',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),acceptance_SHA256=sha(p)))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--output',required=True)
    root=Path(parser.parse_args().output)
    try: main(root)
    except BaseException:
        import traceback
        save(root/'provenance/failure_v1.json',{'traceback':traceback.format_exc()})
        raise
