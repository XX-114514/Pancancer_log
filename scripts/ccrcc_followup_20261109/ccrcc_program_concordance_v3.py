"""Exploratory identity sensitivity, defined after original96 QC review.

Intersect original RNA rules with marker-supported RNA clusters that agree in
all three seeds. No position, ACKR1 or spatial P value determines membership.
Preserve original96 and adjust the two complete views together (192 slots).
"""
from pathlib import Path
from datetime import datetime, timezone
import sys
import json
import hashlib
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import ccrcc_followup_v1 as base
from resumable_public_zip_v1 import atomic_json, exclusive_lock

r = Path(sys.argv[1]).resolve(); o = r / 'program_concordance_v3'


def consensus(path, ref_clusters, name):
    obj = ad.read_h5ad(path, backed='r')
    ix = obj.obs.source_index_PRIVATE.to_numpy(int)
    reference = obj.obs.seed_0.astype(str).to_numpy()
    keep = np.isin(reference, ref_clusters)
    rows = []
    for seed in [1, 2]:
        target = obj.obs['seed_' + str(seed)].astype(str).to_numpy()
        best_labels = []
        for cluster in ref_clusters:
            a = reference == cluster
            candidates = []
            for label in sorted(set(target), key=int):
                b = target == label
                overlap = int((a & b).sum()); union = int((a | b).sum())
                candidates.append((overlap / union, label, overlap, int(b.sum())))
            best = max(candidates, key=lambda z: (z[0], -int(z[1])))
            best_labels.append(best[1])
            rows.append({'program': name, 'reference_seed0_cluster': cluster, 'other_seed': seed,
                         'best_Jaccard_cluster': best[1], 'Jaccard': best[0],
                         'reference_cells': int(a.sum()), 'other_cluster_cells': best[3], 'overlap': best[2]})
        keep &= np.isin(target, best_labels)
    obj.file.close()
    return ix[keep], rows


with exclusive_lock(r / 'program_concordance_v3.lock'):
    o.mkdir(exist_ok=False)
    definitions = [
        ('B', 'identity_refinement_v2/B_ASC/RNA_only_refinement_PRIVATE.h5ad', ['0']),
        ('T', 'identity_refinement_v2/T_NK/RNA_only_refinement_PRIVATE.h5ad', ['0', '1']),
        ('NK', 'identity_refinement_v2/T_NK/RNA_only_refinement_PRIVATE.h5ad', ['2']),
        ('Blood_EC', 'followup_v1/clusters_vascular_stromal/RNA_only_cluster_object_PRIVATE.h5ad', ['0'])]
    atomic_json(o / 'design_before_v3_spatial_effects.json', {
        'UTC': datetime.now(timezone.utc).isoformat(),
        'original96_effects_and_identity_review_already_seen': True,
        'exploratory_sensitivity_not_independent_confirmation': True,
        'membership': 'Original QC/rule identity AND marker-supported RNA program in all3 seeds; max-Jaccard cluster matching',
        'programs': {n: {'parent_object': p, 'seed0_clusters': z} for n, p, z in definitions},
        'program_evidence': {'B': 'MS4A1/BANK1/CD19; CD20 protein',
                             'T': 'CD3/TRAC/CD2; CD3E protein',
                             'NK': 'GNLY/KLRD1/KLRC1/KLRF1/NKG7; lower within-antibody CD3E than T',
                             'Blood_EC': 'PODXL/PECAM1/CD34/ADGRL4/CD93/VWF/CLEC14A; CD31 protein'},
        'spatial_parameters': 'Unchanged complete96; concatenate original96+view96 and use joint192 BH/BY',
        'unmatched_candidates': 'Other in this analysis view; not declared wrong, removed, or doublets',
        'QC_positions_and_lymphatic_mural_rules': 'Preserved',
        'denominator': 'allQC_nonEC means all other QC cells excluding accepted EC and preserved lymphatic candidates',
        'ACKR1_missing_subtype_or_location_used_for_identity': False,
        'V4_or_source_labels_overwritten': False})
    ids, c, xy, xr, rn, xp, pn, scale, total, ng, a = base.load(r)
    pos, count, masks, conflict, depth, qc, original = base.identities(a, total, ng)
    label = np.full(len(ids), 'Other', dtype=object)
    for typ in ['Lymphatic', 'Mural']:
        label[qc & (original == typ)] = typ
    rows, match_rows = [], []
    for typ, path, clusters in definitions:
        selected, matches = consensus(r / path, clusters, typ); match_rows += matches
        eligible = np.zeros(len(ids), dtype=bool); eligible[selected] = True
        selected_mask = eligible & qc & (original == typ)
        label[selected_mask] = typ
        rows.append({'identity': typ, 'original_rule_candidates': int((qc & (original == typ)).sum()),
                     'program_concordant_candidates': int(selected_mask.sum()),
                     'minimum_reference_cluster_Jaccard': min(v['Jaccard'] for v in matches)})
    pd.DataFrame(rows).to_csv(o / 'identity_concordance_summary.tsv', sep='\t', index=False)
    pd.DataFrame(match_rows).to_csv(o / 'cluster_specific_stability.tsv', sep='\t', index=False)
    c['identity_candidate'] = label; c['QC_keep'] = qc; c['RNA_total'] = total; c['RNA_genes'] = ng
    c.to_parquet(o / 'identity_candidates_PRIVATE.parquet')
    spatial = base.spatial_tests(r, xy, total, label, qc, pos, ids, o)
    original_table = pd.read_csv(r / 'followup_v1/spatial_complete96.tsv', sep='\t')
    strict_table = pd.read_csv(o / 'spatial_complete96.tsv', sep='\t')
    original_table['view'] = 'original_marker_candidates'
    strict_table['view'] = 'RNA_program_concordant_sensitivity'
    combined = pd.concat([original_table, strict_table], ignore_index=True)
    assert len(combined) == 192
    combined['BH_q_joint192'] = base.fdr(combined.p_one_sided_conditional)
    combined['BY_q_joint192'] = base.fdr(combined.p_one_sided_conditional, by=True)
    combined.to_csv(o / 'complete192_two_views.tsv', sep='\t', index=False, na_rep='NA')
    rec = {'UTC': datetime.now(timezone.utc).isoformat(), 'status': 'sensitivity_computed_audit_pending',
           'identity': rows, 'spatial': spatial, 'joint_hypotheses': 192,
           'joint_BH_lt05': int((combined.BH_q_joint192 < .05).sum()),
           'joint_BY_lt05': int((combined.BY_q_joint192 < .05).sum()),
           'original96_source_SHA_unchanged': hashlib.sha256((r / 'followup_v1/spatial_complete96.tsv').read_bytes()).hexdigest(),
           'no_new_fine_label_freeze': True, 'clinical_or_communication_models_added': False}
    atomic_json(o / 'complete.json', rec)
    print(json.dumps(rec), flush=True)
