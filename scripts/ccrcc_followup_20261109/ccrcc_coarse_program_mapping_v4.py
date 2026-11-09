"""Exploratory identity sensitivity, defined after original96 QC review.

Intersect original RNA rules with marker-supported coarse RNA programs that
agree in all three seeds, allowing each class to occupy multiple state clusters. No position, ACKR1 or spatial P value determines membership.
Preserve original96 and adjust the three explored views together (288 slots).
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

r = Path(sys.argv[1]).resolve(); o = r / 'program_mapping_v4'


def consensus(path, ref_clusters, name):
    obj = ad.read_h5ad(path, backed='r')
    ix = obj.obs.source_index_PRIVATE.to_numpy(int)
    reference = obj.obs.seed_0.astype(str).to_numpy()
    reference_program = np.isin(reference, ref_clusters)
    keep = reference_program.copy()
    rows = []
    for seed in [1, 2]:
        target = obj.obs['seed_' + str(seed)].astype(str).to_numpy()
        accepted = []
        for cluster in sorted(set(target), key=int):
            z = target == cluster
            overlap = int((reference_program & z).sum())
            purity = overlap / z.sum()
            selected = purity > .5
            if selected:
                accepted.append(cluster)
            rows.append({'program': name, 'other_seed': seed, 'target_cluster': cluster,
                         'target_cluster_cells': int(z.sum()), 'reference_program_overlap': overlap,
                         'program_purity': float(purity), 'mapped_to_program': bool(selected)})
        mapped_program = np.isin(target, accepted)
        intersection = int((reference_program & mapped_program).sum())
        union = int((reference_program | mapped_program).sum())
        rows.append({'program': name, 'other_seed': seed, 'target_cluster': 'PROGRAM_UNION',
                     'target_cluster_cells': int(mapped_program.sum()), 'reference_program_overlap': intersection,
                     'program_purity': intersection / max(1, int(mapped_program.sum())),
                     'mapped_to_program': True, 'Jaccard': intersection / union})
        keep &= mapped_program
    obj.file.close()
    return ix[keep], rows


with exclusive_lock(r / 'program_mapping_v4.lock'):
    o.mkdir(exist_ok=False)
    definitions = [
        ('B', 'identity_refinement_v2/B_ASC/RNA_only_refinement_PRIVATE.h5ad', ['0']),
        ('T', 'identity_refinement_v2/T_NK/RNA_only_refinement_PRIVATE.h5ad', ['0', '1']),
        ('NK', 'identity_refinement_v2/T_NK/RNA_only_refinement_PRIVATE.h5ad', ['2']),
        ('Blood_EC', 'followup_v1/clusters_vascular_stromal/RNA_only_cluster_object_PRIVATE.h5ad', ['0'])]
    atomic_json(o / 'design_before_v3_spatial_effects.json', {
        'UTC': datetime.now(timezone.utc).isoformat(),
        'original96_and_v3_diagnostic_effects_already_seen': True,
        'correction': 'v3 exact cluster-overlap intersection unnecessarily dropped some T states; v4 maps the coarse T union across any number of seed clusters',
        'exploratory_sensitivity_not_independent_confirmation': True,
        'membership': 'Original QC/rule identity AND marker-supported coarse program in all3 seeds; each target seed cluster maps by majority overlap (>0.5) with seed0 PROGRAM UNION, allowing one-to-many state splits',
        'programs': {n: {'parent_object': p, 'seed0_clusters': z} for n, p, z in definitions},
        'program_evidence': {'B': 'MS4A1/BANK1/CD19; CD20 protein',
                             'T': 'CD3/TRAC/CD2; CD3E protein',
                             'NK': 'GNLY/KLRD1/KLRC1/KLRF1/NKG7; lower within-antibody CD3E than T',
                             'Blood_EC': 'PODXL/PECAM1/CD34/ADGRL4/CD93/VWF/CLEC14A; CD31 protein'},
        'spatial_parameters': 'Unchanged complete96; original96 + v3 strict-cluster diagnostic96 + corrected coarse-program96; joint288 BH/BY',
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
                     'minimum_reference_cluster_Jaccard': min(v['Jaccard'] for v in matches if v['target_cluster'] == 'PROGRAM_UNION')})
    pd.DataFrame(rows).to_csv(o / 'identity_concordance_summary.tsv', sep='\t', index=False)
    pd.DataFrame(match_rows).to_csv(o / 'coarse_program_mapping.tsv', sep='\t', index=False)
    c['identity_candidate'] = label; c['QC_keep'] = qc; c['RNA_total'] = total; c['RNA_genes'] = ng
    c.to_parquet(o / 'identity_candidates_PRIVATE.parquet')
    spatial = base.spatial_tests(r, xy, total, label, qc, pos, ids, o)
    original_table = pd.read_csv(r / 'followup_v1/spatial_complete96.tsv', sep='\t')
    strict_table = pd.read_csv(o / 'spatial_complete96.tsv', sep='\t')
    original_table['view'] = 'original_marker_candidates'
    diagnostic = pd.read_csv(r / 'program_concordance_v3/spatial_complete96.tsv', sep='\t')
    diagnostic['view'] = 'v3_strict_cluster_diagnostic'
    strict_table['view'] = 'v4_coarse_program_concordant'
    combined = pd.concat([original_table, diagnostic, strict_table], ignore_index=True)
    assert len(combined) == 288
    combined['BH_q_joint288'] = base.fdr(combined.p_one_sided_conditional)
    combined['BY_q_joint288'] = base.fdr(combined.p_one_sided_conditional, by=True)
    combined.to_csv(o / 'complete288_three_views.tsv', sep='\t', index=False, na_rep='NA')
    rec = {'UTC': datetime.now(timezone.utc).isoformat(), 'status': 'sensitivity_computed_audit_pending',
           'identity': rows, 'spatial': spatial, 'joint_hypotheses': 288,
           'joint_BH_lt05': int((combined.BH_q_joint288 < .05).sum()),
           'joint_BY_lt05': int((combined.BY_q_joint288 < .05).sum()),
           'original96_source_SHA_unchanged': hashlib.sha256((r / 'followup_v1/spatial_complete96.tsv').read_bytes()).hexdigest(),
           'no_new_fine_label_freeze': True, 'clinical_or_communication_models_added': False}
    atomic_json(o / 'complete.json', rec)
    print(json.dumps(rec), flush=True)
