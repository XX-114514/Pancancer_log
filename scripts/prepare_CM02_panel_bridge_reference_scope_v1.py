"""Freeze source-native, reviewed-case reference selection before model effects."""
import argparse
import hashlib
import json
import shutil
import subprocess
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem

CORE = ['FOXP3','IL2RA','CTLA4','IKZF2','TNFRSF18','GZMB','PRF1','CTSW','NKG7','CD8A','CD8B','CD3D','CD3E','TRAC']
AXES = ['TGFB1','TGFBR1','TGFBR2','CCL5','CCR5']
MAJORS = dict(T='T', NK='NK', Myeloid='Myeloid', Neutrophil='Myeloid', Mast='Myeloid',
              B='B_Plasma', Plasma='B_Plasma', Epithelial='Nonimmune', Fibroblast='Nonimmune',
              Endothelial='Nonimmune', Pericyte_SMC='Nonimmune', Erythroid='Nonimmune',
              Glial='Nonimmune', Neuron='Nonimmune', Oligodendrocyte='Nonimmune', Melanocytic='Nonimmune')


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def budget(root, additional):
    fs = subprocess.check_output(['df','-P',str(root)], text=True).splitlines()[-1].split()[0]
    rows = [s.split() for s in subprocess.check_output(['quota','-g','-v','-w','USER002'],text=True).splitlines() if s.split() and s.split()[0] == fs]
    assert len(rows) == 1
    remaining = (int(rows[0][3])-int(rows[0][1].rstrip('*'))) * 1024
    assert min(remaining,shutil.disk_usage(root).free) - additional >= 450 * 1024**3
    return dict(group_remaining_bytes=remaining, reserve_bytes=450*1024**3,
                estimated_additional_bytes=additional, reserve_change_explicitly_user_authorized=True)


def main(a):
    resource = budget(a.root, 16*1024**2)
    z = zarr.open_group(str(a.reference), mode='r')
    obs = read_elem(z['obs'])
    var = read_elem(z['var'])
    sources = list(map(str, read_elem(z['uns/source_order'])))
    measured = np.asarray(read_elem(z['varm/measured_in_source']), bool)
    official = json.loads(a.panel.read_text())['MERFISH_genes']
    assert len(official) == len(set(official)) == 479
    # Reviewed identity flag does not erase an unresolved official-link status.
    okay = obs.patient_identity_verified.astype(str).str.lower().isin(['true','1'])
    okay &= obs.patient_id_source_status.astype(str).eq('exact_original_library_official_metadata_verified')
    okay &= ~obs.annotation_exclude_from_primary.astype(str).str.lower().isin(['true','1'])
    okay &= obs.major.astype(str).isin(MAJORS)
    okay &= ~obs.patient_id.astype(str).str.upper().isin(['','NAN','NONE','NA','UNKNOWN','__SOURCE_VALUE_MISSING__'])
    selected = obs.loc[okay, ['dataset','patient_id','library_id','major','cancer_type','timepoint','assay','annotation_status']].copy()
    selected['original_row'] = np.flatnonzero(okay)
    selected['lineage'] = selected.major.astype(str).map(MAJORS)
    selected['patient_block'] = selected.dataset.astype(str) + '::' + selected.patient_id.astype(str)
    rng = np.random.default_rng(20261020)
    parts = []
    coverage = []
    for (source, patient, lineage), frame in selected.groupby(['dataset','patient_id','lineage'], observed=True, sort=True):
        cap = 2000 if lineage == 'T' else 1000
        ids = rng.permutation(len(frame))[:min(cap,len(frame))]
        chosen = frame.iloc[ids].copy()
        chosen['lineage_evaluation_row'] = np.arange(len(chosen)) < 1000
        parts.append(chosen)
        coverage.append(dict(source=source, patient_block=str(source)+'::'+str(patient), lineage=lineage,
                             all_eligible_cells=len(frame), selected_cells=len(chosen)))
    selected = pd.concat(parts).sort_values('original_row').reset_index(drop=True)
    assert selected.original_row.is_unique
    used_sources = sorted(selected.dataset.astype(str).unique())
    maps = {}
    for source in used_sources:
        lookup = defaultdict(set)
        si = sources.index(source)
        for col in np.flatnonzero(measured[:,si]):
            if str(var.iloc[col].feature_identity_status) == 'ambiguous_source_feature':
                continue
            labels = {str(var.index[col]), str(var.iloc[col].gene_symbol)}
            if str(var.iloc[col].feature_identity_status) == 'stable_Ensembl_ID':
                labels.update(str(var.iloc[col].source_gene_aliases).split('|'))
            for label in labels - {'','nan','None'}:
                lookup[label].add(int(col))
        mapping = {g: next(iter(lookup[g])) if len(lookup[g]) == 1 else -1 for g in set(official+CORE+AXES)}
        collisions = {v for v in mapping.values() if v >= 0 and list(mapping.values()).count(v)>1}
        mapping = {g: (-1 if v in collisions else v) for g,v in mapping.items()}
        maps[source] = mapping
    shared = sorted(g for g in official if all(maps[s][g]>=0 for s in used_sources))
    assert len(shared) >= 50
    lineage_genes = sorted(set(shared)-set(AXES))
    programme_genes = sorted(set(shared)-set(CORE)-set(AXES))
    estimated = int(len(selected)*len(shared)*4*3 + len(selected)*256 + 2*1024**3)
    actual_budget = budget(a.root, estimated)
    out = a.root/'reference_scope_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    selected.to_csv(out/'private_selected_reference_rows.tsv.gz',sep='\t',index=False,compression='gzip')
    pd.DataFrame(coverage).to_csv(out/'private_all_eligible_case_lineage_selection.tsv',sep='\t',index=False)
    contract = dict(created_utc=datetime.now(timezone.utc).isoformat(),
                    producer_sha256=sha(Path(__file__)), panel_sha256=sha(a.panel),
                    reference_shape=list(z['layers/counts'].attrs['shape']),
                    reference_release_sha256=sha(a.reference.parent.parent/'manifests/immune_authority_release_public.json'),
                    sample_selection='seed20261020;max1000/case/lineage;T2000/case for programmes;all reviewed official-linked cases;no spatial/outcomes',
                    known_major_labels_are_reference_proposals_not_independent_biological_truth=True,
                    sources=used_sources, source_native_gene_columns=maps,
                    shared_panel_genes=shared, lineage_feature_genes=lineage_genes,
                    programme_feature_genes=programme_genes, core14=CORE, LR_excluded_genes=AXES,
                    UCell_reference_signatures=dict(regulatory=CORE[:5],cytotoxic=CORE[5:11]),
                    normalization='shared measured biological panel CP10K/log1p;StandardScaler training fold only',
                    rank_cap=1500, reserve=actual_budget,
                    private_selected_rows_sha256=sha(out/'private_selected_reference_rows.tsv.gz'),
                    hyperparameters=dict(logistic_C=[.1,1.,10.],ridge_alpha=[.1,1.,10.,100.],inner_case_folds=3,seed=20261020),
                    eligibility=dict(min_evaluation_cells=100,min_evaluation_cases=5,min_studies=3,
                                     lineage_macroF1=.8,adopted_class_precision=.9,programme_median_rho=.5,programme_each_study_rho=.3,
                                     maximum_probability_threshold=.8),
                    no_model_fit_or_spatial_effect=True)
    (out/'private_frozen_reference_contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    public = dict(status='REVIEWED_CASE_REFERENCE_SELECTION_AND_NATIVE_GENE_SCOPE_FROZEN_NO_RNA_MODEL',
                  finished_utc=datetime.now(timezone.utc).isoformat(),selected_cells=len(selected),
                  sources=len(used_sources),private_contract_sha256=sha(out/'private_frozen_reference_contract.json'),
                  producer_sha256=sha(Path(__file__)),resource=actual_budget,
                  class_labels_are_proposals_not_biological_gold_standard=True)
    (out/'scope_technical_public.json').write_text(json.dumps(public,indent=2)+'\n')
    print(json.dumps(public),flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for field in ['root','reference','panel']:
        p.add_argument('--'+field,type=Path,required=True)
    main(p.parse_args())
