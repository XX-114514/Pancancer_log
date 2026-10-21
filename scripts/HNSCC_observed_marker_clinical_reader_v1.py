"""Read-only reviewed RNA/clinical bridge. Returned individual arrays stay in RAM.

The caller injects pinned native_counts and evidence; no file output or network
downloads are performed here. All original source libraries are checked.
"""
from collections import defaultdict
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
from scipy.sparse import csr_matrix


def read_source(root):
    root = Path(root)
    table = root/'HNSCC_author_clinical_table_user_20261014_v1/private_key_and_response_audit_v1/all57_private_author_table_join.tsv'
    clinical = root/'HNSCC_frozen_MCP_clinical_preparation_20261014_v1/private_patient_clinical_definition.tsv'
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    assert sha(table) == '9942706870518c20c87a45fed44caff6bec22a1a512b6f194dc17d406df5cd22'
    assert sha(clinical) == 'ce3c924ead635e51d322669e9b3504b5abee38a30044394eb9c4f7778f949db6'
    tab = pd.read_csv(table, sep='\t', dtype=str, keep_default_na=False)
    cl = pd.read_csv(clinical, sep='\t', dtype=str, keep_default_na=False)
    assert tab.library_id.is_unique and cl.patient_id.is_unique
    cl = cl.set_index('patient_id')
    bridge = tab.set_index('library_id')
    atlas = root.parents[1]/'runs/20261013_immunotherapy_atlas_v1'
    z = zarr.open_group(str(atlas/'outputs/immune_validation_authority_v1.zarr'), mode='r')
    ds = np.asarray(read_elem(z['obs/dataset'])).astype(str)
    selected = np.flatnonzero(ds == 'GSE286827')
    library = pd.Series(np.asarray(read_elem(z['obs/library_id'])).astype(str)[selected])
    phase = pd.Series(np.asarray(read_elem(z['obs/timepoint'])).astype(str)[selected])
    assert set(library) == set(bridge.index)
    assert np.array_equal(phase, library.map(bridge.timepoint_registry))
    patient = library.map(bridge.patient_registry)
    assert patient.notna().all() and set(patient) <= set(cl.index)
    assert (tab.groupby(['patient_registry','timepoint_registry']).size() == 1).all()
    for _, a in tab.iterrows():
        c = cl.loc[a.patient_registry]
        assert float(a['Tumor.regression']) == float(c.regression)
        assert a.neo_type == c.treatment and a['HPV/p16'] == c.HPV
        assert 0 <= float(c.regression) <= 100
    var = read_elem(z['var'])
    sources = list(map(str, read_elem(z['uns/source_order'])))
    measured = np.asarray(read_elem(z['varm/measured_in_source']), bool)[:, sources.index('GSE286827')]
    human = measured & var.feature_identity_status.astype(str).eq('stable_Ensembl_ID').to_numpy()
    lookup = defaultdict(set)
    for col in np.flatnonzero(human):
        r = var.iloc[col]
        for value in {str(var.index[col]), str(r.gene_symbol)} | set(str(r.source_gene_aliases).split('|')):
            lookup[value].add(int(col))
    columns = []
    for gene in evidence.GENES:
        assert len(lookup[gene]) == 1
        columns.append(next(iter(lookup[gene])))
    assert len(set(columns)) == len(columns)
    axis_sha = hashlib.sha256(json.dumps(dict(zip(evidence.GENES,columns)),sort_keys=True).encode()).hexdigest()
    assert axis_sha == '0ec80d4ee65c4b425ad440e9bf31ab0ffbfd945fa8f288b1974279f2b7bf0ded'
    mt = np.flatnonzero(human & var.gene_symbol.astype(str).str.startswith('MT-').to_numpy())
    assert len(mt) == 13 and not set(mt) & set(columns)
    counts, total, detected = native_counts(z, selected, columns+mt.tolist(), human, measured, lambda v:None)
    mito = counts[:,len(columns):].sum(axis=1)
    counts = counts[:,:len(columns)]
    assert (mito <= total).all()
    result = []
    coverage = []
    for qc in ['human_RNA_QC500_250_MT20','all_retained_authority_cells']:
        good = ((total >= 500) & (detected >= 250) & (mito <= .2*total) & (total > 0)) if qc.startswith('human') else np.ones(len(total),bool)
        flags = evidence.vector_evidence(counts, good, 'dual_regulatory_evidence')
        np.testing.assert_array_equal(flags,evidence.scalar_evidence(counts,good,'dual_regulatory_evidence'))
        np.testing.assert_array_equal(flags[:,1:5].sum(axis=1),good.astype(int))
        for ph in ['pre','post','post2']:
            cells = phase.eq(ph).to_numpy()
            librows = tab.loc[tab.timepoint_registry.eq(ph)]
            emitted = 0
            for _, entry in librows.iterrows():
                take = library.eq(entry.library_id).to_numpy()
                assert np.any(take) and (cells[take]).all()
                observed = flags[take].sum(axis=0).astype(int)
                quality = []
                for k in [1,2]:
                    parent = take & flags[:,k]
                    n = int(parent.sum())
                    quality.append([float(np.mean(np.log1p(total[parent]))) if n else None,
                                    float(np.mean(np.log1p(detected[parent]))) if n else None,
                                    float(np.log1p(n)), n/int(observed[0]) if observed[0] else None])
                c = cl.loc[entry.patient_registry]
                assert c.treatment in ['D','D+T'] and c.HPV in ['positive','negative']
                # Order, not identity, is sufficient for within-phase inference.
                result.append(dict(phase=ph,QC=qc,metrics=observed.tolist(),quality=quality,
                    regression=float(c.regression),clinical=[int(c.treatment=='D+T'),int(c.HPV=='positive')]))
                emitted += 1
            coverage.append(dict(phase=ph,QC=qc,reviewed_patients=emitted,
                retained_cell_observations=int(cells.sum()),QC_cell_observations=int((cells&good).sum())))
    return dict(rows=result,source_coverage=coverage,source='GSE286827',
        source_observations=len(selected),unique_reviewed_patients=len(set(patient)),
        clinical_table_SHA256=sha(table),clinical_definition_SHA256=sha(clinical),marker_column_map_SHA256=axis_sha,
        authority_release_SHA256=sha(atlas/'manifests/immune_authority_release_public.json'),
        scalar_cell_rules_and_raw_CSR_marker_UMI_checks=True,
        all_original_source_library_keys_matched=True,no_files_written=True)
