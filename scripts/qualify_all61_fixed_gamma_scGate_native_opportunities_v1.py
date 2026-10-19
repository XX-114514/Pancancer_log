"""All existing source axes for frozen pilot models; metadata only."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd
from external_functional_coverage_inputs_v1 import read_inputs
from NSCLC_external_lineage_inputs_v2 import load
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT']); guard = resource(root, 1024**2)
    fixed = root / 'gamma_scGate_native_model_gene_qualification_20261019_v1/private_complete_fixed_two_model_native_qualification.json'
    assert sha(fixed) == '1f746d64ff2eeb601f9b8bfd42802b4a4b96f218dbff1fb044c697258ead7778'
    runtime = root / 'gamma_scGate_official_runtime_qualification_20261019_v1'
    assert sha(runtime / 'private_terminal.json') == '8b69085d7dbe3be4b7342d58511066b2804eef4f7618a362112f4a1e24379058'
    j = json.loads(fixed.read_text()); reference = {m['collection']: {s['name']: set(s['genes']) for s in m['signatures']} for m in j['models']}
    own = set(pd.read_csv(root / 'GSE243013_genes.csv.gz', keep_default_na=False).geneSymbol.str.strip().str.upper())
    full = set.union(*reference['HiTME'].values()); pbmc = set.union(*reference['PBMC'].values())
    nkfile = runtime / 'scGate_models-41a45cd3f8bb5f5a7daf21ec276f6a726f6ee0d4/human/HiTME/NK_scGate_Model.tsv'
    nk = pd.read_csv(nkfile, sep='\t', keep_default_na=False)
    models = dict(HiTME_full=full, PBMC_exact=pbmc, HiTME_native_subset=full & own, NK_HiTME_exact=set.union(*(reference['HiTME'][name] for name in nk.name)))
    for genes in reference['HiTME'].values(): assert len(genes & own) >= 1 and len(genes & own) / len(genes) >= .8
    data = load(); axes = read_inputs(root)
    out = root / 'all61_fixed_gamma_scGate_native_opportunities_20261019_v1'; assert not out.exists(); out.mkdir(mode=0o700)
    save = lambda name, value: (out / name).write_text(json.dumps(value, indent=2) + '\n')
    save('plan_before_all_source_model_axis_qualification.json', dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), fixed_model_scope_sha256=sha(fixed), runtime_sha256=sha(runtime / 'private_terminal.json'), NK_definition_sha256=sha(nkfile), resource_guard=guard, input_metadata_hashes=data['hashes'], input_axis_hashes=axes['input_sha256'], models={k: sorted(v) for k, v in models.items()}, scope='All61 source-arm physical axes and every existing material/capture/phase context. The HiTME measured subset is frozen by original NSCLC measurability, never re-trimmed per new source. Complete-model gates require all genes uniquely measurable with no column collision. Physical eligibility is not RNA identity or clinical acceptance.'))
    records = []
    for (arm, dataset), axis in sorted(axes['source_axes'].items()):
        for name, genes in models.items():
            matches = {g: axis['lookup'].get(g, set()) for g in genes}
            unique = {g: next(iter(v)) for g, v in matches.items() if len(v) == 1}
            collisions = len(unique) - len(set(unique.values()))
            records.append(dict(arm=arm, dataset=dataset, model=name, fixed_features=len(genes), uniquely_mapped_features=len(unique), unmapped_features=sum(not v for v in matches.values()), ambiguous_features=sum(len(v) > 1 for v in matches.values()), column_collisions=collisions, all_fixed_features_uniquely_measurable=len(unique) == len(genes) and collisions == 0, required_missing_or_ambiguous_features=json.dumps(sorted(g for g, v in matches.items() if len(v) != 1)), source_dependent_feature_deletion=False, RNA_identity_or_clinical_accepted=False))
    table = pd.DataFrame(records); assert len(table) == 244 and not table.duplicated(['arm', 'dataset', 'model']).any()
    table.to_csv(out / 'private_complete61_four_fixed_model_physical_opportunities.tsv', sep='\t', index=False)
    lib = data['lib'].copy()
    lymphoid = {'T cell', 'T/NK lymphoid', 'Lymphoid cell', 'T', 'NK', 'T_NK_boundary'}
    c = data['counts']; parents = c[c.major.isin(lymphoid)].groupby(['arm', 'dataset', 'library_id']).n_cells.sum()
    lib['frozen_lymphoid_parent_cells'] = parents.reindex(pd.MultiIndex.from_frame(lib[['arm', 'dataset', 'library_id']]), fill_value=0).to_numpy(int)
    context = ['arm', 'dataset', 'cancer_label', 'material_class', 'phase_class', 'capture_scope', 'ex_vivo']
    contexts = []
    for key, g in lib.groupby(context, sort=True, dropna=False):
        reviewed = g.reviewed_patient.astype(bool)
        contexts.append(dict(zip(context, key), libraries=len(g), canonical_native_cells=int(g.annotated_scope_cells.sum()), frozen_lymphoid_parent_cells=int(g.frozen_lymphoid_parent_cells.sum()), reviewed_patient_blocks=int(g.loc[reviewed, 'patient_block'].nunique()), unreviewed_library_units=int((~reviewed).sum()), candidate_gamma_cells_not_yet_measured=True, clinical_or_identity_adoption=False))
    contexts = pd.DataFrame(contexts)
    assert contexts.libraries.sum() == len(lib) == 1690 and contexts.canonical_native_cells.sum() == lib.annotated_scope_cells.sum()
    contexts.to_csv(out / 'private_complete_source_material_capture_phase_contexts.tsv', sep='\t', index=False)
    save('private_all_fixed_model_definitions_and_input_hashes.json', dict(models={k: sorted(v) for k, v in models.items()}, input_metadata_hashes=data['hashes'], input_axis_hashes=axes['input_sha256']))
    assert sum(p.stat().st_size for p in out.iterdir() if p.is_file()) < 1024**2
    receipt = dict(status='COMPLETE_ALL61_FIXED_MODEL_PHYSICAL_AXES_AND_ALL_MATERIAL_CONTEXTS_NO_RNA_OR_CLINICAL_ADOPTION', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), source_arm_model_positions=244, libraries=1690, resource_guard=guard, files={p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()}, RNA_read=False, clinical_or_identity_adoption=False, scientific_outputs_public=False)
    save('private_terminal.json', receipt)
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
