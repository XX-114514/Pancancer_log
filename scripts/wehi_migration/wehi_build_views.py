#!/usr/bin/env python3
"""WEHI: index migrated runs by analysis line and build symlink views; data are not moved.

Lines follow the project's three active analyses plus the shared atlas foundation:
  0_atlas_foundation      source qualification, native archive, QC, public ingest
  1_TME_annotation        new-TME annotation, integration, clustering, engine comparisons
  1_TME_CM                new-TME CoVarNet schemes, comparisons and robustness
  2_immune_atlas_CM       immune atlas CM02, clinical cohorts, earlier communication work
  3_state_programmes      state scoring, state CM, conditional state models, increments
  9_other                 anything not matched
The first matching rule wins. Within a family (name without date, version and numeric suffix)
the newest date/version is flagged, so the latest step of each family is easy to find.

Usage: wehi_build_views.py PC_ROOT [--spatial DIR]
Writes index/runs_by_line.tsv and rebuilds views/ (only symlinks live there).
"""
import argparse
import os
import re
import shutil

RULES = [
    # Collections that feed the new-TME atlas despite 'immune' in their name.
    ('0_atlas_foundation', r'existing_immune_ecology'),
    ('3_state_programmes', r'TME_conditional_|TME_state_|TME_increment|TME_identity_(state|conditional)'
     r'|TME_MOFA|TME_native_UCell|TME_frozen_state|TME_gene_overlap|TME_original_and_conditional'
     r'|predefined_state|TNK_state_pilot|DIALOGUE|signature_proxy|fixed_CD8_state'
     r'|effective_T_cytotoxic'),
    ('2_immune_atlas_CM', r'CM02|immune|immunotherapy|ICI_|ici_|ICB_|clinical_(review|preflight)'
     r'|TNBC|GSE120575|Treg_|CCL5|TNFSF9|cellpair|majorpair|Tcell_pair|chemokine|tumor_local'
     r'|stromal_receiver|auto_communication|patient_source|three_source_author|GSE305078'),
    ('1_TME_CM', r'TME_CM_|TME_all_CM|TME_block_rank'),
    ('1_TME_annotation', r'TME_|standard_matched|full_context_standard'),
    ('0_atlas_foundation', r'full_|source_|native_|guarded_|joined_native|library_qualification'
     r'|public_preprocessed|clinical_scrna_expansion|doublet_ambient|archive_|complete_source'
     r'|common_reference|fine_|pDC|disputed_|bounded_projection|representative_|safe_group'
     r'|original_h5|atlas_census|GSE16|existing_immune_ecology|phase[01]'),
]
NAME_RE = re.compile(r'^(\d{8})_(.+?)(?:_v(\d+))?$')

CURRENT = [
    ('immune_atlas_v1.zarr',
     'Project_v5/runs/20261013_immunotherapy_atlas_v1/outputs/immune_validation_authority_v1.zarr'),
    ('core_fullgene_v1.zarr',
     'Project_v3/annotation_rebuild/20260923_core_anndata_zarr_v1/objects/pancancer_core_fullgene_v1.zarr'),
    ('public_preprocessed_union_v1.zarr',
     'Project_v5/runs/20261012_public_preprocessed_zarr_v2/outputs/public_preprocessed_union_v1.zarr'),
    ('clinical_expansion_native_objects',
     'Project_v5/runs/20261010_clinical_scrna_expansion_v1/objects'),
    ('clinical_expansion_metadata', 'Project_v5/runs/20261010_clinical_scrna_expansion_v1/metadata'),
    ('TME_joined_native_atlas_delivery', 'Project_v5/runs/20261026_joined_native_atlas_delivery_v1'),
    ('frozen', 'Project_v5/frozen'),
    ('celltypist_models', 'resources/celltypist'),
    ('R_ExperimentHub', 'resources/R_ExperimentHub'),
    ('env_specs', 'envs/specs'),
]


def classify(stem):
    for line, pattern in RULES:
        if re.search(pattern, stem):
            return line
    return '9_other'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('pc_root')
    ap.add_argument('--spatial')
    a = ap.parse_args()
    pc = os.path.abspath(a.pc_root)
    runs_dir = os.path.join(pc, 'Project_v5', 'runs')
    rows = []
    for run in sorted(os.listdir(runs_dir)):
        m = NAME_RE.match(run)
        date, stem, ver = (m.group(1), m.group(2), int(m.group(3) or 0)) if m else ('', run, 0)
        family = re.sub(r'_\d{1,3}$', '', stem)
        rows.append({'run': run, 'date': date, 'version': ver, 'family': family,
                     'line': classify(stem)})
    newest = {}
    for r in rows:
        key = (r['line'], r['family'])
        if key not in newest or (r['date'], r['version'], r['run']) > newest[key]:
            newest[key] = (r['date'], r['version'], r['run'])

    os.makedirs(os.path.join(pc, 'index'), exist_ok=True)
    with open(os.path.join(pc, 'index', 'runs_by_line.tsv'), 'w') as f:
        f.write('line\tfamily\tdate\tversion\tlatest_in_family\trun\n')
        for r in sorted(rows, key=lambda r: (r['line'], r['family'], r['date'], r['version'])):
            latest = newest[(r['line'], r['family'])][2] == r['run']
            f.write('%s\t%s\t%s\t%d\t%s\t%s\n' % (r['line'], r['family'], r['date'], r['version'],
                                                  'yes' if latest else '', r['run']))

    views = os.path.join(pc, 'views')
    if os.path.isdir(views):
        for root, dirs, files in os.walk(views):
            for n in files + dirs:
                p = os.path.join(root, n)
                if not os.path.islink(p) and os.path.isfile(p):
                    raise SystemExit('views/ holds a real file, refusing to rebuild: ' + p)
        shutil.rmtree(views)
    for r in rows:
        for sub in ('all', 'latest') if newest[(r['line'], r['family'])][2] == r['run'] else ('all',):
            d = os.path.join(views, r['line'], sub)
            os.makedirs(d, exist_ok=True)
            os.symlink(os.path.relpath(os.path.join(runs_dir, r['run']), d), os.path.join(d, r['run']))
    d = os.path.join(views, 'current')
    os.makedirs(d, exist_ok=True)
    for name, rel in CURRENT:
        target = os.path.join(pc, rel)
        if os.path.exists(target):
            os.symlink(os.path.relpath(target, d), os.path.join(d, name))
    if a.spatial:
        os.symlink(os.path.abspath(a.spatial), os.path.join(views, '4_WEHI_spatial'))

    counts = {}
    for r in rows:
        counts[r['line']] = counts.get(r['line'], 0) + 1
    for line in sorted(counts):
        n_latest = sum(1 for k in newest if k[0] == line)
        print('%-20s %5d runs, %4d families' % (line, counts[line], n_latest))


if __name__ == '__main__':
    main()
