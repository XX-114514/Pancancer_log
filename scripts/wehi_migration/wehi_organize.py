#!/usr/bin/env python3
"""WEHI: move verified migration items into one tree that mirrors the old server layout.

Old run configs hold absolute old-server paths; mirroring the layout keeps them resolvable by
one prefix substitution (see PATH_MAP.tsv written next to the tree). Moves are renames inside
one filesystem: nothing is copied or rewritten. Every move is logged to
migration/relocation.tsv so wehi_verify_relocated.py can re-check all files at their new place.

Usage: wehi_organize.py DATA_ROOT STAGE_ROOT PC_ROOT [--execute]
Without --execute only prints the plan. Refuses to overwrite any existing target.
"""
import argparse
import os
import shutil
import sys

ZARRS = [
    ('immune_validation_authority_v1.zarr',
     'Project_v5/runs/20261013_immunotherapy_atlas_v1/outputs/immune_validation_authority_v1.zarr'),
    ('public_preprocessed_union_v1.zarr',
     'Project_v5/runs/20261012_public_preprocessed_zarr_v2/outputs/public_preprocessed_union_v1.zarr'),
    ('pancancer_core_fullgene_v1.zarr',
     'Project_v3/annotation_rebuild/20260923_core_anndata_zarr_v1/objects/pancancer_core_fullgene_v1.zarr'),
]
ITEMS = ['m0_records', 'immune_authority_zarr', 'core_fullgene_zarr', 'public_union_zarr',
         'm2_working_layer', 'm4a_inputs', 'm4b_runs', 'm4c_env']
RECEIPT_FILES = ['receipt.tsv', 'meta.tsv', 'files.sha256', 'chunks.sha256', 'files.wehi.sha256',
                 'fetch.log']


def plan(data, pc):
    moves = []

    def add(src, dst):
        moves.append((os.path.join(data, src), os.path.join(pc, dst)))

    def each(src_dir, dst_dir):
        for name in sorted(os.listdir(os.path.join(data, src_dir))):
            add(os.path.join(src_dir, name), os.path.join(dst_dir, name))

    m0 = 'm0_tree/Project_v5_records'
    for name in sorted(os.listdir(os.path.join(data, m0))):
        if name == 'runs':
            # Config/provenance snapshots of runs that now exist in full; kept as provenance.
            add(os.path.join(m0, 'runs'), 'migration/m0_run_records')
        else:
            add(os.path.join(m0, name), os.path.join('Project_v5', name))
    each('m0_tree/extra', 'resources/wheels')
    for runs in ('m2_tree/Project_v5_runs', 'm4/m4a_inputs_tree/Project_v5_runs',
                 'm4/m4b_runs_tree/Project_v5_runs'):
        each(runs, 'Project_v5/runs')
    # The clinical expansion run arrives partially (referenced objects and non-raw parts).
    add('m4/m4a_inputs_tree/clinical_scrna_expansion_v1',
        'Project_v5/runs/20261010_clinical_scrna_expansion_v1')
    each('m4/m4a_inputs_tree/resources', 'resources')
    each('m2_tree/extra', 'Project_v3/annotation_rebuild')
    each('m4/m4c_env_tree', 'envs/specs')
    for src, dst in ZARRS:
        add(src, dst)
    for f in ('closure.runs.tsv', 'selected_runs.txt', 'symlinks.tsv', 'm2_exclude.txt'):
        add(os.path.join('m2_tree', f), os.path.join('migration/items/m2_working_layer', f))
    for item, tree in (('m4a_inputs', 'm4/m4a_inputs_tree'), ('m4b_runs', 'm4/m4b_runs_tree')):
        for f in ('source_list.tsv', 'symlinks.tsv'):
            add(os.path.join(tree, f), os.path.join('migration/items', item, f))
    return moves


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('data_root')
    ap.add_argument('stage_root')
    ap.add_argument('pc_root')
    ap.add_argument('--execute', action='store_true')
    a = ap.parse_args()
    data, stage, pc = (os.path.abspath(p) for p in (a.data_root, a.stage_root, a.pc_root))
    moves = plan(data, pc)

    problems = []
    targets = set()
    for src, dst in moves:
        if not os.path.lexists(src):
            problems.append('missing source: ' + src)
        if os.path.lexists(dst) or dst in targets:
            problems.append('target exists: ' + dst)
        targets.add(dst)
    for item in ITEMS:
        if not os.path.isfile(os.path.join(stage, item, 'receipt.tsv')):
            problems.append('missing receipt: ' + item)
    if problems:
        print('\n'.join(problems), file=sys.stderr)
        sys.exit(1)
    print('%d moves planned' % len(moves))
    if not a.execute:
        for src, dst in moves[:20]:
            print('%s -> %s' % (os.path.relpath(src, data), os.path.relpath(dst, pc)))
        return

    os.makedirs(os.path.join(pc, 'migration'), exist_ok=True)
    log = open(os.path.join(pc, 'migration', 'relocation.tsv'), 'w')
    log.write('old_path\tnew_path\n')
    for src, dst in moves:
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        os.rename(src, dst)
        log.write('%s\t%s\n' % (src, dst))
        log.flush()
    log.close()
    for item in ITEMS:
        out = os.path.join(pc, 'migration', 'items', item)
        os.makedirs(out, exist_ok=True)
        for f in RECEIPT_FILES:
            p = os.path.join(stage, item, f)
            if os.path.isfile(p):
                shutil.copy2(p, out)
    print('done; relocation log in migration/relocation.tsv')


if __name__ == '__main__':
    main()
