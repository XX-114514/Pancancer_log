#!/usr/bin/env python3
"""Collect an explicit allowlist of small, aggregate server evidence, never cell rows.

The caller supplies existing roots. Output must be a NEW directory. File hashes refer
to original bytes; published copies replace physical roots with logical placeholders.
Directory presence and producer completion reports do not establish scientific validity.
"""
import argparse
import collections
import csv
import datetime as dt
import hashlib
import json
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for name in ('pc-root', 'workspace', 'source-repo', 'output'):
        ap.add_argument('--' + name, required=True, type=Path)
    a = ap.parse_args()
    a.output.mkdir(parents=True, exist_ok=False)
    captured = dt.datetime.now(dt.timezone.utc).isoformat()
    roots = [(str(a.source_repo), '${CODE_ROOT}'),
             (str(a.workspace), '${WORKSPACE}'), (str(a.pc_root), '${PC_ROOT}'),
             (str(a.pc_root.parent), '${SCRATCH_ROOT}'),
             (str(a.workspace.parent), '${HOME_ROOT}')]
    evidence = []

    def logical(s):
        for physical, alias in roots:
            s = s.replace(physical, alias)
        return s

    def table(name, rows, columns):
        with (a.output / name).open('w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=columns, delimiter='\t', lineterminator='\n')
            w.writeheader()
            for row in rows:
                row = dict(row)
                if not row.get(columns[-1]):
                    row[columns[-1]] = 'NA'
                w.writerow(row)

    def record(p, role, dest='not_published'):
        if not p.is_file():
            return False
        size = p.stat().st_size
        if size > 5 * 1024 * 1024:
            raise ValueError('Evidence file unexpectedly large: ' + logical(str(p)))
        raw = p.read_bytes()
        evidence.append(dict(logical_path=logical(str(p)), role=role, bytes=size,
                             sha256=hashlib.sha256(raw).hexdigest(),
                             mtime_utc=dt.datetime.fromtimestamp(p.stat().st_mtime, dt.timezone.utc).isoformat(),
                             captured_utc=captured, published_copy=dest, published_sha256=''))
        return True

    def copy(p, dest, role='aggregate_result'):
        if record(p, role, 'tables/' + dest):
            out = a.output / 'tables' / dest
            out.parent.mkdir(parents=True, exist_ok=True)
            content = logical(p.read_text())
            if p.suffix == '.tsv':
                content = '\n'.join(line + 'NA' if line.endswith('\t') else line for line in content.splitlines()) + '\n'
            out.write_text(content)
            evidence[-1]['published_sha256'] = hashlib.sha256(out.read_bytes()).hexdigest()

    inv = a.workspace / 'investigations'
    comp = inv / '20261103_V5_V7_recompute'
    for p in sorted((comp / 'tables').glob('*')):
        if p.suffix in ('.tsv', '.json'):
            copy(p, 'v5_v7/' + p.name)
    old = inv / '20261103_0906_blisa'
    for name in ('artifact_integrity.tsv', 'CM01_CM07_CM08_members.tsv',
                 'old_new_cm_alignment.tsv', 'focus_axes_definition.tsv',
                 'focus_axes_replication_summary.tsv', 'panel_axis_coverage.tsv', 'EVIDENCE_MATRIX.tsv'):
        copy(old / 'tables' / name, 'historical_0906/' + name)
    nmf = a.pc_root / 'Project_v5/runs/20261104_V7_NMF_benchmark_v1'
    for name in ('method_summary.tsv', 'ref_cm_matches.tsv', 'seed_stability.tsv',
                 'cross_method_member_cos.tsv', 'unmatched_factors.tsv', 'source_eta2_runtime.tsv'):
        copy(nmf / 'compare' / name, 'nmf/' + name)
    cm = a.pc_root / 'Project_v5/runs/20261103_V7_CM_blisa_v1'
    for branch in ('primary_strict', 'legacy_all_display', 'primary_strict_K9only'):
        for name in ('CoVarNet_execution_report_v1.json', 'rank_survey_v1.tsv'):
            copy(cm / 'outputs' / ('covarnet_' + branch) / name, 'v7_cm/' + branch + '/' + name)
    for p in sorted((cm / 'results').glob('legacy_v7_v5_cm_alignment_seed_*.tsv')):
        copy(p, 'v7_cm/' + p.name)
    for p in sorted((comp / 'v7_cm_blisa/tables').glob('v7_v5_cm_alignment_seed_*.tsv')):
        copy(p, 'v7_cm/' + p.name)
    for p in sorted((cm / 'results').glob('cm_blisa_*/CM_summary_across_units.tsv')):
        copy(p, 'v7_spatial/' + p.parent.name + '.tsv')
    for p in sorted((cm / 'xenium').glob('panel_*_heldout_F1.tsv')):
        copy(p, 'v7_spatial/' + p.name)
    for name in ('graph_metrics.tsv', 'scib_metrics.tsv', 'times_cpu.json', 'times_scvi.json', 'times_scanorama.json'):
        copy(a.pc_root / 'tme_v1/batch_benchmark_v1' / name, 'batch_benchmark/' + name)
    for name in ('tme_summary.json', 'integration_metrics.json', 'malignancy_by_source.tsv'):
        copy(a.pc_root / 'tme_v1/run_v1/report' / name, 'tme_v1/' + name)
    for p in sorted((a.source_repo / 'inventories/tme_v1').glob('*')):
        # Existing reviewed repository tables remain authoritative; only hash them.
        if p.is_file():
            record(p, 'existing_repository_inventory')

    pilots = {}
    for rel in ('20261103_celltypist72_pilot/GSE138709_1000/summary.json',
                '20261103_marker_routed_annotation/GSE138709_v1/summary.json',
                '20261103_marker_routed_annotation/GSE138709_v2/summary.json'):
        p = inv / rel
        if record(p, 'source_level_pilot_summary', 'annotation_pilots.json'):
            pilots[rel] = json.loads(logical(p.read_text()))
    (a.output / 'annotation_pilots.json').write_text(json.dumps(pilots, indent=2) + '\n')
    for row in evidence:
        if row['published_copy'] == 'annotation_pilots.json':
            row['published_sha256'] = hashlib.sha256((a.output / 'annotation_pilots.json').read_bytes()).hexdigest()
    p = inv / '20261103_marker_routed_annotation/submissions/20261103_preflight_v2/plan.json'
    if record(p, 'historical_preflight', 'annotation_preflight.json'):
        plan = json.loads(p.read_text())
        public = {k: plan[k] for k in ('created_utc', 'submitted', 'source_counts', 'code_hashes') if k in plan}
        public['current_completion_not_inferred'] = True
        (a.output / 'annotation_preflight.json').write_text(json.dumps(public, indent=2) + '\n')
        evidence[-1]['published_sha256'] = hashlib.sha256((a.output / 'annotation_preflight.json').read_bytes()).hexdigest()
    copy(a.workspace / 'xenium_GSE274673/lessons/02_qc/summary.json', 'xenium_pilot/qc_summary.json')

    spatial = []
    sections = sorted(p for p in (cm / 'xenium').iterdir() if p.is_dir())
    for size in (25, 50, 100):
        complete = 0
        for section in sections:
            p = section / ('blisa_bin' + str(size)) / 'summary.tsv'
            if record(p, 'section_completion_metadata'):
                complete += 1
        spatial.append(dict(bin_um=size, sections_discovered=len(sections), summaries_present=complete,
                            scheduler_exit_verified=False, biological_validation=False))
    table('spatial_completion.tsv', spatial, list(spatial[0]))

    tasks = []
    for p in sorted((nmf / 'methods').glob('*/seed_*/info.json')):
        if record(p, 'method_runtime_receipt'):
            j = json.loads(p.read_text())
            tasks.append({k: j.get(k, '') for k in ('method', 'seed', 'K', 'input', 'package', 'seconds', 'starts_ok')})
    if tasks:
        table('nmf_tasks.tsv', tasks, list(tasks[0]))
    runroot = a.pc_root / 'Project_v5/runs'
    dirs = sorted(p.name for p in runroot.iterdir() if p.is_dir())
    counts = collections.Counter(n[:8] for n in dirs)
    table('run_directory_counts.tsv', [dict(name_date=k, directories=v, origin='presence_only_not_execution_date')
          for k, v in sorted(counts.items())], ['name_date', 'directories', 'origin'])
    table('recent_run_directories.tsv', [dict(run_id=n, execution_status='not_inferred_from_directory')
          for n in dirs if n[:8].isdigit() and n[:8] >= '20261101'], ['run_id', 'execution_status'])
    # Logs stay private. Export a narrow set of JSON progress fields only.
    p = a.pc_root / 'tme_v1/run_v2/integrate-tme-v1-30660430.out'
    progress = []
    if record(p, 'producer_progress_log'):
        for line in p.read_text().splitlines():
            try:
                j = json.loads(line)
            except ValueError:
                continue
            if j.get('step') in ('matrix', 'pca_done', 'harmony_vars', 'harmony_done', 'done', 'write_done'):
                progress.append({k: j[k] for k in ('t', 'step', 'shape', 'vars') if k in j})
    (a.output / 'tme_run_v2_progress.json').write_text(json.dumps(dict(events=progress,
        final_summary_present=(a.pc_root / 'tme_v1/run_v2/report/tme_summary.json').is_file(),
        scheduler_status='unverified'), indent=2) + '\n')
    # Original reports are hashes only: they contain physical paths and superseded states.
    for p in sorted(inv.glob('*/README.md')):
        record(p, 'historical_working_report')
    for p in [comp / 'V5_V7_COMPARISON_CN.md', comp / 'v7_cm_blisa/V7_CM_BLISA_RESULTS_CN.md',
              inv / '20261104_V7_NMF_benchmark/REPORT_CN.md']:
        record(p, 'historical_working_report')
    table('source_evidence.tsv', evidence, ['logical_path', 'role', 'bytes', 'sha256',
          'mtime_utc', 'captured_utc', 'published_copy', 'published_sha256'])
    (a.output / 'snapshot.json').write_text(json.dumps(dict(captured_utc=captured,
        captured_utc8=dt.datetime.now(dt.timezone(dt.timedelta(hours=8))).isoformat(),
        evidence_files=len(evidence), historical_run_directories=len(dirs),
        privacy='allowlisted aggregates; no cell rows, clinical rows, matrices or full logs',
        scope='named project roots only; no raw data recomputation; presence is not acceptance'), indent=2) + '\n')
    print(json.dumps(dict(evidence_files=len(evidence), run_directories=len(dirs), nmf_tasks=len(tasks))))


if __name__ == '__main__':
    main()
