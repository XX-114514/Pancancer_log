"""Run every predeclared official fit after full native input acceptance."""
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    cap = 128 * 1024**2;budget = resource(root, cap)
    inputs = root / 'immune_RNA_MCP_native_expression_inputs_20261018_v1'
    audit = root / 'immune_RNA_MCP_native_expression_input_audit_20261018_v1/private_terminal.json'
    assert sha(audit) == '58385bca98d751470d7f7026b0e6ff4ca7f1c584dacbfa6b7eb5ac644f7b5d4b'
    proof = json.loads(audit.read_text())
    assert sha(inputs / 'private_terminal.json') == proof['producer_terminal_sha256']
    for n, v in json.loads((inputs / 'private_terminal.json').read_text())['files'].items(): assert sha(inputs / n) == v['sha256']
    runtime = root.parent / '20261014_DIALOGUE_method_pilot_v1'
    archive = runtime / 'sources/DIALOGUE-9c146cc.tar.gz'
    assert sha(archive) == '2f2a8a6dfa2f47f2e9cce553310164d38c862a692e963b3e702bf04d7709e511'
    rscript = Path(os.environ.get('PANCANCER_DIALOGUE_RSCRIPT', str(Path(sys.executable).parents[2] / 'scrna_r/bin/Rscript')))
    assert rscript.is_file()
    worker = Path(__file__).with_name('run_immune_RNA_MCP_official_fit_v1.R');assert worker.is_file()
    out = root / 'immune_RNA_MCP_official_models_20261018_v1';assert not out.exists();out.mkdir(mode=0o700)
    fits = [dict(k=2, seed=s, source_mode='source_adjusted', role='primary' if s == 20261018 else 'seed_stability') for s in [20261018, 20261019, 20261020]]
    fits += [dict(k=3, seed=20261018, source_mode='source_adjusted', role='rank_sensitivity'),
             dict(k=2, seed=20261018, source_mode='source_permissive', role='source_context_sensitivity')]
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), worker_sha256=sha(worker), budget=budget,
        input_acceptance_sha256=sha(audit), input_terminal_sha256=sha(inputs / 'private_terminal.json'),
        official_archive_sha256=sha(archive), fits=fits, primary_anchor='K2 seed20261018 source_adjusted',
        primary_plan='Previously frozen four primary/stability/rank fits retained. Before any MCP result, add one separate source-permissive sensitivity reflecting user acceptance of source/cancer-associated ecology; no replacement of the primary model by the more favourable fit.',
        controls='Full native RNA/genes, original view fraction and original view cell count. Source indicator additionally in primary/stability/rank fits. No clinical phenotype supplied.',
        biological_scope='RNA-defined FOXP3-supported-only, effector-supported-only and original broad Myeloid; no formal cell identity change. Cases include mixed sampling context.',
        heldout='Entire TNBC source is excluded from feature training/PCA/DIALOGUE fitting. Identity review remains sample-case provisional, not untouched independent CM validation. Heldout scoring is a later frozen projection step.',
        storage='Official full.version TRUE, keep all original official intermediate/full/reduced outputs. Old real/synthetic output footprints reviewed before this cap. No original model/data or previous failed files deleted.',
        failures='Every planned fit retained in status, including errors/empty scores; no parameter/seed search on clinical or communication P.',
        output_cap_bytes=cap, stability_and_biological_or_clinical_validation='PENDING')
    (out / 'plan_before_any_official_MCP_fit.json').write_text(json.dumps(plan, indent=2) + '\n')
    statuses = []
    for item in fits:
        fit_id = 'K%d_seed%d_%s' % (item['k'], item['seed'], item['source_mode'])
        destination = out / fit_id;assert not destination.exists()
        used = sum(p.stat().st_size for p in out.rglob('*') if p.is_file())
        assert used < cap
        resource(root, cap - used)
        status = dict(item, fit_id=fit_id, started_UTC=datetime.now(timezone.utc).isoformat(), status='RUNNING')
        (out / 'private_current_fit.json').write_text(json.dumps(status, indent=2) + '\n')
        env = os.environ.copy();env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
        command = [str(rscript), str(worker), str(inputs), str(destination), str(runtime), str(item['k']), str(item['seed']), item['source_mode']]
        with (out / (fit_id + '_stdout_stderr.log')).open('xb') as log:
            try:
                result = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=6 * 3600)
                return_code = result.returncode
            except subprocess.TimeoutExpired:
                return_code = 124
        status['finished_UTC'] = datetime.now(timezone.utc).isoformat();status['exit_code'] = return_code
        summary = destination / 'private_fit_execution_summary.tsv'
        if return_code == 0 and summary.exists():
            row = pd.read_csv(summary, sep='\t').iloc[0]
            status['status'] = str(row.status);status['nonempty_gene_slots'] = int(row.nonempty_gene_slots)
        else:
            status['status'] = 'FAILED_REVIEW_REQUIRED_PRESERVE_ORIGINAL_OUTPUTS'
        statuses.append(status)
        (out / 'private_current_fit.json').write_text(json.dumps(status, indent=2) + '\n')
        (out / 'private_all_planned_fit_status.json').write_text(json.dumps(statuses, indent=2) + '\n')
    files = [p for p in out.rglob('*') if p.is_file()];assert sum(p.stat().st_size for p in files) < cap
    terminal = dict(status='ALL_FIVE_PREDECLARED_OFFICIAL_RNA_MCP_FITS_ATTEMPTED_VALIDATION_PENDING',
        completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
        files={str(p.relative_to(out)): dict(sha256=sha(p), bytes=p.stat().st_size) for p in files},
        successful_executions=sum(s['exit_code'] == 0 for s in statuses),
        source_holdout_or_clinical_or_mechanistic_validation=False)
    (out / 'private_terminal.json').write_text(json.dumps(terminal, indent=2) + '\n')
    print(json.dumps(dict(status=terminal['status'], completed_UTC=terminal['completed_UTC'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
