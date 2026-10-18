"""Run the fixed read-only official-output reconstruction after all fit attempts."""
import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    cap = 64 * 1024**2
    budget = resource(root, cap)
    model = root / 'immune_RNA_MCP_official_models_20261018_v1'
    terminal = model / 'private_terminal.json'
    assert terminal.is_file(), 'Existing five-fit family still pending; do not restart it.'
    statuses = json.loads((model / 'private_all_planned_fit_status.json').read_text())
    assert len(statuses) == 5 and all('exit_code' in s for s in statuses)
    contract = root / 'immune_RNA_MCP_postfit_contract_20261018_v1/plan_before_learned_outputs_or_projection.json'
    assert sha(contract) == '1172d887cffb7c3b0c28a21583154fbe3718a1f13d79edb85500343ead5e2b09'
    cp = json.loads(contract.read_text())
    assert sha(model / 'plan_before_any_official_MCP_fit.json') == cp['official_fit_plan_sha256']
    inputs = root / 'immune_RNA_MCP_native_expression_inputs_20261018_v1'
    assert sha(root / 'immune_RNA_MCP_native_expression_input_audit_20261018_v1/private_terminal.json') == '58385bca98d751470d7f7026b0e6ff4ca7f1c584dacbfa6b7eb5ac644f7b5d4b'
    for name, record in json.loads((inputs / 'private_terminal.json').read_text())['files'].items():
        assert sha(inputs / name) == record['sha256']
    source = Path(__file__).with_suffix('.R')
    assert sha(source) == 'f7738b571a238b5fa34ce34aa370081def84da159fa70b6e26397a87962df2fd'
    model_files = {str(p.relative_to(model)): sha(p) for p in sorted(model.rglob('*.rds'))}
    output = root / 'immune_RNA_MCP_frozen_projection_20261018_v1'
    assert not output.exists()
    output.mkdir(mode=0o700)
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(),
        source_sha256=sha(__file__), R_reader_sha256=sha(source), budget=budget,
        postfit_contract_sha256=sha(contract), original_model_terminal_sha256=sha(terminal),
        original_model_files_sha256=model_files,
        scope='Every planned model and all33 view/programme slots, including missing/no-gene/reconstruction failure. Explicit official sanitized-name map; no model/feature learning, no clinical outcomes.',
        source_offset='Unknown new-source additive effects remain unestimated. No assigning TNBC to a discovery source; only within-heldout-source comparisons are supported.',
        independent_acceptance='PENDING')
    (output / 'plan_before_official_output_read.json').write_text(json.dumps(plan, indent=2) + '\n')
    rscript = os.environ['PANCANCER_DIALOGUE_RSCRIPT']
    with (output / 'private_R_reader_stdout_stderr.log').open('xb') as log:
        result = subprocess.run([rscript, str(source), str(inputs), str(model),
            str(output / 'private_R_projection')], stdout=log, stderr=subprocess.STDOUT,
            timeout=3600)
    assert all(sha(model / name) == digest for name, digest in model_files.items())
    assert sha(terminal) == plan['original_model_terminal_sha256']
    used = sum(p.stat().st_size for p in output.rglob('*') if p.is_file())
    assert used <= cap
    manifest = {str(p.relative_to(output)): sha(p) for p in sorted(output.rglob('*')) if p.is_file()}
    record = dict(status='READER_COMPLETE_INDEPENDENT_ACCEPTANCE_PENDING' if result.returncode == 0 else 'READER_FAILED_PRESERVED_NO_MODEL_RESTART',
        completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
        exit_code=result.returncode, output_bytes=used, files=manifest,
        biological_clinical_or_source_association_validation=False)
    (output / 'private_terminal.json').write_text(json.dumps(record, indent=2) + '\n')
    print(record['status'])
    assert result.returncode == 0


if __name__ == '__main__':
    main()
