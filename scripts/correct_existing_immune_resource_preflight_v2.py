"""Match group quota to the actual output filesystem; retain original scope."""
import hashlib
import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    scope = root / 'existing_immune_fine_native_measurement_scope_20261016_v1'
    old = json.loads((scope / 'private_terminal.json').read_text())
    for name, digest in old['files'].items():
        assert sha(scope / name) == digest
    original = json.loads((scope / 'plan_and_resource_preflight.json').read_text())
    out = root / 'existing_immune_native_resource_preflight_20261016_v2'
    assert not out.exists()
    out.mkdir(mode=0o700)
    df = subprocess.run(['df', '-P', str(root.resolve())], capture_output=True, text=True, timeout=20)
    quota = subprocess.run(['quota', '-g', '-v', '-w', 'USER002'], capture_output=True, text=True, timeout=20)
    assert df.returncode == quota.returncode == 0
    (out / 'private_output_filesystem_df.txt').write_text(df.stdout)
    (out / 'private_exact_group_quota.txt').write_text(quota.stdout)
    filesystem = df.stdout.splitlines()[-1].split()[0]
    rows = [line.split() for line in quota.stdout.splitlines() if line.split() and line.split()[0] == filesystem]
    assert len(rows) == 1, 'Actual output filesystem quota not uniquely resolved'
    row = rows[0]
    used, limit = int(row[1].rstrip('*')), int(row[3])
    assert limit > 0 and used < limit
    remaining = (limit - used) * 1024
    RAM = original['available_ram_bytes'] > original['dense_workspace_bytes'] * 2.5 + 4*1024**3
    needed = original['quota_reserve_decimal_bytes'] + original['future_output_cap_bytes'] + original['quota_extra_headroom_bytes']
    plan = dict(observed_UTC=datetime.now(timezone.utc).isoformat(), code_sha256=sha(__file__), original_scope_terminal_sha256=sha(scope / 'private_terminal.json'), original_resource_plan_sha256=sha(scope / 'plan_and_resource_preflight.json'),
        correction='V1 regex selected the first quota filesystem, which can have zero/unlimited hard limit. V2 matches df of the real output path. Original scope, RAM estimate, output cap, reserve and headroom unchanged.',
        exact_group_quota_remaining_bytes=remaining, exact_required_remaining_bytes=needed, group_quota_gate=remaining > needed, RAM_gate=RAM, resource_gate=RAM and remaining > needed,
        nodes=original['nodes'], native_units=original['native_units'], selected_original_cells=original['selected_original_cells'], no_scope_or_labels_or_expression_changed=True)
    p = out / 'corrected_resource_preflight.json'
    p.write_text(json.dumps(plan, indent=2) + '\n')
    result = dict(status='ACTUAL_OUTPUT_FILESYSTEM_QUOTA_RECONCILED_ORIGINAL_FINE_SCOPE_AND_RESERVE_UNCHANGED', finished_UTC=datetime.now(timezone.utc).isoformat(), producer_sha256=sha(__file__), resource_gate=plan['resource_gate'], files={p.name: sha(p) for p in out.iterdir() if p.is_file()}, no_expression_or_models_computed=True)
    (out / 'private_terminal.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
