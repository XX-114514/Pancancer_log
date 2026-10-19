"""Fresh metadata axes and public stable-ID authority; no expression reads."""
import hashlib
import json
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
from anndata._io.specs import read_elem


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    authority = json.load(sys.stdin)
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    source = Path(os.environ['PANCANCER_INTERFACE_PRODUCTION'])
    out = Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert out.parent.parent == Path('/dev/shm') and not out.exists()
    term_path = source / 'private_terminal.json'
    term = json.loads(term_path.read_text())
    assert term['status'] == 'COMPLETE_ALL_FIXED_MEASUREMENT_ROLES_AND_TWO_UNESTIMABLE_QUESTIONS_NO_EFFECTS'
    for name, proof in term['files'].items():
        p = source / name
        assert p.stat().st_size == proof['bytes'] and sha(p) == proof['sha256']
    plan = json.loads((source / 'private_plan_before_any_current_source_effect.json').read_text())
    definition = plan['private_interface_definition']
    symbols = definition['sender_components'] + [definition['receiver_component']]
    assert authority['symbols'] == symbols and len(authority['stable_ids']) == 3
    contract_path = root.parents[2] / 'Project_v5/runs/20261013_immunotherapy_atlas_v1/config/source_contract_v2_private.json'
    assert sha(contract_path) == '8d7f54a385389d0ecbc825155df727ed065cae2ed1c3f3ed584c1a7896647e71'
    contract = json.loads(contract_path.read_text())
    rec = next(r for r in contract['records'] if r['dataset'] == 'ICB_Bassez')
    positive = set(np.load(root / 'Bassez_fixed_function_readiness_20261014_v1/private_source_positive_confirmed_columns.npy').tolist())
    production = json.loads((source / 'private_all3_fixed_interface_measurement_qualification.json').read_text())
    rows = []
    axes = []
    with h5py.File(rec['path'], 'r') as h:
        for axis in ['var', 'raw/var']:
            v = read_elem(h[axis])
            keys = [str(k).split('.')[0] for k in v.index]
            fields = [k for k in ['gene_symbol', 'feature_name', 'gene_name_original'] if k in v]
            matrix = [[str(n)] + [str(v.iloc[i][k]) for k in fields] for i, n in enumerate(v.index)]
            for role, (symbol, stable) in enumerate(zip(symbols, authority['stable_ids'])):
                stable_positions = {i for i, key in enumerate(keys) if key == stable}
                exact_positions = {i for i, values in enumerate(matrix) if symbol in values}
                norm = lambda s: re.sub('[^A-Za-z0-9]', '', s).casefold()
                lexical = {i for i, values in enumerate(matrix) if norm(symbol) in {norm(s) for s in values}}
                union = stable_positions | exact_positions
                present = bool(len(union) == 1 and next(iter(union)) in positive)
                row = dict(axis=axis, role=role, stable_ID_physical_columns=len(stable_positions),
                           exact_symbol_columns=len(exact_positions), normalized_lexical_candidates=len(lexical),
                           all_candidate_physical_columns=len(union), source_positive_confirmed=present,
                           lexical_mapping_not_adopted=True)
                rows.append(row)
                if axis == 'raw/var':
                    prior = production['roles'][role]
                    assert len(union) == prior['resolved_physical_columns']
                    assert present == prior['unique_source_positive_confirmed']
                    assert len(exact_positions) == prior['direct_columns']
            axes.append(dict(axis=axis, features=len(v), complete_feature_axis_sha256=hashlib.sha256(
                json.dumps(matrix, ensure_ascii=False, separators=(',', ':')).encode()).hexdigest()))
    assert axes[0]['complete_feature_axis_sha256'] == axes[1]['complete_feature_axis_sha256']
    result = json.loads((source / 'private_complete2_source_association_results.json').read_text())
    assert [r['slot'] for r in result] == ['BIC01', 'BIC02']
    for r in result:
        assert r['status'] == 'UNESTIMABLE_FIXED_INTERFACE_MEASUREMENT'
        assert all(r[k] is None for k in ['effect', 'P', 'BH2', 'BY2'])
        assert r['internal_planned_P_for_correction'] == 1 and r['missing_not_negative'] and r['no_endpoint_replacement']
    old = source.parent.parent / 'pancancer-Bassez-fixed-interface-20261020-v1'
    old_receipt = json.loads((old / 'launch_receipt.json').read_text())
    current_receipt = json.loads((source.parent / 'launch_receipt.json').read_text())
    assert old_receipt['private_config_sha256'] == current_receipt['private_config_sha256']
    failure = (old / 'private_job_stdout_stderr.log').read_text()
    assert 'AssertionError' in re.findall(r'\b([A-Za-z]+Error)(?::|$)', failure, re.M)
    assert not (old / 'private_output/private_terminal.json').exists()
    out.mkdir(mode=0o700)
    detail = dict(observed_UTC=datetime.now(timezone.utc).isoformat(), axes=axes, roles=rows,
                  authority_provenance=authority['provenance'], original_interface_config_unchanged=True,
                  v1_failure_preserved=True, v1_source_sha256=old_receipt['source_sha256'],
                  all_expression_and_score_arrays_not_read=True, missing_not_zero=True,
                  no_clinical_effect_or_identity_claim=True)
    (out / 'private_source_schema_only_measurement_acceptance.json').write_text(json.dumps(detail, indent=2) + '\n')
    terminal = dict(status='ACCEPTED_ALL_FIXED_ROLES_BOTH_AXES_AND_PLANNED_UNESTIMABLE_SLOTS',
                    completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
                    producer_terminal_sha256=sha(term_path), native_axes_checked=2, interface_role_positions=6,
                    planned_question_slots_checked=2, no_expression_or_statistical_inference_executed=True,
                    scientific_source_detail_sha256=sha(out / 'private_source_schema_only_measurement_acceptance.json'))
    (out / 'private_terminal.json').write_text(json.dumps(terminal, indent=2) + '\n')
    print(json.dumps(dict(status=terminal['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
