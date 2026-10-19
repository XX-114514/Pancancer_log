"""All-source label and endpoint recount using independent text/token arithmetic."""
import json
import os
import re
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC_external_lineage_inputs_v2 import load, cell_arms
from NSCLC_external_lineage_role_rules_v1 import author_role

CHANNELS = ['author_exact_VCAN_mac', 'author_gamma_with_mac',
            'author_gamma_with_myeloid_parent', 'author_marker_with_myeloid_parent',
            'frozen_marker_grouping']


def classify(value):
    text = str(value).casefold()
    for letter, word in [('γ', 'gamma '), ('δ', 'delta '), ('α', 'alpha '), ('β', 'beta ')]:
        text = text.replace(letter, word)
    tokens = re.findall(r'[^\W_]+', text)
    compact = ''.join(tokens)
    explicit = 'gammadelta' in compact or any(v in {'gd', 'gdt', 'tgd', 'tcrgd', 'tcrgdt'} for v in tokens)
    # Preserve exact complete-word semantics for separated g/d words.
    original = re.sub(r'[_\-]+', ' ', text)
    joined_gamma = re.search(r'gamma\s*delta', original) is not None
    explicit = joined_gamma or any(v in {'gd', 'gdt', 'tgd', 'tcrgd', 'tcrgdt'} for v in tokens)
    mixed = explicit and (any(v in {'nk', 'nkt'} for v in tokens) or re.search(r'alpha\s*beta', original) is not None)
    marker = any(re.fullmatch(r'(?:trdc|trdv|trg|trgc)[0-9]*', v) for v in tokens)
    marker = marker or any(re.fullmatch(r'vdelta[0-9]+', ''.join(tokens[i:i + width])) for i in range(len(tokens)) for width in [1, 2, 3])
    return 2 if mixed else 1 if explicit else 3 if marker else 0


def VCAN(value):
    return 'vcan' in re.findall(r'[^\W_]+', str(value).casefold())


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root, 2 * 1024**2)
    out = root / 'NSCLC_external_endpoint_metadata_audit_20261019_v1'
    assert not out.exists(); out.mkdir(mode=0o700)
    src = root / 'NSCLC_external_endpoint_metadata_qualification_20261019_v1'
    terminal = json.loads((src / 'private_terminal.json').read_text())
    for name, d in terminal['files'].items(): assert sha(src / name) == d['sha256']
    data = load(); lib = data['lib']
    assert data['hashes'] == terminal['source_metadata_hashes']
    expected = np.zeros((len(lib), 5, 2), np.int64)
    diag = np.zeros((len(lib), 4), np.int64)
    rows = pd.read_csv(src / 'private_all_author_label_rule_classifications.tsv.gz', sep='\t', keep_default_na=False).set_index(['arm', 'author_label'])
    oldproofs = json.loads((src / 'private_source_metadata_read_proofs.json').read_text())
    for armnumber, chunk in enumerate(cell_arms(data)):
        assert oldproofs[armnumber]['author_values_sha256'] == chunk['author_values_sha256']
        assert oldproofs[armnumber]['cells'] == chunk['cell_count']
        labels = list(chunk['author'].categories)
        categories = np.array([classify(v) for v in labels] + [0], np.int8)
        vc = np.array([VCAN(v) for v in labels] + [False], bool)
        macro = np.array([author_role(v) == 2 for v in labels] + [False], bool)
        for i, label in enumerate(labels):
            row = rows.loc[(chunk['arm'], str(label))]
            assert row.gamma_text_class == categories[i]
            assert row.VCAN_named == vc[i]
            assert row.unambiguous_macrophage_text == macro[i]
        codes = chunk['author'].codes
        c = categories[codes]; v = vc[codes]; mac = macro[codes]
        node = chunk['node']; majors = chunk['major']
        T = np.array([m in {'T', 'T cell', 'T/NK lymphoid', 'Lymphoid cell'} for m in majors])[node.codes]
        M = np.array([m in {'Macrophage', 'Myeloid'} for m in majors])[node.codes]
        author_gamma = (c == 1) & T
        author_marker = np.isin(c, [1, 3]) & T
        raw_mac = mac & M
        catalog = data['catalogs'][chunk['arm']].set_index('node_id')
        frozen = catalog.loc[np.asarray(node.categories).astype(str), 'label'].astype(str)
        ng = np.array([classify(s) in {1, 3} for s in frozen])[node.codes] & T
        nv = np.array([VCAN(s) for s in frozen])[node.codes] & M
        pairs = [(raw_mac & v, author_gamma), (raw_mac, author_gamma), (M, author_gamma), (M, author_marker), (nv, ng)]
        frame = pd.DataFrame({'library_position': chunk['library_index']})
        for j, pair in enumerate(pairs):
            for e, value in enumerate(pair): frame['e%d_%d' % (j, e)] = value
        flags = [codes < 0, c == 2, (c == 1) & ~T, (c == 3) & ~T]
        for j, flag in enumerate(flags): frame['d%d' % j] = flag
        summed = frame.groupby('library_position', sort=False).sum().reindex(range(len(lib)), fill_value=0)
        for j in range(5):
            for e in range(2): expected[:, j, e] += summed['e%d_%d' % (j, e)].to_numpy(np.int64)
        for j in range(4): diag[:, j] += summed['d%d' % j].to_numpy(np.int64)
        (out / 'private_progress.json').write_text(json.dumps(dict(stage='INDEPENDENT_ALL_CELL_AUTHOR_METADATA_RECOUNT', completed_arm=chunk['arm'], observed_UTC=datetime.now(timezone.utc).isoformat())) + '\n')
    observed = pd.read_csv(src / 'private_all_library_endpoint_counts.tsv.gz', sep='\t', keep_default_na=False)
    assert len(observed) == 5 * len(lib)
    for j, name in enumerate(CHANNELS):
        f = observed[observed.channel.eq(name)].reset_index(drop=True)
        assert f[['arm', 'dataset', 'library_id']].equals(lib[['arm', 'dataset', 'library_id']])
        np.testing.assert_array_equal(f[['sender_cells', 'receiver_cells']], expected[:, j])
        np.testing.assert_array_equal(f.both_endpoints_atleast20, (expected[:, j] >= 20).all(1))
    diagrows = pd.read_csv(src / 'private_all_library_label_diagnostics.tsv.gz', sep='\t', keep_default_na=False)
    np.testing.assert_array_equal(diagrows[['missing_author_label', 'mixed_gamma_label', 'explicit_gamma_parent_discordant', 'marker_gamma_parent_discordant']], diag)
    for filename, columns in [('private_complete_source_endpoint_summary.tsv', ['arm', 'dataset', 'channel']), ('private_complete_source_stratified_endpoint_summary.tsv', ['arm', 'dataset', 'channel', 'cancer_label', 'phase_class', 'material_class', 'therapy_keyword_flag'])]:
        recorded = pd.read_csv(src / filename, sep='\t', keep_default_na=False).set_index(columns)
        groups = observed.groupby(columns, sort=True, dropna=False)
        assert len(recorded) == len(groups)
        for key, group in groups:
            row = recorded.loc[key]
            available = group.both_endpoints_atleast20
            reviewed = group.reviewed_patient.astype(bool)
            integers = dict(libraries=len(group), libraries_receiver_atleast20=int(group.receiver_cells.ge(20).sum()), libraries_both_endpoints_atleast20=int(available.sum()), sender_cells=int(group.sender_cells.sum()), receiver_cells=int(group.receiver_cells.sum()), reviewed_patient_blocks=group.loc[reviewed, 'patient_block'].nunique(), reviewed_patient_blocks_both_endpoints=group.loc[reviewed & available, 'patient_block'].nunique(), unknown_identity_eligible_libraries=int((~reviewed & available).sum()))
            for column, value in integers.items(): assert row[column] == value
            assert not row.any_response_field_in_current_qualification_table
    receipt = dict(status='ACCEPTED_ALL_SOURCE_TEXT_CLASSES_NATIVE_METADATA_COUNTS_ALL_LIBRARIES_SOURCE_AND_STRATIFIED_SUMMARIES', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), producer_terminal_sha256=sha(src / 'private_terminal.json'), source_metadata_hashes=data['hashes'], all_author_categories_checked=len(rows), libraries_checked=len(lib), library_channel_positions_checked=len(observed), expression_matrix_read=False, new_clinical_adoption=False, P_values_computed=False, resource_guard=guard, scientific_outputs_public=False)
    (out / 'private_terminal.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
