"""Prespecified nonsequence annotation syntax, not TCR productivity."""
import numpy as np
import pandas as pd

KEYS = ['sampleID', 'cellID']
GENES = ['TRA_v_gene', 'TRA_j_gene', 'TRA_c_gene', 'TRB_v_gene', 'TRB_j_gene', 'TRB_c_gene']
SELECTED = KEYS + ['sub_cell_type'] + GENES
PREFIX = ['TRAV', 'TRAJ', 'TRAC', 'TRBV', 'TRBJ', 'TRBC']
MISSING = {'', 'NA', 'NAN', 'NONE', 'NULL', '<NA>', 'NO_DATA'}


def classify(frame):
    valid = np.zeros(len(frame), np.uint8)
    unknown = valid.copy()
    for j, (field, prefix) in enumerate(zip(GENES, PREFIX)):
        text = frame[field].str.strip().str.upper()
        absent = text.isin(MISSING).to_numpy()
        present = text.str.startswith(prefix).to_numpy() & ~absent
        valid |= present.astype(np.uint8) << j
        unknown |= (~absent & ~present).astype(np.uint8) << j
    alpha = (valid & 7) == 7
    beta = (valid & 56) == 56
    # Complete-chain and paired-chain flags come from a single original row.
    flags = alpha.astype(np.uint8) | (beta.astype(np.uint8) << 1) | ((alpha & beta).astype(np.uint8) << 2)
    flags |= (((valid & 6) == 6).astype(np.uint8) << 3)
    flags |= (((valid & 48) == 48).astype(np.uint8) << 4)
    return valid, unknown, flags


def fixture():
    rows = [
        ['01', 'BC-1', 'gamma_ref', 'TRAV1', 'TRAJ1', 'TRAC', 'TRBV1', 'TRBJ1', 'TRBC1'],
        ['1', 'BC-1', 'other_ref', '', 'NA', 'None', 'NULL', '<NA>', 'no_data'],
        ['01', 'BC-2', 'gamma_ref', 'TRAV29DV5', '', '', '', '', ''],
        ['01', 'BC-2', 'gamma_ref', '', 'TRAJ2', 'TRAC', '', '', ''],
        ['01', 'BC-3', 'other_ref', 'unresolved', '', '', '', '', ''],
    ]
    frame = pd.DataFrame(rows, columns=SELECTED)
    valid, unknown, flags = classify(frame)
    assert valid.tolist() == [63, 0, 1, 6, 0]
    assert unknown.tolist() == [0, 0, 0, 0, 1]
    assert flags.tolist() == [31, 0, 0, 8, 0]
    assert not (flags[2] | flags[3]) & 1  # Never stitch a V-only row to J/C-only.
    keys = pd.MultiIndex.from_frame(frame[KEYS])
    assert keys[0] != keys[1]
    return dict(status='PASS_LEADING_ZERO_MISSING_FORMAT_AND_NO_CHAIN_STITCHING_FIXTURE')
