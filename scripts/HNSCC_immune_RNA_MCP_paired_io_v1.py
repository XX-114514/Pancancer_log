"""Immutable all12 paired inputs and source-bound dependency waits."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from NSCLC_immune_RNA_MCP_cohort_io_v2 import wait


def read(root):
    source=root/'HNSCC_immune_RNA_MCP_paired_inputs_20261019_v1';ap=root/'HNSCC_immune_RNA_MCP_paired_input_audit_20261019_v1/private_terminal.json'
    a=json.loads(ap.read_text());assert a['status']=='ACCEPTED_ALL12_PAIRED_INPUTS_ALL312_NODE_PREDICTIONS_AND_ALL_COMPLETE_DESIGNS' and a['source_sha256']=='4cdf9265b9066127e981977837fab525e499dd7c67627b9104ed9501562ced12'
    assert sha(source/'private_terminal.json')==a['producer_terminal_sha256']
    t=json.loads((source/'private_terminal.json').read_text())
    for n,v in t['files'].items():assert sha(source/n)==v['sha256']
    planpath=root/'HNSCC_immune_RNA_MCP_paired_dynamics_contract_20261019_v1/plan_before_any_paired_difference_or_effect.json'
    assert sha(planpath)=='8c03ca8786dcf0de2a872e45cbd527b4eecb3181cc39c71297946a024dd91e0a';plan=json.loads(planpath.read_text())
    assert sha(Path(__file__).with_name('HNSCC_immune_RNA_MCP_paired_math_v1.py'))=='747db0b58146f7662df954f73e6c60f07f02a25b84ede726bd8815232e45d12b'
    assert sha(Path(__file__).with_name('HNSCC_immune_RNA_MCP_paired_audit_math_v1.py'))=='f6647cdb87b0481846455506ab5d8ae643ae19725613f74619b080ddb35d0015'
    scope=pd.read_csv(source/'private_all12_input_design_scope.tsv',sep='\t').sort_values('slot');assert len(scope)==12
    arrays={}
    for r in scope.itertuples(index=False):
        kind='shift' if r.kind=='paired_median_shift' else 'cochange'
        with np.load(source/('private_slot%02d_%s_inputs.npz'%(r.slot,kind)),allow_pickle=False) as z:arrays[r.slot]={k:z[k] for k in z.files}
    return source,ap,plan,scope,arrays
