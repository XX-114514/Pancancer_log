"""Immutable new-family inputs and unique outcome-free design mapping."""
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from NSCLC_immune_RNA_MCP_cohort_io_v2 import wait


def read(root):
    directory=root/'immune_RNA_MCP_new_source_family_inputs_20261019_v2';audit=root/'immune_RNA_MCP_new_source_family_input_audit_20261019_v1/private_terminal.json'
    assert sha(audit)=='c6edcbc68aed60916d9dcfc12cae0b81b97a44a23109a36bf0520411e376b087';assert sha(directory/'private_terminal.json')==json.loads(audit.read_text())['producer_terminal_sha256']
    term=json.loads((directory/'private_terminal.json').read_text())
    for name,value in term['files'].items():assert sha(directory/name)==value['sha256']
    plan=json.loads((directory/'plan_before_any_new_source_or_clinical_effect.json').read_text());assert plan['math_sha256']==sha(Path(__file__).with_name('immune_RNA_MCP_new_source_math_v1.py'))
    registry=pd.read_csv(directory/'private_all_registered_models_and_pre_effect_eligibility.tsv',sep='\t',keep_default_na=False)
    arrays={}
    for row in registry.itertuples(index=False):
        with np.load(directory/row.array_file,allow_pickle=False) as a:arrays[row.model_id]={key:a[key] for key in a.files}
    return directory,audit,plan,registry,arrays


def unique_designs(registry,arrays):
    unique={};mapping={}
    for row in registry.itertuples(index=False):
        if row.family=='source_marginal' or row.status!='ESTIMABLE':continue
        a=arrays[row.model_id];h=hashlib.sha256(row.kind.encode())
        for field in ['cases','x','C','names']:
            h.update(field.encode());h.update(str(a[field].dtype).encode());h.update(str(a[field].shape).encode());h.update(a[field].tobytes())
        key=h.hexdigest()
        if key not in unique:unique[key]=row
        mapping[row.model_id]=key
    return list(unique.items()),mapping
