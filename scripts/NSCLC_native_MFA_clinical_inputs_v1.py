"""Frozen input and storage utilities; no clinical effect estimation."""
import json
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha, VIEWS

COMMON=['age_numeric','male','smoker','LUSC','stage_ordinal','cycles_numeric','chemo_binary','targeted_binary']


def matched_quota(root):
    df=subprocess.run(['df','-P',str(root.resolve())],capture_output=True,text=True,timeout=20)
    quota=subprocess.run(['quota','-g','-v','-w','USER002'],capture_output=True,text=True,timeout=20)
    assert df.returncode==quota.returncode==0
    filesystem=df.stdout.splitlines()[-1].split()[0]
    rows=[s.split() for s in quota.stdout.splitlines() if s.split() and s.split()[0]==filesystem]
    assert len(rows)==1
    remaining=(int(rows[0][3])-int(rows[0][1].rstrip('*')))*1024
    assert min(remaining,shutil.disk_usage(root).free)>500_000_000_000+1073741824+268435456
    return dict(group_remaining_bytes=remaining,reserve_bytes=500_000_000_000,buffer_bytes=1073741824,cap_bytes=268435456)


def read():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);pre=root/'NSCLC_native_MFA_clinical_inputs_20261016_v1'
    terminal=json.loads((pre/'private_terminal.json').read_text())
    assert terminal['status'].startswith('COMPLETE_')
    for name,record in terminal['files'].items():assert sha(pre/name)==record['sha256']
    plan=json.loads((pre/'plan_before_any_native_factor_clinical_effect.json').read_text())
    kernel=Path(__file__).with_name('HNSCC_author_state_clinical_kernel_v2.py')
    assert sha(kernel)==plan['kernel_sha256']
    scope=pd.read_csv(pre/'private_all40_primary_models.tsv',sep='\t')
    pairs=pd.read_csv(pre/'private_all60_primary_pairs.tsv',sep='\t')
    assert len(scope)==40 and len(pairs)==60 and scope.record.is_unique and pairs.slot.is_unique
    arrays={}
    for partition in ['discovery','internal_holdout','all_frozen_four_view_cases']:
        for vi,view in enumerate(VIEWS):
            with np.load(pre/('private_'+partition+'_'+view+'.npz'),allow_pickle=False) as z:
                arrays[(partition,view)]={k:z[k] for k in z.files}
    return root,pre,plan,scope,pairs,arrays
