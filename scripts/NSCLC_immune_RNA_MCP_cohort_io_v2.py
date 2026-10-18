"""Source-bound waits and immutable inputs for the fixed RNA cohort family."""
import json
import time
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha


def wait(root,out,script,terminal,allowed):
    launchpath=root/'HNSCC_receiver_followup_job_control_20261016_v2'/script/'launch_receipt.json'
    launch=json.loads(launchpath.read_text());assert launch['source_sha256']==sha(Path(__file__).with_name(script+'.py'))
    deadline=time.monotonic()+21600
    while not terminal.exists():
        proc=Path('/proc')/str(launch['pid'])/'stat';assert proc.exists(),'ORIGINAL_DEPENDENCY_MISSING_NO_RESTART'
        s=proc.read_text().rsplit(')',1)[1].split();assert s[19]==str(launch['proc_start_ticks']) and s[0]!='Z';assert time.monotonic()<deadline
        (out/'private_progress.json').write_text(json.dumps(dict(stage='WAITING_SAME_LIVE_DEPENDENCY',script=script,pid=launch['pid'],startticks=s[19],observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n');time.sleep(15)
    t=json.loads(terminal.read_text());assert t['status'] in allowed and t['source_sha256']==launch['source_sha256']
    return t


def read(root):
    source=root/'NSCLC_immune_RNA_MCP_cohort_inputs_20261019_v1';audit=root/'NSCLC_immune_RNA_MCP_cohort_inputs_audit_20261019_v2/private_terminal.json'
    assert json.loads(audit.read_text())['status']=='ACCEPTED_ALL6_BUNDLES_18_INPUT_DESIGNS_AND_ALL_EXCLUSIONS'
    t=json.loads((source/'private_terminal.json').read_text())
    assert json.loads(audit.read_text())['producer_terminal_sha256']==sha(source/'private_terminal.json')
    for name,item in t['files'].items():assert sha(source/name)==item['sha256']
    plan=json.loads((source/'plan_before_any_source_coordination_or_clinical_effect.json').read_text())
    assert plan['math_sha256']==sha(Path(__file__).with_name('NSCLC_immune_RNA_MCP_cohort_math_v1.py'))
    scope=pd.read_csv(source/'private_all18_design_scope.tsv',sep='\t');assert len(scope)==18
    arrays={}
    for row in scope.itertuples(index=False):
        path=source/('private_slot%02d_%s_inputs.npz'%(row.slot,row.component))
        with np.load(path,allow_pickle=False) as z:arrays[(row.slot,row.component)]={k:z[k] for k in z.files}
    return source,audit,plan,scope,arrays
