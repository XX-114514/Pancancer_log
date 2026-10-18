"""All6 cochange fixed-design finite-null gate, never true RNA effects."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from HNSCC_author_state_clinical_kernel_v2 import fit
from HNSCC_immune_RNA_MCP_paired_math_v1 import design,wild_values
from HNSCC_immune_RNA_MCP_paired_io_v1 import wait,read


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,8*1024**2)
    out=root/'HNSCC_immune_RNA_MCP_paired_calibration_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_fixed_design_null_calibration.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,scope='All6 cochange designs x3fixed error models x256 synthetic data,999 wild each; no true association/paired shift estimated.',seed=2026101904,gate='All available scenario rejection rates at.05<=.10; any failure stops real effects. Finite working-model alarm, not universal typeI guarantee.'),indent=2)+'\n')
    audit=root/'HNSCC_immune_RNA_MCP_paired_input_audit_20261019_v1/private_terminal.json';wait(root,out,'audit_HNSCC_immune_RNA_MCP_paired_inputs_v1',audit,{'ACCEPTED_ALL12_PAIRED_INPUTS_ALL312_NODE_PREDICTIONS_AND_ALL_COMPLETE_DESIGNS'})
    assert sha(root/'HNSCC_immune_RNA_MCP_paired_math_gate_20261019_v1/private_terminal.json')=='0c8a80c69549ea41806f4323913baa05135dd908d85b2bfb5ba2b2773a61aa4c'
    _,_,_,scope,arrays=read(root);rng=np.random.default_rng(2026101904);rows=[];summary=[]
    for item in scope[scope.kind.eq('paired_cochange')].itertuples(index=False):
        a=arrays[item.slot]
        if item.status!='ESTIMABLE':
            for e in ['normal','t3_unit_variance','first_quality_heteroskedastic_normal']:summary.append(dict(slot=item.slot,error_model=e,status='ORIGINAL_DESIGN_UNAVAILABLE',simulations=0,rejections=0,rejection_rate=np.nan))
            continue
        m,status=design(a['x'],a['y'],a['C'],a['names']);assert status=='ESTIMABLE';D=m['D'];n=len(D);mu=D[:,:-1]@np.linspace(-.2,.2,D.shape[1]-1)
        q=a['C'][:,2];q=(q-q.mean())/q.std(ddof=0)
        for ei,e in enumerate(['normal','t3_unit_variance','first_quality_heteroskedastic_normal']):
            reject=0
            for sim in range(256):
                noise=rng.standard_t(3,n)/np.sqrt(3) if e=='t3_unit_variance' else rng.normal(size=n)
                if e=='first_quality_heteroskedastic_normal':noise*=np.exp(.6*q)
                g=fit(mu+noise,D);s=2026101904+item.slot*100000+ei*10000+sim;values=wild_values(g,999,s);count=int((np.abs(values)>=abs(g['t'])-1e-12).sum());reject+=(1+count)/1000<.05
                rows.append(dict(slot=item.slot,error_model=e,simulation=sim,wild_seed=s,beta=g['beta'],HC3_t=g['t'],wild_exceedances=count))
            summary.append(dict(slot=item.slot,error_model=e,status='CALIBRATED',simulations=256,rejections=int(reject),rejection_rate=float(reject/256)))
            (out/'private_progress.json').write_text(json.dumps(dict(phase='SYNTHETIC_NULL_ONLY',slot=item.slot,error_model=e,simulations=256,observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n')
    pd.DataFrame(rows).to_csv(out/'private_all_null_statistics.tsv',sep='\t',index=False);table=pd.DataFrame(summary);assert len(table)==18;table.to_csv(out/'private_all18_null_scenarios.tsv',sep='\t',index=False)
    tested=table[table.status.eq('CALIBRATED')];passed=len(tested)>0 and tested.rejection_rate.le(.10).all()
    terminal=dict(status='PASS_ALL_PAIRED_FIXED_DESIGN_NULL_DIAGNOSTICS' if passed else 'FAIL_NULL_DIAGNOSTIC_NO_TRUE_EFFECTS',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),input_acceptance_sha256=sha(audit),simulated_data_sets=len(rows),all_scenarios=18,actual_effects_computed=False,files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(out.iterdir()) if p.is_file()})
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(dict(status=terminal['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
