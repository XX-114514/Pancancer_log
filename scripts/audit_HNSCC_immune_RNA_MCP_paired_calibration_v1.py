"""All synthetic null data/statistics/exceedances independently replayed via QR."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from HNSCC_immune_RNA_MCP_paired_math_v1 import design
from HNSCC_immune_RNA_MCP_paired_audit_math_v1 import reference_model,reference_wild
from HNSCC_immune_RNA_MCP_paired_io_v1 import wait,read


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,4*1024**2)
    source=root/'HNSCC_immune_RNA_MCP_paired_calibration_20261019_v1';out=root/'HNSCC_immune_RNA_MCP_paired_calibration_audit_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_complete_null_replay.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,scope='All18 null scenarios, every synthetic HC3 beta/t and999 wild values regenerated for each data set; no actual paired/source effects.'),indent=2)+'\n')
    terminal=wait(root,out,'calibrate_HNSCC_immune_RNA_MCP_paired_dynamics_v1',source/'private_terminal.json',{'PASS_ALL_PAIRED_FIXED_DESIGN_NULL_DIAGNOSTICS','FAIL_NULL_DIAGNOSTIC_NO_TRUE_EFFECTS'})
    for n,v in terminal['files'].items():assert sha(source/n)==v['sha256']
    _,_,_,scope,arrays=read(root);stats=pd.read_csv(source/'private_all_null_statistics.tsv',sep='\t').set_index(['slot','error_model','simulation']);summary=pd.read_csv(source/'private_all18_null_scenarios.tsv',sep='\t').set_index(['slot','error_model']);assert len(summary)==18
    rng=np.random.default_rng(2026101904);checked=0;maxerror=0.;passed=True
    for item in scope[scope.kind.eq('paired_cochange')].itertuples(index=False):
        if item.status!='ESTIMABLE':
            for e in ['normal','t3_unit_variance','first_quality_heteroskedastic_normal']:assert summary.loc[(item.slot,e),'status']=='ORIGINAL_DESIGN_UNAVAILABLE'
            continue
        a=arrays[item.slot];m,status=design(a['x'],a['y'],a['C'],a['names']);ref,rs=reference_model(a['x'],a['y'],a['C'],a['names']);assert status==rs=='ESTIMABLE'
        # Recreate the frozen SVD synthetic mean only; estimation is independent QR.
        D=m['D'];n=len(D);mu=D[:,:-1]@np.linspace(-.2,.2,D.shape[1]-1);q=a['C'][:,2];q=(q-q.mean())/q.std(ddof=0)
        for ei,e in enumerate(['normal','t3_unit_variance','first_quality_heteroskedastic_normal']):
            reject=0
            for sim in range(256):
                noise=rng.standard_t(3,n)/np.sqrt(3) if e=='t3_unit_variance' else rng.normal(size=n)
                if e=='first_quality_heteroskedastic_normal':noise*=np.exp(.6*q)
                y=mu+noise;s=2026101904+item.slot*100000+ei*10000+sim;g=sm.OLS(y,ref['D']).fit(cov_type='HC3');values=reference_wild(y,ref['D'],999,s);count=int((np.abs(values)>=abs(g.tvalues[-1])-1e-12).sum());r=stats.loc[(item.slot,e,sim)]
                err=max(abs(r.beta-g.params[-1]),abs(r.HC3_t-g.tvalues[-1]));assert err<1e-8 and r.wild_seed==s and r.wild_exceedances==count;maxerror=max(maxerror,float(err));reject+=(1+count)/1000<.05;checked+=1
            r=summary.loc[(item.slot,e)];assert r.simulations==256 and r.rejections==reject and abs(r.rejection_rate-reject/256)<1e-12;passed=passed and reject/256<=.10
            (out/'private_progress.json').write_text(json.dumps(dict(phase='INDEPENDENT_NULL_ONLY_REPLAY',slot=item.slot,error_model=e,synthetic_data_checked=checked,observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n')
    assert checked==len(stats)==terminal['simulated_data_sets'];assert (terminal['status']=='PASS_ALL_PAIRED_FIXED_DESIGN_NULL_DIAGNOSTICS')==passed
    status='ACCEPTED_ALL_PAIRED_NULL_SIMULATIONS_AND_PASS_GATE' if passed else 'ACCEPTED_NULL_GATE_FAILURE_NO_TRUE_EFFECTS'
    t=dict(status=status,completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),producer_terminal_sha256=sha(source/'private_terminal.json'),synthetic_data_sets_checked=checked,complete_wild_statistics_regenerated=checked*999,maximum_statistic_error=maxerror,actual_effects_computed=False)
    (out/'private_terminal.json').write_text(json.dumps(t,indent=2)+'\n');print(json.dumps(dict(status=status,terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
