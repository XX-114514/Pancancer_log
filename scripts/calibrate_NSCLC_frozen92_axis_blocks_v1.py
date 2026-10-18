"""Three preselected real designs, synthetic null outcomes only, before effects."""
import json
import os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC_frozen92_axis_io_v2 import prepare
from NSCLC_frozen_axis_block_math_v2 import evaluate,wild,confidence,seed
from NSCLC_frozen_axis_block_audit_math_v1 import prepare_axis,evaluate as independent_evaluate,wild as independent_wild


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,2*1024**2);out=root/'NSCLC_frozen92_axis_block_null_calibration_20261019_v1';assert not out.exists();out.mkdir(mode=0o700);src=root/'NSCLC_frozen92_axis_source_inputs_20261019_v2';audit=root/'NSCLC_frozen92_axis_source_input_audit_20261019_v1/private_terminal.json';t=json.loads(audit.read_text());assert t['status'].startswith('ACCEPTED_ALL92_NATIVE_SCORES') and t['producer_terminal_sha256']==sha(src/'private_terminal.json')
    save=lambda name,data:(out/name).write_text(json.dumps(data,indent=2)+'\n')
    a=dict(np.load(src/'private_all129_all92_native_axis_and_design_inputs.npz'));coverage=pd.read_csv(src/'private_all92_source_design_coverage.tsv',sep='\t');eligible=coverage[coverage.status.eq('ESTIMABLE')].sort_values(['author_records','axis_id'],kind='stable');chosen=eligible.iloc[[0,len(eligible)//2,len(eligible)-1]].axis_id.tolist() if len(eligible) else []
    save('plan_before_null_simulations.json',dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),input_audit_sha256=sha(audit),selected_designs=chosen,selection='first median last by original eligibleN then frozenaxisid, never by true effect',repetitions=200,each_stream_draws=999,null_noise='sigma=(.5+.35abs(zage)+.25response); sigma*(.6profileblockGaussian+.8individualGaussian). Nuisance mean fixed in full nuisance span, true response coefficientzero.',gates='Every design rejectrate<=.10 and fullmodel95CIcoverage>=.90. Stop before trueeffects if fails, no seed or gate retries. Finite diagnostic not exact coverage guarantee.',resource_guard=guard,true_clinical_effects_read=False))
    rows=[];maximum=0.
    for selection,axis in enumerate(chosen):
        j=a['axis_ids'].tolist().index(axis);model,g,status,valid=prepare(a,j);_,independent,why,_=prepare_axis(a,j);assert g is not None and independent is not None;age=model['controls'][:,0];z=(age-age.mean())/(age.std() or 1.);sigma=.5+.35*np.abs(z)+.25*model['x'];mu=model['D'][:,:-1]@np.linspace(.1,.3,model['nuisance_rank']);rng=np.random.default_rng(seed(axis,'calibration_outcomes_'+str(selection)));P=[];cover=[]
        for repeat in range(200):
            y=mu+sigma*(.6*rng.normal(size=g['groups'])[g['membership']]+.8*rng.normal(size=len(mu)));fit=evaluate(g,y);restricted=wild(g,y,999,seed(axis,'calibration_null_'+str(selection)+'_'+str(repeat)),True);full=wild(g,y,999,seed(axis,'calibration_full_'+str(selection)+'_'+str(repeat)),False);p=(1+int((np.abs(restricted)>=abs(fit['t'])-1e-12).sum()))/1000;interval=confidence(fit,full);covered=bool(interval[0]<=0<=interval[1]);P.append(p);cover.append(covered);rows.append(dict(selection=selection,axis_id=axis,repeat=repeat,null_P=p,interval_covers_zero=covered))
            other=independent_evaluate(independent,y);nr=independent_wild(independent,y,999,seed(axis,'calibration_null_'+str(selection)+'_'+str(repeat)),True);nf=independent_wild(independent,y,999,seed(axis,'calibration_full_'+str(selection)+'_'+str(repeat)),False);np.testing.assert_allclose([fit['effect'],fit['se'],fit['t']],[other['effect'],other['se'],other['t']],atol=1e-8,rtol=1e-8);np.testing.assert_allclose(restricted,nr,atol=1e-8,rtol=1e-8);np.testing.assert_allclose(full,nf,atol=1e-8,rtol=1e-8);maximum=max(maximum,float(np.max(np.abs(restricted-nr))),float(np.max(np.abs(full-nf))))
        save('private_progress.json',dict(phase='FINITE_SYNTHETIC_NULL_OUTCOMES_ONLY',completed_design=selection+1,total_designs=len(chosen),observed_UTC=datetime.now(timezone.utc).isoformat()))
    frame=pd.DataFrame(rows,columns=['selection','axis_id','repeat','null_P','interval_covers_zero']);frame.to_csv(out/'private_all600_synthetic_null_experiments.tsv',sep='\t',index=False);summary=[]
    for (selection,axis),f in frame.groupby(['selection','axis_id'],sort=False):summary.append(dict(selection=int(selection),axis_id=axis,repetitions=len(f),rejection_rate=float((f.null_P<=.05).mean()),coverage=float(f.interval_covers_zero.mean())))
    good=all(x['rejection_rate']<=.10 and x['coverage']>=.90 for x in summary);save('private_all3_null_design_summaries.json',summary)
    receipt=dict(status='ACCEPTED_PREDECLARED_NULL_DIAGNOSTICS_BEFORE_ALL92_TRUE_EFFECTS' if good else 'FAILED_PREDECLARED_NULL_DIAGNOSTICS_TRUE_EFFECTS_MUST_NOT_RUN',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),input_audit_sha256=sha(audit),selected_designs=chosen,synthetic_experiments=len(rows),true_clinical_effects_read=False,files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()},scientific_outputs_public=False)
    receipt.update(all_synthetic_random_statistics_independently_replayed=True,maximum_independent_statistic_difference=maximum)
    assert sum(v['bytes'] for v in receipt['files'].values())+sum(p.stat().st_size for p in src.iterdir() if p.is_file())<2*1024**2;save('private_terminal.json',receipt);print(json.dumps(dict(status=receipt['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
