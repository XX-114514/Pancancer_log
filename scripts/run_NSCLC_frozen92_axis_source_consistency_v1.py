"""All92 predeclared within-cohort clinical slots; no endpoint/seed replacement."""
import json
import os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import t as student_t
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC_immune_RNA_MCP_cohort_math_v1 import adjust
from NSCLC_frozen_axis_block_math_v2 import evaluate,wild,confidence,seed
from NSCLC_frozen92_axis_io_v2 import prepare,SENSITIVITIES


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,32*1024**2);out=root/'NSCLC_frozen92_axis_source_consistency_20261019_v1';assert not out.exists();out.mkdir(mode=0o700);src=root/'NSCLC_frozen92_axis_source_inputs_20261019_v2';audit=root/'NSCLC_frozen92_axis_source_input_audit_20261019_v1/private_terminal.json';cal=root/'NSCLC_frozen92_axis_block_null_calibration_20261019_v1/private_terminal.json';at=json.loads(audit.read_text());ct=json.loads(cal.read_text());assert at['status'].startswith('ACCEPTED_ALL92_NATIVE_SCORES') and at['producer_terminal_sha256']==sha(src/'private_terminal.json');assert ct['status']=='ACCEPTED_PREDECLARED_NULL_DIAGNOSTICS_BEFORE_ALL92_TRUE_EFFECTS' and ct['input_audit_sha256']==sha(audit)
    st=json.loads((src/'private_terminal.json').read_text())
    for name,d in st['files'].items():assert sha(src/name)==d['sha256']
    save=lambda n,d:(out/n).write_text(json.dumps(d,indent=2)+'\n')
    save('plan_before_all92_true_axis_effects.json',dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),input_contract_sha256=sha(src/'plan_before_all92_source_axis_readouts_and_effects.json'),input_audit_sha256=sha(audit),null_diagnostic_sha256=sha(cal),resource_guard=guard,scope='Complete92 predeclared slots, original129 separate label records; main19999WCRB,19999WCUB,primaryBH/BY92includingfailures. All276fixedsensitivityeffects and184histologydescriptors. No effects were used to select current axis list or change gates.',fixed_design_intervals_not_external_validation=True))
    a=dict(np.load(src/'private_all129_all92_native_axis_and_design_inputs.npz'));coverage=pd.read_csv(src/'private_all92_source_design_coverage.tsv',sep='\t').set_index('axis_id');results=[];sensitivity=[];histology=[]
    for j,axis in enumerate(a['axis_ids']):
        m,g,status,valid=prepare(a,j);row=coverage.loc[axis].to_dict();row.update(axis_id=axis,status=status,effect=np.nan,CI_low=np.nan,CI_high=np.nan,wild_P=1.,wild_extreme=0,CV3_se=np.nan,CV3_t=np.nan,analytic_CV3_P=1.)
        if g is not None:
            try:fit=evaluate(g,m['y'])
            except ValueError as e:row['status']=str(e);g=None
        if g is not None:
            null=wild(g,m['y'],19999,seed(axis,'primary_WCRB'),True);full=wild(g,m['y'],19999,seed(axis,'interval_WCUB'),False);interval=confidence(fit,full);extreme=int((np.abs(null)>=abs(fit['t'])-1e-12).sum());row.update(effect=fit['effect'],CI_low=float(interval[0]),CI_high=float(interval[1]),wild_P=(1+extreme)/20000,wild_extreme=extreme,CV3_se=fit['se'],CV3_t=fit['t'],analytic_CV3_P=float(2*student_t.sf(abs(fit['t']),g['groups']-1)),max_block_leverage=g['max_block_leverage'],residual_df=m['df'],nuisance_rank=m['nuisance_rank']);np.savez_compressed(out/('private_'+axis+'_all_inference_streams.npz'),restricted_null_WCRB_t=null,full_model_WCUB_t=full,leave_profile_block_effect=fit['leave_block_effect'])
        results.append(row)
        for mode in SENSITIVITIES:
            sm,sg,why,_=prepare(a,j,mode);r=dict(axis_id=axis,sensitivity=mode,status=why,effect=np.nan,CV3_t_CI_low=np.nan,CV3_t_CI_high=np.nan)
            if sg is not None:
                try:f=evaluate(sg,sm['y'])
                except ValueError as e:r['status']=str(e);sg=None
            if sg is not None:
                critical=student_t.ppf(.975,sg['groups']-1);r.update(effect=f['effect'],CV3_t_CI_low=f['effect']-critical*f['se'],CV3_t_CI_high=f['effect']+critical*f['se'])
            sensitivity.append(r)
        for cancer in ['LUAD','LUSC']:
            take=valid&(a['cancer_type']==cancer);x=a['response'][take];y=a['standardized_axis_score'][take,j];n0=int((x==0).sum());n1=int((x==1).sum());effect=float(y[x==1].mean()-y[x==0].mean()) if min(n0,n1)>0 else np.nan;histology.append(dict(axis_id=axis,cancer_type=cancer,author_records=int(take.sum()),nonMPR=n0,pCR_or_MPR=n1,unadjusted_effect=effect,status='DESCRIPTIVE_ONLY_NO_HISTOLOGY_P_OR_REPLICATION_CLAIM' if min(n0,n1)>=10 else 'UNDER10_PER_GROUP_OR_MISSING_DESCRIPTIVE_ONLY'))
        save('private_progress.json',dict(phase='ALL92_PREDECLARED_SOURCE_CLINICAL_INFERENCE',completed_slots=j+1,total_slots=92,observed_UTC=datetime.now(timezone.utc).isoformat()))
    r=pd.DataFrame(results);assert len(r)==92;r['BH92']=adjust(r.wild_P);r['BY92']=adjust(r.wild_P,by=True);r['secondary_analytic_CV3_BH92']=adjust(r.analytic_CV3_P);r['secondary_analytic_CV3_BY92']=adjust(r.analytic_CV3_P,by=True);r.to_csv(out/'private_all92_source_consistency_results.tsv',sep='\t',index=False);pd.DataFrame(sensitivity).to_csv(out/'private_all276_fixed_sensitivity_effects.tsv',sep='\t',index=False);pd.DataFrame(histology).to_csv(out/'private_all184_histology_descriptors.tsv',sep='\t',index=False)
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()};assert sum(v['bytes'] for v in files.values())<32*1024**2;receipt=dict(status='COMPLETE_ALL92_SOURCE_CONSISTENCY_MAIN_AND_FIXED_SENSITIVITIES_FULL_AUDIT_PENDING',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),input_audit_sha256=sha(audit),null_calibration_sha256=sha(cal),files=files,main_slots=92,sensitivity_slots=276,histology_descriptors=184,original_clinical87_not_retested=True,scientific_outputs_public=False)
    save('private_terminal.json',receipt);print(json.dumps(dict(status=receipt['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
