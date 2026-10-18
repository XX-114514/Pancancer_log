"""Independent full deleted-QR replay of all92 inferential streams and outputs."""
import json
import os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import t as student_t
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC_immune_RNA_MCP_cohort_io_v2 import wait
from NSCLC_frozen_axis_block_audit_math_v1 import prepare_axis,evaluate,wild
from NSCLC_frozen_axis_block_math_v2 import seed


def correction(p,by=False):
    values=np.asarray(p,float);order=np.argsort(values,kind='stable');m=len(values);q=values[order]*m/np.arange(1,m+1)
    if by:q=q*np.sum(1/np.arange(1,m+1))
    out=np.empty(m);out[order]=np.minimum(1,np.minimum.accumulate(q[::-1])[::-1]);return out


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,2*1024**2);out=root/'NSCLC_frozen92_axis_source_result_audit_20261019_v1';assert not out.exists();out.mkdir(mode=0o700);src=root/'NSCLC_frozen92_axis_source_consistency_20261019_v1';terminal=wait(root,out,'run_NSCLC_frozen92_axis_source_consistency_v1',src/'private_terminal.json',{'COMPLETE_ALL92_SOURCE_CONSISTENCY_MAIN_AND_FIXED_SENSITIVITIES_FULL_AUDIT_PENDING'})
    for name,d in terminal['files'].items():assert sha(src/name)==d['sha256']
    inp=root/'NSCLC_frozen92_axis_source_inputs_20261019_v2';proof=root/'NSCLC_frozen92_axis_source_input_audit_20261019_v1/private_terminal.json';assert sha(proof)=='b07c736a359df5ec4138c7743c2a978bd7f620e867c59344d06fcc466057ecbd';cal=root/'NSCLC_frozen92_axis_block_null_calibration_20261019_v1/private_terminal.json';ct=json.loads(cal.read_text());assert ct['status']=='ACCEPTED_PREDECLARED_NULL_DIAGNOSTICS_BEFORE_ALL92_TRUE_EFFECTS' and ct['all_synthetic_random_statistics_independently_replayed'] and ct['input_audit_sha256']==sha(proof);assert terminal['input_audit_sha256']==sha(proof) and terminal['null_calibration_sha256']==sha(cal)
    a=dict(np.load(inp/'private_all129_all92_native_axis_and_design_inputs.npz'));r=pd.read_csv(src/'private_all92_source_consistency_results.tsv',sep='\t');s=pd.read_csv(src/'private_all276_fixed_sensitivity_effects.tsv',sep='\t').set_index(['axis_id','sensitivity']);h=pd.read_csv(src/'private_all184_histology_descriptors.tsv',sep='\t').set_index(['axis_id','cancer_type']);assert r.axis_id.tolist()==a['axis_ids'].tolist() and len(s)==276 and len(h)==184;P=[];analytic=[];maximum=0.;streams=0
    for j,axis in enumerate(a['axis_ids']):
        m,g,status,valid=prepare_axis(a,j);row=r.iloc[j]
        if g is not None:
            try:f=evaluate(g,m['y'])
            except ValueError as e:status=str(e);g=None
        assert row.status==status;assert row.author_records==int(valid.sum()) and row.profile_blocks==len(set(a['profile_blocks'][valid]))
        if g is None:
            assert row.wild_P==1 and row.analytic_CV3_P==1 and np.isnan(row.effect) and not (src/('private_'+axis+'_all_inference_streams.npz')).exists();P.append(1.);analytic.append(1.)
        else:
            f=evaluate(g,m['y']);null=wild(g,m['y'],19999,seed(axis,'primary_WCRB'),True);full=wild(g,m['y'],19999,seed(axis,'interval_WCUB'),False);saved=dict(np.load(src/('private_'+axis+'_all_inference_streams.npz')));np.testing.assert_allclose(saved['restricted_null_WCRB_t'],null,atol=1e-8,rtol=1e-8);np.testing.assert_allclose(saved['full_model_WCUB_t'],full,atol=1e-8,rtol=1e-8);np.testing.assert_allclose(saved['leave_profile_block_effect'],f['leave_block_effect'],atol=1e-9,rtol=1e-9);maximum=max(maximum,float(np.max(np.abs(saved['restricted_null_WCRB_t']-null))),float(np.max(np.abs(saved['full_model_WCUB_t']-full))));extreme=int((np.abs(null)>=abs(f['t'])-1e-12).sum());p=(1+extreme)/20000;lo,hi=np.quantile(full,[.025,.975]);interval=[f['effect']-hi*f['se'],f['effect']-lo*f['se']];ap=float(2*student_t.sf(abs(f['t']),g['groups']-1));np.testing.assert_allclose([row.effect,row.CV3_se,row.CV3_t,row.CI_low,row.CI_high,row.wild_P,row.analytic_CV3_P],[f['effect'],f['se'],f['t'],*interval,p,ap],atol=1e-8,rtol=1e-8);assert row.wild_extreme==extreme and row.nuisance_rank==m['nuisance_rank'] and row.residual_df==m['df'];P.append(p);analytic.append(ap);streams+=39998
        for mode in ['unadjusted','clinical_RNA_counts_no_fractions','main_plus_total_communication']:
            mm,gg,why,_=prepare_axis(a,j,mode);rr=s.loc[(axis,mode)]
            if gg is not None:
                try:ff=evaluate(gg,mm['y'])
                except ValueError as e:why=str(e);gg=None
            assert rr.status==why
            if gg is None:assert np.isnan(rr.effect)
            else:
                ff=evaluate(gg,mm['y']);critical=student_t.ppf(.975,gg['groups']-1);np.testing.assert_allclose([rr.effect,rr.CV3_t_CI_low,rr.CV3_t_CI_high],[ff['effect'],ff['effect']-critical*ff['se'],ff['effect']+critical*ff['se']],atol=1e-8,rtol=1e-8)
        for cancer in ['LUAD','LUSC']:
            take=valid&(a['cancer_type']==cancer);x=a['response'][take];y=a['standardized_axis_score'][take,j];n0=int((x==0).sum());n1=int((x==1).sum());effect=float(y[x==1].mean()-y[x==0].mean()) if min(n0,n1)>0 else np.nan;rr=h.loc[(axis,cancer)];assert rr.author_records==int(take.sum()) and rr.nonMPR==n0 and rr.pCR_or_MPR==n1;np.testing.assert_allclose(rr.unadjusted_effect,effect,atol=1e-12,rtol=1e-12,equal_nan=True)
        (out/'private_progress.json').write_text(json.dumps(dict(phase='ALL92_FULL_RANDOM_STATISTICS_AND_FIXED_SENSITIVITIES_REPLAY',completed_slots=j+1,total_slots=92,observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n')
    for name,values in [('wild_P',P),('analytic_CV3_P',analytic)]:np.testing.assert_allclose(r[name],values,atol=1e-12,rtol=1e-12)
    for name,values,by in [('BH92',P,False),('BY92',P,True),('secondary_analytic_CV3_BH92',analytic,False),('secondary_analytic_CV3_BY92',analytic,True)]:np.testing.assert_allclose(r[name],correction(values,by),atol=1e-12,rtol=1e-12)
    result=dict(status='ACCEPTED_ALL92_MAIN_ALL_RANDOM_WCRB_WCUB_LEAVEBLOCK_BH_BY_276SENSITIVITIES_184HISTOLOGY_DESCRIPTORS',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),producer_terminal_sha256=sha(src/'private_terminal.json'),input_audit_sha256=sha(proof),resource_guard=guard,main_slots=92,sensitivity_slots=276,histology_descriptors=184,random_statistics_compared=streams,maximum_independent_statistic_difference=maximum,original_clinical87_not_retested=True,scientific_outputs_public=False)
    assert sum(p.stat().st_size for p in out.iterdir() if p.is_file())+sum(p.stat().st_size for p in proof.parent.iterdir() if p.is_file())<2*1024**2;(out/'private_terminal.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(status=result['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
