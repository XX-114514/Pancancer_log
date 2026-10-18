"""Complete prospective 400-slot functional association family and all streams."""
import json
import os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from GSE176021_disjoint_function_math_v1 import observed,permutation_numerators,probability,bootstrap_indices,bootstrap,response_quality_point,adjust
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,96*1024**2)
    ia=root/'GSE176021_disjoint_function_input_audit_20261019_v1/private_terminal.json';assert sha(ia)=='c5ea0671b8fd14bef17033be3d6089df56aeb5a9aea1e15f975b5ca908579b85'
    ma=root/'GSE176021_disjoint_function_math_check_20261019_v1/private_terminal.json';assert sha(ma)=='f05f2b92c4aeb0a663cb37173b70ff842541553462a1632ae472442d4f37aa2a'
    assert sha(Path(__file__).parent/'GSE176021_disjoint_function_math_v1.py')=='2eafd1e53d95e6778706aa4e88cd20d4b1d996c10126a2b78fcecc03d2e7df7e'
    src=root/'GSE176021_disjoint_function_inputs_20261019_v1';t=json.loads((src/'private_terminal.json').read_text())
    for name,v in t['files'].items():assert sha(src/name)==v['sha256']
    out=root/'GSE176021_disjoint_function_associations_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    bound=400*(99999*2+2000*8+10*8+4096)+4*1024**2;assert bound<96*1024**2
    (out/'plan_before_true400_functional_associations.json').write_text(json.dumps({'created_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'input_acceptance_sha256':sha(ia),'math_acceptance_sha256':sha(ma),'resource_guard':guard,'conservative_uncompressed_output_bound_bytes':bound,'prespecified_contract_sha256':sha(src/'plan_before_all400_functional_associations.json'),'scope':'All400 slots plus3200 point-only sensitivities. Exact integer permutation numerators/all2000bootstrap/all10LOO stored for complete replay. No endpoint/control/rank/threshold changes; no old family retry.','scientific_read':'Only after complete independent replay; all individual/learned/scientific tables private.'},indent=2)+'\n')
    reg=pd.read_csv(src/'private_all400_registered_functional_questions.tsv',sep='\t');a=dict(np.load(src/'private_all400_case_inputs.npz'));idx=bootstrap_indices();rows=[];sens=[];random_total=0
    np.savez_compressed(out/'private_common_case_and_bootstrap_arrays.npz',cases=a['cases'],bootstrap_indices=idx)
    for i,r in enumerate(reg.itertuples(index=False)):
        x,y=a['x'][i],a['y'][i];effect,numerator,den=observed(x,y);record=r._asdict();record.update(effect=effect,observed_numerator=numerator,rank_denominator=den,status='ESTIMABLE' if den else 'CONSTANT_INPUT_NO_PRIMARY_ASSOCIATION',P=np.nan,MC_upper_P=np.nan,exceedances=np.nan)
        if den:
            null=permutation_numerators(x,y,r.model_id);p,u,k=probability(null,numerator);record.update(P=p,MC_upper_P=u,exceedances=k);bs=bootstrap(x,y,idx);loo=np.array([observed(np.delete(x,j),np.delete(y,j))[0] for j in range(10)]);random_total+=len(null)
        else:null=np.empty(0,np.int16);bs=np.full(2000,np.nan);loo=np.full(10,np.nan)
        good=np.isfinite(bs);ci=np.quantile(bs[good],[.025,.975]) if good.sum()>=1900 else [np.nan,np.nan]
        record.update(CI_low=float(ci[0]),CI_high=float(ci[1]),valid_bootstrap=int(good.sum()),LOO_min=float(np.nanmin(loo)) if np.isfinite(loo).any() else np.nan,LOO_max=float(np.nanmax(loo)) if np.isfinite(loo).any() else np.nan,LOO_same_sign=bool(np.isfinite(loo).all() and (loo*effect>0).all()))
        rows.append(record);np.savez_compressed(out/('private_'+r.model_id+'_complete_stream.npz'),permutation_numerators=null,bootstrap_effects=bs,leave_one_effects=loo)
        for j,mode in enumerate(a['modes']):
            if mode=='response_quality_partial_two':e,status=response_quality_point(x,y,a['response_quality'][i])
            else:e,_,d=observed(a['sensitivity_x'][i,j],a['sensitivity_y'][i,j]);status='COMPLETE_POINT_ONLY_NO_P' if d else 'CONSTANT_SENSITIVITY_NO_POINT'
            sens.append({'model_id':r.model_id,'mode':mode,'effect':e,'status':status,'case_units':10,'no_new_P':True})
        if (i+1)%25==0:(out/'private_progress.json').write_text(json.dumps({'phase':'FULL400_FIXED_FUNCTIONAL_FAMILY','positions_completed':i+1,'positions':400,'observed_UTC':datetime.now(timezone.utc).isoformat()})+'\n')
    result=pd.DataFrame(rows);result['BH']=adjust(result.P.fillna(1));result['BY']=adjust(result.P.fillna(1),True);result['MC_upper_BH']=adjust(result.MC_upper_P.fillna(1));result.loc[result.P.isna(),['BH','BY','MC_upper_BH']]=np.nan
    result['CI_excludes_zero']=result.CI_low*result.CI_high>0;result['priority_functional_association']=result.BH.lt(.05)&result.MC_upper_BH.lt(.05)&result.CI_excludes_zero&result.LOO_same_sign
    assert len(result)==400 and len(sens)==3200
    result.to_csv(out/'private_all400_primary_functional_results.tsv',sep='\t',index=False);pd.DataFrame(sens).to_csv(out/'private_all3200_point_sensitivity_results.tsv',sep='\t',index=False)
    files={p.name:{'sha256':sha(p),'bytes':p.stat().st_size} for p in out.iterdir() if p.is_file()};actual=sum(v['bytes'] for v in files.values());assert actual<bound
    receipt={'status':'COMPLETE_ALL400_FUNCTIONAL_ASSOCIATIONS3200_SENSITIVITIES_AND_FULL_STREAMS_AUDIT_PENDING','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'files':files,'primary_positions':400,'point_sensitivities':3200,'permutation_statistics':random_total,'bootstrap_positions':800000,'leave_one_positions':4000,'actual_bytes':actual,'no_MPR_or_receptor_binding_claim':True}
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'status':receipt['status'],'terminal_sha256':sha(out/'private_terminal.json')}))


if __name__=='__main__':main()
