"""Complete prespecified six paired-tissue tests and point-only follow-up."""
import json
import os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from GSE176021_paired_tissue_math_v1 import wilcoxon_signed,rho,permutations,mc_probability,bootstrap_indices,bootstrap_effects,adjust
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def point(x,y,cochange):return rho(x,y) if cochange else float(np.median(y))


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,8*1024**2);src=root/'GSE176021_paired_tissue_inputs_20261019_v1';ia=root/'GSE176021_paired_tissue_input_audit_20261019_v1/private_terminal.json';assert sha(ia)=='1998d626c57e0165ab4b66573170d16a904c600e442823b86dfb37dd39ea5f89';ma=root/'GSE176021_paired_tissue_math_check_20261019_v1/private_terminal.json';assert sha(ma)=='2edb76d1215c1b88e7563e19382088b84c26bfcff61f00b87de33110a9437b04'
    t=json.loads((src/'private_terminal.json').read_text())
    for f,v in t['files'].items():assert sha(src/f)==v['sha256']
    out=root/'GSE176021_paired_tissue_associations_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_all6_true_tissue_tests.json').write_text(json.dumps({'created_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'input_acceptance_sha256':sha(ia),'math_acceptance_sha256':sha(ma),'resource_guard':guard,'main':'All6 prospective questions, same10case units, metadata-preeligible universe6 includes unestimable P1. Exact Wilcoxon ranks vs case-permuted cochange; no endpoint/scope/control changes.','all_random_streams':'Every exact sign statistic,2x19999 cochange permutations,6x2000whole-pairbootstrap,complete sharedindices/caseorder/LOO stored for independent replay.','followup':'All54 predefined point sensitivities, all30 scope/tissue/program coactivity positions with LOO only. No new coactivity P, no MPR. Raw case rows and learned data server-private.'},indent=2)+'\n')
    registry=pd.read_csv(src/'private_all6_registered_primary_tests.tsv',sep='\t');sensitivity=pd.read_csv(src/'private_all_prespecified_case_sensitivity_inputs.tsv',sep='\t');case=pd.read_csv(src/'private_all_scope_case_tissue_values.tsv',sep='\t');ix=bootstrap_indices(10);rows=[];sens=[];exactcount=0;randomcount=0
    for r in registry.itertuples(index=False):
        a=np.load(src/('private_'+r.model_id+'_inputs.npz'));x=a['x'];y=a['y'];cases=a['cases'];co=r.outcome=='RE_tissue_cochange';assert len(cases)==10;effect=point(x,y,co);record=r._asdict();record.update(effect=effect,P=np.nan,MC_upper_P=np.nan,sign_P=np.nan,status='ESTIMABLE',n_pairs=10);stream={'cases':cases,'bootstrap_indices':ix}
        if co and not np.isfinite(effect):record.update(status='CONSTANT_COCHANGE_TARGET');bs=np.full(2000,np.nan);loo=np.full(10,np.nan)
        else:
            if co:
                vals=permutations(x,y,r.model_id);p,u,k=mc_probability(vals,effect);record.update(P=p,MC_upper_P=u,statistic=effect,exceedances=k,exact=False);randomcount+=len(vals)
            else:
                z=wilcoxon_signed(y);vals=z['statistics'];record.update(P=z['P'],MC_upper_P=z['P'],sign_P=z['sign_P'],statistic=z['statistic'],nonzero_differences=z['nonzero'],exact=z['exact']);exactcount+=len(vals)
            stream['null_statistics']=vals;bs=bootstrap_effects(x,y,ix,co);loo=np.array([point(np.delete(x,j),np.delete(y,j),co) for j in range(10)])
        finite=np.isfinite(bs);ci=np.quantile(bs[finite],[.025,.975]) if finite.sum()>=1900 else [np.nan,np.nan];record.update(CI_low=float(ci[0]),CI_high=float(ci[1]),valid_bootstrap=int(finite.sum()),LOO_min=float(np.nanmin(loo)) if np.isfinite(loo).any() else np.nan,LOO_max=float(np.nanmax(loo)) if np.isfinite(loo).any() else np.nan,LOO_same_sign=bool(np.isfinite(loo).all() and (loo*effect>0).all()));rows.append(record);stream.update(bootstrap_effects=bs,leave_one_effects=loo);np.savez_compressed(out/('private_'+r.model_id+'_complete_stream.npz'),**stream)
        for mode,d in sensitivity[sensitivity.model_id.eq(r.model_id)].groupby('mode',sort=True):
            good=np.isfinite(d.x)&np.isfinite(d.y);complete=bool(good.all());value=point(d.x.to_numpy(),d.y.to_numpy(),co) if complete else np.nan;sens.append({'model_id':r.model_id,'mode':mode,'effect':value,'registered_case_rows':len(d),'usable_case_rows':int(good.sum()),'primary_N':10,'status':'COMPLETE_POINT_ONLY_NO_P' if complete else 'INCOMPLETE_PRIMARY_CASE_UNIVERSE_NO_POINT','no_new_P':True})
    result=pd.DataFrame(rows);P=result.P.fillna(1).to_numpy();result['BH']=adjust(P);result['BY']=adjust(P,True);result['MC_upper_BH']=adjust(result.MC_upper_P.fillna(1));result.loc[result.P.isna(),['BH','BY','MC_upper_BH']]=np.nan;result['sign_BH']=np.nan;mask=result.sign_P.notna();result.loc[mask,'sign_BH']=adjust(result.loc[mask,'sign_P']);result['CI_excludes_zero']=result.CI_low*result.CI_high>0;result['priority_association']=result.BH.lt(.05)&result.MC_upper_BH.lt(.05)&result.CI_excludes_zero&result.LOO_same_sign;assert len(result)==6 and len(sens)==54
    follow=[];source=pd.read_csv(root/'GSE176021_author_lineage_RNA_scope_20261019_v2/private_all_source_materials.tsv',sep='\t');primary=case[case['mode'].eq('primary')]
    for sc in ['legacy_capture_anchor','author_confirmed_capture','all_author_lineage']:
        for tissue in sorted(source.tissue.unique()):
            for program in ['MCP1','MCP2']:
                d=primary[primary.scope.eq(sc)&primary.tissue.eq(tissue)&primary.program.eq(program)];a=d[d['view'].eq('R_FOXP3_capture')].set_index('case_block').value;b=d[d['view'].eq('E_effector_capture')].set_index('case_block').value;keys=sorted(set(a.index)&set(b.index));x=a.reindex(keys).to_numpy();y=b.reindex(keys).to_numpy();n=len(keys);e=rho(x,y) if n>=3 else np.nan;loo=np.array([rho(np.delete(x,j),np.delete(y,j)) for j in range(n)]) if n>=4 else np.array([]);follow.append({'scope':sc,'tissue':tissue,'program':program,'case_units':n,'rho':e,'LOO_min':float(np.nanmin(loo)) if np.isfinite(loo).any() else np.nan,'LOO_max':float(np.nanmax(loo)) if np.isfinite(loo).any() else np.nan,'status':'POINT_ONLY_NO_NEW_INFERENCE' if n>=10 else 'BELOW10_DESCRIPTIVE_ONLY' if n>=3 else 'INSUFFICIENT_CASES_NO_POINT','no_P':True})
    result.to_csv(out/'private_all6_primary_paired_tissue_results.tsv',sep='\t',index=False);pd.DataFrame(sens).to_csv(out/'private_all54_point_sensitivity_results.tsv',sep='\t',index=False);pd.DataFrame(follow).to_csv(out/'private_all30_coactivity_context_points_only.tsv',sep='\t',index=False)
    files={p.name:{'sha256':sha(p),'bytes':p.stat().st_size} for p in out.iterdir() if p.is_file()};assert sum(v['bytes'] for v in files.values())<8*1024**2
    d={'status':'COMPLETE_ALL6_TISSUE_TESTS54_SENSITIVITIES30_CONTEXT_POINTS_AND_FULL_STREAMS_AUDIT_PENDING','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'files':files,'exact_sign_statistics':exactcount,'permutation_statistics':randomcount,'bootstrap_effects':12000,'no_MPR_or_old_coactivity_P':True}
    (out/'private_terminal.json').write_text(json.dumps(d,indent=2)+'\n');print(json.dumps({'status':d['status'],'terminal_sha256':sha(out/'private_terminal.json')}))


if __name__=='__main__':main()
