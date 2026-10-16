"""All4800 paired/3200 single slots, no observed clinical coefficients."""
import hashlib,itertools,json,os,subprocess
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import rankdata
from HNSCC_function_change_clinical_inputs_v1 import read,clinical_complete,raw_design,sha,NAMES


def transform(x,y):
    if len(y)==0:return None,None,'NO_CLINICALLY_COMPLETE_PATIENT'
    if not np.isfinite(x).all() or not np.isfinite(y).all():return None,None,'NONFINITE_REQUIRED_INPUT'
    d=x.copy()
    for j in [1,2,3,4,5,6,7,10]:
        col=rankdata(d[:,j]);sd=col.std()
        if sd<1e-12:return None,None,'CONSTANT_REQUIRED_CONTINUOUS_VARIABLE'
        d[:,j]=(col-col.mean())/sd
    for j in [8,9]:
        if np.ptp(d[:,j])==0:return None,None,'CONSTANT_REQUIRED_TREATMENT_OR_HPV'
    yy=rankdata(y);sd=yy.std()
    if sd<1e-12:return None,None,'CONSTANT_CONTINUOUS_PATHOLOGY'
    yy=(yy-yy.mean())/sd
    if len(y)<20 or len(y)-11<10:return d,yy,'INSUFFICIENT_PATIENTS_OR_DF'
    s=np.linalg.svd(d,compute_uv=False);r1=int(np.sum(s>s[0]*1e-12));r2=int(np.sum(s>s[0]*np.finfo(float).eps*max(d.shape)))
    if r1!=r2 or r1!=11:return d,yy,'NONFULL_OR_AMBIGUOUS_REQUIRED_DESIGN'
    if s[0]/s[-1]>=1e8:return d,yy,'ILL_CONDITIONED_REQUIRED_DESIGN'
    inv=np.linalg.pinv(d,rcond=1e-12);h=np.einsum('ij,ji->i',d,inv)
    if h.max()>=.999999:return d,yy,'EXTREME_LEVERAGE'
    return d,yy,'READY_NO_OBSERVED_CLINICAL_EFFECT'


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);data=read(r);out=r/'HNSCC_function_change_clinical_preflight_20261016_v1';assert not out.exists();fs=subprocess.check_output(['df','-P',str(r.resolve())],text=True).splitlines()[-1].split()[0];q=subprocess.check_output(['quota','-g','-v','-w','USER002'],text=True);qr=[s.split() for s in q.splitlines() if s.split() and s.split()[0]==fs];assert len(qr)==1;free=(int(qr[0][3])-int(qr[0][1].rstrip('*')))*1024;assert free>500_000_000_000+1_073_741_824+268_435_456
    out.mkdir(mode=0o700);(out/'private_designs').mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),helper_sha256=sha(Path(__file__).with_name('HNSCC_function_change_clinical_inputs_v1.py')),input_sha256=data['input_sha256'],scope='All16 masks x50 functions x6 pairs=4800 primary slots and x4 pools=3200 secondary single slots. No new clinical effects read.',design_columns=NAMES,response='Rank-standardized original continuous pathological regression; one reviewed patient per row. Each endpoint own6 quality changes, own baseline functional RNA, treatment and HPV. No counterpart function as covariate; shared/marginal clinical association not unique incremental information.',transforms='Within exact analysis patients, average ranks/population-SD standardize continuous variables. Intercept and reviewed binary treatment/HPV kept. Fixed scale for all later resampling, no automatic covariate deletion.',gates='n>=20,residualdf>=10,full unambiguous rank,condition<1e8,leverage<.999999; entire primary pair requires both endpoints ready. Missing or failed slots retained, no P yet.',future_inference='Formal HC3/restricted-HC2 wild calibration and patient-bootstrap protocol still to freeze before real effects. Exact numeric designs deduplicated by patient keys,D,Y and strata; aliases retained, no independent-evidence multiplication.',quota_free_bytes=free,output_cap_bytes=268_435_456,new_download=False)
    (out/'plan_before_any_clinical_coefficients.json').write_text(json.dumps(plan,indent=2)+'\n');pools=sorted(data['axis'].pool.unique());models={};model_rows=[];scope=[];link=[];designs={}
    def prepare(mask,ref,pool,patients,slot):
        raw,y,strata=raw_design(data,mask,ref,pool,patients);d,yy,status=transform(raw,y)
        if status!='READY_NO_OBSERVED_CLINICAL_EFFECT':return dict(model_id='',status=status,n=len(patients),residual_df=len(patients)-11)
        payload=json.dumps(patients,separators=(',',':')).encode()+d.astype('<f8').tobytes()+yy.astype('<f8').tobytes()+strata.astype('<i8').tobytes();identifier=hashlib.sha256(payload).hexdigest()
        if identifier not in models:
            models[identifier]=len(models);path=out/'private_designs'/(identifier+'.npz');np.savez_compressed(path,D=d,Y=yy,patients=np.array(patients),strata=strata);inverse=np.linalg.pinv(d,rcond=1e-12);h=np.einsum('ij,ji->i',d,inverse);model_rows.append(dict(model_id=identifier,design_sha256=sha(path),n=len(patients),residual_df=len(patients)-11,condition=float(np.linalg.cond(d)),max_leverage=float(h.max())));designs[identifier]=(d,yy)
        link.append(dict(slot=slot,pool=pool,model_id=identifier));return dict(model_id=identifier,status=status,n=len(patients),residual_df=len(patients)-11)
    for mask in data['masks'].mask_id:
        for ref in data['references'].reference:
            for pool in pools:
                key=(mask,ref,pool);patients=clinical_complete(data,sorted(data['groups'][key].index));slot='single|'+str(mask)+'|'+ref+'|'+pool;result=prepare(mask,ref,pool,patients,slot);scope.append(dict(slot=slot,family='single3200',mask_id=mask,reference=ref,pool_A=pool,pool_B='',capture_patients=len(data['groups'][key]),clinical_complete_patients=len(patients),A_model=result['model_id'],B_model='',status=result['status']))
            for a,b in itertools.combinations(pools,2):
                original=sorted(set(data['groups'][(mask,ref,a)].index)&set(data['groups'][(mask,ref,b)].index));patients=clinical_complete(data,original);slot='pair|'+str(mask)+'|'+ref+'|'+a+'|'+b;A=prepare(mask,ref,a,patients,slot);B=prepare(mask,ref,b,patients,slot);status='READY_NO_OBSERVED_CLINICAL_EFFECT' if A['status']==B['status']=='READY_NO_OBSERVED_CLINICAL_EFFECT' else 'A_'+A['status']+'__B_'+B['status'];scope.append(dict(slot=slot,family='pair4800',mask_id=mask,reference=ref,pool_A=a,pool_B=b,capture_patients=len(original),clinical_complete_patients=len(patients),A_model=A['model_id'],B_model=B['model_id'],status=status))
    table=pd.DataFrame(scope);assert len(table)==8000 and table.slot.is_unique and table.family.value_counts().to_dict()=={'pair4800':4800,'single3200':3200};model=pd.DataFrame(model_rows).sort_values('model_id');pd.DataFrame(link).to_csv(out/'private_all_slot_endpoint_model_links.tsv',sep='\t',index=False);table.to_csv(out/'private_all8000_clinical_scope.tsv',sep='\t',index=False);model.to_csv(out/'private_unique_numeric_designs.tsv',sep='\t',index=False)
    # Composition turnover is metadata, not an additional adjusted variable or effect.
    counts=pd.read_csv(data['paths']['node_counts'],sep='\t',dtype={'library_id':str,'node_id':str});counts=counts[counts.arm.eq('archived_immune_atlas')&counts.dataset.eq('GSE286827')];defs=pd.read_csv(data['paths']['pool_definitions'],sep='\t',dtype=str);defs=defs[defs.arm.eq('archived_immune_atlas')&defs.pool.isin(pools)];join=counts.merge(defs[['pool','node_id']],on='node_id',how='inner',validate='many_to_one');axis=data['axis'];composition=[]
    for pool in pools:
        pivot=join[join.pool.eq(pool)].pivot_table(index='library_id',columns='node_id',values='n_cells',aggfunc='sum',fill_value=0)
        for patient,g in axis[axis.pool.eq(pool)].groupby('patient_block'):
            pre=g[g.timepoint.eq('pre')];post=g[g.timepoint.eq('post')]
            if len(pre)!=1 or len(post)!=1:continue
            pp=pivot.loc[pre.library_id.iloc[0]].to_numpy(float);qq=pivot.loc[post.library_id.iloc[0]].to_numpy(float);assert pp.sum()==pre.n_cells.iloc[0] and qq.sum()==post.n_cells.iloc[0]
            composition.append(dict(patient_block=patient,pool=pool,pre_cells=int(pp.sum()),post_cells=int(qq.sum()),frozen_fine_nodes=int(len(pp)),composition_total_variation=float(np.abs(pp/pp.sum()-qq/qq.sum()).sum()/2),state_change_not_separated=True))
    pd.DataFrame(composition).to_csv(out/'private_all_patient_pool_fine_composition_turnover.tsv',sep='\t',index=False)
    assert sum(p.stat().st_size for p in out.rglob('*') if p.is_file())<plan['output_cap_bytes'];terminal=dict(status='COMPLETE_ALL8000_CLINICAL_SCOPE_AND_UNIQUE_FULL_DESIGNS_NO_OBSERVED_EFFECTS_AUDIT_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),scope_rows=len(table),unique_numeric_designs=len(model),files={p.name:sha(p) for p in out.iterdir() if p.is_file()},new_observed_clinical_coefficients_P=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(terminal))


if __name__=='__main__':main()
