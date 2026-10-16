"""Freeze complete author-state proportion-change clinical scope before effects."""
import hashlib,json,os,subprocess
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import rankdata

COLUMNS=['intercept','baseline_within_parent_fraction','log_pre_parent_cells','log_post_parent_cells','treatment_DplusT','HPV_positive','delta_within_parent_fraction']


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def transform(x,y):
    d=x.copy();yy=y.copy()
    if not np.isfinite(d).all() or not np.isfinite(y).all():return d,yy,'NONFINITE_REQUIRED_INPUT'
    for j in [1,2,3,6]:
        ranks=rankdata(d[:,j]);sd=ranks.std()
        if sd<1e-12:return d,yy,'CONSTANT_REQUIRED_CONTINUOUS_VARIABLE'
        d[:,j]=(ranks-ranks.mean())/sd
    if any(np.ptp(d[:,j])==0 for j in [4,5]):return d,yy,'CONSTANT_REQUIRED_TREATMENT_OR_HPV'
    ranks=rankdata(yy)
    if ranks.std()<1e-12:return d,yy,'CONSTANT_PATHOLOGY'
    yy=(ranks-ranks.mean())/ranks.std()
    if len(y)<20 or len(y)-7<10:return d,yy,'INSUFFICIENT_PATIENTS_OR_DF'
    s=np.linalg.svd(d,compute_uv=False);rank=int((s>s[0]*1e-12).sum());machine=int((s>s[0]*np.finfo(float).eps*max(d.shape)).sum())
    if rank!=machine or rank!=7:return d,yy,'NONFULL_OR_AMBIGUOUS_REQUIRED_DESIGN'
    if s[0]/s[-1]>=1e8:return d,yy,'ILL_CONDITIONED_REQUIRED_DESIGN'
    h=np.einsum('ij,ji->i',d,np.linalg.pinv(d,rcond=1e-12))
    if h.max()>=.999999:return d,yy,'EXTREME_LEVERAGE'
    return d,yy,'READY_NO_OBSERVED_CLINICAL_EFFECT'


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);source=r/'HNSCC_author_state_count_layer_preflight_20261016_v1';proof=source/'private_independent_author_state_count_audit.json';assert json.loads(proof.read_text())['status'].startswith('PASS_ALL57_');assert sha(source/'private_terminal.json')==json.loads(proof.read_text())['producer_terminal_sha256']
    clinical=r/'HNSCC_frozen_MCP_clinical_preparation_20261014_v1/private_patient_clinical_definition.tsv';c=pd.read_csv(clinical,sep='\t',keep_default_na=False,dtype={'patient_id':str}).set_index('patient_id');c['tx']=c.treatment.map({'D':0.,'D+T':1.});c['hpv']=c.HPV.str.strip().str.lower().map({'negative':0.,'positive':1.,'neg':0.,'pos':1.,'-':0.,'+':1.});c['regression']=pd.to_numeric(c.regression);assert np.isfinite(c[['tx','hpv','regression']].to_numpy()).all()
    ll=pd.read_csv(source/'private_all1995_library_author_state_fractions.tsv',sep='\t',dtype={'patient_id':str});defs=pd.read_csv(source/'private_all35_author_immune_state_fraction_definitions.tsv',sep='\t');op=pd.read_csv(source/'private_all241_state_pair_opportunities.tsv',sep='\t');features={};measure=[]
    for row in defs.itertuples():
        rows=[]
        for person,g in ll[ll.feature_id.eq(row.feature_id)].groupby('patient_id'):
            a=g[g.timepoint.eq('pre')];b=g[g.timepoint.eq('post')]
            if len(a)!=1 or len(b)!=1 or min(a.parent_cells.iloc[0],b.parent_cells.iloc[0])<20:continue
            a=a.iloc[0];b=b.iloc[0];assert person in c.index
            rows.append(dict(patient_id=person,feature_id=row.feature_id,parent=row.parent,pre_library=a.library_id,post_library=b.library_id,baseline_fraction=a.fraction,post_fraction=b.fraction,delta_fraction=b.fraction-a.fraction,pre_parent_cells=int(a.parent_cells),post_parent_cells=int(b.parent_cells),pre_state_cells=int(a.state_cells),post_state_cells=int(b.state_cells)))
        g=pd.DataFrame(rows).set_index('patient_id');assert g.index.is_unique;features[row.feature_id]=g;measure.extend(g.reset_index().to_dict('records'))
    out=r/'HNSCC_author_state_clinical_preflight_20261016_v1';assert not out.exists();fs=subprocess.check_output(['df','-P',str(r.resolve())],text=True).splitlines()[-1].split()[0];quota=subprocess.check_output(['quota','-g','-v','-w','USER002'],text=True);qr=[s.split() for s in quota.splitlines() if s.split() and s.split()[0]==fs];assert len(qr)==1;free=(int(qr[0][3])-int(qr[0][1].rstrip('*')))*1024;assert free>500_000_000_000+1_073_741_824+268_435_456
    out.mkdir(mode=0o700);(out/'private_designs').mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),count_audit_sha256=sha(proof),clinical_sha256=sha(clinical),kernel_sha256=sha(Path(__file__).with_name('HNSCC_function_clinical_kernel_v1.py')),question='All published author immune state within-parent fraction changes and continuous pathological regression; shared marginal clinical association across different parents. A cohort count question, not RNA function, absolute abundance, per-cell state reassignment, exact old CM validation, or direct LR communication.',prior_exploration='This cohort has prior programme, LR and function screening. Author-state layer was specified before inspecting latest functional clinical outcomes; present formal inference is exploratory, not untouched confirmation.',scope={'primary_cross_major_pair':221,'auxiliary_broad_single':26,'supplementary_nested_pair':20,'supplementary_nested_single':9},design_columns=COLUMNS,measurement='Raw author state/explicit parent counts, exact unique pre/post reviewed patient, both parent counts>=20. True zero state counts retained. Count QC/annotation version differs from atlas; no expression measurement or state-cell-count threshold.',controls='Own baseline fraction, own log1p pre/post parent capture, reviewed therapy and HPV. Counts have no RNA quality scores; none are fabricated from atlas. Rank/population-SD standardized continuous variables within exact cases; binary unchanged; no automatic covariate deletion.',gates='At least20 patients and10 residual df; full7-column unambiguous rank, condition<1e8, maxh<.999999. Pair same exact patient set at both ends. Unavailable slots retained in whole family.',formal_inference={'wild_draws':199999,'bootstrap_draws':2000,'bootstrap':'Patients within therapy x HPV strata; original rank scale fixed. CI only>=1900 valid full-rank/condition draws. Every LOO descriptive. Reuse audited generic HC3 target/restricted-HC2 Rademacher wild kernel; calibration required before observed inference.','pair_test':'max(P_A,P_B), no counterpart fraction control; shared marginal clinical information, not unique increment. All four frozen families corrected separately BH; BY and simultaneous .99 MC-upper BH sensitivity. No choosing the favourable family or changing scope after effects.','compositional_limits':'Fractions describe relative allocation within a parent. Closure can cause same-parent inverse patterns; primary pairs cross three distinct broad parents. Nested Treg/TAM fractions are supplementary denominators, not independent replication. Clinical two-end association need not prove direct coordination or communication.','calibration':'Three distinct ready designs chosen by smallest n/highest/median leverage without coefficients; normal,t3 unit variance,heteroskedastic normal,1024 each,1999 wild, every null rejection rate<=.075. Additional strong-target power descriptive only.'},no_new_clinical_effects_computed=True,quota_free_bytes=free,output_cap_bytes=268_435_456,reference='https://doi.org/10.1038/s41467-021-27150-6 ; compositional issue reference, not a claim that this rank-regression is scCODA')
    (out/'plan_before_author_state_clinical_effects.json').write_text(json.dumps(plan,indent=2)+'\n');models={};scope=[];links=[]
    def prepare(feature,people,slot):
        g=features[feature].loc[people];cc=c.loc[people];x=np.column_stack([np.ones(len(g)),g.baseline_fraction,np.log1p(g.pre_parent_cells),np.log1p(g.post_parent_cells),cc.tx,cc.hpv,g.delta_fraction]);y=cc.regression.to_numpy(float);strata=(2*cc.tx+cc.hpv).to_numpy(int);d,yy,status=transform(x,y)
        if status!='READY_NO_OBSERVED_CLINICAL_EFFECT':return '',status
        payload=json.dumps(people,separators=(',',':')).encode()+d.astype('<f8').tobytes()+yy.astype('<f8').tobytes()+strata.astype('<i8').tobytes();identifier=hashlib.sha256(payload).hexdigest()
        if identifier not in models:
            p=out/'private_designs'/(identifier+'.npz');np.savez_compressed(p,D=d,Y=yy,patients=np.array(people),strata=strata);h=np.einsum('ij,ji->i',d,np.linalg.pinv(d,rcond=1e-12));models[identifier]=dict(model_id=identifier,design_sha256=sha(p),n=len(people),residual_df=len(people)-7,condition=float(np.linalg.cond(d)),max_leverage=float(h.max()))
        links.append(dict(slot=slot,feature_id=feature,model_id=identifier));return identifier,status
    for row in defs.itertuples():
        family='auxiliary_broad_single' if row.role=='primary_broad' else 'supplementary_nested_single';people=sorted(features[row.feature_id].index);slot=family+'|'+row.feature_id;identifier,status=prepare(row.feature_id,people,slot);scope.append(dict(slot=slot,family=family,feature_A=row.feature_id,feature_B='',clinical_complete_patients=len(people),A_model=identifier,B_model='',status=status))
    for row in op.itertuples():
        family='primary_cross_major_pair' if row.role=='primary_cross_major' else 'supplementary_nested_pair';people=sorted(set(features[row.feature_A].index)&set(features[row.feature_B].index));assert len(people)==row.reviewed_paired_patients;slot=family+'|'+row.feature_A+'|'+row.feature_B;a,sa=prepare(row.feature_A,people,slot);b,sb=prepare(row.feature_B,people,slot);status=sa if sa==sb=='READY_NO_OBSERVED_CLINICAL_EFFECT' else 'A_'+sa+'__B_'+sb;scope.append(dict(slot=slot,family=family,feature_A=row.feature_A,feature_B=row.feature_B,clinical_complete_patients=len(people),A_model=a,B_model=b,status=status))
    scope=pd.DataFrame(scope);assert len(scope)==276 and scope.slot.is_unique and scope.family.value_counts().to_dict()==plan['scope'];scope.to_csv(out/'private_all276_state_clinical_scope.tsv',sep='\t',index=False);pd.DataFrame(models.values()).sort_values('model_id').to_csv(out/'private_unique_numeric_designs.tsv',sep='\t',index=False);pd.DataFrame(links).to_csv(out/'private_all_slot_endpoint_model_links.tsv',sep='\t',index=False);pd.DataFrame(measure).to_csv(out/'private_all_author_state_paired_measurements.tsv',sep='\t',index=False)
    terminal=dict(status='COMPLETE_ALL276_AUTHOR_STATE_CLINICAL_DESIGNS_NO_OBSERVED_EFFECTS_AUDIT_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),unique_models=len(models),files={p.name:sha(p) for p in out.iterdir() if p.is_file()},new_clinical_effects=False);(out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(terminal))


if __name__=='__main__':main()
