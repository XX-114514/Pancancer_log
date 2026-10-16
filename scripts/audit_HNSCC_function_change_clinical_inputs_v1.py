"""Independently reconstruct every logical clinical design and turnover value."""
import hashlib,itertools,json,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from HNSCC_function_change_clinical_inputs_v1 import read,sha,QC,NAMES


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=r/'HNSCC_function_change_clinical_preflight_20261016_v1';rc=json.loads((r/'HNSCC_receiver_followup_job_control_20261016_v2/preflight_HNSCC_function_change_clinical_v1/launch_receipt.json').read_text());deadline=time.monotonic()+14400
    while not (out/'private_terminal.json').exists():
        p=Path('/proc')/str(rc['pid'])/'stat';assert p.exists();a=p.read_text().rsplit(')',1)[1].split();assert a[0]!='Z' and a[19]==str(rc['proc_start_ticks']) and time.monotonic()<deadline;time.sleep(5)
    terminal=json.loads((out/'private_terminal.json').read_text());plan=json.loads((out/'plan_before_any_clinical_coefficients.json').read_text());assert rc['source_sha256']==terminal['code_sha256']==plan['code_sha256']
    for name,h in terminal['files'].items():assert sha(out/name)==h
    data=read(r);assert data['input_sha256']==plan['input_sha256'] and sha(Path(__file__).with_name('HNSCC_function_change_clinical_inputs_v1.py'))==plan['helper_sha256']
    models=pd.read_csv(out/'private_unique_numeric_designs.tsv',sep='\t').set_index('model_id');scope=pd.read_csv(out/'private_all8000_clinical_scope.tsv',sep='\t',keep_default_na=False);links=pd.read_csv(out/'private_all_slot_endpoint_model_links.tsv',sep='\t');assert len(scope)==8000 and scope.slot.is_unique and models.index.is_unique
    cache={}
    for identifier,row in models.iterrows():
        p=out/'private_designs'/(identifier+'.npz');assert sha(p)==row.design_sha256;archive=np.load(p,allow_pickle=False);assert set(archive.files)=={'D','Y','patients','strata'};d=archive['D'];y=archive['Y'];patients=archive['patients'].astype(str).tolist();strata=archive['strata'];payload=json.dumps(patients,separators=(',',':')).encode()+d.astype('<f8').tobytes()+y.astype('<f8').tobytes()+strata.astype('<i8').tobytes();assert hashlib.sha256(payload).hexdigest()==identifier
        assert len(y)==row.n and len(y)-11==row.residual_df;inv=np.linalg.pinv(d,rcond=1e-12);assert np.isclose(np.linalg.cond(d),row.condition,rtol=1e-12) and np.isclose(np.sum(d*inv.T,axis=1).max(),row.max_leverage,rtol=1e-12);cache[identifier]=(d,y,patients,strata)
    axis=data['axis'];clinical=data['clinical'];mapping=axis[['patient_block','patient_id']].drop_duplicates().set_index('patient_block').patient_id;expected_slots=set();expected_links=set();seen=set();pools=sorted(axis.pool.unique())
    for mask in data['masks'].mask_id:
        for ref in data['references'].reference:
            for pool in pools:expected_slots.add('single|'+str(mask)+'|'+ref+'|'+pool)
            for a,b in itertools.combinations(pools,2):expected_slots.add('pair|'+str(mask)+'|'+ref+'|'+a+'|'+b)
    assert set(scope.slot)==expected_slots
    for row in scope.itertuples():
        ga=data['groups'][(row.mask_id,row.reference,row.pool_A)];original=set(ga.index)
        if row.family=='pair4800':original &= set(data['groups'][(row.mask_id,row.reference,row.pool_B)].index)
        assert row.capture_patients==len(original);people=[]
        for person in sorted(original):
            c=clinical.loc[mapping.loc[person]]
            if np.isfinite([c.regression,c.HPV_numeric,c.treatment_numeric]).all():people.append(person)
        assert row.clinical_complete_patients==len(people);statuses=[]
        for pool,identifier in [(row.pool_A,row.A_model)]+([(row.pool_B,row.B_model)] if row.family=='pair4800' else []):
            g=data['groups'][(row.mask_id,row.reference,pool)].loc[people];pre=axis.loc[g.pre_pool_index];post=axis.loc[g.post_pool_index];c=clinical.loc[pre.patient_id.to_numpy()];assert pre.patient_block.tolist()==post.patient_block.tolist()==people
            delta=post[QC].to_numpy(float)-pre[QC].to_numpy(float);delta[:,0]=np.log1p(post.n_cells.to_numpy(float))-np.log1p(pre.n_cells.to_numpy(float));baseline=data['scores']['ULM_score'][data['mask_index'][row.mask_id],g.pre_pool_index.to_numpy(int),data['reference_index'][row.reference]];raw=np.column_stack([np.ones(len(people)),delta,baseline,c.treatment_numeric,c.HPV_numeric,g.rounded_difference]);y=c.regression.to_numpy(float);d=raw.copy();status='READY_NO_OBSERVED_CLINICAL_EFFECT'
            if not len(people):status='NO_CLINICALLY_COMPLETE_PATIENT'
            elif not np.isfinite(d).all() or not np.isfinite(y).all():status='NONFINITE_REQUIRED_INPUT'
            else:
                for j in [1,2,3,4,5,6,7,10]:
                    values=pd.Series(d[:,j]).rank().to_numpy();sd=values.std()
                    if sd<1e-12:status='CONSTANT_REQUIRED_CONTINUOUS_VARIABLE';break
                    d[:,j]=(values-values.mean())/sd
                if status.startswith('READY') and (np.ptp(d[:,8])==0 or np.ptp(d[:,9])==0):status='CONSTANT_REQUIRED_TREATMENT_OR_HPV'
                yrank=pd.Series(y).rank().to_numpy();sd=yrank.std()
                if status.startswith('READY') and sd<1e-12:status='CONSTANT_CONTINUOUS_PATHOLOGY'
                if status.startswith('READY'):
                    y=(yrank-yrank.mean())/sd
                    if len(people)<20 or len(people)-11<10:status='INSUFFICIENT_PATIENTS_OR_DF'
                    else:
                        s=np.linalg.svd(d,compute_uv=False);r1=int(np.sum(s>s[0]*1e-12));r2=int(np.sum(s>s[0]*np.finfo(float).eps*max(d.shape)))
                        if r1!=r2 or r1!=11:status='NONFULL_OR_AMBIGUOUS_REQUIRED_DESIGN'
                        elif s[0]/s[-1]>=1e8:status='ILL_CONDITIONED_REQUIRED_DESIGN'
                        elif np.sum(d*np.linalg.pinv(d,rcond=1e-12).T,axis=1).max()>=.999999:status='EXTREME_LEVERAGE'
            statuses.append(status)
            if status.startswith('READY'):
                assert identifier in cache;dd,yy,pp,ss=cache[identifier];assert pp==people and np.allclose(dd,d,rtol=0,atol=2e-14) and np.allclose(yy,y,rtol=0,atol=2e-14) and np.array_equal(ss,(2*c.treatment_numeric+c.HPV_numeric).to_numpy(int));expected_links.add((row.slot,pool,identifier));seen.add(identifier)
            else:assert identifier==''
        expected_status=statuses[0] if row.family=='single3200' else 'READY_NO_OBSERVED_CLINICAL_EFFECT' if all(s.startswith('READY') for s in statuses) else 'A_'+statuses[0]+'__B_'+statuses[1]
        assert row.status==expected_status
    assert seen==set(models.index) and expected_links==set(map(tuple,links[['slot','pool','model_id']].to_numpy())) and not links.duplicated(['slot','pool']).any()
    # Independent per-node fraction dictionaries include uncaptured nodes as zeros.
    counts=pd.read_csv(data['paths']['node_counts'],sep='\t',dtype={'library_id':str,'node_id':str});counts=counts[counts.arm.eq('archived_immune_atlas')&counts.dataset.eq('GSE286827')];defs=pd.read_csv(data['paths']['pool_definitions'],sep='\t',dtype=str);defs=defs[defs.arm.eq('archived_immune_atlas')&defs.pool.isin(pools)];lookup={(x.library_id,x.node_id):int(x.n_cells) for x in counts.itertuples()};turnover=pd.read_csv(out/'private_all_patient_pool_fine_composition_turnover.tsv',sep='\t',dtype={'patient_block':str});assert not turnover.duplicated(['patient_block','pool']).any()
    for row in turnover.itertuples():
        g=axis[axis.patient_block.eq(row.patient_block)&axis.pool.eq(row.pool)];pre=g[g.timepoint.eq('pre')].iloc[0];post=g[g.timepoint.eq('post')].iloc[0];nodes=set(defs.loc[defs.pool.eq(row.pool),'node_id']);p=np.array([lookup.get((pre.library_id,n),0) for n in sorted(nodes)]);q=np.array([lookup.get((post.library_id,n),0) for n in sorted(nodes)]);assert p.sum()==row.pre_cells and q.sum()==row.post_cells;assert np.isclose(np.abs(p/p.sum()-q/q.sum()).sum()/2,row.composition_total_variation,atol=1e-14) and row.state_change_not_separated
    proof=dict(status='PASS_ALL8000_CLINICAL_SCOPE_ALL_NUMERIC_DESIGNS_EXACT_PATIENTS_NINE_NUISANCES_AND_COMPOSITION_METADATA_NO_EFFECTS',finished_UTC=datetime.now(timezone.utc).isoformat(),audit_source_sha256=sha(__file__),producer_terminal_sha256=sha(out/'private_terminal.json'),scope_rows=len(scope),unique_numeric_designs=len(models),endpoint_links=len(links),turnover_rows_checked=len(turnover),observed_clinical_coefficients_or_P_computed=False)
    p=out/'private_independent_clinical_input_audit.json';assert not p.exists();p.write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof))


if __name__=='__main__':main()
