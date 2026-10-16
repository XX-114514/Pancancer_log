"""Independent full-scope native inputs and geometry audit, without effects."""
import hashlib,json,math,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.linalg import qr,svdvals
from scipy.stats import rankdata

COMMON=['age_numeric','male','smoker','LUSC','stage_ordinal','cycles_numeric','chemo_binary','targeted_binary']
UNIT=['fraction','logN','mean_logRNA','mean_loggenes']

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()

def same(a,b):assert np.allclose(a,b,rtol=1e-9,atol=1e-9,equal_nan=True),(np.shape(a),np.shape(b))

def rank(A):
    s=svdvals(A);return int(np.sum(s>max(A.shape)*np.finfo(float).eps*s[0])) if len(s) else 0

def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);O=R/'NSCLC_automatic_ligand_response_preflight_20261016_v1'
    launch=json.loads((R/'HNSCC_receiver_followup_job_control_20261016_v2/preflight_NSCLC_automatic_ligand_response_v1/launch_receipt.json').read_text());deadline=time.monotonic()+21600
    while not (O/'private_terminal.json').exists():
        p=Path('/proc')/str(launch['pid'])/'stat';assert p.exists();st=p.read_text().rsplit(')',1)[1].split()
        assert st[19]==str(launch['proc_start_ticks']) and st[0]!='Z' and time.monotonic()<deadline;time.sleep(15)
    result_path=O/'private_independent_audit.json';assert not result_path.exists()
    terminal=json.loads((O/'private_terminal.json').read_text());plan=json.loads((O/'plan_before_any_ligand_response_effect.json').read_text())
    lpath=O/'private_all45760_logical_slots.tsv.gz';tpath=O/'private_all_unique_ligand_response_designs.tsv'
    assert sha(lpath)==terminal['logical_scope_sha256'] and sha(tpath)==terminal['unique_scope_sha256']
    logical=pd.read_csv(lpath,sep='\t',keep_default_na=False);table=pd.read_csv(tpath,sep='\t');assert len(logical)==45760 and table.record.is_unique
    P=R/'NSCLC_all_frozen_receiver_prior_ranking_20261016_v1';ppath=P/'private_all960_prospective_prior_slots.tsv';assert sha(ppath)==plan['input_sha256']['prior_slots']
    prior=pd.read_csv(ppath,sep='\t',keep_default_na=False)
    spath=R/'NSCLC_two_end_coordination_preflight_20261016_v1/private_all572_planned_pair_designs.tsv';assert sha(spath)==plan['input_sha256']['original_pair_scope']
    parents=pd.read_csv(spath,sep='\t').set_index('record');assert len(parents)==572
    expected=[]
    for name,p in parents.iterrows():
        for receiver in ['typeA','typeB']:
            sig=f'{p.context}_{receiver}_{p.program}';sub=prior[prior.signature.eq(sig)]
            assert len(sub)==40
            expected.extend((name,receiver,sig,r.programme_sign,int(r.prospective_position),r.ligand) for r in sub.itertuples())
    columns=['parent_pair','receiver_type','receiver_signature','programme_sign','prior_position','ligand']
    assert set(logical[columns].itertuples(index=False,name=None))==set(expected) and len(logical)==len(expected)
    nonempty=logical[logical.ligand.ne('')];assert set(nonempty.hypothesis)==set(table.record) and logical.loc[logical.ligand.eq(''),'hypothesis'].eq('').all()
    assert not table.duplicated(['parent_pair','receiver_type','ligand']).any()
    for _,f in nonempty.groupby(['parent_pair','receiver_type','ligand']):assert f.hypothesis.nunique()==1
    provenance={}
    for r in nonempty.itertuples():
        entry=provenance.setdefault(r.hypothesis,dict(ligands=set(),signs=set()));entry['ligands'].add(r.ligand);entry['signs'].add(r.programme_sign)
    I=R/'NSCLC_prior_joint_receptor_RNA_20261016_v1';audit=json.loads((I/'private_independent_audit.json').read_text())
    assert sha(I/'private_independent_audit.json')==plan['input_sha256']['joint_RNA_audit']
    epath=I/'private_all_frozen_interfaces_joint_RNA_qualification.tsv';assert sha(epath)==audit['qualification_sha256']
    edges=pd.read_csv(epath,sep='\t');edges=edges[edges.native_joint1pct_qualified]
    joint=pd.read_csv(I/'private_all_native_unit_receptor_joint_counts.tsv.gz',sep='\t').set_index(['receptor','unit_index'])
    N=R/'NSCLC_full_native_gene_units_20261016_v1';genes=pd.Index(pd.read_csv(N/'private_native_gene_axis.tsv',sep='\t').gene);M=np.load(N/'mean_logCP10K.npy',mmap_mode='r')
    coordplan=json.loads((R/'NSCLC_two_end_coordination_preflight_20261016_v1/plan_before_any_two_end_effect.json').read_text())
    param_path=R/'NSCLC_context_frozen_parameters_20261016_v1/private_all53_frozen_gene_parameters.tsv';assert sha(param_path)==coordplan['input_sha256']['parameters']
    params=pd.read_csv(param_path,sep='\t');S=R/'NSCLC_context_native_programme_projection_20261016_v1/private_native_author_subtype_programme_scores.tsv.gz'
    assert sha(S)==coordplan['input_sha256']['scores'];scores=pd.read_csv(S,sep='\t')
    cp=R/'NSCLC_frozen_signature_clinical_preflight_20261016_v1/private_all243_clinical_eligibility.tsv';assert sha(cp)==coordplan['input_sha256']['case_eligibility']
    c=pd.read_csv(cp,sep='\t').set_index('sampleID');groups={};masks={};design_hashes=[];checked=0;maxerr=0.
    for row in table.itertuples():
        key=(row.parent_pair,row.receiver_type)
        if key not in groups:
            par=parents.loc[row.parent_pair];rt=row.receiver_type;st='typeB' if rt=='typeA' else 'typeA'
            rn=par.author_subtype_A if rt=='typeA' else par.author_subtype_B;sn=par.author_subtype_B if st=='typeB' else par.author_subtype_A
            sf=scores[scores.context.eq(par.context)&scores.program.eq(par.program)&scores.type.eq(st)&scores.sub_cell_type.eq(sn)].set_index('sampleID')
            rf=scores[scores.context.eq(par.context)&scores.program.eq(par.program)&scores.type.eq(rt)&scores.sub_cell_type.eq(rn)].set_index('sampleID')
            cases=sorted(set(sf.index)&set(rf.index));sf=sf.loc[cases];rf=rf.loc[cases];cc=c.loc[cases]
            take=sf.n_cells.ge(20).to_numpy()&rf.n_cells.ge(20).to_numpy()&cc.clinical_eligible.to_numpy(bool)
            groups[key]=(par,st,rn,sn,sf.loc[take],rf.loc[take],cc.loc[take])
        par,st,rn,sn,sf,rf,cc=groups[key];n=len(rf)
        assert row.context==par.context and row.program==par.program and row.sender_author_subtype==sn and row.receiver_author_subtype==rn and row.matched_native_cases==n
        origin=provenance[row.record];assert origin['ligands']=={row.ligand} and set(str(row.prior_signs).split('|'))==origin['signs']
        er=edges[edges.signature.eq(row.receiver_signature)&edges.ligand.eq(row.ligand)];receptors=sorted(set(er.receptor))
        actual='' if pd.isna(row.qualified_receptors) else row.qualified_receptors;assert actual=='|'.join(receptors)
        status='READY_FOR_INDEPENDENT_DESIGN_AUDIT';D=None
        if row.ligand not in genes:status='LITERAL_LIGAND_NATIVE_GENE_NOT_MEASURED'
        elif not receptors:status='NO_FROZEN_REFERENCE_AND_JOINT_RNA_QUALIFIED_INTERFACE'
        else:
            mk=(row.context,row.program,row.receiver_type,row.ligand)
            if mk not in masks:
                p=params[params.context.eq(row.context)&params.program.eq(row.program)&params.type.eq(row.receiver_type)]
                s=params[params.context.eq(row.context)&params.program.eq(row.program)&params.type.eq(st)]
                remove=set(s.gene)|{row.ligand}|{g for r in receptors for g in r.split('_')};remaining=[r for r in p.itertuples() if r.gene not in remove]
                total=math.fsum(abs(r.weight) for r in p.itertuples());left=math.fsum(abs(r.weight) for r in remaining);present=[r for r in remaining if r.gene in genes]
                retained=left/total;measured=math.fsum(abs(r.weight) for r in present)/left if left else 0.
                valid=retained>=.8 and measured>=.95 and len(present)>=20;masks[mk]=(valid,retained,measured,present)
            valid,retained,measured,present=masks[mk];same([retained,measured],[row.target_retained_original_L1,row.target_native_remaining_L1])
            if not valid:status='RECEIVER_MASK_CHANGES_PROGRAMME_OR_NATIVE_COVERAGE_TOO_MUCH'
            elif n<50:status='INSUFFICIENT_MATCHED_CASES'
            else:
                si=sf.unit_index.to_numpy(int);ui=rf.unit_index.to_numpy(int);column=genes.get_loc(row.ligand)
                x=np.array(M[si,column]);own=np.array(M[ui,column]);y=np.zeros(n)
                for g in present:y+=(np.array(M[ui,genes.get_loc(g.gene)])-g.gene_mean)*g.weight/g.gene_sd
                fractions=[joint.loc[[(r,int(u)) for u in ui],'joint_positive_cells'].to_numpy(float)/rf.n_cells.to_numpy(float) for r in receptors]
                rec=np.maximum.reduce(fractions)
                cov=np.column_stack([cc[COMMON].to_numpy(float),sf[UNIT].to_numpy(float),rf[UNIT].to_numpy(float),own,rec]);sd=cov.std(0)
                Z=np.column_stack([np.ones(n),(cov-cov.mean(0))/np.where(sd>1e-12,sd,1)])
                xr=rankdata(x);xr=(xr-xr.mean())/xr.std() if xr.std()>0 else np.zeros(n);D=np.column_stack([Z,xr]);rnk=rank(D);ss=svdvals(D);condition=float(ss[0]/ss[-1])
                assert row.design_rank==rnk and row.residual_df==n-rnk and row.sender_RNA_positive_cases==int(np.sum(x>0))
                if rnk==20 and condition<1e8:same(condition,row.condition)
                if np.sum(x>0)<10:status='INSUFFICIENT_CASES_WITH_NATIVE_SENDER_LIGAND_RNA'
                elif (sd<=1e-12).any() or rank(Z)!=19 or rnk!=20:status='CONSTANT_OR_COLLINEAR_REQUIRED_DESIGN'
                elif n-rnk<10 or condition>=1e8:status='LOW_DF_OR_ILL_CONDITIONED'
                elif np.ptp(y)<=1e-12:status='NO_MASKED_RECEIVER_SIGNATURE_VARIATION'
                elif np.max(np.sum(qr(D,mode='economic')[0]**2,axis=1))>=.999999:status='EXTREME_LEVERAGE'
                if status=='READY_FOR_INDEPENDENT_DESIGN_AUDIT':
                    path=O/'private_frozen_designs'/(row.record+'.npz');z=np.load(path,allow_pickle=False)
                    assert z['cases'].tolist()==rf.index.tolist()==sf.index.tolist();assert np.array_equal(z['sender_units'],si) and np.array_equal(z['receiver_units'],ui)
                    same(z['X'],x);same(z['Y'],y);maxerr=max(maxerr,float(np.max(np.abs(z['Y']-y))));same(z['Z'],Z);same(z['covariates'],cov)
                    same(z['sender_programme'],sf.raw_weighted_score);same(z['receiver_full_programme'],rf.raw_weighted_score);same(z['response'],cc.response)
                    design_hashes.append(dict(record=row.record,sha256=sha(path)))
        assert row.status==status,(row.record,row.status,status)
        checked+=1
        if checked%250==0:(O/'private_audit_progress.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),checked_unique_hypotheses=checked,total=len(table))))
    assert {p.stem for p in (O/'private_frozen_designs').glob('*.npz')}=={x['record'] for x in design_hashes}
    hashes=O/'private_all_ready_design_SHA256.tsv';pd.DataFrame(design_hashes,columns=['record','sha256']).to_csv(hashes,sep='\t',index=False)
    result=dict(status='PASS_ALL45760_LOGICAL_SLOTS_ALL_UNIQUE_HYPOTHESES_NATIVE_INPUTS_MASKS_FULL_CONTROLS_AND_DESIGN_GATES_NO_EFFECTS',
        finished_UTC=datetime.now(timezone.utc).isoformat(),audit_code_sha256=sha(__file__),producer_terminal_sha256=sha(O/'private_terminal.json'),
        logical_slots=len(logical),unique_hypotheses=checked,ready_designs=len(design_hashes),maximum_receiver_gene_sum_error=maxerr,
        ready_design_hash_table_sha256=sha(hashes),unique_scope_sha256=sha(tpath),new_effects_or_P_BH=False)
    result_path.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))

if __name__=='__main__':main()
