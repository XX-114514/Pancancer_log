"""Freeze all six same-case bundles after complete native-RNA acceptance."""
import itertools
import json
import os
import re
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from preflight_NSCLC_native_MFA_clinical_coverage_v1 import COMMON
from NSCLC_immune_RNA_MCP_cohort_math_v1 import design


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,16*1024**2)
    source=root/'NSCLC_immune_RNA_MCP_native_cells_20261019_v1'
    audit=root/'NSCLC_immune_RNA_MCP_native_cells_audit_20261019_v1/private_terminal.json'
    launchpath=root/'HNSCC_receiver_followup_job_control_20261016_v2/audit_NSCLC_immune_RNA_MCP_native_cells_v1/launch_receipt.json'
    launch=json.loads(launchpath.read_text());assert launch['pid']==23105 and str(launch['proc_start_ticks'])=='314948353'
    assert launch['source_sha256']=='15f8251af7743e7ebda6e6be370729b1813986171629407a90ce6eb4276d0507'
    clinicalroot=root/'NSCLC_immune_RNA_MCP_clinical_scope_20261019_v1'
    assert sha(root/'NSCLC_immune_RNA_MCP_clinical_scope_audit_20261019_v1/private_terminal.json')=='78d83326e28986f614de71d3fe3883424e7a3fc850642673f119c9a9af0d78d4'
    out=root/'NSCLC_immune_RNA_MCP_cohort_inputs_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    views=['R_FOXP3_capture','E_effector_capture','Myeloid']
    slots=[dict(slot=i,program=program,view_A=a,view_B=b) for i,(program,(a,b)) in enumerate(itertools.product(['MCP1','MCP2'],itertools.combinations(views,2)))]
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),math_sha256=sha(Path(__file__).with_name('NSCLC_immune_RNA_MCP_cohort_math_v1.py')),resource_guard=guard,
        native_auditor_launch_sha256=sha(launchpath),slots=slots,
        selection='All243 author materials, strict explicitICI and complete original8 clinical controls/pathology; both RNA views>=20; then drug-identity categories>=10 within same pair. No programme/outcome-effect selection.',
        agent_identity='Only remove explicitly deposited (Cnumber) or (Cnumber-Cnumber) suffix; keep different identities/combination sequences. Original cycles remain separately adjusted; Yes stays unspecified-agent recordedICI category.',
        scope='6 fixed bundles x(coordination,clinical_A,clinical_B), all same author cases per bundle. No actual joint effects computed in this input stage.',
        primary='IUT=max(Pcoord,PA,PB), all6 BH/BY and simultaneous18 MC upper P; fixed19999 HC3 wild,2000 common-stratum casebootstrap,>=1900 valid for CI; seed2026101901.',
        calibration='Before true effects, each eligible design256 normal/t3/age-heteroskedastic null simulations,999 HC3 wild each,all null rejection rates<=.10;seed2026101902. Working-model diagnostic, not universal guarantee.',
        scale='Per endpoint frozen original38-discovery case SD(ddof0). Unknown source constants absorbed by intercept; no cross-view min/product score.',
        controls='Original8 clinical numeric + drug identity full dummy span; clinical own4quality, coordination both4quality with rankednumeric control. Full SVD span, no automatic covariate deletion.',
        gates='>=50 same-case units,>=10 each pathology group, full nuisance span stable/machine-rank agreement, identifiable target,df>=10,condition<1e8,maxleverage<.999999.',
        point_sensitivities=['unadjusted','clinical_background_only','author_subtype_proportions','all_captured_mean','coordination_plus_response','LUAD_only','LUSC_only'],
        clinical_independent_patient_or_direct_interaction_claim=False,no_new_download=True)
    (out/'plan_before_any_source_coordination_or_clinical_effect.json').write_text(json.dumps(plan,indent=2)+'\n')
    deadline=time.monotonic()+21600
    while not audit.exists():
        proc=Path('/proc')/str(launch['pid'])/'stat';assert proc.exists(),'ORIGINAL_NATIVE_AUDITOR_MISSING_NO_RESTART'
        s=proc.read_text().rsplit(')',1)[1].split();assert s[19]==str(launch['proc_start_ticks']) and s[0]!='Z'
        assert time.monotonic()<deadline
        (out/'private_progress.json').write_text(json.dumps(dict(stage='WAITING_SAME_LIVE_COMPLETE_NATIVE_AUDITOR',pid=launch['pid'],startticks=s[19],observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n');time.sleep(15)
    proof=json.loads(audit.read_text());assert proof['status']=='ACCEPTED_ALL_RAW_NATIVE_COUNTS_RNA_GATES_HASH_SELECTIONS_AND_FULL_FROZEN_PROGRAMMES'
    assert proof['source_sha256']==launch['source_sha256'] and proof['producer_terminal_sha256']==sha(source/'private_terminal.json')
    for name,item in json.loads((source/'private_terminal.json').read_text())['files'].items():assert sha(source/name)==item['sha256']
    for name,item in json.loads((clinicalroot/'private_terminal.json').read_text())['files'].items():assert sha(clinicalroot/name)==item['sha256']
    clinical=pd.read_csv(clinicalroot/'private_all243_current_RNA_clinical_scope.tsv',sep='\t',keep_default_na=False).set_index('sampleID',drop=False)
    assert len(clinical)==243 and clinical.index.is_unique
    clinical['agent_identity']=clinical['anti-PD1_therapy'].map(lambda s:re.sub(r'\(C\d+(?:-C\d+)?\)','',str(s)).strip())
    clinical[['sampleID','anti-PD1_therapy','agent_identity','cycles']].to_csv(out/'private_all243_deposited_agent_identity_bridge.tsv',sep='\t',index=False)
    cases=pd.read_csv(source/'private_all_material_programme_scores.tsv',sep='\t')
    cells=pd.read_csv(source/'private_selected_native_cells.tsv.gz',sep='\t')
    native=root/'NSCLC_full_native_gene_units_20261016_v1'
    unit=pd.read_csv(native/'private_all_native_units_and_quality.tsv',sep='\t').set_index('unit_index')
    mapping=np.load(native/'private_cell_unit_index.npy',mmap_mode='r');cells['author_subtype']=unit.sub_cell_type.loc[mapping[cells.native_row.to_numpy(int)]].to_numpy()
    fractions={v:pd.crosstab(cells.loc[cells.view.eq(v),'sampleID'],cells.loc[cells.view.eq(v),'author_subtype'],normalize='index') for v in views}
    projection=root/'immune_RNA_MCP_frozen_projection_20261018_v1'
    terminal=json.loads((projection/'private_terminal.json').read_text())
    name='private_R_projection/private_all_case_mean_scores.tsv';assert sha(projection/name)==terminal['files'][name]
    old=pd.read_csv(projection/name,sep='\t');old=old[old.fit_id.eq('K2_seed20261018_source_adjusted')&old.partition.eq('discovery')]
    scales={}
    for (v,p),g in old.groupby(['view','program']):
        assert len(g)==38 and g['sample'].is_unique;scales[(v,p)]=float(g.score.std(ddof=0));assert scales[(v,p)]>0
    pd.DataFrame([dict(view=v,program=p,SD=sd,discovery_cases=38) for (v,p),sd in scales.items()]).to_csv(out/'private_all6_frozen_discovery_case_scales.tsv',sep='\t',index=False)
    scope=[];exclusions=[]
    for item in slots:
        a=cases[cases.view.eq(item['view_A'])&cases.program.eq(item['program'])].set_index('sampleID')
        b=cases[cases.view.eq(item['view_B'])&cases.program.eq(item['program'])].set_index('sampleID')
        initial=clinical.index[clinical.prospective_RNA_clinical_eligible.astype(str).str.lower().eq('true')]
        pair=sorted(set(initial)&set(a.index)&set(b.index))
        counts=clinical.loc[pair,'agent_identity'].value_counts();keep=set(counts[counts>=10].index)
        selected=[s for s in pair if clinical.loc[s,'agent_identity'] in keep]
        for s in clinical.index:
            reason='ELIGIBLE' if s in selected else 'INELIGIBLE_ORIGINAL_CLINICAL_ROLE_CONTROLS_OUTCOME' if s not in initial else 'ONE_OR_BOTH_RNA_VIEWS_UNCAPTURED' if s not in pair else 'DRUG_IDENTITY_GROUP_FEWER_THAN10'
            exclusions.append(dict(slot=item['slot'],sampleID=s,reason=reason))
        c=clinical.loc[selected];A=a.loc[selected];B=b.loc[selected]
        response=c.response.to_numpy(float);n0=int((response==0).sum());n1=int((response==1).sum())
        agent=pd.get_dummies(c.agent_identity,dtype=float);agent_names=['agent::'+str(s) for s in agent.columns]
        clinicalC=c[COMMON].to_numpy(float)
        qa=A[['cellQ','geneQ','fraction','logN']].to_numpy(float);qb=B[['cellQ','geneQ','fraction','logN']].to_numpy(float)
        fa=fractions[item['view_A']].reindex(selected,fill_value=0);fb=fractions[item['view_B']].reindex(selected,fill_value=0)
        strata=(c.cancer_type.astype(str)+'|'+c.agent_identity).to_numpy(str)
        for component in ['coordination','clinical_A','clinical_B']:
            both=component=='coordination';side='A' if component=='clinical_A' else 'B'
            own=qa if side=='A' else qb
            quality=np.c_[qa,qb] if both else own
            qnames=['quality::'+v+'::'+f for v in ([item['view_A'],item['view_B']] if both else [item['view_A'] if side=='A' else item['view_B']]) for f in ['cellQ','geneQ','fraction','logN']]
            C=np.c_[clinicalC,quality,agent.to_numpy()];names=np.array(COMMON+qnames+agent_names)
            fine=np.c_[fa.to_numpy(),fb.to_numpy()] if both else fa.to_numpy() if side=='A' else fb.to_numpy()
            fine_names=['fine::'+v+'::'+str(label) for v,f in ([(item['view_A'],fa),(item['view_B'],fb)] if both else [(item['view_A'],fa)] if side=='A' else [(item['view_B'],fb)]) for label in f.columns]
            fineC=np.c_[C,fine];fine_names=np.array(list(names)+fine_names)
            if both:x=A.score.to_numpy();y=B.score.to_numpy();allx=A.all_captured_mean_sensitivity.to_numpy();ally=B.all_captured_mean_sensitivity.to_numpy();kind='coordination'
            else:
                endpoint=A if side=='A' else B;v=item['view_A'] if side=='A' else item['view_B'];sd=scales[(v,item['program'])]
                x=response;y=endpoint.score.to_numpy()/sd;allx=response;ally=endpoint.all_captured_mean_sensitivity.to_numpy()/sd;kind='clinical'
            model,status=design(x,y,C,names,kind)
            if min(n0,n1)<10:model=None;status='FEWER_THAN10_IN_ONE_PATHOLOGY_GROUP'
            rec=dict(**item,component=component,kind=kind,n_cases=len(selected),response0=n0,response1=n1,drug_identity_groups=len(keep),eligible_before_drug_group=len(pair),status=status,
                nuisance_rank=model['nuisance_rank'] if model else -1,residual_df=model['df'] if model else -1,max_leverage=model['max_leverage'] if model else np.nan,condition=model['condition'] if model else np.nan)
            scope.append(rec)
            np.savez_compressed(out/('private_slot%02d_%s_inputs.npz'%(item['slot'],component)),x=x,y=y,C=C,names=names,fineC=fineC,fine_names=fine_names,
                cases=np.asarray(selected,str),strata=strata,cancer=c.cancer_type.to_numpy(str),response=response,allcapture_x=allx,allcapture_y=ally,kind=np.asarray(kind))
    pd.DataFrame(scope).to_csv(out/'private_all18_design_scope.tsv',sep='\t',index=False)
    pd.DataFrame(exclusions).to_csv(out/'private_all1458_slot_case_exclusions.tsv',sep='\t',index=False)
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(out.iterdir()) if p.is_file()};assert sum(v['bytes'] for v in files.values())<16*1024**2
    t=dict(status='COMPLETE_ALL6_SAME_CASE_BUNDLES_AND18_FULL_DESIGNS_INPUT_AUDIT_PENDING',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),native_acceptance_sha256=sha(audit),files=files,new_actual_effects_computed=False)
    (out/'private_terminal.json').write_text(json.dumps(t,indent=2)+'\n')
    print(json.dumps(dict(status=t['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
