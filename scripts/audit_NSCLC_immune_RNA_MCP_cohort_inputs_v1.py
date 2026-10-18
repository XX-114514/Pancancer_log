"""Rebuild all same-case cohort inputs before association or calibration."""
import itertools
import json
import os
import re
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from preflight_NSCLC_native_MFA_clinical_coverage_v1 import COMMON
from NSCLC_immune_RNA_MCP_cohort_audit_math_v1 import prepare


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,4*1024**2)
    source=root/'NSCLC_immune_RNA_MCP_cohort_inputs_20261019_v1'
    launchpath=root/'HNSCC_receiver_followup_job_control_20261016_v2/prepare_NSCLC_immune_RNA_MCP_cohort_inputs_v1/launch_receipt.json'
    launch=json.loads(launchpath.read_text());assert launch['source_sha256']=='5a21cf0c337b335c57156e22842490c18aa861af8b0f524a6789b64020d6562b'
    out=root/'NSCLC_immune_RNA_MCP_cohort_inputs_audit_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_independent_all_input_replay.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,
        scope='All6 pair cohorts,18 models/complete clinical and RNA control values,all agent normalizations/exclusions,original training scales and all fine-label frequencies. Independent dictionaries/counters and pandas ranks/pivotedQR; no target effects.'),indent=2)+'\n')
    deadline=time.monotonic()+21600
    while not (source/'private_terminal.json').exists():
        p=Path('/proc')/str(launch['pid'])/'stat';assert p.exists(),'ORIGINAL_INPUT_PRODUCER_MISSING_NO_RESTART'
        s=p.read_text().rsplit(')',1)[1].split();assert s[19]==str(launch['proc_start_ticks']) and s[0]!='Z';assert time.monotonic()<deadline
        (out/'private_progress.json').write_text(json.dumps(dict(stage='WAITING_SAME_LIVE_INPUT_PRODUCER',pid=launch['pid'],startticks=s[19],observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n');time.sleep(15)
    terminal=json.loads((source/'private_terminal.json').read_text());assert terminal['source_sha256']==launch['source_sha256']
    for name,item in terminal['files'].items():assert sha(source/name)==item['sha256']
    clinical=pd.read_csv(root/'NSCLC_immune_RNA_MCP_clinical_scope_20261019_v1/private_all243_current_RNA_clinical_scope.tsv',sep='\t',keep_default_na=False).set_index('sampleID')
    agent={case:re.sub(r'\(C\d+(?:-C\d+)?\)','',str(s)).strip() for case,s in clinical['anti-PD1_therapy'].items()}
    bridge=pd.read_csv(source/'private_all243_deposited_agent_identity_bridge.tsv',sep='\t',keep_default_na=False).set_index('sampleID')
    for case in clinical.index:assert bridge.loc[case,'agent_identity']==agent[case] and bridge.loc[case,'anti-PD1_therapy']==clinical.loc[case,'anti-PD1_therapy']
    native=root/'NSCLC_immune_RNA_MCP_native_cells_20261019_v1';scores=pd.read_csv(native/'private_all_material_programme_scores.tsv',sep='\t').set_index(['sampleID','view','program'])
    capture=pd.read_csv(native/'private_all729_material_view_capture.tsv',sep='\t').set_index(['sampleID','view'])
    cells=pd.read_csv(native/'private_selected_native_cells.tsv.gz',sep='\t');units=pd.read_csv(root/'NSCLC_full_native_gene_units_20261016_v1/private_all_native_units_and_quality.tsv',sep='\t').set_index('unit_index')
    um=np.load(root/'NSCLC_full_native_gene_units_20261016_v1/private_cell_unit_index.npy',mmap_mode='r');labels=units.sub_cell_type.loc[um[cells.native_row.to_numpy(int)]].to_numpy(str)
    counts=Counter((case,view,label) for case,view,label in zip(cells.sampleID,cells.view,labels));viewlabels={v:sorted(set(labels[cells.view.eq(v)])) for v in cells.view.unique()}
    sizes=Counter(zip(cells.sampleID,cells.view))
    old=pd.read_csv(root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection/private_all_case_mean_scores.tsv',sep='\t')
    old=old[old.fit_id.eq('K2_seed20261018_source_adjusted')&old.partition.eq('discovery')]
    scales={}
    for (view,program),g in old.groupby(['view','program']):
        assert len(g)==38;v=g.score.to_numpy(float);scales[(view,program)]=float(np.sqrt(np.sum((v-v.mean())**2)/38))
    savedscales=pd.read_csv(source/'private_all6_frozen_discovery_case_scales.tsv',sep='\t').set_index(['view','program'])
    for key,value in scales.items():assert abs(savedscales.loc[key,'SD']-value)<1e-12
    scope=pd.read_csv(source/'private_all18_design_scope.tsv',sep='\t').set_index(['slot','component']);ex=pd.read_csv(source/'private_all1458_slot_case_exclusions.tsv',sep='\t').set_index(['slot','sampleID'])
    assert len(scope)==18 and len(ex)==1458
    views=['R_FOXP3_capture','E_effector_capture','Myeloid'];slots=list(itertools.product(['MCP1','MCP2'],itertools.combinations(views,2)));rows_checked=0
    eligible=set(clinical.index[clinical.prospective_RNA_clinical_eligible.astype(str).str.lower().eq('true')])
    for slot,(program,(va,vb)) in enumerate(slots):
        pre=[case for case in sorted(eligible) if capture.loc[(case,va),'full_view_cells']>=20 and capture.loc[(case,vb),'full_view_cells']>=20]
        frequency=Counter(agent[c] for c in pre);selected=[c for c in pre if frequency[agent[c]]>=10]
        for case in clinical.index:
            expected='ELIGIBLE' if case in selected else 'INELIGIBLE_ORIGINAL_CLINICAL_ROLE_CONTROLS_OUTCOME' if case not in eligible else 'ONE_OR_BOTH_RNA_VIEWS_UNCAPTURED' if case not in pre else 'DRUG_IDENTITY_GROUP_FEWER_THAN10'
            assert ex.loc[(slot,case),'reason']==expected
        c=clinical.loc[selected];agent_names=sorted({agent[case] for case in selected});ac=np.array([[int(agent[case]==a) for a in agent_names] for case in selected],float).reshape(len(selected),len(agent_names))
        strata=np.array([clinical.loc[case,'cancer_type']+'|'+agent[case] for case in selected]);resp=c.response.to_numpy(float)
        qa=np.array([[scores.loc[(case,va,program),f] for f in ['cellQ','geneQ','fraction','logN']] for case in selected]).reshape(len(selected),4);qb=np.array([[scores.loc[(case,vb,program),f] for f in ['cellQ','geneQ','fraction','logN']] for case in selected]).reshape(len(selected),4)
        f_a=np.array([[counts[(case,va,label)]/sizes[(case,va)] for label in viewlabels[va]] for case in selected]).reshape(len(selected),len(viewlabels[va]));f_b=np.array([[counts[(case,vb,label)]/sizes[(case,vb)] for label in viewlabels[vb]] for case in selected]).reshape(len(selected),len(viewlabels[vb]))
        for comp in ['coordination','clinical_A','clinical_B']:
            both=comp=='coordination';a=comp=='clinical_A';q=np.c_[qa,qb] if both else qa if a else qb
            qv=[va,vb] if both else [va] if a else [vb]
            names=COMMON+['quality::'+v+'::'+f for v in qv for f in ['cellQ','geneQ','fraction','logN']]+['agent::'+s for s in agent_names]
            C=np.c_[c[COMMON].to_numpy(float),q,ac]
            fv=[(va,f_a),(vb,f_b)] if both else [(va,f_a)] if a else [(vb,f_b)]
            fine=np.c_[C,*[f for _,f in fv]];fn=names+['fine::'+v+'::'+label for v,_ in fv for label in viewlabels[v]]
            yview=vb if both or not a else va
            if both:x=np.array([scores.loc[(case,va,program),'score'] for case in selected]);kind='coordination'
            else:x=resp;kind='clinical'
            y=np.array([scores.loc[(case,yview,program),'score'] for case in selected]);all_y=np.array([scores.loc[(case,yview,program),'all_captured_mean_sensitivity'] for case in selected])
            if not both:y/=scales[(yview,program)];all_y/=scales[(yview,program)]
            all_x=np.array([scores.loc[(case,va,program),'all_captured_mean_sensitivity'] for case in selected]) if both else resp
            z=np.load(source/('private_slot%02d_%s_inputs.npz'%(slot,comp)),allow_pickle=False)
            assert z['cases'].tolist()==selected and z['names'].tolist()==names and z['fine_names'].tolist()==fn and z['strata'].tolist()==strata.tolist()
            for name,value in [('x',x),('y',y),('C',C),('fineC',fine),('allcapture_x',all_x),('allcapture_y',all_y),('response',resp)]:assert np.allclose(z[name],value,atol=1e-12,rtol=0)
            assert z['cancer'].tolist()==c.cancer_type.tolist() and str(z['kind'])==kind
            m,status=prepare(x,y,C,names,kind)
            if min(np.count_nonzero(resp==0),np.count_nonzero(resp==1))<10:m=None;status='FEWER_THAN10_IN_ONE_PATHOLOGY_GROUP'
            r=scope.loc[(slot,comp)];assert r.status==status and r.n_cases==len(selected) and r.response0==np.count_nonzero(resp==0) and r.response1==np.count_nonzero(resp==1)
            if m:
                assert r.nuisance_rank==m['nuisance_rank'] and r.residual_df==m['df'];assert abs(r.condition-m['condition'])<1e-7 and abs(r.max_leverage-m['max_leverage'])<1e-9
            rows_checked+=len(selected)
    accepted=dict(status='ACCEPTED_ALL6_BUNDLES_18_INPUT_DESIGNS_AND_ALL_EXCLUSIONS',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),producer_terminal_sha256=sha(source/'private_terminal.json'),
        models=18,pair_bundles=6,exclusion_positions=1458,complete_model_case_rows_checked=rows_checked,real_association_effects_computed=False)
    (out/'private_terminal.json').write_text(json.dumps(accepted,indent=2)+'\n');print(json.dumps(dict(status=accepted['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
