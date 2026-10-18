"""Rebuild every paired vector and node prediction with dictionary/Counter logic."""
import json
import os
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from HNSCC_immune_RNA_MCP_paired_audit_math_v1 import reference_model


def close(a,b):assert np.allclose(a,b,rtol=0,atol=1e-9,equal_nan=True)


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,4*1024**2)
    source=root/'HNSCC_immune_RNA_MCP_paired_inputs_20261019_v1';launchpath=root/'HNSCC_receiver_followup_job_control_20261016_v2/prepare_HNSCC_immune_RNA_MCP_paired_inputs_v1/launch_receipt.json'
    launch=json.loads(launchpath.read_text());assert launch['source_sha256']=='950be6b420d69ebf2adb2deb4e886ea07fcec0bd0957de28ab18fa5e70d6c3e8'
    out=root/'HNSCC_immune_RNA_MCP_paired_input_audit_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_all_paired_input_replay.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,producer_launch_sha256=sha(launchpath),scope='All12 exact paired vectors/quality/complete-case masks, all312 original-node predictions and all6 full rank designs; independent author static-field dictionary and node-count Counter, pandas/QR reference. No real effects/P.'),indent=2)+'\n')
    deadline=time.monotonic()+21600
    while not (source/'private_terminal.json').exists():
        p=Path('/proc')/str(launch['pid'])/'stat';assert p.exists(),'PRODUCER_MISSING_WITHOUT_TERMINAL'
        st=p.read_text().rsplit(')',1)[1].split();assert st[19]==str(launch['proc_start_ticks']) and st[0]!='Z';assert time.monotonic()<deadline
        (out/'private_progress.json').write_text(json.dumps(dict(phase='WAITING_SAME_PAIRED_INPUT_PRODUCER',observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n');time.sleep(15)
    t=json.loads((source/'private_terminal.json').read_text());assert t['status']=='COMPLETE_ALL12_PAIRED_INPUTS_FULL_DESIGNS_AND_COMPOSITION_SENSITIVITY_NO_EFFECTS'
    for n,v in t['files'].items():assert sha(source/n)==v['sha256']
    contract=root/'HNSCC_immune_RNA_MCP_paired_dynamics_contract_20261019_v1';planpath=contract/'plan_before_any_paired_difference_or_effect.json';assert sha(planpath)=='8c03ca8786dcf0de2a872e45cbd527b4eecb3181cc39c71297946a024dd91e0a';plan=json.loads(planpath.read_text())
    assert sha(root/'HNSCC_immune_RNA_MCP_paired_native_measurement_audit_20261019_v1/private_terminal.json')=='f3c916ecb8f4d40b535ca0c4486a0ae8a825dbff1d362ae6b64e721b5c981640'
    meta=pd.read_csv(contract/'private_all26_static_identity_design_controls.tsv',sep='\t');patients=meta.patient_registry.to_numpy(str);assert len(patients)==len(set(patients))==26
    author=pd.read_csv(root/'HNSCC_author_clinical_table_user_20261014_v1/private_key_and_response_audit_v1/all57_private_author_table_join.tsv',sep='\t',keep_default_na=False,usecols=['patient_registry','library_id','timepoint_registry','neo_type','HPV/p16'])
    static=[]
    for r in meta.itertuples(index=False):
        a=author[author.patient_registry.eq(r.patient_registry)];assert a.neo_type.nunique()==a['HPV/p16'].nunique()==1
        pre=a[a.timepoint_registry.eq('pre')];post=a[a.timepoint_registry.eq('post')];assert len(pre)==len(post)==1
        assert r.pre_sampleID=='GSE286827::'+pre.library_id.iloc[0] and r.post_sampleID=='GSE286827::'+post.library_id.iloc[0]
        treatment=pre.neo_type.iloc[0];hpv=pre['HPV/p16'].iloc[0];values=[float(treatment=='D+T') if treatment in ['D','D+T'] else np.nan,float(hpv=='positive') if hpv in ['negative','positive'] else np.nan];static.append(values)
        close([r.treatment_DplusT,r.HPV_positive],values);assert r.complete_static_controls==bool(np.isfinite(values).all())
    static=np.array(static)
    assay=root/'HNSCC_immune_RNA_MCP_paired_native_measurement_20261019_v1';scores=pd.read_csv(assay/'private_all624_material_programme_measurements.tsv',sep='\t')
    smap={(r.patient_registry,r.phase,r.view,r.program,r.variant):r for r in scores.itertuples(index=False)};assert len(smap)==624
    scales={(r.view,r.program):r.discovery_case_sd for r in pd.read_csv(assay/'private_all6_discovery_case_scales.tsv',sep='\t').itertuples(index=False)}
    fields=['cellQ','geneQ','fraction','logN'];views=['R_FOXP3_capture','E_effector_capture','Myeloid']
    old=root/'existing_immune_fine_native_measurement_scope_20261016_v1';unit=pd.read_csv(old/'private_all_frozen_T_myeloid_node_units.tsv',sep='\t',dtype={'node_id':str}).set_index('unit_index').node_id.to_dict();mapping=np.load(old/'private_original_cell_to_fine_node_unit.npy',mmap_mode='r')
    counts={}
    for r in pd.read_csv(assay/'private_selected_cells_and_all_scores.tsv.gz',sep='\t').itertuples(index=False):counts.setdefault((int(r.material_index),views[int(r.view_index)]),Counter())[unit[int(mapping[int(r.global_row)])]]+=1
    coefficients=pd.read_csv(root/'immune_RNA_MCP_original_node_decomposition_20261018_v1/private_frozen_discovery_node_coefficients.tsv',sep='\t',dtype={'node_id':str});coef={(r.view,r.program,r.node_id):r.coefficient for r in coefficients.itertuples(index=False)}
    mats=pd.read_csv(assay/'private_all52_paired_materials.tsv',sep='\t');pred={}
    for m,row in enumerate(mats.itertuples(index=False)):
        for view in views:
            c=counts[(m,view)];N=sum(c.values())
            for p in ['MCP1','MCP2']:
                unseen=sum(n for node,n in c.items() if (view,p,node) not in coef);value=sum(n*coef[(view,p,node)] for node,n in c.items())/N if unseen==0 else np.nan;pred[(row.patient_registry,row.phase,view,p)]=(value,unseen,N)
    savedcomp=pd.read_csv(source/'private_all312_frozen_node_predictions.tsv',sep='\t');assert len(pred)==len(savedcomp)==312
    for r in savedcomp.itertuples(index=False):
        val,unseen,N=pred[(r.patient_registry,r.phase,r.view,r.program)];close(r.composition_prediction,val);assert r.unseen_selected_cells==unseen and r.selected_cells==N
    def vector(v,p,variant='full',field='score'):
        return np.array([getattr(smap[(patient,'post',v,p,variant)],field)-getattr(smap[(patient,'pre',v,p,variant)],field) for patient in patients])
    scope=pd.read_csv(source/'private_all12_input_design_scope.tsv',sep='\t').set_index('slot');assert len(scope)==12;slots=0;models=0
    for item in plan['slots']:
        i=item['slot'];p=item['program']
        if item['kind']=='paired_median_shift':
            v=item['view'];sc=scales[(v,p)];d=vector(v,p)/sc;z=np.load(source/('private_slot%02d_shift_inputs.npz'%i));assert np.array_equal(z['patients'],patients)
            for phase in ['pre','post']:close(z[phase],np.array([smap[(patient,phase,v,p,'full')].score/sc for patient in patients]))
            close(z['delta'],d);close(z['scale'],sc);close(z['all_captured_delta'],vector(v,p,field='all_captured_mean_sensitivity')/sc);close(z['fixed_stress_delta'],vector(v,p,variant='fixed_stress_complement')/sc)
            nonzero=int((np.abs(d)>1e-10).sum());assert scope.loc[i,'nonzero_pairs']==nonzero and scope.loc[i,'n_patients']==26
            status='ESTIMABLE' if len(d)>=20 and nonzero>=10 else 'INSUFFICIENT_PAIRED_OR_NONZERO_COVERAGE';assert scope.loc[i,'status']==status
        else:
            a=item['view_A'];b=item['view_B'];x=vector(a,p)/scales[(a,p)];y=vector(b,p)/scales[(b,p)];qa=np.column_stack([vector(a,p,field=f) for f in fields]);qb=np.column_stack([vector(b,p,field=f) for f in fields]);C=np.column_stack([static,qa,qb]);good=np.isfinite(np.column_stack([x,y,C])).all(axis=1)
            z=np.load(source/('private_slot%02d_cochange_inputs.npz'%i));assert np.array_equal(z['patients'],patients[good]);close(z['x'],x[good]);close(z['y'],y[good]);close(z['C'],C[good]);assert z['names'].tolist()==plan['cochange']['controls']
            for label,variant,field in [('all_captured','full','all_captured_mean_sensitivity'),('fixed_stress','fixed_stress_complement','score')]:
                close(z[label+'_x'],(vector(a,p,variant=variant,field=field)/scales[(a,p)])[good]);close(z[label+'_y'],(vector(b,p,variant=variant,field=field)/scales[(b,p)])[good])
            cp=[];un=[]
            for view in [a,b]:
                cp.append(np.array([pred[(patient,'post',view,p)][0]-pred[(patient,'pre',view,p)][0] for patient in patients]));un.append(np.array([pred[(patient,'post',view,p)][1]+pred[(patient,'pre',view,p)][1] for patient in patients]))
            close(z['C_plus_composition'],np.column_stack([C,*cp])[good]);assert np.array_equal(z['composition_unseen'],np.column_stack(un)[good])
            m,status=reference_model(x[good],y[good],C[good],z['names']);assert scope.loc[i,'status']==status and scope.loc[i,'n_patients']==int(good.sum())
            if m is not None:
                assert scope.loc[i,'nuisance_rank']==m['nuisance_rank'] and scope.loc[i,'residual_df']==m['residual_df'];close(scope.loc[i,'condition'],m['condition']);models+=1
        slots+=1
    t=dict(status='ACCEPTED_ALL12_PAIRED_INPUTS_ALL312_NODE_PREDICTIONS_AND_ALL_COMPLETE_DESIGNS',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),producer_terminal_sha256=sha(source/'private_terminal.json'),slots_replayed=slots,cochange_designs_estimable=models,all26_static_controls_independently_replayed=True,no_true_effect_P_or_CI_computed=True)
    (out/'private_terminal.json').write_text(json.dumps(t,indent=2)+'\n');print(json.dumps(dict(status=t['status'],terminal_sha256=sha(out/'private_terminal.json'),estimable_cochange_designs=models)))


if __name__=='__main__':main()
