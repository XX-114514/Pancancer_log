"""All frozen paired vectors/designs, no aggregate effects or hypothesis tests."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from HNSCC_immune_RNA_MCP_paired_math_v1 import design


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,16*1024**2)
    assay=root/'HNSCC_immune_RNA_MCP_paired_native_measurement_20261019_v1';accepted=root/'HNSCC_immune_RNA_MCP_paired_native_measurement_audit_20261019_v1/private_terminal.json'
    assert sha(accepted)=='f3c916ecb8f4d40b535ca0c4486a0ae8a825dbff1d362ae6b64e721b5c981640'
    at=json.loads(accepted.read_text());assert sha(assay/'private_terminal.json')==at['producer_terminal_sha256']
    for n,v in json.loads((assay/'private_terminal.json').read_text())['files'].items():assert sha(assay/n)==v['sha256']
    contract=root/'HNSCC_immune_RNA_MCP_paired_dynamics_contract_20261019_v1'
    planpath=contract/'plan_before_any_paired_difference_or_effect.json';assert sha(planpath)=='8c03ca8786dcf0de2a872e45cbd527b4eecb3181cc39c71297946a024dd91e0a'
    plan=json.loads(planpath.read_text());assert len(plan['slots'])==12
    for n,v in json.loads((contract/'private_terminal.json').read_text())['files'].items():assert sha(contract/n)==v['sha256']
    meta=pd.read_csv(contract/'private_all26_static_identity_design_controls.tsv',sep='\t');assert meta.patient_registry.is_unique and len(meta)==26
    patients=meta.patient_registry.to_numpy(str)
    out=root/'HNSCC_immune_RNA_MCP_paired_inputs_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_all_paired_input_vectors.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,assay_acceptance_sha256=sha(accepted),frozen_contract_sha256=sha(planpath),scope='All12 frozen input slots, all identities/deltas/quality/designs and fixed original-node prediction sensitivity. No true effects/P/CI or direction counts computed.'),indent=2)+'\n')
    scores=pd.read_csv(assay/'private_all624_material_programme_measurements.tsv',sep='\t');assert len(scores)==624
    scales=pd.read_csv(assay/'private_all6_discovery_case_scales.tsv',sep='\t').set_index(['view','program']).discovery_case_sd
    node=root/'immune_RNA_MCP_original_node_decomposition_20261018_v1';audit=root/'immune_RNA_MCP_original_node_decomposition_audit_20261018_v1/private_terminal.json';assert sha(audit)=='f109cecf33be6ad8a937bce8c50d567aa1426f1a1e1ea05a0c20335ce4e123d7'
    nt=json.loads((node/'private_terminal.json').read_text());filename='private_frozen_discovery_node_coefficients.tsv';assert sha(node/filename)==nt['files'][filename]
    coefficients=pd.read_csv(node/filename,sep='\t',dtype={'node_id':str})
    oldscope=root/'existing_immune_fine_native_measurement_scope_20261016_v1';units=pd.read_csv(oldscope/'private_all_frozen_T_myeloid_node_units.tsv',sep='\t',dtype={'node_id':str});mapping=np.load(oldscope/'private_original_cell_to_fine_node_unit.npy',mmap_mode='r')
    selected=pd.read_csv(assay/'private_selected_cells_and_all_scores.tsv.gz',sep='\t');assert (mapping[selected.global_row.to_numpy(int)]>=0).all()
    selected['node_id']=units.set_index('unit_index').node_id.loc[mapping[selected.global_row.to_numpy(int)]].to_numpy()
    materials=pd.read_csv(assay/'private_all52_paired_materials.tsv',sep='\t').set_index('sampleID')
    composition=[];views=['R_FOXP3_capture','E_effector_capture','Myeloid'];fields=['cellQ','geneQ','fraction','logN']
    for vi,view in enumerate(views):
        for program in ['MCP1','MCP2']:
            lookup=coefficients[coefficients.view.eq(view)&coefficients.program.eq(program)].set_index('node_id').coefficient
            for m,g in selected[selected.view_index.eq(vi)].groupby('material_index'):
                values=g.node_id.map(lookup);base=materials[materials.index.isin(scores[scores.material_index.eq(m)].sampleID.unique())];assert len(base)==1
                composition.append(dict(patient_registry=base.iloc[0].patient_registry,phase=base.iloc[0].phase,view=view,program=program,composition_prediction=float(values.mean()) if values.notna().all() else np.nan,unseen_selected_cells=int(values.isna().sum()),selected_cells=len(g)))
    comp=pd.DataFrame(composition);assert len(comp)==312;comp.to_csv(out/'private_all312_frozen_node_predictions.tsv',sep='\t',index=False)
    frame_cache={}
    def vectors(view,program):
        key=(view,program)
        if key in frame_cache:return frame_cache[key]
        f=scores[scores.view.eq(view)&scores.program.eq(program)];sc=float(scales.loc[key]);assert sc>0
        r={}
        for variant in ['full','fixed_stress_complement']:
            g=f[f.variant.eq(variant)].set_index(['patient_registry','phase']);assert len(g)==52 and g.index.is_unique
            pre=g.loc[pd.MultiIndex.from_arrays([patients,np.repeat('pre',26)])];post=g.loc[pd.MultiIndex.from_arrays([patients,np.repeat('post',26)])]
            assert np.array_equal(pre.sampleID,meta.pre_sampleID) and np.array_equal(post.sampleID,meta.post_sampleID)
            r[variant]=(post.score.to_numpy()-pre.score.to_numpy())/sc
            if variant=='full':
                r['pre']=pre.score.to_numpy()/sc;r['post']=post.score.to_numpy()/sc
                r['all_captured']=(post.all_captured_mean_sensitivity.to_numpy()-pre.all_captured_mean_sensitivity.to_numpy())/sc
                r['quality']=post[fields].to_numpy(float)-pre[fields].to_numpy(float)
        g=comp[comp.view.eq(view)&comp.program.eq(program)].set_index(['patient_registry','phase'])
        pre=g.loc[pd.MultiIndex.from_arrays([patients,np.repeat('pre',26)])];post=g.loc[pd.MultiIndex.from_arrays([patients,np.repeat('post',26)])]
        r['composition']=post.composition_prediction.to_numpy()-pre.composition_prediction.to_numpy();r['unseen']=post.unseen_selected_cells.to_numpy()+pre.unseen_selected_cells.to_numpy();r['scale']=sc;frame_cache[key]=r;return r
    scope_rows=[]
    for item in plan['slots']:
        slot=item['slot'];kind=item['kind']
        if kind=='paired_median_shift':
            a=vectors(item['view'],item['program']);delta=a['full'];nonzero=int((np.abs(delta)>1e-10).sum());status='ESTIMABLE' if len(delta)>=20 and nonzero>=10 else 'INSUFFICIENT_PAIRED_OR_NONZERO_COVERAGE'
            np.savez_compressed(out/('private_slot%02d_shift_inputs.npz'%slot),patients=patients,pre=a['pre'],post=a['post'],delta=delta,all_captured_delta=a['all_captured'],fixed_stress_delta=a['fixed_stress_complement'],scale=a['scale'])
            scope_rows.append(dict(slot=slot,kind=kind,program=item['program'],view=item['view'],n_patients=len(delta),status=status,nonzero_pairs=nonzero))
        else:
            a=vectors(item['view_A'],item['program']);b=vectors(item['view_B'],item['program']);C=np.c_[meta[['treatment_DplusT','HPV_positive']].to_numpy(float),a['quality'],b['quality']]
            names=plan['cochange']['controls'];complete=meta.complete_static_controls.to_numpy(bool)&np.isfinite(np.c_[a['full'],b['full'],C]).all(axis=1);x=a['full'][complete];y=b['full'][complete];cc=C[complete];m,status=design(x,y,cc,names)
            np.savez_compressed(out/('private_slot%02d_cochange_inputs.npz'%slot),patients=patients[complete],x=x,y=y,C=cc,names=np.asarray(names),all_captured_x=a['all_captured'][complete],all_captured_y=b['all_captured'][complete],fixed_stress_x=a['fixed_stress_complement'][complete],fixed_stress_y=b['fixed_stress_complement'][complete],C_plus_composition=np.c_[cc,a['composition'][complete],b['composition'][complete]],composition_unseen=np.c_[a['unseen'][complete],b['unseen'][complete]])
            scope_rows.append(dict(slot=slot,kind=kind,program=item['program'],view_A=item['view_A'],view_B=item['view_B'],n_patients=len(x),status=status,nuisance_rank=m['nuisance_rank'] if m else np.nan,residual_df=m['residual_df'] if m else np.nan,condition=m['condition'] if m else np.nan))
    pd.DataFrame(scope_rows).to_csv(out/'private_all12_input_design_scope.tsv',sep='\t',index=False)
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(out.iterdir()) if p.is_file()}
    t=dict(status='COMPLETE_ALL12_PAIRED_INPUTS_FULL_DESIGNS_AND_COMPOSITION_SENSITIVITY_NO_EFFECTS',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),files=files)
    (out/'private_terminal.json').write_text(json.dumps(t,indent=2)+'\n');print(json.dumps(dict(status=t['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
