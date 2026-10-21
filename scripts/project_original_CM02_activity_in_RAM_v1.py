"""Read-only original full-basis projection, no learned/individual export.

Reuses the exact original R parity contract; no new NMF, scale fit, outcome
test or restricted seven-member refit is performed. The returned patient
mapping belongs only to the caller's RAM.
"""
import ast
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
import pandas as pd
from scipy.optimize import nnls,lsq_linear

SCOPE='T_whole_TME_primary_tumor'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def project(root,atlas,tab,prior_source,prior_source_sha):
    assert sha(prior_source)==prior_source_sha
    statements=ast.parse(Path(prior_source).read_text()).body
    candidates=[s.value.value for s in statements if isinstance(s,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='R_CHECK' for t in s.targets)]
    assert len(candidates)==1 and isinstance(candidates[0],str)
    old=root/'HNSCC_new_CM_frozen_model_export_20261016_v2'/SCOPE;scope=atlas/'CM'/SCOPE
    files={name:old/name for name in ['A_private.tsv','H_private.tsv','F_private.tsv','norm_private.tsv','scaling_private.tsv']}
    files.update(frequency=scope/'inputs/frequency_private.tsv',metadata=scope/'inputs/sample_metadata_private.tsv',
        original_model=scope/'discovery_v3/NMF_discovery.rds',
        prior_baseline=root/'HNSCC_new_CM_baseline_preparation_20261016_v1/private_patient_baseline_component_scores.tsv')
    previous=json.loads((old.parent/'private_terminal.json').read_text());accepted=next(r for r in previous['scopes'] if r['scope']==SCOPE)
    assert sha(files['A_private.tsv'])==accepted['effective_basis_sha256']=='8fa2e933d2fd109da7c0d2e4a1a88bd1c31edaf3c93fb2ab09ed012792295da0'
    assert sha(files['norm_private.tsv'])==accepted['normalization_sha256']=='916efbace913881d2d3f9b683e0dda17356baaabf2c24edbb8a56f62aef740af'
    launch=json.loads((atlas/'manifests/CM_T_whole_TME_primary_tumor_recovery_v2_started_private.json').read_text());command=launch['command']
    assert command[-1]=='discovery_v3' and Path(command[2]).resolve()==scope.resolve()
    assert sha(command[1])=='4cf9af42279924f071ec1cffdbf49f1467cd4693aa85732bad6de9ad08485feb'
    helper=Path(command[3])/'R/CoVarNet.R'
    checked=subprocess.run([command[0],'--vanilla','-e',candidates[0],str(scope),str(old),str(helper)],capture_output=True,text=True,timeout=120)
    assert checked.returncode==0,'original_model_R_parity_failed'
    contract=json.loads(checked.stdout)
    A,H,F,norm,scale=[pd.read_csv(files[name],sep='\t',index_col=0) for name in ['A_private.tsv','H_private.tsv','F_private.tsv','norm_private.tsv','scaling_private.tsv']]
    fq=pd.read_csv(files['frequency'],sep='\t',index_col=0)
    meta=pd.read_csv(files['metadata'],sep='\t',dtype=str,keep_default_na=False).set_index('sampleID')
    assert all(A.index.equals(x.index) for x in [fq,norm,scale])
    assert list(A.columns)==list(H.index) and list(norm.columns)==list(fq.columns[1:]) and set(norm.columns)==set(meta.index)
    np.testing.assert_allclose((fq.iloc[:,1:].to_numpy()-scale.minimum.to_numpy()[:,None])/scale['range'].to_numpy()[:,None],norm,atol=1e-12,rtol=1e-12)
    a=A.to_numpy();assert np.linalg.matrix_rank(a)==a.shape[1] and np.linalg.cond(a)<=1e8
    for j in range(H.shape[1]):np.testing.assert_allclose(nnls(a,F.iloc[:,j].to_numpy())[0],H.iloc[:,j],atol=1e-6,rtol=1e-7)
    source=meta.loc[meta.dataset.eq('GSE286827')].copy();source['library']=source.index.str.removeprefix('GSE286827::')
    assert set(source.library)==set(tab.library_id)
    registry=tab.set_index('library_id');source['patient']=source.library.map(registry.patient_registry)
    assert source.library.map(registry.timepoint_registry).to_numpy().tolist()==source.timepoint.tolist()
    assert source.groupby('independent_unit').patient.nunique().max()==1 and source.groupby('patient').independent_unit.nunique().max()==1
    legacy=pd.read_csv(files['prior_baseline'],sep='\t');legacy=legacy.loc[legacy.scope.eq(SCOPE)].set_index(['patient_id','component'])
    activity={};residual={};j=list(A.columns).index('CM02');precount=0
    for (case,phase),g in source.groupby(['patient','timepoint'],sort=True):
        y=norm.loc[:,g.index].mean(axis=1).to_numpy();h,error=nnls(a,y,maxiter=10000)
        other=lsq_linear(a,y,bounds=(0,np.inf),method='bvls',tol=1e-12,max_iter=10000);assert other.success
        np.testing.assert_allclose(h,other.x,atol=1e-7,rtol=1e-7)
        gradient=a.T@(a@h-y);tol=1e-7*(1+np.linalg.norm(a,2)*np.linalg.norm(y))
        assert gradient.min()>=-tol and np.max(np.abs(h*gradient))<=tol*(1+np.linalg.norm(h))
        if phase=='pre':
            np.testing.assert_allclose(h,[legacy.loc[(case,c),'activity'] for c in A.columns],atol=1e-7,rtol=1e-7);precount+=1
        activity[case,phase]=float(h[j]);residual.setdefault(phase,[]).append(float(error/max(np.linalg.norm(y),1e-12)))
    receipt=dict(original_R_parity=contract,input_SHA256={k:sha(p) for k,p in files.items()},
        original_helper_SHA256=sha(helper),all_original_H_NNLS_replayed=True,all_phase_NNLS_BVLS_KKT_passed=True,
        prior_baseline_exact_replay_cases=precount,projection_relative_residual_quantiles={ph:np.quantile(v,[0,.25,.5,.75,1]).tolist() for ph,v in residual.items()},
        original_full_basis_minmax_used=True,no_new_NMF_or_individual_or_learned_arrays_written=True)
    return activity,receipt
