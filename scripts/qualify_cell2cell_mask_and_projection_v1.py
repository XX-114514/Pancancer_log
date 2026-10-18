"""Synthetic-only official CP, missing-mask and fixed-factor projection gate."""
import json
import os
import sys
from datetime import datetime,timezone
from pathlib import Path


def worker():
    import importlib.metadata as md
    import inspect
    import numpy as np
    import torch
    from scipy.optimize import nnls
    from packaging.requirements import Requirement
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=Path(os.environ['PANCANCER_SYNTHETIC_OUT'])
    sites=[root/('Tensor_cell2cell_GPU_runtime_20261019_'+v)/'site_packages' for v in ['v3','v2']];sys.path[:0]=list(map(str,sites))
    import cell2cell as c2c
    import tensorly as tl
    from tensorly.decomposition import non_negative_parafac
    expected={'cell2cell':'0.8.0','tensorly':'0.9.0','gseapy':'1.0.3','requests':'2.32.5','numpy':'1.26.4'}
    for name,version in expected.items():assert md.version(name)==version
    seen=set();todo=['cell2cell','tensorly','gseapy','requests']
    while todo:
        name=todo.pop()
        if name.lower() in seen:continue
        seen.add(name.lower())
        for text in md.distribution(name).requires or []:
            r=Requirement(text)
            if r.marker and not r.marker.evaluate({'extra':''}):continue
            v=md.version(r.name);assert not r.specifier or r.specifier.contains(v,prereleases=True);todo.append(r.name)
    assert torch.cuda.is_available();tl.set_backend('pytorch');torch.backends.cuda.matmul.allow_tf32=False;torch.backends.cudnn.allow_tf32=False
    rng=np.random.default_rng(2026101911);shape=(12,9,4,3);factors=[rng.uniform(.4,1.4,(n,2)) for n in shape]
    factors[1][:4,1]=0;factors[1][4:,0]=0;factors[1][-1,:]=0
    x=np.einsum('ak,bk,ck,dk->abcd',*factors);mask=(rng.random(shape)>.18).astype(float);mask[:,8,:,:]=1
    assert np.count_nonzero((x==0)&(mask==1))>0
    def fit(device,fill,iterations):
        t=tl.tensor(np.where(mask>0,x,fill),dtype=torch.float64,device=device);m=tl.tensor(mask,dtype=torch.float64,device=device)
        cp=non_negative_parafac(t,rank=2,init='random',random_state=2026101912,n_iter_max=iterations,tol=0,mask=m,normalize_factors=False)
        return cp,tl.to_numpy(tl.cp_to_tensor(cp))
    cp,recon=fit('cuda:0',0,1200);other,otherrecon=fit('cuda:0',37,1200)
    np.testing.assert_allclose(recon,otherrecon,atol=1e-9,rtol=1e-9)
    cpu,cpurecon=fit('cpu',0,1200);np.testing.assert_allclose(recon,cpurecon,atol=1e-7,rtol=1e-7)
    rel=float(np.linalg.norm((x-recon)*mask)/np.linalg.norm(x*mask));assert rel<.005
    order_names=[[str(i) for i in range(n)] for n in shape]
    tensor=c2c.tensor.PreBuiltTensor(tl.tensor(np.where(mask>0,x,0),dtype=torch.float64,device='cuda:0'),order_names=order_names,mask=tl.tensor(mask,dtype=torch.float64,device='cuda:0'))
    tensor.compute_tensor_factorization(rank=2,init='random',random_state=2026101912,runs=1,normalize_loadings=False,var_ordered_factors=False,n_iter_max=1200,tol=0,verbose=False)
    np.testing.assert_allclose(recon,tl.to_numpy(tl.cp_to_tensor(tensor.tl_object)),atol=1e-9,rtol=1e-9)
    np.testing.assert_array_equal(tl.to_numpy(tensor.mask),mask)
    assert int(((tl.to_numpy(tensor.tensor)==0)&(tl.to_numpy(tensor.mask)==1)).sum())==int(((x==0)&(mask==1)).sum())
    known=rng.uniform(.3,1.3,(5,2));design=np.einsum('lk,sk,rk->lsrk',*factors[1:]).reshape(-1,2)
    heldout=np.einsum('ak,bk,ck,dk->abcd',known,*factors[1:]);projection=[];normal_equation_error=0.
    for i,y in enumerate(heldout):
        observed=rng.random(y.shape)>.23;observed[-1,:,:]=True;D=design[observed.ravel()];v=y[observed];assert np.linalg.matrix_rank(D)==2
        solution,_=nnls(D,v);normal=np.linalg.solve(D.T@D,D.T@v);np.testing.assert_allclose(solution,known[i],atol=1e-10,rtol=0);np.testing.assert_allclose(solution,normal,atol=1e-10,rtol=0)
        # Missing placeholders cannot enter the observed projection.
        bad=np.where(observed,y,999.);again,_=nnls(D,bad[observed]);np.testing.assert_array_equal(solution,again)
        normal_equation_error=max(normal_equation_error,float(np.max(np.abs(solution-normal))));projection.append(solution.tolist())
    receipt={'status':'ACCEPTED_OFFICIAL_CELL2CELL_TENSORLY_CPU_GPU_MASK_TRUE_ZERO_AND_HELDOUT_PROJECTION_SYNTHETIC_GATE','completed_UTC':datetime.now(timezone.utc).isoformat(),'versions':{n:md.version(n) for n in ['cell2cell','tensorly','gseapy','requests','numpy','scipy','torch']},'dependency_distributions_checked':len(seen),'synthetic_shape':shape,'known_rank':2,'observed_true_zero_count':int(((x==0)&(mask==1)).sum()),'observed_relative_reconstruction_error':rel,'placeholder_invariance_max_error':float(np.max(np.abs(recon-otherrecon))),'CPU_GPU_reconstruction_max_error':float(np.max(np.abs(recon-cpurecon))),'fixed_noncontext_projection_normal_equation_max_error':normal_equation_error,'APIs':{'PreBuiltTensor':str(inspect.signature(c2c.tensor.PreBuiltTensor)),'compute_tensor_factorization':str(inspect.signature(c2c.tensor.PreBuiltTensor.compute_tensor_factorization)),'non_negative_parafac':str(inspect.signature(non_negative_parafac))},'tensorly_core_identical_to_cell2cell_wrapper_fixture':True,'runtime_layers_sha256':{v:__import__('hashlib').sha256((root/('Tensor_cell2cell_GPU_runtime_20261019_'+v)/'private_official_binary_dependency_manifest.json').read_bytes()).hexdigest() for v in ['v3','v2']},'new_real_expression_or_clinical_effects':False}
    (out/'private_worker_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


def main():
    import subprocess
    from NSCLC_native_MFA_interpretation_inputs_v1 import sha
    from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,2*1024**2);out=root/'Tensor_cell2cell_synthetic_gate_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    plan={'created_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'resource_guard':guard,'gate':'No additional software installation. Existing official V2+V3 isolated layers. V3 imports succeeded; its receipt failed only at wrong PrebuiltTensor capitalization; actual official API is PreBuiltTensor. Retain all failed records. Known rank2/true zero/missing placeholder invariance/CPU-GPU and official wrapper-engine parity; independent known heldout NNLS/normal-equation projection. No real input or clinical outcome read.','solver_contract':'Official Tensorly multiplicative-update nonnegative CP with random init, fixed iteration count/tol0 avoids missing-placeholder dependence through early-stopping error normalization; observed reconstruction selects starts. The Cell2cell wrapper calls this exact engine, demonstrated by fixture.'}
    (out/'plan_before_synthetic_gate.json').write_text(json.dumps(plan,indent=2)+'\n')
    python=Path(sys.executable).resolve().parents[2]/'scvi/bin/python';env=os.environ.copy();env.update(PANCANCER_SYNTHETIC_OUT=str(out),PYTHONDONTWRITEBYTECODE='1',NUMBA_CACHE_DIR=str(out/'numba_cache'),MPLCONFIGDIR=str(out/'matplotlib_cache'),MPLBACKEND='Agg')
    with (out/'private_worker_stdout.log').open('xb') as so,(out/'private_worker_stderr.log').open('xb') as se:r=subprocess.run([str(python),'-u',str(Path(__file__).resolve()),'--worker'],env=env,stdout=so,stderr=se,timeout=300)
    assert r.returncode==0,'SYNTHETIC_WORKER_FAILED';receipt=json.loads((out/'private_worker_receipt.json').read_text());receipt.update(source_sha256=sha(__file__),resource_guard=guard,worker_receipt_sha256=sha(out/'private_worker_receipt.json'))
    assert sum(p.stat().st_size for p in out.rglob('*') if p.is_file())<2*1024**2
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'status':receipt['status'],'terminal_sha256':sha(out/'private_terminal.json')}))


if __name__=='__main__':
    if sys.argv[1:]==['--worker']:worker()
    else:main()
