"""Independent bounded KL-deviance check of all fixed v1 projection fixtures."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import argparse,json,hashlib
from datetime import datetime,timezone
import numpy as np
from scipy.optimize import least_squares

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run',required=True);a=ap.parse_args();r=Path(a.run);out=r/'public/projection_deviance_audit_v2.json';assert not out.exists();rows=[]
    for p in sorted((r/'private').glob('activities_*.npz')):
        d=np.load(p);W=d['W'];V=d['V'];H=d['H'];scale=W.sum(axis=0);A=W/scale
        for j in range(min(3,V.shape[1])):
            v=V[:,j];target=H[:,j]*scale
            def residual(h):
                z=np.maximum(np.einsum('ij,j->i',A,h),1e-15);dev=z.copy();nz=v>0;t=z[nz]/v[nz]-1
                dev[nz]=v[nz]*(t-np.log1p(t));return np.sign(z-v)*np.sqrt(2*np.maximum(dev,0))
            def objective(h):
                z=np.maximum(np.sum(A*h[None,:],axis=1),1e-15);nz=v>0
                return float(np.sum(z-v)+np.sum(v[nz]*(np.log(v[nz])-np.log(z[nz]))))
            best=None
            for fac in [.7,1.3]:
                initial=np.full(A.shape[1],max(v.sum()/A.shape[1],1e-5)*fac)
                fit=least_squares(residual,initial,bounds=(np.full(A.shape[1],1e-12),np.full(A.shape[1],np.inf)),method='trf',ftol=1e-12,xtol=1e-12,gtol=1e-12,max_nfev=10000)
                value=objective(fit.x)
                if best is None or value<best[0]:best=(value,fit)
            prod=objective(target);gap=prod-best[0];grad=np.array([np.sum(A[:,k]*(1-v/np.maximum(A@target,1e-15))) for k in range(A.shape[1])]);kkt=float(np.max(abs(np.where(target>1e-7,grad,np.minimum(grad,0)))))
            rows.append(dict(fixture_file_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),fixture_index=j,production_objective=prod,audit_objective=best[0],objective_gap=float(gap),production_KKT=kkt,audit_success=bool(best[1].success),pass_objective=bool(abs(gap)<1e-5*(1+abs(best[0])) and kkt<1e-3)))
    result=dict(schema='KL_deviance_independent_projection_fixture_audit_v2',observed_UTC=datetime.now(timezone.utc).isoformat(),files=len(list((r/'private').glob('activities_*.npz'))),fixtures=len(rows),passed=sum(x['pass_objective'] for x in rows),status='PASS_ALL_SELECTED_PROJECTION_FIXTURES' if all(x['pass_objective'] for x in rows) else 'FAIL_PRESERVE_UNMATCHED_FIXTURES',rows=rows,scientific_inputs_or_activities_changed=False,independent_biological_validation=False)
    out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='rows'}),flush=True)

if __name__=='__main__':main()
