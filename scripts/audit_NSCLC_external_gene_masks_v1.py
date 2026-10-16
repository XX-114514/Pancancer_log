"""Augmented QR/explicit column omission replay of all native mask calibrations."""
import json
from collections import Counter
from datetime import datetime,timezone

import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular

from NSCLC_native_MFA_interpretation_inputs_v1 import load,sha,VIEWS
from NSCLC_native_MFA_clinical_inputs_v1 import matched_quota
from external_functional_coverage_inputs_v1 import read_inputs


def main():
    data=load();root=data['root'];src=root/'NSCLC_native_MFA_external_gene_mask_calibration_20261016_v1'
    terminal=json.loads((src/'private_terminal.json').read_text())
    for name,item in terminal['files'].items():assert sha(src/name)==item['sha256']
    proof=src/'private_independent_native_mask_audit.json';assert not proof.exists();matched_quota(root)
    table=pd.read_csv(src/'private_all2440_source_mask_calibration.tsv',sep='\t').set_index(['arm','dataset','view','factor'])
    with np.load(src/'private_all244_masks_and507520_native_partial_scores.npz',allow_pickle=False) as z:scores=z['scores'];masks=z['masks'];keys=z['keys'].astype(str)
    axes=read_inputs(root)['source_axes'];physical=pd.read_csv(root/'NSCLC_native_MFA_external_gene_qualification_20261016_v1/private_all2440_native_factor_source_gene_coverage.tsv',sep='\t').set_index(['arm','dataset','view','factor'])
    assert len(table)==2440 and len(keys)==244 and table.index.is_unique
    checked=0;maximum_score_error=0.;maximum_metric_error=0.
    for vi,view in enumerate(VIEWS):
        p=data['prepared'][vi];w=data['model']['weights_'+str(vi)];sw=np.sqrt(data['model']['precision_'+str(vi)])
        q,r=np.linalg.qr(np.vstack([sw[:,None]*w,np.eye(10)]),mode='reduced');coef=solve_triangular(r,q[:len(w)].T*sw)
        full=data['model']['own_view_scores_'+str(vi)];training=data['train'];sd=full[training].std(axis=0,ddof=1)
        genes=data['genes'][p['selected_native_columns']]
        for (arm,dataset),axis in sorted(axes.items()):
            target=[next(iter(axis['lookup'][g])) if g in axis['lookup'] and len(axis['lookup'][g])==1 else None for g in genes]
            count=Counter(v for v in target if v is not None);valid=np.array([v is not None and count[v]==1 for v in target])
            idx=np.where(keys==arm+'|'+dataset+'|'+view)[0];assert len(idx)==1;idx=int(idx[0]);assert np.array_equal(valid,masks[idx])
            partial=p['X'][:,valid]@coef[:,valid].T
            error=float(np.max(abs(partial-scores[idx])));assert error<1e-9;maximum_score_error=max(maximum_score_error,error);checked+=partial.size
            for k in range(10):
                key=(arm,dataset,view,'NSCLC_MFA_20261016_F%02d'%(k+1));row=table.loc[key];passed=bool(physical.loc[key,'physical_coverage_pass'])
                for name,use in [('discovery',training),('internal_holdout',~training)]:
                    x=full[use,k];y=partial[use,k];dx=x-x.mean();dy=y-y.mean();vary=x.std()>1e-12 and y.std()>1e-12 and sd[k]>1e-12
                    pearson=float(dx@dy/np.sqrt((dx@dx)*(dy@dy))) if vary else np.nan
                    rx=pd.Series(x).rank().to_numpy();ry=pd.Series(y).rank().to_numpy();rx-=rx.mean();ry-=ry.mean()
                    spearman=float(rx@ry/np.sqrt((rx@rx)*(ry@ry))) if vary else np.nan
                    delta=y-x;rmse=float(np.sqrt(delta@delta/len(delta))/sd[k]) if sd[k]>1e-12 else np.nan
                    values=[pearson,spearman,rmse,float(delta.mean()/sd[k]),float(y.std()/x.std())]
                    observed=[row[name+'_'+field] for field in ['Pearson','Spearman','RMSE_SD','mean_shift_SD','SD_ratio']]
                    assert np.allclose(values,observed,atol=1e-8,rtol=1e-8,equal_nan=True)
                    good=np.isfinite(values)&np.isfinite(observed)
                    if good.any():maximum_metric_error=max(maximum_metric_error,float(np.max(abs(np.array(values)[good]-np.array(observed)[good]))))
                    passed=passed and pearson>=.95 and rmse<=.25
                assert bool(row.native_mask_calibration_pass)==bool(passed)
    assert checked==507520
    report=dict(status='PASS_ALL244_NATIVE_SOURCE_MASKS_507520_PARTIAL_SCORES_AND2440_CALIBRATION_GATES_INDEPENDENT_QR',finished_UTC=datetime.now(timezone.utc).isoformat(),
        source_sha256=sha(__file__),producer_terminal_sha256=sha(src/'private_terminal.json'),partial_score_values_checked=checked,
        maximum_score_error=maximum_score_error,maximum_metric_error=maximum_metric_error,external_expression_or_clinical_validation=False)
    proof.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],finished_UTC=report['finished_UTC'],audit_sha256=sha(proof))))


if __name__=='__main__':main()
