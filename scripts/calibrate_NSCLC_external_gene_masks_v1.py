"""All source masks on original native cases, no external expression or effects."""
import json
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from NSCLC_native_MFA_interpretation_inputs_v1 import load,sha,VIEWS
from NSCLC_native_MFA_clinical_inputs_v1 import matched_quota
from external_functional_coverage_inputs_v1 import read_inputs


def main():
    data=load();root=data['root'];genesrc=root/'NSCLC_native_MFA_external_gene_qualification_20261016_v1'
    proof=genesrc/'private_independent_gene_qualification_audit.json';assert sha(proof)=='096d5166f1016fac39ca4feea462a3399eb7b120c3adfa1e8166d20f74c09a22'
    out=root/'NSCLC_native_MFA_external_gene_mask_calibration_20261016_v1';assert not out.exists();matched_quota(root);out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),primary_model_sha256=sha(data['model_path']),physical_gene_audit_sha256=sha(proof),
        question='For all61 source-arm x4 views x10 factors, does the exact available-gene partial frozen projection preserve the original native score?',
        masks='Prior unique physical mapping, same symbol/alias/collision rules. Missing gene contributions are omitted from a partial signature; not expression-zero imputation. All failed physical masks also retained.',
        measurement='Original prepared native X, original fixed own-view projection A. Full scores must exactly replay saved own-view scores. Partial score=X @ (A * available_mask). No refit, offset correction, weight rescaling or new gene selection.',
        partitions='Original discovery145 and internal_holdout63 separately; no clinical outcome used.',
        gates='Physical gene gate unchanged; in BOTH native partitions require Pearson>=.95 and RMSE<=.25 original DISCOVERY score SD(ddof1), fixed before calibration. Spearman, signed mean shift and variance ratio also descriptive. Every2440 slot retained.',
        limits='Engineering measurement gate on original native cases, not external expression comparability, lineage/patient validation or biological replication. RNA count/TMM/QC definition remains an external-stage requirement.',
        external_expression_or_P=False,new_data_download=False,old_MFA_or_CM_changed=False)
    (out/'plan_before_any_native_gene_mask_calibration.json').write_text(json.dumps(plan,indent=2)+'\n')
    axes=read_inputs(root)['source_axes'];physical=pd.read_csv(genesrc/'private_all2440_native_factor_source_gene_coverage.tsv',sep='\t').set_index(['arm','dataset','view','factor'])
    records=[];masks=[];scores=[];keys=[]
    for vi,view in enumerate(VIEWS):
        p=data['prepared'][vi];w=data['model']['weights_'+str(vi)];tau=data['model']['precision_'+str(vi)]
        a=np.linalg.solve(np.eye(10)+w.T@(tau[:,None]*w),w.T*tau)
        full=p['X']@a.T;saved=data['model']['own_view_scores_'+str(vi)];assert np.allclose(full,saved,atol=1e-10,rtol=1e-10)
        sd=full[data['train']].std(axis=0,ddof=1);gene=data['genes'][p['selected_native_columns']]
        for (arm,dataset),axis in sorted(axes.items()):
            columns=[next(iter(axis['lookup'][g])) if g in axis['lookup'] and len(axis['lookup'][g])==1 else -1 for g in gene]
            count=pd.Series([c for c in columns if c>=0]).value_counts().to_dict()
            mask=np.array([c>=0 and count[c]==1 for c in columns]);partial=p['X']@(a*mask).T
            mask_index=len(masks);masks.append(mask);scores.append(partial);keys.append(arm+'|'+dataset+'|'+view)
            for k in range(10):
                factor='NSCLC_MFA_20261016_F%02d'%(k+1);row=dict(arm=arm,dataset=dataset,view=view,factor=factor,mask_index=mask_index,measured_features=int(mask.sum()),physical_gene_gate=bool(physical.loc[(arm,dataset,view,factor),'physical_coverage_pass']))
                for name,take in [('discovery',data['train']),('internal_holdout',~data['train'])]:
                    x=full[take,k];y=partial[take,k];delta=y-x
                    varying=x.std()>1e-12 and y.std()>1e-12 and sd[k]>1e-12
                    row[name+'_Pearson']=float(np.corrcoef(x,y)[0,1]) if varying else np.nan
                    row[name+'_Spearman']=float(spearmanr(x,y).statistic) if varying else np.nan
                    row[name+'_RMSE_SD']=float(np.sqrt(np.mean(delta**2))/sd[k]) if sd[k]>1e-12 else np.nan
                    row[name+'_mean_shift_SD']=float(delta.mean()/sd[k]) if sd[k]>1e-12 else np.nan
                    row[name+'_SD_ratio']=float(y.std()/x.std()) if x.std()>1e-12 else np.nan
                row['native_mask_calibration_pass']=bool(row['physical_gene_gate'] and all(row[name+'_Pearson']>=.95 and row[name+'_RMSE_SD']<=.25 for name in ['discovery','internal_holdout']))
                records.append(row)
    assert len(records)==2440 and len(masks)==244
    pd.DataFrame(records).to_csv(out/'private_all2440_source_mask_calibration.tsv',sep='\t',index=False)
    np.savez_compressed(out/'private_all244_masks_and507520_native_partial_scores.npz',masks=np.array(masks),scores=np.array(scores),keys=np.array(keys))
    terminal=dict(status='COMPLETE_ALL244_SOURCE_VIEW_MASKS_AND2440_FACTOR_NATIVE_CALIBRATIONS_AUDIT_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),
        files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()},external_expression_or_clinical_effects=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'],finished_UTC=terminal['finished_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
