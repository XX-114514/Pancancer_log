"""Independent Python TMM and augmented-QR replay of every native universe."""
import json
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular
from scipy.stats import spearmanr

from NSCLC_native_MFA_interpretation_inputs_v1 import load,sha,VIEWS
from NSCLC_native_MFA_clinical_inputs_v1 import matched_quota
from audit_NSCLC_native_MFA_inputs_v1 import tmm
from external_functional_coverage_inputs_v1 import read_inputs


def main():
    data=load();root=data['root'];matched_quota(root)
    src=root/'NSCLC_native_MFA_external_RNA_universe_calibration_20261016_v1'
    terminal=src/'private_terminal.json';t=json.loads(terminal.read_text())
    assert t['status'].startswith('COMPLETE_ALL61_SOURCE_NATIVE_RNA_UNIVERSES_')
    for name,item in t['files'].items():assert sha(src/name)==item['sha256']
    with np.load(src/'private_all_universe_native_scores.npz',allow_pickle=False) as z:
        masks=z['masks'];saved=z['scores']
    assert masks.shape==(48,31831) and saved.shape==(4,48,208,10) and masks[0].all()
    mask_table=pd.read_csv(src/'private_native_universe_masks.tsv.gz',sep='\t',index_col=0)
    assert mask_table.index.astype(str).tolist()==data['genes'].tolist() and np.array_equal(mask_table.to_numpy().T.astype(bool),masks)
    source_map=pd.read_csv(src/'private_all61_source_native_universe_map.tsv',sep='\t').set_index(['arm','dataset'])
    axes=read_inputs(root)['source_axes'];assert len(source_map)==len(axes)==61
    for key,axis in axes.items():
        cols=[next(iter(axis['lookup'][g])) if g in axis['lookup'] and len(axis['lookup'][g])==1 else None for g in data['genes']]
        count=Counter(c for c in cols if c is not None)
        mask=np.array([c is not None and count[c]==1 for c in cols])
        row=source_map.loc[key];assert np.array_equal(mask,masks[int(row.mask_index)]) and row.measured_native_universe_genes==int(mask.sum())
    table=pd.read_csv(src/'private_all2440_native_RNA_universe_calibrations.tsv',sep='\t').set_index(['arm','dataset','view','factor'])
    prior=pd.read_csv(root/'NSCLC_native_MFA_external_gene_mask_calibration_20261016_v1/private_all2440_source_mask_calibration.tsv',sep='\t').set_index(['arm','dataset','view','factor'])
    assert len(table)==2440 and table.index.is_unique
    cases=data['cases'].sampleID.astype(str).tolist();train=data['train'];inp=root/'NSCLC_native_MFA_inputs_20261016_v1'
    maximum_tmm=0.;maximum_score=0.;checked=0;gates_checked=0
    checkpoint=root/'HNSCC_receiver_followup_job_control_20261016_v2'/Path(__file__).stem/'private_audit_checkpoint.json'
    for vi,view in enumerate(VIEWS):
        counts=pd.read_csv(inp/('private_'+view+'_pseudobulk.tsv.gz'),sep='\t',index_col=0)
        assert counts.columns.tolist()==cases and counts.index.astype(str).tolist()==data['genes'].tolist()
        raw=counts.to_numpy().T;assert np.array_equal(raw,raw.astype(np.int64))
        norm=pd.read_csv(src/('private_'+view+'_all_universe_TMM.tsv.gz'),sep='\t',dtype={'sampleID':str})
        assert len(norm)==48*208 and not norm.duplicated(['mask_id','sampleID']).any()
        p=data['prepared'][vi];selected=p['selected_native_columns'];w=data['model']['weights_'+str(vi)]
        sw=np.sqrt(data['model']['precision_'+str(vi)])
        q,r=np.linalg.qr(np.vstack([sw[:,None]*w,np.eye(10)]),mode='reduced')
        a=solve_triangular(r,q[:len(w)].T*sw)
        correction=p['design']@p['regression_coefficients'];full=data['model']['own_view_scores_'+str(vi)]
        sd=full[train].std(axis=0,ddof=1);rms=float(p['training_view_RMS']);all_scores=[]
        for mi,mask in enumerate(masks):
            rr=raw[:,mask];reference=rr[train].sum(axis=0)
            factors=np.array([tmm(row,reference) for row in rr]);centered=factors/np.exp(np.log(factors[train]).mean())
            libraries=rr.sum(axis=1);effective=libraries*centered
            observed=norm[norm.mask_id.eq('mask_%03d'%mi)].set_index('sampleID').loc[cases]
            delta=float(np.max(abs(factors-observed.raw_factor.to_numpy())));maximum_tmm=max(maximum_tmm,delta)
            assert np.allclose(factors,observed.raw_factor,rtol=1e-10,atol=1e-12),delta
            assert np.array_equal(libraries,observed.library_size)
            assert np.allclose(centered,observed.factor,rtol=1e-10,atol=1e-12)
            assert np.allclose(effective,observed.effective_library_size,rtol=1e-10,atol=1e-6)
            take=mask[selected]
            value=np.log2(1+raw[:,selected[take]].astype(float)*1e6/effective[:,None])
            score=((value-correction[:,take])/rms)@a[:,take].T
            error=float(np.max(abs(score-saved[vi,mi])));maximum_score=max(maximum_score,error)
            assert np.allclose(score,saved[vi,mi],atol=1e-8,rtol=1e-9),error
            if mi==0:assert np.allclose(score,full,atol=1e-8,rtol=1e-9)
            all_scores.append(score);checked+=score.size
            checkpoint.write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),view=view,completed_masks=mi+1,total_masks=len(masks)))+'\n')
        for key,row in source_map.iterrows():
            score=all_scores[int(row.mask_index)]
            for k in range(10):
                factor='NSCLC_MFA_20261016_F%02d'%(k+1);index=(*key,view,factor);observed=table.loc[index];old=prior.loc[index]
                assert bool(observed.original_physical_gate)==bool(old.physical_gene_gate)
                assert bool(observed.original_prepared_X_mask_gate)==bool(old.native_mask_calibration_pass)
                passed=bool(old.native_mask_calibration_pass)
                for name,take in [('discovery',train),('internal_holdout',~train)]:
                    x=full[take,k];y=score[take,k];dx=x-x.mean();dy=y-y.mean();delta=y-x
                    pearson=float(dx@dy/np.sqrt((dx@dx)*(dy@dy))) if y.std()>1e-12 else np.nan
                    rank=float(spearmanr(x,y).statistic) if y.std()>1e-12 else np.nan
                    rmse=float(np.sqrt(delta@delta/len(delta))/sd[k]);shift=float(delta.mean()/sd[k])
                    assert np.allclose([pearson,rank,rmse,shift],[observed[name+'_'+v] for v in ['Pearson','Spearman','RMSE_SD','mean_shift_SD']],atol=1e-8,rtol=1e-8,equal_nan=True)
                    passed=passed and pearson>=.95 and rmse<=.25
                assert bool(passed)==bool(observed.native_RNA_universe_gate);gates_checked+=1
    assert checked==399360 and gates_checked==2440
    proof=src/'private_independent_native_RNA_universe_audit.json';assert not proof.exists()
    report=dict(status='PASS_ALL48_NATIVE_COUNT_UNIVERSES39936_PYTHON_TMM399360_SCORE_VALUES_AND2440_GATES',
        finished_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),producer_terminal_sha256=sha(terminal),
        python_TMM_formula_source_sha256=sha(Path(__file__).with_name('audit_NSCLC_native_MFA_inputs_v1.py')),
        maximum_raw_TMM_factor_error=maximum_tmm,maximum_score_error=maximum_score,
        external_expression_or_clinical_validation=False,original_quality_held_fixed_only_gene_universe_perturbed=True)
    proof.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],finished_UTC=report['finished_UTC'],audit_sha256=sha(proof))))


if __name__=='__main__':main()
