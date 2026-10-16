"""All external TMM/QR gene indices, weighted least-squares QC and statuses."""
import hashlib
import json
import os
import shutil
import subprocess
import time
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular

from NSCLC_native_MFA_interpretation_inputs_v1 import load,sha,VIEWS
from diagnose_NSCLC_native_RNA_TMM_rounding_v1 import calculate

CONTEXT=['arm','dataset','cancer_label','sample_class','material_class','capture_scope','ex_vivo']


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    src=root/'NSCLC_external_fixed_gene_programme_assay_20261016_v1'
    receipt=json.loads((root/'HNSCC_receiver_followup_job_control_20261016_v2/prepare_NSCLC_external_fixed_gene_programme_assay_v1/launch_receipt.json').read_text())
    deadline=time.monotonic()+7200
    while not (src/'private_terminal.json').exists():
        proc=Path('/proc')/str(receipt['pid'])/'stat';assert proc.exists(),'Owned assay producer missing without terminal'
        f=proc.read_text().rsplit(')',1)[1].split();assert f[0]!='Z' and f[19]==str(receipt['proc_start_ticks'])
        assert time.monotonic()<deadline,'Audit observation expired; do not infer producer stopped'
        time.sleep(3)
    fs=subprocess.check_output(['df','-P',str(root.resolve())],text=True).splitlines()[-1].split()[0]
    q=subprocess.check_output(['quota','-g','-v','-w','USER002'],text=True)
    rows=[s.split() for s in q.splitlines() if s.split() and s.split()[0]==fs];assert len(rows)==1
    free=(int(rows[0][3])-int(rows[0][1].rstrip('*')))*1024
    assert min(free,shutil.disk_usage(root).free)>500_000_000_000+1073741824+1048576
    tpath=src/'private_terminal.json';terminal=json.loads(tpath.read_text())
    assert terminal['source_sha256']==receipt['source_sha256'] and terminal['status'].startswith('COMPLETE_ALL61_SOURCE_FIXED_GENE_PROGRAMME_')
    for name,item in terminal['files'].items():assert sha(src/name)==item['sha256']
    data=load();inp=root/'NSCLC_native_MFA_inputs_20261016_v1';cal=root/'NSCLC_native_MFA_external_RNA_universe_calibration_20261016_v1'
    count=root/'NSCLC_native_MFA_external_frozen_role_count_inputs_20261016_v2'
    assert sha(count/'private_independent_native_count_input_audit.json')=='7bd028cc1c89ac7324ac25c65042b5b184cb377162a38f9544caf7bb24048c04'
    units=pd.read_csv(src/'private_all6760_programme_units_metadata_quality_TMM.tsv.gz',sep='\t',keep_default_na=False)
    original_units=pd.read_csv(count/'private_all6760_native_count_units_and_quality.tsv.gz',sep='\t',keep_default_na=False)
    assert len(units)==6760
    for col in original_units.columns:assert np.array_equal(units[col],original_units[col]),col
    with np.load(src/'private_all6760_fixed_programme_indices.npz',allow_pickle=False) as z:saved={k:z[k] for k in z.files}
    quality=np.column_stack([np.log1p(units.n_cells),units[['mean_log_full_source_RNA','mean_log_full_source_detected_genes','symbol_pct_mito_RNA','symbol_pct_ribo_RNA']].apply(pd.to_numeric,errors='coerce').to_numpy()])
    assert np.allclose(quality,saved['quality'],equal_nan=True)
    measurable=units.n_cells.ge(20).to_numpy() & units.full_source_RNA_sum.gt(0).to_numpy() & units.common_native_RNA_sum.gt(0).to_numpy() & np.isfinite(quality).all(axis=1)
    assert np.array_equal(measurable,units.RNA_measurement_available)
    assert np.isfinite(saved['raw'][measurable]).all() and np.isnan(saved['raw'][~measurable]).all()
    cases=data['cases'].sampleID.astype(str).tolist();train=data['train'];references=[];centers=[];coef=[]
    for vi,view in enumerate(VIEWS):
        native=pd.read_csv(inp/('private_'+view+'_pseudobulk.tsv.gz'),sep='\t',index_col=0);assert native.columns.tolist()==cases
        references.append(native.to_numpy()[:,train].sum(axis=1))
        norm=pd.read_csv(cal/('private_'+view+'_all_universe_TMM.tsv.gz'),sep='\t',dtype={'sampleID':str})
        centers.append([float(np.exp(np.log(norm[norm.mask_id.eq('mask_%03d'%m)].set_index('sampleID').loc[cases].raw_factor.to_numpy()[train]).mean())) for m in range(48)])
        w=data['model']['weights_'+str(vi)];sw=np.sqrt(data['model']['precision_'+str(vi)])
        q,r=np.linalg.qr(np.vstack([sw[:,None]*w,np.eye(10)]),mode='reduced')
        coef.append(solve_triangular(r,q[:len(w)].T*sw)/float(data['prepared'][vi]['training_view_RMS']))
    maps=pd.read_csv(cal/'private_all61_source_native_universe_map.tsv',sep='\t').set_index(['arm','dataset'])
    gates=pd.read_csv(cal/'private_all2440_native_RNA_universe_calibrations.tsv',sep='\t').set_index(['arm','dataset','view','factor'])
    manifest=pd.read_csv(count/'private_all61_count_archive_manifest.tsv',sep='\t');raw=np.full((6760,10),np.nan);gate=np.zeros((6760,10),bool)
    factors=saved['factors'].astype(str).tolist();assert factors==['NSCLC_MFA_20261016_F%02d'%k for k in range(1,11)]
    max_tmm=0.;max_raw=0.;normalizations=0
    checkpoint=root/'HNSCC_receiver_followup_job_control_20261016_v2'/Path(__file__).stem/'private_audit_checkpoint.json'
    for entry in manifest.itertuples():
        with np.load(count/entry.archive,allow_pickle=False) as z:counts=z['counts'];columns=z['common_native_columns'];rows=z['global_unit_indices']
        inverse=np.full(31831,-1,int);inverse[columns]=np.arange(len(columns));mi=int(maps.loc[(entry.arm,entry.dataset),'mask_index'])
        for vi,view in enumerate(VIEWS):
            all_ids=rows[rows%4==vi]
            for k,factor in enumerate(factors):gate[all_ids,k]=bool(gates.loc[(entry.arm,entry.dataset,view,factor),'native_RNA_universe_gate'])
            take=(rows%4==vi)&measurable[rows];ids=rows[take]
            if not len(ids):continue
            obs=counts[take];ref=references[vi][columns]
            f=np.array([calculate(row,ref,mode=3)[0] for row in obs]);effective=obs.sum(axis=1)*f/centers[vi][mi]
            observed=units.loc[ids,['raw_TMM_factor','native_discovery_centered_TMM_factor','effective_common_library_size']].astype(float).to_numpy()
            max_tmm=max(max_tmm,float(np.max(abs(f-observed[:,0]))))
            assert np.allclose(f,observed[:,0],rtol=1e-10,atol=1e-12)
            assert np.allclose(f/centers[vi][mi],observed[:,1],rtol=1e-10,atol=1e-12)
            assert np.allclose(effective,observed[:,2],rtol=1e-10,atol=1e-6)
            selected=data['prepared'][vi]['selected_native_columns'];present=inverse[selected]>=0
            value=np.log2(1+obs[:,inverse[selected[present]]].astype(float)*1e6/effective[:,None])
            raw[ids]=value@coef[vi][:,present].T;max_raw=max(max_raw,float(np.max(abs(raw[ids]-saved['raw'][ids]))));normalizations+=len(ids)
        checkpoint.write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),stage='independent_TMM_and_projection',completed_source_index=int(entry.source_index)))+'\n')
    assert np.allclose(raw,saved['raw'],atol=1e-8,rtol=1e-9,equal_nan=True) and np.array_equal(gate,saved['source_gene_gate'])
    fit=pd.read_csv(src/'private_all_source_quality_context_fit_support.tsv',sep='\t').set_index(['context_id','view'])
    assert fit.index.is_unique
    params={(r['context_id'],r['view']):r for r in json.loads((src/'private_source_quality_projection_parameters.json').read_text())}
    residual=np.full_like(raw,np.nan);supported=np.zeros(6760,bool);variable=np.zeros_like(gate)
    fits_checked=0;max_residual=0.
    for levels,g in units.groupby(CONTEXT,sort=True,dropna=False):
        context='QC_'+hashlib.sha256(json.dumps([str(v) for v in levels],ensure_ascii=False,separators=(',',':')).encode()).hexdigest()[:16]
        assert g.quality_context_id.eq(context).all()
        for view in VIEWS:
            ids=g.index[(g['view']==view)&g.RNA_measurement_available].to_numpy(int);observed=fit.loc[(context,view)];fits_checked+=1
            assert observed.libraries==len(ids)
            if not len(ids):assert observed.status=='NO_MEASURABLE_LIBRARY_POOLS';continue
            block=units.loc[ids,'resampling_block'].astype(str);weight=1/block.map(block.value_counts()).to_numpy(float)
            q=quality[ids];mean=(weight[:,None]*q).sum(axis=0)/weight.sum();sd=np.sqrt((weight[:,None]*(q-mean)**2).sum(axis=0)/weight.sum())
            design=np.column_stack([np.ones(len(ids)),(q-mean)/np.where(sd>0,sd,1)])
            sw=np.sqrt(weight);beta,_,rank,_=np.linalg.lstsq(sw[:,None]*design,sw[:,None]*raw[ids],rcond=1e-10)
            good=block.nunique()>=10 and len(ids)-rank>=5
            assert observed.resampling_blocks==block.nunique() and observed.quality_rank==rank
            assert observed.status==('SOURCE_QUALITY_PROJECTION_SUPPORTED' if good else 'INSUFFICIENT_SOURCE_QUALITY_FIT_SUPPORT')
            if not good:continue
            p=params[context,view];assert np.allclose(p['quality_mean'],mean) and np.allclose(p['quality_SD'],sd)
            residual[ids]=raw[ids]-design@beta;supported[ids]=True;variable[ids]=np.average(residual[ids]**2,weights=weight,axis=0)>1e-20
            max_residual=max(max_residual,float(np.max(abs(residual[ids]-saved['source_QC'][ids]))))
    assert np.allclose(residual,saved['source_QC'],atol=1e-8,rtol=1e-8,equal_nan=True)
    assert np.array_equal(supported,saved['quality_fit_supported']) and np.array_equal(variable,saved['variable_programme'])
    status=pd.read_csv(src/'private_all67600_programme_measurement_status.tsv.gz',sep='\t');assert len(status)==67600 and not status.duplicated(['unit_index','factor']).any()
    for k,factor in enumerate(factors):
        observed=status[status.factor.eq(factor)].set_index('unit_index').loc[range(6760)]
        reason=np.full(6760,'SOURCE_QC_PROGRAMME_INDEX_AVAILABLE',object)
        reason[~variable[:,k]]='PROGRAMME_RESIDUAL_CONSTANT_NOT_TESTABLE';reason[~supported]='INSUFFICIENT_SOURCE_QUALITY_FIT_SUPPORT'
        reason[~gate[:,k]]='SOURCE_GENE_MEASUREMENT_GATE_FAILED';reason[~measurable]='SOURCE_OR_COMMON_RNA_OR_QUALITY_UNMEASURABLE'
        reason[units.n_cells.to_numpy()<20]='CELL_POOL_BELOW20';reason[units.n_cells.to_numpy()==0]='CELL_POOL_ABSENT'
        assert np.array_equal(reason,observed.status)
        assert np.array_equal(gate[:,k],observed.source_gene_gate) and np.array_equal(supported,observed.source_quality_fit_supported) and np.array_equal(variable[:,k],observed.variable_programme)
    proof=src/'private_independent_fixed_programme_assay_audit.json';assert not proof.exists()
    report=dict(status='PASS_ALL61_EXTERNAL_TMM_QR_PROGRAMME_INDICES_SOURCE_QC_MODELS_AND67600_STATUSES',
        finished_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),producer_terminal_sha256=sha(tpath),
        external_library_pool_normalizations_checked=normalizations,source_quality_context_fits_checked=fits_checked,
        maximum_TMM_error=max_tmm,maximum_raw_index_error=max_raw,maximum_source_QC_index_error=max_residual,
        exact_native_factor_identity_or_clinical_validation=False,new_P=False)
    proof.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],finished_UTC=report['finished_UTC'],audit_sha256=sha(proof))))


if __name__=='__main__':main()
