"""Fixed lung gene vectors, native TMM and source-own quality projection.

Explicit assay extension, not exact original latent-factor replication. No
external outcome, association, feature selection or clinical model is used.
"""
import hashlib
import io
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_external_lineage_inputs_v2 import load as metadata_inputs
from NSCLC_native_MFA_interpretation_inputs_v1 import load as native_inputs,sha,VIEWS

ASSAY='external_nativeCount_frozenA_sourceQC_v1'
CONTEXT=['arm','dataset','cancer_label','sample_class','material_class','capture_scope','ex_vivo']
QFIELDS=['logN','mean_logRNA','mean_loggenes','symbol_pct_mito','symbol_pct_ribo']
R_CODE=r'''
args<-commandArgs(trailingOnly=TRUE)
.libPaths(c(file.path(args[1],c('library_v7','library_v5','library_v2','library')),.libPaths()))
stopifnot(as.character(packageVersion('edgeR'))=='4.0.16')
ng<-as.integer(args[2]);nu<-as.integer(args[3]);center<-as.numeric(args[4])
con<-file('stdin','rb');x<-readBin(con,what='double',n=ng*(nu+1),size=8,endian='little');close(con)
stopifnot(length(x)==ng*(nu+1));x<-matrix(x,nrow=ng,ncol=nu+1)
ref<-x[,1];refsize<-sum(ref);stopifnot(refsize>0,center>0,all(x>=0),all(is.finite(x)))
f<-vapply(seq_len(nu),function(i){
  obs<-x[,i+1];lib<-sum(obs);stopifnot(lib>0)
  v<-edgeR::calcNormFactors(cbind(ref,obs),lib.size=c(refsize,lib),refColumn=1,
      method='TMM',logratioTrim=.3,sumTrim=.05,doWeighting=TRUE,Acutoff=-1e10)
  v[2]/v[1]
},numeric(1))
lib<-colSums(x[,-1,drop=FALSE])
write.table(data.frame(position=seq_len(nu),raw_factor=f,factor=f/center,
  library_size=lib,effective_library_size=lib*f/center),stdout(),sep='\t',quote=FALSE,row.names=FALSE)
'''


def resource(root,cap=16*1024**2):
    fs=subprocess.check_output(['df','-P',str(root.resolve())],text=True).splitlines()[-1].split()[0]
    q=subprocess.check_output(['quota','-g','-v','-w','USER002'],text=True)
    rows=[s.split() for s in q.splitlines() if s.split() and s.split()[0]==fs];assert len(rows)==1
    free=(int(rows[0][3])-int(rows[0][1].rstrip('*')))*1024
    assert min(free,shutil.disk_usage(root).free)>500_000_000_000+1073741824+cap
    return dict(group_remaining_bytes=free,reserve_bytes=500_000_000_000,buffer_bytes=1073741824,output_cap_bytes=cap)


def normalize(reference,counts,center,rscript,helper,runtime,log):
    assert reference.ndim==1 and counts.ndim==2 and counts.shape[1]==len(reference)
    assert (counts.sum(axis=1)>0).all() and counts.max()<2**53 and reference.max()<2**53
    data=np.column_stack([reference,counts.T]).astype('<f8',order='F').tobytes(order='F')
    result=subprocess.run([str(rscript),str(helper),str(runtime),str(len(reference)),str(len(counts)),format(center,'.17g')],
        input=data,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True,timeout=300)
    with log.open('ab') as h:h.write(result.stderr)
    table=pd.read_csv(io.StringIO(result.stdout.decode()),sep='\t')
    assert table.position.tolist()==list(range(1,len(counts)+1))
    assert np.array_equal(table.library_size,counts.sum(axis=1))
    assert np.isfinite(table.effective_library_size).all() and table.effective_library_size.gt(0).all()
    return table


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget=resource(root)
    scope=root/'NSCLC_native_MFA_external_quality_domain_review_20261016_v1/PRIVATE_NEXT_FIXED_GENE_PROGRAMME_ASSAY_SCOPE.json'
    assert sha(scope)=='55df467a5e42db1b0a5fcb84573ab84672b6d5e9fdfff62e72109936a50f96f9'
    src=root/'NSCLC_native_MFA_external_frozen_role_count_inputs_20261016_v2'
    count_proof=src/'private_independent_native_count_input_audit.json'
    assert sha(count_proof)=='7bd028cc1c89ac7324ac25c65042b5b184cb377162a38f9544caf7bb24048c04'
    tpath=src/'private_terminal.json';terminal=json.loads(tpath.read_text())
    assert json.loads(count_proof.read_text())['producer_terminal_sha256']==sha(tpath)
    for name,item in terminal['files'].items():assert sha(src/name)==item['sha256']
    cal=root/'NSCLC_native_MFA_external_RNA_universe_calibration_20261016_v1'
    proof=cal/'private_independent_native_RNA_universe_audit.json'
    assert sha(proof)=='f9822656c7be74fa7953b04c8807b5da718bc2084bf3a946e7943504ed31d5b3'
    assert json.loads(proof.read_text())['producer_terminal_sha256']==sha(cal/'private_terminal.json')
    data=native_inputs();meta=metadata_inputs();assert data['root']==meta['root']==root
    out=root/'NSCLC_external_fixed_gene_programme_assay_20261016_v1';assert not out.exists();out.mkdir(mode=0o700)
    helper=out/'frozen_binary_reference_TMM.R';helper.write_text(R_CODE)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),assay_id=ASSAY,
        frozen_scope_sha256=sha(scope),count_audit_sha256=sha(count_proof),RNA_universe_audit_sha256=sha(proof),budget=budget,
        model_sha256=sha(data['model_path']),scope='All61 source-arm/6760 frozen compatible pool positions/all10 fixed factors; no outcome or association selection.',
        original_factor_replication=False,source_quality_adaptation=True,clinical_or_association_model=False,
        quality='Own5 observed quality variables plus intercept. Full weighted SVD span threshold maxSingular*1e-10, weighted source centering/SD. Patient block totalweight1; unknown library fallback retained. Require>=10blocks and nLibraries-rank>=5.',
        context=CONTEXT,phase_or_treatment_centered=False,parent_fraction_bridged_or_used=False,
        binary_check='Before external indices, original four views xfirst2 frozen cases against accepted TMM factors/effective sizes; reference and geometric constant original145 only.',
        no_new_data_download_or_annotation_changes=True)
    (out/'plan_before_any_external_fixed_programme_measurement.json').write_text(json.dumps(plan,indent=2)+'\n')
    units=pd.read_csv(src/'private_all6760_native_count_units_and_quality.tsv.gz',sep='\t')
    original_columns=units.columns.tolist();lib=meta['lib']
    fields=['arm','dataset','library_id','patient_block','reviewed_patient','resampling_block','cancer_label','tissue','sample_class','timepoint','treatment_arm','capture_scope','ex_vivo','material_class','phase_class']
    units=units.merge(lib[fields],on=['arm','dataset','library_id'],how='left',sort=False,validate='many_to_one')
    assert units.unit_index.tolist()==list(range(6760)) and units.resampling_block.notna().all()
    quality=np.column_stack([np.log1p(units.n_cells),units[['mean_log_full_source_RNA','mean_log_full_source_detected_genes','symbol_pct_mito_RNA','symbol_pct_ribo_RNA']].to_numpy(float)])
    measurable=units.n_cells.ge(20).to_numpy() & units.full_source_RNA_sum.gt(0).to_numpy() & units.common_native_RNA_sum.gt(0).to_numpy() & np.isfinite(quality).all(axis=1)
    raw=np.full((6760,10),np.nan);residual=raw.copy();normalizers=np.full((6760,3),np.nan);source_gate=np.zeros((6760,10),bool)
    factors=['NSCLC_MFA_20261016_F%02d'%k for k in range(1,11)]
    gates=pd.read_csv(cal/'private_all2440_native_RNA_universe_calibrations.tsv',sep='\t').set_index(['arm','dataset','view','factor'])
    source_masks=pd.read_csv(cal/'private_all61_source_native_universe_map.tsv',sep='\t').set_index(['arm','dataset'])
    inp=root/'NSCLC_native_MFA_inputs_20261016_v1';cases=data['cases'].sampleID.astype(str).tolist();train=data['train']
    references=[];projections=[];centers=[];checks=[]
    rscript=Path(sys.executable).resolve().parents[2]/'scrna_r/bin/Rscript';runtime=root.parents[2]/'Project_v5/runs/20261014_DIALOGUE_method_pilot_v1'
    stderr=out/'private_binary_TMM_stderr.log'
    for vi,view in enumerate(VIEWS):
        counts=pd.read_csv(inp/('private_'+view+'_pseudobulk.tsv.gz'),sep='\t',index_col=0)
        assert counts.columns.tolist()==cases and counts.index.astype(str).tolist()==data['genes'].tolist()
        ref=counts.to_numpy()[:,train].sum(axis=1);references.append(ref)
        norms=pd.read_csv(cal/('private_'+view+'_all_universe_TMM.tsv.gz'),sep='\t',dtype={'sampleID':str})
        ctr=[]
        for mi in range(48):
            n=norms[norms.mask_id.eq('mask_%03d'%mi)].set_index('sampleID').loc[cases]
            ctr.append(float(np.exp(np.log(n.raw_factor.to_numpy()[train]).mean())))
        centers.append(ctr)
        baseline=pd.read_csv(inp/('private_'+view+'_TMM.tsv'),sep='\t',dtype={'sampleID':str}).set_index('sampleID').loc[cases]
        result=normalize(ref,counts.iloc[:,:2].to_numpy().T,ctr[0],rscript,helper,runtime,stderr)
        assert np.allclose(result.raw_factor,baseline.raw_factor.iloc[:2],rtol=1e-10,atol=1e-12)
        assert np.allclose(result.effective_library_size,baseline.effective_library_size.iloc[:2],rtol=1e-10,atol=1e-6)
        checks.append(dict(view=view,fixed_original_cases_checked=2,passed=True))
        w=data['model']['weights_'+str(vi)];tau=data['model']['precision_'+str(vi)]
        a=np.linalg.solve(np.eye(10)+w.T@(tau[:,None]*w),w.T*tau)/float(data['prepared'][vi]['training_view_RMS'])
        projections.append(a)
    (out/'private_all8_original_binary_TMM_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
    manifest=pd.read_csv(src/'private_all61_count_archive_manifest.tsv',sep='\t')
    for record in manifest.itertuples():
        path=src/record.archive;assert sha(path)==record.sha256
        with np.load(path,allow_pickle=False) as z:counts=z['counts'];columns=z['common_native_columns'];rows=z['global_unit_indices']
        physical=np.full(31831,-1,int);physical[columns]=np.arange(len(columns));mi=int(source_masks.loc[(record.arm,record.dataset),'mask_index'])
        for vi,view in enumerate(VIEWS):
            all_ids=rows[rows%4==vi]
            for k,factor in enumerate(factors):source_gate[all_ids,k]=bool(gates.loc[(record.arm,record.dataset,view,factor),'native_RNA_universe_gate'])
            take=(rows%4==vi)&measurable[rows];ids=rows[take]
            if not len(ids):continue
            obs=counts[take];norm=normalize(references[vi][columns],obs,centers[vi][mi],rscript,helper,runtime,stderr)
            selected=data['prepared'][vi]['selected_native_columns'];available=physical[selected]>=0
            value=np.log2(1+obs[:,physical[selected[available]]].astype(float)*1e6/norm.effective_library_size.to_numpy()[:,None])
            raw[ids]=value@projections[vi][:,available].T
            normalizers[ids]=norm[['raw_factor','factor','effective_library_size']].to_numpy(float)
        (out/'private_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),stage='fixed_gene_projection',completed_source_index=int(record.source_index),total_sources=61))+'\n')
    assert np.isfinite(raw[measurable]).all() and np.isnan(raw[~measurable]).all()
    units['assay_id']=ASSAY;units['RNA_measurement_available']=measurable
    for j,name in enumerate(QFIELDS):units['source_quality_'+name]=quality[:,j]
    for j,name in enumerate(['raw_TMM_factor','native_discovery_centered_TMM_factor','effective_common_library_size']):units[name]=normalizers[:,j]
    units['quality_context_id']='';qc_pass=np.zeros(6760,bool);variable=np.zeros((6760,10),bool);fits=[]
    parameters=[]
    for levels,g in units.groupby(CONTEXT,sort=True,dropna=False):
        key=json.dumps([str(v) for v in levels],ensure_ascii=False,separators=(',',':'))
        context='QC_'+hashlib.sha256(key.encode()).hexdigest()[:16]
        units.loc[g.index,'quality_context_id']=context
        for vi,view in enumerate(VIEWS):
            ids=g.index[(g['view']==view)&g.RNA_measurement_available].to_numpy(int)
            if not len(ids):
                fits.append(dict(context_id=context,view=view,libraries=0,resampling_blocks=0,quality_rank=0,status='NO_MEASURABLE_LIBRARY_POOLS'))
                continue
            block=units.loc[ids,'resampling_block'].astype(str);weights=1/block.map(block.value_counts()).to_numpy(float)
            q=quality[ids];mean=np.average(q,weights=weights,axis=0);sd=np.sqrt(np.average((q-mean)**2,weights=weights,axis=0))
            design=np.column_stack([np.ones(len(ids)),(q-mean)/np.where(sd>0,sd,1)])
            sw=np.sqrt(weights);u,s,vt=np.linalg.svd(sw[:,None]*design,full_matrices=False);rank=int((s>s[0]*1e-10).sum())
            blocks=block.nunique();supported=blocks>=10 and len(ids)-rank>=5
            status='SOURCE_QUALITY_PROJECTION_SUPPORTED' if supported else 'INSUFFICIENT_SOURCE_QUALITY_FIT_SUPPORT'
            fits.append(dict(context_id=context,view=view,libraries=len(ids),resampling_blocks=int(blocks),quality_rank=rank,status=status))
            if not supported:continue
            b=vt[:rank].T@((u[:,:rank].T@(sw[:,None]*raw[ids]))/s[:rank,None])
            residual[ids]=raw[ids]-design@b;qc_pass[ids]=True
            var=np.average(residual[ids]**2,weights=weights,axis=0);variable[ids]=var>1e-20
            parameters.append(dict(context_id=context,view=view,quality_mean=mean.tolist(),quality_SD=sd.tolist(),rank=rank,coefficients=b.tolist(),blocks=int(blocks),libraries=len(ids)))
    statuses=[]
    for i in range(6760):
        for k,factor in enumerate(factors):
            reason='SOURCE_QC_PROGRAMME_INDEX_AVAILABLE'
            if units.n_cells.iloc[i]==0:reason='CELL_POOL_ABSENT'
            elif units.n_cells.iloc[i]<20:reason='CELL_POOL_BELOW20'
            elif not measurable[i]:reason='SOURCE_OR_COMMON_RNA_OR_QUALITY_UNMEASURABLE'
            elif not source_gate[i,k]:reason='SOURCE_GENE_MEASUREMENT_GATE_FAILED'
            elif not qc_pass[i]:reason='INSUFFICIENT_SOURCE_QUALITY_FIT_SUPPORT'
            elif not variable[i,k]:reason='PROGRAMME_RESIDUAL_CONSTANT_NOT_TESTABLE'
            statuses.append(dict(unit_index=i,factor=factor,source_gene_gate=bool(source_gate[i,k]),source_quality_fit_supported=bool(qc_pass[i]),variable_programme=bool(variable[i,k]),status=reason))
    np.savez_compressed(out/'private_all6760_fixed_programme_indices.npz',raw=raw,source_QC=residual,source_gene_gate=source_gate,quality_fit_supported=qc_pass,variable_programme=variable,factors=np.array(factors),quality=quality)
    units.to_csv(out/'private_all6760_programme_units_metadata_quality_TMM.tsv.gz',sep='\t',index=False,compression='gzip')
    pd.DataFrame(statuses).to_csv(out/'private_all67600_programme_measurement_status.tsv.gz',sep='\t',index=False,compression='gzip')
    pd.DataFrame(fits).to_csv(out/'private_all_source_quality_context_fit_support.tsv',sep='\t',index=False)
    (out/'private_source_quality_projection_parameters.json').write_text(json.dumps(parameters,separators=(',',':'))+'\n')
    terminal=dict(status='COMPLETE_ALL61_SOURCE_FIXED_GENE_PROGRAMME_ASSAY_INDICES_AUDIT_PENDING_NO_ASSOCIATIONS',
        finished_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),assay_id=ASSAY,
        files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()},
        exact_original_latent_factor_or_identity_replication=False,new_P_or_clinical_model=False)
    assert sum(item['bytes'] for item in terminal['files'].values())<budget['output_cap_bytes']
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'],finished_UTC=terminal['finished_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
