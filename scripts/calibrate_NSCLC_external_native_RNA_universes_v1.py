"""Whole-count gene-universe/TMM perturbation on original frozen cases.

No external RNA, clinical label, factor refit or learned feature export. Native
quality variables are held fixed to isolate ONLY gene-universe normalization.
"""
import hashlib
import json
import os
import subprocess
import sys
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from NSCLC_native_MFA_interpretation_inputs_v1 import load,sha,VIEWS
from NSCLC_native_MFA_clinical_inputs_v1 import matched_quota
from external_functional_coverage_inputs_v1 import read_inputs

R_MASK_TMM=r'''
args<-commandArgs(trailingOnly=TRUE)
.libPaths(c(file.path(args[4],c('library_v7','library_v5','library_v2','library')),.libPaths()))
stopifnot(as.character(packageVersion('edgeR'))=='4.0.16')
counts<-as.matrix(read.delim(gzfile(args[1]),row.names=1,check.names=FALSE))
axis<-read.delim(args[2],check.names=FALSE)
masks<-as.matrix(read.delim(gzfile(args[3]),row.names=1,check.names=FALSE))
stopifnot(identical(rownames(counts),rownames(masks)),identical(colnames(counts),as.character(axis$sampleID)))
train<-axis$partition=='discovery'
con<-gzfile(args[5],'wt')
writeLines('mask_id\tsampleID\traw_factor\tfactor\tlibrary_size\teffective_library_size',con)
for(j in seq_len(ncol(masks))){
  take<-masks[,j]==1
  if(sum(take)<2)stop('Insufficient universe for TMM')
  cc<-counts[take,,drop=FALSE]
  ref<-rowSums(cc[,train,drop=FALSE]);refsize<-sum(ref);lib<-colSums(cc)
  stopifnot(all(lib>0),refsize>0)
  factor<-vapply(seq_len(ncol(cc)),function(i){
    v<-edgeR::calcNormFactors(cbind(ref,cc[,i]),lib.size=c(refsize,lib[i]),method='TMM',refColumn=1,
        logratioTrim=.3,sumTrim=.05,doWeighting=TRUE,Acutoff=-1e10)
    v[2]/v[1]
  },numeric(1))
  normalizer<-exp(mean(log(factor[train])))
  result<-data.frame(mask_id=colnames(masks)[j],sampleID=colnames(cc),raw_factor=factor,
    factor=factor/normalizer,library_size=lib,effective_library_size=lib*factor/normalizer)
  write.table(result,con,row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t')
  cat(sprintf('Completed_native_mask %d/%d\n',j,ncol(masks)));flush.console()
}
close(con)
'''


def main():
    data=load();root=data['root'];matched_quota(root)
    question=root/'NSCLC_native_MFA_external_opportunity_review_20261016_v1/PRIVATE_NEXT_NATIVE_RNA_UNIVERSE_CALIBRATION_AND_IDENTITY_SCOPE.json'
    assert sha(question)=='cf6be4eec2749592d37d730db65847005bb8811633275f6266f8923a514ae481'
    out=root/'NSCLC_native_MFA_external_RNA_universe_calibration_20261016_v1'
    assert not out.exists();out.mkdir(mode=0o700)
    inp=root/'NSCLC_native_MFA_inputs_20261016_v1';axes=read_inputs(root)['source_axes']
    genes=data['genes'];train=data['train'];cases=data['cases'].sampleID.astype(str).tolist()
    masks=[np.ones(len(genes),bool)];mask_keys={hashlib.sha256(masks[0].tobytes()).hexdigest():0};source_map=[]
    for (arm,source),axis in sorted(axes.items()):
        positions=[next(iter(axis['lookup'][g])) if g in axis['lookup'] and len(axis['lookup'][g])==1 else None for g in genes]
        multiplicity=Counter(c for c in positions if c is not None)
        mask=np.array([c is not None and multiplicity[c]==1 for c in positions],bool)
        digest=hashlib.sha256(mask.tobytes()).hexdigest()
        if digest not in mask_keys:mask_keys[digest]=len(masks);masks.append(mask)
        source_map.append(dict(arm=arm,dataset=source,mask_index=mask_keys[digest],measured_native_universe_genes=int(mask.sum()),mask_sha256=digest))
    masks=np.array(masks);assert len(source_map)==61
    pd.DataFrame(masks.T.astype(int),index=genes,columns=['mask_%03d'%i for i in range(len(masks))]).to_csv(out/'private_native_universe_masks.tsv.gz',sep='\t',compression='gzip',index_label='gene')
    pd.DataFrame(source_map).to_csv(out/'private_all61_source_native_universe_map.tsv',sep='\t',index=False)
    (out/'frozen_all_universe_TMM.R').write_text(R_MASK_TMM)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),next_question_sha256=sha(question),
        model_sha256=sha(data['model_path']),source_mask_sha256=sha(out/'private_native_universe_masks.tsv.gz'),
        scope='All61 sources x4 views x10 factors, original208 cases and original145/63 split, plus full-universe positive replay.',
        normalization='Restrict original count pseudobulk AND discovery pooled TMM reference to uniquely mapped common native genes. Official original edgeR4.0.16 parameters; geometric centering on145 discovery only.',
        scoring='Use original per-gene quality regression coefficients, original quality design, original RMS and fixed own-view projection. Omit unavailable selected features. No normalization/weight/quality refit or phenotype used.',
        gates='Original physical and prepared-X mask gates retained; BOTH original partitions Pearson>=.95 and RMSE<=.25 full-score discovery SD. Full-universe TMM and scores must replay saved outputs.',
        limits='Native quality variables unchanged to isolate count-universe effects ONLY. This does not certify external platform/quality/identity, nor is it external expression replication.',
        new_external_expression_or_P=False,new_download=False)
    (out/'plan_before_native_count_universe_calibration.json').write_text(json.dumps(plan,indent=2)+'\n')
    prior=pd.read_csv(root/'NSCLC_native_MFA_external_gene_mask_calibration_20261016_v1/private_all2440_source_mask_calibration.tsv',sep='\t').set_index(['arm','dataset','view','factor'])
    runtime=root.parents[2]/'Project_v5/runs/20261014_DIALOGUE_method_pilot_v1'
    rscript=Path(sys.executable).resolve().parents[2]/'scrna_r/bin/Rscript'
    records=[];all_scores=[];normalizer_max=0.;baseline_max=0.
    for vi,view in enumerate(VIEWS):
        counts_path=inp/('private_'+view+'_pseudobulk.tsv.gz')
        original_terminal=json.loads((inp/'private_terminal.json').read_text());assert sha(counts_path)==original_terminal['files'][counts_path.name]
        norm_path=out/('private_'+view+'_all_universe_TMM.tsv.gz')
        with (out/('private_'+view+'_normalization.log')).open('x') as log:
            subprocess.run([str(rscript),str(out/'frozen_all_universe_TMM.R'),str(counts_path),str(inp/'private_case_axis.tsv'),str(out/'private_native_universe_masks.tsv.gz'),str(runtime),str(norm_path)],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=7200)
        norm=pd.read_csv(norm_path,sep='\t',dtype={'sampleID':str})
        assert len(norm)==len(masks)*len(cases)
        original_norm=pd.read_csv(inp/('private_'+view+'_TMM.tsv'),sep='\t',dtype={'sampleID':str}).set_index('sampleID').loc[cases]
        base_norm=norm[norm.mask_id.eq('mask_000')].set_index('sampleID').loc[cases]
        for col in ['raw_factor','factor','library_size','effective_library_size']:
            assert np.allclose(base_norm[col],original_norm[col],rtol=1e-10,atol=1e-8)
            normalizer_max=max(normalizer_max,float(np.max(abs(base_norm[col]-original_norm[col]))))
        p=data['prepared'][vi];selected=p['selected_native_columns'];w=data['model']['weights_'+str(vi)];tau=data['model']['precision_'+str(vi)]
        coef=np.linalg.solve(np.eye(10)+w.T@(tau[:,None]*w),w.T*tau)
        fixed_quality=p['design']@p['regression_coefficients'];rms=float(p['training_view_RMS'])
        raw=pd.read_csv(counts_path,sep='\t',index_col=0);assert raw.columns.tolist()==cases and raw.index.astype(str).tolist()==genes.tolist()
        raw_selected=raw.iloc[selected].to_numpy().T;full=data['model']['own_view_scores_'+str(vi)];sd=full[train].std(axis=0,ddof=1)
        view_scores=[]
        for mi,mask in enumerate(masks):
            normalization=norm[norm.mask_id.eq('mask_%03d'%mi)].set_index('sampleID').loc[cases]
            scale=normalization.effective_library_size.to_numpy(float)
            assert np.isfinite(scale).all() and (scale>0).all()
            take=mask[selected]
            value=np.log2(1+raw_selected[:,take].astype(float)*1e6/scale[:,None])
            score=((value-fixed_quality[:,take])/rms)@coef[:,take].T
            view_scores.append(score)
            if mi==0:
                error=float(np.max(abs(score-full)));assert error<1e-8;baseline_max=max(baseline_max,error)
        all_scores.append(np.array(view_scores))
        for src in source_map:
            partial=view_scores[src['mask_index']]
            for k in range(10):
                factor='NSCLC_MFA_20261016_F%02d'%(k+1);old=prior.loc[(src['arm'],src['dataset'],view,factor)]
                row=dict(src,view=view,factor=factor,original_physical_gate=bool(old.physical_gene_gate),original_prepared_X_mask_gate=bool(old.native_mask_calibration_pass))
                for name,take in [('discovery',train),('internal_holdout',~train)]:
                    x=full[take,k];y=partial[take,k];delta=y-x
                    row[name+'_Pearson']=float(np.corrcoef(x,y)[0,1]) if y.std()>1e-12 else np.nan
                    row[name+'_Spearman']=float(spearmanr(x,y).statistic) if y.std()>1e-12 else np.nan
                    row[name+'_RMSE_SD']=float(np.sqrt(np.mean(delta**2))/sd[k])
                    row[name+'_mean_shift_SD']=float(delta.mean()/sd[k])
                row['native_RNA_universe_gate']=bool(row['original_prepared_X_mask_gate'] and all(row[name+'_Pearson']>=.95 and row[name+'_RMSE_SD']<=.25 for name in ['discovery','internal_holdout']))
                records.append(row)
        (out/'private_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),completed_view=view,completed_view_count=vi+1))+'\n')
    assert len(records)==2440
    pd.DataFrame(records).to_csv(out/'private_all2440_native_RNA_universe_calibrations.tsv',sep='\t',index=False)
    np.savez_compressed(out/'private_all_universe_native_scores.npz',scores=np.array(all_scores),masks=masks)
    terminal=dict(status='COMPLETE_ALL61_SOURCE_NATIVE_RNA_UNIVERSES_TMM_AND2440_FACTOR_CALIBRATIONS_AUDIT_PENDING',
        finished_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),full_universe_max_score_error=baseline_max,
        full_universe_max_TMM_numeric_error=normalizer_max,files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()},
        external_expression_or_clinical_validation=False)
    assert sum(v['bytes'] for v in terminal['files'].values())<268435456
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'],finished_UTC=terminal['finished_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
