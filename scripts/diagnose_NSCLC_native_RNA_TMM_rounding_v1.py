"""Numerical-only diagnosis of observed Python/R TMM disagreement, no effects."""
import json
import math
import os
import subprocess
import sys
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from NSCLC_native_MFA_clinical_inputs_v1 import sha,matched_quota


def calculate(obs,ref,mode=0):
    no=float(obs.sum());nr=float(ref.sum());valid=(obs>0)&(ref>0)
    o=obs[valid].astype(float);r=ref[valid].astype(float)
    scalar=lambda x:np.fromiter((math.log2(float(v)) for v in x),dtype=float,count=len(x))
    lr=(scalar if mode in [1,3] else np.log2)((o/no)/(r/nr))
    logger=scalar if mode in [2,3] else np.log2
    avg=(logger(o/no)+logger(r/nr))/2
    variance=(no-o)/no/o+(nr-r)/nr/r
    good=np.isfinite(lr)&np.isfinite(avg)&(avg>-1e10);lr=lr[good];avg=avg[good];variance=variance[good]
    if not len(lr) or np.max(abs(lr))<1e-6:return 1.,None
    n=len(lr);rl=rankdata(lr);ra=rankdata(avg)
    keep=(rl>=np.floor(n*.3)+1)&(rl<=n-np.floor(n*.3))&(ra>=np.floor(n*.05)+1)&(ra<=n-np.floor(n*.05))
    weight=1/variance[keep];f=np.sum(lr[keep]*weight)/np.sum(weight)
    return (2**float(f) if np.isfinite(f) else 1.),dict(lr=lr,average=avg,variance=variance,keep=keep)


R_DIAG=r'''
args<-commandArgs(trailingOnly=TRUE)
.libPaths(c(file.path(args[5],c('library_v7','library_v5','library_v2','library')),.libPaths()))
stopifnot(as.character(packageVersion('edgeR'))=='4.0.16')
counts<-as.matrix(read.delim(gzfile(args[1]),row.names=1,check.names=FALSE))
axis<-read.delim(args[2],check.names=FALSE)
masks<-as.matrix(read.delim(gzfile(args[3]),row.names=1,check.names=FALSE))
cc<-counts[masks[,38]==1,,drop=FALSE];ref<-rowSums(cc[,axis$partition=='discovery',drop=FALSE]);obs<-cc[,as.integer(args[4])]
no<-sum(obs);nr<-sum(ref);positive<-obs>0 & ref>0;o<-obs[positive];r<-ref[positive]
lr<-log2((o/no)/(r/nr));a<-(log2(o/no)+log2(r/nr))/2;v<-(no-o)/no/o+(nr-r)/nr/r
ok<-is.finite(lr)&is.finite(a)&a> -1e10;lr<-lr[ok];a<-a[ok];v<-v[ok];n<-length(lr)
keep<-rank(lr)>=floor(n*.3)+1 & rank(lr)<=n-floor(n*.3) & rank(a)>=floor(n*.05)+1 & rank(a)<=n-floor(n*.05)
write.table(data.frame(logratio=lr,average=a,variance=v,keep=keep),args[6],row.names=FALSE,quote=FALSE,sep='\t')
writeLines(deparse(get('.calcFactorTMM',asNamespace('edgeR'))),args[7])
'''


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);matched_quota(root)
    out=root/'NSCLC_native_MFA_external_RNA_TMM_rounding_diagnosis_20261016_v1';assert not out.exists();out.mkdir(mode=0o700)
    src=root/'NSCLC_native_MFA_external_RNA_universe_calibration_20261016_v1';inp=root/'NSCLC_native_MFA_inputs_20261016_v1'
    masks=pd.read_csv(src/'private_native_universe_masks.tsv.gz',sep='\t',index_col=0).to_numpy(bool)
    raw=pd.read_csv(inp/'private_author_Treg_pseudobulk.tsv.gz',sep='\t',index_col=0)
    counts=raw.to_numpy().T[:,masks[:,37]];axis=pd.read_csv(inp/'private_case_axis.tsv',sep='\t',dtype={'sampleID':str});ref=counts[axis.partition.eq('discovery')].sum(axis=0)
    norms=pd.read_csv(src/'private_author_Treg_all_universe_TMM.tsv.gz',sep='\t',dtype={'sampleID':str})
    norms=norms[norms.mask_id.eq('mask_037')].set_index('sampleID').loc[axis.sampleID]
    variants=[];results=[]
    for mode in range(4):
        factor=np.array([calculate(row,ref,mode)[0] for row in counts]);results.append(factor)
        error=np.abs(factor-norms.raw_factor.to_numpy())
        variants.append(dict(mode=mode,scalar_ratio_log2=mode in [1,3],scalar_abundance_log2=mode in [2,3],max_raw_factor_error=float(error.max()),discrepant_cases=int((error>1e-10).sum())))
    bad=np.flatnonzero(abs(results[0]-norms.raw_factor.to_numpy())>1e-10);assert len(bad)>0
    script=out/'private_R_exact_trimming_diagnosis.R';script.write_text(R_DIAG)
    rscript=Path(sys.executable).resolve().parents[2]/'scrna_r/bin/Rscript';runtime=root.parents[2]/'Project_v5/runs/20261014_DIALOGUE_method_pilot_v1'
    table=out/'private_first_numeric_disagreement_R_trimming.tsv';official=out/'private_official_edgeR_calcFactorTMM_source.R'
    with (out/'private_R_diagnosis.log').open('x') as log:
        subprocess.run([str(rscript),str(script),str(inp/'private_author_Treg_pseudobulk.tsv.gz'),str(inp/'private_case_axis.tsv'),str(src/'private_native_universe_masks.tsv.gz'),str(int(bad[0])+1),str(runtime),str(table),str(official)],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=120)
    exact=pd.read_csv(table,sep='\t');comparisons=[]
    for mode in range(4):
        factor,parts=calculate(counts[int(bad[0])],ref,mode)
        comparisons.append(dict(mode=mode,trimming_membership_differences=int((parts['keep']!=exact.keep).sum()),
            maximum_logratio_difference=float(np.max(abs(parts['lr']-exact.logratio))),maximum_average_difference=float(np.max(abs(parts['average']-exact.average)))))
    report=dict(status='COMPLETE_OBSERVED_TMM_DISAGREEMENT_NUMERICAL_DIAGNOSIS_NO_THRESHOLD_CHANGE',
        finished_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),variants=variants,first_numeric_case_comparisons=comparisons,
        private_R_source_sha256=sha(official),no_genes_cases_or_factors_selected_for_scientific_effect=True,
        original_R_production_or_physical_gene_gates_changed=False)
    p=out/'private_numerical_diagnosis.json';p.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))


if __name__=='__main__':main()
