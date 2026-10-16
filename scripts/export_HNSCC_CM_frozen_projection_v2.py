"""Export original scaling and nsNMF W/S without learning new parameters."""
import hashlib, json, os, subprocess
from datetime import datetime, timezone
from pathlib import Path
import numpy as np, pandas as pd

SCOPES = ['whole_TME', 'T_lineage_biopsy', 'whole_TME_primary_tumor', 'T_whole_TME_primary_tumor']

R_CODE = r'''
args<-commandArgs(trailingOnly=TRUE)
atlas<-args[[1]];out<-args[[2]];source_file<-args[[3]]
suppressPackageStartupMessages(library(NMF))
expressions<-parse(source_file)
hit<-vapply(expressions,function(x)is.call(x)&&is.symbol(x[[1]])&&as.character(x[[1]])%in%c('<-','=')&&is.symbol(x[[2]])&&as.character(x[[2]])=='freq_normalize',logical(1))
stopifnot(sum(hit)==1L);eval(expressions[[which(hit)]])
scopes<-readLines(file.path(out,'scopes.txt'))
idx<-jsonlite::fromJSON(file.path(atlas,'manifests','CM_recovery_v3_public.json'))$results
for(scope in scopes){
 folder<-file.path(atlas,'CM',scope)
 accepted<-file.path(folder,idx$discovery_attempt[match(scope,idx$scope)])
 od<-file.path(out,scope);dir.create(od,mode='0700')
 fq<-read.delim(file.path(folder,'inputs','frequency_private.tsv'),check.names=FALSE,row.names=1)
 norm<-freq_normalize(fq,normalize='minmax')
 saved<-readRDS(file.path(accepted,'NMF_discovery.rds'))
 W<-as.matrix(basis(saved));H<-as.matrix(coef(saved));model<-fit(saved)
 S<-if(is(model,'NMFns'))smoothing(model)else diag(ncol(W))
 dimnames(S)<-list(colnames(W),colnames(W))
 A<-W%*%S;F<-as.matrix(fitted(saved))
 oldW<-as.matrix(read.delim(file.path(accepted,'W_private.tsv'),row.names=1,check.names=FALSE))
 oldH<-as.matrix(read.delim(file.path(accepted,'H_patient_block_private.tsv'),row.names=1,check.names=FALSE))
 stopifnot(identical(dimnames(W),dimnames(oldW)),identical(dimnames(H),dimnames(oldH)))
 stopifnot(max(abs(W-oldW))<1e-9,max(abs(H-oldH))<1e-9,max(abs(A%*%H-F))<1e-8)
 raw<-as.matrix(fq[,-1]);lo<-apply(raw,1,min);span<-apply(raw,1,max)-lo
 stopifnot(all(span>0),max(abs((raw-lo)/span-norm))<1e-12)
 for(name in c('W','S','A','H','F','norm'))write.table(get(name),file.path(od,paste0(name,'_private.tsv')),sep='\t',quote=FALSE,col.names=NA)
 write.table(data.frame(node=rownames(raw),minimum=lo,range=span),file.path(od,'scaling_private.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
 jsonlite::write_json(list(scope=scope,model_class=class(model),NMF_version=as.character(packageVersion('NMF')),
  original_W_H_exact=TRUE,original_fitted_equals_W_S_H=TRUE,original_minmax_equals_python_formula=TRUE,
  smoothing_off_diagonal_nonzero=any(S[row(S)!=col(S)]>0)),file.path(od,'R_model_contract.json'),auto_unbox=TRUE,pretty=TRUE)
}
'''


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    project=Path(os.environ['PANCANCER_PROJECT_ROOT']);root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    atlas=project/'Project_v5/runs/20261013_immunotherapy_atlas_v1'
    prior=root/'HNSCC_new_CM_clinical_scope_20261016_v1/private_terminal.json'
    scope_audit=json.loads(prior.read_text());assert scope_audit['status']=='ALL15_NEW_ATLAS_CM_CLINICAL_SCOPE_AUDITED_NO_EFFECT_TEST'
    selected=[s for s in scope_audit['scopes'] if s['baseline_reviewed_author_patients']>=20]
    assert {s['scope'] for s in selected}==set(SCOPES) and sum(s['rank'] for s in selected)==30
    out=root/'HNSCC_new_CM_frozen_model_export_20261016_v2';assert not out.exists();out.mkdir()
    commands=[]
    for path in (atlas/'manifests').glob('CM_whole_TME*started_private.json'):
        command=json.loads(path.read_text()).get('command',[])
        if len(command)>3 and str(command[1]).endswith('run_immune_CM.R'):
            commands.append(command)
    assert commands
    Rscript=commands[0][0];source=Path(commands[0][3])/'R/CoVarNet.R';assert source.exists()
    assert all((Path(c[3])/'R/CoVarNet.R').read_bytes()==source.read_bytes() for c in commands)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(Path(__file__)),source_CoVarNet_R_sha256=sha(source),
        prior_scope_audit_sha256=sha(prior),scopes=SCOPES,primary_slots=30,
        primary_projection='NNLS onto original effective nsNMF basis A=W S, verified against original fitted model. W/S fixed, no parameter refit.',
        clarification_before_effects='Raw W alone omits original nsNMF smoothing; keep model-faithful coefficients as primary and raw-W coefficients only as descriptive sensitivity.',
        scaling='Original full-discovery per-node min/range, not baseline-refitted; clinical use remains internal exploratory.',
        no_clinical_outcome_read_or_model_fit=True)
    (out/'plan_before_exports_or_clinical_effects.json').write_text(json.dumps(plan,indent=2)+'\n')
    (out/'scopes.txt').write_text('\n'.join(SCOPES)+'\n');script=out/'export_original_R_model.R';script.write_text(R_CODE)
    env=os.environ.copy();env['OMP_NUM_THREADS']='2';env['OPENBLAS_NUM_THREADS']='2';env['MKL_NUM_THREADS']='2'
    run=subprocess.run([str(Rscript),str(script),str(atlas),str(out),str(source)],capture_output=True,text=True,env=env,timeout=120)
    (out/'R_export_stdout.log').write_text(run.stdout);(out/'R_export_stderr.log').write_text(run.stderr)
    assert run.returncode==0,run.stderr[-2500:]
    summaries=[]
    for scope in SCOPES:
        od=out/scope
        W=pd.read_csv(od/'W_private.tsv',sep='\t',index_col=0);S=pd.read_csv(od/'S_private.tsv',sep='\t',index_col=0)
        A=pd.read_csv(od/'A_private.tsv',sep='\t',index_col=0);H=pd.read_csv(od/'H_private.tsv',sep='\t',index_col=0)
        F=pd.read_csv(od/'F_private.tsv',sep='\t',index_col=0);norm=pd.read_csv(od/'norm_private.tsv',sep='\t',index_col=0)
        frequency=pd.read_csv(atlas/'CM'/scope/'inputs/frequency_private.tsv',sep='\t',index_col=0)
        scale=pd.read_csv(od/'scaling_private.tsv',sep='\t',index_col=0)
        assert W.index.equals(frequency.index) and A.index.equals(W.index)
        assert list(W.columns)==list(H.index)==list(A.columns)==list(S.index)==list(S.columns)
        assert np.allclose(W.to_numpy()@S.to_numpy(),A,atol=1e-12,rtol=1e-12)
        assert np.allclose(A.to_numpy()@H.to_numpy(),F,atol=1e-8,rtol=1e-10)
        expected=(frequency.iloc[:,1:].to_numpy()-scale.minimum.to_numpy()[:,None])/scale['range'].to_numpy()[:,None]
        assert np.allclose(expected,norm,atol=1e-12,rtol=1e-12)
        contract=json.loads((od/'R_model_contract.json').read_text())
        summaries.append(dict(scope=scope,rank=W.shape[1],nodes=W.shape[0],original_model_contract=contract,
            effective_basis_rank=int(np.linalg.matrix_rank(A)),effective_basis_condition=float(np.linalg.cond(A)),
            effective_basis_sha256=sha(od/'A_private.tsv'),normalization_sha256=sha(od/'norm_private.tsv')))
    report=dict(status='ALL4_ORIGINAL_MINMAX_AND_NSNMF_W_S_H_PARITY_EXPORTED_NO_CLINICAL_EFFECTS',
        finished_UTC=datetime.now(timezone.utc).isoformat(),scopes=summaries,original_models_changed=False)
    (out/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()
