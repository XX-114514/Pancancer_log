"""Training-only pseudobulk normalization, technical regression and feature selection."""
import gzip,hashlib,json,os,re,subprocess,sys
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd

R_NORMALIZE=r'''
args<-commandArgs(trailingOnly=TRUE)
.libPaths(c(file.path(args[3],c('library_v7','library_v5','library_v2','library')),.libPaths()))
stopifnot(as.character(packageVersion('edgeR'))=='4.0.16')
counts<-as.matrix(read.delim(gzfile(args[1]),row.names=1,check.names=FALSE))
axis<-read.delim(args[2],check.names=FALSE)
stopifnot(identical(colnames(counts),as.character(axis$sampleID)))
train<-axis$partition=='discovery'
ref<-rowSums(counts[,train,drop=FALSE]);refsize<-sum(ref)
lib<-colSums(counts)
factor<-vapply(seq_len(ncol(counts)),function(i){
  v<-edgeR::calcNormFactors(cbind(ref,counts[,i]),lib.size=c(refsize,lib[i]),
       method='TMM',refColumn=1,logratioTrim=.3,sumTrim=.05,doWeighting=TRUE,Acutoff=-1e10)
  v[2]/v[1]
},numeric(1))
normalizer<-exp(mean(log(factor[train])))
write.table(data.frame(sampleID=colnames(counts),raw_factor=factor,
    factor=factor/normalizer,library_size=lib,effective_library_size=lib*factor/normalizer),
    args[4],row.names=FALSE,quote=FALSE,sep='\t')
'''


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as stream:
        for block in iter(lambda:stream.read(16777216),b''):h.update(block)
    return h.hexdigest()


def quota_free(root):
    fs=subprocess.check_output(['df','-P',str(root.resolve())],text=True).splitlines()[-1].split()[0]
    text=subprocess.check_output(['quota','-g','-v','-w','USER002'],text=True)
    rows=[s.split() for s in text.splitlines() if s.split() and s.split()[0]==fs];assert len(rows)==1
    free=(int(rows[0][3])-int(rows[0][1].rstrip('*')))*1024
    assert free>500_000_000_000+1_073_741_824+268_435_456
    return free


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);qual=root/'NSCLC_native_multicell_discovery_preflight_20261016_v1'
    proof=json.loads((qual/'private_terminal.json').read_text());assert proof['eligible_for_next_input_phase']
    for name,h in proof['files'].items():assert sha(qual/name)==h
    previous=json.loads((qual/'plan_before_any_native_discovery_effect.json').read_text())
    roles=pd.read_csv(qual/'private_all243_author_case_therapy_partition.tsv',sep='\t',keep_default_na=False,dtype={'sampleID':str})
    axis=roles[roles.discovery_role.isin(['discovery','internal_holdout'])].sort_values('sampleID').copy();axis=axis.rename(columns={'discovery_role':'partition'})
    cases=axis.sampleID.tolist();train=axis.partition.eq('discovery').to_numpy();assert len(cases)==208 and train.sum()==145 and (~train).sum()==63
    known=['Tislelizumab','Pembrolizumab','Camrelizumab','Sintilimab','Nivolumab','Toripalimab']
    axis['therapy_group']=[next(iter(found)) if len(found:={x for x in known if x in s})==1 else 'MULTI_OR_UNSPECIFIED' for s in axis.therapy_original]
    axis['holdout_block']=axis.cancer_type+'::'+axis.therapy_group
    native=root/'NSCLC_full_native_gene_units_20261016_v1';native_proof=root/'NSCLC_full_native_gene_unit_audit_20261016_v1/private_terminal.json'
    accepted=json.loads(native_proof.read_text());assert accepted['producer_terminal_sha256']==sha(native/'private_terminal.json')
    units=pd.read_csv(native/'private_all_native_units_and_quality.tsv',sep='\t',dtype={'sampleID':str});genes=pd.read_csv(native/'private_native_gene_axis.tsv',sep='\t',keep_default_na=False);assert genes.native_column.tolist()==list(range(31831)) and genes.gene.is_unique
    raw=np.load(native/'sum_raw_counts.npy',mmap_mode='r');assert raw.shape==(11701,31831) and raw.dtype==np.int64
    output=root/'NSCLC_native_MFA_inputs_20261016_v1';assert not output.exists();free=quota_free(root);output.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),qualification_sha256=sha(qual/'private_terminal.json'),native_audit_sha256=sha(native_proof),primary_views=previous['primary_views'],scope='All208 previously frozen author cases and four disjoint views; 145 discovery/63 internal holdout unchanged. No clinical outcomes or previous programme gene lists read.',normalization='Sum native integer counts over all author subtypes in each view. Whole31831-gene library sizes. Official edgeR4.0.16 TMM per case against a discovery-only pooled-count reference; factor ratios from public two-column API, geometric centring over discovery only. log2(1+1e6*count/effective_library_size).',technical_controls='View-own log1p cell count, cell-weighted mean log1p full native RNA and detected genes, symbol-prefix MT and RPL/RPS percent RNA; own fraction within original broad parent for CD8/Treg/Macrophage. B/Plasma equals entire B parent, so its structurally constant parent fraction is absent by definition. Train-only mean/SD and full-rank OLS; no automatic deletion or clinical/histology regressor.',features='At least100 raw counts in ceil(.25*145) discovery cases; exclude symbol MT-/RPL/RPS, positive discovery technical-residual variance. Top min(2000,eligible) residual variance per view, ties gene symbol, require>=500. Holdout never selects genes/centres/coefficients. One training RMS per view, no per-gene variance rescaling.',model_before_any_fit={'implementation':'Official mofapy2==0.7.2 private runtime','factors':10,'likelihood':'gaussian','spikeslab_weights':False,'ard_weights':True,'ard_factors':False,'spikeslab_factors':False,'dropR2':None,'max_iterations':1000,'convergence':'medium','primary_seed':20261017,'stability_seeds':[20261018,20261019],'selection':'Primary seed fixed; additional seeds only stability with discovery-only factor alignment, never pick best clinical/holdout fit. No k/penalty/endpoint search. All10 slots retained.'},holdout_before_any_fit={'score':'Each view separately: noise-weighted fixed-loading ridge argmin ||sqrt(tau)*(y-Wz)||²+||z||². Only own view expression, never posterior score using other views. Frozen model/input parameters. This is a defined ridge score, not claimed exact variational posterior.','family':60,'slots':'10 fixed factors x6 unordered view pairs, all63 identical cases; absent/constant/unestimable slots retained p=1 internally.','point':'Pearson of average ranks centred within fixed histology x therapy-group blocks; blocks derived only from official treatment/histology. Also unblocked Spearman descriptive. Singleton blocks contribute zero centred ranks. Need>=30 cases in blocks of size>=2, nonzero residual variance.','permutations':19999,'permutation':'Shuffle B-case assignments only within frozen blocks. Plus-one two-sided P, no precision escalation. Conditional case-exchangeability assumption; author cases are not newly reviewed independent patients.','bootstrap':2000,'interval':'Resample author cases within frozen blocks, recalculate ranks and block means each draw; retain failures, CI only>=1900valid. LOO descriptive.','multiplicity':'Complete60 BH; BY and joint.99 Clopper-Pearson Monte Carlo upper-bound BH separately. Seeds fixed by factor/pair names. No post-fit slot selection.','limits':'Internal source consistency, not independent patient validation or physical LR interaction. Subtype mixing/clinical treatment effects remain distinct from technical control; complete subtype-composition diagnostics follow without replacing main fit.'},quota_free_bytes=free,output_cap_bytes=268_435_456,new_expression_download=False,original_identity_or_CM_changed=False)
    (output/'plan_before_native_expression_and_MFA_fit.json').write_text(json.dumps(plan,indent=2)+'\n');axis.to_csv(output/'private_case_axis.tsv',sep='\t',index=False);(output/'frozen_TMM_normalize.R').write_text(R_NORMALIZE)
    rscript=Path(sys.executable).resolve().parents[2]/'scrna_r/bin/Rscript';runtime=Path(os.environ['PANCANCER_PROJECT_ROOT'])/'Project_v5/runs/20261014_DIALOGUE_method_pilot_v1';summary=[];quality_rows=[]
    for view in plan['primary_views']:
        labels=previous['view_author_labels'][view];parent='T/NK cell' if view in ['author_CD8','author_Treg'] else 'Myeloid cell' if view=='author_Macrophage' else 'B cell'
        counts=np.zeros((len(cases),len(genes)),np.int64);quality=[]
        for i,case in enumerate(cases):
            g=units[units.sampleID.eq(case)&units.sub_cell_type.isin(labels)];indices=g.unit_index.to_numpy(int);n=int(g.n_cells.sum());assert n>=20
            counts[i]=np.asarray(raw[indices]).sum(axis=0,dtype=np.int64);total=int(counts[i].sum());assert total==int(g.native_total_RNA.sum()) and total>0
            parent_n=int(units.loc[units.sampleID.eq(case)&units.major_cell_type.eq(parent),'n_cells'].sum());assert parent_n>=n
            mt=genes.gene.str.startswith('MT-').to_numpy();rb=genes.gene.str.match(r'^(RPL|RPS)').to_numpy()
            q=[np.log1p(n),float(np.average(g.mean_logRNA,weights=g.n_cells)),float(np.average(g.mean_loggenes,weights=g.n_cells)),100*float(counts[i,mt].sum())/total,100*float(counts[i,rb].sum())/total]
            names=['logN','mean_logRNA','mean_loggenes','symbol_pct_mito','symbol_pct_ribo']
            if view!='author_B_Plasma':q.append(n/parent_n);names.append('fraction_within_broad_parent')
            quality.append(q);quality_rows.append(dict(sampleID=case,view=view,n_cells=n,parent_cells=parent_n,full_native_RNA=total,**dict(zip(names,q))))
        file=output/('private_'+view+'_pseudobulk.tsv.gz');pd.DataFrame(counts.T,index=genes.gene,columns=cases).to_csv(file,sep='\t',compression='gzip',index_label='gene')
        tmm=output/('private_'+view+'_TMM.tsv')
        with (output/('private_'+view+'_R_normalization.log')).open('x') as log:subprocess.run([str(rscript),str(output/'frozen_TMM_normalize.R'),str(file),str(output/'private_case_axis.tsv'),str(runtime),str(tmm)],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=300)
        norm=pd.read_csv(tmm,sep='\t',dtype={'sampleID':str});assert norm.sampleID.tolist()==cases and np.isfinite(norm.effective_library_size).all() and (norm.effective_library_size>0).all()
        value=np.log2(1+counts.astype(float)*1e6/norm.effective_library_size.to_numpy()[:,None]);q=np.array(quality);mean=q[train].mean(axis=0);sd=q[train].std(axis=0);assert (sd>1e-10).all();design=np.column_stack([np.ones(len(q)),(q-mean)/sd]);assert np.linalg.matrix_rank(design[train])==design.shape[1] and np.linalg.cond(design[train])<1e8
        coefficients=np.linalg.lstsq(design[train],value[train],rcond=None)[0];residual=value-design@coefficients;var=residual[train].var(axis=0);detected=(counts[train]>=100).sum(axis=0);eligible=(detected>=int(np.ceil(.25*train.sum())))&~(mt|rb)&(var>1e-12);indices=np.flatnonzero(eligible);assert len(indices)>=500;selected=indices[np.lexsort((genes.gene.to_numpy()[indices],-var[indices]))[:2000]];rms=float(np.sqrt(np.mean(np.square(residual[train][:,selected]))));assert rms>1e-10;x=residual[:,selected]/rms;assert np.isfinite(x).all()
        np.savez_compressed(output/('private_'+view+'_prepared.npz'),X=x,selected_native_columns=selected,selected_genes=genes.gene.to_numpy(str)[selected],regression_coefficients=coefficients[:,selected],quality= q,quality_mean=mean,quality_sd=sd,design=design,training_view_RMS=np.array(rms),cases=np.array(cases),train=train,covariate_names=np.array(names))
        pd.DataFrame(dict(native_column=genes.native_column,discovery_ge100_cases=detected,technical_excluded=mt|rb,technical_residual_variance=var,eligible=eligible,selected=np.isin(np.arange(len(genes)),selected))).to_csv(output/('private_'+view+'_all_native_feature_selection.tsv.gz'),sep='\t',index=False)
        summary.append(dict(view=view,discovery_cases=int(train.sum()),internal_holdout_cases=int((~train).sum()),selected_features=len(selected),eligible_features=int(eligible.sum()),technical_covariates=len(names),design_rank=design.shape[1],training_view_RMS=rms,holdout_used_for_selection=False))
        assert sum(p.stat().st_size for p in output.iterdir() if p.is_file())<plan['output_cap_bytes']
    pd.DataFrame(summary).to_csv(output/'private_all_four_view_input_summary.tsv',sep='\t',index=False);pd.DataFrame(quality_rows).to_csv(output/'private_all832_case_view_quality.tsv',sep='\t',index=False)
    terminal=dict(status='COMPLETE_FOUR_NATIVE_PSEUDOBULK_VIEWS_TRAIN_ONLY_TMM_QC_REGRESSION_HVG_NO_MODEL_FIT_AUDIT_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),files={p.name:sha(p) for p in output.iterdir() if p.is_file()},expression_views_prepared=True,model_or_clinical_fit=False);(output/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(terminal))


if __name__=='__main__':main()
