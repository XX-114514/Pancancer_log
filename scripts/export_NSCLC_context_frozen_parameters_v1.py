"""Export supported frozen coefficients and source-blind partial calibration.

No refit, new cell scores, clinical outcomes or identity assignments.
"""
import hashlib,json,os,subprocess,sys
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import spearmanr


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


RCODE=r'''
args<-commandArgs(trailingOnly=TRUE);projection<-args[1];components<-args[2];out<-args[3]
s<-read.delim(file.path(projection,'public_slot_scope.tsv'),check.names=FALSE)
s<-s[s$status=='FROZEN_RECONSTRUCTION_PASS',];rows<-list();checks<-list()
for(i in seq_len(nrow(s))){
 q<-s[i,];key<-paste(q$context,q$type,q$program,sep='_')
 w<-readRDS(file.path(projection,paste0(key,'_weights.rds')))
 b<-readRDS(file.path(components,paste0(key,'_gene_coefficients.rds')))
 stopifnot(identical(as.character(w$genes),as.character(b$genes)),max(abs(w$weight-b$weights))<1e-12)
 stopifnot(identical(b$quality_columns,c('(Intercept)','cellQ','geneQ','fraction','logN')))
 stopifnot(nrow(b$gene_quality_coefficients)==5,ncol(b$gene_quality_coefficients)==length(w$genes))
 reconstructed<-as.numeric(b$gene_quality_coefficients%*%w$weight)
 err<-max(abs(reconstructed-as.numeric(w$quality_coefficients)))
 stopifnot(is.finite(err),err<1e-6,all(w$gene_sd>0),!anyDuplicated(w$genes))
 a<-data.frame(context=q$context,type=q$type,program=q$program,gene=w$genes,weight=w$weight,gene_mean=w$gene_mean,gene_sd=w$gene_sd,
   correction_intercept=b$gene_quality_coefficients[1,],correction_cellQ=b$gene_quality_coefficients[2,],
   correction_geneQ=b$gene_quality_coefficients[3,],correction_fraction=b$gene_quality_coefficients[4,],correction_logN=b$gene_quality_coefficients[5,],check.names=FALSE)
 rows[[key]]<-a;checks[[key]]<-data.frame(context=q$context,type=q$type,program=q$program,genes=length(w$genes),maximum_weighted_quality_coefficient_error=err)
}
write.table(do.call(rbind,rows),file.path(out,'private_all53_frozen_gene_parameters.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
write.table(do.call(rbind,checks),file.path(out,'private_all53_parameter_parity.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
'''


def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);V=Path(json.loads((R/'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root'])
    F=V/'state_after_DIALOGUE_frozen_projection_20261016_v1';C=V/'state_after_DIALOGUE_gene_components_20261016_v1'
    O=R/'NSCLC_context_frozen_parameters_20261016_v1';assert not O.exists();O.mkdir(mode=0o700)
    fi=json.loads((F/'inventory.json').read_text());ci=json.loads((C/'inventory.json').read_text())
    scope=pd.read_csv(F/'public_slot_scope.tsv',sep='\t');supported=scope[scope.status.eq('FROZEN_RECONSTRUCTION_PASS')];assert len(scope)==56 and len(supported)==53
    for row in supported.itertuples():
        key=f'{row.context}_{row.type}_{row.program}'
        for root,name,inv in [(F,key+'_weights.rds',fi),(C,key+'_gene_coefficients.rds',ci)]:assert sha(root/name)==inv['files'][name]
    script=O/'private_export_frozen_RDS.R';script.write_text(RCODE)
    rb=Path(os.environ.get('PANCANCER_DIALOGUE_RSCRIPT',str(Path(sys.executable).parents[2]/'scrna_r/bin/Rscript')));assert rb.is_file()
    (O/'plan_before_coefficient_export.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),R_code_sha256=sha(script),
        model_scope='All53 supported slots, all3 unsupported slots retained separately; no clinical or LR-effect selection.',
        expression='Freeze original gene mean/SD, signed weight, every gene-level five-term correction, and training quality mean/SD.',
        gate='Keep original95% native absolute-weight rule; no scaling remaining weights or filling absent genes.',
        calibration='For all original blocks, compare full versus NSCLC-observed-gene contribution sums separately by original discovery/holdout partition; descriptive, no P or outcome tests.',
        limitations='Partial-score calibration is within discovery source, not NSCLC biological replication, source identity matching or training-domain acceptance.'),indent=2)+'\n')
    subprocess.run([str(rb),'--vanilla',str(script),str(F),str(C),str(O)],check=True,timeout=300)
    params=pd.read_csv(O/'private_all53_frozen_gene_parameters.tsv',sep='\t');assert not params.duplicated(['context','type','program','gene']).any()
    parity=pd.read_csv(O/'private_all53_parameter_parity.tsv',sep='\t');assert len(parity)==53 and parity.maximum_weighted_quality_coefficient_error.lt(1e-6).all()
    native=set(pd.read_csv(R/'NSCLC_full_native_gene_units_20261016_v1/private_native_gene_axis.tsv',sep='\t').gene)
    B=pd.read_csv(C/'private_all_block_gene_components.tsv.gz',sep='\t');quality=[];calibration=[]
    for row in scope.itertuples():
        root=V/'state_after_DIALOGUE_inputs_20261016_v3'/f'context_{row.context:04d}'
        inv=json.loads((root/'private_input_inventory.json').read_text());name=inv['endpoints'][0 if row.type=='typeA' else 1]
        if (row.context,row.type) not in [(v['context'],v['type']) for v in quality]:
            q=np.load(root/(row.type+'_frozen_transform.npz'),allow_pickle=False)
            quality.append(dict(context=row.context,type=row.type,original_endpoint=name,mean_logRNA=float(q['quality_mean'][0]),mean_loggenes=float(q['quality_mean'][1]),sd_logRNA=float(q['quality_sd'][0]),sd_loggenes=float(q['quality_sd'][1])))
        if row.status!='FROZEN_RECONSTRUCTION_PASS':
            calibration.append(dict(context=row.context,type=row.type,program=row.program,status=row.status));continue
        p=params[params.context.eq(row.context)&params.type.eq(row.type)&params.program.eq(row.program)].set_index('gene')
        f=B[B.context.eq(row.context)&B.type.eq(row.type)&B.program.eq(row.program)].copy()
        assert set(f.gene)==set(p.index)
        observed=p.index.isin(native);fraction=float(np.abs(p.weight[observed]).sum()/np.abs(p.weight).sum())
        assert np.allclose(f.weight,p.weight.reindex(f.gene).to_numpy(),atol=1e-12,rtol=1e-12)
        f['full_contribution']=f.contribution*f.weight;f['native_observed_contribution']=f.full_contribution.where(f.gene.isin(native),0.)
        for partition,d in f.groupby('partition',sort=True):
            means=d.groupby('sample')[['full_contribution','native_observed_contribution']].sum()
            x=means.full_contribution.to_numpy();y=means.native_observed_contribution.to_numpy()
            calibration.append(dict(context=row.context,type=row.type,program=row.program,status='SUPPORTED_FIXED_PARTIAL_CALIBRATION',partition=partition,
                blocks=len(means),native_absolute_weight=fraction,native95_gate=bool(fraction>=.95),
                full_partial_spearman=float(spearmanr(x,y).statistic),RMSE=float(np.sqrt(np.mean((x-y)**2))),
                full_IQR=float(np.quantile(x,.75)-np.quantile(x,.25)),partial_IQR=float(np.quantile(y,.75)-np.quantile(y,.25))))
    pd.DataFrame(quality).to_csv(O/'private_all28_training_quality_transforms.tsv',sep='\t',index=False)
    pd.DataFrame(calibration).to_csv(O/'private_all56_slot_partial_score_calibration.tsv',sep='\t',index=False,na_rep='NA')
    report=dict(status='ALL53_FROZEN_GENE_PARAMETERS_28_QUALITY_TRANSFORMS_AND56_SLOT_PARTIAL_CALIBRATION_COMPLETE_NO_NSCLC_SCORES',
        finished_UTC=datetime.now(timezone.utc).isoformat(),supported_slots=len(parity),quality_transforms=len(quality),
        source_projection_inventory_sha256=sha(F/'inventory.json'),source_component_inventory_sha256=sha(C/'inventory.json'),
        files={p.name:sha(p) for p in O.glob('*.tsv')},clinical_outcomes_used=False,fine_identity_mapping_validated=False,coefficients_refitted=False)
    (O/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()
