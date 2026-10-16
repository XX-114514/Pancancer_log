"""Frozen programme projection within explicitly compatible author lineages.

All56 parent slots accounted. No clinical outcome read or fine identity transfer.
Apply each author subtype separately, preserving subtype/quality/composition.
"""
import hashlib,json,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);V=Path(json.loads((R/'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root'])
    N=R/'NSCLC_full_native_gene_units_20261016_v1';E=R/'NSCLC_context_frozen_parameters_20261016_v1'
    O=R/'NSCLC_context_native_programme_projection_20261016_v1';assert not O.exists();O.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),
        model_scope='All56 original context/endpoint/programme slots, preserving3 unsupported slots and all native-weight failures.',
        lineage_scope='Curated ancestry-compatible Macrophage (all author Mphi_ subtypes), CD8 (all author CD8T_ subtypes), Treg (all author CD4T_Treg_ subtypes). Other lineages remain explicitly unmapped, never projected onto immune cells.',
        identity='Ancestry-compatible expression measurement, not transfer/acceptance of any original fine subcluster identity or formal state.',
        units='Every original author case/subtype unit; no old selected200 or clinical outcome restrictions. n<20 scores descriptive and ineligible for later inference.',
        denominator='Macrophage fraction within all author macrophages; T subtype fraction within all author CD4T_, CD8T_ and T_gdT_, excluding NK/ILC. Preserve original major (Macrophage or T cell).',
        projection='Per gene (native mean per-cell log1pCP10K minus frozen gene mean)/frozen SD, subtract original gene-level intercept/cellQ/geneQ/fraction/logN correction, multiply unchanged signed frozen weight; sum only genuinely measured genes, no rescaling or zero imputation.',
        quality='True full-native mean(log1pRNA/genes), standardized by original discovery mean/SD; original subtype-cell-count and within-major correction retained.',
        coverage='Original >=95% native absolute weight rule. Report complete measurement separately. Never claim full programme when partial.',
        support='Report original discovery unit-level ranges of all4 covariates and native unit out-of-range flags; diagnostic, not effect-selected exclusions.',
        prerequisites='Independent all-output/old-overlap and64-unit/all-gene raw audit PASS, hash-verified frozen parameter export and source clinical scope reviewed.',
        clinical='No clinical outcomes merged or association tested. Freeze later case/exposure/outcome/family contracts before effects.')
    (O/'plan_before_native_programme_projection.json').write_text(json.dumps(plan,indent=2)+'\n')
    A=R/'NSCLC_full_native_gene_unit_audit_20261016_v1';j=R/'HNSCC_receiver_followup_job_control_20261016_v2/audit_NSCLC_full_native_gene_units_v1/launch_receipt.json'
    launch=json.loads(j.read_text());deadline=time.monotonic()+10800
    while not (A/'private_terminal.json').exists():
        p=Path('/proc')/str(launch['pid'])/'stat';assert p.exists()
        st=p.read_text().rsplit(')',1)[1].split();assert st[19]==str(launch['proc_start_ticks']) and st[0]!='Z'
        assert time.monotonic()<deadline;time.sleep(15)
    audit=json.loads((A/'private_terminal.json').read_text());assert audit['status']=='PASS_ALL_OUTPUTS_OLD_OVERLAP_AND_PRESPECIFIED64_RAW_UNITS_ALL_NATIVE_GENES'
    et=json.loads((E/'private_terminal.json').read_text())
    for name,h in et['files'].items():assert sha(E/name)==h
    assert json.loads((R/'NSCLC_all_case_clinical_scope_20261016_v1/private_terminal.json').read_text())['status'].startswith('ALL243_CASE_METADATA')
    params=pd.read_csv(E/'private_all53_frozen_gene_parameters.tsv',sep='\t');qt=pd.read_csv(E/'private_all28_training_quality_transforms.tsv',sep='\t').set_index(['context','type'])
    scope=pd.read_csv(V/'state_after_DIALOGUE_frozen_projection_20261016_v1/public_slot_scope.tsv',sep='\t')
    units=pd.read_csv(N/'private_all_native_units_and_quality.tsv',sep='\t',keep_default_na=False)
    old_units=pd.read_csv(R/'all_CM_state_native_input_preflight_20261014_v1/private_all_frozen_sample_node_units.tsv.gz',sep='\t',keep_default_na=False).set_index('unit_index',drop=False)
    major_totals=old_units.groupby(['sample_uid','major']).n_cells.sum()
    genes=pd.read_csv(N/'private_native_gene_axis.tsv',sep='\t').gene
    M=np.load(N/'mean_logCP10K.npy',mmap_mode='r');assert len(genes)==M.shape[1] and genes.is_unique
    fam={'Macrophage':units.sub_cell_type.str.startswith('Mφ_'), 'CD8':units.sub_cell_type.str.startswith('CD8T_'), 'Treg':units.sub_cell_type.str.startswith('CD4T_Treg_')}
    assert [units.loc[fam[k],'sub_cell_type'].nunique() for k in ['Macrophage','CD8','Treg']]==[11,10,3]
    isT=units.sub_cell_type.str.startswith(('CD4T_','CD8T_','T_gdT_'))
    denominators={'Macrophage':units[fam['Macrophage']].groupby('sampleID').n_cells.sum(),'T cell':units[isT].groupby('sampleID').n_cells.sum()}
    output=[];slots=[];domains=[]
    correction_columns=['correction_intercept','correction_cellQ','correction_geneQ','correction_fraction','correction_logN']
    for s in scope.itertuples():
        tr=qt.loc[(s.context,s.type)];endpoint=tr.original_endpoint
        family='Macrophage' if endpoint.startswith('Macrophage_') else 'CD8' if endpoint.startswith('CD8_T_cell_') else 'Treg' if endpoint.startswith('Regulatory_T_cell_') else None
        row=dict(context=s.context,type=s.type,program=s.program,original_endpoint=endpoint,source_family=family,status='PENDING',original_status=s.status,
            native_weight_fraction=np.nan,complete_native=False,author_subtypes=0,all_native_units=0,native_units_ge20=0,fine_identity_match_validated=False)
        if s.status!='FROZEN_RECONSTRUCTION_PASS':row['status']=s.status;slots.append(row);continue
        if family is None:row['status']='NO_CURATED_COMPATIBLE_NATIVE_LINEAGE';slots.append(row);continue
        p=params[params.context.eq(s.context)&params.type.eq(s.type)&params.program.eq(s.program)].copy()
        ix=pd.Index(genes).get_indexer(p.gene);present=ix>=0;fraction=float(np.abs(p.weight.to_numpy()[present]).sum()/np.abs(p.weight).sum())
        row.update(native_weight_fraction=fraction,complete_native=bool(present.all()))
        if fraction<.95:row['status']='NATIVE_WEIGHT_BELOW_FROZEN95_RULE';slots.append(row);continue
        f=units[fam[family]].copy();f['original_endpoint']=endpoint;f['source_family']=family;f['context']=s.context;f['type']=s.type;f['program']=s.program
        uindex=f.unit_index.to_numpy(int);major='Macrophage' if family=='Macrophage' else 'T cell'
        f['fraction']=f.n_cells.to_numpy()/denominators[major].reindex(f.sampleID).to_numpy()
        f['logN']=np.log1p(f.n_cells);f['cellQ']=(f.mean_logRNA-tr.mean_logRNA)/tr.sd_logRNA;f['geneQ']=(f.mean_loggenes-tr.mean_loggenes)/tr.sd_loggenes
        D=np.column_stack([np.ones(len(f)),f[['cellQ','geneQ','fraction','logN']].to_numpy()])
        active=p.loc[present];columns=ix[present];w=active.weight.to_numpy();mu=active.gene_mean.to_numpy();sd=active.gene_sd.to_numpy()
        x=(np.asarray(M[np.ix_(uindex,columns)])-mu)/sd;co=active[correction_columns].to_numpy().T
        raw=x@w;correction=D@(co@w);score=raw-correction
        # Independent algebraic arrangement, no reconstruction through producer code.
        other=np.sum((x-D@co)*w,axis=1);assert np.allclose(score,other,atol=1e-7,rtol=1e-8)
        f['score']=score;f['raw_weighted_score']=raw;f['frozen_quality_composition_correction']=correction
        f['native_weight_fraction']=fraction;f['complete_native_measurement']=bool(present.all());f['later_inference_cell_count_eligible']=f.n_cells.ge(20)
        root=V/'state_after_DIALOGUE_inputs_20261016_v3'/f'context_{s.context:04d}'
        md=pd.read_csv(root/(s.type+'_discovery_metadata.tsv'),sep='\t')
        train=md.groupby('unit_index')[['cellQ','geneQ']].mean();base=old_units.loc[train.index]
        assert base.subCluster.eq(endpoint).all() and base.major.eq(major).all()
        train['fraction']=base.n_cells.to_numpy()/np.array([major_totals.loc[(a,b)] for a,b in zip(base.sample_uid,base.major)])
        train['logN']=np.log1p(base.n_cells.to_numpy())
        for name in ['cellQ','geneQ','fraction','logN']:
            lo=float(train[name].min());hi=float(train[name].max());outside=(f[name]<lo)|(f[name]>hi)
            f[name+'_outside_training_unit_range']=outside
            domains.append(dict(context=s.context,type=s.type,program=s.program,original_endpoint=endpoint,source_family=family,covariate=name,
                discovery_unit_minimum=lo,discovery_unit_maximum=hi,native_unit_minimum=float(f[name].min()),native_unit_maximum=float(f[name].max()),
                native_units=len(f),native_units_outside=int(outside.sum()),native_ge20_units=int(f.n_cells.ge(20).sum()),native_ge20_outside=int((outside&f.n_cells.ge(20)).sum())))
        f['any_covariate_outside_training_unit_range']=f[[n+'_outside_training_unit_range' for n in ['cellQ','geneQ','fraction','logN']]].any(axis=1)
        output.append(f);row.update(status='PARTIAL_FROZEN_NATIVE_MEASUREMENT_NOT_FINE_IDENTITY_ACCEPTANCE',author_subtypes=int(f.sub_cell_type.nunique()),all_native_units=len(f),native_units_ge20=int(f.n_cells.ge(20).sum()))
        slots.append(row)
    scores=pd.concat(output,ignore_index=True);slots=pd.DataFrame(slots);assert len(slots)==56
    assert not scores.duplicated(['context','type','program','unit_index']).any() and np.isfinite(scores.score).all()
    scores.to_csv(O/'private_native_author_subtype_programme_scores.tsv.gz',sep='\t',index=False)
    slots.to_csv(O/'private_all56_projection_scope.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(domains).to_csv(O/'private_all_projected_slots_training_domain.tsv',sep='\t',index=False)
    report=dict(status='ALL56_CONTEXT_SLOTS_ACCOUNTED_FROZEN_AUTHOR_LINEAGE_SCORES_AND_TRAINING_DOMAIN_COMPLETE_NO_CLINICAL_TESTS',finished_UTC=datetime.now(timezone.utc).isoformat(),
        scope_status_counts=slots.status.value_counts().to_dict(),projected_slot_count=len(output),native_unit_score_records=len(scores),
        full_native_audit_sha256=sha(A/'private_terminal.json'),frozen_export_terminal_sha256=sha(E/'private_terminal.json'),
        files={p.name:sha(p) for p in O.iterdir() if p.name.endswith(('.tsv','.tsv.gz'))},gene_weights_or_corrections_refitted=False,
        clinical_outcomes_used=False,fine_identity_transfer_accepted=False,out_of_domain_exclusions_applied=False)
    (O/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()
