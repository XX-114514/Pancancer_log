"""Freeze all measured signatures/reference masks and native units before functional RNA."""
import hashlib,json,os,subprocess
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);native=r/'NSCLC_full_native_gene_units_20261016_v1';projection=r/'NSCLC_context_native_programme_projection_20261016_v1'
    na=r/'NSCLC_full_native_gene_unit_audit_20261016_v1/private_terminal.json';pa=projection/'private_independent_audit.json'
    assert json.loads(na.read_text())['status'].startswith('PASS_ALL_OUTPUTS_OLD_OVERLAP') and json.loads(pa.read_text())['status'].startswith('PASS_ALL56')
    old=json.loads((native/'private_terminal.json').read_text());assert sha(native/'private_terminal.json')==json.loads(na.read_text())['producer_terminal_sha256']
    assert sha(projection/'private_terminal.json')==json.loads(pa.read_text())['source_terminal_sha256']
    paths=dict(matrix=native/'mean_logCP10K.npy',gene_axis=native/'private_native_gene_axis.tsv',scores=projection/'private_native_author_subtype_programme_scores.tsv.gz',
               projection_scope=projection/'private_all56_projection_scope.tsv',parameters=r/'NSCLC_context_frozen_parameters_20261016_v1/private_all53_frozen_gene_parameters.tsv',
               clinical=r/'NSCLC_frozen_signature_clinical_preflight_20261016_v1/private_all243_clinical_eligibility.tsv',clinical_scope=r/'NSCLC_frozen_signature_clinical_preflight_20261016_v1/private_all_planned_signature_subtype_designs.tsv',
               official_ULM=r/'fullgene_receiver_target_code_v1/collectri_reference_v1/official_decoupler192_ULM_source.txt',native_audit=na,projection_audit=pa,
               decision=r/'NSCLC_post_joint_clinical_decision_20261016_v1/PRIVATE_DECISION_AND_NEXT_SCOPE.json',decision_audit=r/'NSCLC_post_joint_clinical_decision_20261016_v1/private_decision_scope_audit.json')
    refroot=r/'CRC_fine_pair_DIALOGUE_code_20261014_v1/hallmark50_public_reference_v3';paths['reference_provenance']=refroot/'reference_provenance.json'
    digests={k:sha(p) for k,p in paths.items()}
    assert digests['matrix']==old['outputs']['mean_logCP10K.npy']['sha256'] and digests['gene_axis']==old['outputs']['private_native_gene_axis.tsv']['sha256']
    previous=json.loads((r/'NSCLC_two_end_coordination_preflight_20261016_v1/plan_before_any_two_end_effect.json').read_text())
    assert digests['scores']==previous['input_sha256']['scores'] and digests['parameters']==previous['input_sha256']['parameters']
    assert digests['official_ULM']=='17954dbf35dbb10f358a330d271588867fc54be11552997726992353157664a1'
    out=r/'NSCLC_disjoint_hallmark_preflight_20261016_v1';assert not out.exists();out.mkdir(mode=0o700)
    df=subprocess.run(['df','-P',str(r.resolve())],capture_output=True,text=True,timeout=20);quota=subprocess.run(['quota','-g','-v','-w','USER002'],capture_output=True,text=True,timeout=20)
    assert df.returncode==quota.returncode==0;fs=df.stdout.splitlines()[-1].split()[0];qr=[x.split() for x in quota.stdout.splitlines() if x.split() and x.split()[0]==fs];assert len(qr)==1
    remaining=(int(qr[0][3])-int(qr[0][1].rstrip('*')))*1024;assert remaining>500_000_000_000+1_073_741_824+268_435_456
    mem={x.split(':')[0]:x.split(':')[1].strip() for x in Path('/proc/meminfo').read_text().splitlines()};ram=int(mem['MemAvailable'].split()[0])*1024;assert ram>2*1024**3
    s=pd.read_csv(paths['scores'],sep='\t');genes=pd.Index(pd.read_csv(paths['gene_axis'],sep='\t').gene);params=pd.read_csv(paths['parameters'],sep='\t')
    slots=s[['context','type','program','source_family']].drop_duplicates().sort_values(['context','type','program']);assert len(slots)==24 and genes.is_unique
    pscope=pd.read_csv(paths['projection_scope'],sep='\t');expected=pscope[pscope.status.str.startswith('PARTIAL_FROZEN')]
    assert set(map(tuple,slots[['context','type','program']].to_numpy()))==set(map(tuple,expected[['context','type','program']].to_numpy()))
    assert np.load(paths['matrix'],mmap_mode='r').shape==(11701,len(genes))
    prov=json.loads(paths['reference_provenance'].read_text());assert len(prov['records'])==50;refs={};refhash={}
    for z in prov['records']:
        p=refroot/(z['name']+'.json');assert sha(p)==z['sha256'];refs[z['name']]=set(json.loads(p.read_text())[z['name']]['geneSymbols']);refhash[z['name']]=sha(p)
    references=sorted(refs);assert len(references)==50
    estimated_output=len(s)*50*(8*3+2)+16*1024**2;assert estimated_output<268_435_456
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),kernel_sha256=sha(Path(__file__).with_name('disjoint_hallmark_ULM_kernel_v1.py')),
              input_sha256=digests,reference_sha256=refhash,scope='All24 actually measured frozen signatures, all50 references and every original native unit including n<20; all212 original author-subtype clinical slots x50 descriptions retained.',
              mask='Exclude union of BOTH original partner learned genes from targets AND complete native background. Exact native physical symbols, no fuzzy mapping or refit.',
              gates='Remaining reference>=95% measured, >=20 targets, >=500 non-target background; reference-retained fraction reported, fragment not whole pathway.',
              measurement='Pinned officialULM1.9.2, target+1/rest0, per-cell log1pCP10K mean; gene-pattern t is NOT case P or protein activity. No new P/BH.',
              QA='Also store target/non-target means and actual positive-RNA target count; zero foreground can induce background-driven scores and is not functional response evidence.',
              inference='Initially no-new-P full functional descriptions under original complete clinical/quality controls; no old clinical/P/BH endpoint changes.',
              units=len(s),unique_native_units=int(s.unit_index.nunique()),signature_slots=24,functional_slots=1200,description_slots=10600,batch_units=128,
              raw_independent_replay_scope='Previously all output/old overlap +64 prespecified original units; this stage audits ALL new scores against the fixed cache, not all raw cells again.',
              quota_free_bytes=remaining,RAM_available_bytes=ram,quota_reserve_decimal_bytes=500_000_000_000,output_cap_bytes=268_435_456,estimated_output_bytes=estimated_output)
    (out/'plan_before_functional_RNA_measurements.json').write_text(json.dumps(plan,indent=2)+'\n')
    masks={};coverage=[]
    for row in slots.itertuples():
        signature=f'{row.context}_{row.type}_{row.program}';mask_id=f'C{row.context}_{row.program}'
        p=params[params.context.eq(row.context)&params.program.eq(row.program)];assert len(p)>0;forbidden=set(p.gene)
        forbidden_native=np.flatnonzero(genes.isin(forbidden));bg=np.flatnonzero(~genes.isin(forbidden))
        if mask_id not in masks:masks[mask_id]=dict(context=int(row.context),program=row.program,forbidden_native_columns=forbidden_native.tolist(),targets={})
        for j,name in enumerate(references):
            remaining_set=refs[name]-forbidden;ix=genes.get_indexer(sorted(remaining_set));target=sorted(int(x) for x in ix if x>=0);fraction=len(target)/len(remaining_set) if remaining_set else 0.
            status='MEASURABLE_INPUT_NO_RNA_SCORE' if fraction>=.95 and len(target)>=20 and len(bg)-len(target)>=500 else 'INSUFFICIENT_NATIVE_REMAINING_TARGET_OR_BACKGROUND_COVERAGE'
            masks[mask_id]['targets'][name]=target
            coverage.append(dict(readout_id=f'HMRNA{len(coverage)+1:04d}',signature=signature,context=row.context,type=row.type,program=row.program,source_family=row.source_family,mask_id=mask_id,reference=name,reference_column=j,
                                 original_reference_genes=len(refs[name]),reference_genes_after_training_exclusion=len(remaining_set),native_measured_targets=len(target),remaining_target_coverage=fraction,
                                 original_reference_fraction_retained=len(remaining_set)/len(refs[name]),native_background_genes=len(bg),native_nontarget_background_genes=len(bg)-len(target),status=status))
    scope=pd.DataFrame(coverage);assert len(scope)==1200 and not scope.duplicated(['signature','reference']).any()
    scope.to_csv(out/'private_all1200_functional_scope.tsv',sep='\t',index=False)
    (out/'private_physical_column_masks.json').write_text(json.dumps(masks,separators=(',',':'))+'\n')
    axis=s[['sampleID','sub_cell_type','unit_index','n_cells','context','type','program','source_family']].copy();axis.insert(0,'score_row',np.arange(len(s)))
    assert not axis.duplicated(['context','type','program','unit_index']).any();axis.to_csv(out/'private_score_row_axis.tsv.gz',sep='\t',index=False,compression='gzip')
    pd.DataFrame(dict(reference_column=np.arange(50),reference=references)).to_csv(out/'private_reference_axis.tsv',sep='\t',index=False)
    clin=pd.read_csv(paths['clinical_scope'],sep='\t');assert len(clin)==212
    desc=[]
    for z in clin.itertuples():
        for name in references:desc.append(dict(description_id=f'HMDF{len(desc)+1:05d}',original_clinical_record=z.record,context=z.context,type=z.type,program=z.program,sub_cell_type=z.sub_cell_type,reference=name,old_eligible_cases=z.eligible_cases))
    pd.DataFrame(desc).to_csv(out/'private_all10600_description_scope.tsv',sep='\t',index=False)
    terminal=dict(status='COMPLETE_ALL24_SIGNATURE_ALL1200_FUNCTION_ALL10600_DESCRIPTION_SCOPE_NO_RNA_SCORES',finished_UTC=datetime.now(timezone.utc).isoformat(),
                  signature_slots=24,functional_scope=1200,native_score_rows=len(s),coverage_status_counts=scope.status.value_counts().to_dict(),files={p.name:sha(p) for p in out.iterdir() if p.is_file()},new_RNA_scores_computed=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(terminal))


if __name__=='__main__':main()
