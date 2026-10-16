"""Effect-blind physical coverage of all original functions and programme slots."""
import hashlib,json,os,time,subprocess
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from external_functional_coverage_inputs_v1 import read_inputs


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def available(genes,lookup,forbidden=frozenset()):
    mapping={g:next(iter(lookup[g])) for g in genes if len(lookup.get(g,set()))==1};hist=Counter(mapping.values())
    return {g:c for g,c in mapping.items() if hist[c]==1 and c not in forbidden}


def one(mask,reference,axis,params):
    forbidden_genes=set(params.loc[params.context.eq(mask.context)&params.program.eq(mask.program),'gene']);forbidden=set().union(*(axis['lookup'].get(g,set()) for g in forbidden_genes)) if forbidden_genes else set();background=axis['native']-forbidden;remaining=reference-forbidden_genes;targets=available(remaining,axis['lookup'],forbidden);physical=set(targets.values())
    assert physical<=background;fraction=len(physical)/len(remaining) if remaining else 0.;gate=fraction>=.95 and len(physical)>=20 and len(background)-len(physical)>=500
    return dict(original_reference_genes=len(reference),reference_genes_after_training_exclusion=len(remaining),original_reference_fraction_retained=len(remaining)/len(reference),native_unique_targets=len(physical),remaining_reference_native_fraction=fraction,native_background_genes=len(background),native_nontarget_background_genes=len(background)-len(physical),unavailable_or_collision_or_masked_targets=len(remaining)-len(physical),technical_measurement_gate=gate)


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);op=r/'complete_external_lineage_opportunities_20261016_v2/private_independent_opportunity_audit.json';job=r/'HNSCC_receiver_followup_job_control_20261016_v2/audit_complete_external_lineage_opportunities_v2';rc=json.loads((job/'launch_receipt.json').read_text());deadline=time.monotonic()+14400
    while not op.exists():
        p=Path('/proc')/str(rc['pid'])/'stat';assert p.exists();s=p.read_text().rsplit(')',1)[1].split();assert s[19]==str(rc['proc_start_ticks']) and s[0]!='Z' and time.monotonic()<deadline;time.sleep(5)
    assert json.loads(op.read_text())['status'].startswith('PASS_ALL61_SOURCE_ARM')
    data=read_inputs(r);out=r/'complete_external_functional_gene_coverage_20261016_v1';assert not out.exists()
    df=subprocess.check_output(['df','-P',str(r.resolve())],text=True,timeout=20);fs=df.splitlines()[-1].split()[0];quota=subprocess.check_output(['quota','-g','-v','-w','USER002'],text=True,timeout=20);qr=[x.split() for x in quota.splitlines() if x.split() and x.split()[0]==fs];assert len(qr)==1;free=(int(qr[0][3])-int(qr[0][1].rstrip('*')))*1024;assert free>500_000_000_000+1_073_741_824+268_435_456
    out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),input_helper_sha256=sha(Path(__file__).with_name('external_functional_coverage_inputs_v1.py')),opportunity_audit_sha256=sha(op),input_sha256=data['input_sha256'],reference_sha256=data['reference_sha256'],legacy_gene_axis_sha256=data['legacy_gene_axis_sha256'],
        scope='All61 source-arm entries x all16 fixed masks x all50 functions; every original56 programme slot, including3 unsupported slots, retained per source. No effect or lineage/clinical eligibility selection.',
        policy='Exclude BOTH partners ALL learned gene names/possible physical columns from remaining reference and complete unambiguous native background; include native genes absent from NSCLC. Unique exact symbols/canonical/stable-ID aliases only; multiple matches and same-column target collisions unavailable, never choose by weights or effects.',
        gates='Function>=95% remaining reference native,>=20 unique targets,>=500 nontarget background. Frozen programme>=95% original absolute weight with collisions removed, no reweight/refit. Original retained-reference fraction and gene counts separate.',
        limits='Technical header coverage only, not RNA expression, exact NSCLC endpoint mapping, design rank, independent patients or clinical validation. Legacy axes have no new alias inference; unresolved mapping not biological absence.',quota_free_bytes=free,output_cap_bytes=268_435_456,new_expression_scores_P_or_download=False)
    (out/'plan_before_external_gene_coverage.json').write_text(json.dumps(plan,indent=2)+'\n')
    params=data['parameters'];refs=data['references'];masks=data['masks'];original=data['original_scope'].drop_duplicates(['mask_id','reference']).set_index(['mask_id','reference']);calibration=[]
    for mask in masks.itertuples():
        for name,genes in sorted(refs.items()):
            got=one(mask,genes,data['NSCLC_calibration'],params);old=original.loc[(mask.mask_id,name)]
            assert got['native_unique_targets']==old.native_measured_targets and got['native_background_genes']==old.native_background_genes and got['native_nontarget_background_genes']==old.native_nontarget_background_genes and got['reference_genes_after_training_exclusion']==old.reference_genes_after_training_exclusion
            calibration.append(dict(mask_id=mask.mask_id,reference=name,**got))
    assert len(calibration)==800;pd.DataFrame(calibration).to_csv(out/'private_all800_original_NSCLC_coverage_calibration.tsv',sep='\t',index=False)
    functions=[];programmes=[];axes=[]
    for (arm,source),axis in data['source_axes'].items():
        axes.append(dict(arm=arm,dataset=source,physical_native_features=axis['physical_features'],unambiguous_background_features=len(axis['native']),mapping_policy=axis['mapping'],clinical_or_independence_accepted=False))
        for mask in masks.itertuples():
            for name,genes in sorted(refs.items()):functions.append(dict(arm=arm,dataset=source,mask_id=mask.mask_id,context=mask.context,program=mask.program,reference=name,**one(mask,genes,axis,params)))
        for slot in data['programmes'].itertuples():
            p=params[params.context.eq(slot.context)&params.type.eq(slot.type)&params.program.eq(slot.program)];row=dict(arm=arm,dataset=source,context=slot.context,type=slot.type,program=slot.program,original_endpoint=slot.original_endpoint,NSCLC_curated_family=slot.source_family,upstream_status=slot.original_status,frozen_genes=len(p),native_unique_genes=0,native_absolute_weight_fraction=np.nan,programme95_gate=False,status='UPSTREAM_NO_FROZEN_COEFFICIENTS')
            if len(p):
                mapped=available(set(p.gene),axis['lookup']);use=p.gene.isin(mapped);fraction=float(p.weight.abs().to_numpy()[use].sum()/p.weight.abs().sum());row.update(native_unique_genes=len(mapped),native_absolute_weight_fraction=fraction,programme95_gate=fraction>=.95,status='TECHNICALLY_COVERED_NO_RNA_OR_LINEAGE_ACCEPTANCE' if fraction>=.95 else 'INCOMPLETE_OR_UNRESOLVED_NATIVE_WEIGHT_NO_ZERO_IMPUTATION')
            programmes.append(row)
    assert len(functions)==61*800 and len(programmes)==61*56
    for name,rows in [('private_all48800_functional_native_coverage.tsv',functions),('private_all3416_programme_native_coverage.tsv',programmes),('private_all61_physical_source_axis_summaries.tsv',axes)]:pd.DataFrame(rows).to_csv(out/name,sep='\t',index=False,na_rep='NA')
    assert sum(p.stat().st_size for p in out.iterdir())<plan['output_cap_bytes']
    terminal=dict(status='COMPLETE_ALL61_SOURCE_ARM_ALL48800_FUNCTION_ALL3416_PROGRAMME_COVERAGE_AND800_CALIBRATION_AUDIT_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),source_arm_entries=len(axes),functional_rows=len(functions),programme_rows=len(programmes),files={p.name:sha(p) for p in out.iterdir() if p.is_file()},no_expression_scores_P_or_download=True)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(terminal))


if __name__=='__main__':main()
