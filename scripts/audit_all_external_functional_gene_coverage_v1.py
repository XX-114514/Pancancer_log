"""Independently replay all physical coverage metadata, without effect selection."""
import hashlib,json,math,os,time
from collections import defaultdict
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from external_functional_coverage_inputs_v1 import read_inputs


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def resolve(genes,axis,blocked=frozenset()):
    # Invert candidates to physical-column owners; remove every collision owner.
    owners=defaultdict(list)
    for gene in sorted(genes):
        hits=axis['lookup'].get(gene,set())
        if len(hits)==1:
            column,=hits
            if column not in blocked:owners[column].append(gene)
    return {names[0]:column for column,names in owners.items() if len(names)==1}


def functional(mask,reference,axis,parameters):
    genes=set(parameters.loc[(parameters.context==mask.context)&(parameters.program==mask.program),'gene'])
    blocked=set()
    for gene in genes:blocked.update(axis['lookup'].get(gene,set()))
    retained=reference.difference(genes);targets=resolve(retained,axis,blocked)
    background=axis['native'].difference(blocked);assert set(targets.values()).issubset(background)
    n=len(targets);fraction=n/len(retained) if retained else 0.
    return {'original_reference_genes':len(reference),'reference_genes_after_training_exclusion':len(retained),
        'original_reference_fraction_retained':len(retained)/len(reference),'native_unique_targets':n,
        'remaining_reference_native_fraction':fraction,'native_background_genes':len(background),
        'native_nontarget_background_genes':len(background)-n,'unavailable_or_collision_or_masked_targets':len(retained)-n,
        'technical_measurement_gate':fraction>=.95 and n>=20 and len(background)-n>=500}


def compare(row,wanted):
    for key,value in wanted.items():
        got=getattr(row,key)
        if isinstance(value,float):assert np.isclose(got,value,rtol=0,atol=2e-14),(key,got,value)
        else:assert got==value,(key,got,value)


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=r/'complete_external_functional_gene_coverage_20261016_v1'
    rc=json.loads((r/'HNSCC_receiver_followup_job_control_20261016_v2/preflight_all_external_functional_gene_coverage_v1/launch_receipt.json').read_text());deadline=time.monotonic()+14400
    while not (out/'private_terminal.json').exists():
        p=Path('/proc')/str(rc['pid'])/'stat';assert p.exists();fields=p.read_text().rsplit(')',1)[1].split()
        assert fields[0]!='Z' and fields[19]==str(rc['proc_start_ticks']) and time.monotonic()<deadline;time.sleep(5)
    terminal=json.loads((out/'private_terminal.json').read_text());plan=json.loads((out/'plan_before_external_gene_coverage.json').read_text())
    assert terminal['code_sha256']==rc['source_sha256']==plan['code_sha256']
    for name,digest in terminal['files'].items():assert sha(out/name)==digest
    data=read_inputs(r);assert data['input_sha256']==plan['input_sha256'] and data['reference_sha256']==plan['reference_sha256'] and data['legacy_gene_axis_sha256']==plan['legacy_gene_axis_sha256']
    assert sha(Path(__file__).with_name('external_functional_coverage_inputs_v1.py'))==plan['input_helper_sha256']
    masks={row.mask_id:row for row in data['masks'].itertuples()};parameters=data['parameters'];refs=data['references']
    calibration=pd.read_csv(out/'private_all800_original_NSCLC_coverage_calibration.tsv',sep='\t');assert len(calibration)==800 and not calibration.duplicated(['mask_id','reference']).any()
    original=data['original_scope'].drop_duplicates(['mask_id','reference']).set_index(['mask_id','reference'])
    for row in calibration.itertuples():
        want=functional(masks[row.mask_id],refs[row.reference],data['NSCLC_calibration'],parameters);compare(row,want);old=original.loc[(row.mask_id,row.reference)]
        assert want['native_unique_targets']==old.native_measured_targets and want['native_background_genes']==old.native_background_genes and want['native_nontarget_background_genes']==old.native_nontarget_background_genes and want['reference_genes_after_training_exclusion']==old.reference_genes_after_training_exclusion
    functions=pd.read_csv(out/'private_all48800_functional_native_coverage.tsv',sep='\t');keys=['arm','dataset','mask_id','reference'];assert len(functions)==48800 and not functions.duplicated(keys).any()
    expected={(arm,source,mask,name) for arm,source in data['source_axes'] for mask in masks for name in refs};assert set(map(tuple,functions[keys].to_numpy()))==expected
    failures=0
    for row in functions.itertuples():
        mask=masks[row.mask_id];assert row.context==mask.context and row.program==mask.program
        want=functional(mask,refs[row.reference],data['source_axes'][(row.arm,row.dataset)],parameters);compare(row,want);failures+=not want['technical_measurement_gate']
    programmes=pd.read_csv(out/'private_all3416_programme_native_coverage.tsv',sep='\t',keep_default_na=False);pkeys=['arm','dataset','context','type','program'];assert len(programmes)==3416 and not programmes.duplicated(pkeys).any()
    slots={(s.context,s.type,s.program):s for s in data['programmes'].itertuples()};assert set(map(tuple,programmes[pkeys].to_numpy()))=={(arm,source,*slot) for arm,source in data['source_axes'] for slot in slots}
    unsupported=0
    for row in programmes.itertuples():
        slot=slots[(row.context,row.type,row.program)];assert row.original_endpoint==slot.original_endpoint and row.NSCLC_curated_family==slot.source_family and row.upstream_status==slot.original_status
        sub=parameters[(parameters.context==row.context)&(parameters.type==row.type)&(parameters.program==row.program)];assert row.frozen_genes==len(sub)
        if not len(sub):
            assert row.native_unique_genes==0 and row.native_absolute_weight_fraction=='NA' and not row.programme95_gate and row.status=='UPSTREAM_NO_FROZEN_COEFFICIENTS';unsupported+=1;continue
        mapped=resolve(set(sub.gene),data['source_axes'][(row.arm,row.dataset)]);numerator=math.fsum(abs(float(p.weight)) for p in sub.itertuples() if p.gene in mapped);denominator=math.fsum(abs(float(p.weight)) for p in sub.itertuples());fraction=numerator/denominator
        assert row.native_unique_genes==len(mapped) and abs(float(row.native_absolute_weight_fraction)-fraction)<2e-14 and row.programme95_gate==(fraction>=.95)
        assert row.status==('TECHNICALLY_COVERED_NO_RNA_OR_LINEAGE_ACCEPTANCE' if fraction>=.95 else 'INCOMPLETE_OR_UNRESOLVED_NATIVE_WEIGHT_NO_ZERO_IMPUTATION')
    axes=pd.read_csv(out/'private_all61_physical_source_axis_summaries.tsv',sep='\t');assert len(axes)==61 and not axes.duplicated(['arm','dataset']).any()
    for row in axes.itertuples():
        axis=data['source_axes'][(row.arm,row.dataset)];assert row.physical_native_features==axis['physical_features'] and row.unambiguous_background_features==len(axis['native']) and row.mapping_policy==axis['mapping'] and not row.clinical_or_independence_accepted
    # Explicit multihit, learned alias and collision fixtures, using no real genes.
    fixture={'lookup':{'a':{1},'alias':{1},'b':{2},'ambiguous':{3,4}},'native':{1,2,3,4}}
    assert resolve({'a','alias','b','ambiguous'},fixture)=={'b':2} and resolve({'a','b'},fixture,{1})=={'b':2}
    assert unsupported==61*3
    proof=dict(status='PASS_ALL48800_FUNCTION_ALL3416_PROGRAMME_ALL800_NSCLC_CALIBRATION_AND61_NATIVE_AXES',finished_UTC=datetime.now(timezone.utc).isoformat(),audit_source_sha256=sha(__file__),producer_terminal_sha256=sha(out/'private_terminal.json'),input_helper_sha256=plan['input_helper_sha256'],functional_rows_checked=len(functions),programme_rows_checked=len(programmes),calibration_rows_checked=len(calibration),unsupported_slots_retained=unsupported,collision_multihit_learned_alias_cases_checked=True,technical_failure_rows_present=bool(failures),new_expression_scores_P_or_clinical_acceptance=False)
    p=out/'private_independent_gene_coverage_audit.json';assert not p.exists();p.write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof))


if __name__=='__main__':main()
