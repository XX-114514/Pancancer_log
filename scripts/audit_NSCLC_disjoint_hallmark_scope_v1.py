"""Independent complete native-column and unit scope audit without RNA scores."""
import hashlib,json,os,time
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);o=r/'NSCLC_disjoint_hallmark_preflight_20261016_v1';job=r/'HNSCC_receiver_followup_job_control_20261016_v2/preflight_NSCLC_disjoint_hallmark_v1'
    receipt=json.loads((job/'launch_receipt.json').read_text());deadline=time.monotonic()+7200
    while not (o/'private_terminal.json').exists():
        p=Path('/proc')/str(receipt['pid'])/'stat';assert p.exists();s=p.read_text().rsplit(')',1)[1].split()
        assert s[19]==str(receipt['proc_start_ticks']) and s[0]!='Z' and time.monotonic()<deadline;time.sleep(5)
    terminal=json.loads((o/'private_terminal.json').read_text())
    for n,h in terminal['files'].items():assert sha(o/n)==h
    plan=json.loads((o/'plan_before_functional_RNA_measurements.json').read_text());native=r/'NSCLC_full_native_gene_units_20261016_v1'
    inputs=dict(matrix=native/'mean_logCP10K.npy',gene_axis=native/'private_native_gene_axis.tsv',scores=r/'NSCLC_context_native_programme_projection_20261016_v1/private_native_author_subtype_programme_scores.tsv.gz',
                parameters=r/'NSCLC_context_frozen_parameters_20261016_v1/private_all53_frozen_gene_parameters.tsv',clinical_scope=r/'NSCLC_frozen_signature_clinical_preflight_20261016_v1/private_all_planned_signature_subtype_designs.tsv',
                clinical=r/'NSCLC_frozen_signature_clinical_preflight_20261016_v1/private_all243_clinical_eligibility.tsv')
    for key,p in inputs.items():assert sha(p)==plan['input_sha256'][key]
    s=pd.read_csv(inputs['scores'],sep='\t');axis=pd.read_csv(o/'private_score_row_axis.tsv.gz',sep='\t')
    assert axis.score_row.tolist()==list(range(len(s))) and len(s)==plan['units']
    pd.testing.assert_frame_equal(axis[s[axis.columns[1:]].columns],s[axis.columns[1:]],check_dtype=False)
    assert not axis.duplicated(['context','type','program','unit_index']).any()
    g=pd.read_csv(inputs['gene_axis'],sep='\t');assert g.gene.is_unique;lookup=dict(zip(g.gene,g.native_column))
    assert np.load(inputs['matrix'],mmap_mode='r').shape==(11701,len(g))
    params=pd.read_csv(inputs['parameters'],sep='\t');sc=pd.read_csv(o/'private_all1200_functional_scope.tsv',sep='\t');masks=json.loads((o/'private_physical_column_masks.json').read_text())
    refs=r/'CRC_fine_pair_DIALOGUE_code_20261014_v1/hallmark50_public_reference_v3';ref_axis=pd.read_csv(o/'private_reference_axis.tsv',sep='\t');assert len(ref_axis)==50
    assert ref_axis.reference.tolist()==sorted(plan['reference_sha256']) and ref_axis.reference_column.tolist()==list(range(50))
    assert len(sc)==1200 and sc.signature.nunique()==24 and not sc.duplicated(['signature','reference']).any();verified=0
    for row in sc.itertuples():
        forbidden=set(params.loc[params.context.eq(row.context)&params.program.eq(row.program),'gene']);assert len(forbidden)>0
        expected_forbidden=sorted(lookup[x] for x in forbidden if x in lookup);mask=masks[row.mask_id];assert expected_forbidden==mask['forbidden_native_columns']
        path=refs/(row.reference+'.json');assert sha(path)==plan['reference_sha256'][row.reference]
        original=set(json.loads(path.read_text())[row.reference]['geneSymbols']);retained=original.difference(forbidden);observed=retained.intersection(lookup)
        physical=sorted(lookup[x] for x in observed);assert physical==mask['targets'][row.reference]
        assert not set(physical).intersection(expected_forbidden)
        coverage=len(observed)/len(retained) if retained else 0.;background=len(g)-len(expected_forbidden)
        assert row.original_reference_genes==len(original) and row.reference_genes_after_training_exclusion==len(retained) and row.native_measured_targets==len(observed)
        assert np.isclose(coverage,row.remaining_target_coverage) and np.isclose(len(retained)/len(original),row.original_reference_fraction_retained)
        assert background==row.native_background_genes and background-len(observed)==row.native_nontarget_background_genes
        eligible=coverage>=.95 and len(observed)>=20 and background-len(observed)>=500
        assert (row.status=='MEASURABLE_INPUT_NO_RNA_SCORE')==eligible
        assert row.reference_column==int(ref_axis.set_index('reference').loc[row.reference,'reference_column'])
        assert row.signature==f'{row.context}_{row.type}_{row.program}' and ((s.context==row.context)&s.type.eq(row.type)&s.program.eq(row.program)).any();verified+=1
    c=pd.read_csv(inputs['clinical_scope'],sep='\t');desc=pd.read_csv(o/'private_all10600_description_scope.tsv',sep='\t');assert len(c)==212 and len(desc)==10600
    assert not desc.duplicated(['context','type','program','sub_cell_type','reference']).any()
    assert set(map(tuple,desc[['context','type','program','sub_cell_type','original_clinical_record']].drop_duplicates().to_numpy()))==set(map(tuple,c[['context','type','program','sub_cell_type','record']].to_numpy()))
    assert desc.groupby('original_clinical_record').size().eq(50).all()
    result=dict(status='PASS_ALL1200_MASKS_NATIVE_COVERAGE_ALL_ORIGINAL_UNITS_AND10600_DESCRIPTION_KEYS_NO_RNA_SCORES',finished_UTC=datetime.now(timezone.utc).isoformat(),
                audit_source_sha256=sha(__file__),producer_terminal_sha256=sha(o/'private_terminal.json'),functional_slots_checked=verified,native_rows_checked=len(axis),description_slots_checked=len(desc),new_functional_scores=False)
    p=o/'private_independent_scope_audit.json';assert not p.exists();p.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
