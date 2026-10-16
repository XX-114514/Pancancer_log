"""Audit all author cases and metadata conflicts, without outcome associations."""
import hashlib,json,os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd

FIELDS=['gender','age','smoking_history','cancer_type','pre_treatment_staging','anti-PD1_therapy','chemotherapy','targeted_therapy','cycles','pathological_response','pathological_response_rate','radiological_response']


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(16777216),b''):h.update(b)
    return h.hexdigest()


def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);O=R/'NSCLC_all_case_clinical_scope_20261016_v1'
    assert not O.exists();O.mkdir(mode=0o700)
    path=R/'GSE243013_NSCLC_immune_scRNA_metadata.csv.gz'
    (O/'plan_before_metadata_review.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),
        metadata_sha256=sha(path),scope='All source cases and author51 subtypes; no old selected200 restriction, no expression or effect selection.',
        missingness='Retain author literal values; any within-case field conflict quarantined, not resolved by majority vote.',
        outcomes='No effect tests or numerical response threshold inferred; source author category semantics need primary paper review.',
        patients='Case IDs are not by themselves reviewed independent-patient linkage.',
        sampling='Published post-treatment surgical immune-only compartment, not baseline prediction or full-TME fractions.'),indent=2)+'\n')
    m=pd.read_csv(path,usecols=['sampleID','major_cell_type','sub_cell_type']+FIELDS,keep_default_na=False,dtype=str)
    assert len(m)==1254749 and m.sampleID.nunique()==243 and m.sub_cell_type.nunique()==51
    assert m.groupby('sub_cell_type').major_cell_type.nunique().eq(1).all()
    rows=[];conflicts=[]
    for case,f in m.groupby('sampleID',sort=True):
        row=dict(sampleID=case,n_cells=len(f))
        for name in FIELDS:
            counts=f[name].value_counts(dropna=False);values=counts.index.tolist()
            row[name]=values[0] if len(values)==1 else None
            row[name+'_within_case_consistent']=len(values)==1
            if len(values)>1:conflicts.append(dict(sampleID=case,field=name,literal_values=json.dumps(values),cell_counts=json.dumps(counts.tolist())))
        rows.append(row)
    cases=pd.DataFrame(rows);cases.to_csv(O/'private_all243_author_case_fields.tsv',sep='\t',index=False,na_rep='CONFLICT')
    pd.DataFrame(conflicts,columns=['sampleID','field','literal_values','cell_counts']).to_csv(O/'private_within_case_metadata_conflicts.tsv',sep='\t',index=False)
    literal=[]
    for name in FIELDS:
        c=cases[name].fillna('WITHIN_CASE_CONFLICT').value_counts(dropna=False)
        for value,n in c.items():literal.append(dict(field=name,literal_value=value,case_count=int(n)))
    literal=pd.DataFrame(literal);literal.to_csv(O/'private_all_case_literal_field_counts.tsv',sep='\t',index=False)
    n=m.groupby(['sampleID','major_cell_type','sub_cell_type']).size().reset_index(name='n_cells')
    n.to_csv(O/'private_author_case_subtype_cell_counts.tsv',sep='\t',index=False)
    capture=[]
    for (major,sub),f in n.groupby(['major_cell_type','sub_cell_type'],sort=True):
        capture.append(dict(major=major,author_subtype=sub,native_cells=int(f.n_cells.sum()),cases_any=len(f),cases_ge20=int(f.n_cells.ge(20).sum()),
            smallest_nonzero_unit=int(f.n_cells.min()),median_nonzero_unit=float(f.n_cells.median())))
    capture=pd.DataFrame(capture);capture.to_csv(O/'private_all51_author_subtype_capture.tsv',sep='\t',index=False)
    assert capture.native_cells.sum()==len(m)
    ct=pd.crosstab([cases.cancer_type,cases['anti-PD1_therapy'],cases.chemotherapy],cases.pathological_response,dropna=False)
    ct.to_csv(O/'private_histology_exposure_author_response_crosstab.tsv',sep='\t')
    rates=cases[['pathological_response','pathological_response_rate']].copy()
    rates['numeric_rate']=pd.to_numeric(rates.pathological_response_rate,errors='coerce')
    rate_summary=rates.groupby('pathological_response',dropna=False).agg(cases=('pathological_response_rate','size'),numeric_cases=('numeric_rate','count'),minimum=('numeric_rate','min'),maximum=('numeric_rate','max'))
    rate_summary.to_csv(O/'private_author_response_rate_ranges.tsv',sep='\t')
    report=dict(status='ALL243_CASE_METADATA_AND51_AUTHOR_SUBTYPE_CAPTURE_REVIEWED_NO_ASSOCIATIONS_OR_PATIENT_LINKAGE_CLAIM',
        finished_UTC=datetime.now(timezone.utc).isoformat(),native_cells=len(m),author_cases=len(cases),author_subtypes=len(capture),
        within_case_conflict_records=len(conflicts),conflicts_by_field=dict(pd.Series([x['field'] for x in conflicts],dtype=str).value_counts().astype(int)),
        source_metadata_sha256=sha(path),source_paper='https://pubmed.ncbi.nlm.nih.gov/40147443/',
        files={p.name:sha(p) for p in O.glob('*.tsv')},original_labels_changed=False,old200_eligibility_reused=False,
        response_rate_meaning_assumed=False,patient_identity_reviewed=False,clinical_effects_tested=False)
    (O/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report));print('LITERAL_FIELDS',literal[~literal.field.isin(['age','pathological_response_rate'])].to_json(orient='records'))
    print('RATE_RANGES',rate_summary.reset_index().to_json(orient='records'))
    print('TREG_CAPTURE',capture[capture.author_subtype.str.contains('Treg')].to_json(orient='records'))


if __name__=='__main__':main()
