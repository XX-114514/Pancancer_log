"""Small official clinical-covariate provenance; records stay in RAM.

Use only exact patient fields/keys and recognized clinical columns. Coded
headers are suppressed. No survival endpoints, expression or outcome models.
"""
from pathlib import Path
from datetime import datetime,timezone
from collections import defaultdict
import csv,hashlib,importlib.util,io,json,math,os,re,traceback,zipfile
import pandas as pd


def sha(b):return hashlib.sha256(b).hexdigest()


def load(name,path,expected):
    assert sha(path.read_bytes())==expected
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def field_role(value):
    key=re.sub('[^a-z0-9]','',value.lower())
    roles={'age':{'age','ageyears','ageatdiagnosis','ageattreatment'},'sex':{'sex','gender'},
        'stage':{'stage','tnmstage','clinicalstage'},'site':{'site','tumorsite','tumoursite'}}
    return next((role for role,keys in roles.items() if key in keys),None)



def author_patient_keys(patients,prefix):
    # Literal patient branch of the SHA-pinned author function; no eval/R.
    assert re.fullmatch('[A-Za-z_]{1,16}',prefix)
    mapping={p:prefix+(p if 'F' in p else p+'F') for p in patients}
    assert len(set(mapping.values()))==len(mapping)
    return mapping


def main():
    out=Path(os.environ['PANCANCER_RAM_OUTPUT']);assert not out.exists();out.mkdir(mode=0o700)
    def save(n,v):(out/n).write_text(json.dumps(v,indent=2)+'\n')
    save('plan_before_clinical_covariate_metadata.json',dict(frozen_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(Path(__file__).read_bytes()),
        selected_member='data/input/clidata/mel_puch_cli.csv',maximum_complete_member_bytes=65536,
        scope='Recognized age/sex/stage/site fields and exact overlap with explicit imaging patient values; no effect/endpoint/survival models.',
        row_key_and_transposed_explicit_column_key_orientations_both_checked=True,
        no_prefix_or_numeric_identity_inference=True,author_defined_literal_join_format_reviewed_before_use=True,no_coded_headers_rows_or_individual_arrays_written=True,
        no_expression_or_learned_arrays_selected=True))
    helper=Path(os.environ['PANCANCER_IMC_TRANSPORT_HELPER'])
    source=load('published_archive_transport',helper,'9fb200d98e321584fde8bca6cb92f72151ac33ec00d64e04d93e7e2b33e870e2')
    transport=load('bounded_archive_transport',Path(os.environ['PANCANCER_IMC_BOUNDED_TRANSPORT']),'709f9b5ced875fea663557bad936f32d004e12d5daf291048b7de476064b0d68')
    d=json.loads(source.fetch(source.BASE+'api/records/6838169',2*1024**2));f=next(x for x in d['files'] if x['key']=='Melanoma_code_data.zip')
    assert f['size']==668227028 and f['checksum']=='md5:da74b7228763b76dcd2d9830c1a10cd9'
    rf=transport.RangeFile(source,f['links']['self'],f['size'])
    with zipfile.ZipFile(rf) as z:
        a=z.read('data/input/hc_result/hc_result_sample_anno_invasive.csv')
        assert sha(a)=='59c9f04e4e12cc7484f1cc262b96b77c7a5eae4efff7d2510fbc23921b6f55dc'
        bridge=[dict(ROI=r['ROI'],patient=r['patient']) for r in csv.DictReader(io.StringIO(a.decode('utf-8-sig')))]
        patients={r['patient'] for r in bridge};assert all(patients)
        info=z.getinfo('data/input/clidata/mel_puch_cli.csv');assert info.file_size<=65536
        b=z.read(info);rows=list(csv.DictReader(io.StringIO(b.decode('utf-8-sig'))));columns=list(rows[0])
        df=pd.read_csv(io.BytesIO(b),dtype=str,keep_default_na=False)
        assert len(rows)==len(df)
        for i,c in enumerate(columns):assert [r[c] for r in rows]==df.iloc[:,i].tolist()
        opportunities=[]
        for key in columns:
            values=[r[key] for r in rows]
            if len(values)!=len(set(values)):continue
            records={r[key]:r for r in rows};overlap=patients&set(records)
            if not overlap:continue
            profiles=[]
            for col in columns:
                role=field_role(col)
                if role is None:continue
                selected=[records[p][col] for p in overlap]
                profile=dict(role=role,source_field=col,matched_patient_units=len(overlap),
                    nonempty_values=sum(v!='' for v in selected),distinct_values=len(set(selected)))
                if role=='age':
                    valid=[]
                    for v in selected:
                        try:
                            value=float(v)
                            if math.isfinite(value):valid.append(value)
                        except ValueError:pass
                    independent=[v for v in pd.to_numeric(pd.Series(selected),errors='coerce').dropna().tolist() if math.isfinite(v)]
                    assert sorted(valid)==sorted(independent)
                    profile.update(numeric_values=len(valid),adult_year_range_18_to120_values=sum(18<=v<=120 for v in valid),
                        age_units_not_inferred_from_range_alone=True)
                profiles.append(profile)
            opportunities.append(dict(orientation='patient_rows',key_field=key if key in ['', 'patient','patient_id','Patient','PatientID','sample','Sample'] else 'noncontract_key_field_name_suppressed',
                exact_matched_imaging_patients=len(overlap),all_imaging_patients_covered=overlap==patients,covariate_profiles=profiles))

        functions=z.read('Rscript/functions.R');figures=z.read('Rscript/generate_Figs.Rmd')
        assert sha(functions)=='bfb6f99c80947643655aa5d4d8ef310e5f0a1c24aac86326482cb32c123565ee'
        assert sha(figures)=='b213b123c189457f7870d042094f85415f9dfef9bfc52e296393bfcc9cabf4af'
        fn=functions.decode('utf-8-sig');fig=figures.decode('utf-8-sig')
        assert "if (stringr::str_detect(each, 'F'))" in fn and "a = paste0(each,'F')" in fn
        assert "if (IDtype == 'patient')" in fn
        lines=fig.splitlines()
        assert 'mel_puch_cli.csv' in lines[363] and 'row.names = 1' in lines[363]
        assert "renamePID(exp.sample$patient,IDtype = 'patient')" in lines[374]
        # Read the one literal constant used by the author for this clinical join.
        prefix_match=re.fullmatch(r'\s*rownames\(exp.sample\) <- paste0\([\'\"]([A-Za-z_]{1,16})[\'\"],rownames\(exp.sample\)\)\s*',lines[375])
        assert prefix_match is not None
        prefix=prefix_match.group(1);normalized=author_patient_keys(patients,prefix)
        first=columns[0];records={r[first]:r for r in rows};assert len(records)==len(rows)
        overlap={p for p,v in normalized.items() if v in records};profiles=[]
        for col in columns:
            role=field_role(col)
            if role is None:continue
            selected=[records[normalized[p]][col] for p in overlap]
            profile=dict(role=role,source_field=col,matched_patient_units=len(overlap),
                nonempty_values=sum(v!='' for v in selected),distinct_values=len(set(selected)))
            if role=='age':
                valid=[]
                for v in selected:
                    try:
                        number=float(v)
                        if math.isfinite(number):valid.append(number)
                    except ValueError:pass
                other=[v for v in pd.to_numeric(pd.Series(selected),errors='coerce').dropna().tolist() if math.isfinite(v)]
                assert sorted(valid)==sorted(other)
                profile.update(numeric_values=len(valid),adult_year_range_18_to120_values=sum(18<=v<=120 for v in valid),age_units_not_inferred_from_range_alone=True)
            profiles.append(profile)
        opportunities.append(dict(orientation='author_declared_first_column_clinical_index',
            exact_matched_imaging_patients=len(overlap),all_imaging_patients_covered=overlap==patients,
            author_patient_formatter='Preserve names containing F; otherwise append F; then prepend the author literal clinical-join prefix.',
            generic_author_prefix=prefix,source_author_functions_SHA256=sha(functions),source_author_figures_SHA256=sha(figures),
            explicit_author_clinical_first_index_and_patient_branch_checked=True,normalized_imaging_keys_unique=True,
            no_case_ids_or_heuristically_chosen_prefixes=True,covariate_profiles=profiles))

        # Some public clinical tables store coded patients as column headers.
        overlap=patients&set(columns)
        if overlap:
            trait_key=columns[0];profiles=[]
            for row in rows:
                role=field_role(row[trait_key])
                if role is None:continue
                selected=[row[p] for p in overlap]
                profiles.append(dict(role=role,source_field=row[trait_key],matched_patient_units=len(overlap),
                    nonempty_values=sum(v!='' for v in selected),distinct_values=len(set(selected))))
            opportunities.append(dict(orientation='explicit_patient_columns',exact_matched_imaging_patients=len(overlap),
                all_imaging_patients_covered=overlap==patients,covariate_profiles=profiles))
    final=json.loads(source.fetch(source.BASE+'api/records/6838169',2*1024**2));other=next(x for x in final['files'] if x['key']==f['key'])
    assert all(f.get(k)==other.get(k) for k in ['id','size','checksum','links']);rf.close()
    save('private_complete_clinical_covariate_key_qualification.json',dict(status='COMPLETE_SMALL_CLINICAL_COVARIATE_LITERAL_KEY_QUALIFICATION',
        completed_UTC=datetime.now(timezone.utc).isoformat(),producer_SHA256=sha(Path(__file__).read_bytes()),
        complete_member_SHA256=sha(b),complete_CRC_checked=True,member_bytes=len(b),source_rows=len(rows),
        safe_recognized_column_fields=[c for c in columns if field_role(c) is not None],other_column_count=sum(field_role(c) is None for c in columns),
        exact_key_opportunities=opportunities,CSV_Pandas_full_string_value_replay=True,
        no_coded_headers_or_record_values_written=True,no_survival_or_response_association_computed=True,
        bulk_patient_cohort_not_claimed_independent_imaging_validation=True))
    save('private_terminal.json',dict(status='COMPLETE_SMALL_CLINICAL_COVARIATE_LITERAL_KEY_QUALIFICATION',completed_UTC=datetime.now(timezone.utc).isoformat(),
        source_SHA256=sha(Path(__file__).read_bytes()),files={p.name:dict(SHA256=sha(p.read_bytes()),bytes=p.stat().st_size) for p in out.iterdir()}))
    print(json.dumps(dict(status='COMPLETE_SMALL_CLINICAL_COVARIATE_LITERAL_KEY_QUALIFICATION',terminal_SHA256=sha((out/'private_terminal.json').read_bytes()))))


if __name__=='__main__':
    try:main()
    except Exception as error:
        out=Path(os.environ['PANCANCER_RAM_OUTPUT']);out.mkdir(mode=0o700,parents=True,exist_ok=True)
        value=dict(status='FAILED_NOT_ACCEPTED',observed_UTC=datetime.now(timezone.utc).isoformat(),error_class=type(error).__name__,
            frames=[dict(file=Path(x.filename).name,line=x.lineno,function=x.name) for x in traceback.extract_tb(error.__traceback__)],no_record_values_in_error=True)
        (out/'private_failed_terminal.json').write_text(json.dumps(value,indent=2)+'\n')
        print(json.dumps(dict(status=value['status'],error_class=value['error_class'])))
        raise SystemExit(1)
