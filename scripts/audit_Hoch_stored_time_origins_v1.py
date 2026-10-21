"""Two declared calendar anchors, days only; no outcome association or repair.

Both names are explicit in the public author preparation code. The audit
reports agreement, never selects a winning clinical endpoint or edits data.
All image/block/patient records remain process RAM, outputs source aggregates.
"""
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import argparse,csv,hashlib,importlib.util,io,json,os,zipfile
import numpy as np
import qualify_Hoch_block_landmark_scope_v1 as scope


def audit(rows):
    groups=defaultdict(list)
    for row in rows:
        r={k:scope.canonical(row.get(k)) for k in scope.FIELDS}
        if all(r[k] for k in ['PatientID','BlockID','ImageNumber']):groups[r['PatientID'],r['BlockID']].append(r)
    contexts=[];conflicts=0
    for records in groups.values():
        values={f:{r[f] for r in records if r[f]} for f in scope.FIELDS if f not in ['PatientID','BlockID','ImageNumber']}
        if any(len(v)>1 for v in values.values()):conflicts+=1;continue
        contexts.append({f:next(iter(v)) if v else '' for f,v in values.items()})
    result=[]
    for endfield,storedfield in [('Date_death','Time_to_death_or_last_PET'),('Date_progression','Time_to_progression_or_last_PET')]:
        for anchorfield in ['Date_surgery','Start_treatment']:
            counts=Counter();differences=[];ratios=[]
            for r in contexts:
                start=scope.dates.explicit_date(r[anchorfield]);end=scope.dates.explicit_date(r[endfield] or r['Last_PET'])
                if start is None or end is None or not r[storedfield]:counts['UNAVAILABLE_DATE_OR_STORED_VALUE']+=1;continue
                try:stored=float(r[storedfield])
                except ValueError:counts['UNPARSEABLE_STORED_VALUE']+=1;continue
                if not np.isfinite(stored):counts['UNPARSEABLE_STORED_VALUE']+=1;continue
                expected=(end-start).days;delta=stored-expected;differences.append(delta)
                if expected:ratios.append(stored/expected)
                counts['EXACT_DAYS_AGREEMENT' if abs(delta)<1e-8 else 'DAYS_MISMATCH']+=1
                if abs(stored-round(stored))<1e-8:counts['STORED_VALUE_IS_INTEGER']+=1
            result.append(dict(endpoint_field=endfield,stored_duration_field=storedfield,declared_calendar_anchor=anchorfield,
                audit_counts=dict(counts),stored_minus_recomputed_days_min_q25_median_q75_max=np.quantile(differences,[0,.25,.5,.75,1]).tolist() if len(differences)>=5 else None,
                ratio_to_recomputed_days_min_q25_median_q75_max=np.quantile(ratios,[0,.25,.5,.75,1]).tolist() if len(ratios)>=5 else None))
    return dict(status='COMPLETE_DECLARED_CALENDAR_ANCHOR_AUDIT_NO_ENDPOINT_REPAIR',accepted_block_contexts=len(contexts),within_block_conflicts=conflicts,
        rows=result,no_new_unit_guesses_or_stored_value_replacement=True,no_clinical_tests_or_cell_objects=True,no_individual_values_exported=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--scripts-root',type=Path,required=True);a=p.parse_args();assert not a.output.exists()
    path=a.scripts_root/'melanoma_IMC_existing_clash_transport_v1.py';assert hashlib.sha256(path.read_bytes()).hexdigest()=='546b491e20e6797e5e255bfa226fa568a83dca08eb01334a27128fb0587755ac'
    spec=importlib.util.spec_from_file_location('bounded_calendar_audit_reader',path);t=importlib.util.module_from_spec(spec);spec.loader.exec_module(t);t.PROXY_PORT=int(os.environ['PANCANCER_EXISTING_CLASH_PORT'])
    rf=t.RangeFile(None,'https://zenodo.org/api/records/5994136/files/data_for_analysis.zip/content',7123744479);rf.max_chunk=64*1024
    def directory(z):return hashlib.sha256(json.dumps([(x.filename,x.file_size,x.compress_size,x.CRC) for x in z.infolist()],separators=(',',':')).encode()).hexdigest()
    try:
        with zipfile.ZipFile(rf) as z:
            assert directory(z)==scope.DIRECTORY_SHA
            raw=z.read('data_for_analysis_zenodo/protein/clinical_data_protein.csv');assert len(raw)<256*1024 and hashlib.sha256(raw).hexdigest()==scope.CLINICAL_SHA
            reader=csv.DictReader(io.StringIO(raw.decode('utf-8-sig')));assert set(scope.FIELDS)<=set(reader.fieldnames)
            rows=[{k:r[k] for k in scope.FIELDS} for r in reader];result=audit(rows)
        with zipfile.ZipFile(rf) as replay:assert directory(replay)==scope.DIRECTORY_SHA
        result.update(completed_UTC=datetime.now(timezone.utc).isoformat(),clinical_member_SHA256=scope.CLINICAL_SHA,fresh_directory_and_full_small_member_CRC=True,wire_bytes=rf.transferred,producer_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    finally:rf.close()


if __name__=='__main__':main()
