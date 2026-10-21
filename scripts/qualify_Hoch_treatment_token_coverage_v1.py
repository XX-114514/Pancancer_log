"""Whitelist-token input qualification only, not regimen or clinical acceptance.

Unknown free text is suppressed. Counts are metadata rows, not independent
patients. Drug abbreviations/brands are reported as tokens, without inferring
concurrent/sequential schedules or treatment eligibility.
"""
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import argparse
import csv
import hashlib
import importlib.util
import io
import json
import os
import re
import zipfile

TOKENS=['nivolumab','nivo','opdivo','pembrolizumab','pembro','keytruda',
        'ipilimumab','ipi','yervoy','atezolizumab','atezo','durvalumab',
        'durva','avelumab','pd1','pd-1','ctla4','ctla-4','immunotherapy']
DIRECTORY_SHA='dff220e35d38da2b60c42950190d71b2259f12aae3c21ad4bc6c05ad5c5d2893'
CLINICAL_SHA='ac4a02e4d0e23504aa980b21585de74c4d7390d2e6e0ec2399303c5c29ca5b25'


def summarize(rows):
    summaries=[]
    for field in ['Drug_after_surgery','Treatment_after_surgery','Name_last_drug_used','Last_sys_treatment_before_surgery']:
        counts=Counter()
        for row in rows:
            text=row.get(field,'').strip().lower()
            found=[token for token in TOKENS if re.search(r'(?<![a-z0-9])'+re.escape(token)+r'(?![a-z0-9])',text)]
            counts['|'.join(found) if found else 'NO_WHITELIST_TOKEN_OR_MISSING']+=1
        summaries.append(dict(field=field,source_metadata_row_token_patterns=dict(counts)))
    return dict(status='COMPLETE_TREATMENT_TOKEN_INPUT_QUALIFICATION_ONLY',field_coverage=summaries,
        counts_are_metadata_rows_not_independent_patients=True,unknown_free_text_suppressed=True,
        tokens_not_equated_with_regimen_or_schedule=True,no_date_patient_outcome_or_cell_rows_exported=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--scripts-root',type=Path,required=True);a=p.parse_args()
    assert not a.output.exists()
    path=a.scripts_root/'melanoma_IMC_existing_clash_transport_v1.py'
    assert hashlib.sha256(path.read_bytes()).hexdigest()=='546b491e20e6797e5e255bfa226fa568a83dca08eb01334a27128fb0587755ac'
    spec=importlib.util.spec_from_file_location('token_small_metadata_transport',path)
    transport=importlib.util.module_from_spec(spec);spec.loader.exec_module(transport)
    transport.PROXY_PORT=int(os.environ['PANCANCER_EXISTING_CLASH_PORT'])
    rf=transport.RangeFile(None,'https://zenodo.org/api/records/5994136/files/data_for_analysis.zip/content',7123744479);rf.max_chunk=64*1024
    def directory(z):return hashlib.sha256(json.dumps([(x.filename,x.file_size,x.compress_size,x.CRC) for x in z.infolist()],separators=(',',':')).encode()).hexdigest()
    try:
        with zipfile.ZipFile(rf) as z:
            assert directory(z)==DIRECTORY_SHA
            raw=z.read('data_for_analysis_zenodo/protein/clinical_data_protein.csv');assert hashlib.sha256(raw).hexdigest()==CLINICAL_SHA
            selected=[{key:row.get(key,'') for key in ['Drug_after_surgery','Treatment_after_surgery','Name_last_drug_used','Last_sys_treatment_before_surgery']}
                for row in csv.DictReader(io.StringIO(raw.decode('utf-8-sig')))]
            result=summarize(selected)
        with zipfile.ZipFile(rf) as replay:assert directory(replay)==DIRECTORY_SHA
        result.update(completed_UTC=datetime.now(timezone.utc).isoformat(),clinical_member_SHA256=CLINICAL_SHA,
            full_small_member_CRC_and_fresh_directory_replay=True,wire_bytes=rf.transferred,
            producer_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        a.output.write_text(json.dumps(result,indent=2)+'\n')
    finally:rf.close()


if __name__=='__main__':main()
