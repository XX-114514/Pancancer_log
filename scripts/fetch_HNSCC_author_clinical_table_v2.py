"""Fetch Table S1 via documented Europe PMC OA supplementary service.
The publisher transport failure is preserved. No access-control bypass.
Only Table S1 is extracted; no expression matrices or patient exports.
"""
from datetime import datetime,timezone
from pathlib import Path
import hashlib,io,json,os,urllib.request,zipfile
import pandas as pd


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=root/'HNSCC_author_clinical_table_20261014_v2'
    assert not out.exists();out.mkdir()
    url='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12629806/supplementaryFiles'
    try:
        with urllib.request.urlopen(url,timeout=45) as response:
            b=response.read(64*1024*1024+1);final_url=response.url
        assert b.startswith(b'PK') and len(b)<=64*1024*1024,'Expected bounded public supplement ZIP'
        with zipfile.ZipFile(io.BytesIO(b)) as archive:
            names=[n for n in archive.namelist() if Path(n).name.lower() in ['mmc2.xlsx','1-s2.0-s2666379125004811-mmc2.xlsx']]
            assert len(names)==1,'Unambiguous author Table S1 required'
            data=archive.read(names[0]);assert data.startswith(b'PK') and len(data)<2_000_000
    except Exception as error:
        (out/'failed_fetch.json').write_text(json.dumps(dict(status='PUBLIC_OA_METADATA_FETCH_FAILED_NO_BYPASS',
          observed_UTC=datetime.now(timezone.utc).isoformat(),source_url=url,error_type=type(error).__name__,error=str(error)),indent=2)+'\n')
        print('PUBLIC_OA_METADATA_FETCH_FAILED_NO_BYPASS');return
    p=out/'author_Table_S1_original.xlsx';p.write_bytes(data)
    workbook=pd.ExcelFile(p);sheets=[]
    for name in workbook.sheet_names:
        frame=pd.read_excel(workbook,sheet_name=name,header=None)
        frame.to_csv(out/(name.replace('/','_')+'_private_raw_cells.tsv'),sep='\t',index=False,header=False)
        sheets.append(dict(sheet=name,rows=len(frame),columns=len(frame.columns)))
    provenance=dict(status='AUTHOR_CLINICAL_TABLE_FETCHED_AND_READ_ONLY_INSPECTED',observed_UTC=datetime.now(timezone.utc).isoformat(),
      source_url=url,final_url=final_url,api_documentation='https://europepmc.org/RestfulWebService',
      doi='10.1016/j.xcrm.2025.102408',table_definition='Table S1: clinical metadata and analyzed cell counts',
      bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),archive_sha256=hashlib.sha256(b).hexdigest(),sheets=sheets,
      original_workbook_unchanged=True,patient_rows_not_exported=True)
    (out/'private_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print('AUTHOR_CLINICAL_TABLE_FETCHED_AND_READ_ONLY_INSPECTED_RESULTS_REMAIN_SERVER')


if __name__=='__main__':main()
