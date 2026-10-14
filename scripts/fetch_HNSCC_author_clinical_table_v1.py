"""Fetch only the publisher's small public clinical Table S1 to server.
Read-only workbook inspection; no workbook edits or patient row exports.
"""
from datetime import datetime,timezone
from pathlib import Path
import hashlib,json,os,urllib.request
import pandas as pd


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=root/'HNSCC_author_clinical_table_20261014_v1'
    assert not out.exists();out.mkdir()
    url='https://ars.els-cdn.com/content/image/1-s2.0-S2666379125004811-mmc2.xlsx'
    request=urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0 (scientific metadata retrieval)'})
    try:
        with urllib.request.urlopen(request,timeout=30) as response:
            b=response.read(2_000_001);final_url=response.url;content_type=response.headers.get('Content-Type')
        assert b.startswith(b'PK') and len(b)<=2_000_000,'Expected small XLSX, not login/HTML or expression data'
    except Exception as error:
        status=dict(status='PUBLIC_PUBLISHER_METADATA_FETCH_FAILED_NO_BYPASS',observed_UTC=datetime.now(timezone.utc).isoformat(),
                    source_url=url,error_type=type(error).__name__,error=str(error))
        (out/'failed_fetch.json').write_text(json.dumps(status,indent=2)+'\n')
        print('PUBLIC_PUBLISHER_METADATA_FETCH_FAILED_NO_BYPASS');return
    p=out/'author_Table_S1_original.xlsx';p.write_bytes(b)
    workbook=pd.ExcelFile(p);sheets=[]
    for name in workbook.sheet_names:
        frame=pd.read_excel(workbook,sheet_name=name,header=None)
        # Headers/definitions are retained privately; no row values emitted.
        frame.to_csv(out/(name.replace('/','_')+'_private_raw_cells.tsv'),sep='\t',index=False,header=False)
        sheets.append(dict(sheet=name,rows=len(frame),columns=len(frame.columns)))
    provenance=dict(status='AUTHOR_CLINICAL_TABLE_FETCHED_AND_READ_ONLY_INSPECTED',observed_UTC=datetime.now(timezone.utc).isoformat(),
                    source_url=url,final_url=final_url,doi='10.1016/j.xcrm.2025.102408',
                    table_definition='Table S1: clinical metadata and analyzed cell counts',
                    bytes=len(b),content_type=content_type,sha256=hashlib.sha256(b).hexdigest(),sheets=sheets,
                    original_workbook_unchanged=True,patient_rows_not_exported=True)
    (out/'private_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print('AUTHOR_CLINICAL_TABLE_FETCHED_AND_READ_ONLY_INSPECTED_RESULTS_REMAIN_SERVER')


if __name__=='__main__':main()
