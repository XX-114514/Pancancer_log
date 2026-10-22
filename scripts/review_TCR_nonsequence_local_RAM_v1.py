"""Local-only public annotation review. No SSH or human-record export.

Only six nonsequence columns enter the parser. Only source summaries are saved.
The complete compressed file is never sent to another host or saved to disk.
"""
from pathlib import Path
from datetime import datetime,timezone
from urllib.request import Request,ProxyHandler,build_opener
import hashlib
import io
import json
import os
import resource
import numpy as np
import pandas as pd
import qualify_NSCLC_reported_TCR_expansion_RAM_v1 as fields


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def now():return datetime.now(timezone.utc).isoformat()
def save(path,obj):
    assert not path.exists();path.write_text(json.dumps(obj,indent=2)+'\n')


def summarize(frame,plan):
    before=len(frame);frame=frame.drop_duplicates(fields.SELECTED).copy()
    exact_duplicates=before-len(frame);ambiguous=frame.duplicated(['sampleID','cellID'],keep=False)
    ambiguous_rows=int(ambiguous.sum());frame=frame.loc[~ambiguous].copy().reset_index(drop=True)
    number,expanded,known,profile=fields.expansion_qualification(frame)
    allowed=plan['author_lineage_names'];labels=allowed['Treg']+allowed['CD8T']+allowed['CD4T_nonTreg']
    axes=[];pairs=[]
    for axis in ['sub_cell_type','T_new_name']:
        groups={label:[label] for label in labels}
        groups.update(Treg_POOL=allowed['Treg'],CD8_POOL=allowed['CD8T'],other_CD4_POOL=allowed['CD4T_nonTreg'])
        grouped={}
        cases=sorted(frame.sampleID.unique());codes=pd.Index(cases).get_indexer(frame.sampleID)
        for group,names in groups.items():
            use=frame[axis].isin(names).to_numpy();called=use&known
            counts=np.bincount(codes[called],minlength=len(cases));positive=np.bincount(codes[called],weights=expanded[called],minlength=len(cases)).astype(int)
            grouped[group]=counts
            oracle=frame.loc[called,['sampleID']].assign(expanded=expanded[called].astype(int)).groupby('sampleID').expanded.sum().reindex(cases,fill_value=0).to_numpy(int)
            np.testing.assert_array_equal(positive,oracle)
            axes.append(dict(annotation_axis=axis,author_group=group,source_records=int(use.sum()),known_expansion_records=int(called.sum()),
                expanded_records=int(positive.sum()) if profile['operational_reported_expansion_qualification'] else None,
                materials_any=int((counts>0).sum()),materials_ge20=int((counts>=20).sum()),
                annotation_not_linked_to_original_RNA_or_clinical=True))
        for r in ['Treg_POOL']+allowed['Treg']:
            for e in allowed['CD8T']:
                pairs.append(dict(annotation_axis=axis,R_group=r,E_group=e,
                    materials_R_E_and_otherCD4_TCR_ge20=int(((grouped[r]>=20)&(grouped[e]>=20)&(grouped['other_CD4_POOL']>=20)).sum()),
                    no_RNA_QC_clinical_or_patient_qualification=True))
    return dict(status='COMPLETE_LOCAL_ONLY_NONSEQUENCE_TCR_STRUCTURE_NO_RNA_OR_CLINICAL_LINKAGE',source='GSE243013',
        original_source_rows=before,exact_selected_field_duplicates_collapsed=exact_duplicates,ambiguous_cell_key_rows_quarantined=ambiguous_rows,
        source_materials=int(frame.sampleID.nunique()),semantics=profile,all_fixed_author_group_capture=axes,all_fixed_pair_capture=pairs,
        both_annotation_axes_reported_without_selecting_favorable_mapping=True,
        no_CDR3_original_clonotype_or_receptor_gene_columns_selected=True,
        no_raw_file_identifier_or_individual_array_disk_or_network_export=True,
        TRG_TRD_not_assayed_by_these_columns_not_gamma_delta_validation=True,
        no_hypothesis_test_independent_patient_or_original_CM_validation=True),frame


def main():
    task=Path(os.environ['PANCANCER_LOCAL_TCR_TASK']);code=Path(__file__).parent;frozen=task/'private_frozen_local_measurement_v1.json';plan=json.loads(frozen.read_text());item=plan['source']
    for name,digest in plan['code_SHA256'].items():assert sha(code/name)==digest
    fields.synthetic();started=now();phase='OFFICIAL_HTTP_HEAD'
    try:
        opener=build_opener(ProxyHandler({'https':os.environ['PANCANCER_METADATA_HTTP_PROXY']}))
        with opener.open(Request(item['url'],method='HEAD'),timeout=35) as response:
            assert response.status==200 and response.url==item['url'] and int(response.headers['Content-Length'])==item['bytes']
            assert response.headers.get('ETag')==item['etag'] and response.headers.get('Last-Modified')==item['last_modified']
        phase='OFFICIAL_PUBLIC_GET_LOCAL_RAM_ONLY';buffer=io.BytesIO();digest=hashlib.sha256()
        with opener.open(Request(item['url']),timeout=45) as response:
            assert response.status==200 and response.url==item['url'] and int(response.headers['Content-Length'])==item['bytes']
            assert response.headers.get('ETag')==item['etag'] and response.headers.get('Last-Modified')==item['last_modified']
            while True:
                block=response.read(1024**2)
                if not block:break
                buffer.write(block);digest.update(block);assert buffer.tell()<=item['bytes']
        assert buffer.tell()==item['bytes'];phase='LOCAL_NONSEQUENCE_COLUMNS'
        frame=fields.parse_selected(buffer.getvalue(),plan);buffer.close();summary,frame=summarize(frame,plan)
        assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss<plan['local_RSS_cap_bytes']
        save(task/'private_complete_local_measurement_v1.json',summary)
        save(task/'private_local_terminal_v1.json',dict(status=summary['status'],started_UTC=started,completed_UTC=now(),
            fetched_public_object_SHA256=digest.hexdigest(),frozen_contract_SHA256=sha(frozen),summary_SHA256=sha(task/'private_complete_local_measurement_v1.json'),
            no_remote_process_or_record_transmission=True,no_raw_original_file_sequence_or_individual_records_saved=True,
            full_public_gzip_CRC_size_headers_and_private_nonsequence_field_fixture_accepted=True))
        del frame
        print(json.dumps(dict(status='COMPLETE_LOCAL_ONLY_TCR_PREFLIGHT_RAM_RELEASED',operational_reported_expansion_qualified=summary['semantics']['operational_reported_expansion_qualification'])),flush=True)
    except Exception as error:
        save(task/'private_local_failed_terminal_v1.json',dict(status='LOCAL_PREFLIGHT_NOT_ACCEPTED',observed_UTC=now(),phase=phase,
            error_type=type(error).__name__,HTTP_code=getattr(error,'code',None),no_raw_errors_or_partial_records_emitted=True))
        raise RuntimeError('LOCAL_PREFLIGHT_NOT_ACCEPTED') from None


if __name__=='__main__':main()
