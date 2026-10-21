"""Inspect generic author code for explicit registry-link semantics; never run it.

Individual literals are masked. Only generic code/schema and source category
counts are written. This does not infer identities from coded prefixes.
"""
from pathlib import Path
from datetime import datetime,timezone
from collections import Counter
import csv,hashlib,importlib.util,io,json,os,re,traceback,zipfile


def sha(b):return hashlib.sha256(b).hexdigest()


def module(name,path,expected):
    assert sha(path.read_bytes())==expected
    spec=importlib.util.spec_from_file_location(name,path)
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result)
    return result


def main():
    helper=Path(os.environ['PANCANCER_IMC_TRANSPORT_HELPER'])
    source=module('published_archive_transport',helper,'9fb200d98e321584fde8bca6cb92f72151ac33ec00d64e04d93e7e2b33e870e2')
    retry_path=helper.parent.parent/'pancancer-CM02-melanoma-ROI-registry-20261021-v2'/'qualify_melanoma_IMC_ROI_registry_v2.py'
    retry=module('bounded_archive_transport',retry_path,'190b81b9c771d06e1b3b63506a77c13ef29202f96553b456973796cec5062910')
    out=Path(os.environ['PANCANCER_RAM_OUTPUT']);assert not out.exists();out.mkdir(mode=0o700)
    def save(n,v):(out/n).write_text(json.dumps(v,indent=2)+'\n')
    save('plan_before_author_code_review.json',dict(frozen_UTC=datetime.now(timezone.utc).isoformat(),
        producer_SHA256=sha(Path(__file__).read_bytes()),scope='Four complete generic author code files and the complete small ROI registry; no cell data.',
        no_author_code_execution=True,no_case_or_ROI_literals_or_records_written=True,
        no_inferred_patient_rules_or_new_effects=True))
    raw=source.fetch(source.BASE+'api/records/6838169',2*1024**2);record=json.loads(raw)
    f=next(x for x in record['files'] if x['key']=='Melanoma_code_data.zip')
    assert f['size']==668227028 and f['checksum']=='md5:da74b7228763b76dcd2d9830c1a10cd9'
    rf=retry.BoundedRetryRangeFile(source,f['links']['self'],f['size'])
    safe={'sample','sample_id','sampleROI','sampleROI2','ROI','roi','patient','patient_id','patientID','region','response','Position','position','x','y','X','Y','CT','IM','R','NR','invasive','core','central','margin','tumor','invasive margin','core tumor','Treg','CD4T','CD8T','Treg cells','CD4 T cells','CD8 T cells'}
    quote=re.compile(r'''(["'])(?:\\.|(?!\1).)*?\1''')
    def redact(m):
        s=m.group()[1:-1]
        return m.group() if s in safe or re.fullmatch(r'[_ ,;:|.\-]{1,4}',s) else '"[masked_literal]"'
    proofs=[]
    with zipfile.ZipFile(rf) as z:
        for name in ['Rscript/Spatial_analysis.Rmd','Rscript/functions.R','Rscript/generate_Figs.Rmd','Rscript/Clustering.Rmd']:
            info=z.getinfo(name);assert info.file_size<1024**2
            b=z.read(info);text=b.decode('utf-8-sig');lines=text.splitlines();chosen=[]
            for i,line in enumerate(lines):
                if re.search(r'patient|sampleROI|Position|separate\(|colnames\(|read\.csv|centroid|positive|decoy',line,re.I):
                    chosen.append(dict(line=i+1,masked_generic_code=quote.sub(redact,line)[:600]))
            proofs.append(dict(generic_path=name,complete_member_SHA256=sha(b),CRC_checked=True,selected_masked_lines=chosen))
        info=z.getinfo('data/input/meta/ROI_meta_afterQC.csv');b=z.read(info)
        assert sha(b)=='d0fa2438b7b91157a7848654fb1f23e40c9dbcca8020f7f5ff0566ce7acede21'
        rows=list(csv.DictReader(io.StringIO(b.decode('utf-8-sig'))))
        categories={}
        allowed_region={'invasive','core','central','margin','tumor','invasive margin','core tumor','invasive_margin','core_tumor','CT','IM'}
        counts=Counter(r['region'] for r in rows)
        if all(x.lower() in allowed_region or x in allowed_region for x in counts):categories=dict(counts)
        key_profile=dict(field='sampleROI2',nonempty=sum(bool(r['sampleROI2']) for r in rows),unique=len({r['sampleROI2'] for r in rows}),
            delimiter_counts={sep:dict(Counter(r['sampleROI2'].count(sep) for r in rows)) for sep in ['_','-','/']},
            no_values_or_prefixes_exported=True)
    final=json.loads(source.fetch(source.BASE+'api/records/6838169',2*1024**2))
    other=next(x for x in final['files'] if x['key']==f['key'])
    assert all(f.get(k)==other.get(k) for k in ['id','size','checksum','links'])
    save('private_complete_masked_link_code_review.json',dict(status='COMPLETE_MASKED_AUTHOR_LINK_CODE_AND_ROI_SCHEMA_REVIEW',
        completed_UTC=datetime.now(timezone.utc).isoformat(),producer_SHA256=sha(Path(__file__).read_bytes()),
        generic_author_code=proofs,ROI_registry_SHA256=sha(b),region_source_category_counts=categories,
        ROI_key_schema_profile=key_profile,published_file_metadata_unchanged=True,
        no_case_rule_or_patient_acceptance_inferred=True,no_cell_rows_or_clinical_effects_read=True,
        requested_range_bytes=rf.requested,successful_range_bytes=rf.transferred,failed_transport_classes=rf.failures))
    save('private_terminal.json',dict(status='COMPLETE_MASKED_AUTHOR_LINK_CODE_AND_ROI_SCHEMA_REVIEW',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(Path(__file__).read_bytes()),
        files={p.name:dict(SHA256=sha(p.read_bytes()),bytes=p.stat().st_size) for p in out.iterdir()}))
    print(json.dumps(dict(status='COMPLETE_MASKED_AUTHOR_LINK_CODE_AND_ROI_SCHEMA_REVIEW',terminal_SHA256=sha((out/'private_terminal.json').read_bytes()))))


if __name__=='__main__':
    try:main()
    except Exception as error:
        out=Path(os.environ['PANCANCER_RAM_OUTPUT']);out.mkdir(mode=0o700,parents=True,exist_ok=True)
        payload=dict(status='FAILED_NOT_ACCEPTED',observed_UTC=datetime.now(timezone.utc).isoformat(),error_class=type(error).__name__,
            frames=[dict(file=Path(x.filename).name,line=x.lineno,function=x.name) for x in traceback.extract_tb(error.__traceback__)],no_individual_values_in_errors=True)
        (out/'private_failed_terminal.json').write_text(json.dumps(payload,indent=2)+'\n')
        print(json.dumps(dict(status=payload['status'],error_class=payload['error_class'])))
        raise SystemExit(1)
