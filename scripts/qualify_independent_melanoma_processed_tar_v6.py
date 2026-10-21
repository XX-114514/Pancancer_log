"""Read a fixed public processed archive in RAM; qualify schema, never outcomes.

Author code is untrusted source material and is never executed. No archive,
individual table, image, sequence or expression array is saved. A separately
hashed own callback may use the accepted in-memory archive after qualification.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import csv
import gzip
import bz2
import lzma
import hashlib
import importlib.util
import io
import json
from pathlib import Path, PurePosixPath
import re
import resource
import sys
import tarfile
import time
import traceback

SIZE = 45548888
MD5 = 'f03adede5b7b81c11569608393aad5df'
URL = 'https://zenodo.org/api/records/5903179/files/Moldoveanu_2022_CyTOF_melanoma.tar.gz/content'
SAFE_FIELDS = {'sample', 'sample_id', 'cell_id', 'CellID', 'cell_type', 'celltype',
    'cellType', 'x', 'y', 'X', 'Y', 'Position', 'Ki67', 'Ki-67', 'CD3', 'CD4',
    'CD8', 'CD8a', 'FOXP3', 'FoxP3', 'PD1', 'PD-1', 'CTLA4', 'LAG3', 'TCF7',
    'CD45RO', 'CD45RA', 'CD45', 'Tumor', 'Cell', 'Sample', 'Patient', 'Type'}
SAFE_TABLE = re.compile(r'^(?:ST[1-9]_[A-Za-z_]+\.(?:txt|tsv)|(?:non-)?ICI_meanIntensity\.tsv\.gz)$')


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def load(path, expected, name):
    assert sha(path.read_bytes()) == expected
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def source_summary_only(value):
    if isinstance(value, dict):
        assert not set(value).intersection({'patient', 'patient_id', 'sample_id', 'cell_id',
            'barcode', 'donor', 'matrix', 'weights', 'patient_scores'})
        for item in value.values(): source_summary_only(item)
    elif isinstance(value, list):
        assert len(value) <= 1000
        if value and all(isinstance(x, (int, float)) for x in value): assert len(value) <= 5
        for item in value: source_summary_only(item)
    else:
        assert value is None or isinstance(value, (str, int, float, bool))


def qualify(raw, expected_size=SIZE, expected_md5=MD5):
    assert len(raw) == expected_size and hashlib.md5(raw).hexdigest() == expected_md5
    with gzip.GzipFile(fileobj=io.BytesIO(raw)) as compressed:
        unpacked = compressed.read(256*1024**2+1)
    assert len(unpacked) <= 256*1024**2
    archive = tarfile.open(fileobj=io.BytesIO(unpacked), mode='r:')
    entries = archive.getmembers()
    assert 0 < len(entries) <= 1000 and len({entry.name for entry in entries}) == len(entries)
    for entry in entries:
        path = PurePosixPath(entry.name)
        assert not path.is_absolute() and '..' not in path.parts
        assert entry.isfile() or entry.isdir()
        assert 0 <= entry.size <= 200*1024**2
    tables, code_hashes = [], {}
    for entry in entries:
        name = PurePosixPath(entry.name).name
        if not entry.isfile(): continue
        if name.endswith('.R') and re.fullmatch(r'[A-Za-z0-9_]+\.R', name):
            assert entry.size <= 1024**2
            code_hashes[name] = sha(archive.extractfile(entry).read())
        if not SAFE_TABLE.fullmatch(name): continue
        stream = archive.extractfile(entry)
        signature = stream.read(8); stream.seek(0)
        codec = 'plain_text_candidate'
        try:
            if signature.startswith(b'\x1f\x8b'):
                codec = 'gzip'; stream = gzip.GzipFile(fileobj=stream)
            elif signature.startswith(b'BZh'):
                codec = 'bzip2'; stream = bz2.BZ2File(stream)
            elif signature.startswith(b'\xfd7zXZ\x00'):
                codec = 'xz'; stream = lzma.LZMAFile(stream)
            header = stream.readline(65537)
            assert len(header) <= 65536 and header.endswith(b'\n')
            fields = next(csv.reader([header.decode('utf-8-sig')], delimiter='\t'))
            assert len(fields)>1 and set(fields)&SAFE_FIELDS
            tables.append(dict(status='TABLE_HEADER_ACCEPTED_NOT_ROWS',
                generic_table_name=name, member_bytes=entry.size, detected_codec=codec,
                field_count=len(fields), recognized_generic_fields=sorted(set(fields)&SAFE_FIELDS),
                header_SHA256=sha(header), header_only_not_rows=True))
        except (UnicodeError, OSError, EOFError, AssertionError, csv.Error) as error:
            tables.append(dict(status='TABLE_HEADER_NOT_ACCEPTED',generic_table_name=name,
                member_bytes=entry.size,detected_codec=codec,error_class=type(error).__name__))
        finally:
            stream.close()
    result = dict(status='COMPLETE_ARCHIVE_INTEGRITY_WITH_TABLE_HEADER_QUALIFICATION_NOT_CELL_ACCEPTANCE',
        archive_SHA256=sha(raw), archive_MD5_verified=True, full_gzip_CRC_checked=True,
        tar_directory_SHA256=sha(json.dumps([(x.name,x.size,x.type.decode('ascii')) for x in entries],separators=(',',':')).encode()),
        member_count=len(entries), uncompressed_bytes=len(unpacked), generic_tables=tables,
        generic_author_code_SHA256=code_hashes, no_author_code_executed=True,
        no_table_rows_or_outcomes_parsed=True, no_arrays_or_payloads_written=True)
    source_summary_only(result)
    return result, archive


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--scripts-root', type=Path, required=True)
    parser.add_argument('--pipeline-request',type=Path)
    args = parser.parse_args()
    args.output.mkdir(mode=0o700, parents=False, exist_ok=False)
    def save(name, value):
        source_summary_only(value)
        (args.output/name).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
    phase='archive_input'
    try:
        transport = load(args.scripts_root/'melanoma_IMC_verified_dns_transport_v2.py',
            '12ffe4ac9c1c37efefbe42e8153928473da3cee8997a2086118d46ba81d7c25c','public_strict_ranges')
        rf = transport.RangeFile(None, URL, SIZE, checkpoint=lambda x: save('private_transport_progress.json', x))
        chunks = []
        while rf.tell() < SIZE:
            chunks.append(rf.read(min(8*1024**2, SIZE-rf.tell())))
        raw = b''.join(chunks)
        chunks.clear();rf.close()
        assert sha(raw)=='d11cd7a5a1672ff002e40658f8d22bea0898733386d24197d817de9adca56112'
        result, archive = qualify(raw)
        del raw
        result.update(completed_UTC=datetime.now(timezone.utc).isoformat(), producer_SHA256=sha(Path(__file__).read_bytes()))
        assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024) <= 2*1024**3
        save('private_complete_archive_schema.json', result)
        scripts = args.output.parent/'scripts'
        scripts.mkdir(mode=0o700,exist_ok=True)
        if args.pipeline_request is not None:
            phase='own_pipeline'
            request=args.pipeline_request
            assert request.resolve().parent==args.output.parent.resolve() and not request.is_symlink()
            protocol=json.loads(request.read_text())
            assert set(protocol)=={'step','code_file','code_SHA256','parameters'}
            step=protocol['step'];assert re.fullmatch('[A-Za-z0-9_]{1,60}',step)
            filename=protocol['code_file'];assert Path(filename).name==filename and filename.endswith('.py')
            path=scripts/filename;assert path.resolve().parent==scripts.resolve() and not path.is_symlink()
            own=load(path,protocol['code_SHA256'],'fixed_'+step)
            save('private_RAM_session_status.json',dict(status='RUNNING_FROZEN_SINGLE_PIPELINE',step=step,
                request_SHA256=sha(request.read_bytes()),observed_UTC=datetime.now(timezone.utc).isoformat()))
            answer=own.analyze(dict(archive=archive,qualification=result),protocol['parameters'])
            source_summary_only(answer)
            assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)<=2*1024**3
            save('private_single_pipeline_source_summary.json',answer)
            archive.close()
            save('private_RAM_session_status.json',dict(status='RAM_SESSION_CLOSED_AFTER_SINGLE_PIPELINE',
                scientific_status=answer['status'],observed_UTC=datetime.now(timezone.utc).isoformat()))
            return
        started=time.monotonic();done=set()
        save('private_RAM_session_status.json',dict(status='READY_WITH_ARCHIVE_RAM_FOR_FROZEN_OWN_QUALIFICATION', maximum_hold_seconds=45*60, observed_UTC=datetime.now(timezone.utc).isoformat()))
        while time.monotonic()-started < 45*60 and len(done)<3:
            for request in sorted(scripts.glob('frozen_step_request_*.json')):
                if request.name in done: continue
                try:
                    protocol=json.loads(request.read_text())
                    assert set(protocol)=={'step','code_file','code_SHA256','parameters'}
                    step=protocol['step'];assert re.fullmatch('[A-Za-z0-9_]{1,60}',step)
                    filename=protocol['code_file'];assert Path(filename).name==filename and filename.endswith('.py')
                    path=scripts/filename;assert path.resolve().parent==scripts.resolve() and not path.is_symlink()
                    own=load(path,protocol['code_SHA256'],'frozen_'+step)
                    save('private_RAM_session_status.json',dict(status='RUNNING_FROZEN_OWN_STEP',step=step,request_SHA256=sha(request.read_bytes()),observed_UTC=datetime.now(timezone.utc).isoformat()))
                    answer=own.analyze(dict(archive=archive,qualification=result),protocol['parameters'])
                    source_summary_only(answer)
                    assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)<=2*1024**3
                    save('private_frozen_step_'+step+'_source_summary.json',answer)
                except Exception as error:
                    save('private_step_'+request.stem+'_not_accepted.json',dict(status='FROZEN_STEP_NOT_ACCEPTED',error_class=type(error).__name__,frames=[dict(file=Path(x.filename).name,line=x.lineno,function=x.name) for x in traceback.extract_tb(error.__traceback__)]))
                done.add(request.name)
                save('private_RAM_session_status.json',dict(status='READY_WITH_ARCHIVE_RAM_FOR_FROZEN_OWN_QUALIFICATION',completed_steps=sorted(done),observed_UTC=datetime.now(timezone.utc).isoformat()))
            time.sleep(5)
        archive.close()
        save('private_RAM_session_status.json',dict(status='RAM_SESSION_CLOSED',completed_steps=sorted(done),observed_UTC=datetime.now(timezone.utc).isoformat()))
    except Exception as error:
        save('private_failed_terminal.json',dict(status='ARCHIVE_SCOPE_NOT_ACCEPTED' if phase=='archive_input' else 'OWN_PIPELINE_NOT_ACCEPTED_ARCHIVE_INPUT_REMAINS_SEPARATE',phase=phase,error_class=type(error).__name__,frames=[dict(file=Path(x.filename).name,line=x.lineno,function=x.name) for x in traceback.extract_tb(error.__traceback__)],observed_UTC=datetime.now(timezone.utc).isoformat(),no_individual_values_in_error=True))
        if 'archive' in locals():archive.close()
        save('private_RAM_session_status.json',dict(status='RAM_SESSION_CLOSED_AFTER_FAILURE',phase=phase,observed_UTC=datetime.now(timezone.utc).isoformat()))
        raise SystemExit(1)


if __name__ == '__main__': main()
