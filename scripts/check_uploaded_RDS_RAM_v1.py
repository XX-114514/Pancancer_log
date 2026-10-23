"""Meaningful artificial-only input corruption and concurrent-upload checks."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import tempfile
from read_uploaded_RDS_RAM_v1 import UploadedInputError, read_verified_upload
from pipe_native_RDS_RAM_v1 import pipe_native_RDS


def check(Rscript=None):
    root = Path(tempfile.mkdtemp(prefix='ARTIFICIAL_upload_input_checks_'))
    payload = bytes(range(256))*16
    md5, sha = hashlib.md5(payload).hexdigest(), hashlib.sha256(payload).hexdigest()
    path = root/'ARTIFICIAL_VALID.bin'; path.write_bytes(payload)
    defaults = dict(expected_bytes=len(payload), expected_MD5=md5, expected_SHA256=sha,
                    maximum_input_bytes=8192, maximum_seconds=30, chunk_bytes=512)
    buffer, m, s, proof = read_verified_upload(path, **defaults)
    assert buffer.getvalue() == payload and m.hexdigest() == md5 and s.hexdigest() == sha
    buffer.close()
    assert proof['bytes_read'] == len(payload) and proof['no_network_requests'] is True
    checks = ['COMPLETE_HASHED_READ_ONLY_BUFFER']

    def reject(p, reason, **changes):
        config = dict(defaults); config.update(changes)
        try:
            read_verified_upload(p, **config)
        except UploadedInputError as error:
            assert error.code_name == reason, (error.code_name, reason)
            checks.append(reason)
        else:
            raise AssertionError('ARTIFICIAL_BAD_INPUT_ACCEPTED')

    short = root/'ARTIFICIAL_PARTIAL.bin';short.write_bytes(payload[:100])
    reject(short, 'UPLOAD_LENGTH_NOT_COMPLETE')
    reject(path, 'WHOLE_FILE_IDENTITY_MISMATCH', expected_MD5='0'*32)
    reject(path, 'WHOLE_FILE_IDENTITY_MISMATCH', expected_SHA256='0'*64)
    reject(path, 'INVALID_INPUT_RESOURCE_BUDGET', maximum_input_bytes=16)
    reject(path, 'INPUT_READ_TIME_BUDGET', maximum_seconds=1e-12)
    link = root/'ARTIFICIAL_SYMLINK';link.symlink_to(path)
    reject(link, 'INPUT_NOT_OWNED_REGULAR_FILE')
    reject(root, 'INPUT_NOT_OWNED_REGULAR_FILE')
    grow = root/'ARTIFICIAL_GROW.bin';grow.write_bytes(payload)
    def grow_once(count):
        if count == 512:
            with grow.open('ab') as stream:stream.write(b'X')
    reject(grow, 'INPUT_GREW_DURING_READ', progress=grow_once)
    changed = root/'ARTIFICIAL_CHANGE.bin';changed.write_bytes(payload)
    def modify_once(count):
        if count == 512:
            with changed.open('r+b') as stream:stream.write(b'Y')
    reject(changed, 'INPUT_CHANGED_DURING_READ', progress=modify_once)
    replaced = root/'ARTIFICIAL_REPLACE.bin';replaced.write_bytes(payload)
    def replace_once(count):
        if count == 512:
            replaced.rename(root/'ARTIFICIAL_REPLACED_ORIGINAL_RETAINED.bin')
            replaced.write_bytes(payload)
    reject(replaced, 'INPUT_CHANGED_DURING_READ', progress=replace_once)
    if Rscript:
        generated = subprocess.run([Rscript,'--vanilla','-e',
            'z<-file("/dev/stdout","wb");saveRDS(list(ARTIFICIAL_ONLY=TRUE,x=1:32),z,compress="gzip");close(z)'],
            stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=30,check=True).stdout
        assert generated
        rds = root/'ARTIFICIAL_ONLY.rds';rds.write_bytes(generated)
        buffer, _, _, _ = read_verified_upload(rds, expected_bytes=len(generated),
            expected_MD5=hashlib.md5(generated).hexdigest(),expected_SHA256=hashlib.sha256(generated).hexdigest(),
            maximum_input_bytes=8192,maximum_seconds=30)
        child, receipt = pipe_native_RDS(buffer,[Rscript,'--vanilla','-e',
            'z<-file("stdin","rb");x<-readRDS(z);close(z);stopifnot(isTRUE(x$ARTIFICIAL_ONLY),identical(x$x,1:32));cat("{\\"status\\":\\"ARTIFICIAL_RDS_PIPE_PASS\\"}")'],
            maximum_expanded_bytes=8192,timeout=30)
        buffer.close()
        assert child.returncode == 0 and json.loads(child.stdout)['status'] == 'ARTIFICIAL_RDS_PIPE_PASS'
        assert receipt['no_native_or_expanded_payload_files'] is True
        checks.append('READ_ONLY_EXISTING_ARTIFICIAL_RDS_TO_RAM_TO_ORIGINAL_R_PIPE')
    result = dict(status='PASS_ARTIFICIAL_UPLOAD_INPUT_GATES',checks=checks,
                  artificial_fixture_directory_retained=True,no_real_native_or_clinical_file_read=True)
    print(json.dumps(result))
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--Rscript')
    check(parser.parse_args().Rscript)
