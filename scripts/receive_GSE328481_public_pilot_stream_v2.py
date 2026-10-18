"""Receive one effect-blind public processed pilot through a byte-preserving relay.

No local matrix persistence required; selection and quota checked on server.
"""
import hashlib
import json
import sys
from datetime import datetime,timezone
from pathlib import Path

from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC_native_MFA_interpretation_inputs_v1 import sha


def main():
    root=Path(sys.argv[1]);budget=resource(root,128*1024**2)
    src=root/'GSE328481_public_processed_spatial_metadata_20261018_v2'
    t=json.loads((src/'private_terminal.json').read_text())
    for name,item in t['files'].items():assert sha(src/name)==item['sha256']
    files=json.loads((src/'public_all_processed_file_inventory.json').read_text())
    assert len(files)==11 and all(x['status']==200 and int(x['content_length'])>0 and x['filename'].endswith('.h5ad') for x in files)
    chosen=min(files,key=lambda x:(int(x['content_length']),x['url']));size=int(chosen['content_length']);assert size+8*1024**2<128*1024**2
    out=root/'GSE328481_public_pilot_stream_transfer_20261018_v2';assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),budget=budget,
        metadata_terminal_sha256=sha(src/'private_terminal.json'),chosen=chosen,
        reason='Whole-file GET remains slow despite local byte relay, while a4MiB Range request completed in2.90s. Use four bounded exact contiguous ranges with Last-Modified/ETag, whole SHA and size checks; no new proxy tunnel or port.',
        scope='Same minimum-size public processed pilot selected before expression. Local client keeps only bounded network buffers and no local expression file. Server stores source once in this new directory, original partial retained.',
        assertions='Exact official HEAD byte count, HDF5 signature, streaming SHA and independent full-file SHA readback. Local client SHA must match server. Not official expected hash or biological qualification.')
    (out/'plan_before_public_byte_relay.json').write_text(json.dumps(plan,indent=2)+'\n')
    target=out/chosen['filename'];n=0;digest=hashlib.sha256();last=-1
    with target.open('xb') as f:
        while True:
            b=sys.stdin.buffer.read(1024**2)
            if not b:break
            n+=len(b);assert n<=size;f.write(b);digest.update(b)
            if n//(8*1024**2)!=last:
                last=n//(8*1024**2);(out/'private_transfer_progress.json').write_text(json.dumps(dict(bytes_received=n,expected_bytes=size,observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n')
    assert n==size and target.stat().st_size==size and sha(target)==digest.hexdigest()
    with target.open('rb') as f:assert f.read(8)==b'\x89HDF\r\n\x1a\n'
    report=dict(status='COMPLETE_PUBLIC_PROCESSED_PILOT_BYTE_TRANSFER_ONLY',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),chosen=chosen,bytes=n,download_sha256=digest.hexdigest(),
        files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()},expression_or_clonotype_values_interpreted=False)
    (out/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(status=report['status'],completed_UTC=report['completed_UTC'],download_sha256=report['download_sha256'],terminal_sha256=sha(out/'private_terminal.json'))))

if __name__=='__main__':main()
