"""Bounded processed-TCR CSV header, without data rows or sequences."""
import csv
import hashlib
import io
import json
import os
import urllib.request
import zlib
from datetime import datetime, timezone
from pathlib import Path
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource

URL = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE243nnn/GSE243013/suppl/GSE243013_T_with_TCR_annotation.csv.gz'


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root, 64 * 1024)
    out = root / 'public_GSE243013_TCR_header_qualification_20261019_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    cap = 256 * 1024
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), public_URL=URL, compressed_payload_cap_bytes=cap, scope='Only first decompressed CSV header, no data rows, full file, sequences, patient or cell values', resource_guard=guard)
    (out / 'plan_before_public_TCR_header_request.json').write_text(json.dumps(plan, indent=2) + '\n')
    payload = bytearray()
    decoded = bytearray()
    unzip = zlib.decompressobj(16 + zlib.MAX_WBITS)
    with urllib.request.urlopen(urllib.request.Request(URL, headers={'Range': 'bytes=0-262143'}), timeout=30) as response:
        assert response.status == 206 and response.headers.get('Content-Range') == 'bytes 0-262143/15940803'
        for chunk in iter(lambda: response.read(8192), b''):
            assert len(payload) + len(chunk) <= cap
            payload.extend(chunk)
            decoded.extend(unzip.decompress(chunk, 1024 * 1024))
            if b'\n' in decoded:
                break
        assert b'\n' in decoded
        header = bytes(decoded).split(b'\n', 1)[0].decode('utf-8-sig')
        fields = next(csv.reader(io.StringIO(header)))
        assert fields == ['sampleID', 'cellID', 'sub_cell_type', 'TRA_v_gene', 'TRA_j_gene', 'TRA_c_gene', 'TRA_cdr3', 'TRB_v_gene', 'TRB_j_gene', 'TRB_c_gene', 'TRB_cdr3', 'clonotype', 'expansion', 'clonotype_number', 'T_new_name']
        proof = dict(public_URL=URL, observed_UTC=datetime.now(timezone.utc).isoformat(), HTTP_status=response.status, range_response=response.headers['Content-Range'], requested_compressed_byte_cap=cap, compressed_bytes_read=len(payload), partial_compressed_sha256=hashlib.sha256(payload).hexdigest(), CSV_header_fields=fields, explicit_TRG_TRD_columns_present=False, alpha_beta_named_columns_present=True, full_table_content_or_productivity_not_inspected=True, complete_file_downloaded=False, whole_file_SHA_or_CRC_verified=False, data_rows_or_sequences_exported=False)
    (out / 'private_header_qualification.json').write_text(json.dumps(proof, indent=2) + '\n')
    receipt = dict(status='QUALIFIED_PUBLIC_PROCESSED_TCR_HEADER_ONLY_EXPLICIT_ALPHA_BETA_NOT_GAMMA_CONTENT_VERIFIED', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), header_qualification_sha256=sha(out / 'private_header_qualification.json'), resource_guard=guard, scientific_identity_not_accepted=True)
    (out / 'private_terminal.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
