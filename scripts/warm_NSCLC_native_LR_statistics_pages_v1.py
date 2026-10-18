"""Sequential read-only cache warming with full original-array SHA verification."""
import hashlib
import json
import os
from datetime import datetime,timezone
from pathlib import Path
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,1024**2);out=root/'NSCLC_native_LR_statistics_page_warm_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    native=root/'NSCLC_full_native_gene_units_20261016_v1';t=json.loads((native/'private_terminal.json').read_text());assert t['status']=='ALL_NATIVE_NSCLC_FULL_GENE_UNIT_SUMMARIES_COMPLETE_INDEPENDENT_RAW_AUDIT_PENDING'
    memory={s.split(':')[0]:s.split(':')[1].strip() for s in Path('/proc/meminfo').read_text().splitlines()};assert int(memory['MemAvailable'].split()[0])*1024>16*1024**3
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,scope='Read existing sum_raw_counts/positive_cell_counts/mean_logCP10K arrays sequentially,64MiBbuffer only; no matrix copy/write, no source download, no production/audit restart. Original sparse-column mmap comparison has observed diskpage waits despite abundant available memory. SHA checks bind each complete original array while OS page cache warms.',source_terminal_sha256=sha(native/'private_terminal.json'))
    (out/'plan_before_readonly_page_warm.json').write_text(json.dumps(plan,indent=2)+'\n');records=[]
    for name in ['sum_raw_counts.npy','positive_cell_counts.npy','mean_logCP10K.npy']:
        p=native/name;expected=t['outputs'][name];assert p.stat().st_size==expected['bytes'];h=hashlib.sha256();read=0
        with p.open('rb') as f:
            while True:
                b=f.read(64*1024**2)
                if not b:break
                h.update(b);read+=len(b)
                (out/'private_progress.json').write_text(json.dumps(dict(phase='SEQUENTIAL_EXISTING_ARRAY_READ',array=name,bytes_read=read,array_bytes=expected['bytes'],observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n')
        assert read==expected['bytes'] and h.hexdigest()==expected['sha256'];records.append(dict(array=name,bytes=read,sha256=h.hexdigest()))
    receipt=dict(status='COMPLETE_READONLY_SEQUENTIAL_PAGE_WARM_AND_ALL3_ORIGINAL_ARRAY_SHA_PARITY',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,arrays=records,no_existing_artifact_changed=True,no_inference_executed=True)
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'status':receipt['status'],'terminal_sha256':sha(out/'private_terminal.json')}))


if __name__=='__main__':main()
