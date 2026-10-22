"""Full synthetic HTTP/RDS/native analysis/audit/exit path; no real network or records."""
from pathlib import Path
from datetime import datetime, timezone
from email.message import Message
from contextlib import redirect_stdout
from types import SimpleNamespace
import base64
import hashlib
import io
import json
import os
import subprocess
import sys
import tempfile


def main():
    R = Path(sys.executable).resolve().parents[2]/'scrna_r/bin/Rscript'
    pc = TEST_PHENOTYPE_CONFIG
    maker = DERIVE_SOURCE + r'''
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(jsonlite))
pc<-fromJSON(Sys.getenv("PANCANCER_LN_TEST_PC"))
meta<-do.call(rbind,lapply(seq_len(7),function(p)do.call(rbind,lapply(c("N1","N2"),function(nd)data.frame(patient=paste0("FAKE_P_",p),node=nd,image=paste0("FAKE_IMAGE_",p,"_",nd),fine=rep(c("D4","D8","B","M"),c(80,80,40,40)),j=c(1:80,1:80,1:40,1:40),area=50+(1:240)%%11,stringsAsFactors=FALSE)))))
a<-matrix(1,length(pc$physical_panel),nrow(meta),dimnames=list(pc$physical_panel,NULL))
a["CD4",meta$fine=="D4"]<-10;a["CD8a",meta$fine=="D8"]<-10;a["CD19",meta$fine=="B"]<-10;a[c("CD14","CD68","CD163"),meta$fine=="M"]<-10;a["CD45RO",meta$fine%in%c("D4","D8")]<-10
for(m in c("FoxP3","CTLA4"))a[m,meta$fine=="D4"]<-1+meta$j[meta$fine=="D4"]/5
for(m in c("PD1","TIM3"))a[m,meta$fine=="D8"]<-1+meta$j[meta$fine=="D8"]/5
for(m in c("CCR7","CD127"))a[m,meta$fine=="D8"]<-1+(81-meta$j[meta$fine=="D8"])/5
a["MKI67",]<-1+(meta$j*13)%%80/5
xy<-cbind(meta$j*3+ifelse(meta$fine=="D8",.1,0),ifelse(meta$fine=="D4",0,ifelse(meta$fine=="D8",.2,ifelse(meta$fine=="B",30,60))))
x<-SingleCellExperiment(assays=list(counts=a,asin=asinh(a)),colData=S4Vectors::DataFrame(meta))
SingleCellExperiment::int_colData(x)$spatialCoords<-xy
ph<-derive_LN_hierarchical_protein_states_RAM(x,pc,"counts")
counts<-setNames(lapply(ph$summary$state_records,function(r)r$cell_rows),vapply(ph$summary$state_records,function(r)r$state,character(1)))
cat(toJSON(list(raw=base64_enc(serialize(x,NULL)),expected_state_counts=counts,expected_parent_groups=ph$summary$resolved_parent_groups,cells=ncol(x)),auto_unbox=TRUE,null="null"))
'''
    env = dict(os.environ, PANCANCER_LN_TEST_PC=json.dumps(pc))
    fixture_root=Path(tempfile.mkdtemp(prefix='pancancer_LN_SPEC_SYNTHETIC_fixture_',dir='/dev/shm'))
    fixture_code=fixture_root/'make_synthetic_fixture.R';fixture_code.write_text(maker)
    made = subprocess.run([str(R), '--vanilla', str(fixture_code)], capture_output=True, env=env, timeout=120)
    if made.returncode:
        raise RuntimeError('SYNTHETIC_SPECIFICITY_FIXTURE_FAILED')
    fixture = json.loads(made.stdout)
    payload = base64.b64decode(fixture.pop('raw'))
    record_id, key = 'SYNTHETIC_NO_NETWORK', 'synthetic.rds'
    metadata_url = 'https://zenodo.org/api/records/'+record_id
    content_url = 'https://zenodo.org/records/'+record_id+'/files/'+key+'?download=1'
    metadata = json.dumps(dict(id=record_id, metadata=dict(access_right='open'),
        files=[dict(key=key,size=len(payload),checksum='md5:'+hashlib.md5(payload).hexdigest())])).encode()
    requested, controls = [], []

    class Response(io.BytesIO):
        def __init__(self, value, url, span=None, declared_length=None):
            super().__init__(value);self.url=url;self.status=206 if span else 200
            self.headers=Message();self.headers['content-length']=str(len(value) if declared_length is None else declared_length);self.headers['content-type']='application/octet-stream'
            if span:self.headers['content-range']=span

    class FakeOpener:
        def open(self, request, timeout):
            span=request.get_header('Range');requested.append((request.full_url,span))
            if request.full_url==metadata_url:
                assert span is None
                return Response(metadata,metadata_url)
            assert request.full_url==content_url and span.startswith('bytes=')
            start,end=map(int,span[6:].split('-'))
            body=payload[start:end+1]
            if len(body)>32:body=body[:-31]
            return Response(body,content_url,f'bytes {start}-{end}/{len(payload)}',declared_length=end-start+1)

    def control_dir(**kwargs):
        root=Path(tempfile.mkdtemp(prefix='pancancer_LN_SPEC_SYNTHETIC_',dir='/dev/shm'));controls.append(root);return str(root)

    spatial=dict(image_field='image',patient_field='patient',node_field='node',k=1,
        minimum_parent=20,minimum_class=2,minimum_mixed_blocks=2,minimum_movable_class=2,
        grids=[2,4,8],primary_grid=4,permutations=19,bootstrap_replicates=2000,minimum_cases=5,seed=20261023)
    analysis=dict(phenotype_config=pc,spatial_config=spatial,expected_cells=fixture['cells'],
        expected_state_counts=fixture['expected_state_counts'],expected_parent_groups=fixture['expected_parent_groups'],
        family_contract_SHA256=hashlib.sha256(b'SYNTHETIC_SECONDARY_FAMILY').hexdigest())
    plan=dict(record=record_id,expected_file_bytes=len(payload),expected_file_MD5=hashlib.md5(payload).hexdigest(),
        expected_file_SHA256=hashlib.sha256(payload).hexdigest(),file_key_SHA256=hashlib.sha256(key.encode()).hexdigest(),
        R_code_SHA256=hashlib.sha256(R_SOURCE.encode()).hexdigest(),pipe_code_SHA256=hashlib.sha256(PIPE_SOURCE.encode()).hexdigest(),
        range_code_SHA256=hashlib.sha256(RANGE_SOURCE.encode()).hexdigest(),range_core_code_SHA256=hashlib.sha256(RANGE_CORE_SOURCE.encode()).hexdigest(),
        parent_address_cap_bytes=4*1024**3,inherited_address_hard_cap_bytes=16*1024**3,R_address_cap_bytes=12*1024**3,
        minimum_available_memory_bytes=24*1024**3,maximum_download_seconds=180,maximum_R_seconds=240,maximum_expanded_bytes=64*1024**2,
        previous_source_pid=0,previous_source_start_ticks=0,previous_R_pid=0,previous_R_start_ticks=0,
        range_chunk_bytes=4096,range_workers=4,range_request_timeout_seconds=10,maximum_range_seconds=30,
        range_maximum_attempts=4,range_retry_backoff_seconds=0,range_maximum_attempt_seconds=90,analysis_config=analysis)
    ns=dict(__name__='synthetic_secondary_reader',R_SOURCE=R_SOURCE,PIPE_SOURCE=PIPE_SOURCE,
        RANGE_SOURCE=RANGE_SOURCE,RANGE_CORE_SOURCE=RANGE_CORE_SOURCE,RANGE_OPENER_FACTORY=lambda:FakeOpener())
    exec(READER_SOURCE,ns);ns['build_opener']=lambda *a,**kw:FakeOpener();ns['tempfile']=SimpleNamespace(mkdtemp=control_dir)
    original=sys.stdin;captured=io.StringIO()
    try:
        sys.stdin=io.StringIO(json.dumps(plan))
        with redirect_stdout(captured):ns['main']()
    finally:sys.stdin=original
    frames=[json.loads(line) for line in captured.getvalue().splitlines()];final=frames[-1]
    if final['kind']!='COMPLETE_PUBLIC_PROTEIN_SECONDARY_SPECIFICITY_RAM_RELEASED':
        print(json.dumps(dict(status='FULL_SYNTHETIC_NATIVE_WORKFLOW_NOT_ACCEPTED',safe_frames=frames[-2:])));raise SystemExit(2)
    assert final['source_R_exit_code']==0 and final['source_payload_RAM_released'] and final['no_dynamic_followup_nodes']
    assert final['summary']['audit']['status']=='PASS_FULL_SECONDARY_SPECIFICITY_NUMERIC_AUDIT'
    assert len(final['summary']['production']['primary_rows'])==15
    spans={f'bytes={s}-{min(s+4095,len(payload)-1)}' for s in range(0,len(payload),4096)}
    assert requested[0]==(metadata_url,None) and all(u==content_url for u,v in requested[1:])
    actual={span for u,span in requested[1:]}
    expected=set(spans)
    for span in spans:
        start,end=map(int,span[6:].split('-'))
        if end-start+1>32:expected.add(f'bytes={end-30}-{end}')
    assert actual==expected and len(requested)==1+len(expected)
    assert final['range_assembly_proof']['logical_ranges_requiring_continuation']>0
    allowed={'private_verified_public_input_identity.json','native_reader.R','private_specificity_phenotype_ready.json',
        'private_specificity_production_started.json','private_specificity_production_complete.json','private_specificity_audit_complete.json'}
    assert len(controls)==1 and {f.name for f in controls[0].iterdir()}<=allowed
    assert 'FAKE_P_' not in captured.getvalue() and 'FAKE_IMAGE_' not in captured.getvalue()
    print(json.dumps(dict(status='PASS_FULL_SYNTHETIC_SUFFIX_HTTP_NATIVE_SPECIFICITY_AUDIT_AND_EXIT',
        completed_UTC=datetime.now(timezone.utc).isoformat(),actual_source_R_exit_code=0,
        native_entrypoint_executed=True,frozen_phenotype_and_full_secondary_family_replayed=True,
        independent_numeric_audit=final['summary']['audit'],exact_range_coverage_MD5_and_SHA256_checked=True,
        only_generic_code_and_aggregate_control_files_created=True,no_real_network_or_human_records_read=True)))


if __name__=='__main__':main()
