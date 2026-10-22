"""Source-runtime synthetic native codec check. No real data or file writes."""
from pathlib import Path
from datetime import datetime, timezone
import base64
import bz2
import gzip
import hashlib
import io
import json
import lzma
import subprocess
import sys


def main():
    R = Path(sys.executable).resolve().parents[2]/'scrna_r/bin/Rscript'
    create = '''
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(jsonlite))
invisible(setClass("CodecSyntheticLeaf",contains="SingleCellExperiment",slots=c(spatialCoords="matrix")))
m<-matrix(seq_len(60),6,10,dimnames=list(c("CD3","CD4","CD8a","FOXP3","PD1","Ki67"),NULL))
s<-SingleCellExperiment(assays=list(counts=m),colData=S4Vectors::DataFrame(patient_id=rep(c("FAKE_A","FAKE_B"),each=5)))
x<-new("CodecSyntheticLeaf",s,spatialCoords=matrix(seq_len(20),10,2))
cat(toJSON(list(raw=base64_enc(serialize(x,NULL))),auto_unbox=TRUE))
'''
    made = subprocess.run([str(R), '--vanilla', '-e', create], capture_output=True, timeout=90)
    if made.returncode:
        print(json.dumps(dict(status='SYNTHETIC_CODEC_GENERATOR_FAILED', diagnostic=made.stderr.decode()[-2000:])))
        raise SystemExit(2)
    raw = base64.b64decode(json.loads(made.stdout)['raw'])
    reader = R_SOURCE + chr(10) + '''
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(jsonlite))
x<-readRDS(file("stdin","rb"))
cat(toJSON(qualify_structure(x),auto_unbox=TRUE,null="null",digits=16))
'''
    receipts, baseline = [], None
    for name, value in [('none', raw), ('gzip', gzip.compress(raw)), ('bzip2', bz2.compress(raw)), ('xz', lzma.compress(raw))]:
        ns = {}
        exec(PIPE_SOURCE, ns)
        result, technical = ns['pipe_native_RDS'](io.BytesIO(value), [str(R), '--vanilla', '-e', reader], maximum_expanded_bytes=64*1024**2, timeout=90)
        if result.returncode:
            print(json.dumps(dict(status='SYNTHETIC_NATIVE_CODEC_FAILED', codec=name, source_R_exit_code=result.returncode)))
            raise SystemExit(2)
        q = json.loads(result.stdout)
        assert q['source_cells'] == 10 and 'FAKE_A' not in result.stdout.decode()
        if baseline is None:
            baseline = q
        else:
            assert q == baseline
        assert technical['codec'] == name and technical['expanded_stream_bytes'] == len(raw)
        receipts.append(dict(codec=name, identical_qualified_structure=True, no_individual_values_emitted=True))
    print(json.dumps(dict(status='PASS_FOUR_NATIVE_PIPE_CODECS_IDENTICAL_SOURCE_ONLY_STRUCTURE',
        observed_UTC=datetime.now(timezone.utc).isoformat(), R_schema_SHA256=hashlib.sha256(R_SOURCE.encode()).hexdigest(), pipe_code_SHA256=hashlib.sha256(PIPE_SOURCE.encode()).hexdigest(),
        formats=receipts, no_human_data_or_source_requests=True)))


if __name__ == '__main__':
    main()
