"""Bounded, task-local official scGate/UCell runtime; synthetic cells only.

Fetch immutable public software/model archives, without any expression data.
No shared package library is modified and no dependency update is performed.
"""
import hashlib
import io
import json
import os
import subprocess
import sys
import tarfile
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource

SOURCES = [
    ('UCell', 'f778e091f19e28d3d05bde134663fd6b2f938267'),
    ('scGate', '14dcede8f60298fddd2c47aa139e3edffbb31386'),
    ('scGate_models', '41a45cd3f8bb5f5a7daf21ec276f6a726f6ee0d4'),
]

CHECK = r'''
args<-commandArgs(TRUE)
.libPaths(c(args[1],file.path(args[2],c('library_v7','library_v5','library_v2','library')),.libPaths()))
for(p in c('Seurat','dplyr','patchwork','ggridges','colorspace','reshape2','ggplot2','BiocParallel','data.table','Matrix','BiocNeighbors','SingleCellExperiment','SummarizedExperiment'))stopifnot(requireNamespace(p,quietly=TRUE))
cat('EXISTING_IMPORTS_ACCEPTED\n')
'''

SMOKE = r'''
args<-commandArgs(TRUE)
.libPaths(c(args[1],file.path(args[2],c('library_v7','library_v5','library_v2','library')),.libPaths()))
suppressPackageStartupMessages(library(scGate))
stopifnot(as.character(packageVersion('scGate'))=='1.7.3',as.character(packageVersion('UCell'))=='2.17.1')
master<-file.path(args[3],'human','HiTME','master_table.tsv')
model<-scGate::load_scGate_model(file.path(args[3],'human','HiTME','gdT_scGate_Model.tsv'),master.table=master)
stopifnot(all(nchar(model$signature)>0))
tab<-read.delim(master,check.names=FALSE)
g<-unique(sub('[+-]$','',unlist(strsplit(tab$signature,';',fixed=TRUE))))
g<-unique(c(g,paste0('SYNTHETIC_BACKGROUND_',seq_len(500))))
set.seed(20261019)
x<-matrix(rpois(length(g)*240,.05),nrow=length(g),dimnames=list(g,paste0('synthetic_',seq_len(240))))
sets<-list(c('PTPRC','LAPTM5','SRGN','CXCR4','CD52','LCK','CD3D','CD3E','CD3G','CD2','TRDC','TRGC1','TRGC2','TRDV1','TRDV2'),c('PTPRC','LCK','KLRD1','NKG7','NCR1','FCGR3A'),c('PTPRC','LCK','CD3D','CD3E','CD3G','CD2','TRAC','TRBC1','TRBC2','CD8A','CD8B'),c('PTPRC','SPI1','LYZ','CSF1R','MSR1','MAFB','CD300E','ITGAX','CD68'))
for(i in seq_along(sets))x[intersect(sets[[i]],rownames(x)),((i-1)*60+1):(i*60)]<-x[intersect(sets[[i]],rownames(x)),((i-1)*60+1):(i*60)]+20
obj<-Seurat::CreateSeuratObject(counts=Matrix::Matrix(x,sparse=TRUE),min.cells=0,min.features=0)
obj<-Seurat::NormalizeData(obj,verbose=FALSE)
obj<-scGate::scGate(obj,model=model,seed=123,pos.thr=.2,neg.thr=.2,maxRank=1500,nfeatures=2000,pca.dim=30,k.param=30,ncores=1,progressbar=FALSE,verbose=FALSE)
stopifnot(ncol(obj)==240,all(obj$is.pure%in%c('Pure','Impure')))
cat('OFFICIAL_SYNTHETIC_RUNTIME_ACCEPTED\n')
cat('scGate=',as.character(packageVersion('scGate')),';UCell=',as.character(packageVersion('UCell')),'\n',sep='')
'''


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT']); guard = resource(root, 12 * 1024**2)
    decision = root / 'NSCLC_gamma_delta_RNA_candidate_review_20261019_v1/private_candidate_extension_decision.json'
    assert sha(decision) == '14c2ae9566e5d8d87b62ad9cee126360e60d21aafc9a241f59a02c171ad382d4'
    assert json.loads(decision.read_text())['continue_fixed_candidate_only_extension'] is False
    out = root / 'gamma_scGate_official_runtime_qualification_20261019_v1'; assert not out.exists(); out.mkdir(mode=0o700)
    lib = out / 'library'; lib.mkdir(mode=0o700)
    runtime = root.parents[2] / 'Project_v5/runs/20261014_DIALOGUE_method_pilot_v1'
    R = Path(sys.executable).resolve().parents[2] / 'scrna_r/bin/R'
    Rscript = R.with_name('Rscript'); assert R.is_file() and Rscript.is_file()
    save = lambda name, value: (out / name).write_text(json.dumps(value, indent=2) + '\n')
    save('plan_before_official_runtime_preparation.json', dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), resource_guard=guard, source_commits=SOURCES, scope='Software and public predefined models only. No expression download, real-cell classification or clinical analysis. Task-local package library; do not install/update shared dependencies.', prerequisite_scientific_review_sha256=sha(decision), archive_total_cap_bytes=4 * 1024**2, extracted_total_cap_bytes=4 * 1024**2, full_runtime_cap_bytes=12 * 1024**2, official_defaults='pos/neg0.2,maxRank1500,k30,nfeatures2000,pca30,seed123. Synthetic run is compatibility only, not biological classifier validation.'))
    check = out / 'check_existing_imports.R'; check.write_text(CHECK)
    env = os.environ.copy(); env['R_LIBS_USER'] = ':'.join([str(lib)] + [str(runtime / v) for v in ['library_v7', 'library_v5', 'library_v2', 'library']])
    log = out / 'private_software_preparation.log'
    with log.open('xb') as f:
        p = subprocess.run([str(Rscript), '--vanilla', str(check), str(lib), str(runtime)], stdout=f, stderr=f, env=env, timeout=90)
        assert p.returncode == 0, 'EXISTING_REQUIRED_DEPENDENCY_UNAVAILABLE_NO_INSTALL'
    total_download = 0; extracted = 0; records = []; roots = {}
    for repo, commit in SOURCES:
        url = 'https://codeload.github.com/carmonalab/' + repo + '/tar.gz/' + commit
        with urllib.request.urlopen(url, timeout=45) as response: blob = response.read(4 * 1024**2 + 1)
        total_download += len(blob); assert total_download < 4 * 1024**2
        archive = out / (repo + '_immutable_public_source.tar.gz'); archive.write_bytes(blob)
        with tarfile.open(fileobj=io.BytesIO(blob), mode='r:gz') as t:
            members = t.getmembers(); extracted += sum(m.size for m in members)
            assert extracted < 4 * 1024**2
            assert all(not Path(m.name).is_absolute() and '..' not in Path(m.name).parts and not m.issym() and not m.islnk() for m in members)
            top = {Path(m.name).parts[0] for m in members}; assert len(top) == 1
            t.extractall(out)
        source = out / next(iter(top)); roots[repo] = source
        records.append(dict(repository='carmonalab/' + repo, commit=commit, archive_sha256=hashlib.sha256(blob).hexdigest(), archive_bytes=len(blob)))
        if repo in ['UCell', 'scGate']:
            with log.open('ab') as f:
                p = subprocess.run([str(R), 'CMD', 'INSTALL', '--no-multiarch', '--no-docs', '--library=' + str(lib), str(source)], stdout=f, stderr=f, env=env, timeout=180)
                assert p.returncode == 0, 'TASK_LOCAL_OFFICIAL_PACKAGE_INSTALL_FAILED'
        assert sum(p.stat().st_size for p in out.rglob('*') if p.is_file()) < 12 * 1024**2
        save('private_progress.json', dict(stage='TASK_LOCAL_SOFTWARE_ONLY', completed_repository=repo, observed_UTC=datetime.now(timezone.utc).isoformat()))
    smoke = out / 'official_synthetic_compatibility.R'; smoke.write_text(SMOKE)
    with log.open('ab') as f:
        p = subprocess.run([str(Rscript), '--vanilla', str(smoke), str(lib), str(runtime), str(roots['scGate_models'])], stdout=f, stderr=f, env=env, timeout=240)
        assert p.returncode == 0, 'OFFICIAL_SYNTHETIC_COMPATIBILITY_FAILED_NO_REAL_CELLS'
    assert sum(p.stat().st_size for p in out.rglob('*') if p.is_file()) < 12 * 1024**2
    save('private_public_software_archive_provenance.json', records)
    receipt = dict(status='ACCEPTED_TASK_LOCAL_OFFICIAL_SCGATE_UCELL_AND_PREDEFINED_MODELS_SYNTHETIC_ONLY', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), public_software_archives=records, scGate='1.7.3', UCell='2.17.1', predefined_model_commit=SOURCES[-1][1], synthetic_cells=240, real_cell_RNA_scored=False, shared_dependencies_modified=False, expression_data_downloaded=False, resource_guard=guard, output_bytes=sum(p.stat().st_size for p in out.rglob('*') if p.is_file()))
    save('private_terminal.json', receipt)
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
