"""Outcome-blind fixed official scGate pilot; native full-gene arrays RAM only."""
import gzip
import io
import json
import os
import struct
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from gamma_delta_RNA_rules_v1 import MARKERS, flags
from prepare_NSCLC_full_native_gene_units_v1 import TrackedReader, chunks
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource

R_CODE = r'''
args<-commandArgs(TRUE)
.libPaths(c(args[1],file.path(args[2],c('library_v7','library_v5','library_v2','library')),.libPaths()))
stopifnot(as.character(packageVersion('scGate'))=='1.7.3',as.character(packageVersion('UCell'))=='2.17.1')
suppressPackageStartupMessages(library(scGate))
con<-file('stdin','rb');d<-readBin(con,'integer',n=3,size=4,endian='little');stopifnot(length(d)==3)
i<-readBin(con,'integer',n=d[3],size=4,endian='little')+1L
j<-readBin(con,'integer',n=d[3],size=4,endian='little')+1L
v<-readBin(con,'double',n=d[3],size=8,endian='little')
expectedRNA<-readBin(con,'double',n=d[2],size=8,endian='little')
expectedGenes<-readBin(con,'integer',n=d[2],size=4,endian='little');close(con)
stopifnot(length(i)==d[3],length(j)==d[3],length(v)==d[3],all(is.finite(v)),all(v>0),all(v==floor(v)))
genes<-toupper(trimws(read.csv(args[4],check.names=FALSE,stringsAsFactors=FALSE,na.strings=character(0))$geneSymbol))
stopifnot(length(genes)==d[1],!anyDuplicated(genes))
x<-Matrix::sparseMatrix(i=i,j=j,x=v,dims=d[1:2],dimnames=list(genes,paste0('anonymous_cell_',seq_len(d[2]))))
stopifnot(identical(as.numeric(Matrix::colSums(x)),as.numeric(expectedRNA)),identical(as.integer(Matrix::colSums(x>0)),as.integer(expectedGenes)))
rm(i,j,v);gc()
obj<-Seurat::CreateSeuratObject(counts=x,min.cells=0,min.features=0)
stopifnot(nrow(obj)==length(genes),!anyDuplicated(rownames(obj)))
obj<-Seurat::NormalizeData(obj,normalization.method='LogNormalize',scale.factor=10000,verbose=FALSE)
loadmodel<-function(collection,name){scGate::load_scGate_model(file.path(args[3],'human',collection,paste0(name,'_scGate_Model.tsv')),master.table=file.path(args[3],'human',collection,'master_table.tsv'))}
full<-loadmodel('HiTME','gdT');pbmc<-loadmodel('PBMC','gdT');nk<-loadmodel('HiTME','NK')
present<-function(m){all(sub('[+-]$','',unlist(strsplit(m$signature,';',fixed=TRUE)))%in%genes)}
stopifnot(present(pbmc),present(nk),!present(full))
subset<-full
for(k in seq_len(nrow(subset))){
  old<-strsplit(subset$signature[k],';',fixed=TRUE)[[1]]
  keep<-sub('[+-]$','',old)%in%genes
  stopifnot(sum(keep)>=1,mean(keep)>=.8)
  subset$signature[k]<-paste(old[keep],collapse=';')
}
stopifnot(present(subset))
models<-list(PBMC_exact=pbmc,HiTME_native_subset=subset,NK_HiTME_exact=nk)
answer<-data.frame(position=seq_len(ncol(obj))-1L)
answer$HiTME_full<-'UNMEASURABLE'
for(name in names(models)){
  result<-scGate::scGate(obj,model=models[[name]],seed=123,pos.thr=.2,neg.thr=.2,maxRank=1500,nfeatures=2000,pca.dim=30,k.param=30,min.cells=30,ncores=1,return.CellOntology=FALSE,progressbar=FALSE,verbose=FALSE)
  stopifnot(ncol(result)==d[2],all(result$is.pure%in%c('Pure','Impure')))
  answer[[name]]<-as.character(result$is.pure)
  rm(result);gc()
}
cat('BEGIN_CLASSIFICATION_TSV\n')
write.table(answer,stdout(),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
'''


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT']); guard = resource(root, 4 * 1024**2)
    runtime = root / 'gamma_scGate_official_runtime_qualification_20261019_v1'
    assert sha(runtime / 'private_terminal.json') == '8b69085d7dbe3be4b7342d58511066b2804eef4f7618a362112f4a1e24379058'
    gene_scope = root / 'gamma_scGate_native_model_gene_qualification_20261019_v1/private_complete_fixed_two_model_native_qualification.json'
    assert sha(gene_scope) == '1f746d64ff2eeb601f9b8bfd42802b4a4b96f218dbff1fb044c697258ead7778'
    native = root / 'NSCLC_full_native_gene_units_20261016_v1'
    unit = np.load(native / 'private_cell_unit_index.npy', allow_pickle=False)
    RNA = np.load(native / 'private_full_RNA.npy', allow_pickle=False)
    detected = np.load(native / 'private_detected_genes.npy', allow_pickle=False)
    u = pd.read_csv(native / 'private_all_native_units_and_quality.tsv', sep='\t', keep_default_na=False).sort_values('unit_index')
    assert len(unit) == len(RNA) == len(detected) == 1254749 and len(u) == 11701
    selection = np.sort(np.random.default_rng(2026101901).choice(len(unit), 20000, replace=False)).astype(np.int32)
    nnnz = int(detected[selection].sum())
    memory = int(next(s.split()[1] for s in Path('/proc/meminfo').read_text().splitlines() if s.startswith('MemAvailable:'))) * 1024
    assert memory > nnnz * 96 + 10 * 1024**3
    out = root / 'NSCLC_gamma_scGate_fixed_pilot_20261019_v1'; assert not out.exists(); out.mkdir(mode=0o700)
    save = lambda name, value: (out / name).write_text(json.dumps(value, indent=2) + '\n')
    buffer = io.BytesIO(); np.savez_compressed(buffer, native_row=selection)
    (out / 'private_outcome_free_uniform_native_selection.npz').write_bytes(buffer.getvalue())
    source = root / 'count_matrix_download_v1/download_state.json'; download = json.loads(source.read_text()); matrix = Path(download['path'])
    gene_path = root / 'GSE243013_genes.csv.gz'; genes = pd.read_csv(gene_path, keep_default_na=False).geneSymbol.str.strip().str.upper().to_numpy(str)
    assert len(genes) == len(set(genes)) == 31831 and (genes != '').all() and matrix.stat().st_size == download['bytes']
    models = runtime / 'scGate_models-41a45cd3f8bb5f5a7daf21ec276f6a726f6ee0d4'
    assert models.is_dir() and (runtime / 'library/scGate').is_dir() and (runtime / 'library/UCell').is_dir()
    save('plan_before_real_fixed_scGate_pilot_RNA.json', dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), R_source_sha256=__import__('hashlib').sha256(R_CODE.encode()).hexdigest(), resource_guard=guard, runtime_acceptance_sha256=sha(runtime / 'private_terminal.json'), native_model_scope_sha256=sha(gene_scope), selection_sha256=sha(out / 'private_outcome_free_uniform_native_selection.npz'), seed=2026101901, cells=20000, native_gene_universe=len(genes), selected_nonzero_bound=nnnz, RAM_guard_bytes=nnnz * 96 + 10 * 1024**3, RAM_available_bytes=memory, source_matrix_sha256=download['sha256'], models=['HiTME_full_unmeasurable', 'PBMC_exact', 'HiTME_native_measured_subset', 'NK_HiTME_exact'], all_model_slots_retained=True, scope='Uniform cells from all native source cells, natural prevalence; no author class, patient, treatment, outcome, LR or function chosen. Author reference used only after scoring. Official defaults fixed and no threshold search. Full HiTME not run because of physical missing marker; explicit measured-subset adaptation and PBMC tissue-domain restriction retained.', full_RNA_or_score_matrix_persisted=False, clinical_or_formal_identity_adoption=False))
    inverse = np.full(len(unit), -1, np.int32); inverse[selection] = np.arange(len(selection))
    totals = np.zeros(len(unit), np.int64); ngenes = np.zeros(len(unit), np.int32)
    parts = []; entries = 0; last = -1; checkpoint = -1; start = time.monotonic()
    with matrix.open('rb') as binary:
        tracked = TrackedReader(binary)
        with gzip.GzipFile(fileobj=tracked, mode='rb') as stream:
            assert stream.readline().startswith(b'%%MatrixMarket matrix coordinate real general')
            line = stream.readline()
            while line.startswith(b'%'): line = stream.readline()
            dims = tuple(map(int, line.split())); assert dims == (len(unit), len(genes), 2010550708)
            for a in chunks(stream):
                assert np.isfinite(a).all() and np.array_equal(a, np.floor(a)) and (a[:, 2] > 0).all() and a[:, 2].max() < 2**31
                row = a[:, 0].astype(np.int64) - 1; col = a[:, 1].astype(np.int64) - 1; val = a[:, 2].astype(np.int32)
                coord = row * len(genes) + col; assert coord[0] > last and (np.diff(coord) > 0).all(); last = int(coord[-1])
                lo = int(row[0]); hi = int(row[-1]) + 1
                totals[lo:hi] += np.bincount(row - lo, weights=val, minlength=hi - lo).astype(np.int64)
                ngenes[lo:hi] += np.bincount(row - lo, minlength=hi - lo).astype(np.int32)
                take = inverse[row] >= 0
                parts.append((inverse[row[take]], col[take].astype(np.int32), val[take])); entries += len(a)
                if entries // 100000000 != checkpoint:
                    checkpoint = entries // 100000000
                    save('private_progress.json', dict(stage='REPRESENTATIVE_FULL_NATIVE_RNA_EXTRACTION_RAM_ONLY', entries=entries, expected_entries=dims[2], elapsed_seconds=time.monotonic() - start, observed_UTC=datetime.now(timezone.utc).isoformat()))
    assert entries == dims[2] and tracked.bytes == download['bytes'] and tracked.digest.hexdigest() == download['sha256']
    np.testing.assert_array_equal(totals, RNA); np.testing.assert_array_equal(ngenes, detected)
    rr = np.concatenate([p[0] for p in parts]); cc = np.concatenate([p[1] for p in parts]); vv = np.concatenate([p[2] for p in parts]); del parts
    x = coo_matrix((vv, (rr, cc)), shape=(len(selection), len(genes))).tocsr(); del rr, cc, vv
    assert x.nnz == nnnz and x.has_canonical_format
    hard = dict(np.load(root / 'NSCLC_gamma_delta_RNA_candidates_20261019_v1/private_candidate_identity_sidecar.npz', allow_pickle=False))
    oldbits = np.zeros(len(unit), np.uint8); oldbits[hard['native_row']] = hard['rule_bits']
    np.testing.assert_array_equal(flags(x[:, pd.Index(genes).get_indexer(MARKERS)].toarray()), oldbits[selection])
    coords = x.tocoo()
    payload = struct.pack('<iii', len(genes), len(selection), x.nnz) + coords.col.astype('<i4').tobytes() + coords.row.astype('<i4').tobytes() + coords.data.astype('<f8').tobytes() + RNA[selection].astype('<f8').tobytes() + detected[selection].astype('<i4').tobytes()
    worker = out / 'official_fixed_pilot_worker.R'; worker.write_text(R_CODE)
    rt = root.parents[2] / 'Project_v5/runs/20261014_DIALOGUE_method_pilot_v1'
    models = runtime / 'scGate_models-41a45cd3f8bb5f5a7daf21ec276f6a726f6ee0d4'
    Rscript = Path(sys.executable).resolve().parents[2] / 'scrna_r/bin/Rscript'
    save('private_progress.json', dict(stage='OFFICIAL_FIXED_SCGATE_THREE_MODELS_RAM_ONLY', observed_UTC=datetime.now(timezone.utc).isoformat()))
    with (out / 'private_official_R_worker.log').open('xb') as log:
        p = subprocess.run([str(Rscript), '--vanilla', str(worker), str(runtime / 'library'), str(rt), str(models), str(gene_path)], input=payload, stdout=subprocess.PIPE, stderr=log, timeout=1800)
    assert p.returncode == 0, 'OFFICIAL_REAL_PILOT_FAILED_NO_MODEL_RETUNING'
    marker = b'BEGIN_CLASSIFICATION_TSV\n'; assert p.stdout.count(marker) == 1
    text = p.stdout.split(marker)[1]; assert len(text) < 4 * 1024**2
    classes = pd.read_csv(io.BytesIO(text), sep='\t', keep_default_na=False)
    assert classes.position.tolist() == list(range(20000))
    names = ['HiTME_full', 'PBMC_exact', 'HiTME_native_subset', 'NK_HiTME_exact']
    assert classes.HiTME_full.eq('UNMEASURABLE').all()
    for name in names[1:]: assert classes[name].isin(['Pure', 'Impure']).all()
    metadata = u.iloc[unit[selection]].reset_index(drop=True)
    records = []
    for name in names:
        for (major, subtype), group in metadata.groupby(['major_cell_type', 'sub_cell_type'], sort=True):
            values = classes.loc[group.index, name]
            records.append(dict(source='GSE243013', model=name, major_cell_type=major, sub_cell_type=subtype, pilot_cells=len(group), pure_cells=int(values.eq('Pure').sum()) if name != 'HiTME_full' else None, model_status='NATIVE_GENE_UNMEASURABLE' if name == 'HiTME_full' else 'FIXED_MODEL_EXECUTED_NOT_FORMAL_IDENTITY'))
    pd.DataFrame(records).to_csv(out / 'private_complete_author_type_model_classification_summary.tsv', sep='\t', index=False)
    buffer = io.BytesIO(); np.savez_compressed(buffer, native_row=selection, models=np.array(names), classes=classes[names].to_numpy(str))
    (out / 'private_complete_pilot_candidate_classification_sidecar.npz').write_bytes(buffer.getvalue())
    save('private_full_input_and_output_transfer_proofs.json', dict(raw_matrix_sha256=tracked.digest.hexdigest(), raw_entries=entries, all_native_RNA_and_detection_exact=True, selected_all_gene_nonzeros=x.nnz, binary_input_sha256=__import__('hashlib').sha256(payload).hexdigest(), worker_stdout_sha256=__import__('hashlib').sha256(p.stdout).hexdigest(), classification_TSV_sha256=__import__('hashlib').sha256(text).hexdigest(), all_selected_old_hard_rule_bits_exact=True, anonymous_R_cells=True, clinical_or_author_labels_supplied_to_R=False, no_full_RNA_or_score_matrix_persisted=True))
    files = {f.name: dict(sha256=sha(f), bytes=f.stat().st_size) for f in out.iterdir() if f.is_file()}
    assert sum(v['bytes'] for v in files.values()) < 3 * 1024**2
    receipt = dict(status='COMPLETE_REPRESENTATIVE_FIXED_SCGATE_PILOT_ALL_MODEL_SLOTS_STRUCTURAL_REVIEW_PENDING', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), files=files, real_pilot_cells=20000, all_model_slots=4, executed_models=3, clinical_or_formal_identity_adoption=False, scientific_outputs_public=False)
    save('private_terminal.json', receipt)
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
