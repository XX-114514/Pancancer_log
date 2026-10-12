#!/usr/bin/env python3
"""Preserve small private filtering/audit evidence beside the final Zarr.

Run once after the full assembly audit passes. No matrices, R exchange files,
credentials or unrelated inputs are selected; nothing is removed.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil

from merge_public_preprocessed_zarr import effective_capacity, RESERVE


def digest(path):
    checksum=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):checksum.update(block)
    return checksum.hexdigest()


def preserve(preprocess,collection,final_root):
    terminal=json.loads((final_root/'manifests/assembly_terminal.json').read_text())
    if not terminal['all_checks_pass'] or not terminal['status'].startswith('completed_audited_'):
        raise ValueError('Require full final Zarr audit before preserving terminal evidence')
    original=json.loads((preprocess/'manifests/independent_preparation_audit.json').read_text())
    if not original['all_checks_pass'] or original['n_datasets']!=9:raise ValueError('Require original full audit')
    target_zarr=final_root/'outputs/public_preprocessed_union_v1.zarr'
    if target_zarr.is_symlink():raise ValueError('Final Zarr is an external link')
    for directory,folders,files in os.walk(target_zarr):
        if any((Path(directory)/name).is_symlink() for name in folders+files):
            raise ValueError('Final Zarr contains external links')
    selected=[]
    def add(path,relative,required=True):
        if not path.is_file():
            if required:raise FileNotFoundError(path)
            return
        selected.append((path,relative))
    for name in ['execution_start.json','execution_terminal.json','independent_preparation_audit.json',
                 'operator_storage_hold_v1.json','operator_storage_resume_v1.json']:
        add(preprocess/'manifests'/name,Path('execution')/name)
    for name in ['preprocess_doublet_ambient.py','run_ambient_preparation.py',
                 'decontx_library_worker.R','audit_ambient_preparation.py']:
        add(preprocess/'scripts'/name,Path('code')/name)
    for source in original['datasets']:
        name=source['dataset'];folder=preprocess/'objects'/name
        add(folder/'dataset_summary.json',Path('sources')/name/'dataset_summary.json')
        if (folder/'native_QC_input_contract.json').exists():
            for filename in ['native_QC_input_contract.json','native_QC_ledger_private.tsv.gz']:
                add(folder/filename,Path('sources')/name/filename)
        else:
            for filename in ['accepted_'+name+'.json',name+'_cell_qc_and_filter_reason.tsv.gz']:
                add(collection/'manifests'/filename,Path('sources')/name/filename)
    for name in ['baseline_manifest.json','baseline_fingerprint_complete.json','fingerprint_gene_vocabulary.json']:
        add(collection/'manifests'/name,Path('baseline')/name)
    for name in ['identifier_audit_repair_start_v4.json','identifier_audit_repair_terminal_v4.json',
                 'independent_preparation_audit_identifier_v4.json']:
        add(preprocess/'manifests'/name,Path('audit_repair')/name,required=False)
    for name in ['independent_preparation_audit.json','execution_terminal.json','audit_ambient_preparation.py']:
        add(preprocess/'config/identifier_audit_repair_v4'/name,Path('audit_repair/original_failed_v3')/name,required=False)
    total=sum(path.stat().st_size for path,_ in selected)
    if total>1_000_000_000:raise ValueError('Small-evidence contract exceeds 1GB; preserve attempt, review explicitly')
    free,capacity=effective_capacity(final_root)
    if free<total+RESERVE:raise RuntimeError('Cannot preserve evidence and final-volume 500GB reserve')
    destination=final_root/'provenance';destination.mkdir(mode=0o700)
    rows=[]
    for source,relative in selected:
        target=destination/relative;target.parent.mkdir(parents=True,exist_ok=True,mode=0o700)
        if target.exists():raise ValueError('Existing evidence attempt preserved')
        checksum=digest(source);shutil.copyfile(source,target)
        if digest(target)!=checksum:raise ValueError('Provenance copy checksum mismatch')
        rows.append(dict(relative_path=str(relative),bytes=source.stat().st_size,sha256=checksum))
    private=dict(observed_utc=datetime.now(timezone.utc).isoformat(),all_checks_pass=True,
                 n_files=len(rows),total_bytes=total,files=rows,temporary_inputs_deleted=False)
    with (destination/'private_evidence_index.json').open('x') as stream:json.dump(private,stream,indent=2)
    summary={k:v for k,v in private.items() if k!='files'}
    summary.update(schema='public_preprocessing_durable_evidence_v1',source_final_AC=False,
        final_reserve_bytes=RESERVE,copy_script_sha256=digest(Path(__file__)),
        final_audit_sha256=digest(final_root/'manifests/independent_final_audit.json'),
        standalone_Zarr_arrays_no_external_links=True)
    with (final_root/'manifests/durable_evidence_summary.json').open('x') as stream:json.dump(summary,stream,indent=2)
    return summary


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    for name in ['preprocess-root','collection','final-root']:parser.add_argument('--'+name,required=True)
    args=parser.parse_args()
    print(json.dumps(preserve(Path(args.preprocess_root),Path(args.collection),Path(args.final_root))))
