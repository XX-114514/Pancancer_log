"""Reconstruct every author label by independent dictionary keys and quoted R export."""
import io
import json
import os
import subprocess
import sys
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from immune_RNA_MCP_new_source_io_v1 import wait


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,8*1024**2)
    source=root/'GSE176021_exact_original_author_annotation_sidecar_20261019_v2';out=root/'GSE176021_exact_original_author_annotation_audit_20261019_v2';assert not out.exists();out.mkdir(mode=0o700)
    terminal=wait(root,out,'restore_GSE176021_exact_author_annotation_sidecar_v2',source/'private_terminal.json',{'COMPLETE_EXACT_AUTHOR_LABEL_SIDECAR_ALL_RAW_AND_ATLAS_ROWS_AUDIT_PENDING'})
    for name,v in terminal['files'].items():assert sha(source/name)==v['sha256']
    original=root.parents[2]/'Project_v2/Downstream/0716_Data/raw_geo_supplements/GSE176021'
    with h5py.File(original/'GSE176021_rawdata.h5ad','r') as h:
        obs=h['obs'];rawids=np.asarray(read_elem(obs[str(obs.attrs['_index'])])).astype(str);barcode=np.asarray(read_elem(obs['barcode'])).astype(str);imid=np.asarray(read_elem(obs['imid'])).astype(str);patient=np.asarray(read_elem(obs['patient_id'])).astype(str)
    raw_keys={lib+'_'+bc:i for i,(lib,bc) in enumerate(zip(imid,barcode))};assert len(raw_keys)==902463
    atlas=root.parents[1]/'runs/20261013_immunotherapy_atlas_v1';a=zarr.open_group(str(atlas/'outputs/immune_validation_authority_v1.zarr'),mode='r');global_rows=np.flatnonzero(np.asarray(read_elem(a['obs/dataset'])).astype(str)=='GSE176021');rows=np.asarray(read_elem(a['obs/source_row_position']))[global_rows].astype(np.int64);ids=np.asarray(read_elem(a['obs/parent_cell_id'])).astype(str)[global_rows];assert np.array_equal(rawids[rows],ids) and np.array_equal(rows,np.load(atlas/'manifests/GSE176021_source_rows.npy'))
    tissue=np.asarray(read_elem(a['obs/tissue'])).astype(str)[global_rows];saved=np.load(source/'private_all_atlas_cells_exact_author_bridge.npz');assert np.array_equal(saved['global_atlas_row'],global_rows) and np.array_equal(saved['original_raw_row'],rows)
    reference=json.loads((source/'private_author_label_lookup_and_coverage.json').read_text());all_codes={};metadata_counts={}
    for branch,expected,expected_matched in [('CD3',560916,552499),('CD8',235851,233249)]:
        path=original/('GSE176021_'+branch+'_annotations.rds.gz');assert sha(path)==reference[branch]['source_sha256'];R=Path(sys.executable).resolve().parents[2]/'scrna_r/bin/Rscript'
        script='p<-commandArgs(TRUE)[1];con<-pipe(paste("gzip -dc",shQuote(p)),"rb");d<-readRDS(gzcon(con));stopifnot(identical(rownames(d),as.character(d$barcode)));write.csv(d[,c("barcode","imid","CellType")],stdout(),row.names=FALSE);close(con)'
        r=subprocess.run([str(R),'--vanilla','-e',script,str(path)],capture_output=True,timeout=120);assert r.returncode==0
        author=pd.read_csv(io.BytesIO(r.stdout),keep_default_na=False,dtype=str);assert len(author)==expected and author.barcode.is_unique
        hits=np.array([raw_keys.get(x,-1) for x in author.barcode],np.int64);valid=hits>=0;assert valid.sum()==expected_matched and np.unique(hits[valid]).size==valid.sum()
        assert np.array_equal(patient[hits[valid]],author.imid.to_numpy(str)[valid])
        rowmap={int(rawpos):int(authorpos) for authorpos,rawpos in enumerate(hits) if rawpos>=0};current=np.array([rowmap.get(int(x),-1) for x in rows],np.int32);assert np.array_equal(current,saved[branch+'_author_row']) and np.array_equal(np.flatnonzero(~valid),saved[branch+'_unmatched_author_rows'])
        labels=reference[branch]['labels'];assert labels==sorted(set(author.CellType));code_by_name={name:i for i,name in enumerate(labels)};code=np.array([code_by_name[author.CellType.iloc[j]] if j>=0 else -1 for j in current],np.int16);assert np.array_equal(code,saved[branch+'_label_code']);all_codes[branch]=code
        assert reference[branch]['atlas_labelled']==int((code>=0).sum()) and reference[branch]['atlas_unlabelled']==int((code<0).sum());metadata_counts[branch]={'exact_raw_matches':int(valid.sum()),'unmatched_author_rows':int((~valid).sum()),'current_atlas_labelled':int((code>=0).sum())}
    coverage=pd.read_csv(source/'private_all_source_tissue_author_reference_coverage.tsv',sep='\t',keep_default_na=False)
    for r in coverage.itertuples(index=False):
        codes=all_codes[r.branch];mask=tissue==r.tissue;label=-1 if r.author_label=='UNMAPPED_OR_OUTSIDE_AUTHOR_ANNOTATION' else reference[r.branch]['labels'].index(r.author_label)
        assert r.labelled_cells==int((mask&(codes==label)).sum()) and r.all_atlas_source_cells_in_tissue==int(mask.sum())
    native=root/'immune_RNA_MCP_all_atlas_native_measurement_20261019_v2';pools=pd.read_csv(native/'private_all_capture_pools.tsv',sep='\t');pools=pools[pools.sampleID.str.startswith('GSE176021::')];member=np.load(native/'private_all_captured_membership.npz');poolids=member['pool_index'];gr=member['global_row'];byglobal={int(x):i for i,x in enumerate(global_rows)};counts=Counter();sample_views={}
    for p in pools.itertuples(index=False):
        original_rows=gr[poolids==p.pool_index];ix=np.array([byglobal[int(x)] for x in original_rows]);assert len(set(tissue[ix]))==1;t=tissue[ix[0]];sample_views.setdefault(p.sampleID,{})[getattr(p,'view')]=set(map(int,original_rows))
        for branch in ['CD3','CD8']:
            counts.update((t,getattr(p,'view'),branch,reference[branch]['labels'][c] if c>=0 else 'UNMAPPED_OR_OUTSIDE_AUTHOR_ANNOTATION') for c in all_codes[branch][ix])
    cv=pd.read_csv(source/'private_all_captured_view_author_reference_coverage.tsv',sep='\t',keep_default_na=False)
    for r in cv.itertuples(index=False):assert r.cells==counts[(r.tissue,r.view,r.branch,r.author_label)]
    assert int(cv.cells.sum())==sum(counts.values())
    overlap=[]
    for views in sample_views.values():
        R=views.get('R_FOXP3_capture',set());E=views.get('E_effector_capture',set());overlap.append({'R_cells':len(R),'E_cells':len(E),'shared_cells':len(R&E),'R_only_cells':len(R-E),'E_only_cells':len(E-R)})
    overlap_summary={k:sum(x[k] for x in overlap) for k in ['R_cells','E_cells','shared_cells','R_only_cells','E_only_cells']}
    receipt={'status':'ACCEPTED_ALL_RAW_AUTHOR_KEYS_PATIENTS_CURRENT_ATLAS_LABELS_AND_FULL_COVERAGE','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'producer_terminal_sha256':sha(source/'private_terminal.json'),'resource_guard':guard,'metadata_counts':metadata_counts,'all_current_source_rows':868494,'all_original_author_rows':[560916,235851],'capture_view_overlap':overlap_summary,'clinical_labels_or_effects':False,'authority_mutated':False,'label_role':'Original author reference values, not formal current subtype/state adoption; original correction remains a reference qualification caveat.'}
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'status':receipt['status'],'terminal_sha256':sha(out/'private_terminal.json')}))


if __name__=='__main__':main()
