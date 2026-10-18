"""Restore author labels by exact library-namespace/barcode, preserving authorities."""
import io
import json
import os
import subprocess
import sys
from datetime import datetime,timezone
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def original_frames(root):
    original=root.parents[2]/'Project_v2/Downstream/0716_Data/raw_geo_supplements/GSE176021'
    with h5py.File(original/'GSE176021_rawdata.h5ad','r') as h:
        o=h['obs'];idx=np.asarray(read_elem(o[str(o.attrs['_index'])])).astype(str)
        raw=pd.DataFrame({f:np.asarray(read_elem(o[f])).astype(str) for f in ['barcode','imid','patient_id']},index=pd.Index(idx,name='raw_index'))
    assert len(raw)==902463 and raw.index.is_unique
    return original,raw


def author_frame(path):
    R=Path(sys.executable).resolve().parents[2]/'scrna_r/bin/Rscript'
    script='p<-commandArgs(TRUE)[1];con<-pipe(paste("gzip -dc",shQuote(p)),"rb");x<-readRDS(gzcon(con));x$author_row_name<-rownames(x);write.table(x,stdout(),sep="\\t",quote=FALSE,row.names=FALSE,na="NA");close(con)'
    r=subprocess.run([str(R),'--vanilla','-e',script,str(path)],capture_output=True,timeout=120)
    assert r.returncode==0 and len(r.stdout)<256*1024**2,'AUTHOR_METADATA_READ_FAILED'
    return pd.read_csv(io.BytesIO(r.stdout),sep='\t',keep_default_na=False,dtype={'barcode':str,'imid':str,'CellType':str,'author_row_name':str})


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,16*1024**2)
    pre=root/'GSE176021_original_author_annotation_key_qualification_20261019_v1/private_terminal.json';assert pre.exists()
    out=root/'GSE176021_exact_original_author_annotation_sidecar_20261019_v2';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_any_label_sidecar.json').write_text(json.dumps({'created_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'resource_guard':guard,'preflight_sha256':sha(pre),'scope':'Both entire CD3/CD8 original annotation tables; exact original raw imid + underscore + barcode. Raw/author keys must be unique, patient exact. All unmatched retained; no fuzzy/library alias or barcode-only rescue. Current atlas exact original source row and cell-ID bridge.','outputs':'Private metadata-only sidecar plus complete aggregate coverage, no original authority mutation or new scores/P. Author reference labels, no formal state adoption. Figure erratum DOI10.1038/s41586-021-03893-6 tracked, no inferred relabel.'},indent=2)+'\n')
    original,raw=original_frames(root);key=pd.Index(raw.imid+'_'+raw.barcode);assert key.is_unique
    atlas=root.parents[1]/'runs/20261013_immunotherapy_atlas_v1';a=zarr.open_group(str(atlas/'outputs/immune_validation_authority_v1.zarr'),mode='r');ds=np.asarray(read_elem(a['obs/dataset'])).astype(str);global_rows=np.flatnonzero(ds=='GSE176021');source_rows=np.load(atlas/'manifests/GSE176021_source_rows.npy');native_rows=np.asarray(read_elem(a['obs/source_row_position']))[global_rows].astype(np.int64);orig=np.asarray(read_elem(a['obs/parent_cell_id'])).astype(str)[global_rows]
    assert len(global_rows)==868494 and np.array_equal(source_rows,native_rows) and np.array_equal(raw.index.to_numpy(str)[source_rows],orig)
    tissue=np.asarray(read_elem(a['obs/tissue'])).astype(str)[global_rows];labels={};arrays={'global_atlas_row':global_rows.astype(np.int32),'original_raw_row':source_rows.astype(np.int32)};stats=[]
    for branch,nexpected,nmatched in [('CD3',560916,552499),('CD8',235851,233249)]:
        path=original/('GSE176021_'+branch+'_annotations.rds.gz');d=author_frame(path);assert len(d)==nexpected and d.barcode.is_unique and np.array_equal(d.barcode,d.author_row_name)
        loc=key.get_indexer(d.barcode);matched=loc>=0;assert matched.sum()==nmatched and len(np.unique(loc[matched]))==nmatched
        assert np.array_equal(raw.patient_id.to_numpy(str)[loc[matched]],d.imid.to_numpy(str)[matched])
        pos=np.full(len(raw),-1,np.int32);pos[loc[matched]]=np.flatnonzero(matched);current=pos[source_rows];v=current>=0
        names=sorted(d.CellType.unique());lookup=pd.Index(names);code=lookup.get_indexer(d.CellType).astype(np.int16);current_code=np.full(len(current),-1,np.int16);current_code[v]=code[current[v]]
        arrays[branch+'_author_row']=current;arrays[branch+'_label_code']=current_code;arrays[branch+'_unmatched_author_rows']=np.flatnonzero(~matched).astype(np.int32)
        labels[branch]={'labels':names,'source_sha256':sha(path),'rows':len(d),'exact_raw_matches':int(matched.sum()),'unmatched_original_author_rows':int((~matched).sum()),'atlas_labelled':int(v.sum()),'atlas_unlabelled':int((~v).sum()),'values_status':'Original author reference values; not formal current state labels. Preserve original values; figure erratum tracked separately.'}
        for t in sorted(set(tissue)):
            mask=tissue==t
            for c,name in enumerate(names):stats.append({'branch':branch,'tissue':t,'author_label':name,'labelled_cells':int((mask&(current_code==c)).sum()),'all_atlas_source_cells_in_tissue':int(mask.sum())})
            stats.append({'branch':branch,'tissue':t,'author_label':'UNMAPPED_OR_OUTSIDE_AUTHOR_ANNOTATION','labelled_cells':int((mask&~v).sum()),'all_atlas_source_cells_in_tissue':int(mask.sum())})
    np.savez_compressed(out/'private_all_atlas_cells_exact_author_bridge.npz',**arrays)
    (out/'private_author_label_lookup_and_coverage.json').write_text(json.dumps(labels,indent=2)+'\n');pd.DataFrame(stats).to_csv(out/'private_all_source_tissue_author_reference_coverage.tsv',sep='\t',index=False)
    native=root/'immune_RNA_MCP_all_atlas_native_measurement_20261019_v2';member=np.load(native/'private_all_captured_membership.npz');pools=pd.read_csv(native/'private_all_capture_pools.tsv',sep='\t');pools=pools[pools.sampleID.str.startswith('GSE176021::')];pool_idx=member['pool_index'];gr=member['global_row'];mapping=pd.Index(global_rows);rows=[]
    for p in pools.itertuples(index=False):
        ids=gr[pool_idx==p.pool_index];ix=mapping.get_indexer(ids);assert (ix>=0).all();ts=set(tissue[ix]);assert len(ts)==1
        for branch in ['CD3','CD8']:
            codes=arrays[branch+'_label_code'][ix]
            for c,name in enumerate(labels[branch]['labels']):rows.append({'tissue':next(iter(ts)),'view':getattr(p,'view'),'branch':branch,'author_label':name,'cells':int((codes==c).sum())})
            rows.append({'tissue':next(iter(ts)),'view':getattr(p,'view'),'branch':branch,'author_label':'UNMAPPED_OR_OUTSIDE_AUTHOR_ANNOTATION','cells':int((codes<0).sum())})
    coverage=pd.DataFrame(rows).groupby(['tissue','view','branch','author_label'],sort=True).cells.sum().reset_index();coverage.to_csv(out/'private_all_captured_view_author_reference_coverage.tsv',sep='\t',index=False)
    files={p.name:{'sha256':sha(p),'bytes':p.stat().st_size} for p in out.iterdir() if p.is_file()};assert sum(v['bytes'] for v in files.values())<16*1024**2
    result={'status':'COMPLETE_EXACT_AUTHOR_LABEL_SIDECAR_ALL_RAW_AND_ATLAS_ROWS_AUDIT_PENDING','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'resource_guard':guard,'source_atlas_cells':868494,'raw_cells':902463,'files':files,'authority_mutated':False,'scores_or_clinical_associations':False}
    (out/'private_terminal.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':result['status'],'terminal_sha256':sha(out/'private_terminal.json')}))


if __name__=='__main__':main()
