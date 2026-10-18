"""Read double-gzip author metadata and qualify exact cell-key bridges, no scores."""
import hashlib
import io
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from anndata._io.specs import read_elem
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,8*1024**2)
    original=root.parents[2]/'Project_v2/Downstream/0716_Data/raw_geo_supplements/GSE176021'
    assert original.is_dir()
    out=root/'GSE176021_original_author_annotation_key_qualification_20261019_v1'
    assert not out.exists();out.mkdir(mode=0o700)
    plan={'created_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'resource_guard':guard,'scope':'All original author CD3/CD8 rows and raw H5AD metadata only. No expression, programme score, clinical label or association read. No authority mutation. Exact keys only, no fuzzy or barcode-only cross-library rescue.','format':'External gzip wraps original compressed RDS; external gzip stream then R gzcon reads inner RDS.','original_paper':'10.1038/s41586-021-03752-4','author_correction':'10.1038/s41586-021-03893-6; figure label swaps not automatically applied to 2024 cached annotation values.'}
    (out/'plan_before_metadata_join.json').write_text(json.dumps(plan,indent=2)+'\n')
    with h5py.File(original/'GSE176021_rawdata.h5ad','r') as h:
        o=h['obs'];idx=np.asarray(read_elem(o[str(o.attrs['_index'])])).astype(str)
        fields=['original_cell_id','barcode','imid','gsm_id','patient_id','sample_label','cd3_CellType','cd8_CellType']
        raw=pd.DataFrame({f:np.asarray(read_elem(o[f])).astype(str) for f in fields},index=pd.Index(idx,name='raw_index'))
    assert len(raw)==902463 and raw.index.is_unique
    missing={'','nan','na','none','unknown','not_available','not_assessed','unassigned','<na>'}
    empty_counts={f:int((~raw[f].str.lower().isin(missing)).sum()) for f in ['cd3_CellType','cd8_CellType']}
    R=Path(sys.executable).resolve().parents[2]/'scrna_r/bin/Rscript';assert R.is_file()
    rcode='p<-commandArgs(TRUE)[1];con<-pipe(paste("gzip -dc",shQuote(p)),"rb");x<-readRDS(gzcon(con));stopifnot(is.data.frame(x));x$author_row_name<-rownames(x);write.table(x,stdout(),sep="\\t",quote=FALSE,row.names=FALSE,na="NA");close(con)'
    results=[];candidates=[];streams=[]
    for branch,expected in [('CD3',560916),('CD8',235851)]:
        path=original/('GSE176021_'+branch+'_annotations.rds.gz')
        r=subprocess.run([str(R),'--vanilla','-e',rcode,str(path)],capture_output=True,timeout=120)
        assert r.returncode==0 and len(r.stdout)<256*1024**2,'R_METADATA_DECODE_FAILED'
        d=pd.read_csv(io.BytesIO(r.stdout),sep='\t',keep_default_na=False,dtype={'barcode':str,'imid':str,'CellType':str,'author_row_name':str})
        assert len(d)==expected and set(d)=={'nCount_RNA','nFeature_RNA','barcode','imid','CellType','UMAP_1','UMAP_2','author_row_name'}
        label_missing=d.CellType.str.lower().isin(missing)
        audit={'branch':branch,'rows':len(d),'source_sha256':sha(path),'author_row_unique':bool(d.author_row_name.is_unique),'barcode_unique':bool(d.barcode.is_unique),'composite_unique':bool((d.imid+'|'+d.barcode).is_unique),'valid_label_rows':int((~label_missing).sum()),'label_counts':d.CellType.value_counts().to_dict(),'old_raw_nonmissing':empty_counts['cd3_CellType' if branch=='CD3' else 'cd8_CellType']}
        for rk,rv in [('raw_index',raw.index.to_numpy(str)),('original_cell_id',raw.original_cell_id.to_numpy(str)),('barcode',raw.barcode.to_numpy(str)),('imid_barcode',(raw.imid+'|'+raw.barcode).to_numpy(str))]:
            for ak,av in [('author_row_name',d.author_row_name.to_numpy(str)),('barcode',d.barcode.to_numpy(str)),('imid_barcode',(d.imid+'|'+d.barcode).to_numpy(str))]:
                unique=bool(pd.Index(av).is_unique);raw_unique=bool(pd.Index(rv).is_unique)
                if unique:
                    loc=pd.Index(av).get_indexer(rv);good=loc>=0
                    imid_agree=bool(np.array_equal(raw.imid.to_numpy(str)[good],d.imid.to_numpy(str)[loc[good]]))
                    matches=int(good.sum());author_matched=int(np.unique(loc[good]).size)
                else:matches=int(pd.Index(rv).isin(av).sum());author_matched=None;imid_agree=False
                candidates.append({'branch':branch,'raw_key':rk,'author_key':ak,'raw_key_unique':raw_unique,'author_key_unique':unique,'matched_raw_rows':matches,'matched_author_rows':author_matched,'exact_imid_agreement':imid_agree})
        results.append(audit)
        streams.append({'branch':branch,'export_payload_sha256':hashlib.sha256(r.stdout).hexdigest(),'export_bytes':len(r.stdout)})
    pd.DataFrame(candidates).to_csv(out/'private_all_exact_key_candidates.tsv',sep='\t',index=False)
    (out/'private_all_author_schema_counts_and_decode_hashes.json').write_text(json.dumps({'raw_rows':len(raw),'raw_empty_counts':empty_counts,'branches':results,'R_export_payloads':streams},indent=2)+'\n')
    files={p.name:{'sha256':sha(p),'bytes':p.stat().st_size} for p in out.iterdir() if p.is_file()}
    assert sum(v['bytes'] for v in files.values())<8*1024**2
    receipt={'status':'COMPLETE_DOUBLE_GZIP_AUTHOR_METADATA_AND_ALL_EXACT_KEY_CANDIDATES_QUALIFIED_NO_MAPPING_ADOPTED','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'raw_rows':len(raw),'author_rows':[560916,235851],'files':files,'expression_scores_clinical_effects_read':False,'authority_modified':False}
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'terminal_sha256':sha(out/'private_terminal.json')}))


if __name__=='__main__':main()
