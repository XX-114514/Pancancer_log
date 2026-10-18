"""Reconstruct original source coordinates after complete RNA/row-ID replay.

New small sidecars only; legacy matrices/annotations and archives remain unchanged.
"""
import json
import os
import tarfile
from datetime import datetime,timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from anndata._io.specs import read_elem

from diagnose_GSE299886_original_metadata_identity_v1 import readcsv,strings
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def dense_rows(group,start,stop,nf):
    assert str(group.attrs['encoding-type'])=='csr_matrix'
    ptr=np.asarray(group['indptr'][start:stop+1],np.int64);lo,hi=int(ptr[0]),int(ptr[-1])
    return csr_matrix((np.asarray(group['data'][lo:hi]),np.asarray(group['indices'][lo:hi]),ptr-lo),shape=(stop-start,nf)).toarray()


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget=resource(root,64*1024**2)
    diag=root/'GSE299886_original_metadata_identity_diagnostic_20261018_v1'
    terminal=json.loads((diag/'private_terminal.json').read_text())
    for n,item in terminal['files'].items():assert sha(diag/n)==item['sha256']
    joins=pd.read_csv(diag/'private_all_original_ID_join_diagnostics.tsv',sep='\t')
    chosen=joins[joins.basis.eq('EXACT')&joins.complete_unique_bijection].copy();assert len(chosen)==3 and chosen.sample_id.is_unique
    headers=json.loads((diag/'private_original_source_member_headers.json').read_text())
    schemas=json.loads((root/'CM02_existing_targeted_spatial_schema_20261018_v2/private_all_eligible_targeted_input_schema.json').read_text());byid={r['sample_id']:r for r in schemas}
    out=root/'GSE299886_original_coordinate_sidecars_20261018_v1';assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),budget=budget,original_ID_diagnostic_sha256=sha(diag/'private_terminal.json'),
        scope='Every3 original coverage-eligible MERFISH assay, exact string IDs only. Recompute every retained original RNA value in source row order against BOTH existing X/counts; require raw-count IDs equal existing barcode row by row. Then exact unique metadata EntityID join.',
        repair='New original-ID/metadata-row/center_xy/FOV/volume/RNA vectors with current matrix SHA. No original matrix/obs/annotation rewrite, no deletion or new download.',
        scale='Source center_x/y coordinates, physical unit unverified; old um_assumed is not evidence. FOV remains technical, no independent patient/core claim.',
        no_spatial_effect_or_new_state_labels=True)
    (out/'plan_before_original_coordinate_reconstruction.json').write_text(json.dumps(plan,indent=2)+'\n')
    reports=[]
    for number,item in enumerate(chosen.itertuples()):
        h=next(r for r in headers if r['sample_id']==item.sample_id);schema=byid[item.sample_id];path=Path(schema['analysis_file']);assert sha(path)==schema['input_sha256']
        archive=Path(h['archive']);before=dict(bytes=archive.stat().st_size,mtime_ns=archive.stat().st_mtime_ns)
        with h5py.File(path,'r') as saved,tarfile.open(archive,'r:') as tar:
            n,nf=schema['shape'];obs=saved['obs'];barcode=strings(read_elem(obs['barcode']))
            index=saved['var'].attrs['_index'];index=index.decode() if isinstance(index,bytes) else str(index)
            features=strings(read_elem(saved['var'][index])).tolist()
            raw_columns=h['count_header'];columns=features
            if not set(columns)<=set(raw_columns) and 'gene_symbol' in saved['var']:columns=strings(read_elem(saved['var']['gene_symbol'])).tolist()
            assert len(columns)==nf and len(set(columns))==nf and set(columns)<=set(raw_columns)
            meta=readcsv(tar,h['metadata_member'],dtype=str,keep_default_na=False)
            ids=meta[item.metadata_id_column].to_numpy(str);assert len(set(ids))==len(ids)==n and (ids!='').all()
            indexed=meta.set_index(item.metadata_id_column,drop=False);row_index={v:i for i,v in enumerate(ids)}
            RNA=np.zeros(n,np.int64);processed=0;max_error=0.;source_ids=[];RNA_values=0
            f=tar.extractfile(h['count_member']);assert f is not None
            import gzip
            stream=gzip.GzipFile(fileobj=f) if h['count_member'].endswith('.gz') else f
            for chunk in pd.read_csv(stream,dtype={item.count_id_column:str},keep_default_na=False,chunksize=2048):
                stop=processed+len(chunk);cid=chunk[item.count_id_column].to_numpy(str)
                assert np.array_equal(cid,barcode[processed:stop]),'Original source RNA IDs differ from existing matrix row axis'
                values=chunk[columns].to_numpy(float);assert np.isfinite(values).all() and (values>=0).all() and (values==np.floor(values)).all()
                for key in ['X','layers/counts']:
                    assert key in saved
                    observed=dense_rows(saved[key],processed,stop,nf)
                    assert np.array_equal(values,observed),'Original RNA parity failed; do not attach coordinates'
                    max_error=max(max_error,float(np.max(np.abs(values-observed))))
                RNA[processed:stop]=values.sum(axis=1,dtype=np.float64).astype(np.int64);RNA_values+=values.size
                source_ids.extend(cid);processed=stop
            assert processed==n and len(set(source_ids))==n
            source_ids=np.array(source_ids,str);matched=indexed.loc[source_ids]
            xy=matched[['center_x','center_y']].apply(pd.to_numeric,errors='coerce').to_numpy(float)
            assert xy.shape==(n,2) and np.isfinite(xy).all()
            fov=matched['fov'].to_numpy(str);valid_fov=~np.isin(fov,['','nan','NaN','None'])
            categories,codes=np.unique(fov,return_inverse=True)
            volume=pd.to_numeric(matched['volume'],errors='coerce').to_numpy(float) if 'volume' in matched else np.full(n,np.nan)
            metadata_RNA=pd.to_numeric(matched['transcript_count'],errors='coerce').to_numpy(float) if 'transcript_count' in matched else np.full(n,np.nan)
            mdrows=np.array([row_index[v] for v in source_ids],np.int64)
            file=out/('private_source_%02d_original_coordinates.npz'%number)
            np.savez_compressed(file,source_cell_id=source_ids,existing_row_index=np.arange(n),source_metadata_row_index=mdrows,xy_source=xy,FOV_code=codes.astype(np.int32),FOV_categories=categories,FOV_known=valid_fov,source_volume=volume,deposited_panel_total=RNA,source_metadata_transcript_count=metadata_RNA,existing_matrix_sha256=np.array(schema['input_sha256']))
            with np.load(file,allow_pickle=False) as z:
                assert np.array_equal(z['source_cell_id'],barcode) and np.array_equal(z['xy_source'],xy) and np.array_equal(z['deposited_panel_total'],RNA)
            reports.append(dict(sample_id=item.sample_id,output=file.name,rows=n,features=nf,original_retained_RNA_values=RNA_values,X_and_counts_exact_original_RNA_replayed=True,maximum_RNA_difference=max_error,source_and_existing_barcode_row_order_exact=True,source_metadata_unique_bijection=True,all_source_xy_finite=True,unique_xy=int(len(np.unique(xy,axis=0))),FOV_known_rows=int(valid_fov.sum()),FOV_categories=len(categories),source_volume_finite=int(np.isfinite(volume).sum()),metadata_transcript_count_available=int(np.isfinite(metadata_RNA).sum()),metadata_transcript_total_mismatches=int((np.isfinite(metadata_RNA)&(metadata_RNA!=RNA)).sum()),matrix_SHA256=schema['input_sha256'],source_archive_stat=before,source_members=dict(counts=h['count_member'],metadata=h['metadata_member']),physical_coordinate_unit='UNVERIFIED_SOURCE_CENTER_X_Y',patient_core_identity='UNREVIEWED'))
        assert before==dict(bytes=archive.stat().st_size,mtime_ns=archive.stat().st_mtime_ns) and sha(path)==schema['input_sha256']
    (out/'private_all_source_coordinate_reconstruction_receipt.json').write_text(json.dumps(reports,indent=2)+'\n')
    files=[p for p in out.iterdir() if p.is_file()];assert sum(p.stat().st_size for p in files)<64*1024**2
    t=dict(status='COMPLETE_ALL3_ORIGINAL_COORDINATE_SIDECARS_FULL_RNA_AND_ID_REPLAY_NO_SPATIAL_TEST',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),sources=len(reports),RNA_values_replayed=sum(r['original_retained_RNA_values'] for r in reports),files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in files},original_inputs_or_annotations_modified=False,unit_patient_or_cell_state_validation=False)
    (out/'private_terminal.json').write_text(json.dumps(t,indent=2)+'\n');print(json.dumps(dict(status=t['status'],completed_UTC=t['completed_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))

if __name__=='__main__':main()
