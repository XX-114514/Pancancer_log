"""ALL archive axes/invariants plus96 prespecified pools, all common raw genes.

Fresh sparse COO/bin-count aggregation differs from production group-matrix
projection. Not an all6760-pool independent full-gene raw replay.
"""
import json
import os
import shutil
import subprocess
import time
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
from scipy.sparse import csr_matrix

from NSCLC_external_lineage_inputs_v2 import load,cell_arms,sha,VIEWS
from external_functional_coverage_inputs_v1 import read_inputs


def quota(root):
    fs=subprocess.check_output(['df','-P',str(root.resolve())],text=True).splitlines()[-1].split()[0]
    text=subprocess.check_output(['quota','-g','-v','-w','USER002'],text=True)
    rows=[s.split() for s in text.splitlines() if s.split() and s.split()[0]==fs];assert len(rows)==1
    free=(int(rows[0][3])-int(rows[0][1].rstrip('*')))*1024
    assert min(free,shutil.disk_usage(root).free)>500_000_000_000+1073741824+1048576


def main():
    # Wait on this exact owned producer, never on an unverified lock file.
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    producer=r/'HNSCC_receiver_followup_job_control_20261016_v2/extract_NSCLC_external_frozen_role_count_inputs_v2'
    receipt=json.loads((producer/'launch_receipt.json').read_text());deadline=time.monotonic()+7200
    pending=r/'NSCLC_native_MFA_external_frozen_role_count_inputs_20261016_v2/private_terminal.json'
    while not pending.exists():
        proc=Path('/proc')/str(receipt['pid'])/'stat'
        assert proc.exists(),'Owned count producer missing without terminal; no restart'
        state=proc.read_text().rsplit(')',1)[1].split()
        assert state[0]!='Z' and state[19]==str(receipt['proc_start_ticks'])
        assert time.monotonic()<deadline,'Audit observation window expired; producer not declared stopped'
        time.sleep(3)
    data=load();root=data['root'];quota(root);lib=data['lib'];nu=len(lib)*4
    out=root/'NSCLC_native_MFA_external_frozen_role_count_inputs_20261016_v2'
    tpath=out/'private_terminal.json';terminal=json.loads(tpath.read_text())
    assert terminal['status'].startswith('COMPLETE_ALL61_SOURCE6760_FROZEN_ROLE_NATIVE_COUNT_INPUTS_')
    for name,item in terminal['files'].items():assert sha(out/name)==item['sha256']
    plan=json.loads((out/'plan_before_all_external_native_count_extraction.json').read_text())
    qpath=out/'private_all6760_native_count_units_and_quality.tsv.gz'
    units=pd.read_csv(qpath,sep='\t',keep_default_na=False);assert len(units)==nu and units.unit_index.tolist()==list(range(nu))
    selection=pd.read_csv(out/'private_prespecified96_raw_count_audit_units.tsv.gz',sep='\t',keep_default_na=False)
    assert len(selection)==96 and selection.unit_index.is_unique
    assert sha(out/'private_prespecified96_raw_count_audit_units.tsv.gz')==plan['predetermined_raw_audit_units_sha256']
    audit_indices=selection.unit_index.to_numpy(int);unit_to_audit=np.full(nu,-1,int);unit_to_audit[audit_indices]=np.arange(96)
    for r in selection.itertuples():
        u=units.iloc[r.unit_index]
        assert (r.arm,r.dataset,r.library_id,getattr(r,'view'),r.n_cells)==(u.arm,u.dataset,u.library_id,u['view'],u.n_cells)
    native=pd.read_csv(out/'private_complete_native_gene_axis.tsv',sep='\t',keep_default_na=False)
    original=pd.read_csv(root/'NSCLC_full_native_gene_units_20261016_v1/private_native_gene_axis.tsv',sep='\t',keep_default_na=False)
    assert native.gene.tolist()==original.gene.tolist() and native.native_column.tolist()==list(range(31831))
    genes=native.gene.to_numpy(str);axes=read_inputs(root)['source_axes'];infos={};source_for_audit=np.full(96,-1,int)
    manifests=pd.read_csv(out/'private_all61_count_archive_manifest.tsv',sep='\t');assert len(manifests)==61
    expected_counts={};fresh_counts={};archive_entries=0
    for row in manifests.itertuples():
        key=(row.arm,row.dataset);axis=axes[key]
        target=[next(iter(axis['lookup'][g])) if g in axis['lookup'] and len(axis['lookup'][g])==1 else None for g in genes]
        multiplicity=Counter(c for c in target if c is not None)
        available=np.array([c is not None and multiplicity[c]==1 for c in target]);nc=np.flatnonzero(available)
        physical=np.array([target[i] for i in nc],int)
        ix=lib.index[lib.arm.eq(row.arm)&lib.dataset.eq(row.dataset)].to_numpy(int)
        allunits=(ix[:,None]*4+np.arange(4)).reshape(-1)
        path=out/row.archive;assert sha(path)==row.sha256
        with np.load(path,allow_pickle=False) as z:
            counts=z['counts'];columns=z['common_native_columns'];globalids=z['global_unit_indices']
        assert counts.dtype==np.int64 and counts.shape==(len(allunits),len(nc)) and (counts>=0).all()
        assert np.array_equal(columns,nc) and np.array_equal(globalids,allunits)
        assert np.array_equal(counts.sum(axis=1),units.iloc[allunits].common_native_RNA_sum.to_numpy())
        audit_rows=np.flatnonzero(np.isin(audit_indices,allunits));source_for_audit[audit_rows]=row.source_index
        mapglobal={u:i for i,u in enumerate(allunits)}
        expected_counts[row.source_index]=counts[[mapglobal[u] for u in audit_indices[audit_rows]]]
        fresh_counts[row.source_index]=np.zeros((len(audit_rows),len(nc)),np.int64)
        inv=np.full(96,-1,int);inv[audit_rows]=np.arange(len(audit_rows))
        infos[row.source_index]=dict(arm=row.arm,dataset=row.dataset,physical=physical,native=nc,audit_rows=audit_rows,inverse=inv)
        archive_entries+=counts.size
    assert (source_for_audit>=0).all() and len(set(zip(selection.arm,selection.dataset)))==61
    # Independently compare complete fresh-marker control archive to prior audit.
    control=root/'NSCLC_native_MFA_external_role_marker_profiles_20261016_v2'
    assert sha(control/'private_independent_all_role_marker_audit.json')=='b40688206ea0abfa9ffd0bb795e85e871d59c745db13646416a3831f32d3c1b4'
    baseline=pd.read_csv(control/'private_all6760_frozen_role_units.tsv',sep='\t',keep_default_na=False)
    for field in baseline.columns:assert np.array_equal(units[field],baseline[field]),field
    with np.load(root/'NSCLC_native_MFA_external_lineage_qualification_20261016_v2/private_all1690_three_channel_role_counts.npz',allow_pickle=False) as z:
        parents=z['parent_counts'].reshape(-1)
    assert np.array_equal(parents,units.diagnostic_parent_cells)
    fractions=np.divide(units.n_cells.to_numpy(),parents,out=np.full(nu,np.nan),where=parents>0)
    assert np.allclose(pd.to_numeric(units.diagnostic_fraction_within_parent,errors='coerce'),fractions,equal_nan=True)
    with np.load(control/'private_fixed_marker_pool_sums.npz',allow_pickle=False) as z:
        expected_n=z['n_cells'];expected_log=z['sum_logCP10K'];expected_positive=z['positive_counts']
    assert np.array_equal(units.n_cells,expected_n)
    with np.load(out/'private_all44_fresh_raw_marker_controls.npz',allow_pickle=False) as z:
        assert np.array_equal(z['n_cells'],expected_n) and np.array_equal(z['positive_counts'],expected_positive)
        assert np.allclose(z['sum_logCP10K'],expected_log,atol=1e-6,rtol=1e-9)
    quality=np.zeros((96,12));oldbase=root.parents[2]/'Project_v3/annotation_rebuild/20261010_uniform_expression_annotation_v1'
    checkpoint=root/'HNSCC_receiver_followup_job_control_20261016_v2'/Path(__file__).stem/'private_audit_checkpoint.json'
    for arm in cell_arms(data):
        original_arm=arm['arm']=='frozen_original_annotation'
        z=zarr.open_group(str(data['old_zarr'] if original_arm else data['atlas']),mode='r');raw=z['X'] if original_arm else z['layers']['counts']
        nr,nf=raw.attrs['shape'];destination=np.full(nr,-1,np.int16);role=arm['roles'][0];ok=role>=0
        global_units=arm['library_index'][ok]*4+role[ok]
        rr=np.arange(len(role))
        if original_arm:
            with np.load(data['old']/'private_canonical_row_selection.npz',allow_pickle=False) as f:rr=f['global_rows']
        destination[rr[ok]]=unit_to_audit[global_units]
        if not original_arm:
            order=np.asarray(read_elem(z['uns']['source_order'])).astype(str).tolist()
            measured=np.asarray(read_elem(z['varm']['measured_in_source']),bool);symbol=pd.Series(np.asarray(read_elem(z['var']['gene_symbol'])).astype(str))
        for info in infos.values():
            if info['arm']!=arm['arm']:continue
            colmap=np.full(nf,-1,int)
            if original_arm:
                g=pd.read_csv(oldbase/'by_gse'/info['dataset']/'genes.tsv',sep='\t',keep_default_na=False);cols=g.core_column.to_numpy(int)
                colmap[cols[info['physical']]]=np.arange(len(info['native']))
                present=np.zeros(nf,bool);present[cols]=True;mt=np.zeros(nf,bool);rb=np.zeros(nf,bool)
                mt[cols]=g.gene.str.startswith('MT-').to_numpy();rb[cols]=g.gene.str.match(r'^(RPL|RPS)').to_numpy()
            else:
                colmap[info['physical']]=np.arange(len(info['native']));present=measured[:,order.index(info['dataset'])]
                mt=symbol.str.startswith('MT-').to_numpy();rb=symbol.str.match(r'^(RPL|RPS)').to_numpy()
            info.update(colmap=colmap,present=present,mito=mt,ribo=rb)
        for start in range(0,nr,4096):
            stop=min(nr,start+4096);take=destination[start:stop]>=0
            if not take.any():continue
            ptr=np.asarray(raw['indptr'][start:stop+1],np.int64);lo,hi=int(ptr[0]),int(ptr[-1])
            matrix=csr_matrix((np.asarray(raw['data'][lo:hi]),np.asarray(raw['indices'][lo:hi]),ptr-lo),shape=(stop-start,nf))[take].copy()
            matrix.sum_duplicates();matrix.eliminate_zeros();matrix.sort_indices()
            ids=destination[start:stop][take];source=source_for_audit[ids]
            for si in np.unique(source):
                sub=source==si;m=matrix[sub];audit_ids=ids[sub];co=m.tocoo();info=infos[int(si)]
                assert info['present'][co.col].all()
                rna=np.bincount(co.row,weights=co.data,minlength=len(audit_ids));detected=np.bincount(co.row,minlength=len(audit_ids))
                col=info['colmap'][co.col];keep=col>=0
                common=np.bincount(co.row[keep],weights=co.data[keep],minlength=len(audit_ids));common_detected=np.bincount(co.row[keep],minlength=len(audit_ids))
                local=info['inverse'][audit_ids];nc=len(info['native']);packed=local[co.row[keep]]*nc+col[keep]
                values=np.bincount(packed,weights=co.data[keep],minlength=len(info['audit_rows'])*nc)
                assert np.array_equal(values,values.astype(np.int64))
                fresh_counts[int(si)]+=values.reshape(fresh_counts[int(si)].shape).astype(np.int64)
                mito=np.bincount(co.row,weights=co.data*info['mito'][co.col],minlength=len(audit_ids))
                ribo=np.bincount(co.row,weights=co.data*info['ribo'][co.col],minlength=len(audit_ids))
                items=[np.ones(len(audit_ids)),rna,common,np.log1p(rna),np.log1p(detected),np.log1p(common),np.log1p(common_detected),mito,ribo,rna==0,detected,common_detected]
                for j,value in enumerate(items):quality[:,j]+=np.bincount(audit_ids,weights=value,minlength=96)
            if start//4096%32==0 or stop==nr:checkpoint.write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),arm=arm['arm'],rows_read=stop,total_rows=nr,audit_cells_seen=int(quality[:,0].sum())))+'\n')
        for info in infos.values():
            if info['arm']==arm['arm']:
                for field in ['colmap','present','mito','ribo']:info.pop(field)
    genes_checked=0
    for si in infos:
        assert np.array_equal(fresh_counts[si],expected_counts[si]),infos[si]['dataset']
        genes_checked+=fresh_counts[si].size
    qa=units.iloc[audit_indices]
    for column,j in [('n_cells',0),('full_source_RNA_sum',1),('common_native_RNA_sum',2),('zero_full_RNA_cells',9)]:assert np.array_equal(qa[column].to_numpy(),quality[:,j])
    fields=[('mean_log_full_source_RNA',3),('mean_log_full_source_detected_genes',4),('mean_log_common_native_RNA',5),('mean_log_common_native_detected_genes',6),('mean_full_source_detected_genes',10),('mean_common_native_detected_genes',11)]
    for field,j in fields:
        expected=np.divide(quality[:,j],quality[:,0],out=np.full(96,np.nan),where=quality[:,0]>0)
        assert np.allclose(pd.to_numeric(qa[field],errors='coerce'),expected,atol=1e-10,rtol=1e-10,equal_nan=True)
    for field,j in [('symbol_pct_mito_RNA',7),('symbol_pct_ribo_RNA',8)]:
        expected=np.divide(100*quality[:,j],quality[:,1],out=np.full(96,np.nan),where=quality[:,1]>0)
        assert np.allclose(pd.to_numeric(qa[field],errors='coerce'),expected,atol=1e-10,rtol=1e-10,equal_nan=True)
    proof=out/'private_independent_native_count_input_audit.json';assert not proof.exists()
    report=dict(status='PASS_ALL61_ARCHIVES6760_UNIT_AXES_CONTROLS_AND96_PRESPECIFIED_RAW_POOLS_ALL_COMMON_GENES',
        finished_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),producer_terminal_sha256=sha(tpath),
        all_archive_integer_entries_checked=archive_entries,prespecified_pools_raw_checked=96,raw_common_gene_values_checked=genes_checked,
        raw_cells_in96_pools=int(quality[:,0].sum()),all6760_full_gene_raw_independent_replay=False,
        external_scores_clinical_or_identity_promotion=False)
    proof.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],finished_UTC=report['finished_UTC'],audit_sha256=sha(proof))))


if __name__=='__main__':main()
