"""Exact native full-gene counts for frozen cells; no state fitting or promotion.
Private H5ADs retain cell IDs, original nodes, samples and inherited blocks.
"""
from pathlib import Path
from datetime import datetime, timezone
import os, json, hashlib, gc
import numpy as np
import pandas as pd
import anndata as ad
import zarr
from scipy.sparse import csr_matrix, vstack, load_npz


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()


def balanced_positions(obs, cap, seed):
    rng=np.random.default_rng(seed)
    queues=[rng.permutation(g.index.to_numpy()) for _,g in obs.groupby('analysis_block',sort=True)]
    selected=[];offset=0
    while len(selected)<min(cap,len(obs)):
        advanced=False
        for q in queues:
            if offset<len(q) and len(selected)<cap:
                selected.append(int(q[offset]));advanced=True
        if not advanced:break
        offset+=1
    return sorted(selected)


def main():
    started=datetime.now(timezone.utc).isoformat()
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);A=Path(os.environ['PANCANCER_ANNOTATION_ROOT'])
    C=Path(os.environ['PANCANCER_STATE_INPUT_PLAN']);plan=json.loads(C.read_text())
    assert sha(__file__)==plan['producer_sha256']
    U=A/'20261010_uniform_expression_annotation_v1';P=R/'immune_same_name_family_source_review_v1'
    N=P/'mac_T_IFNG_receptor_single_cell_v1';O=R/'state_singlecell_inputs_20261014_v1'
    assert not O.exists()
    for name,h in plan['frozen_input_sha256'].items():
        assert sha(N/name)==h,name
    original_audit=json.loads((N/'independent_audit.json').read_text())
    assert original_audit['cells']==240533
    contract=json.loads((U/'input_contract.json').read_text());assert sha(U/'input_contract.json')==plan['input_contract_sha256']
    g=zarr.open_group(contract['core_zarr'],mode='r');X=g['X'];assert X.attrs['encoding-type']=='csr_matrix'
    z=np.load(N/'private_selected_cell_gene_bits.npz',allow_pickle=False)
    global_rows=z['global_rows'];unit_index=z['unit_index'];bits=z['gene_bits']
    units=pd.read_csv(N/'private_all1832_sample_node_RNA_counts.tsv',sep='\t')
    assert len(global_rows)==240533 and len(units)==1832 and (np.diff(global_rows)>0).all()
    assert np.array_equal(units.unit_index,np.arange(1832))
    block_map={};review_map={}
    for cid in plan['contexts']:
        p=P/'major_receiver_TF_preflight_v3/private_context_arrays'/(cid+'.npz')
        assert sha(p)==plan['context_sha256'][cid]
        q=np.load(p,allow_pickle=False)
        for sample,block,review in zip(q['samples'],q['blocks'],q['identity_reviewed']):
            sample=str(sample);block=str(block)
            assert sample not in block_map or block_map[sample]==block
            block_map[sample]=block;review_map[sample]=bool(review)
    assert set(units.sample_uid).issubset(block_map) and units.sample_uid.nunique()==212
    gene_tables={}
    for sid in plan['sources']:
        root=U/'by_gse'/sid;manifest=json.loads((root/'manifest.json').read_text())
        assert sha(root/'manifest.json')==plan['source_manifest_sha256'][sid]
        for name in ['genes.tsv','units.tsv','positive_counts.npz','sum_log1p_CP10K.npz']:
            assert sha(root/name)==manifest['output_sha256'][name]
        gene_tables[sid]=pd.read_csv(root/'genes.tsv',sep='\t')
        assert gene_tables[sid].gene.is_unique and gene_tables[sid].core_column.is_unique
    common=sorted(set.intersection(*[set(v.gene) for v in gene_tables.values()]));assert len(common)==16625
    O.mkdir();(O/'private_native_sources').mkdir()
    summaries=[];training=[];training_obs=[];trained_rows=[];assembly_values_checked=0
    for source_number,sid in enumerate(plan['sources']):
        cell_mask=units.gse_id.to_numpy()[unit_index]==sid
        rows=global_rows[cell_mask];source_unit=unit_index[cell_mask];gt=gene_tables[sid]
        cols=gt.core_column.to_numpy();gene_names=gt.gene.to_numpy(dtype=str)
        frame=units.iloc[source_unit][['gse_id','sample_uid','node_id','subCluster','major','unit_index']].reset_index(drop=True)
        frame['core_position']=rows;frame['analysis_block']=frame.sample_uid.map(block_map)
        frame['patient_identity_reviewed']=frame.sample_uid.map(review_map)
        frame['state_discovery_source']=sid!=plan['heldout_source']
        index_key=g['obs'].attrs['_index'];frame.index=pd.Index(g['obs'][index_key].get_orthogonal_selection((rows,)).astype(str),name='cell_id')
        assert frame.index.is_unique
        chunks=[]
        for begin in range(0,len(rows),plan['cell_chunk']):
            chosen=rows[begin:begin+plan['cell_chunk']];lo=int(chosen[0]);hi=int(chosen[-1])+1
            pointer=X['indptr'][lo:hi+1].astype(np.int64);a,b=int(pointer[0]),int(pointer[-1])
            full=csr_matrix((X['data'][a:b],X['indices'][a:b],pointer-a),shape=(hi-lo,contract['raw_shape'][1]))
            full=full[chosen-lo];native=full[:,cols].tocsr()
            assert np.isfinite(native.data).all() and (native.data>=0).all() and np.array_equal(native.data,np.floor(native.data))
            assert np.array_equal(np.asarray(native.sum(axis=1)).ravel(),np.asarray(full.sum(axis=1)).ravel())
            assert native.data.max(initial=0)<np.iinfo(np.int32).max
            chunks.append(native.astype(np.int32))
        raw=vstack(chunks,format='csr');del chunks
        total=np.asarray(raw.sum(axis=1)).ravel();det=np.asarray((raw>0).sum(axis=1)).ravel();assert (total>0).all()
        frame['full_RNA_total']=total;frame['full_detected_genes']=det
        cp=raw.astype(np.float64);cp.data*=np.repeat(10000/total,np.diff(cp.indptr));cp.data=np.log1p(cp.data)
        native_units=pd.read_csv(U/'by_gse'/sid/'units.tsv',sep='\t')
        original=pd.read_csv(U/'all_node_number_catalog.tsv',sep='\t')
        native_units['node_id']=original.node_id.to_numpy()[native_units.node_index]
        exact=native_units[native_units.candidate_eligible].set_index(['sample_uid','node_id'])
        assert exact.index.is_unique
        local_units=units[units.gse_id.eq(sid)];lookup={u:i for i,u in enumerate(local_units.unit_index)}
        local_codes=np.asarray([lookup[u] for u in source_unit])
        incidence=csr_matrix((np.ones(len(rows),np.int64),(local_codes,np.arange(len(rows)))),shape=(len(local_units),len(rows)))
        fresh_positive=(incidence@(raw>0).astype(np.int32)).toarray()
        fresh_log=(incidence@cp).toarray();old_rows=[int(exact.loc[(r.sample_uid,r.node_id)].unit_row) for r in local_units.itertuples()]
        expected_positive=load_npz(U/'by_gse'/sid/'positive_counts.npz')[old_rows].toarray()
        expected_log=load_npz(U/'by_gse'/sid/'sum_log1p_CP10K.npz')[old_rows].toarray()
        assert np.array_equal(fresh_positive,expected_positive)
        assert np.allclose(fresh_log,expected_log,rtol=1e-9,atol=1e-8)
        assert np.array_equal(np.bincount(local_codes,minlength=len(local_units)),local_units.n_cells)
        assembly_values_checked+=fresh_positive.size+fresh_log.size
        for bit,gene in enumerate(['IFNG','IFNGR1','IFNGR2']):
            hit=np.asarray(raw[:,np.flatnonzero(gene_names==gene)].sum(axis=1)).ravel()>0
            assert np.array_equal(hit,(bits[cell_mask]&(1<<bit))>0)
        selected=[]
        if sid!=plan['heldout_source']:
            for major_number,major in enumerate(plan['majors']):
                subset=frame.reset_index().loc[frame.major.eq(major).to_numpy()]
                selected.extend(balanced_positions(subset,plan['max_training_cells_per_source_major'],plan['seed']+source_number*10+major_number))
        selected=sorted(selected);frame['state_training_selected']=False
        if selected:
            frame.iloc[selected,frame.columns.get_loc('state_training_selected')]=True
            common_index=pd.Index(gene_names).get_indexer(common);assert (common_index>=0).all()
            training.append(raw[selected][:,common_index]);training_obs.append(frame.iloc[selected].copy());trained_rows.extend(rows[selected].tolist())
        obj=ad.AnnData(X=raw,obs=frame,var=pd.DataFrame(index=pd.Index(gene_names,name='gene')))
        obj.uns['source_gene_axis_preserved']=True;obj.uns['canonical_core_manifest_sha256']=contract['core_manifest_sha256']
        obj.write_h5ad(O/'private_native_sources'/(sid+'.h5ad'),compression='gzip',compression_opts=4)
        summaries.append(dict(gse_id=sid,cells=len(rows),genes=len(gene_names),libraries=frame.sample_uid.nunique(),analysis_blocks=frame.analysis_block.nunique(),T_cells=int(frame.major.eq('T cell').sum()),Macrophage_cells=int(frame.major.eq('Macrophage').sum()),training_cells=len(selected),heldout_source=sid==plan['heldout_source'],all_gene_positive_and_log_aggregate_checks=fresh_positive.size+fresh_log.size))
        print(json.dumps(summaries[-1]),flush=True)
        del obj,raw,cp,fresh_positive,fresh_log,expected_positive,expected_log;gc.collect()
    train=ad.AnnData(X=vstack(training,format='csr'),obs=pd.concat(training_obs),var=pd.DataFrame(index=pd.Index(common,name='gene')))
    assert train.obs_names.is_unique and len(trained_rows)==train.n_obs and not train.obs.gse_id.eq(plan['heldout_source']).any()
    train.uns['normalization_contract']='Use obs full_RNA_total native denominator, not common-axis sum; raw counts for cNMF plus externally prepared full-native TPM'
    train.write_h5ad(O/'private_training_common_counts.h5ad',compression='gzip',compression_opts=4)
    pd.DataFrame(summaries).to_csv(O/'source_scope_summary.tsv',sep='\t',index=False)
    assert sum(s['cells'] for s in summaries)==240533 and sum(s['libraries'] for s in summaries)==424
    inventory=dict(status='COMPLETE_EXACT_SINGLECELL_INPUTS_NOT_STATE_FIT',started_UTC=started,finished_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(C),cells=240533,sample_node_units=1832,source_sample_major_units=424,common_genes=len(common),training_cells=train.n_obs,source_summaries=summaries,all_aggregate_values_checked=assembly_values_checked,heldout_source=plan['heldout_source'],heldout_is_internal_previously_used_source=True,formal_states_assigned=False,scientific_release=False,outputs={str(p.relative_to(O)):sha(p) for p in O.rglob('*') if p.is_file()})
    (O/'inventory.json').write_text(json.dumps(inventory,indent=2)+'\n');print(json.dumps({k:v for k,v in inventory.items() if k not in ['outputs','source_summaries']}),flush=True)


if __name__=='__main__':main()
