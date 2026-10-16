"""Replay every native panel cell and all fourteen case/summary records.

Use sparse gene-wise nonzero positions for receptor/lineage intersections.
Missing native marker sets remain unavailable rather than biological zeros.
"""
import hashlib,json,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd,zarr
from scipy.sparse import csr_matrix


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(16777216),b''):h.update(b)
    return h.hexdigest()


def eq(a,b):
    assert np.isclose(a,b,atol=1e-10,rtol=1e-10,equal_nan=True),(a,b)


def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);P=Path(os.environ['PANCANCER_PROJECT_ROOT'])
    N=R/'DLG_LR_native_receptor_panels_20261016_v3';O=R/'DLG_LR_native_receptor_panel_audit_20261016_v1'
    assert not O.exists();O.mkdir(mode=0o700)
    receipt=R/'HNSCC_receiver_followup_job_control_20261016_v2/prepare_DLG_LR_native_receptor_panels_v3/launch_receipt.json'
    launch=json.loads(receipt.read_text());deadline=time.monotonic()+3600
    while not (N/'private_terminal.json').exists():
        p=Path('/proc')/str(launch['pid'])/'stat';assert p.exists()
        x=p.read_text().rsplit(')',1)[1].split();assert x[19]==str(launch['proc_start_ticks']) and x[0]!='Z'
        assert time.monotonic()<deadline;time.sleep(10)
    terminal=json.loads((N/'private_terminal.json').read_text());plan=json.loads((N/'plan_before_native_RNA_checks.json').read_text())
    assert terminal['status']=='ALL14_FROZEN_LR_NATIVE_RNA_AND_LINEAGE_PANELS_COMPLETE_NO_NEW_EFFECT_TESTS'
    assert sha(N/'all14_private_native_receptor_summary.tsv')==terminal['summary_sha256']
    assert sha(N/'private_all14_native_case_receptor_checks.tsv')==terminal['case_table_sha256']
    (O/'plan_before_independent_panel_audit.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),
        producer_terminal_sha256=sha(N/'private_terminal.json'),scope='Every selected native cell, every panel gene, all original units and all14 candidate case records.',
        independent_method='Canonical CSR row batches and sparse per-gene positions, separate unit aggregation and set intersections; no producer helper imports.'),indent=2)+'\n')
    V=Path(json.loads((R/'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root'])
    original=pd.read_csv(V/'state_after_DIALOGUE_LR_response_20261016_v1/all_estimable_exploratory_LR_receiver_associations.tsv',sep='\t')
    original=original[original.BH_q.lt(.05)];assert len(original)==14
    F=R/'all_CM_state_native_input_preflight_20261014_v1'
    units=pd.read_csv(F/'private_all_frozen_sample_node_units.tsv.gz',sep='\t',keep_default_na=False).set_index('unit_index',drop=False)
    rowmap=np.load(F/'private_canonical_row_selection.npz',allow_pickle=False)
    AR=P/'Project_v3/annotation_rebuild/20261010_uniform_expression_annotation_v1'
    c=json.loads((AR/'input_contract.json').read_text());X=zarr.open_group(c['core_zarr'],mode='r')['X'];ng=c['raw_shape'][1]
    ex=pd.read_csv(N/'private_all_native_unit_panel_expression.tsv',sep='\t').set_index(['unit_index','gene'])
    detail=pd.read_csv(N/'private_all14_native_case_receptor_checks.tsv',sep='\t')
    summary=pd.read_csv(N/'all14_private_native_receptor_summary.tsv',sep='\t')
    keys=['context','source','target','MCP','ligand','receptor'];ss=summary.set_index(keys)
    assert set(ss.index)==set(original.set_index(keys).index) and len(ss)==14
    cache={};allcells=0;maxerr=0.
    for info in terminal['source_records']:
        source=info['gse_id'];d=N/source
        panel=np.load(d/'private_native_panel_counts.npy',mmap_mode='r');RNA=np.load(d/'private_full_RNA.npy',mmap_mode='r');det=np.load(d/'private_detected_genes.npy',mmap_mode='r')
        assert sha(d/'private_native_panel_counts.npy')==info['panel_sha256']
        md=pd.read_csv(d/'private_all_native_cell_scope.tsv.gz',sep='\t');genes=pd.read_csv(d/'private_panel_genes.tsv',header=None)[0].to_numpy(str)
        measurement=pd.read_csv(d/'private_panel_gene_measurement.tsv',sep='\t').set_index('gene')
        gt=pd.read_csv(AR/'by_gse'/source/'genes.tsv',sep='\t');assert sha(AR/'by_gse'/source/'genes.tsv')==info['gene_axis_sha256']
        gmap=gt.groupby('gene').core_column.agg(lambda x:set(x));cols=[]
        for g in genes:assert len(gmap[g])==1;cols.append(next(iter(gmap[g])))
        for r in measurement.itertuples():
            expected=r.Index in gmap and len(gmap[r.Index])==1;assert bool(r.source_measured_unique)==expected
        rr=md.global_row.to_numpy(int);uid=md.unit_index.to_numpy(int)
        assert md.panel_row.tolist()==list(range(len(md))) and len(rr)==len(set(rr)) and np.all(np.diff(rr)>0)
        expected_rows=rowmap['global_rows'][np.isin(rowmap['unit_index'],np.unique(uid))]
        assert np.array_equal(rr,np.sort(expected_rows))
        assert np.array_equal(uid,rowmap['unit_index'][np.searchsorted(rowmap['global_rows'],rr)])
        assert len(rr)==int(units.loc[np.unique(uid),'n_cells'].sum())==info['native_cells']
        for start in range(0,len(rr),768):
            end=min(start+768,len(rr));take=rr[start:end];lo=int(take.min());hi=int(take.max())+1
            ip=np.asarray(X['indptr'][lo:hi+1],np.int64);a,b=int(ip[0]),int(ip[-1])
            raw=csr_matrix((X['data'][a:b],X['indices'][a:b],ip-a),shape=(hi-lo,ng))[take-lo].tocsr();raw.sum_duplicates();raw.eliminate_zeros()
            assert np.array_equal(np.asarray(raw.sum(axis=1)).ravel(),RNA[start:end])
            assert np.array_equal(np.diff(raw.indptr),det[start:end])
            for j,col in enumerate(cols):assert np.array_equal(raw[:,col].toarray().ravel(),panel[start:end,j])
        for unit in np.unique(uid):
            use=np.flatnonzero(uid==unit);a=csr_matrix(panel[use]);G={g:j for j,g in enumerate(genes)}
            qc=RNA[use];dg=det[use];n=len(use);assert n==int(units.loc[unit,'n_cells'])
            norm=a.multiply(10000./qc[:,None]).tocsr();ln=norm.copy();ln.data=np.log1p(ln.data)
            means=np.asarray(norm.sum(axis=0)).ravel()/n;logs=np.asarray(ln.sum(axis=0)).ravel()/n;positive=np.asarray((a>0).sum(axis=0)).ravel()
            for j,g in enumerate(genes):
                r=ex.loc[(unit,g)];assert r.measured_unique and int(r.positive_cells)==positive[j] and r.n_cells==n
                eq(r.mean_CP10K,means[j]);eq(r.mean_logCP10K,logs[j]);eq(r.mean_logRNA,np.log1p(qc).mean());eq(r.mean_loggenes,np.log1p(dg).mean())
                maxerr=max(maxerr,abs(r.mean_CP10K-means[j]),abs(r.mean_logCP10K-logs[j]))
            for g in set(measurement.index)-set(genes):
                r=ex.loc[(unit,g)];assert not r.measured_unique and pd.isna(r.mean_CP10K) and pd.isna(r.positive_cells)
            positives={g:set(a[:,j].nonzero()[0]) for g,j in G.items()}
            positives2={g:set((a[:,j]>=2).nonzero()[0]) for g,j in G.items()}
            flags={}
            for name,markers in plan['markers'].items():
                if set(markers)<=set(genes):
                    hit={i:sum(i in positives[g] for g in markers) for i in set().union(*(positives[g] for g in markers))}
                    flags[name]={i for i,v in hit.items() if v>=(1 if name=='Treg_related' else 2)}
            if 'T_lineage' in flags and 'Treg_related' in flags:flags['T_lineage_and_Treg_RNA']=flags['T_lineage']&flags['Treg_related']
            cache[int(unit)]=dict(n=n,genes=set(genes),positive=positives,positive2=positives2,flags=flags)
        allcells+=len(md)
        (O/'private_progress.json').write_text(json.dumps(dict(source=source,all_native_cells_audited=allcells,observed_UTC=datetime.now(timezone.utc).isoformat())))
    for r in detail.itertuples():
        s=cache[r.source_unit];t=cache[r.target_unit];ls=r.ligand.split('_');rs=r.receptor.split('_')
        assert r.source_cells==s['n'] and r.target_cells==t['n']
        lp=set.intersection(*(s['positive'][g] for g in ls));rp=set.intersection(*(t['positive'][g] for g in rs));rp2=set.intersection(*(t['positive2'][g] for g in rs))
        assert r.ligand_positive_cells==len(lp) and r.receptor_joint1_cells==len(rp) and r.receptor_joint2_cells==len(rp2)
        eq(r.receptor_joint1_fraction,len(rp)/t['n']);eq(r.receptor_joint2_fraction,len(rp2)/t['n'])
        if r.ligand=='MIF':eq(r.receptor_joint_CD44_cells,len(rp&t['positive']['CD44']))
        for name in list(plan['markers'])+['T_lineage_and_Treg_RNA']:
            a=getattr(r,'target_'+name+'_cells',np.nan);b=getattr(r,'receptor_joint1_and_'+name+'_cells',np.nan)
            if name in t['flags']:eq(a,len(t['flags'][name]));eq(b,len(rp&t['flags'][name]))
            else:assert pd.isna(a) and pd.isna(b)
    for key,f in detail.groupby(keys,sort=True):
        r=ss.loc[key];assert len(f)==r.matched_blocks and len(set(f.block))==len(f)
        orig=original.set_index(keys).loc[key];assert int(r.original_inference_observed_blocks)==int(orig.blocks)
        assert f.source_cells.sum()==r.native_source_cells and f.target_cells.sum()==r.native_target_cells
        eq(r.pooled_joint1_fraction,f.receptor_joint1_cells.sum()/f.target_cells.sum());eq(r.pooled_joint2_fraction,f.receptor_joint2_cells.sum()/f.target_cells.sum())
        assert r.blocks_with_any_joint_RNA==f.receptor_joint1_cells.gt(0).sum()
        assert r.blocks_with_atleast1pct_joint_RNA==(f.receptor_joint1_cells*100>=f.target_cells).sum()
        for name in list(plan['markers'])+['T_lineage_and_Treg_RNA']:
            col='receptor_joint1_and_'+name+'_cells';available=col in f and f[col].notna().all()
            assert bool(r['joint1_'+name+'_measured'])==available
            expected=f[col].sum()/f.receptor_joint1_cells.sum() if available and f.receptor_joint1_cells.sum()>0 else np.nan
            eq(r['joint1_'+name+'_fraction'],expected)
    report=dict(status='PASS_ALL_NATIVE_PANEL_CELLS_ALL_UNIT_MEANS_ALL586_CASE_CHECKS_ALL14_SUMMARIES_MISSINGNESS',
        finished_UTC=datetime.now(timezone.utc).isoformat(),producer_terminal_sha256=sha(N/'private_terminal.json'),
        native_cells=allcells,native_units=len(cache),case_records=len(detail),candidate_records=len(ss),maximum_mean_error=maxerr,
        full_native_canonical_panel_replayed=True,missing_lineage_sets_preserved=True,surface_complex_verified=False,new_hypothesis_tests=False)
    (O/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()
