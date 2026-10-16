"""All-unit joint RNA audit and independent raw replay of a fixed cell prefix.

This is RNA compatibility, never protein-complex or causal signaling proof.
"""
import gzip,hashlib,json,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()

def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);O=R/'NSCLC_prior_joint_receptor_RNA_20261016_v1'
    launch=json.loads((R/'HNSCC_receiver_followup_job_control_20261016_v2/prepare_NSCLC_prior_joint_receptor_RNA_v1/launch_receipt.json').read_text());deadline=time.monotonic()+14400
    while not (O/'private_terminal.json').exists():
        p=Path('/proc')/str(launch['pid'])/'stat';assert p.exists();st=p.read_text().rsplit(')',1)[1].split()
        assert st[19]==str(launch['proc_start_ticks']) and st[0]!='Z' and time.monotonic()<deadline;time.sleep(15)
    output=O/'private_independent_audit.json';assert not output.exists()
    terminal=json.loads((O/'private_terminal.json').read_text());plan=json.loads((O/'plan_before_native_joint_RNA.json').read_text())
    for f,h in terminal['files'].items():assert sha(O/f)==h
    N=R/'NSCLC_full_native_gene_units_20261016_v1';native=pd.Index(pd.read_csv(N/'private_native_gene_axis.tsv',sep='\t').gene)
    units=pd.read_csv(N/'private_all_native_units_and_quality.tsv',sep='\t');unit=np.load(N/'private_cell_unit_index.npy',mmap_mode='r');n=len(unit)
    assert sha(N/'private_native_gene_axis.tsv')==plan['native_axis_sha256'] and sha(N/'private_cell_unit_index.npy')==plan['cell_unit_sha256']
    assert np.array_equal(np.bincount(unit,minlength=len(units)),units.n_cells)
    axis=pd.read_csv(O/'private_native_panel_axis.tsv',sep='\t');assert axis.gene.is_unique and axis.native_column.tolist()==native.get_indexer(axis.gene).tolist()
    a=np.load(O/'private_cell_panel_positive_packed.npy');b=np.load(O/'private_cell_panel_positive2_packed.npy')
    assert a.shape==b.shape==(n,(len(axis)+7)//8) and a.dtype==b.dtype==np.uint8
    a=np.unpackbits(a,axis=1,bitorder='little',count=len(axis)).astype(bool);b=np.unpackbits(b,axis=1,bitorder='little',count=len(axis)).astype(bool);assert (b<=a).all()
    positive=np.load(N/'positive_cell_counts.npy',mmap_mode='r');positions={}
    for j,row in enumerate(axis.itertuples()):
        ix=np.flatnonzero(a[:,j]);positions[row.gene]=set(ix.tolist())
        counts=pd.Series(unit[ix]).value_counts().reindex(range(len(units)),fill_value=0).to_numpy(np.int64)
        assert np.array_equal(counts,positive[:,int(row.native_column)])
    # Direct raw coordinates for a predetermined, outcome-free bounded sample.
    cells=np.load(O/'private_prespecified_prefix64_cells.npy');assert cells.tolist()==np.sort(np.random.default_rng(2026101606).choice(min(4096,n),64,replace=False)).tolist()
    lookup={int(c):j for j,c in enumerate(cells)};column={int(c):j for j,c in enumerate(axis.native_column)}
    expected=np.zeros((len(cells),len(axis)),np.int64)
    download=json.loads((R/'count_matrix_download_v1/download_state.json').read_text());assert terminal['source_compressed_sha256']==download['sha256']==plan['matrix_expected_sha256']
    with gzip.open(download['path'],'rt') as stream:
        assert stream.readline().startswith('%%MatrixMarket matrix coordinate');line=stream.readline()
        while line.startswith('%'):line=stream.readline()
        nr,ng,nnz=map(int,line.split());assert nr==n and ng==len(native)
        reader=pd.read_csv(stream,sep=r'\s+',header=None,names=['row','col','count'],dtype=np.int64,engine='c',chunksize=1000000)
        for f in reader:
            selected=f[f.row.isin(cells+1)&f.col.isin(axis.native_column.to_numpy()+1)]
            for r in selected.itertuples(index=False):expected[lookup[int(r.row)-1],column[int(r.col)-1]]=int(r.count)
            if int(f.row.iloc[-1])>int(cells.max())+1:break
        reader.close()
    assert np.array_equal(expected>0,a[cells]) and np.array_equal(expected>=2,b[cells])
    table=pd.read_csv(O/'private_all_native_unit_receptor_joint_counts.tsv.gz',sep='\t');assert not table.duplicated(['receptor','unit_index']).any()
    receptor_names=sorted(table.receptor.unique());assert len(table)==len(receptor_names)*len(units)
    second={g:set(np.flatnonzero(b[:,j]).tolist()) for j,g in enumerate(axis.gene)}
    for receptor in receptor_names:
        f=table[table.receptor.eq(receptor)].set_index('unit_index').loc[range(len(units))];parts=receptor.split('_')
        if len(parts)==1:
            assert np.array_equal(f.joint_positive_cells,positive[:,native.get_loc(parts[0])]) and f.joint_positive2_cells.isna().all()
        else:
            ix=np.fromiter(set.intersection(*(positions[g] for g in parts)),np.int64)
            ix2=np.fromiter(set.intersection(*(second[g] for g in parts)),np.int64)
            for rows,col in [(ix,'joint_positive_cells'),(ix2,'joint_positive2_cells')]:
                counts=pd.Series(unit[rows]).value_counts().reindex(range(len(units)),fill_value=0).to_numpy(np.int64)
                assert np.array_equal(counts,f[col])
        assert (f.joint_positive_cells.to_numpy()<=units.n_cells.to_numpy()).all()
    summary=pd.read_csv(O/'private_all_lineage_receptor_joint_coverage.tsv',sep='\t');assert len(summary)==3*len(receptor_names)
    for family,prefix in {'Macrophage':'M\u03c6_','CD8':'CD8T_','Treg':'CD4T_Treg_'}.items():
        us=units[units.sub_cell_type.str.startswith(prefix)];total=int(us.n_cells.sum());f=table[table.unit_index.isin(us.unit_index)]
        for row in summary[summary.source_family.eq(family)].itertuples():
            s=f[f.receptor.eq(row.receptor)];one=int(s.joint_positive_cells.sum())
            assert row.native_cells==total and row.joint_positive_cells==one and bool(row.pooled_joint1pct_gate)==(one*100>=total)
            assert np.isclose(row.joint_positive_fraction,one/total)
            if '_' in row.receptor:
                two=int(s.joint_positive2_cells.sum());assert row.joint_positive2_cells==two and np.isclose(row.joint_positive2_fraction,two/total)
            else:assert pd.isna(row.joint_positive2_cells) and pd.isna(row.joint_positive2_fraction)
    P=R/'NSCLC_all_frozen_receiver_prior_ranking_20261016_v1';ep=P/'private_all_shortlist_reference_interfaces.tsv';assert sha(ep)==plan['interface_table_sha256']
    interfaces=pd.read_csv(ep,sep='\t');cover=summary.set_index(['source_family','receptor']);qualified=[]
    for row in interfaces.itertuples(index=False):
        d=row._asdict();old=bool(row.reference_and_marginal_compatible);key=(row.source_family,row.receptor)
        if old:
            c=cover.loc[key];gate=bool(c.pooled_joint1pct_gate)
            d.update(native_pooled_joint_positive_fraction=float(c.joint_positive_fraction),native_joint1pct_qualified=gate,status='REFERENCE_AND_NATIVE_JOINT_RNA_COMPATIBLE' if gate else 'NATIVE_SAME_CELL_RNA_COVERAGE_INSUFFICIENT')
        else:d.update(native_pooled_joint_positive_fraction=None,native_joint1pct_qualified=False,status='ORIGINAL_REFERENCE_OR_MARGINAL_GATE_FAILED')
        d['protein_complex_assembly_verified']=False;qualified.append(d)
    result_path=O/'private_all_frozen_interfaces_joint_RNA_qualification.tsv';pd.DataFrame(qualified).to_csv(result_path,sep='\t',index=False,na_rep='NA')
    result=dict(status='PASS_ALL_UNIT_GENE_MARGINALS_ALL_RECEPTOR_SET_INTERSECTIONS_ALL_POOLS_AND_FIXED64_RAW_CELL_REPLAY',
        finished_UTC=datetime.now(timezone.utc).isoformat(),audit_code_sha256=sha(__file__),producer_terminal_sha256=sha(O/'private_terminal.json'),
        raw_replay_cells=len(cells),raw_replay_sample_is_prefix_bounded=True,all_cells_independently_raw_replayed=False,
        panel_genes=len(axis),native_units=len(units),receptors=len(receptor_names),qualification_sha256=sha(result_path),
        new_LR_or_clinical_tests=False,protein_complex_or_causal_axis_validated=False)
    output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))

if __name__=='__main__':main()
