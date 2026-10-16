"""Independent native-input gate: all outputs/old overlap and 64 raw units.

The raw replay covers every native gene in 64 outcome-free units, not all units.
Producer full-cell conservation and this independently parsed replay are separate.
"""
import gzip,hashlib,json,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(16777216),b''):h.update(b)
    return h.hexdigest()


def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    I=R/'NSCLC_full_native_gene_units_20261016_v1'
    O=R/'NSCLC_full_native_gene_unit_audit_20261016_v1';assert not O.exists();O.mkdir(mode=0o700)
    save=lambda n,d:(O/n).write_text(json.dumps(d,indent=2)+'\n')
    units=pd.read_csv(R/'native_unit_expression_v1/unit_registry.tsv',sep='\t',keep_default_na=False).sort_values('unit_index')
    selected=[]
    for _,f in units.groupby('sub_cell_type',sort=True):
        a=f.sort_values(['n_cells','unit_index']);selected.append(int(a.iloc[len(a)//2].unit_index))
    selected=list(dict.fromkeys(selected+[int(units.loc[units.n_cells.idxmin(),'unit_index']),int(units.loc[units.n_cells.idxmax(),'unit_index'])]))
    rng=np.random.default_rng(2026101601)
    for x in rng.permutation(units.unit_index.to_numpy(int)):
        if int(x) not in selected:selected.append(int(x))
        if len(selected)==64:break
    selected=np.array(sorted(selected),int);assert len(selected)==64
    assert set(units.set_index('unit_index').loc[selected,'sub_cell_type'])==set(units.sub_cell_type)
    save('plan_before_independent_raw_audit.json',dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),
        selection='One median-size unit per original author state, global smallest/largest, then fixed-seed random unused units to64; no clinical data used.',
        selection_seed=2026101601,units=64,raw_gene_scope='Every native gene in each selected unit; independently parsed full matrix stream with pandas C whitespace parser.',
        outputs='All file checksums and array axes/shapes; every old provided gene in every unit; all unit/cell conservation.',
        restriction='Does not claim independent raw replay of all11701 units or clinical/identity acceptance.'))
    pd.DataFrame(dict(unit_index=selected)).to_csv(O/'private_prespecified_64_units.tsv',sep='\t',index=False)
    receipt=R/'HNSCC_receiver_followup_job_control_20261016_v2/prepare_NSCLC_full_native_gene_units_v1/launch_receipt.json'
    launch=json.loads(receipt.read_text());deadline=time.monotonic()+21600
    while not (I/'private_terminal.json').exists():
        proc=Path('/proc')/str(launch['pid'])/'stat'
        assert proc.exists(),'Producer handle missing without terminal; do not restart blindly'
        st=proc.read_text().rsplit(')',1)[1].split();assert st[19]==str(launch['proc_start_ticks']) and st[0]!='Z'
        assert time.monotonic()<deadline
        save('private_progress.json',dict(phase='WAIT_FOR_VERIFIED_LIVE_NATIVE_PRODUCER',observed_UTC=datetime.now(timezone.utc).isoformat()))
        time.sleep(20)
    terminal=json.loads((I/'private_terminal.json').read_text())
    assert terminal['status']=='ALL_NATIVE_NSCLC_FULL_GENE_UNIT_SUMMARIES_COMPLETE_INDEPENDENT_RAW_AUDIT_PENDING'
    for n,v in terminal['outputs'].items():assert sha(I/n)==v['sha256'],n
    u=pd.read_csv(I/'private_all_native_units_and_quality.tsv',sep='\t',keep_default_na=False)
    g=pd.read_csv(I/'private_native_gene_axis.tsv',sep='\t',keep_default_na=False)
    assert u[['sampleID','sub_cell_type','unit_index','n_cells']].equals(units[['sampleID','sub_cell_type','unit_index','n_cells']])
    assert g.native_column.tolist()==list(range(len(g))) and g.gene.is_unique and g.source_measured.all()
    md=pd.read_csv(R/'GSE243013_NSCLC_immune_scRNA_metadata.csv.gz',usecols=['cellID','sampleID','sub_cell_type','total_counts'],keep_default_na=False)
    assert pd.MultiIndex.from_frame(md[['sampleID','cellID']]).is_unique
    ids=pd.MultiIndex.from_frame(units[['sampleID','sub_cell_type']]).get_indexer(pd.MultiIndex.from_frame(md[['sampleID','sub_cell_type']]))
    nc=units.n_cells.to_numpy(int);nu=len(u);ng=len(g)
    assert np.array_equal(np.bincount(ids,minlength=nu),nc)
    assert np.array_equal(np.load(I/'private_cell_unit_index.npy'),ids)
    RNA=np.load(I/'private_full_RNA.npy',mmap_mode='r');det=np.load(I/'private_detected_genes.npy',mmap_mode='r')
    assert np.array_equal(RNA,md.total_counts.to_numpy(np.int64))
    assert (det>0).all() and (det<=ng).all()
    sums=np.load(I/'sum_raw_counts.npy',mmap_mode='r');pos=np.load(I/'positive_cell_counts.npy',mmap_mode='r')
    cp=np.load(I/'mean_CP10K.npy',mmap_mode='r');logcp=np.load(I/'mean_logCP10K.npy',mmap_mode='r')
    assert all(a.shape==(nu,ng) for a in [sums,pos,cp,logcp])
    assert sums.dtype==np.int64 and pos.dtype==np.int32 and cp.dtype==logcp.dtype==np.float64
    assert np.array_equal(sums.sum(axis=1),np.bincount(ids,weights=RNA,minlength=nu).astype(np.int64))
    assert np.array_equal(pos.sum(axis=1,dtype=np.int64),np.bincount(ids,weights=det,minlength=nu).astype(np.int64))
    assert np.allclose(cp.sum(axis=1),10000,atol=1e-6,rtol=0)
    for a in [sums,pos,cp,logcp]:assert np.isfinite(a).all() and (a>=0).all()
    assert (pos<=nc[:,None]).all()
    for name,v in [('mean_logRNA',np.log1p(RNA)),('mean_loggenes',np.log1p(det))]:
        assert np.allclose(u[name],np.bincount(ids,weights=v,minlength=nu)/nc,atol=1e-12,rtol=1e-12)
    old_errors={}
    with np.load(R/'native_unit_expression_v1/native_unit_expression.npz',allow_pickle=False) as old:
        oc=np.flatnonzero(old['provided']);cols=pd.Index(g.gene).get_indexer(old['genes'][oc]);assert (cols>=0).all()
        for name,array in [('mean_CP10K',cp),('mean_logCP10K',logcp),('detected_fraction',pos)]:
            expected=old[name];maximum=0.
            for start in range(0,len(oc),128):
                c=cols[start:start+128];j=oc[start:start+128];actual=np.asarray(array[:,c],float)
                if name=='detected_fraction':actual/=nc[:,None]
                e=np.max(np.abs(actual-expected[:,j]));maximum=max(maximum,float(e))
                assert np.allclose(actual,expected[:,j],atol=1e-8 if name!='detected_fraction' else 1e-12,rtol=1e-9 if name!='detected_fraction' else 0)
            old_errors[name]=maximum
    lookup=np.full(nu,-1,int);lookup[selected]=np.arange(len(selected));cell_unit=lookup[ids]
    rowsum=np.zeros(len(md),np.int64);rowdet=np.zeros(len(md),np.int32)
    raw=np.zeros((64,ng),np.int64);positive=np.zeros((64,ng),np.int64)
    normalized=np.zeros((64,ng));logged=np.zeros((64,ng));entries=0
    ds=json.loads((R/'count_matrix_download_v1/download_state.json').read_text());matrix=Path(ds['path'])
    assert sha(matrix)==ds['sha256']==terminal['matrix_sha256']
    with gzip.open(matrix,'rt',encoding='ascii') as f:
        first=f.readline();assert first.strip()=='%%MatrixMarket matrix coordinate real general'
        line=f.readline()
        while line.startswith('%'):line=f.readline()
        dims=tuple(map(int,line.split()));assert dims==(len(md),ng,terminal['nonzero_coordinates'])
        reader=pd.read_csv(f,sep=r'\s+',names=['r','c','v'],dtype=np.float64,chunksize=1000000,engine='c')
        for chunk in reader:
            a=chunk.to_numpy();assert np.isfinite(a).all() and np.array_equal(a,np.floor(a))
            r=a[:,0].astype(np.int64)-1;c=a[:,1].astype(np.int64)-1;v=a[:,2].astype(np.int64)
            assert (v>0).all() and r.min()>=0 and r.max()<len(md) and c.min()>=0 and c.max()<ng
            pick=cell_unit[r]>=0;ri=r[pick];ci=c[pick];vi=v[pick];ui=cell_unit[ri]
            if len(ri):
                lo=int(ri.min());hi=int(ri.max())+1
                rowsum[lo:hi]+=np.bincount(ri-lo,weights=vi,minlength=hi-lo).astype(np.int64)
                rowdet[lo:hi]+=np.bincount(ri-lo,minlength=hi-lo).astype(np.int32)
                ix=ui*ng+ci;vv=vi*10000./RNA[ri]
                raw+=np.bincount(ix,weights=vi,minlength=64*ng).reshape(64,ng).astype(np.int64)
                positive+=np.bincount(ix,minlength=64*ng).reshape(64,ng)
                normalized+=np.bincount(ix,weights=vv,minlength=64*ng).reshape(64,ng)
                logged+=np.bincount(ix,weights=np.log1p(vv),minlength=64*ng).reshape(64,ng)
            entries+=len(a)
            if entries%25000000==0:
                save('private_progress.json',dict(phase='INDEPENDENT_RAW_REPLAY_64_UNITS_ALL_GENES',entries=entries,expected_entries=dims[2],observed_UTC=datetime.now(timezone.utc).isoformat()))
    assert entries==dims[2]
    cells=cell_unit>=0;assert np.array_equal(rowsum[cells],RNA[cells]) and np.array_equal(rowdet[cells],det[cells])
    assert np.array_equal(raw,sums[selected]) and np.array_equal(positive,pos[selected])
    normalized/=nc[selected,None];logged/=nc[selected,None]
    assert np.allclose(normalized,cp[selected],atol=1e-8,rtol=1e-10)
    assert np.allclose(logged,logcp[selected],atol=1e-10,rtol=1e-10)
    report=dict(status='PASS_ALL_OUTPUTS_OLD_OVERLAP_AND_PRESPECIFIED64_RAW_UNITS_ALL_NATIVE_GENES',finished_UTC=datetime.now(timezone.utc).isoformat(),
        producer_terminal_sha256=sha(I/'private_terminal.json'),all_output_checksums_verified=True,all_unit_conservation_verified=True,
        old_overlap_genes=len(oc),old_overlap_max_errors=old_errors,independent_raw_units=64,independent_raw_cells=int(cells.sum()),
        independent_raw_genes=ng,all_author_states_covered=True,raw_stream_entries=entries,selected_raw_sums_and_positive_exact=True,
        maximum_CP10K_error=float(np.max(np.abs(normalized-cp[selected]))),maximum_logCP10K_error=float(np.max(np.abs(logged-logcp[selected]))),
        all_units_independently_raw_replayed=False,clinical_effects_tested=False,new_download=False)
    save('private_terminal.json',report);print(json.dumps(report))


if __name__=='__main__':main()
