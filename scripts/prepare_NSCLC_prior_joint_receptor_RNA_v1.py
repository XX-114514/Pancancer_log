"""Native same-cell RNA for every frozen reference-compatible prior receptor.

Reuses the existing compressed count matrix. No new download or LR effects.
Only multicomponent receptors require another native coordinate read.
"""
import gzip,hashlib,json,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from prepare_NSCLC_full_native_gene_units_v1 import TrackedReader,chunks,sha

def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);P=R/'NSCLC_all_frozen_receiver_prior_ranking_20261016_v1'
    launch=json.loads((R/'HNSCC_receiver_followup_job_control_20261016_v2/audit_NSCLC_all_frozen_receiver_priors_v1/launch_receipt.json').read_text());deadline=time.monotonic()+7200
    while not (P/'private_independent_audit.json').exists():
        s=Path('/proc')/str(launch['pid'])/'stat';assert s.exists();st=s.read_text().rsplit(')',1)[1].split()
        assert st[19]==str(launch['proc_start_ticks']) and st[0]!='Z' and time.monotonic()<deadline;time.sleep(10)
    assert json.loads((P/'private_independent_audit.json').read_text())['status'].startswith('PASS_ALL24_BACKGROUNDS_')
    assert json.loads((R/'NSCLC_full_native_gene_unit_audit_20261016_v1/private_terminal.json').read_text())['status'].startswith('PASS_ALL_OUTPUTS_')
    O=R/'NSCLC_prior_joint_receptor_RNA_20261016_v1';assert not O.exists();O.mkdir(mode=0o700)
    ep=P/'private_all_shortlist_reference_interfaces.tsv';edges=pd.read_csv(ep,sep='\t');selected=edges[edges.reference_and_marginal_compatible]
    receptors=sorted(set(selected.receptor));multi=[r for r in receptors if '_' in r];panel=sorted({g for r in multi for g in r.split('_')})
    N=R/'NSCLC_full_native_gene_units_20261016_v1';native=pd.Index(pd.read_csv(N/'private_native_gene_axis.tsv',sep='\t').gene)
    units=pd.read_csv(N/'private_all_native_units_and_quality.tsv',sep='\t');unit=np.load(N/'private_cell_unit_index.npy',mmap_mode='r');ncells=len(unit)
    assert units.unit_index.tolist()==list(range(len(units))) and np.array_equal(np.bincount(unit,minlength=len(units)),units.n_cells)
    native_columns=native.get_indexer(panel);assert (native_columns>=0).all()
    download=json.loads((R/'count_matrix_download_v1/download_state.json').read_text());matrix=Path(download['path'])
    assert matrix.stat().st_size==download['bytes'];memory=int(next(x for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')).split()[1])*1024
    assert memory>ncells*len(panel)*5+2*1024**3
    audit_cells=np.sort(np.random.default_rng(2026101606).choice(min(4096,ncells),64,replace=False));np.save(O/'private_prespecified_prefix64_cells.npy',audit_cells)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),prior_audit_sha256=sha(P/'private_independent_audit.json'),
        interface_table_sha256=sha(ep),scope='All frozen candidate literal reference-supported receptors whose complete components pass pooled native marginal10% in at least one compatible receiver lineage; no effect selection.',
        no_rank21_replacement=True,native_cells=ncells,native_units=len(units),receptors=len(receptors),multicomponent_receptors=len(multi),panel_genes=len(panel),
        measurement='All component genes positive in the SAME original cell. Native count>=1 main, >=2 each component sensitivity for multicomponent receptors. Singletons reuse independently audited native unit positive counts.',
        reference_gate='Original per-signature reference_and_marginal_compatible remains required. Same-cell pooled1% is a future eligibility gate, not protein assembly or tumor signaling proof.',
        native_axis_sha256=sha(N/'private_native_gene_axis.tsv'),cell_unit_sha256=sha(N/'private_cell_unit_index.npy'),matrix_expected_sha256=download['sha256'],
        audit='All unit/gene marginals, all receptor unit joins and pool summaries; independent pandas C raw coordinate replay of64 predetermined cells within original first4096 cells, not all cells or all author subtypes.',
        no_LR_or_clinical_effects=True)
    (O/'plan_before_native_joint_RNA.json').write_text(json.dumps(plan,indent=2)+'\n')
    pd.DataFrame(dict(gene=panel,native_column=native_columns)).to_csv(O/'private_native_panel_axis.tsv',sep='\t',index=False)
    one=np.zeros((ncells,len(panel)),bool);two=one.copy();lookup=np.full(len(native),-1,int);lookup[native_columns]=np.arange(len(panel));seen=0;previous=-1
    with matrix.open('rb') as raw:
        tracked=TrackedReader(raw)
        with gzip.GzipFile(fileobj=tracked,mode='rb') as stream:
            header=stream.readline();assert header.startswith(b'%%MatrixMarket matrix coordinate')
            line=stream.readline()
            while line.startswith(b'%'):line=stream.readline()
            nr,ng,nnz=map(int,line.split());assert nr==ncells and ng==len(native)
            for block in chunks(stream):
                assert np.isfinite(block).all() and np.array_equal(block,np.floor(block)) and (block[:,2]>0).all()
                row=block[:,0].astype(np.int64)-1;col=block[:,1].astype(np.int64)-1;value=block[:,2]
                key=row*ng+col;assert key[0]>previous and (np.diff(key)>0).all();previous=int(key[-1]);seen+=len(block)
                loc=lookup[col];take=loc>=0;rr=row[take];cc=loc[take];vv=value[take];one[rr,cc]=True;two[rr,cc]=vv>=2
                if seen//50000000 != (seen-len(block))//50000000:
                    (O/'private_progress.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),coordinates_read=seen,total_coordinates=nnz)))
        assert seen==nnz and tracked.digest.hexdigest()==download['sha256'] and tracked.bytes==download['bytes']
    np.save(O/'private_cell_panel_positive_packed.npy',np.packbits(one,axis=1,bitorder='little'))
    np.save(O/'private_cell_panel_positive2_packed.npy',np.packbits(two,axis=1,bitorder='little'))
    counts=np.load(N/'positive_cell_counts.npy',mmap_mode='r');data=[];axis={g:i for i,g in enumerate(panel)}
    for receptor in receptors:
        components=receptor.split('_')
        if len(components)==1:
            j=native.get_loc(components[0]);a=np.asarray(counts[:,j],np.int64);b=None
        else:
            columns=[axis[g] for g in components]
            a=np.bincount(unit,weights=one[:,columns].all(axis=1),minlength=len(units)).astype(np.int64)
            b=np.bincount(unit,weights=two[:,columns].all(axis=1),minlength=len(units)).astype(np.int64)
        data.extend(dict(receptor=receptor,unit_index=int(u.unit_index),joint_positive_cells=int(a[u.unit_index]),joint_positive2_cells=int(b[u.unit_index]) if b is not None else None) for u in units.itertuples())
    table=pd.DataFrame(data);table.to_csv(O/'private_all_native_unit_receptor_joint_counts.tsv.gz',sep='\t',index=False,na_rep='NA',compression='gzip')
    summaries=[]
    for family,prefix in {'Macrophage':'M\u03c6_','CD8':'CD8T_','Treg':'CD4T_Treg_'}.items():
        uu=units[units.sub_cell_type.str.startswith(prefix)];total=int(uu.n_cells.sum());tt=table[table.unit_index.isin(uu.unit_index)]
        for receptor,f in tt.groupby('receptor',sort=True):
            a=int(f.joint_positive_cells.sum());b=int(f.joint_positive2_cells.sum()) if f.joint_positive2_cells.notna().all() else None
            summaries.append(dict(source_family=family,receptor=receptor,native_cells=total,joint_positive_cells=a,joint_positive_fraction=a/total,pooled_joint1pct_gate=a*100>=total,
                joint_positive2_cells=b,joint_positive2_fraction=b/total if b is not None else None,same_cell_RNA_not_protein_assembly=True))
    pd.DataFrame(summaries).to_csv(O/'private_all_lineage_receptor_joint_coverage.tsv',sep='\t',index=False,na_rep='NA')
    files=['private_native_panel_axis.tsv','private_prespecified_prefix64_cells.npy','private_cell_panel_positive_packed.npy','private_cell_panel_positive2_packed.npy','private_all_native_unit_receptor_joint_counts.tsv.gz','private_all_lineage_receptor_joint_coverage.tsv']
    result=dict(status='ALL_FROZEN_PRIOR_RECEPTOR_NATIVE_JOINT_RNA_COMPLETE_INDEPENDENT_AUDIT_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),
        source_coordinates=seen,source_compressed_sha256=tracked.digest.hexdigest(),native_cells=ncells,receptors=len(receptors),multicomponent_receptors=len(multi),
        files={f:sha(O/f) for f in files},new_LR_or_clinical_effects=False,protein_or_causal_axis_validated=False)
    (O/'private_terminal.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))

if __name__=='__main__':main()
