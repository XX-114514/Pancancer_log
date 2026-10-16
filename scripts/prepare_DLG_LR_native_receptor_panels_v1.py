"""Native RNA/lineage checks for every original corrected-significant LR row.

No new effect tests or cell relabeling. Receptor RNA is not surface protein.
"""
import hashlib,json,os
from collections import defaultdict
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd,zarr
from scipy.sparse import csr_matrix

MARKERS={
    'T_lineage':['TRAC','CD3D','CD3E'],
    'Treg_related':['FOXP3','IL2RA','CTLA4'],
    'myeloid_RNA':['LST1','LYZ','TYROBP'],
    'B_lineage_RNA':['MS4A1','CD79A','CD79B'],
    'epithelial_RNA':['EPCAM','KRT8','KRT18'],
    'endothelial_RNA':['PECAM1','VWF','EMCN'],
}


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(16777216),b''):h.update(block)
    return h.hexdigest()


def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);P=Path(os.environ['PANCANCER_PROJECT_ROOT'])
    V=Path(json.loads((R/'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root'])
    S=V/'state_after_DIALOGUE_LR_response_20261016_v1'
    original_path=S/'all_estimable_exploratory_LR_receiver_associations.tsv'
    t=pd.read_csv(original_path,sep='\t');selected=t[t.BH_q.lt(.05)].sort_values(['context','source','target','MCP','ligand','receptor'])
    assert len(selected)==14
    F=R/'all_CM_state_native_input_preflight_20261014_v1';fi=json.loads((F/'inventory.json').read_text())
    unit_path=F/'private_all_frozen_sample_node_units.tsv.gz';assert sha(unit_path)==fi['private_unit_table_sha256']
    row_path=F/'private_canonical_row_selection.npz';assert sha(row_path)==fi['private_row_selection_sha256']
    u=pd.read_csv(unit_path,sep='\t',keep_default_na=False);assert u.unit_index.is_unique
    u=u.set_index('unit_index',drop=False)
    out=R/'DLG_LR_native_receptor_panels_20261016_v1';assert not out.exists();out.mkdir(mode=0o700)
    genes=sorted({g for r in selected.itertuples() for value in [r.ligand,r.receptor] for g in value.split('_')}|{g for values in MARKERS.values() for g in values}|{'CD44'})
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),
        original_result_sha256=sha(original_path),scope='All fourteen original BH-significant records, no substitution or new inferential p family.',
        populations='Every original native cell in both endpoints and every originally selected context library; full units, not the capped MCP cell draws.',
        markers=MARKERS,lineage_rules='At least2 of3 T genes and at least1 of3 Treg-related genes; other marker sets at least2of3. RNA evidence only, not formal cell identity or doublet calls.',
        receptor='All literal receptor-chain genes in the SAME native cell, report >=1 and >=2 counts; extra CD44 joint RNA descriptive for MIF only, not a new selected axis.',
        detection='Per original case/sample block and pooled native-cell counts; zeros retained, absent or ambiguous native genes missing, no expression-based gene substitutions.',
        expression='All measured panel genes: full-native-RNA mean per-cell log1pCP10K and mean CP10K, for later sender/receiver decomposition; no effect computed here.',
        interpretation='RNA cannot determine surface localization, complex assembly, binding or causality. Myeloid/B RNA co-occurrence flags sensitivity, not automatic contamination or cell exclusion.',
        original_labels_changed=False,original_p_BH_changed=False,new_expression_download=False)
    (out/'plan_before_native_RNA_checks.json').write_text(json.dumps(plan,indent=2)+'\n')
    with np.load(row_path,allow_pickle=False) as z:global_rows=z['global_rows'];unit_ids=z['unit_index']
    assert len(global_rows)==4532933 and len(global_rows)==len(unit_ids)
    endpoints={};needed=set()
    for context in selected.context.unique():
        path=V/'state_after_DIALOGUE_inputs_20261016_v3'/f'context_{int(context):04d}'
        ci=json.loads((path/'private_input_inventory.json').read_text())
        part=pd.read_csv(path/'private_block_library_partition.tsv',sep='\t')
        endpoints[int(context)]=(ci['endpoints'],part)
        needed.update(part.unit_index_A.astype(int));needed.update(part.unit_index_B.astype(int))
    assert needed<=set(u.index)
    chosen=u.loc[sorted(needed)]
    AR=P/'Project_v3/annotation_rebuild/20261010_uniform_expression_annotation_v1'
    contract=json.loads((AR/'input_contract.json').read_text());core=zarr.open_group(contract['core_zarr'],mode='r');X=core['X'];core_genes=contract['raw_shape'][1]
    population={};expression_rows=[];source_records=[]
    for source,source_units in chosen.groupby('gse_id',sort=True):
        folder=out/source;folder.mkdir(mode=0o700)
        gt_path=AR/'by_gse'/source/'genes.tsv';gt=pd.read_csv(gt_path,sep='\t')
        expected=pd.read_csv(F/'all_source_native_scope.tsv',sep='\t').set_index('gse_id').loc[source,'native_gene_axis_sha256']
        assert sha(gt_path)==expected
        lookup=defaultdict(set)
        for r in gt.itertuples():lookup[str(r.gene)].add(int(r.core_column))
        mapped={gene:next(iter(lookup[gene])) for gene in genes if len(lookup[gene])==1}
        present=sorted(mapped);columns=np.array([mapped[g] for g in present],int)
        allowed=np.zeros(core_genes,bool);allowed[gt.core_column.to_numpy(int)]=True
        positions=np.flatnonzero(np.isin(unit_ids,source_units.unit_index.to_numpy(int)))
        order=np.argsort(global_rows[positions]);positions=positions[order];rr=global_rows[positions]
        ids=unit_ids[positions].astype(int);assert len(set(rr))==len(rr)
        assert len(rr)==int(source_units.n_cells.sum())
        panel=np.zeros((len(rr),len(present)),np.int32);RNA=np.zeros(len(rr),np.int64);detected=np.zeros(len(rr),np.int32)
        start=0
        while start<len(rr):
            end=min(start+500,len(rr),int(np.searchsorted(rr,rr[start]+5000,side='left')));assert end>start
            take=rr[start:end];lo=int(take[0]);hi=int(take[-1])+1;ip=np.asarray(X['indptr'][lo:hi+1],np.int64);a,b=int(ip[0]),int(ip[-1])
            raw=csr_matrix((X['data'][a:b],X['indices'][a:b],ip-a),shape=(hi-lo,core_genes))[take-lo].tocsr()
            raw.sum_duplicates();raw.eliminate_zeros();assert (raw.data>0).all() and np.isfinite(raw.data).all()
            assert np.array_equal(raw.data,np.floor(raw.data)) and allowed[raw.indices].all()
            raw_values=raw[:,columns].toarray();assert raw_values.max(initial=0)<np.iinfo(np.int32).max
            panel[start:end]=raw_values.astype(np.int32);RNA[start:end]=np.asarray(raw.sum(axis=1)).ravel().astype(np.int64)
            detected[start:end]=np.diff(raw.indptr);start=end
        assert (RNA>0).all()
        np.save(folder/'private_native_panel_counts.npy',panel,allow_pickle=False)
        np.save(folder/'private_full_RNA.npy',RNA,allow_pickle=False)
        np.save(folder/'private_detected_genes.npy',detected,allow_pickle=False)
        md=pd.DataFrame(dict(panel_row=np.arange(len(rr)),global_row=rr,unit_index=ids))
        md=md.merge(source_units[['unit_index','sample_uid','effective_resampling_block','subCluster','major']],on='unit_index',how='left',validate='many_to_one').sort_values('panel_row')
        md.to_csv(folder/'private_all_native_cell_scope.tsv.gz',sep='\t',index=False)
        pd.DataFrame([dict(gene=g,source_measured_unique=g in mapped,core_column=mapped.get(g)) for g in genes]).to_csv(folder/'private_panel_gene_measurement.tsv',sep='\t',index=False,na_rep='NA')
        pd.Series(present).to_csv(folder/'private_panel_genes.tsv',sep='\t',index=False,header=False)
        gi={gene:i for i,gene in enumerate(present)}
        for uid in source_units.unit_index:
            use=ids==uid;count=int(use.sum());assert count==int(u.loc[uid,'n_cells'])
            raw=panel[use];norm=raw*10000/RNA[use,None]
            flags={}
            for name,markers in MARKERS.items():
                if set(markers)<=set(present):flags[name]=(raw[:,[gi[g] for g in markers]]>0).sum(axis=1)>=(1 if name=='Treg_related' else 2)
            flags['T_lineage_and_Treg_RNA']=flags['T_lineage']&flags['Treg_related'] if 'T_lineage' in flags and 'Treg_related' in flags else None
            population[int(uid)]=dict(source=source,raw=raw,genes=gi,RNA=RNA[use],detected=detected[use],flags=flags)
            for gene in genes:
                exists=gene in gi;j=gi.get(gene)
                expression_rows.append(dict(unit_index=int(uid),gse_id=source,sample_uid=u.loc[uid,'sample_uid'],endpoint=u.loc[uid,'subCluster'],gene=gene,n_cells=count,
                    measured_unique=exists,positive_cells=int((raw[:,j]>0).sum()) if exists else np.nan,
                    mean_CP10K=float(norm[:,j].mean()) if exists else np.nan,mean_logCP10K=float(np.log1p(norm[:,j]).mean()) if exists else np.nan,
                    mean_logRNA=float(np.log1p(RNA[use]).mean()),mean_loggenes=float(np.log1p(detected[use]).mean())))
        source_records.append(dict(gse_id=source,native_cells=len(rr),native_units=len(source_units),measured_panel_genes=len(present),
            panel_sha256=sha(folder/'private_native_panel_counts.npy'),gene_axis_sha256=sha(gt_path)))
        (out/'private_progress.json').write_text(json.dumps(dict(completed_source=source,observed_UTC=datetime.now(timezone.utc).isoformat())))
    detail=[]
    for r in selected.itertuples():
        pair,part=endpoints[int(r.context)];reverse=r.source==pair[1];assert r.source in pair and r.target in pair
        source_col='unit_index_B' if reverse else 'unit_index_A';target_col='unit_index_A' if reverse else 'unit_index_B'
        for row in part.itertuples():
            suid=int(getattr(row,source_col));tuid=int(getattr(row,target_col));s=population[suid];tgt=population[tuid]
            ligands=r.ligand.split('_');receptors=r.receptor.split('_')
            sok=all(g in s['genes'] for g in ligands);rok=all(g in tgt['genes'] for g in receptors)
            sp=np.all(s['raw'][:,[s['genes'][g] for g in ligands]]>=1,axis=1) if sok else None
            rp=np.all(tgt['raw'][:,[tgt['genes'][g] for g in receptors]]>=1,axis=1) if rok else None
            rp2=np.all(tgt['raw'][:,[tgt['genes'][g] for g in receptors]]>=2,axis=1) if rok else None
            rec=dict(context=int(r.context),gse_id=r.gse_id,source=r.source,target=r.target,MCP=r.MCP,ligand=r.ligand,receptor=r.receptor,
                block=row.effective_resampling_block,partition=row.partition,sample_uid=row.sample_uid,source_unit=suid,target_unit=tuid,
                source_cells=len(s['raw']),target_cells=len(tgt['raw']),ligand_components_measured=sok,receptor_components_measured=rok,
                ligand_positive_cells=int(sp.sum()) if sok else np.nan,receptor_joint1_cells=int(rp.sum()) if rok else np.nan,
                receptor_joint2_cells=int(rp2.sum()) if rok else np.nan,receptor_joint1_fraction=float(rp.mean()) if rok else np.nan,
                receptor_joint2_fraction=float(rp2.mean()) if rok else np.nan,
                receptor_joint_CD44_cells=int((rp&(tgt['raw'][:,tgt['genes']['CD44']]>0)).sum()) if rok and 'CD44' in tgt['genes'] and r.ligand=='MIF' else np.nan)
            for name,flag in tgt['flags'].items():
                rec['target_'+name+'_cells']=int(flag.sum()) if flag is not None else np.nan
                rec['receptor_joint1_and_'+name+'_cells']=int((rp&flag).sum()) if rok and flag is not None else np.nan
            detail.append(rec)
    detail=pd.DataFrame(detail);summary=[]
    for key,f in detail.groupby(['context','source','target','MCP','ligand','receptor'],sort=True):
        value=dict(zip(['context','source','target','MCP','ligand','receptor'],key));available=f.receptor_components_measured.all()
        value.update(matched_blocks=len(f),native_source_cells=int(f.source_cells.sum()),native_target_cells=int(f.target_cells.sum()),
            receptor_components_all_measured=available,blocks_with_any_joint_RNA=int(f.receptor_joint1_cells.gt(0).sum()) if available else np.nan,
            blocks_with_atleast1pct_joint_RNA=int((f.receptor_joint1_cells*100>=f.target_cells).sum()) if available else np.nan,
            pooled_joint1_fraction=float(f.receptor_joint1_cells.sum()/f.target_cells.sum()) if available else np.nan,
            pooled_joint2_fraction=float(f.receptor_joint2_cells.sum()/f.target_cells.sum()) if available else np.nan)
        for name in list(MARKERS)+['T_lineage_and_Treg_RNA']:
            column='receptor_joint1_and_'+name+'_cells'
            value['joint1_'+name+'_fraction']=float(f[column].sum()/f.receptor_joint1_cells.sum()) if column in f and available and f.receptor_joint1_cells.sum()>0 else np.nan
        summary.append(value)
    detail.to_csv(out/'private_all14_native_case_receptor_checks.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(summary).to_csv(out/'all14_private_native_receptor_summary.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(expression_rows).to_csv(out/'private_all_native_unit_panel_expression.tsv',sep='\t',index=False,na_rep='NA')
    report=dict(status='ALL14_FROZEN_LR_NATIVE_RNA_AND_LINEAGE_PANELS_COMPLETE_NO_NEW_EFFECT_TESTS',finished_UTC=datetime.now(timezone.utc).isoformat(),
        source_records=source_records,records=len(summary),native_case_records=len(detail),
        summary_sha256=sha(out/'all14_private_native_receptor_summary.tsv'),case_table_sha256=sha(out/'private_all14_native_case_receptor_checks.tsv'),
        source_result_sha256=sha(original_path),current_tumor_surface_complex_verified=False,cell_labels_changed=False)
    (out/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
    print(pd.DataFrame(summary)[pd.DataFrame(summary).ligand.eq('MIF')].to_json(orient='records'))


if __name__=='__main__':main()
