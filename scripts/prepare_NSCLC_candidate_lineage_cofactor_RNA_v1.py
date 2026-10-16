"""All audited primary-screen records: native lineage/cofactor RNA diagnostics.

Public marker references are independent of learned programme gene lists.
No hypothesis retesting, identity relabel, new exclusions or protein claims.
"""
import gzip,json,os,shutil
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from prepare_NSCLC_full_native_gene_units_v1 import TrackedReader,chunks,sha

MARKERS={
    'T_lineage':['CD3D','CD3E','TRAC'],
    'myeloid_lineage':['LST1','TYROBP','AIF1','FCER1G','CD68'],
    'Treg_markers':['FOXP3','IL2RA','CTLA4'],
    'CD8_markers':['CD8A','CD8B'],
    'mesenchymal_markers':['COL1A1','COL1A2','COL3A1','DCN','LUM'],
    'endothelial_markers':['PECAM1','VWF','KDR'],
    'epithelial_markers':['EPCAM','KRT8','KRT18'],
    'platelet_markers':['PF4','PPBP'],
}

def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);S=R/'NSCLC_automatic_ligand_response_results_20261016_v1'
    audit=json.loads((S/'private_independent_audit.json').read_text());assert audit['status'].startswith('PASS_COMPLETE_UNIQUE_FAMILY_')
    source=S/'private_all_unique_ligand_response_results.tsv';assert sha(source)==audit['result_sha256'];table=pd.read_csv(source,sep='\t')
    selected=table[table.BH_q_complete_unique_family.lt(.05)].copy();assert len(selected)>0
    B=R/'NSCLC_automatic_ligand_molecular_review_20261016_v1';rulespath=B/'private_public_literature_component_rules.json';rules=json.loads(rulespath.read_text())['ligand_rules']
    assert set(selected.ligand)<=set(rules)
    O=R/'NSCLC_candidate_lineage_cofactor_RNA_20261016_v1';assert not O.exists();O.mkdir(mode=0o700)
    N=R/'NSCLC_full_native_gene_units_20261016_v1';native=pd.Index(pd.read_csv(N/'private_native_gene_axis.tsv',sep='\t').gene)
    units=pd.read_csv(N/'private_all_native_units_and_quality.tsv',sep='\t');unit=np.load(N/'private_cell_unit_index.npy',mmap_mode='r');ncells=len(unit)
    panel=set(g for group in MARKERS.values() for g in group)
    for row in selected.itertuples():
        panel.add(row.ligand)
        panel.update(g for receptor in row.qualified_receptors.split('|') for g in receptor.split('_'))
        rule=rules[row.ligand]
        panel.update(g for role in ['sender_sets','receiver_sets'] for group in rule[role] for g in group);panel.update(rule['auxiliary_genes'])
    panel=sorted(panel);positions=native.get_indexer(panel);available=positions>=0
    pd.DataFrame(dict(gene=panel,native_column=positions,measured=available)).to_csv(O/'private_native_panel_axis.tsv',sep='\t',index=False)
    bytes_planned=ncells*len(panel);assert shutil.disk_usage(O).free>500*1024**3+bytes_planned*2
    cells=np.sort(np.random.default_rng(2026101613).choice(min(4096,ncells),64,replace=False));np.save(O/'private_prespecified_prefix64_cells.npy',cells)
    download=json.loads((R/'count_matrix_download_v1/download_state.json').read_text());matrix=Path(download['path']);assert matrix.stat().st_size==download['bytes']
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),result_audit_sha256=sha(S/'private_independent_audit.json'),result_sha256=sha(source),
        molecular_rules_sha256=sha(rulespath),selected_records=selected.record.tolist(),scope='ALL primaryBH records from the original complete family, not only biologically favored ligands. All original native cells measured, original matched case-unit sets used for each record.',
        public_marker_reference_sets=MARKERS,expression='Native count>=1 and>=2 capped RNA; a ligand/receptor group requires every component in SAME cell.',
        lineage='Original author labels unchanged. Describe >=2 canonical lineage markers among GOI-positive cells, plus >=2 foreign-type marker coincidences; not an automatic identity/contamination/exclusion decision.',
        cofactors='All literature component groups retained, including diagnostic alternatives NOT tested as replacement hypotheses. RNA co-detection does not prove protein stoichiometry, release, redox, surface localization or signaling.',
        no_learned_programme_gene_lists_read_or_exported=True,no_P_or_BH_or_new_exclusions=True,matrix_expected_sha256=download['sha256'],
        native_axis_sha256=sha(N/'private_native_gene_axis.tsv'),cell_unit_sha256=sha(N/'private_cell_unit_index.npy'))
    (O/'plan_before_native_lineage_cofactor_read.json').write_text(json.dumps(plan,indent=2)+'\n')
    counts=np.zeros((ncells,len(panel)),np.uint8);lookup=np.full(len(native),-1,int);lookup[positions[available]]=np.flatnonzero(available);seen=0;previous=-1
    with matrix.open('rb') as raw:
        tracked=TrackedReader(raw)
        with gzip.GzipFile(fileobj=tracked,mode='rb') as stream:
            assert stream.readline().startswith(b'%%MatrixMarket matrix coordinate');line=stream.readline()
            while line.startswith(b'%'):line=stream.readline()
            nr,ng,nnz=map(int,line.split());assert nr==ncells and ng==len(native)
            for block in chunks(stream):
                assert np.isfinite(block).all() and np.array_equal(block,np.floor(block)) and (block[:,2]>0).all()
                rows=block[:,0].astype(np.int64)-1;cols=block[:,1].astype(np.int64)-1;key=rows*ng+cols
                assert key[0]>previous and (np.diff(key)>0).all();previous=int(key[-1]);seen+=len(block)
                loc=lookup[cols];use=loc>=0;counts[rows[use],loc[use]]=np.minimum(block[use,2],2).astype(np.uint8)
                if seen//50000000!=(seen-len(block))//50000000:
                    (O/'private_progress.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),coordinates=seen,total=nnz)))
        assert seen==nnz and tracked.bytes==download['bytes'] and tracked.digest.hexdigest()==download['sha256']
    np.save(O/'private_native_capped2_panel.npy',counts)
    axis={g:i for i,g in enumerate(panel)};diagnostics=[];summaries=[];cache={};designroot=R/'NSCLC_automatic_ligand_response_preflight_20261016_v1/private_frozen_designs'
    for row in selected.itertuples():
        z=np.load(designroot/(row.record+'.npz'),allow_pickle=False)
        for role in ['sender','receiver']:
            ix=z[role+'_units'].astype(int);n=units.n_cells.to_numpy(int)[ix]
            subtype=getattr(row,role+'_author_subtype');core='myeloid_lineage' if subtype.startswith('M\u03c6_') else 'T_lineage'
            groups=[tuple([row.ligand])] if role=='sender' else [tuple(r.split('_')) for r in row.qualified_receptors.split('|')]
            groups+=list(map(tuple,rules[row.ligand][role+'_sets']));groups=sorted(set(tuple(sorted(g)) for g in groups))
            for group in groups:
                measured=all(g in axis and available[axis[g]] for g in group)
                for threshold in [1,2]:
                    label='|'.join(group);key=(group,threshold,core)
                    if measured and key not in cache:
                        active=(counts[:,[axis[g] for g in group]]>=threshold).all(axis=1)
                        values={'GOI_positive_cells':active}
                        for marker,genes in MARKERS.items():
                            columns=[axis[g] for g in genes if available[axis[g]]]
                            ok=len(columns)>=2
                            flag=(counts[:,columns]>=threshold).sum(axis=1)>=2 if ok else np.zeros(ncells,bool)
                            values[marker+'_co_positive_cells']=active&flag
                        cache[key]={name:np.bincount(unit,weights=v,minlength=len(units)).astype(np.int64) for name,v in values.items()}
                    if measured:
                        measured_values={name:v[ix] for name,v in cache[key].items()};positive=int(measured_values['GOI_positive_cells'].sum())
                        for j,u in enumerate(ix):diagnostics.append(dict(record=row.record,role=role,author_subtype=subtype,group=label,threshold=threshold,unit_index=int(u),native_cells=int(n[j]),measured=True,**{name:int(v[j]) for name,v in measured_values.items()}))
                        sums={name:int(v.sum()) for name,v in measured_values.items()}
                        summaries.append(dict(record=row.record,role=role,author_subtype=subtype,group=label,threshold=threshold,measured=True,native_cells=int(n.sum()),cases=len(ix),positive_cases=int(np.sum(measured_values['GOI_positive_cells']>0)),
                            GOI_positive_fraction=positive/int(n.sum()),core_lineage_marker_set=core,core_lineage_among_GOI_positive=sums[core+'_co_positive_cells']/positive if positive else None,
                            **sums))
                    else:summaries.append(dict(record=row.record,role=role,author_subtype=subtype,group=label,threshold=threshold,measured=False,status='REQUIRED_NATIVE_GENE_NOT_MEASURED',native_cells=int(n.sum()),cases=len(ix)))
    pd.DataFrame(diagnostics).to_csv(O/'private_all_matched_unit_lineage_cofactor_counts.tsv.gz',sep='\t',index=False,compression='gzip')
    pd.DataFrame(summaries).to_csv(O/'private_all_record_lineage_cofactor_summaries.tsv',sep='\t',index=False,na_rep='NA')
    files=['private_native_panel_axis.tsv','private_prespecified_prefix64_cells.npy','private_native_capped2_panel.npy','private_all_matched_unit_lineage_cofactor_counts.tsv.gz','private_all_record_lineage_cofactor_summaries.tsv']
    result=dict(status='ALL_ORIGINAL_PRIMARY_SCREEN_RECORDS_NATIVE_LINEAGE_COFACTOR_DIAGNOSTICS_COMPLETE_AUDIT_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),
        result_sha256=sha(source),native_cells=ncells,source_coordinates=seen,source_compressed_sha256=tracked.digest.hexdigest(),panel_columns=len(panel),records=len(selected),summary_rows=len(summaries),
        files={f:sha(O/f) for f in files},no_new_P_BH_or_identity_changes=True,no_learned_programme_gene_list_export=True,protein_or_receptor_signaling_proved=False)
    (O/'private_terminal.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))

if __name__=='__main__':main()
