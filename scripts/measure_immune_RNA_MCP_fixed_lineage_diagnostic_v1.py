"""Fixed marker RNA diagnosis of all original views, no identity or main-test change."""
import json
import os
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from anndata.io import read_elem
from scipy.sparse import csr_matrix
from scipy.stats import spearmanr

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from test_immune_RNA_MCP_source_coordination_v1 import fit_rank_partial


PANELS={
    'T_core':['CD3D','CD3E','TRAC'],
    'Myeloid_core':['LYZ','LST1','TYROBP','FCER1G','CSF1R','AIF1'],
    'B_core':['CD79A','CD79B','MS4A1','CD19'],
    'epithelial':['EPCAM','KRT8','KRT18','KRT19'],
    'IFN':['STAT1','IRF1','ISG15','IFIT1','IFIT3','MX1'],
    'proliferation':['MKI67','TOP2A','STMN1','TYMS','PCNA'],
    'immediate_early':['FOS','FOSB','JUN','JUNB','JUND','ATF3','EGR1','DUSP1','HSPA1A','HSPA1B','DNAJB1','IER2','IER3','NR4A1','BTG2']}
THRESHOLDS={'T_core':2,'Myeloid_core':3,'B_core':2,'epithelial':2}


def point(x,y,C):
    f,status=fit_rank_partial(x,y,C);return (f['effect'] if f else np.nan),status


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,128*1024**2)
    inputs=root/'immune_RNA_MCP_native_expression_inputs_20261018_v1'
    assert sha(inputs/'private_terminal.json')=='2474355a003d96232b752ad2c77ba92a506f4c68829c7e0adf434045cf03d656'
    for name,item in json.loads((inputs/'private_terminal.json').read_text())['files'].items():assert sha(inputs/name)==item['sha256']
    md=pd.read_csv(inputs/'private_selected_native_cell_depth.tsv',sep='\t').rename(columns={'fit_sample':'sample'})
    assert len(md)==19018 and md.global_row.is_unique and md.global_row.is_monotonic_increasing
    genes=sorted(set(g for panel in PANELS.values() for g in panel));assert len(genes)==43
    out=root/'immune_RNA_MCP_fixed_lineage_diagnostic_20261018_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_all_fixed_marker_RNA_read.json').write_text(json.dumps(dict(
        created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,output_cap_bytes=16*1024**2,
        panels=PANELS,thresholds=THRESHOLDS,scope='All19018 original cells,3 views/3 sources/2 programmes; all fixed lineage and3 state readouts. No outcome or main-effect-based selection.',
        limits='Cross-lineage RNA not a doublet/identity/physical-interaction call. No gene weights or original394 main tests changed. Fixed IFN/cycle genes fit-excluded; stress set overlaps possible.',
        sensitivities='18 original-anchor coordination point descriptions after fixed marker complement; no new P/FDR/CI. Full-view abundance/count controls retained; subset RNA/gene means and fixed-node predicted means recomputed.'),indent=2)+'\n')
    z=zarr.open_group(str(root.parents[1]/'runs/20261013_immunotherapy_atlas_v1/outputs/immune_validation_authority_v1.zarr'),mode='r')
    order=list(np.asarray(read_elem(z['uns/source_order'])).astype(str));measured=np.asarray(read_elem(z['varm/measured_in_source']),bool)
    symbol=np.asarray(read_elem(z['var/gene_symbol'])).astype(str);canonical=np.asarray(read_elem(z['var/canonical_gene_identity'])).astype(str)
    identity=np.asarray(read_elem(z['var/feature_identity_status'])).astype(str);aliases=np.asarray(read_elem(z['var/source_gene_aliases'])).astype(str)
    columns={}
    for dataset in sorted(md.dataset.unique()):
        lookup=defaultdict(set)
        for col in np.flatnonzero(measured[:,order.index(dataset)]):
            if canonical[col].startswith('AMBIGUOUS::') or identity[col]=='ambiguous_source_feature':continue
            names={symbol[col],canonical[col]}
            if identity[col]=='stable_Ensembl_ID':names|=set(aliases[col].split('|'))
            for gene in names&set(genes):lookup[gene].add(int(col))
        columns[dataset]=np.array([next(iter(lookup[g])) if len(lookup[g])==1 else -1 for g in genes],int)
    pd.DataFrame(dict(gene=genes,**columns)).to_csv(out/'private_fixed_marker_source_columns.tsv',sep='\t',index=False)
    raw=np.full((len(md),len(genes)),-1,np.int32);counts=z['layers/counts'];nc,nf=map(int,counts.attrs['shape']);rows=md.global_row.to_numpy(int);seen=0
    for start in range(0,nc,8192):
        stop=min(nc,start+8192);a,b=np.searchsorted(rows,[start,stop])
        if b>a:
            ptr=np.asarray(counts['indptr'][start:stop+1],np.int64)
            x=csr_matrix((np.asarray(counts['data'][int(ptr[0]):int(ptr[-1])]),np.asarray(counts['indices'][int(ptr[0]):int(ptr[-1])]),ptr-ptr[0]),shape=(stop-start,nf))[rows[a:b]-start].copy()
            x.sum_duplicates();x.eliminate_zeros();x.sort_indices()
            assert np.array_equal(np.asarray(x.sum(axis=1)).ravel(),md.fullRNA.iloc[a:b]) and np.array_equal(np.diff(x.indptr),md.full_detected_genes.iloc[a:b])
            for dataset in md.dataset.iloc[a:b].unique():
                take=np.flatnonzero(md.dataset.iloc[a:b].to_numpy()==dataset);valid=np.flatnonzero(columns[dataset]>=0)
                values=x[take][:,columns[dataset][valid]].toarray();assert np.isfinite(values).all() and (values>=0).all() and (values==np.floor(values)).all() and values.max()<2**31
                raw[np.ix_(a+take,valid)]=values.astype(np.int32)
            seen+=b-a
    assert seen==19018
    np.savez_compressed(out/'private_all19018_fixed43_raw_RNA.npz',raw=raw,genes=np.array(genes),global_rows=rows,fullRNA=md.fullRNA.to_numpy(),full_detected_genes=md.full_detected_genes.to_numpy())
    logs=np.log1p(np.maximum(raw,0).astype(float)*10000/md.fullRNA.to_numpy()[:,None]);logs[raw<0]=np.nan
    flags={}
    for panel,names in PANELS.items():
        cols=[genes.index(g) for g in names];valid=(raw[:,cols]>=0).all(axis=1)
        md[panel+'_meanlog']=np.where(valid,logs[:,cols].mean(axis=1),np.nan)
        if panel in THRESHOLDS:
            flags[panel]=np.where(valid,(raw[:,cols]>0).sum(axis=1)>=THRESHOLDS[panel],np.nan)
            md[panel+'_capture']=flags[panel]
    assert md.loc[~md.view.eq('Myeloid'),'T_core_capture'].eq(1).all()
    md['marker_complement']=False
    for view in md.view.unique():
        foreign=['T_core','B_core','epithelial'] if view=='Myeloid' else ['Myeloid_core','B_core','epithelial']
        take=md.view.eq(view);f=md.loc[take,[k+'_capture' for k in foreign]].to_numpy(float)
        md.loc[take,'marker_complement']=np.isfinite(f).all(axis=1)&(f==0).all(axis=1)
    q=[]
    for view in md.view.unique():
        for part in ['discovery','TNBC_holdout']:
            q.append(pd.read_csv(inputs/'private_inputs'/(view+'_'+part+'_metadata.tsv'),sep='\t')[['cell','cellQ','geneQ']])
    md=md.merge(pd.concat(q,ignore_index=True),on='cell',validate='one_to_one')
    old=root/'existing_immune_fine_native_measurement_scope_20261016_v1';units=pd.read_csv(old/'private_all_frozen_T_myeloid_node_units.tsv',sep='\t').set_index('unit_index')
    mapping=np.load(old/'private_original_cell_to_fine_node_unit.npy',mmap_mode='r');md['node_id']=units.node_id.loc[mapping[md.global_row.to_numpy(int)]].to_numpy()
    md.to_csv(out/'private_all_original_cells_fixed_markers.tsv.gz',sep='\t',index=False)
    projection=root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection'
    scores=pd.read_csv(projection/'private_all_cell_score_reconstruction.tsv.gz',sep='\t');scores=scores[scores.fit_id.eq('K2_seed20261018_source_adjusted')]
    cell=scores[['cell','view','program','score']].merge(md,on=['cell','view'],validate='many_to_one');assert len(cell)==38036
    beta=pd.read_csv(root/'immune_RNA_MCP_original_node_decomposition_20261018_v1/private_frozen_discovery_node_coefficients.tsv',sep='\t').set_index(['view','program','node_id']).coefficient
    saved_meta=pd.read_csv(root/'immune_RNA_MCP_source_coordination_20261018_v1/private_all_case_view_controls.tsv',sep='\t').set_index(['sample','view'])
    cases=[]
    for (sample,dataset,view,program),g in cell.groupby(['sample','dataset','view','program'],sort=True):
        clear=g[g.marker_complement];meta=saved_meta.loc[(sample,view)]
        def comp(f):
            b=beta.reindex(pd.MultiIndex.from_arrays([np.repeat(view,len(f)),np.repeat(program,len(f)),f.node_id])).to_numpy()
            return float(b.mean()) if len(b) and np.isfinite(b).all() else np.nan
        row=dict(sample=sample,dataset=dataset,view=view,program=program,selected_cells=len(g),complement_cells=len(clear),score_original=float(g.score.mean()),score_complement=float(clear.score.mean()) if len(clear) else np.nan,
            composition_original=comp(g),composition_complement=comp(clear),fraction=float(meta.fraction),logN=float(meta.logN),cellQ_original=float(g.cellQ.mean()),geneQ_original=float(g.geneQ.mean()),
            cellQ_complement=float(clear.cellQ.mean()) if len(clear) else np.nan,geneQ_complement=float(clear.geneQ.mean()) if len(clear) else np.nan)
        for panel in ['IFN','proliferation','immediate_early']:row[panel]=float(g[panel+'_meanlog'].mean())
        cases.append(row)
    case=pd.DataFrame(cases);assert len(case)==396;case.to_csv(out/'private_all396_original_and_marker_complement_scores.tsv',sep='\t',index=False)
    capture=[]
    for (dataset,view),g in md.groupby(['dataset','view'],sort=True):
        for panel in THRESHOLDS:
            f=g[panel+'_capture'];capture.append(dict(dataset=dataset,view=view,panel=panel,selected_cells=len(g),case_units=g['sample'].nunique(),measurable_cells=int(f.notna().sum()),capture_cells=int(f.eq(1).sum()),descriptive_not_identity=True))
    pd.DataFrame(capture).to_csv(out/'private_all36_source_marker_capture_descriptions.tsv',sep='\t',index=False)
    state=[]
    for (dataset,view,program),g in case.groupby(['dataset','view','program'],sort=True):
        C=g[['fraction','logN','cellQ_original','geneQ_original']].to_numpy(float)
        for panel in ['IFN','proliferation','immediate_early']:
            valid=np.isfinite(g[[panel,'score_original']].to_numpy()).all(axis=1);x=g.loc[valid,panel].to_numpy();y=g.loc[valid,'score_original'].to_numpy()
            effect,status=point(x,y,C[valid]);state.append(dict(dataset=dataset,view=view,program=program,panel=panel,complete_units=len(x),
                raw_rho=float(spearmanr(x,y).statistic) if len(x)>=3 and min(np.ptp(x),np.ptp(y))>0 else np.nan,adjusted_point=effect,status=status,no_P_or_CI=True))
    pd.DataFrame(state).to_csv(out/'private_all54_fixed_state_readout_points.tsv',sep='\t',index=False)
    contract=json.loads((root/'immune_RNA_MCP_postfit_contract_20261018_v1/plan_before_learned_outputs_or_projection.json').read_text());sensitivity=[]
    for slot in contract['primary_source_association']['slots']:
        a=case[case.view.eq(slot['view_A'])&case.program.eq(slot['program'])];b=case[case.view.eq(slot['view_B'])&case.program.eq(slot['program'])]
        joined=a.merge(b,on=['sample','dataset','program'],suffixes=('_A','_B'),validate='one_to_one')
        for dataset,g in joined.groupby('dataset',sort=True):
            old=[c+'_original_'+side for side in ['A','B'] for c in ['cellQ','geneQ','composition']]+[c+'_'+side for side in ['A','B'] for c in ['fraction','logN']]
            new=[c.replace('_original_','_complement_') for c in old]
            cols=['score_original_A','score_original_B','score_complement_A','score_complement_B']+old+new
            valid=(g.complement_cells_A.ge(20)&g.complement_cells_B.ge(20)).to_numpy()&np.isfinite(g[cols].to_numpy(float)).all(axis=1);d=g[valid]
            before,status0=point(d.score_original_A.to_numpy(),d.score_original_B.to_numpy(),d[old].to_numpy(float))
            after,status1=point(d.score_complement_A.to_numpy(),d.score_complement_B.to_numpy(),d[new].to_numpy(float))
            sensitivity.append(dict(dataset=dataset,program=slot['program'],view_A=slot['view_A'],view_B=slot['view_B'],source_units=len(g),same_comparable_units=len(d),
                original_same_cases_point=before,marker_complement_point=after,original_status=status0,complement_status=status1,no_P_or_CI=True))
    pd.DataFrame(sensitivity).to_csv(out/'private_all18_marker_complement_coordination_points.tsv',sep='\t',index=False)
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(out.iterdir()) if p.is_file()};assert sum(f['bytes'] for f in files.values())<16*1024**2
    (out/'private_terminal.json').write_text(json.dumps(dict(status='COMPLETE_ALL_FIXED_LINEAGE_AND_STATE_DIAGNOSTICS_INDEPENDENT_ACCEPTANCE_PENDING',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),files=files,native_cells=19018,fixed_genes=43,source_capture_descriptions=36,state_readout_points=54,coordination_sensitivities=18,
        no_identity_changes_or_new_P=True),indent=2)+'\n')
    print('ALL_FIXED_LINEAGE_AND_STATE_DIAGNOSTIC_POSITIONS_COMPLETE_ACCEPTANCE_PENDING')


if __name__=='__main__':main()
