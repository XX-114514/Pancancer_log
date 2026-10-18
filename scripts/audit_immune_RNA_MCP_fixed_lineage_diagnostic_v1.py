"""Independent full raw43 assay replay and all source diagnostic descriptions."""
import json
import os
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from anndata.io import read_elem, sparse_dataset

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from measure_immune_RNA_MCP_fixed_lineage_diagnostic_v1 import PANELS, THRESHOLDS
from audit_immune_RNA_MCP_stability_and_coordination_v1 import independent_fit, rho, equal


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,8*1024**2)
    source=root/'immune_RNA_MCP_fixed_lineage_diagnostic_20261018_v1'
    assert sha(source/'private_terminal.json')=='0ba53c8accfb21462271154b90f98cd4940f69efc85905d2f10d5a29aebe9c35'
    for name,item in json.loads((source/'private_terminal.json').read_text())['files'].items():assert sha(source/name)==item['sha256']
    out=root/'immune_RNA_MCP_fixed_lineage_diagnostic_audit_20261018_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_raw43_and_all_source_replay.json').write_text(json.dumps(dict(
        created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),guard=guard,
        producer_terminal_sha256=sha(source/'private_terminal.json'),scope='All19018 original cells/all43 native-marker slots from alternate unitwise CSR, all marker/clear flags,396 case scores,36 source capture counts,54 state points,18 complement points. Cache arrays once; no identity or main-P modification.'),indent=2)+'\n')
    with np.load(source/'private_all19018_fixed43_raw_RNA.npz',allow_pickle=False) as blob:raw={k:blob[k] for k in blob.files}
    genes=raw['genes'].tolist();assert genes==sorted(set(g for v in PANELS.values() for g in v)) and len(genes)==43
    inputs=root/'immune_RNA_MCP_native_expression_inputs_20261018_v1'
    md=pd.read_csv(inputs/'private_selected_native_cell_depth.tsv',sep='\t').rename(columns={'fit_sample':'sample'})
    assert np.array_equal(md.global_row,raw['global_rows']) and np.array_equal(md.fullRNA,raw['fullRNA'])
    stored=pd.read_csv(source/'private_fixed_marker_source_columns.tsv',sep='\t').set_index('gene')
    z=zarr.open_group(str(root.parents[1]/'runs/20261013_immunotherapy_atlas_v1/outputs/immune_validation_authority_v1.zarr'),mode='r')
    order=list(np.asarray(read_elem(z['uns/source_order'])).astype(str));measured=np.asarray(read_elem(z['varm/measured_in_source']),bool)
    symbol=np.asarray(read_elem(z['var/gene_symbol'])).astype(str);canonical=np.asarray(read_elem(z['var/canonical_gene_identity'])).astype(str)
    identity=np.asarray(read_elem(z['var/feature_identity_status'])).astype(str);aliases=np.asarray(read_elem(z['var/source_gene_aliases'])).astype(str)
    col={};requested=set(genes)
    for dataset in sorted(md.dataset.unique()):
        labels=defaultdict(set)
        for j in np.flatnonzero(measured[:,order.index(dataset)]):
            if canonical[j].startswith('AMBIGUOUS::') or identity[j]=='ambiguous_source_feature':continue
            names={symbol[j],canonical[j]}
            if identity[j]=='stable_Ensembl_ID':names|=set(aliases[j].split('|'))
            for name in names&requested:labels[name].add(int(j))
        col[dataset]=np.array([next(iter(labels[g])) if len(labels[g])==1 else -1 for g in genes],int)
        assert np.array_equal(col[dataset],stored.loc[genes,dataset])
        for names in PANELS.values():
            values=[col[dataset][genes.index(g)] for g in names];valid=[v for v in values if v>=0]
            assert len(valid)==len(set(valid))
    counts=sparse_dataset(z['layers/counts']);replay=np.full_like(raw['raw'],-1);read=0
    for (_,dataset,_),g in md.groupby(['sample','dataset','view'],sort=True):
        x=counts[g.global_row.to_numpy(int),:].tocsr();x.sum_duplicates();x.eliminate_zeros();x.sort_indices()
        assert np.array_equal(np.asarray(x.sum(axis=1)).ravel(),g.fullRNA) and np.array_equal(np.diff(x.indptr),g.full_detected_genes)
        valid=col[dataset]>=0;replay[np.ix_(g.index.to_numpy(),np.flatnonzero(valid))]=x[:,col[dataset][valid]].toarray();read+=len(g)
    assert read==19018 and np.array_equal(replay,raw['raw'])
    saved_cells=pd.read_csv(source/'private_all_original_cells_fixed_markers.tsv.gz',sep='\t')
    assert saved_cells.cell.tolist()==md.cell.tolist()
    logged=np.log1p(replay.astype(float)*10000/md.fullRNA.to_numpy()[:,None]);logged[replay<0]=np.nan
    for panel,names in PANELS.items():
        where=[genes.index(g) for g in names];valid=(replay[:,where]>=0).all(axis=1)
        values=np.where(valid,logged[:,where].mean(axis=1),np.nan);assert np.allclose(values,saved_cells[panel+'_meanlog'],atol=1e-9,rtol=0,equal_nan=True)
        md[panel+'_meanlog']=values
        if panel in THRESHOLDS:
            value=np.where(valid,(replay[:,where]>0).sum(axis=1)>=THRESHOLDS[panel],np.nan)
            assert np.allclose(value,saved_cells[panel+'_capture'],equal_nan=True);md[panel+'_capture']=value
    md['marker_complement']=False
    for view,g in md.groupby('view'):
        foreign=['T_core','B_core','epithelial'] if view=='Myeloid' else ['Myeloid_core','B_core','epithelial']
        values=g[[k+'_capture' for k in foreign]].to_numpy();md.loc[g.index,'marker_complement']=np.isfinite(values).all(axis=1)&(values==0).all(axis=1)
    assert np.array_equal(md.marker_complement,saved_cells.marker_complement)
    q=[]
    for view in md.view.unique():
        for part in ['discovery','TNBC_holdout']:q.append(pd.read_csv(inputs/'private_inputs'/(view+'_'+part+'_metadata.tsv'),sep='\t')[['cell','cellQ','geneQ']])
    md=md.merge(pd.concat(q,ignore_index=True),on='cell',validate='one_to_one')
    old=root/'existing_immune_fine_native_measurement_scope_20261016_v1';mapping=np.load(old/'private_original_cell_to_fine_node_unit.npy',mmap_mode='r')
    units=pd.read_csv(old/'private_all_frozen_T_myeloid_node_units.tsv',sep='\t').set_index('unit_index')
    md['node_id']=units.node_id.loc[mapping[md.global_row.to_numpy(int)]].to_numpy()
    scores=pd.read_csv(root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection/private_all_cell_score_reconstruction.tsv.gz',sep='\t')
    scores=scores[scores.fit_id.eq('K2_seed20261018_source_adjusted')][['cell','view','program','score']]
    cell=scores.merge(md,on=['cell','view'],validate='many_to_one');assert len(cell)==38036
    coefficient=pd.read_csv(root/'immune_RNA_MCP_original_node_decomposition_20261018_v1/private_frozen_discovery_node_coefficients.tsv',sep='\t')[['view','program','node_id','coefficient']]
    cell=cell.merge(coefficient,on=['view','program','node_id'],how='left',validate='many_to_one')
    meta=pd.read_csv(root/'immune_RNA_MCP_source_coordination_20261018_v1/private_all_case_view_controls.tsv',sep='\t').set_index(['sample','view'])
    saved=pd.read_csv(source/'private_all396_original_and_marker_complement_scores.tsv',sep='\t').set_index(['sample','view','program'])
    rows=[]
    for (sample,dataset,view,program),g in cell.groupby(['sample','dataset','view','program'],sort=True):
        clear=g[g.marker_complement];m=meta.loc[(sample,view)]
        r=dict(sample=sample,dataset=dataset,view=view,program=program,selected_cells=len(g),complement_cells=len(clear),score_original=float(g.score.mean()),score_complement=float(clear.score.mean()) if len(clear) else np.nan,
            composition_original=float(g.coefficient.mean()) if g.coefficient.notna().all() else np.nan,
            composition_complement=float(clear.coefficient.mean()) if len(clear) and clear.coefficient.notna().all() else np.nan,
            fraction=float(m.fraction),logN=float(m.logN),cellQ_original=float(g.cellQ.mean()),geneQ_original=float(g.geneQ.mean()),
            cellQ_complement=float(clear.cellQ.mean()) if len(clear) else np.nan,geneQ_complement=float(clear.geneQ.mean()) if len(clear) else np.nan)
        for panel in ['IFN','proliferation','immediate_early']:r[panel]=float(g[panel+'_meanlog'].mean())
        s=saved.loc[(sample,view,program)]
        for name,value in r.items():
            if name not in ['sample','dataset','view','program']:equal(value,s[name])
        rows.append(r)
    case=pd.DataFrame(rows);assert len(case)==396
    capture=pd.read_csv(source/'private_all36_source_marker_capture_descriptions.tsv',sep='\t')
    for r in capture.itertuples(index=False):
        g=md[md.dataset.eq(r.dataset)&md.view.eq(r.view)];f=g[r.panel+'_capture']
        assert r.selected_cells==len(g) and r.case_units==g['sample'].nunique() and r.measurable_cells==f.notna().sum() and r.capture_cells==f.eq(1).sum()
    state=pd.read_csv(source/'private_all54_fixed_state_readout_points.tsv',sep='\t')
    for r in state.itertuples(index=False):
        g=case[case.dataset.eq(r.dataset)&case.view.eq(r.view)&case.program.eq(r.program)]
        valid=np.isfinite(g[[r.panel,'score_original']].to_numpy()).all(axis=1);g=g[valid]
        x=g[r.panel].to_numpy();y=g.score_original.to_numpy();C=g[['fraction','logN','cellQ_original','geneQ_original']].to_numpy(float)
        fit,status=independent_fit(x,y,C);assert status==r.status and r.complete_units==len(g)
        equal(r.raw_rho,rho(x,y));equal(r.adjusted_point,fit['effect'] if fit else np.nan)
    sensitivity=pd.read_csv(source/'private_all18_marker_complement_coordination_points.tsv',sep='\t')
    for r in sensitivity.itertuples(index=False):
        a=case[case.view.eq(r.view_A)&case.program.eq(r.program)&case.dataset.eq(r.dataset)]
        b=case[case.view.eq(r.view_B)&case.program.eq(r.program)&case.dataset.eq(r.dataset)]
        g=a.merge(b,on=['sample','dataset','program'],suffixes=('_A','_B'),validate='one_to_one')
        oldC=[c+'_original_'+s for s in ['A','B'] for c in ['cellQ','geneQ','composition']]+[c+'_'+s for s in ['A','B'] for c in ['fraction','logN']]
        newC=[c.replace('_original_','_complement_') for c in oldC]
        fields=['score_original_A','score_original_B','score_complement_A','score_complement_B']+oldC+newC
        valid=(g.complement_cells_A.ge(20)&g.complement_cells_B.ge(20)).to_numpy()&np.isfinite(g[fields].to_numpy(float)).all(axis=1);d=g[valid]
        f,s0=independent_fit(d.score_original_A.to_numpy(),d.score_original_B.to_numpy(),d[oldC].to_numpy(float))
        h,s1=independent_fit(d.score_complement_A.to_numpy(),d.score_complement_B.to_numpy(),d[newC].to_numpy(float))
        assert r.original_status==s0 and r.complement_status==s1 and r.same_comparable_units==len(d)
        equal(r.original_same_cases_point,f['effect'] if f else np.nan);equal(r.marker_complement_point,h['effect'] if h else np.nan)
    assert len(capture)==36 and len(state)==54 and len(sensitivity)==18
    (out/'private_terminal.json').write_text(json.dumps(dict(status='ACCEPTED_ALL_RAW43_FIXED_LINEAGE_FLAGS_STATE_AND_COMPLEMENT_DESCRIPTIONS',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),producer_terminal_sha256=sha(source/'private_terminal.json'),
        native_cells_replayed=19018,fixed_gene_slots=19018*43,measured_raw_values=int((replay>=0).sum()),case_scores_replayed=396,
        source_capture_descriptions=36,state_points=54,complement_points=18,no_identity_or_main_P_changes=True),indent=2)+'\n')
    print('ALL_FIXED_LINEAGE_AND_STATE_DIAGNOSTIC_POSITIONS_INDEPENDENTLY_ACCEPTED')


if __name__=='__main__':main()
