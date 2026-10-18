"""Independent all-selected-cell joint-marker/count replay, original-core leads.

Shares accepted annotation I/O; aggregation uses explicit per-metric np.add.at.
No new identity labels or statistical tests. Scientific tables remain private.
"""
import json
import os
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from scipy.sparse import csr_matrix

from NSCLC_external_lineage_inputs_v2 import load,cell_arms
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from external_functional_coverage_inputs_v1 import read_inputs
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from measure_NSCLC_supported_external_source_identity_v1 import OUT,PANELS


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);resource(root,8*1024**2)
    out=root/OUT;terminal=json.loads((out/'private_terminal.json').read_text())
    assert terminal['source_sha256']==sha(Path(__file__).with_name('measure_NSCLC_supported_external_source_identity_v1.py'))
    for name,item in terminal['files'].items():assert sha(out/name)==item['sha256']
    plan=json.loads((out/'plan_before_source_positive_cell_expression.json').read_text());assert plan['panels']==PANELS
    assert all(r['arm']=='frozen_original_annotation' for r in plan['selected_records'])
    table=pd.read_csv(out/'private_all_library_node_source_lineage_counts.tsv.gz',sep='\t')
    marg=pd.read_csv(out/'private_all_selected_gene_measurability_and_marginals.tsv',sep='\t')
    data=load();ref=read_inputs(root);assert data['hashes']==plan['metadata_hashes']
    arm=next(cell_arms(data));assert arm['arm']=='frozen_original_annotation'
    with np.load(data['old']/'private_canonical_row_selection.npz',allow_pickle=False) as z:rows=z['global_rows']
    z=zarr.open_group(str(data['old_zarr']),mode='r');raw=z['X'];nr,nf=map(int,raw.attrs['shape'])
    maximum=0.;selectedcells=0;fullmetrics=0
    for lead in plan['selected_records']:
        g=table[table.model.eq(lead['model'])].reset_index(drop=True)
        lookup={(int(v.library_index),str(v.node_id)):i for i,v in enumerate(g.itertuples())}
        libs=np.array(sorted(g.library_index.unique()),int);take=np.isin(arm['library_index'],libs)
        ids=np.full(nr,-1,np.int32)
        nodes=np.asarray(arm['node'].categories).astype(str)[arm['node'].codes[take]]
        ids[rows[take]]=[lookup[(int(li),ni)] for li,ni in zip(arm['library_index'][take],nodes)]
        genes=[lead['ligand']]+sorted({v for vv in PANELS.values() for v in vv})
        source=ref['source_axes'][(lead['arm'],lead['dataset'])]['lookup']
        mapping=root.parents[2]/'Project_v3/annotation_rebuild/20261010_uniform_expression_annotation_v1/by_gse'/lead['dataset']/'genes.tsv'
        assert sha(mapping)==ref['legacy_gene_axis_sha256'][lead['dataset']]
        columns=pd.read_csv(mapping,sep='\t').core_column.to_numpy(int)
        positions=np.array([columns[next(iter(source[v]))] if v in source and len(source[v])==1 else -1 for v in genes],int)
        validgene=positions>=0;assert len(set(positions[validgene]))==validgene.sum()
        numeric=[c for c in g.columns if c in ['n_cells','source_positive','source_sum_raw','source_sum_logCP10K','full_RNA_sum'] or c.startswith('ligand_ge')]
        sums={k:np.zeros(len(g)) for k in numeric};positive=np.zeros(len(genes),np.int64)
        for start in range(0,nr,1024):
            stop=min(start+1024,nr);keep=ids[start:stop]>=0
            if not keep.any():continue
            ptr=np.asarray(raw['indptr'][start:stop+1],np.int64);lo,hi=int(ptr[0]),int(ptr[-1])
            x=csr_matrix((np.asarray(raw['data'][lo:hi]),np.asarray(raw['indices'][lo:hi]),ptr-lo),shape=(stop-start,nf))[keep].copy()
            x.sum_duplicates();x.eliminate_zeros()
            panel=np.zeros((x.shape[0],len(genes)),np.int64)
            for j in np.flatnonzero(validgene):panel[:,j]=np.asarray(x[:,int(positions[j])].toarray()).ravel()
            positive+=(panel>0).sum(axis=0);dest=ids[start:stop][keep];total=np.asarray(x.sum(axis=1)).ravel();selectedcells+=len(dest)
            sourcecount=panel[:,0];q=dict(n_cells=np.ones(len(dest)),source_positive=sourcecount>=1,source_sum_raw=sourcecount,source_sum_logCP10K=np.log1p(np.divide(sourcecount*10000.,total,out=np.zeros(len(total)),where=total>0)),full_RNA_sum=total)
            checks={}
            for k,p in PANELS.items():
                js=[genes.index(v) for v in p]
                checks[k]=((panel[:,js]>0).sum(axis=1),bool(validgene[js].all()))
            combos={'Tcore_ge1':('T_core',1),'Tcore_ge2':('T_core',2),'CD8_ge1':('CD8',1),'NK_ge1':('NK_cytotoxic',1),'NK_ge2':('NK_cytotoxic',2),'gamma_delta_ge1':('gamma_delta',1),'myeloid_ge1':('myeloid',1),'epithelial_ge1':('epithelial',1)}
            for threshold in [1,2,3]:
                hit=sourcecount>=threshold;q[f'ligand_ge{threshold}']=hit
                for name,(p,n) in combos.items():q[f'ligand_ge{threshold}_{name}']=hit&(checks[p][0]>=n) if checks[p][1] else np.full(len(dest),np.nan)
                ok=checks['T_core'][1] and checks['NK_cytotoxic'][1]
                q[f'ligand_ge{threshold}_Tcore_ge2_AND_NK_ge2']=hit&(checks['T_core'][0]>=2)&(checks['NK_cytotoxic'][0]>=2) if ok else np.full(len(dest),np.nan)
            for name in numeric:np.add.at(sums[name],dest,q[name])
        for name in numeric:
            expected=g[name].to_numpy(float);actual=sums[name]
            assert np.array_equal(np.isnan(actual),np.isnan(expected)),name
            assert np.allclose(actual,expected,atol=1e-9,rtol=1e-10,equal_nan=True),name
            finite=np.isfinite(actual)&np.isfinite(expected)
            if finite.any():maximum=max(maximum,float(np.max(np.abs(actual[finite]-expected[finite]))))
            fullmetrics+=len(g)
        mg=marg[marg.model.eq(lead['model'])].set_index('gene')
        for j,name in enumerate(genes):
            assert bool(mg.loc[name,'measured'])==bool(validgene[j])
            if validgene[j]:assert int(mg.loc[name,'all_selected_cells_positive'])==int(positive[j])
    proof=out/'private_independent_full_selected_cell_identity_audit.json';assert not proof.exists()
    result=dict(status='PASS_ALL_SELECTED_CANONICAL_CELLS_ALL_NODE_METRICS_ALL_THRESHOLDS_INDEPENDENT_AGGREGATION',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),producer_terminal_sha256=sha(out/'private_terminal.json'),models=len(plan['selected_records']),selected_cells_replayed=selectedcells,node_metrics_replayed=fullmetrics,max_absolute_difference=maximum,scope='All selected raw cells and all joint-marker thresholds; accepted annotation/category I/O and marker contract shared. Not independent cell-identity truth, secretion or biological replication.')
    proof.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(status=result['status'],completed_UTC=result['completed_UTC'],audit_sha256=sha(proof))))

if __name__=='__main__':main()
