"""Post-selection, same-cell source/lineage diagnostics for every supported lead.

Frozen labels are retained. Marker co-RNA is not secretion, TCR or doublet proof.
All biological outputs stay in the private analysis directory.
"""
import json
import os
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from scipy.sparse import csr_matrix

from NSCLC_external_lineage_inputs_v2 import load, cell_arms
from external_functional_coverage_inputs_v1 import read_inputs
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource

PANELS = dict(T_core=['CD3D','CD3E','TRAC'], CD8=['CD8A','CD8B'],
              NK_cytotoxic=['NKG7','GNLY','KLRD1'], gamma_delta=['TRDC','TRGC1','TRGC2'],
              myeloid=['LYZ','LST1'], epithelial=['EPCAM','KRT19'])
OUT = 'NSCLC_supported_external_source_identity_20261018_v1'


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget=resource(root,32*1024**2)
    assert int(next(x for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')).split()[1])*1024>3*1024**3
    src=root/'NSCLC_candidate_external_association_results_20261018_v1'
    proof=src/'private_independent_external_candidate_result_audit.json'
    assert sha(proof)=='dc9e287dea77c2044f03f4b7087c3cc1462436f2d913088315e33736dd0a818b'
    terminal=json.loads((src/'private_terminal.json').read_text())
    assert json.loads(proof.read_text())['producer_terminal_sha256']==sha(src/'private_terminal.json')
    name='private_all_estimable_external_effects_and_sensitivities.tsv.gz'
    assert sha(src/name)==terminal['files'][name]['sha256']
    effects=pd.read_csv(src/name,sep='\t')
    scope_path=root/'NSCLC_candidate_external_association_scope_20261018_v1.json'
    assert sha(scope_path)=='d4141c94fe11d2fbccd764fbc2bff503d8a068e7e32132cb6f5c790bb158d752'
    candidates=pd.DataFrame(json.loads(scope_path.read_text())['candidates']).set_index('record')
    cases_path=root/'NSCLC_candidate_external_association_inputs_20261018_v1/private_all_case_measurement_units.tsv.gz'
    case_terminal=json.loads((cases_path.parent/'private_terminal.json').read_text())
    assert sha(cases_path)==case_terminal['files'][cases_path.name]['sha256']
    cases=pd.read_csv(cases_path,sep='\t')
    keep=effects.family.eq('all_units') & effects.partial_rank_rho.gt(0) & effects.BH_q.lt(.05) & effects.MC_upper_BH_q.lt(.05) & effects.CI_low.gt(0) & effects.LOO_all_valid & effects.LOO_same_sign
    selected=effects.loc[keep].copy();assert len(selected)>0
    selected['ligand']=selected.record.map(candidates.ligand)
    data=load();ref=read_inputs(root)
    out=root/OUT;assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),budget=budget,
              accepted_effect_audit_sha256=sha(proof),case_sha256=sha(cases_path),metadata_hashes=data['hashes'],gene_axes=ref['input_sha256'],
              lead_rule='All main-family positive BH<.05, MC-upper-BH<.05, CI lower>0, all LOO valid and same sign. This is explicitly POST-SELECTION explanation, not a new discovery family.',
              panels=PANELS,source_count_thresholds=[1,2,3],
              scope='Every canonical annotated cell in exactly the paired libraries of each supported lead, including all frozen nodes/major lineages. No source relabel or outcome testing.',
              marker_logic='T_core>=1/>=2 genes; CD8>=1; NK/cytotoxic>=1/>=2; gamma_delta/myeloid/epithelial>=1. Report overlaps; NK/cytotoxic is not definitive NK identity. A diagnostic requires every member of its panel uniquely measurable; missing is NaN.',
              limitations='Co-RNA does not establish endogenous synthesis, secretion, protein receptor assembly or eliminate ambient RNA/doublets. Frozen subtype labels are not formal cell states. No new P/CI or main reclassification.',
              selected_records=selected[['model','record','arm','dataset','ligand','active_case_rows_json']].to_dict('records'))
    (out/'plan_before_source_positive_cell_expression.json').write_text(json.dumps(plan,indent=2)+'\n')
    tables=[];genetables=[];replays=[]
    for arm in cell_arms(data):
        todo=selected[selected.arm.eq(arm['arm'])]
        if not len(todo):continue
        original=arm['arm']=='frozen_original_annotation'
        z=zarr.open_group(str(data['old_zarr'] if original else data['atlas']),mode='r')
        raw=z['X'] if original else z['layers']['counts'];nr,nf=map(int,raw.attrs['shape'])
        rr=np.arange(arm['cell_count'])
        if original:
            with np.load(data['old']/'private_canonical_row_selection.npz',allow_pickle=False) as f:rr=f['global_rows']
        nodes=arm['node'];nn=len(nodes.categories)
        for lead in todo.itertuples():
            c=cases.iloc[json.loads(lead.active_case_rows_json)]
            libs=np.unique(np.concatenate([np.array(json.loads(v),int) for v in c.library_indices_json]))
            assert data['lib'].loc[libs,'dataset'].eq(lead.dataset).all()
            take=np.isin(arm['library_index'],libs)
            packed=arm['library_index'][take]*nn+nodes.codes[take]
            keys,inv=np.unique(packed,return_inverse=True);ng=len(keys)
            group=np.full(nr,-1,np.int32);group[rr[take]]=inv
            glib=keys//nn;gnode=keys%nn
            genes=[lead.ligand]+sorted({v for vv in PANELS.values() for v in vv})
            assert len(set(genes))==len(genes)
            axis=ref['source_axes'][(arm['arm'],lead.dataset)]['lookup']
            positions=np.array([next(iter(axis[g])) if g in axis and len(axis[g])==1 else -1 for g in genes],int)
            multiplicity=Counter(positions[positions>=0]);positions=np.array([p if p>=0 and multiplicity[p]==1 else -1 for p in positions])
            if original:
                p=root.parents[2]/'Project_v3/annotation_rebuild/20261010_uniform_expression_annotation_v1/by_gse'/lead.dataset/'genes.tsv'
                assert sha(p)==ref['legacy_gene_axis_sha256'][lead.dataset]
                sourcecols=pd.read_csv(p,sep='\t').core_column.to_numpy(int)
                physical=np.array([sourcecols[p] if p>=0 else -1 for p in positions])
            else:physical=positions
            measured=physical>=0;assert measured[0]
            gi={g:i for i,g in enumerate(genes)};panelcols={k:np.array([gi[g] for g in v],int) for k,v in PANELS.items()}
            metrics=['n_cells','source_positive','source_sum_raw','source_sum_logCP10K','full_RNA_sum']
            for threshold in [1,2,3]:
                metrics += [f'ligand_ge{threshold}']+[f'ligand_ge{threshold}_{k}' for k in ['Tcore_ge1','Tcore_ge2','CD8_ge1','NK_ge1','NK_ge2','Tcore_ge2_AND_NK_ge2','gamma_delta_ge1','myeloid_ge1','epithelial_ge1']]
            values=np.zeros((ng,len(metrics)));marginal=np.zeros((ng,len(genes)),np.int64)
            for start in range(0,nr,4096):
                stop=min(start+4096,nr);valid=group[start:stop]>=0
                if not valid.any():continue
                ptr=np.asarray(raw['indptr'][start:stop+1],np.int64);lo,hi=int(ptr[0]),int(ptr[-1])
                x=csr_matrix((np.asarray(raw['data'][lo:hi]),np.asarray(raw['indices'][lo:hi]),ptr-lo),shape=(stop-start,nf))[valid].copy()
                x.sum_duplicates();x.eliminate_zeros();x.sort_indices()
                assert np.isfinite(x.data).all() and (x.data>0).all() and np.equal(x.data,np.floor(x.data)).all()
                p=np.zeros((x.shape[0],len(genes)),np.int64);p[:,measured]=x[:,physical[measured]].toarray()
                total=np.asarray(x.sum(axis=1)).ravel();assert (p<=total[:,None]).all()
                dest=group[start:stop][valid];counts={k:(p[:,v]>0).sum(axis=1) for k,v in panelcols.items()}
                hits=dict(Tcore_ge1=counts['T_core']>=1,Tcore_ge2=counts['T_core']>=2,CD8_ge1=counts['CD8']>=1,NK_ge1=counts['NK_cytotoxic']>=1,NK_ge2=counts['NK_cytotoxic']>=2,gamma_delta_ge1=counts['gamma_delta']>=1,myeloid_ge1=counts['myeloid']>=1,epithelial_ge1=counts['epithelial']>=1)
                hits['Tcore_ge2_AND_NK_ge2']=hits['Tcore_ge2']&hits['NK_ge2']
                m=[np.ones(len(dest)),p[:,0]>0,p[:,0],np.log1p(p[:,0]*np.divide(10000.,total,out=np.zeros(len(total)),where=total>0)),total]
                for t in [1,2,3]:
                    positive=p[:,0]>=t;m.append(positive)
                    for k in ['Tcore_ge1','Tcore_ge2','CD8_ge1','NK_ge1','NK_ge2','Tcore_ge2_AND_NK_ge2','gamma_delta_ge1','myeloid_ge1','epithelial_ge1']:m.append(positive&hits[k])
                grouping=csr_matrix((np.ones(len(dest)),(dest,np.arange(len(dest)))),shape=(ng,len(dest)))
                values+=grouping@np.column_stack(m);marginal+=(grouping@(p>0)).astype(np.int64)
            table=pd.DataFrame(values,columns=metrics)
            for name,panel in [('Tcore_ge1','T_core'),('Tcore_ge2','T_core'),('CD8_ge1','CD8'),('NK_ge1','NK_cytotoxic'),('NK_ge2','NK_cytotoxic'),('gamma_delta_ge1','gamma_delta'),('myeloid_ge1','myeloid'),('epithelial_ge1','epithelial')]:
                if not measured[panelcols[panel]].all():
                    for t in [1,2,3]:table[f'ligand_ge{t}_{name}']=np.nan
            if not (measured[panelcols['T_core']].all() and measured[panelcols['NK_cytotoxic']].all()):
                for t in [1,2,3]:table[f'ligand_ge{t}_Tcore_ge2_AND_NK_ge2']=np.nan
            cat=data['catalogs'][arm['arm']].set_index('node_id');nodeids=np.asarray(nodes.categories).astype(str)[gnode]
            for k,v in dict(model=lead.model,record=lead.record,ligand=lead.ligand,arm=arm['arm'],dataset=lead.dataset,library_index=glib,node_id=nodeids,major=cat.loc[nodeids,'major'].to_numpy(str),frozen_role=cat.loc[nodeids,'frozen_role'].to_numpy(int)).items():table[k]=v
            expected=data['counts'][data['counts'].arm.eq(arm['arm']) & data['counts'].dataset.eq(lead.dataset)]
            lookup={(r.library_id,r.node_id):r.n_cells for r in expected.itertuples()}
            assert np.array_equal(table.n_cells,[lookup[(data['lib'].iloc[li].library_id,ni)] for li,ni in zip(glib,nodeids)])
            old=root/'NSCLC_candidate_external_exact_RNA_20261018_v1/private_all_frozen_pool_exact_panel_measurements.npz'
            with np.load(old,allow_pickle=False) as z0:
                col=list(z0['genes']).index(lead.ligand)
                for li in libs:
                    for v in range(4):
                        g=table[table.library_index.eq(li)&table.frozen_role.eq(v)];u=li*4+v
                        for a,b in [(g.n_cells.sum(),z0['n_cells'][u]),(g.source_positive.sum(),z0['positive_counts'][u,col]),(g.source_sum_logCP10K.sum(),z0['sum_logCP10K'][u,col])]:assert np.isclose(a,b,atol=1e-10,rtol=1e-10,equal_nan=True)
            for j,gene in enumerate(genes):genetables.append(dict(model=lead.model,dataset=lead.dataset,gene=gene,measured=bool(measured[j]),all_selected_cells_positive=int(marginal[:,j].sum()) if measured[j] else None))
            tables.append(table);replays.append(dict(model=lead.model,selected_libraries=len(libs),selected_cells=int(table.n_cells.sum()),all_nodes=len(table),counts_and_original_ligand_marginals_replayed=True))
        if len(todo)==len(selected):break
    dataout=pd.concat(tables,ignore_index=True)
    dataout.to_csv(out/'private_all_library_node_source_lineage_counts.tsv.gz',sep='\t',index=False,compression='gzip')
    pd.DataFrame(genetables).to_csv(out/'private_all_selected_gene_measurability_and_marginals.tsv',sep='\t',index=False)
    groupcols=['model','record','ligand','arm','dataset','major','frozen_role']
    numeric=[c for c in dataout.columns if c not in groupcols+['library_index','node_id']]
    summary=dataout.groupby(groupcols,dropna=False)[numeric].sum(min_count=1).reset_index()
    summary.to_csv(out/'private_major_level_source_identity_summary.tsv',sep='\t',index=False)
    files=[p for p in out.iterdir() if p.is_file()];assert sum(p.stat().st_size for p in files)<32*1024**2
    term=dict(status='COMPLETE_POST_SELECTION_SOURCE_IDENTITY_DIAGNOSTIC_MARGINALS_CHECKED_NOT_INDEPENDENT_JOINT_REPLAY',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),selected_models=len(selected),replay=replays,files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in files},new_P_CI_or_annotation_changes=False)
    (out/'private_terminal.json').write_text(json.dumps(term,indent=2)+'\n')
    print(json.dumps(dict(status=term['status'],completed_UTC=term['completed_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))

if __name__=='__main__':main()
