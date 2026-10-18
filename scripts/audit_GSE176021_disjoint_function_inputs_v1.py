"""Rebuild every functional input directly from accepted individual score arrays."""
import json
import os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from immune_RNA_MCP_new_source_io_v1 import wait


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,2*1024**2)
    out=root/'GSE176021_disjoint_function_input_audit_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    src=root/'GSE176021_disjoint_function_inputs_20261019_v1'
    term=wait(root,out,'prepare_GSE176021_disjoint_function_inputs_v1',src/'private_terminal.json',{'COMPLETE_ALL400_FUNCTIONAL_INPUTS_AND3200_POINT_SENSITIVITIES_AUDIT_PENDING'})
    for name,v in term['files'].items():assert sha(src/name)==v['sha256']
    assert sha(root/'GSE176021_disjoint_function_native_RNA_audit_20261019_v1/private_terminal.json')=='b07e68106a615cb1f32f3e33f7acac18b5624c4f3cb82aa0180857d768bad21f'
    a=dict(np.load(src/'private_all400_case_inputs.npz'));reg=pd.read_csv(src/'private_all400_registered_functional_questions.tsv',sep='\t')
    inp=root/'GSE176021_paired_tissue_inputs_20261019_v1';cases=np.load(inp/'private_AP01_inputs.npz')['cases'];assert np.array_equal(a['cases'],cases) and np.array_equal(a['model_ids'],reg.model_id)
    scope=root/'GSE176021_author_lineage_RNA_scope_20261019_v2';pairs=pd.read_csv(scope/'private_all_scope_material_pair_positions.tsv',sep='\t');pairs=pairs[pairs.scope.eq('all_author_lineage')&pairs.case_block.isin(cases)&pairs.tissue.isin(['tumor','adjacent normal lung'])&pairs.both_views_eligible];assert len(pairs)==60
    pools=pd.read_csv(scope/'private_all_eligible_pools.tsv',sep='\t');pools=pools[pools.scope.eq('all_author_lineage')&pools.sampleID.isin(pairs.sampleID)].merge(pairs[['sampleID','case_block','tissue']],on='sampleID',validate='many_to_one');assert len(pools)==120
    ss=dict(np.load(root/'GSE176021_disjoint_function_native_RNA_20261019_v1/private_all_selected_cell_functional_scores.npz'));old=dict(np.load(root/'GSE176021_author_lineage_native_RNA_20261019_v1/private_selected_cell_quality_scores.npz'))
    scales=pd.read_csv(inp/'private_original_discovery_case_scales.tsv',sep='\t').set_index(['view','program']).discovery_case_sd
    values=pd.read_csv(inp/'private_all_scope_case_tissue_values.tsv',sep='\t');values=values[values.scope.eq('all_author_lineage')].set_index(['case_block','tissue','view','program','mode']).value
    tissue_functions={};tissue_scores={};tissue_quality={}
    for cc in cases:
        for tissue in ['tumor','adjacent normal lung']:
            for view in ['R_FOXP3_capture','E_effector_capture']:
                group=pools[(pools.case_block==cc)&(pools.tissue==tissue)&(pools['view']==view)];assert len(group)>0
                funcs=[];scores=[];quality=[]
                for p in group.itertuples(index=False):
                    ix=np.flatnonzero(ss['pool_index']==p.pool_index);oi=np.flatnonzero(old['pool_index']==p.pool_index);assert len(ix)==len(oi)==p.selected_cells and np.array_equal(ss['global_row'][ix],np.sort(old['global_row'][oi]))
                    funcs.append(np.c_[ss['ucell'][ix].astype(float).sum(axis=0)/len(ix),ss['mean_logCP10K'][ix].sum(axis=0)/len(ix)])
                    scores.append(old['scores'][oi].sum(axis=0)/len(oi));quality.append(old['quality'][oi,:2].sum(axis=0)/len(oi))
                tissue_functions[(cc,tissue,view)]=sum(funcs)/len(funcs);tissue_scores[(cc,tissue,view)]=sum(scores)/len(scores);tissue_quality[(cc,tissue,view)]=sum(quality)/len(quality)
    refs=ss['references'].tolist();seen=set();modes=['ucell1000','ucell2000','mean_logCP10K','fixed_stress_complement','own_population_background','unresidualized_weighted_RNA','RNA_quality_only_background','response_quality_partial_two'];assert a['modes'].tolist()==modes
    for i,r in enumerate(reg.itertuples(index=False)):
        key=(r.reference,r.program,r.predictor_view,r.response_view);assert key not in seen;seen.add(key)
        assert r.model_id=='CF'+str(i+1).zfill(3) and r.primary_case_units==10 and r.primary_materials==60 and r.family_denominator==400 and r.metadata_preeligible
        assert r.role==('WITHIN_LINEAGE' if r.predictor_view==r.response_view else 'CROSS_LINEAGE')
        rid=refs.index(r.reference);pc=2*(int(r.program[-1])-1);sd=scales.loc[(r.predictor_view,r.program)]
        x=[];y=[];q=[];sensx=[[] for _ in modes];sensy=[[] for _ in modes]
        for cc in cases:
            T,N=(cc,'tumor',r.predictor_view),(cc,'adjacent normal lung',r.predictor_view)
            RT,RN=(cc,'tumor',r.response_view),(cc,'adjacent normal lung',r.response_view)
            xx=(tissue_scores[T][pc]-tissue_scores[N][pc])/sd;yy=tissue_functions[RT][rid,0]-tissue_functions[RN][rid,0]
            x.append(xx);y.append(yy);q.append(tissue_quality[RT]-tissue_quality[RN])
            for j,mode in enumerate(modes):
                xp=xx;yp=yy
                if j<3:yp=tissue_functions[RT][rid,j+1]-tissue_functions[RN][rid,j+1]
                elif j==3:xp=(tissue_scores[T][pc+1]-tissue_scores[N][pc+1])/sd
                elif j<7:xp=(values.loc[(cc,'tumor',r.predictor_view,r.program,mode)]-values.loc[(cc,'adjacent normal lung',r.predictor_view,r.program,mode)])/sd
                sensx[j].append(xp);sensy[j].append(yp)
        for actual,expected in [(a['x'][i],x),(a['y'][i],y),(a['sensitivity_x'][i],sensx),(a['sensitivity_y'][i],sensy),(a['response_quality'][i],q)]:np.testing.assert_allclose(actual,expected,atol=1e-10,rtol=0)
    expected={(ref,pr,pv,rv) for ref in refs for pr in ['MCP1','MCP2'] for pv in ['R_FOXP3_capture','E_effector_capture'] for rv in ['R_FOXP3_capture','E_effector_capture']};assert seen==expected and len(reg)==400
    receipt={'status':'ACCEPTED_ALL400_FUNCTIONAL_INPUTS3200_SENSITIVITIES_COMPLETE_CASE_AND_ROLE_ALIGNMENT','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'producer_terminal_sha256':sha(src/'private_terminal.json'),'resource_guard':guard,'all_primary_questions':400,'all_sensitivity_positions':3200,'whole_case_units':10,'real_functional_associations_executed':False}
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'status':receipt['status'],'terminal_sha256':sha(out/'private_terminal.json')}))


if __name__=='__main__':main()
