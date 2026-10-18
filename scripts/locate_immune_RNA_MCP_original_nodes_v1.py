"""Locate all programme scores in original nodes and full captured-view mixtures."""
import json
import os
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget=resource(root,16*1024**2)
    scope=root/'immune_RNA_MCP_discovery_scope_20261018_v1'
    assert sha(scope/'private_terminal.json')=='0f4833a07bd5cdf5e8f6806d06805777dddf78adad60a6151eca991dfdda6196'
    for name,item in json.loads((scope/'private_terminal.json').read_text())['files'].items():assert sha(scope/name)==item['sha256']
    chosen=pd.read_csv(scope/'private_metadata_only_case_source_partition.tsv',sep='\t',keep_default_na=False)
    chosen=chosen[~chosen.partition.eq('small_source_projection_only')].copy();assert len(chosen)==66
    chosen=chosen[['sampleID','fit_sample','dataset']].rename(columns={'fit_sample':'sample'})
    native=root/'existing_immune_fine_native_measurement_scope_20261016_v1'
    prior=json.loads((native/'private_terminal.json').read_text())
    for name,digest in prior['files'].items():assert sha(native/name)==digest
    units=pd.read_csv(native/'private_all_frozen_T_myeloid_node_units.tsv',sep='\t',keep_default_na=False)
    units['node_id']=units.node_id.astype(str);assert units.unit_index.is_unique
    mapping=np.load(native/'private_original_cell_to_fine_node_unit.npy',mmap_mode='r')
    fields=['node_id','major','subcluster','subtype_candidate','subtype_evidence_status','identity_promoted']
    catalogue=units[fields].drop_duplicates();assert catalogue.node_id.is_unique
    assert catalogue[catalogue.major.eq('T')].shape[0]==16 and catalogue[catalogue.major.eq('Myeloid')].shape[0]==13
    source=root/'immune_CM_fixed_panel_same_cell_RNA_20261018_v1'
    audit=root/'immune_CM_fixed_panel_same_cell_RNA_audit_20261018_v1/private_terminal.json'
    raw_proof=json.loads(audit.read_text());assert sha(source/'private_terminal.json')==raw_proof['producer_terminal_sha256']
    for name,item in json.loads((source/'private_terminal.json').read_text())['files'].items():assert sha(source/name)==item['sha256']
    T=pd.read_csv(source/'private_original_T_unit_scope.tsv',sep='\t',keep_default_na=False);T['node_id']=T.node_id.astype(str)
    z=np.load(source/'private_fixed34_gene_T_counts_and_depth.npz',allow_pickle=False)
    raw=z['raw_gene_counts'];genes=z['genes'].tolist();tunit=z['original_T_unit_index']
    fixed=json.loads((root.parents[1]/'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json').read_text())
    core=(raw[:,[genes.index(g) for g in fixed['panels']['pan_T']]]>0).sum(1)>=2
    reg=((raw[:,[genes.index(g) for g in fixed['panels']['regulatory_associated']]]>0).sum(1)>=2)&(raw[:,genes.index('FOXP3')]>0)
    eff=(raw[:,[genes.index(g) for g in ['GZMB','PRF1','CTSW','NKG7']]]>0).sum(1)>=2
    full=[]
    for view,mask in [('R_FOXP3_capture',core&reg&~eff),('E_effector_capture',core&eff&~reg)]:
        counts=np.bincount(tunit[mask],minlength=len(T));f=T[['dataset','library_id','node_id']].copy();f['full_node_cells']=counts
        f['sampleID']=f.dataset.astype(str)+'::'+f.library_id.astype(str)
        f=f.merge(chosen[['sampleID','sample']],on='sampleID',how='inner',validate='many_to_one');f['view']=view
        full.append(f[['sample','dataset','view','node_id','full_node_cells']])
    f=units[units.major.eq('Myeloid')].copy();f['sampleID']=f.dataset.astype(str)+'::'+f.library_id.astype(str)
    f=f.merge(chosen[['sampleID','sample']],on='sampleID',how='inner',validate='many_to_one')
    f['view']='Myeloid';f['full_node_cells']=f.n_cells
    full.append(f[['sample','dataset','view','node_id','full_node_cells']])
    full=pd.concat(full,ignore_index=True).groupby(['sample','dataset','view','node_id'],sort=True).full_node_cells.sum().reset_index()
    dense=[]
    for view in ['R_FOXP3_capture','E_effector_capture','Myeloid']:
        nodes=catalogue.loc[catalogue.major.eq('Myeloid' if view=='Myeloid' else 'T'),'node_id'].tolist()
        for case in chosen.itertuples(index=False):
            one=full[full['sample'].eq(case.sample)&full.view.eq(view)].set_index('node_id').full_node_cells.reindex(nodes,fill_value=0)
            total=int(one.sum());assert total>=20
            for node,count in one.items():dense.append(dict(sample=case.sample,dataset=case.dataset,view=view,node_id=node,full_node_cells=int(count),full_view_cells=total,fraction=float(count/total)))
    dense=pd.DataFrame(dense);assert len(dense)==2970
    selected=pd.read_csv(root/'immune_RNA_MCP_native_expression_inputs_20261018_v1/private_frozen_original_cell_selection.tsv',sep='\t')
    selected['unit_index']=mapping[selected.global_row.to_numpy(int)];assert selected.unit_index.ge(0).all()
    selected=selected.merge(units[['unit_index']+fields],on='unit_index',validate='many_to_one')
    assert selected.loc[selected.view.eq('Myeloid'),'major'].eq('Myeloid').all()
    assert selected.loc[~selected.view.eq('Myeloid'),'major'].eq('T').all()
    for view in selected.view.unique():
        actual=selected[selected.view.eq(view)].groupby('fit_sample').full_view_cells.first()
        expected=dense[dense.view.eq(view)].groupby('sample').full_view_cells.first()
        assert np.array_equal(actual.sort_index(),expected.reindex(actual.sort_index().index))
    projection=root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection'
    scores=pd.read_csv(projection/'private_all_cell_score_reconstruction.tsv.gz',sep='\t',dtype={'cell':str,'sample':str})
    scores=scores[scores.fit_id.eq('K2_seed20261018_source_adjusted')]
    scores=scores.merge(selected[['cell','view','node_id']],on=['cell','view'],validate='many_to_one')
    assert len(scores)==38036
    scores['within_case_centered_score']=scores.score-scores.groupby(['sample','view','program']).score.transform('mean')
    means=scores.groupby(['sample','dataset','view','program','node_id'],sort=True).agg(
        selected_cells=('score','size'),mean_score=('score','mean'),within_case_centered_score=('within_case_centered_score','mean')).reset_index()
    case_scores=scores.groupby(['sample','dataset','view','program'],sort=True).score.mean().reset_index()
    output=root/'immune_RNA_MCP_original_node_localization_20261018_v1'
    assert not output.exists();output.mkdir(mode=0o700)
    (output/'plan_before_original_node_interpretation.json').write_text(json.dumps(dict(
        created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),budget=budget,
        scope='All16 original T nodes in each disjoint RNA T view, all13 original myeloid nodes; both anchor programmes, all3 sources. 270 source/view/programme/node positions. Full original conditional view node counts and selected-cell identity localization.',
        limits='Scores trained on fixed100-cell caps. Full-view node fractions are separately reconstructed from all native marker counts or all original myeloid unit counts. Within-case centred node means locate participants; composition correlation uses equal case means. No P/clinical labels or identity promotion; broad-view expression can reflect both subtype mixtures and within-node states.'),indent=2)+'\n')
    dense.to_csv(output/'private_all2970_full_view_node_compositions.tsv',sep='\t',index=False)
    means.to_csv(output/'private_selected_case_node_score_means.tsv',sep='\t',index=False)
    rows=[]
    for view in ['R_FOXP3_capture','E_effector_capture','Myeloid']:
        cats=catalogue[catalogue.major.eq('Myeloid' if view=='Myeloid' else 'T')]
        for program in ['MCP1','MCP2']:
            for source_name in sorted(chosen.dataset.unique()):
                for node in cats.itertuples(index=False):
                    m=means[means.view.eq(view)&means.program.eq(program)&means.dataset.eq(source_name)&means.node_id.eq(node.node_id)]
                    f=dense[dense.view.eq(view)&dense.dataset.eq(source_name)&dense.node_id.eq(node.node_id)]
                    cs=case_scores[case_scores.view.eq(view)&case_scores.program.eq(program)&case_scores.dataset.eq(source_name)]
                    pair=f.merge(cs,on=['sample','dataset','view'],validate='one_to_one')
                    corr=float(spearmanr(pair.fraction,pair.score).statistic) if len(pair)>=3 and pair.fraction.nunique()>1 and pair.score.nunique()>1 else np.nan
                    rows.append(dict(dataset=source_name,program=program,view=view,node_id=node.node_id,
                        subcluster=node.subcluster,subtype_candidate=node.subtype_candidate,subtype_evidence_status=node.subtype_evidence_status,
                        identity_promoted=node.identity_promoted,source_units=len(f),full_view_units_with_node=int(f.full_node_cells.gt(0).sum()),
                        selected_units_with_node=len(m),selected_cells=int(m.selected_cells.sum()),
                        mean_within_case_centered_score=float(m.within_case_centered_score.mean()) if len(m) else np.nan,
                        full_node_fraction_vs_own_programme_rho=corr,descriptive_no_P=True))
    table=pd.DataFrame(rows);assert len(table)==270
    table.to_csv(output/'private_all270_source_original_node_localizations.tsv',sep='\t',index=False)
    files={p.name:sha(p) for p in sorted(output.iterdir()) if p.is_file()}
    (output/'private_terminal.json').write_text(json.dumps(dict(status='COMPLETE_ALL_ORIGINAL_NODE_LOCALIZATION_AND_FULL_COMPOSITION_DESCRIPTIONS_REVIEW_PENDING',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),files=files,
        full_case_view_node_records=2970,source_localizations=270,identity_promoted_by_this_analysis=False),indent=2)+'\n')
    print('ALL_ORIGINAL_NODE_AND_FULL_COMPOSITION_DESCRIPTIONS_DONE_REVIEW_PENDING')


if __name__=='__main__':main()
