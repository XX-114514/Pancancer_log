"""All-function biological interpretation, with effect-blind display contexts."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha, VIEWS
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource

SOLID={'RECORDED_PRIMARY_TUMOR','RECORDED_SOLID_TUMOR','RECORDED_SOLID_METASTASIS','TUMOR_CAPTURE_SITE_UNSPECIFIED'}


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget=resource(root,16*1024**2)
    src=root/'NSCLC_external_disjoint_function_coupling_20261016_v1'
    proof=src/'private_independent_functional_coupling_audit.json'
    assert sha(proof)=='4dcc2e3a1378e0d55b99e316c2b29cfdc5b5f69e0bb87f200f89854e74bd24ef'
    assert json.loads(proof.read_text())['producer_terminal_sha256']==sha(src/'private_terminal.json')
    out=root/'NSCLC_external_disjoint_function_review_20261016_v1';assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),budget=budget,
              scope='All50 functions/all10 factors/all four compatible pools. Complete table and all-material descriptions retained; displayed source summaries use reviewed-patient solid-tumor contexts only.',
              display_context='One record per GSE/factor/own pool, selected by pre_treatment, unknown_recorded_time, on_treatment, post_treatment, other phases; lexical phase, largest active cases, lexical slot. Partner-QC context is retained explicitly; no choice by effect or function.',
              summary='Across selected distinct GSE: signed median/min/max ULM, target and nontarget partial ranks, source counts, and same-orientation ULM plus target-mean |rho|>=.3 with non-low support counts. These are descriptive, no P or independent validation.',
              display='All50 references in each four-pool heatmap; colours only if >=2 finite sources, single-source values retained in table. No source/context is an untouched independent validation by default.',
              main_coordination_P_or_clinical_values_loaded=False,
              limitations='Shared samples/denominators, phase/site/sorting and cell-mixture differences remain. Not formal state annotation, exact lung latent-factor replication, pathway activation proof or direct communication.')
    (out/'plan_before_any_functional_effect_review.json').write_text(json.dumps(plan,indent=2)+'\n')
    table=pd.read_csv(src/'private_all_fixed_disjoint_function_couplings.tsv.gz',sep='\t')
    assert len(table)==31500 and not table.duplicated(['record','reference']).any()
    expanded=[]
    for side,other in [('A','B'),('B','A')]:
        x=table.copy();x['own_view']=x['view_'+side];x['partner_view']=x['view_'+other]
        for key in ['ULM_score','target_mean','nontarget_mean']:
            x[key+'_own_partial_rho']=x[key+'_programme'+side+'_function'+side+'_partial_rho']
        x['own_low_foreground_support']=x['low_foreground_support_'+side]
        expanded.append(x)
    full=pd.concat(expanded,ignore_index=True)
    # Keep all materials/families in a separate complete synopsis.
    all_scope=full.groupby(['family','material_class','dataset','phase_key','factor','own_view','reference'],dropna=False).agg(
        pair_contexts=('record','nunique'),descriptive_median_ULM_rho=('ULM_score_own_partial_rho','median'),
        minimum_ULM_rho=('ULM_score_own_partial_rho','min'),maximum_ULM_rho=('ULM_score_own_partial_rho','max')).reset_index()
    all_scope.to_csv(out/'private_all_material_phase_function_synopsis.tsv.gz',sep='\t',index=False,compression='gzip',na_rep='NA')
    eligible=full[full.family.eq('reviewed_patients') & full.material_class.isin(SOLID)].copy()
    priority={'pre_treatment':0,'unknown_recorded_time':1,'on_treatment':2,'post_treatment':3}
    records=eligible.drop_duplicates(['record','own_view']).copy()
    records['phase_priority']=records.phase_key.map(priority).fillna(4)
    records=records.sort_values(['phase_priority','phase_key','active_units','slot'],ascending=[True,True,False,True])
    chosen=records.drop_duplicates(['dataset','factor','own_view'])[['record','own_view']]
    selected=eligible.merge(chosen,on=['record','own_view'],how='inner',validate='many_to_one')
    assert selected.groupby(['dataset','factor','own_view']).reference.nunique().eq(50).all()
    selected.to_csv(out/'private_effect_blind_display_source_function_rows.tsv.gz',sep='\t',index=False,compression='gzip',na_rep='NA')
    summaries=[]
    for levels,g in selected.groupby(['factor','own_view','reference'],sort=True):
        factor,view,reference=levels
        assert g.dataset.is_unique
        row=dict(factor=factor,view=view,reference=reference,selected_sources=len(g),
                 low_support_sources=int(g.own_low_foreground_support.sum()))
        for key in ['ULM_score','target_mean','nontarget_mean']:
            value=g[key+'_own_partial_rho'];row[key+'_finite_sources']=int(value.notna().sum())
            for tag,stat in [('median',value.median()),('minimum',value.min()),('maximum',value.max())]:row[key+'_'+tag+'_rho']=float(stat)
        good=~g.own_low_foreground_support
        row['positive_ULM_and_target_sources']=int((good & g.ULM_score_own_partial_rho.ge(.3) & g.target_mean_own_partial_rho.ge(.3)).sum())
        row['negative_ULM_and_target_sources']=int((good & g.ULM_score_own_partial_rho.le(-.3) & g.target_mean_own_partial_rho.le(-.3)).sum())
        summaries.append(row)
    summary=pd.DataFrame(summaries)
    assert not summary.duplicated(['factor','view','reference']).any()
    summary.to_csv(out/'private_all_source_disjoint_function_summaries.tsv',sep='\t',index=False,na_rep='NA')
    figures=[];names=sorted(table.reference.unique());factors=sorted(table.factor.unique())
    for vi,view in enumerate(VIEWS):
        values=summary[summary.view.eq(view)].pivot(index='reference',columns='factor',values='ULM_score_median_rho').reindex(index=names,columns=factors)
        counts=summary[summary.view.eq(view)].pivot(index='reference',columns='factor',values='ULM_score_finite_sources').reindex(index=names,columns=factors)
        matrix=values.to_numpy().copy();matrix[counts.to_numpy()<2]=np.nan
        fig,axes=plt.subplots(1,2,figsize=(13,15),constrained_layout=True,gridspec_kw={'width_ratios':[3,2]})
        im=axes[0].imshow(matrix,aspect='auto',vmin=-1,vmax=1,cmap='RdBu_r')
        axes[0].set_yticks(np.arange(50),[n.removeprefix('HALLMARK_') for n in names],fontsize=8)
        axes[0].set_xticks(np.arange(10),[f.rsplit('_',1)[-1] for f in factors],fontsize=9)
        axes[0].set_title('Median partial rho: own programme / disjoint function\nColours require at least two measured GSE',fontsize=10)
        fig.colorbar(im,ax=axes[0],shrink=.35,label='Descriptive rho, no P')
        cc=axes[1].imshow(counts.to_numpy(),aspect='auto',vmin=0,vmax=max(4,float(np.nanmax(counts.to_numpy()))),cmap='Blues')
        axes[1].set_yticks(np.arange(50),['']*50)
        axes[1].set_xticks(np.arange(10),[f.rsplit('_',1)[-1] for f in factors],fontsize=9)
        axes[1].set_title('Distinct measured sources\nReviewed-patient solid-tumor contexts',fontsize=10)
        fig.colorbar(cc,ax=axes[1],shrink=.35,label='Source count; not proof of independence')
        fig.suptitle(view+' | all50 fixed functions, genes excluded from all programme views\nPhase/partner selected using metadata only; primary coordination results not read',fontsize=12)
        name=f'private_pool_{vi}_all50_disjoint_function_heatmap.png';fig.savefig(out/name,dpi=145);plt.close(fig);figures.append(name)
    terminal=dict(status='COMPLETE_ALL_DISJOINT_FUNCTION_SOURCE_SYNOPSIS_AND_FOUR_FULL_MAPS_VISUAL_REVIEW_PENDING',
                  finished_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),coupling_audit_sha256=sha(proof),figures=figures,
                  files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()},
                  new_P_or_primary_coordination_or_clinical_results_loaded=False,clinical_or_direct_communication_validation=False)
    assert sum(v['bytes'] for v in terminal['files'].values())<budget['output_cap_bytes']
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'],finished_UTC=terminal['finished_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
