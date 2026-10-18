"""Private complete source-level review and all400 functional forest panels."""
import json
import os
import textwrap
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from immune_RNA_MCP_new_source_io_v1 import wait


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,12*1024**2);out=root/'GSE176021_disjoint_function_complete_review_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    audit=root/'GSE176021_disjoint_function_association_audit_20261019_v1/private_terminal.json'
    wait(root,out,'audit_GSE176021_disjoint_function_family_v1',audit,{'ACCEPTED_ALL400_FUNCTIONAL_RESULTS3200_SENSITIVITIES_COMPLETE_RANDOM_STREAM_REPLAY'})
    src=root/'GSE176021_disjoint_function_associations_20261019_v1'
    (out/'plan_before_complete_private_functional_review.json').write_text(json.dumps({'created_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'resource_guard':guard,'association_acceptance_sha256':sha(audit),'scope':'Read every400 main and3200sensitivity positions. Fixed primary criteria unchanged. Four figures, each all50 references/two response lineages for oneprogramme/predictor lineage; no effect-selected display subset. Scientific values/plots private; only source aggregates authorized for chat.','limits':'Exploratory single-source functional complements, not binding, formalstate or clinicalprediction. All95percent CIs marginal, BH/BY400 and finite bootstrap/LOO/MCcriteria separated. Response quality controls only point sensitivity atN10.'},indent=2)+'\n')
    result=pd.read_csv(src/'private_all400_primary_functional_results.tsv',sep='\t');sensitivity=pd.read_csv(src/'private_all3200_point_sensitivity_results.tsv',sep='\t');assert len(result)==400 and len(sensitivity)==3200
    q=sensitivity[sensitivity['mode'].eq('response_quality_partial_two')].set_index('model_id');assert len(q)==400
    result['response_quality_partial_point']=result.model_id.map(q.effect);result['response_quality_point_status']=result.model_id.map(q.status)
    result['quality_point_retains_sign']=result.effect*result.response_quality_partial_point>0
    result['interpretation_tier']=np.select([result.priority_functional_association & result.role.eq('CROSS_LINEAGE'),result.priority_functional_association & result.role.eq('WITHIN_LINEAGE'),result.BH.lt(.05),result.effect.notna()],['CROSS_LINEAGE_FUNCTIONAL_CANDIDATE','WITHIN_LINEAGE_FUNCTIONAL_CANDIDATE','BH_PASS_STABILITY_INSUFFICIENT','NO_FAMILY_ADJUSTED_ASSOCIATION'],default='NOT_ESTIMABLE')
    summary=[]
    for (program,role),d in result.groupby(['program','role'],sort=True):
        summary.append({'program':program,'role':role,'positions':len(d),'estimable':int(d.effect.notna().sum()),'BH_below05':int(d.BH.lt(.05).sum()),'BY_below05':int(d.BY.lt(.05).sum()),'fixed_priority':int(d.priority_functional_association.sum()),'priority_quality_available':int((d.priority_functional_association & d.response_quality_partial_point.notna()).sum()),'priority_quality_same_sign':int((d.priority_functional_association & d.quality_point_retains_sign).sum())})
    rank=result.sort_values(['priority_functional_association','BH','BY'],ascending=[False,True,True],na_position='last');rank.to_csv(out/'private_complete400_ranked_with_quality_points.tsv',sep='\t',index=False)
    science={'created_UTC':datetime.now(timezone.utc).isoformat(),'source':'GSE176021','cancer':'NSCLC','phase':'current_material_metadata_unknown; source_protocol_neoadjuvant_PD1_resection_not_case-level_adoption','case_units':10,'family_denominator':400,'summary_all_groups':summary,'all400_primary':json.loads(result.to_json(orient='records',double_precision=15)),'all3200_point_sensitivities':json.loads(sensitivity.to_json(orient='records',double_precision=15)),'limits':['All pathways are predictor-gene-excluded complements, not unmodified canonical signatures.','Within/cross lineage associations not independent mechanisms: MCP predictors already cochange.','Original frozen correction does not replace complete current-case quality/abundance control.','Normal reference, shared upstream and treatment environment remain possible.','No MPR case bridge, no binding or clinicalprediction inference.','Family is closed after complete review; no significance-driven pathway/control changes.']}
    (out/'private_complete_source_functional_review.json').write_text(json.dumps(science,indent=2)+'\n')
    references=sorted(result.reference.unique());views=['R_FOXP3_capture','E_effector_capture'];names={'R_FOXP3_capture':'author Treg','E_effector_capture':'author CD8 union'}
    for program in ['MCP1','MCP2']:
        for pv in views:
            fig,axes=plt.subplots(1,2,figsize=(15,16),sharey=True,layout='constrained')
            for ax,rv in zip(axes,views):
                d=result[result.program.eq(program)&result.predictor_view.eq(pv)&result.response_view.eq(rv)].set_index('reference').loc[references];assert len(d)==50
                for j,(_,r) in enumerate(d.iterrows()):
                    if not np.isfinite(r.effect):ax.text(-1.0,j,'not estimable',fontsize=6,color='grey');continue
                    color='#22699e' if r.priority_functional_association else '#b45d22' if r.BH<.05 else '#777777'
                    if np.isfinite(r.CI_low) and np.isfinite(r.CI_high):ax.plot([r.CI_low,r.CI_high],[j,j],color=color,lw=.8,alpha=.8)
                    ax.plot(r.effect,j,'o',ms=4,mec=color,mfc=color if r.priority_functional_association else 'white')
                ax.axvline(0,color='black',lw=.6);ax.set_xlim(-1.05,1.05);ax.set_ylim(49.8,-.8);ax.set_xticks([-1,-.5,0,.5,1]);ax.grid(axis='x',alpha=.15);ax.set_title('Response: '+names[rv]+'\n'+('within lineage' if pv==rv else 'cross lineage'));ax.set_xlabel('Spearman correlation of paired tissue changes')
            axes[0].set_yticks(np.arange(50));axes[0].set_yticklabels(['\n'.join(textwrap.wrap(x.removeprefix('HALLMARK_').replace('_',' '),28)) for x in references],fontsize=6.8)
            fig.suptitle('GSE176021 | '+program+' predictor: '+names[pv]+' | all50 functional complements\n10 shared case pairs; one complete400-test family',fontsize=12)
            fig.supxlabel('Lines: marginal95% bootstrap CI. Filled: all fixed priority criteria. Orange open: BH<0.05 only. Grey: other.\nComplementary functional association does not establish receptor binding or clinical prediction.',fontsize=8)
            stem='private_'+program+'_predictor_'+('R' if pv==views[0] else 'E')+'_all100_functional_forests'
            fig.savefig(out/(stem+'.png'),dpi=180);fig.savefig(out/(stem+'.svg'));plt.close(fig)
    files={p.name:{'sha256':sha(p),'bytes':p.stat().st_size} for p in out.iterdir() if p.is_file()};assert sum(v['bytes'] for v in files.values())<12*1024**2
    receipt={'status':'COMPLETE_ALL400_PRIVATE_FUNCTIONAL_REVIEW_AND4_FULL_FOREST_FIGURES_VISUAL_PENDING','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'association_acceptance_sha256':sha(audit),'files':files,'main_positions':400,'point_sensitivities':3200,'forest_panels':8,'no_individual_scatter_or_gene_weights_exported':True}
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'status':receipt['status'],'terminal_sha256':sha(out/'private_terminal.json')}))


if __name__=='__main__':main()
