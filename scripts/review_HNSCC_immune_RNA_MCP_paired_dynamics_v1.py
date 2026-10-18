"""All12 source forest/review and all6 descriptive cross-lag covariance identities."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from HNSCC_immune_RNA_MCP_paired_io_v1 import read
from HNSCC_immune_RNA_MCP_paired_math_v1 import design


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,8*1024**2)
    audit=root/'HNSCC_immune_RNA_MCP_paired_dynamics_audit_20261019_v1/private_terminal.json';assert sha(audit)=='99c0eae3ee66bb591a3a1a2460ce49ff3438abdcaae8847596ce9d9f9889cdeb'
    source=root/'HNSCC_immune_RNA_MCP_paired_dynamics_20261019_v1';t=json.loads((source/'private_terminal.json').read_text())
    for n,v in t['files'].items():assert sha(source/n)==v['sha256']
    _,_,plan,_,arrays=read(root)
    out=root/'HNSCC_immune_RNA_MCP_paired_source_review_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_all_source_review_and_covariance_diagnostics.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,
        scope='All12 primary results and48points, full forest. Post-result descriptive covariance decomposition for all6 pairs;24 view/phase/quality training-domain summaries. No new P/CI or primary family change.',
        diagnostic='Scores use original discovery SD. Apply SAME frozen primary nuisance projector to pre/post score vectors; Cov(deltaA,deltaB)=Cov(postA,postB)+Cov(preA,preB)-Cov(postA,preB)-Cov(preA,postB). Raw-score covariance identity, NOT a decomposition of primary partial rank rho, nor causal attribution.',
        reason='Some baselines participated in programme learning; delta covariance can inherit preexisting covariance. Treatment/phase/procedure and discovery overlap require separate checks.'),indent=2)+'\n')
    result=pd.read_csv(source/'private_all12_paired_dynamics_associations.tsv',sep='\t');sens=pd.read_csv(source/'private_all48_point_sensitivities.tsv',sep='\t')
    shifts={(r['program'],r['view']):arrays[r['slot']] for r in plan['slots'] if r['kind']=='paired_median_shift'}
    diagnostics=[]
    for item in plan['slots']:
        if item['kind']!='paired_cochange':continue
        a=arrays[item['slot']];m,status=design(a['x'],a['y'],a['C'],a['names']);r=dict(slot=item['slot'],program=item['program'],view_A=item['view_A'],view_B=item['view_B'],n_patients=len(a['x']),status=status)
        if m is not None:
            sa=shifts[(item['program'],item['view_A'])];sb=shifts[(item['program'],item['view_B'])]
            ai=pd.Index(sa['patients']).get_indexer(a['patients']);bi=pd.Index(sb['patients']).get_indexer(a['patients']);assert (ai>=0).all() and (bi>=0).all()
            M=np.eye(len(ai))-m['Q']@m['Q'].T
            apre=M@sa['pre'][ai];apost=M@sa['post'][ai];bpre=M@sb['pre'][bi];bpost=M@sb['post'][bi];n=len(ai)
            def cov(x,y):return float(x@y/(n-1))
            def corr(x,y):return float(x@y/np.linalg.norm(x)/np.linalg.norm(y)) if min(np.linalg.norm(x),np.linalg.norm(y))>1e-12 else np.nan
            c00=cov(apre,bpre);c11=cov(apost,bpost);c10=cov(apost,bpre);c01=cov(apre,bpost);delta=cov(apost-apre,bpost-bpre);identity=c11+c00-c10-c01
            assert abs(delta-identity)<1e-9
            r.update(pre_pre_cov=c00,post_post_cov=c11,post_pre_cov=c10,pre_post_cov=c01,delta_cov=delta,identity_error=abs(delta-identity),same_projector_pre_rho=corr(apre,bpre),same_projector_post_rho=corr(apost,bpost),raw_delta_rho=corr(apost-apre,bpost-bpre),new_P_or_CI=False)
        diagnostics.append(r)
    diag=pd.DataFrame(diagnostics);assert len(diag)==6;diag.to_csv(out/'private_all6_descriptive_covariance_identities.tsv',sep='\t',index=False)
    assay=root/'HNSCC_immune_RNA_MCP_paired_native_measurement_20261019_v1';cell=pd.read_csv(assay/'private_selected_cells_and_all_scores.tsv.gz',sep='\t');material=pd.read_csv(assay/'private_all52_paired_materials.tsv',sep='\t');cell['phase']=material.iloc[cell.material_index.to_numpy(int)].phase.to_numpy();domains=pd.read_csv(assay/'private_all12_training_domain_descriptions.tsv',sep='\t').set_index(['view','covariate']);views=['R_FOXP3_capture','E_effector_capture','Myeloid'];domainrows=[]
    for (vi,phase),g in cell.groupby(['view_index','phase']):
        view=views[vi]
        for f in ['cellQ','geneQ','fraction','logN']:
            d=domains.loc[(view,f)];v=g[f].to_numpy();domainrows.append(dict(view=view,phase=phase,covariate=f,selected_cells=len(g),selected_cells_outside_training=int(((v<d.training_min)|(v>d.training_max)).sum()),selected_cell_median=float(np.median(v)),training_min=float(d.training_min),training_max=float(d.training_max),new_P_or_CI=False))
    dr=pd.DataFrame(domainrows);assert len(dr)==24;dr.to_csv(out/'private_all24_phase_training_domain_descriptions.tsv',sep='\t',index=False)
    display={'R_FOXP3_capture':'R','E_effector_capture':'E','Myeloid':'My'};colors={'MCP1':'#286b97','MCP2':'#b05b24'}
    fig,axes=plt.subplots(1,2,figsize=(13,6.6))
    for ax,kind,title in zip(axes,['paired_median_shift','paired_cochange'],['Paired programme change','Paired programme cochange']):
        d=result[result.kind.eq(kind)].sort_values('slot');ys=np.arange(6)
        labels=[]
        for y,r in zip(ys,d.itertuples(index=False)):
            label=r.program+' '+(display[r.view] if kind=='paired_median_shift' else display[r.view_A]+' / '+display[r.view_B]);labels.append(label+'  q='+f'{r.BH12:.3g}')
            ax.plot([r.CI_low,r.CI_high],[y,y],color=colors[r.program],lw=2);ax.scatter(r.effect,y,s=35,color=colors[r.program],zorder=3)
        ax.set_yticks(ys,labels,fontsize=10);ax.invert_yaxis();ax.axvline(0,color='#777777',ls='--',lw=.8);ax.set_title(title,fontsize=12,pad=12);ax.grid(axis='x',color='#eeeeee');ax.set_axisbelow(True)
        ax.set_xlabel('Median post minus pre\n(original discovery case SD)' if kind=='paired_median_shift' else 'Adjusted partial rank correlation\nof paired programme changes',fontsize=10)
        for side in ['top','right']:ax.spines[side].set_visible(False)
    fig.suptitle('GSE286827 | 26 reviewed paired patients | Frozen 12-question family',fontsize=14,y=.98)
    fig.text(.015,.025,'R: FOXP3-supported T capture; E: effector-supported T capture; My: original myeloid pool. No new formal identities.\nLeft: exact sign tests, distribution-free >=95% median intervals. Right: HC3 wild tests, 95% patient-bootstrap intervals.\nAll 12 tests share BH/BY correction. Baseline/patient overlap with discovery; time/procedure/therapy are confounded.',fontsize=9,va='bottom')
    fig.subplots_adjust(left=.16,right=.985,top=.86,bottom=.23,wspace=.8);fig.savefig(out/'private_all12_paired_forest.png',dpi=180,facecolor='white');plt.close(fig)
    review=dict(reviewed_UTC=datetime.now(timezone.utc).isoformat(),dataset='GSE286827',audit_sha256=sha(audit),all12_primary=result.to_dict(orient='records'),all48_sensitivities=sens.to_dict(orient='records'),all6_covariance_diagnostics=diag.to_dict(orient='records'),all24_phase_domain_descriptions=dr.to_dict(orient='records'),
        interpretation='MCP2 R and E paired shifts pass BH/BY12. R/E cochange passes BH/BY andMCupperBH12, but bootstrap CI spans0: no prespecified priority cochange candidate. Retain source-internal T-associated temporal programme clue; not validated communication, formal cell state, clinical response prediction or causal therapy effect.',
        robustness='All-capture/stress-complement R/E cochange point directions retained; original-node composition point uses25 cases due unseen node. No new sensitivity P or replacement of the main family.',
        next='Use descriptive pre/post/cross-lag covariance and phase-domain summaries to distinguish preexisting discovery covariance and processing context before any mechanistic follow-up. Do not retest alternative P/CI to bypass failed priority gate.')
    (out/'private_complete_paired_scientific_review.json').write_text(json.dumps(review,indent=2)+'\n')
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(out.iterdir()) if p.is_file()}
    t=dict(status='COMPLETE_ALL12_PAIRED_SOURCE_REVIEW_FOREST_AND_DESCRIPTIVE_DIAGNOSTICS_VISUAL_PENDING',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),files=files)
    (out/'private_terminal.json').write_text(json.dumps(t,indent=2)+'\n');print(json.dumps(dict(status=t['status'],terminal_sha256=sha(out/'private_terminal.json'),figure_sha256=sha(out/'private_all12_paired_forest.png'))))


if __name__=='__main__':main()
