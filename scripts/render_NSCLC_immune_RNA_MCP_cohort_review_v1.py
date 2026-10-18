"""Source-aggregate forest; no learned genes or individual records displayed."""
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


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,8*1024**2)
    audit=root/'NSCLC_immune_RNA_MCP_cohort_associations_audit_20261019_v1/private_terminal.json'
    assert sha(audit)=='ad5ca0c34c263577f981656beabcb03a482175521427b13453e5d07799120214'
    source=root/'NSCLC_immune_RNA_MCP_cohort_associations_20261019_v1';t=json.loads((source/'private_terminal.json').read_text())
    for name in ['private_all18_component_associations.tsv','private_all6_joint_cohort_clinical_bundles.tsv','private_all126_point_sensitivities.tsv']:assert sha(source/name)==t['files'][name]['sha256']
    out=root/'NSCLC_immune_RNA_MCP_cohort_source_review_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_all_source_review.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,
        scope='All6 source-aggregate bundles with all18 effects and95CI; all126 sensitivities retained in private complete review. No cases/genes/weights in figure; no primary-family change.',
        limits='Primary joint BH6/BY6 govern claims; nominal components or paired signals cannot replace the joint family. Author-case post-treatment cohort; no independent patient or causal mechanism assertion.'),indent=2)+'\n')
    components=pd.read_csv(source/'private_all18_component_associations.tsv',sep='\t')
    joint=pd.read_csv(source/'private_all6_joint_cohort_clinical_bundles.tsv',sep='\t').sort_values('slot')
    sens=pd.read_csv(source/'private_all126_point_sensitivities.tsv',sep='\t')
    display={'R_FOXP3_capture':'R','E_effector_capture':'E','Myeloid':'My'}
    labels=[r.program+' '+display[r.view_A]+' / '+display[r.view_B]+'  (n='+str(r.n_cases)+')' for r in joint.itertuples(index=False)]
    fig,axes=plt.subplots(1,4,figsize=(14,6.6),gridspec_kw={'width_ratios':[1.1,1,1,.7]},sharey=True)
    colors={'MCP1':'#286b97','MCP2':'#b05b24'};ys=np.arange(6)
    for ax,component,title in zip(axes[:3],['coordination','clinical_A','clinical_B'],['Programme coordination','A: pathology association','B: pathology association']):
        d=components[components.component.eq(component)].set_index('slot').loc[joint.slot]
        for y,r in zip(ys,d.itertuples()):
            color=colors[r.program];ax.plot([r.CI_low,r.CI_high],[y,y],color=color,lw=1.8);ax.scatter(r.effect,y,s=35,color=color,zorder=3)
        ax.axvline(0,color='#777777',ls='--',lw=.8);ax.set_title(title,fontsize=11,pad=14);ax.grid(axis='x',color='#eeeeee');ax.set_axisbelow(True)
        ax.set_xlabel('Partial rank correlation' if component=='coordination' else 'Response 1 minus 0\n(frozen discovery SD)',fontsize=9)
        for side in ['top','right']:ax.spines[side].set_visible(False)
    axes[0].set_yticks(ys,labels,fontsize=10);axes[0].invert_yaxis();axes[3].axis('off');axes[3].set_title('Joint BH / BY',fontsize=10,pad=14)
    for y,r in zip(ys,joint.itertuples(index=False)):axes[3].text(.04,y,f'{r.BH6:.3f} / {r.BY6:.3f}',va='center',fontsize=10,transform=axes[3].get_yaxis_transform())
    fig.suptitle('GSE243013 | Post-neoadjuvant NSCLC | Frozen six-bundle analysis',fontsize=14,y=.98)
    fig.text(.01,.035,'R: FOXP3-supported T capture; E: effector-supported T capture; My: author myeloid pool. These are not new formal cell identities.\nPoints: adjusted effects; bars: conditional 95% case-bootstrap CI. Joint P=max(coordination P, A P, B P); all six planned bundles retained.',fontsize=9,va='bottom')
    fig.subplots_adjust(left=.24,right=.985,top=.87,bottom=.19,wspace=.25)
    fig.savefig(out/'private_all6_joint_forest.png',dpi=180,facecolor='white');plt.close(fig)
    review=dict(reviewed_UTC=datetime.now(timezone.utc).isoformat(),dataset='GSE243013',source_sampling_context='post_neoadjuvant_surgery',unit='author_case',audit_sha256=sha(audit),
        all6_joint=joint.to_dict(orient='records'),all18_components=components.to_dict(orient='records'),all126_point_sensitivities=sens.to_dict(orient='records'),
        primary_joint_BH_below005=int(joint.BH6.lt(.05).sum()),minimum_joint_BH=float(joint.BH6.min()),
        interpretation='No primary joint candidate qualified. R/My MCP1 andMCP2 have directionally coherent nominal coassociation with all CI excluding0/LOO stable, but common-family BH about.134. Retain as exploratory not validated. Myeloid single-end clinical signals must not replace the two-end joint requirement. E-related programme coordination is not itself a clinical association.',
        limits='Fine-author-composition sensitivity has extreme-leverage failures especially R-containing models; no claim to have removed all subtype/state mixing. LUAD-only design often unavailable. Reused endpoints across overlapping case sets are not independent replication.',
        decision='Close fixed6 family without threshold/outcome/endpoint tuning. Qualify existing source/paired-material opportunities before further experiments. Stable-patient RNA background versus coordinated changes is a next biological question, not a retrospective primary-test replacement.')
    (out/'private_complete_family_scientific_review.json').write_text(json.dumps(review,indent=2)+'\n')
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(out.iterdir()) if p.is_file()}
    (out/'private_terminal.json').write_text(json.dumps(dict(status='COMPLETE_ALL6_SOURCE_REVIEW_FOREST_VISUAL_ACCEPTANCE_PENDING',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),files=files),indent=2)+'\n')
    print(json.dumps(dict(terminal_sha256=sha(out/'private_terminal.json'),review_sha256=sha(out/'private_complete_family_scientific_review.json'),figure_sha256=sha(out/'private_all6_joint_forest.png'))))


if __name__=='__main__':main()
