"""All fixed sources and both unit families, legible non-mechanistic forests."""
import hashlib
import itertools
import json
import os
from datetime import datetime,timezone
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha,VIEWS
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def short_source(value):
    if str(value).startswith('GSE'):return str(value)
    text=str(value)
    if text.startswith('CXG_'):return text[:12]
    return text[:14]+'_'+hashlib.sha256(text.encode()).hexdigest()[:4]


def short_phase(value):
    mapped={'pre_treatment':'pre','post_treatment':'post','on_treatment':'on','unknown_recorded_time':'time?','infusion_product':'product'}
    value=str(value)
    return mapped.get(value,value if len(value)<=15 else value[:10]+'_'+hashlib.sha256(value.encode()).hexdigest()[:4])


def material(value):
    mapped={'SOLID_TUMOR_CONTEXTS':'tumor','NORMAL_OR_ADJACENT':'normal','BLOOD':'blood','INFUSION_PRODUCT':'product','BONE_MARROW':'marrow','LYMPH_NODE':'node'}
    return mapped.get(str(value),'other')


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget=resource(root,32*1024**2)
    src=root/'NSCLC_external_programme_coordination_review_20261016_v3'
    proof=src/'private_independent_context_review_audit.json';audit=json.loads(proof.read_text())
    assert audit['status'].startswith('PASS_ALL630_RNA_BACKGROUND_')
    assert audit['review_terminal_sha256']==sha(src/'private_terminal.json')
    terminal=json.loads((src/'private_terminal.json').read_text())
    p=src/'private_all_primary_and_sensitivity_review.tsv';assert sha(p)==terminal['files'][p.name]['sha256']
    full=pd.read_csv(p,sep='\t');assert len(full)==630 and full.record.is_unique
    out=root/'NSCLC_external_coordination_forests_20261018_v2';assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),context_review_audit_sha256=sha(proof),input_table_sha256=sha(p),budget=budget,
              scope='All10 factors x2 existing unit families x6 pairs, all630 fixed contexts; plotting-only, no model/selection/FDR changes.',
              symbols='Star: BH<.05 and CI excludes zero. Hollow: CI unavailable under1900 valid bootstrap gate. Blue: all active units reviewed, grey: sample/library fallback.',
              layout='3 rows x2 columns; source/phase/material short labels with private exact mapping. All text artist bounds checked; repeated GSE contexts not independent replicates.')
    (out/'plan_before_reformatting_existing_forests.json').write_text(json.dumps(plan,indent=2)+'\n')
    source_map={value:short_source(value) for value in full.dataset.unique()}
    assert len(set(source_map.values()))==len(source_map)
    short={VIEWS[0]:'CD8',VIEWS[1]:'Treg',VIEWS[2]:'Mac',VIEWS[3]:'B/Plasma'}
    pairs=list(itertools.combinations(VIEWS,2));figures=[];manifest=[];text_checked=0
    for factor in sorted(full.factor.unique()):
        for family in ['reviewed_patients','all_units']:
            fig,axes=plt.subplots(3,2,figsize=(16,14),constrained_layout=True)
            part=full[full.factor.eq(factor)&full.family.eq(family)]
            for ax,(va,vb) in zip(axes.flat,pairs):
                g=part[part.view_A.eq(va)&part.view_B.eq(vb)].sort_values(['material_group','dataset','phase_key','arm','context_id'])
                labels=[];base=[]
                for row in g.itertuples():base.append(f'{source_map[row.dataset]}  {short_phase(row.phase_key)}  {material(row.material_group)}  n={row.active_units}')
                for j,row in enumerate(g.itertuples()):
                    label=base[j]
                    if base.count(label)>1:label+='  '+hashlib.sha256(row.context_id.encode()).hexdigest()[:4]
                    labels.append(label);known=row.active_units==row.reviewed_units;color='#2166ac' if known else '#777777'
                    has_ci=bool(np.isfinite(row.CI_low) and np.isfinite(row.CI_high))
                    if has_ci:ax.plot([row.CI_low,row.CI_high],[j,j],color=color,linewidth=1.3)
                    supported=has_ci and row.BH_q<.05 and (row.CI_low>0 or row.CI_high<0)
                    if supported:ax.scatter([row.partial_rank_rho],[j],marker='*',color=color,s=85,zorder=3)
                    elif has_ci:ax.scatter([row.partial_rank_rho],[j],marker='o',color=color,s=25,zorder=3)
                    else:ax.scatter([row.partial_rank_rho],[j],marker='o',facecolors='none',edgecolors=color,s=32,zorder=3)
                    manifest.append(dict(factor=factor,family=family,view_A=va,view_B=vb,record=row.record,dataset=row.dataset,phase_key=row.phase_key,material_class=row.material_class,material_group=row.material_group,display_label=label,active_units=row.active_units,CI_available=has_ci,star=supported))
                ax.set_yticks(np.arange(len(g)),labels,fontsize=8);ax.invert_yaxis();ax.axvline(0,color='black',linewidth=.7)
                ax.set_xlim(-1.05,1.05);ax.set_xticks([-1,-.5,0,.5,1]);ax.set_title(short[va]+' / '+short[vb],fontsize=11)
                ax.set_xlabel('Partial rank correlation',fontsize=9);ax.grid(axis='x',alpha=.14)
                if g.empty:ax.text(.5,.5,'No estimable context under fixed gates',ha='center',va='center',transform=ax.transAxes,fontsize=9)
            label='Reviewed-patient subset' if family=='reviewed_patients' else 'All-unit exploratory family'
            fig.suptitle(factor+' | '+label+'\nFixed gene-programme co-variation, not direct communication. Star: BH + CI support; hollow: CI unavailable.\nBlue: reviewed units; grey: library fallback. Multiple contexts of one GSE are not independent sources.',fontsize=11)
            fig.canvas.draw();renderer=fig.canvas.get_renderer();width,height=fig.bbox.width,fig.bbox.height
            for ax in axes.flat:
                for artist in [*ax.get_yticklabels(),*ax.get_xticklabels(),ax.title,ax.xaxis.label]:
                    if not artist.get_text():continue
                    bb=artist.get_window_extent(renderer)
                    assert bb.x0>=-1 and bb.y0>=-1 and bb.x1<=width+1 and bb.y1<=height+1,'Plot text outside figure bounds'
                    text_checked+=1
            name='private_'+factor+'_'+family+'_all_context_forests_v2.png';fig.savefig(out/name,dpi=140,bbox_inches='tight');plt.close(fig);figures.append(name)
    manifest=pd.DataFrame(manifest);assert len(manifest)==630 and manifest.record.is_unique and set(manifest.record)==set(full.record)
    manifest.to_csv(out/'private_all630_exact_plot_labels.tsv.gz',sep='\t',index=False,compression='gzip')
    report=dict(status='COMPLETE_ALL20_FORESTS_ALL630_CONTEXTS_LABEL_BOUNDS_PASS_VISUAL_REVIEW_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),context_review_audit_sha256=sha(proof),figures=figures,artist_bounds_checked=text_checked,
                files={f.name:dict(sha256=sha(f),bytes=f.stat().st_size) for f in out.iterdir() if f.is_file()},statistical_values_changed=False,clinical_or_direct_communication_validation=False)
    assert sum(item['bytes'] for item in report['files'].values())<budget['output_cap_bytes']
    (out/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],finished_UTC=report['finished_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
