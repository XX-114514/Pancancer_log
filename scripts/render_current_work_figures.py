#!/usr/bin/env python3
"""Render current-work figures from published aggregate evidence only.

No new estimates, tests, fitted models or cell-level exports. Run from any cwd.
The separate render_current_work_umap.py renders server-only embeddings.
"""
from pathlib import Path
import hashlib
import json
import platform

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'assets/publication/20261105_current_work'
OUT.mkdir(parents=True, exist_ok=True)
AUD = ROOT / 'inventories/server_analysis/2026-11-05'
BLUE, GREEN, ORANGE, GREY = '#0072B2', '#009E73', '#D55E00', '#8b95a3'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8,
    'axes.titlesize': 10, 'axes.labelsize': 8, 'xtick.labelsize': 7,
    'ytick.labelsize': 7, 'axes.spines.top': False, 'axes.spines.right': False,
    'axes.linewidth': .65, 'pdf.fonttype': 42, 'savefig.dpi': 260})
SOURCES = set()
FIGURES = {}


def table(path, **kw):
    path = Path(path); SOURCES.add(path)
    return pd.read_csv(path, sep='\t', **kw)


def record(path):
    path = Path(path); SOURCES.add(path)
    return json.loads(path.read_text())


def save(fig, name, sources, scope):
    for ext in ['png', 'pdf']:
        fig.savefig(OUT / f'{name}.{ext}', bbox_inches='tight')
    FIGURES[name] = {'scope': scope, 'sources': [str(Path(p).relative_to(ROOT)) for p in sources]}
    plt.close(fig)


def heading(ax, text):
    ax.set_title(text, loc='left', fontweight='bold', pad=12)


def datasets():
    cp = ROOT / 'inventories/integration_v2/composition_by_cancer.tsv'
    bp = ROOT / 'inventories/integration_v2/composition_by_batch.tsv'
    vp = AUD / 'version_crosswalk.tsv'
    c, b, versions = table(cp), table(bp), table(vp)
    assert c.cells.sum() == b.cells.sum() == 4720931 and b.sources.sum() == 72
    fig = plt.figure(figsize=(12, 9)); gs = fig.add_gridspec(2, 2, width_ratios=[1.4, 1], hspace=.6, wspace=.60)
    ax = fig.add_subplot(gs[:, 0]); y = np.arange(len(c))
    cmap = plt.cm.viridis; norm = Normalize(1, 9)
    ax.barh(y, c.cells / 1e6, color=cmap(norm(c.sources)))
    ax.set_yticks(y, c.analysis_main); ax.invert_yaxis(); ax.set_xlim(0, .73)
    for i, row in c.iterrows(): ax.text(row.cells / 1e6 + .009, i, f'{row.sources} src', va='center', fontsize=7)
    ax.set_xlabel('Cells (millions); labels give source records')
    heading(ax, 'a  integration_v2: 28 cancer strata')
    ax = fig.add_subplot(gs[0, 1]); y = np.arange(4)
    ax.barh(y, b.cells / 1e6, color=[BLUE, '#56B4E9', GREEN, ORANGE])
    ax.set_yticks(y, ['Existing set', 'Supplement v1', 'Supplement v2', 'Supplement v3']); ax.invert_yaxis()
    ax.set_xlim(0, 3.35)
    for i, row in b.iterrows(): ax.text(row.cells / 1e6 + .04, i, f'{row.cells:,}\n{row.sources} sources', va='center', fontsize=8)
    ax.set_xlabel('Cells (millions)'); heading(ax, 'b  Expanded data layer: 72 sources')
    ax = fig.add_subplot(gs[1, 1])
    v = versions.set_index('version')
    import re
    expanded = int(re.match(r'(\d+)_cells', v.loc['integration_v2_20261102', 'universe']).group(1))
    tme, eligible = map(int, re.match(r'(\d+)_TME_(\d+)_CoVarNet_eligible', v.loc['TME_current_scope_20261105', 'universe']).groups())
    counts = [expanded, tme, eligible]
    ax.barh(range(3), np.array(counts) / 1e6, color=[GREY, BLUE, GREEN])
    ax.set_yticks(range(3), ['integration_v2', 'Current TME scope', 'CoVarNet-eligible']); ax.invert_yaxis()
    for i, n in enumerate(counts): ax.text(.08, i, f'{n:,}', va='center', color='white', fontweight='bold')
    ax.set_xlim(0, 5); ax.set_xlabel('Cells (millions)')
    heading(ax, 'c  Data layer and current analysis eligibility')
    fig.text(.52, .02, 'Counts do not imply final TME run_v2 acceptance.\nPatient keys and analysis units are not automatically independent patients.', fontsize=8)
    save(fig, 'fig1_data_scope', [cp, bp, vp], 'integration_v2 plus current TME/CoVarNet scope; not 0906 or CT75')


def batch():
    gp, sp = AUD/'tables/batch_benchmark/graph_metrics.tsv', AUD/'tables/batch_benchmark/scib_metrics.tsv'
    g, s = table(gp, index_col=0), table(sp, index_col=0)
    order = ['unintegrated','v1_fulldata_harmony','harmony_source_unit','harmony_source_t2','bbknn_source','harmony_bbknn','scvi_source','scvi_unit','scanorama','combat']
    names = ['Unintegrated','Original full-data Harmony','Harmony source + unit','Harmony source theta2','BBKNN source','Harmony + BBKNN','scVI source','scVI unit','Scanorama','ComBat']
    colors = [GREEN if k == 'harmony_source_t2' else BLUE for k in order]
    fig, axes = plt.subplots(1, 3, figsize=(12, 5), gridspec_kw={'width_ratios':[1,1,1.2]})
    fig.subplots_adjust(wspace=.3, left=.22, right=.99, bottom=.16, top=.78)
    y = np.arange(len(order))
    for ax, field, title, xmax in zip(axes[:2], ['source_iLISI_median','purity_fine'], ['a  Source mixing','b  Marker-label retention'], [8.5,1]):
        vals = g.loc[order, field]
        ax.barh(y, vals, color=colors, height=.6); ax.invert_yaxis(); ax.set_xlim(0,xmax)
        ax.set_yticks(y, names if ax == axes[0] else ['']*len(y))
        for i, val in enumerate(vals): ax.text(val+.02*xmax,i,f'{val:.3f}',va='center',fontsize=7)
        heading(ax,title)
    axes[0].set_xlabel('Graph source iLISI (median)'); axes[1].set_xlabel('Fine-marker neighborhood purity')
    ax=axes[2]
    for i,k in enumerate(order):
        if k in s.index:
            bio=float(s.loc[k,'Bio conservation']); bat=float(s.loc[k,'Batch correction']); total=float(s.loc[k,'Total'])
            assert abs(total-(.6*bio+.4*bat))<1e-5
            ax.barh(i,.6*bio,color=BLUE,height=.6); ax.barh(i,.4*bat,left=.6*bio,color=ORANGE,height=.6)
            ax.text(total+.015,i,f'{total:.3f}',va='center',fontsize=7)
        else: ax.text(.02,i,'Graph-only method' if k in ['bbknn_source','harmony_bbknn'] else 'Not reported',va='center',color=GREY,fontsize=7)
    ax.set_yticks(y,['']*len(y)); ax.set_ylim(9.5,-.5); ax.set_xlim(0,.82)
    heading(ax,'c  scIB embedding score'); ax.set_xlabel('0.6 × biology + 0.4 × batch')
    from matplotlib.patches import Patch
    fig.legend(handles=[Patch(color=BLUE,label='Biology contribution'),Patch(color=ORANGE,label='Batch contribution')],frameon=False,fontsize=7,loc='upper right',bbox_to_anchor=(.98,.93),ncol=2)
    fig.suptitle('Batch removal must preserve biological structure',fontsize=13,y=.98)
    fig.text(.22,.025,'Graph metrics: shared 300k-cell benchmark. scIB: up to 100k confidently labeled cells.\nMarker labels are derived from expression; they are not independent reference truth.',fontsize=8)
    save(fig,'fig3_batch_metrics',[gp,sp],'TME v1 shared benchmark; graph and embedding metrics have different evaluation sets')


def annotation():
    p=AUD/'annotation_pilots.json'; data=record(p)['20261103_marker_routed_annotation/GSE138709_v2/summary.json']
    route=data['route_counts']; status=data['reference_status_counts']
    assert sum(route.values())==sum(status.values())==11381 and not data['accepted_fine_identity']
    rnames={'excluded_malignant':'Malignant: excluded','immune_celltypist':'Immune CellTypist','vascular_celltypist':'Vascular CellTypist','singler_hpca':'SingleR HPCA','unresolved_no_compatible_reference':'No compatible reference'}
    snames={'candidate_immune_parent_supported':'Immune parent supported','candidate_hpca_parent_supported':'HPCA parent supported','candidate_vascular_parent_supported':'Vascular parent supported','reference_parent_conflict_or_pruned':'Conflict / pruned','reference_parent_conflict_or_unresolved':'Conflict / unresolved','unresolved_no_compatible_reference':'No compatible reference'}
    fig,axes=plt.subplots(1,2,figsize=(11,4.6));fig.subplots_adjust(wspace=.9,left=.18,right=.98,bottom=.18)
    for ax,d,names,title,exclude in [(axes[0],route,rnames,'a  Routing of 11,381 source cells',False),(axes[1],status,snames,'b  Reference status: 4,206 TME candidates',True)]:
        rows=[(names[k],v,k) for k,v in d.items() if not(exclude and k=='excluded_malignant')]
        y=np.arange(len(rows));vals=[x[1] for x in rows]
        ax.barh(y,vals,color=[ORANGE if 'conflict' in k or 'unresolved' in k else BLUE for _,_,k in rows]);ax.set_yticks(y,[x[0] for x in rows]);ax.invert_yaxis()
        for i,val in enumerate(vals):ax.text(val+max(vals)*.025,i,f'{val:,}',va='center')
        ax.set_xlim(0,max(vals)*1.28);ax.set_xlabel('Cells');heading(ax,title)
    fig.text(.18,.04,'GSE138709 routed candidate annotation v2. Parent support is not accepted fine identity.\nThe two immune models share a reference family and are not independent validation.',fontsize=8)
    save(fig,'fig4_annotation_routes',[p],'GSE138709 full-source routed annotation v2; not all 72 sources')


def nmf():
    rp=AUD/'tables/v7_cm/primary_strict/rank_survey_v1.tsv';npth=AUD/'nmf_scale_sensitivity.tsv'
    r,n=table(rp),table(npth)
    fig,axes=plt.subplots(2,2,figsize=(12,8));fig.subplots_adjust(hspace=.55,wspace=.65,left=.09,right=.97,bottom=.12)
    for ax,col,title,label in [(axes[0,0],'cophenetic','a  V7 strict rank survey','Cophenetic correlation'),(axes[1,0],'RSS','c  Reconstruction trade-off','Residual sum of squares')]:
        ax.plot(r['rank'],r[col],'-o',color=BLUE,ms=3,lw=1)
        ax.axvline(4,color=GREEN,ls='--',lw=1);ax.axvline(9,color=ORANGE,ls=':',lw=1)
        ax.set_xticks([2,4,6,9,12,16,20]);ax.set_xlabel('Rank K');ax.set_ylabel(label);heading(ax,title)
    axes[0,0].text(.04,.04,'K4: rule-based sensitivity\nK9: fixed historical comparison',transform=axes[0,0].transAxes,fontsize=8)
    names=n.method.str.replace('NMF_nsNMF','nsNMF');y=np.arange(len(n))
    for ax,old,new,title in [(axes[0,1],'equivalent_rescaling_original_cos','equivalent_rescaling_canonical_cos','b  Same fit after equivalent factor rescaling'),(axes[1,1],'original_seed_cos','canonical_seed_cos','d  Between-seed member similarity')]:
        for i,row in n.iterrows():ax.plot([row[old],row[new]],[i,i],color=GREY,lw=1)
        ax.scatter(n[old],y,color=ORANGE,s=20,label='Original matching',zorder=3)
        ax.scatter(n[new],y,color=GREEN,s=20,label='Column-scale canonicalized',zorder=3)
        ax.set_yticks(y,names);ax.invert_yaxis();ax.set_xlim(.60,1.025);ax.set_xlabel('Matched mean cosine');heading(ax,title)
    axes[0,1].legend(loc='lower center',bbox_to_anchor=(.5,1.12),fontsize=7,frameon=False,ncol=2)
    fig.text(.09,.025,'V7 benchmark: K9, five seeds per method. Methods differ in inputs and preprocessing.\nHigh member reproducibility is not evidence of a correct rank, stable patient activity or biological validity.',fontsize=8)
    save(fig,'fig5_cm_identifiability',[rp,npth],'Historical V7 strict rank and nine-method scale sensitivity; no refitting')


def spatial():
    files=[];fig,axes=plt.subplots(1,2,figsize=(11,4.8));fig.subplots_adjust(wspace=.3,left=.07,right=.98,bottom=.19,top=.83)
    for ax,seed in zip(axes,[20261103,20261104]):
        d=[]
        for membership in ['wtop10','network']:
            p=AUD/f'tables/v7_spatial/cm_blisa_primary_strict_K9only_seed{seed}_{membership}_bin50.tsv';files.append(p);d.append(table(p).set_index('cm'))
        cms=list(d[0].index);y=np.arange(len(cms))
        for offset,t,color,label in [(-.17,d[0],BLUE,'W top10'),(.17,d[1],ORANGE,'Author network')]:
            vals=t.reindex(cms);fractions=vals.units_z_pos/vals.n_units
            ax.barh(y+offset,fractions,height=.3,color=color,label=label)
            for i,(_,row) in enumerate(vals.iterrows()):
                ax.text(fractions.iloc[i]+.015,i+offset,f'{int(row.units_z_pos)}/{int(row.n_units)}',va='center',fontsize=6.5)
        ax.axvline(.5,color=GREY,ls='--',lw=.8);ax.set_xlim(0,1.22);ax.set_xticks([0,.5,1]);ax.set_yticks(y,cms);ax.invert_yaxis();ax.set_xlabel('Fraction of analysis units with z > 0')
        heading(ax,f'{"a" if seed==20261103 else "b"}  V7 strict K9, seed {seed}')
    fig.legend(*axes[1].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.53,.99),ncol=2,fontsize=7,frameon=False)
    fig.text(.07,.04,'50 micrometre bins; all nine modules are shown. Analysis units are not all verified independent patients.\nThe 0.5 guide is a sign fraction, not a significance threshold. Member definitions can change the denominator.',fontsize=8)
    save(fig,'fig6_spatial_membership',files,'V7 50um exploratory positivity; independent of CT75 numbering')



def v7_followup():
    folder = ROOT/'docs/reports/20261103_V7_CM_blisa/tables'
    lp, ap = folder/'v7_primary_K9_LOSO_r.tsv', folder/'CM_activity_share_by_sample_class.tsv'
    loso, activity = table(lp, index_col=0), table(ap, index_col=0)
    assert loso.shape == (53, 9) and activity.n_samples.sum() == 1207
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.8), gridspec_kw={'width_ratios':[1, 1.35]})
    fig.subplots_adjust(wspace=.62, left=.06, right=.96, bottom=.2, top=.85)
    rng = np.random.default_rng(20261105)
    ax = axes[0]
    for i, cm in enumerate(loso.columns):
        vals = loso[cm].to_numpy()
        ax.scatter(vals, i+rng.uniform(-.18,.18,len(vals)), s=9, alpha=.6, color=ORANGE if cm=='CM04' else BLUE, linewidths=0)
        ax.text(1.025, i, f'{int((vals>=.8).sum())}/53', va='center', fontsize=7)
    ax.axvline(.8,color=GREY,ls='--',lw=.8)
    ax.set_yticks(range(9),loso.columns);ax.invert_yaxis();ax.set_xlim(0,1.16)
    ax.set_xlabel('Matched W-column Pearson r on shared nodes')
    heading(ax,'a  All 53 leave-one-source-out results')
    ax=axes[1]
    values=activity[loso.columns].to_numpy()
    im=ax.imshow(values,aspect='auto',cmap='YlGnBu',vmin=0,vmax=1)
    labels=['Blood','Bone marrow','Lymph node','Metastatic tumour','Normal / adjacent','Primary / unspecified tumour','Unknown','Unmatched']
    ax.set_yticks(range(8),[f'{label} (n={n})' for label,n in zip(labels,activity.n_samples)])
    ax.set_xticks(range(9),loso.columns,rotation=45,ha='right')
    for i in range(8):
        for j in range(9):
            ax.text(j,i,f'{values[i,j]:.2f}',ha='center',va='center',fontsize=6.5,color='white' if values[i,j]>.55 else 'black')
    heading(ax,'b  Reported CM activity share by material class')
    fig.colorbar(im,ax=ax,fraction=.035,pad=.03,label='Reported activity share')
    fig.text(.06,.035,'V7 strict K9, reference seed 20261103. Each point omits one source; fits are not independent replicates.\nShared-node matching can hide coverage loss. Material classes are historical metadata, not verified patient-level labels.',fontsize=8)
    save(fig,'figS2_v7_source_material',[lp,ap],'Latest fa2d26a aggregate report; shared-node W correlations and material-level activity shares, not causal attribution')

if __name__=='__main__':
    datasets();batch();annotation();nmf();spatial();v7_followup()
    def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
    manifest={'renderer':'scripts/render_current_work_figures.py','renderer_sha256':sha(__file__),
              'python':platform.python_version(),'packages':{'matplotlib':matplotlib.__version__,'numpy':np.__version__,'pandas':pd.__version__},
              'figures':FIGURES,'inputs_sha256':{str(p.relative_to(ROOT)):sha(p) for p in sorted(SOURCES)},
              'outputs_sha256':{p.name:sha(p) for p in sorted(OUT.glob('*')) if p.suffix in ['.png','.pdf']},
              'new_statistical_tests':False,'new_scientific_model_fits':False}
    (OUT/'figure_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'figures':list(FIGURES),'output':str(OUT)}))
