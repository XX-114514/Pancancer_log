"""Source-balanced technical figures; no learned gene lists or individual rows."""
import hashlib,json,os,time,shutil
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    parent=Path('/dev/shm/pancancer-state-proxy-quality-20261020-v2/private_output')
    audit=Path('/dev/shm/pancancer-state-proxy-quality-audit-20261020-v2/private_output')
    out=Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert not out.exists() and out.parent.parent==Path('/dev/shm') and shutil.disk_usage('/dev/shm').free>512*1024**2
    for _ in range(360):
        if (audit/'private_terminal.json').is_file():break
        time.sleep(5)
    else:raise TimeoutError('Wait for same independent auditor, no restart')
    a=json.loads((audit/'private_terminal.json').read_text());assert a['status'].startswith('ACCEPTED_ALL_QUALITY_FIDELITY')
    assert a['producer_terminal_sha256']==sha(parent/'private_terminal.json')
    t=json.loads((parent/'private_terminal.json').read_text())
    for name in ['private_all_source_major_programme_axis_fidelity.tsv.gz','private_complete_programme_major_axis_triage.tsv']:
        assert sha(parent/name)==t['files'][name]['sha256']
    s=pd.read_csv(parent/'private_all_source_major_programme_axis_fidelity.tsv.gz',sep='\t')
    q=pd.read_csv(parent/'private_complete_programme_major_axis_triage.tsv',sep='\t')
    assert len(s)==127636 and len(q)==3256
    axes=sorted(q.axis_id.unique());labels={a:'Panel '+chr(65+i) for i,a in enumerate(axes)}
    out.mkdir(mode=0o700);coverage=[];raw=[]
    for axis in axes:
        g=q[q.axis_id.eq(axis)]
        assert len(g)==814
        coverage.append(dict(panel=labels[axis],axis_sha256=axis,all_frozen_programmes=len(g),
                             any_qualifiable_source=int(g.qualifiable_sources.gt(0).sum()),
                             at_least3_qualifiable_sources=int(g.qualifiable_sources.ge(3).sum()),
                             stable_and_at_least3_sources=int((g.algorithm_stable_candidate&g.qualifiable_sources.ge(3)).sum()),
                             exploratory_triage_pass=int(g.exploratory_proxy_triage_pass.sum()),
                             original_full_signature_spatial_gate_pass=0))
    cov=pd.DataFrame(coverage);cov.to_csv(out/'source_balanced_panel_stage_counts.tsv',sep='\t',index=False)
    for row in q.itertuples(index=False):
        g=s[s.programme_id.eq(row.programme_id)&s.major.eq(row.major)&s.axis_id.eq(row.axis_id)&s.estimable_identity_strata.gt(0)&s.partial_signature_genes.ge(5)]
        raw.append(dict(panel=labels[row.axis_id],axis_sha256=row.axis_id,programme_id=row.programme_id,major=row.major,
                        source_count=len(g),source_median_raw_cell_rho=float(g.raw_block_weighted_rho.median()) if len(g) else np.nan,
                        source_median_adjusted_cell_rho=row.source_median_adjusted_cell_rho,
                        source_median_adjusted_block_rho=row.source_median_adjusted_block_rho,
                        provisional_technical_pass=bool(row.exploratory_proxy_triage_pass)))
    desc=pd.DataFrame(raw);desc.to_csv(out/'private_all3256_source_balanced_descriptors.tsv',sep='\t',index=False)
    majors=q.groupby(['major','axis_id'],sort=True).agg(programme_positions=('programme_id','size'),
            provisional_technical_pass=('exploratory_proxy_triage_pass','sum')).reset_index()
    majors['panel']=majors.axis_id.map(labels);majors.to_csv(out/'major_panel_count_summary.tsv',sep='\t',index=False)
    # Broad source summaries have no programme identifiers, signatures or unit rows.
    sr=s.groupby(['gse_id','major','axis_id'],sort=True).agg(programme_positions=('programme_id','size'),
                  estimable_programme_positions=('estimable_identity_strata',lambda z:int((z>0).sum()))).reset_index()
    sr.to_csv(out/'all_source_major_panel_coverage_counts.tsv',sep='\t',index=False)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none'})
    fig,axs=plt.subplots(2,2,figsize=(10.5,9.5))
    for ax,axis in zip(axs.ravel(),axes):
        g=desc[desc.axis_sha256.eq(axis)&desc.source_count.gt(0)]
        ax.scatter(g.source_median_raw_cell_rho,g.source_median_adjusted_cell_rho,s=12,c='#808b96',alpha=.45,linewidths=0)
        good=g[g.provisional_technical_pass]
        ax.scatter(good.source_median_raw_cell_rho,good.source_median_adjusted_cell_rho,s=20,c='#167c80',linewidths=0)
        ax.plot([-1,1],[-1,1],color='.7',lw=.8);ax.axhline(0,color='.85',lw=.7);ax.axvline(0,color='.85',lw=.7)
        ax.set(xlim=(-1.02,1.02),ylim=(-1.02,1.02),title=f'{labels[axis]} | {len(g)}/814 with reference coverage',
               xlabel='Median raw within-identity correlation',ylabel='Median quality-adjusted correlation')
    fig.suptitle('Partial-panel fidelity to reference genes outside the panel',fontsize=14)
    fig.tight_layout(rect=(0,.12,1,.94))
    fig.text(.035,.025,'Each point is one frozen programme/panel, summarized with equal source weight.\nRNA amount and detected-gene ranks are required covariates; each block receives equal total weight.\nTeal: provisional technical triage only. Cross-K programmes are correlated; no spatial or biological validation.',fontsize=9)
    for ext in ['png','svg']:fig.savefig(out/('source_balanced_quality_fidelity.'+ext),dpi=170)
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(10,5.4));x=np.arange(4);width=.18
    for j,(col,label,color) in enumerate([('any_qualifiable_source','At least one source','#b5c3cf'),('at_least3_qualifiable_sources','At least three sources','#778f9e'),('stable_and_at_least3_sources','Stable and at least three','#b7a17b'),('exploratory_triage_pass','Provisional fidelity triage','#167c80')]):
        ax.bar(x+(j-1.5)*width,cov[col],width,label=label,color=color)
    ax.set_xticks(x,cov.panel);ax.set(ylabel='Frozen programme positions (814 per panel)',title='Coverage and provisional technical triage')
    ax.legend(fontsize=8);fig.tight_layout(rect=(0,.16,1,1))
    fig.text(.04,.025,'Full-signature spatial measurement gates remain unmet for all four panels.\nAt least 50 cells /10 blocks per identity stratum; sparse strata are unavailable, not biological negatives.\nProject thresholds are exploratory measurement triage, not validated state labels or interaction tests.',fontsize=9)
    for ext in ['png','svg']:fig.savefig(out/('source_balanced_fidelity_coverage.'+ext),dpi=170)
    plt.close(fig)
    report=dict(status='COMPLETE_SOURCE_ONLY_QUALITY_FIDELITY_FIGURES_PENDING_ACTUAL_VISUAL_REVIEW',completed_UTC=datetime.now(timezone.utc).isoformat(),
                source_sha256=sha(__file__),producer_terminal_sha256=sha(parent/'private_terminal.json'),audit_terminal_sha256=sha(audit/'private_terminal.json'),
                programme_panel_descriptors=3256,figure_panels=4,source_count_rows=len(sr),
                individual_records_learned_gene_names_weights_or_matrices_output=False,scientific_outputs_public=False,
                full_signature_gate_or_identity_changed=False,
                files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()})
    assert sum(v['bytes'] for v in report['files'].values())<8*1024**2
    (out/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['status','programme_panel_descriptors','figure_panels','source_count_rows']}))


if __name__=='__main__':main()
