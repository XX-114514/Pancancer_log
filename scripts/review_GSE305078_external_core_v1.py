"""Private complete external evidence and bootstrap-geometry review; no new test."""
import hashlib,json,os
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);src=r/'GSE305078_external_ligand_programme_core_results_20261016_v1'
    proof=src/'private_independent_audit.json';audit=json.loads(proof.read_text());assert audit['status'].startswith('PASS_ALL676')
    assert sha(src/'private_terminal.json')==audit['producer_terminal_sha256']
    terminal=json.loads((src/'private_terminal.json').read_text())
    for n,h in terminal['files'].items():assert sha(src/n)==h
    t=pd.read_csv(src/'private_all676_external_core_results.tsv',sep='\t');cfg=json.loads((r/'external_ligand_programme_core_followup_config_20261016_v1.json').read_text())
    out=r/'GSE305078_external_core_private_review_20261016_v1';assert not out.exists();out.mkdir(mode=0o700)
    rows=[]
    for i,row in enumerate(t.itertuples()):
        if not np.isfinite(row.wild_P):continue
        a=np.load(src/'private_resamples'/(row.model_id+'.npz'),allow_pickle=False);b=a['bootstrap_rho'];valid=np.isfinite(b)
        assert valid.sum()==row.bootstrap_valid and np.allclose(np.quantile(b[valid],[.025,.975]),[row.bootstrap_low,row.bootstrap_high])
        rng=np.random.default_rng(cfg['seed']+i*10+1);unique=np.array([len(np.unique(rng.integers(0,row.library_units,size=row.library_units))) for _ in range(2000)])
        rows.append(dict(model_id=row.model_id,library_units=row.library_units,main_design_rank=row.effective_rank,
                         bootstrap_valid=int(valid.sum()),unique_units_minimum=int(unique.min()),unique_units_median=float(np.median(unique)),
                         bootstrap_fraction_unique_le_main_rank_plus2=float(np.mean(unique<=row.effective_rank+2)),
                         finite_bootstrap_fraction_at_abs_rho1=float(np.mean(np.abs(b[valid])>=1-1e-10)),
                         original_CI_width=row.bootstrap_high-row.bootstrap_low,no_draws_removed=True))
    geometry=pd.DataFrame(rows);assert len(geometry)==terminal['estimated_models']
    geometry.to_csv(out/'private_all_estimated_bootstrap_geometry.tsv',sep='\t',index=False)
    coverage=t.groupby(['stage','extension_family','readout','status'],dropna=False).size().reset_index(name='planned_slots')
    coverage.to_csv(out/'private_all676_coverage.tsv',sep='\t',index=False)
    estimates=t[t.wild_P.notna()].sort_values(['ligand','readout','model_id']).copy();estimates.to_csv(out/'private_all_estimated_forest_data.tsv',sep='\t',index=False,na_rep='NA')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9})
    fig,ax=plt.subplots(figsize=(10,9));colors=plt.get_cmap('tab10');ligands=sorted(estimates.ligand.unique())
    for i,row in enumerate(estimates.itertuples()):
        ax.plot([row.bootstrap_low,row.bootstrap_high],[i,i],color=colors(ligands.index(row.ligand)),lw=1,alpha=.75)
        ax.scatter(row.partial_rank_rho,i,color=colors(ligands.index(row.ligand)),marker='s' if row.readout=='masked_programme' else '^',s=30)
        ax.text(1.045,i,f'{row.library_units} / {row.BH_q_complete676:.2f}',va='center',fontsize=8)
    ax.set_yticks(np.arange(len(estimates)),[f'{x.model_id} | {x.ligand} | '+('programme' if x.readout=='masked_programme' else str(x.regulator)) for x in estimates.itertuples()])
    ax.invert_yaxis();ax.axvline(0,color='#888888',ls='--',lw=.8)
    ax.set(xlim=(-1.02,1.02),xlabel='Full core-control partial rank correlation; original 95% case-row bootstrap interval',title='All estimable TNBC lineage-extension models; complete 676-slot family')
    ax.text(1.045,-1,'N / BH q',fontsize=8,va='center');fig.text(.29,.025,'Pretreatment physical libraries; patient independence and exact NSCLC author-state mapping unverified.',fontsize=8)
    fig.subplots_adjust(left=.29,right=.88,bottom=.09,top=.94);fig.savefig(out/'private_all_estimable_forest.png',dpi=150);fig.savefig(out/'private_all_estimable_forest.svg');plt.close(fig)
    report=['# Complete external lineage-extension review','',f'Input independent audit: {sha(proof)}','',
            'This is a frozen source-stage lineage extension, not exact author-state replication. No original family P/q changed and no later node/TF/phase substitution.','',
            '## Complete planned coverage','',coverage.to_string(index=False),'','## All estimable effects, including opposing and imprecise cases','',
            estimates[['model_id','stage','ligand','sender_node','receiver_node','readout','library_units','partial_rank_rho','wild_P','bootstrap_low','bootstrap_high','BH_q_complete676','BY_q_complete676','MC_upper_BH_q_complete676','LOO_direction_stable']].to_string(index=False),'',
            '## Bootstrap geometry','',geometry.to_string(index=False),'',
            'Rows resampled with replacement can leave very few distinct libraries relative to the complete control rank. Near-unit correlations and wide intervals are retained; no resamples were dropped to narrow the interval. Distinct-library count is a diagnostic, not a replacement denominator or rerun rule.','',
            '## Same-core discovery bridge (descriptive, selected parents)','',pd.read_csv(src/'private_all_discovery_core10_bridge_points.tsv',sep='\t').to_string(index=False),'',
            '## Decision','',
            'No robust external support from this fixed family; retain all results and stop same-cohort threshold/TF/endpoint iterations. Opposing programme signs and mixed functional directions constrain transport to these captured fine nodes; different cancer, pre/post-treatment context and unverified exact states prevent a blanket rejection of original NSCLC biology.','',
            'Continue exact context/MCP/cell-end clinical and multicellular-coordination questions without borrowing another programme or treating post-treatment pathology association as prediction. Receptor binding, surface proteins, ligand release and actual signalling remain unverified.']
    (out/'PRIVATE_COMPLETE_EXTERNAL_REVIEW.md').write_text('\n'.join(report)+'\n')
    result=dict(status='COMPLETE_PRIVATE_ALL676_COVERAGE_ALL_ESTIMABLE_FOREST_AND_BOOTSTRAP_DIAGNOSTICS_NO_NEW_INFERENCE_DISPLAY_AUDIT_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),
                code_sha256=sha(__file__),input_audit_sha256=sha(proof),files={p.name:sha(p) for p in out.iterdir() if p.is_file()},scientific_table_or_learned_parameter_public_exported=False)
    (out/'private_terminal.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
