"""All-source descriptive marker review after full numeric audit; no classifier."""
import json
import os
from datetime import datetime,timezone
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from NSCLC_native_MFA_clinical_inputs_v1 import matched_quota,sha,VIEWS

DISPLAY={
    'author_CD8':['CD3D','CD8A','CD8B'],
    'author_Treg':['FOXP3','IL2RA','CTLA4'],
    'author_Macrophage':['LST1','TYROBP','CSF1R'],
    'author_B_Plasma':['CD79A','MS4A1','MZB1'],
}


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);matched_quota(root)
    src=root/'NSCLC_native_MFA_external_role_marker_profiles_20261016_v2'
    proof=src/'private_independent_all_role_marker_audit.json'
    assert sha(proof)=='b40688206ea0abfa9ffd0bb795e85e871d59c745db13646416a3831f32d3c1b4'
    audit=json.loads(proof.read_text());assert audit['producer_terminal_sha256']==sha(src/'private_terminal.json')
    terminal=json.loads((src/'private_terminal.json').read_text())
    for name,item in terminal['files'].items():assert sha(src/name)==item['sha256']
    out=root/'NSCLC_native_MFA_external_role_marker_review_20261016_v1';assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),input_audit_sha256=sha(proof),
        scope='All61 source-arm x4 frozen candidate pools x44 prespecified marker genes. ALL units and missingness retained in upstream; aggregate only n>=20 library pools for this display.',
        summary='Cell-weighted marker means/fractions over eligible libraries; per-library median also. Source aggregates are not independent patients, inferential effects or subtype validation.',
        display=DISPLAY,selection='Same fixed three lineage markers per view shown for EVERY source; entire44-marker summary saved. No significance or outcome selection.',
        interpretation='Sources lacking author fine labels are not assumed devoid of these cells. Marker expression is descriptive, not a pass/fail classifier or new cell-state assignment.',
        no_new_P_or_identity_changes=True)
    (out/'plan_before_complete_marker_profile_review.json').write_text(json.dumps(plan,indent=2)+'\n')
    units=pd.read_csv(src/'private_all6760_frozen_role_units.tsv',sep='\t',keep_default_na=False)
    with np.load(src/'private_fixed_marker_pool_sums.npz',allow_pickle=False) as z:
        sizes=z['n_cells'];logs=z['sum_logCP10K'];positive=z['positive_counts'];available=z['source_measured'];markers=z['markers'].astype(str)
    assert len(markers)==44 and len(units)==6760
    sources=units[['arm','dataset']].drop_duplicates().sort_values(['arm','dataset'])
    records=[]
    for (arm,source,view),g in units.groupby(['arm','dataset','view'],sort=True):
        idx=g.loc[g.n20,'unit_index'].to_numpy(int);n=int(sizes[idx].sum());all_n=int(g.n_cells.sum())
        for j,gene in enumerate(markers):
            measurement=available[g.unit_index.to_numpy(int)//4,j];assert np.all(measurement==measurement[0])
            measured=bool(measurement[0]);status='MEASURED' if measured and len(idx) else 'NO_LIBRARY_POOL_GE20' if measured else 'UNMEASURED_GENE_OR_UNRESOLVED'
            records.append(dict(arm=arm,dataset=source,view=view,gene=gene,all_candidate_cells=all_n,
                eligible_libraries=len(idx),eligible_candidate_cells=n,status=status,
                pooled_mean_logCP10K=float(logs[idx,j].sum()/n) if status=='MEASURED' else np.nan,
                pooled_positive_fraction=float(positive[idx,j].sum()/n) if status=='MEASURED' else np.nan,
                median_library_mean_logCP10K=float(np.median(logs[idx,j]/sizes[idx])) if status=='MEASURED' else np.nan,
                median_library_positive_fraction=float(np.median(positive[idx,j]/sizes[idx])) if status=='MEASURED' else np.nan))
    result=pd.DataFrame(records);assert len(result)==61*4*44
    result.to_csv(out/'private_all10736_source_role_marker_descriptions.tsv.gz',sep='\t',index=False,compression='gzip')
    indexed=result.set_index(['arm','dataset','view','gene']);assert indexed.index.is_unique
    fig,axes=plt.subplots(1,4,figsize=(13,19),sharey=True)
    cmap=plt.get_cmap('viridis');norm=matplotlib.colors.Normalize(0,3)
    for vi,view in enumerate(VIEWS):
        ax=axes[vi]
        for i,r in enumerate(sources.itertuples()):
            for j,gene in enumerate(DISPLAY[view]):
                row=indexed.loc[(r.arm,r.dataset,view,gene)]
                if row.status=='MEASURED':
                    ax.scatter(j,i,s=8+110*row.pooled_positive_fraction,c=[cmap(norm(row.pooled_mean_logCP10K))],edgecolors='#666666',linewidths=.3)
                elif row.status=='NO_LIBRARY_POOL_GE20':ax.scatter(j,i,marker='x',s=18,c='#999999',linewidths=.6)
                else:ax.scatter(j,i,marker='s',s=14,c='#dddddd',linewidths=0)
        ax.set_xticks(range(3),DISPLAY[view],rotation=60,ha='right',fontsize=9)
        ax.set_xlim(-.5,2.5);ax.set_ylim(60.8,-.8);ax.set_title('Frozen '+view.removeprefix('author_')+' candidate',fontsize=9)
        ax.grid(axis='y',alpha=.15)
        ax.set_yticks(range(61),[('I: ' if r.arm=='archived_immune_atlas' else 'O: ')+r.dataset for r in sources.itertuples()],fontsize=7)
    fig.suptitle('All61 sources: fixed marker footprints of four frozen candidate pools\nEligible library pools >=20 cells; dot size=positive-cell fraction; color=pooled mean logCP10K\nX=no eligible library pool; gray square=gene unmeasured/unresolved. Descriptive only; no identity acceptance.',fontsize=10,y=.995)
    fig.tight_layout(rect=[0,0,.94,.955]);cax=fig.add_axes([.945,.25,.012,.5])
    fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,label='mean logCP10K (colors capped at3)')
    fig.savefig(out/'private_all61_four_role_marker_footprints.png',dpi=160);plt.close(fig)
    report=dict(status='COMPLETE_ALL10736_SOURCE_ROLE_MARKER_DESCRIPTIONS_AND_ALL61_SOURCE_FIGURE',
        finished_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),
        files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()},
        no_inference_or_identity_promotion=True)
    (out/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],finished_UTC=report['finished_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
