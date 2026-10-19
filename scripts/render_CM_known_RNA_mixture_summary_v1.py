"""All seventeen targets and two prespecified gene scales, source aggregates only."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    root=Path(os.environ['PANCANCER_RAM_OUTPUT']);out=root/'private_source_figures_v1';assert not out.exists();out.mkdir(mode=0o700)
    t=json.loads((root/'private_terminal.json').read_text())
    p=root/'private_all_source_target_recovery_metrics.json';assert hashlib.sha256(p.read_bytes()).hexdigest()==t['files'][p.name]['SHA256']
    table=pd.DataFrame(json.loads(p.read_text()));plan=json.loads((root/'private_plan_before_expression_and_mixture_effects.json').read_text())['plan']
    targets=plan['reference_classes']+list(plan['marginal_targets'])+['paired_R_only_E_minimum'];sources=plan['prior_reference_plan']['source_scope']
    labels=['T anchored reg only','T anchored reg + effector','T effector only','T reg without anchor','T reg without anchor + effector','T neither capture','T core not captured','NK','Other assigned','Boundary / unresolved','All T','Anchored regulatory total','Broad regulatory total','Effector total','Broad reg + effector overlap','Effector without anchored reg','Paired R-only / E minimum']
    scales=plan['scalings'];fig,axes=plt.subplots(2,2,figsize=(15.5,13.5),sharey=True)
    fig.subplots_adjust(left=.26,right=.86,top=.87,bottom=.14,wspace=.20,hspace=.42)
    matrices=[]
    for regime,value in [('combined_mixed','MAE'),('pure_profile_anchors','designed_zero_max_contribution')]:
        for scaling in scales:
            sub=table[table.regime.eq(regime)&table.scaling.eq(scaling)]
            matrix=sub.pivot(index='target',columns='source',values=value).reindex(index=targets,columns=sources).to_numpy(float)
            matrices.append(matrix)
    for row in range(2):
        vmax=max(float(np.nanmax(matrices[row*2])),float(np.nanmax(matrices[row*2+1])))
        for col,scaling in enumerate(scales):
            ax=axes[row,col];values=matrices[row*2+col]
            im=ax.imshow(np.ma.masked_invalid(values),vmin=0,vmax=max(vmax,.1),cmap='YlOrRd',aspect='auto')
            ax.set_xticks(range(len(sources)),sources,rotation=25,ha='right',fontsize=9)
            ax.set_yticks(range(len(targets)),labels,fontsize=8)
            if col:ax.tick_params(labelleft=False)
            ax.set_title(('Train-mean weighted (primary)' if scaling=='train_mean_weighted' else 'Unweighted (fixed sensitivity)')+'\n'+('Mixed RNA: mean absolute\nfraction error' if row==0 else 'Pure profile: max false RNA fraction\nfor designed-zero target'),fontsize=10)
            for i in range(values.shape[0]):
                for j in range(values.shape[1]):
                    if np.isfinite(values[i,j]):ax.text(j,i,f'{values[i,j]:.2f}',ha='center',va='center',fontsize=7,color='black' if values[i,j]<vmax*.65 else 'white')
            for edge in [6.5,9.5,15.5]:ax.axhline(edge,color='0.45',lw=.6)
        position=axes[row,1].get_position()
        color_axis=fig.add_axes([.89,position.y0,.015,position.height])
        bar=fig.colorbar(im,cax=color_axis)
        bar.set_label('RNA fraction error; lower is better',fontsize=9)
    fig.suptitle('2026-10-14 immune atlas CM (v1): known RNA-mixture source-heldout calibration\nAll 10 components / 17 targets retained; these are simulation diagnostics',fontsize=13)
    fig.text(.025,.025,'Mixed panel combines equal numbers of dense, sparse and T-diluted known-RNA mixtures.\nPure profiles are pooled operational-class RNA, not purified biological cell types; zero means not added to that simulation.\nRNA fraction errors are not original tissue cell proportions, patient outcomes, P/FDR or direct communication validation.',fontsize=9)
    for suffix in ['png','svg']:fig.savefig(out/('private_all17_target_mixture_errors_and_false_contributions.'+suffix),dpi=180)
    plt.close(fig)
    files={p.name:dict(bytes=p.stat().st_size,SHA256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in out.iterdir()}
    result=dict(status='SOURCE_MIXTURE_FIGURE_RENDERED_VISUAL_REVIEW_PENDING',completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),parent_terminal_SHA256=hashlib.sha256((root/'private_terminal.json').read_bytes()).hexdigest(),files=files)
    (out/'private_figure_terminal.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(status=result['status'],files=files)))


if __name__=='__main__':main()
