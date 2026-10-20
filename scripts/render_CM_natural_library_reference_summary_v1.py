"""Complete source-level recovery heatmaps, including unavailable positions."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    root = Path(os.environ['PANCANCER_PARENT_OUTPUT'])
    out = root/'private_source_figures_v1'; assert not out.exists()
    terminal = json.loads((root/'private_terminal.json').read_text())
    p = root/'private_all_target_recovery_metrics.json'
    assert sha(p) == terminal['files'][p.name]['SHA256']
    rows = json.loads(p.read_text())
    plan = json.loads((root/'private_plan_before_native_expression.json').read_text())['plan']
    sources = plan['prior_reference_plan']['source_scope']
    domains = plan['source_domains']; architectures = plan['reference_architectures']
    scales = plan['scalings']; views = plan['coefficient_views']
    core = plan['supported_reference_classes']
    targets = core+['Other_supported_reference_pool','Boundary_supported_reference_pool']+list(plan['marginal_targets'])+['paired_R_only_E_minimum']
    assert len(targets) == 28 and len(set(targets)) == 28
    out.mkdir(mode=0o700)
    for scale in scales:
        for view in views:
            sub = [r for r in rows if r['scaling'] == scale and r['coefficient_view'] == view]
            x = np.full((len(targets),16),np.nan); exists = np.zeros_like(x,dtype=bool)
            labels = []
            for si,source in enumerate(sources):
                for di,domain in enumerate(domains):
                    for ai,arch in enumerate(architectures):
                        col = si*4+di*2+ai
                        labels.append(source+'\n'+('cross' if di==0 else 'within')+' / '+('separate' if ai==0 else 'merged'))
                        for r in sub:
                            if r['source']==source and r['domain']==domain and r['architecture']==arch:
                                j=targets.index(r['target']); exists[j,col]=True
                                if r['MAE'] is not None: x[j,col]=r['MAE']
            assert int(exists.sum()) == len(sub) == 344
            values = np.log10(1+100*x); maximum = max(float(np.nanmax(values)),.1)
            cmap=plt.get_cmap('viridis').copy();cmap.set_bad('#ededed')
            fig,ax=plt.subplots(figsize=(14,13))
            fig.subplots_adjust(left=.34,right=.92,bottom=.22,top=.88)
            im=ax.imshow(values,aspect='auto',interpolation='nearest',cmap=cmap,norm=Normalize(0,maximum))
            for j,i in zip(*np.where(~exists)):
                ax.text(i,j,'–',ha='center',va='center',color='#888888',fontsize=9)
            for j,i in zip(*np.where(exists & np.isnan(x))):
                ax.text(i,j,'NA',ha='center',va='center',fontsize=7)
            ax.set_yticks(np.arange(28),targets,fontsize=9)
            ax.set_xticks(np.arange(16),labels,rotation=65,ha='right',fontsize=8)
            for edge in [3.5,7.5,11.5]: ax.axvline(edge,color='white',lw=2)
            ax.set_title('All fixed targets: source-level mean absolute recovery error\n'+scale+' | '+view,fontsize=12,pad=12)
            cb=fig.colorbar(im,ax=ax,fraction=.025,pad=.02)
            ticks=[v for v in [0,.001,.01,.05,.1,.25,.5,1,2,5,10,100] if np.log10(1+100*v)<=maximum]
            cb.set_ticks(np.log10(1+100*np.asarray(ticks)),[str(v) for v in ticks]);cb.set_label('MAE (nonlinear color scale)',fontsize=9)
            fig.text(.02,.045,'Cross = held-out source; within = held-out conservative claimed-identity component. Libraries are not newly verified patients.\nSeparate and merged references use identical eligible training cells. All 24 captured RNA classes remain in truth denominators.\n– = target is not defined in this architecture; NA = registered target is unestimable. Raw coefficients are not automatically RNA fractions.\nDiagnostic errors, without new P/FDR or clinical/identity acceptance. Complete tables include correlations, calibration and missingness.',fontsize=9)
            name='private_recovery_'+scale+'_'+view
            for ext in ['png','svg']: fig.savefig(out/(name+'.'+ext),dpi=140)
            plt.close(fig)
    files={p.name:dict(SHA256=sha(p),bytes=p.stat().st_size) for p in out.iterdir()}
    result=dict(status='ALL_FOUR_RECOVERY_HEATMAPS_RENDERED_VISUAL_REVIEW_PENDING',completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(__file__),input_SHA256=sha(root/'private_all_target_recovery_metrics.json'),registered_metric_positions=1376,all_positions_including_missing_represented=True,files=files)
    (out/'private_figure_terminal.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'])))


if __name__ == '__main__': main()
