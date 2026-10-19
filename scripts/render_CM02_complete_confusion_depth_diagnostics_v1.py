"""Required complete confusion and depth diagnostics; no refit or new P."""
import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import rankdata

def main(a):
    root=a.root;r=root/'model_validation_v1';out=root/'confusion_depth_diagnostics_v1';assert not out.exists();out.mkdir(mode=0o700)
    predictions=pd.read_csv(r/'private_all_heldout_predictions.tsv.gz',sep='\t',float_precision='round_trip')
    meta=pd.read_csv(root/'reference_native_inputs_v2/private_selected_reference_metadata.tsv.gz',sep='\t',keep_default_na=False)
    raw=np.load(root/'reference_native_inputs_v2/private_panel_counts.npy',mmap_mode='r');depth=raw.sum(1)
    groups=meta.patient_block.astype(str).to_numpy();sources=sorted(meta.dataset.unique());modes=['native','matched_depth_stress']
    classes=['B_Plasma','Myeloid','NK','Nonimmune','T'];plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'svg.fonttype':'none'})
    def save(fig,name):
        fig.savefig(out/(name+'.png'),dpi=400,bbox_inches='tight');fig.savefig(out/(name+'.pdf'),bbox_inches='tight');fig.savefig(out/(name+'.svg'),bbox_inches='tight');fig.savefig(out/(name+'_preview.jpg'),dpi=120,bbox_inches='tight');plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(12,8))
    for i,mode in enumerate(modes):
        for j,source in enumerate(sources):
            matrix=np.load(r/(source+'_'+mode+'_private_confusion.npy'));matrix=matrix/matrix.sum(1,keepdims=True)
            ax=axes[i,j];ax.imshow(matrix,vmin=0,vmax=1,cmap='Blues');ax.set(xticks=range(5),yticks=range(5),xticklabels=classes,yticklabels=classes,xlabel='Predicted',ylabel='Reference proposal',title=source+' / '+('native' if i==0 else 'depth stress'));ax.tick_params(axis='x',labelrotation=45)
            for y in range(5):
                for x in range(5):ax.text(x,y,f'{matrix[y,x]:.2f}',ha='center',va='center',fontsize=7,color='white' if matrix[y,x]>.6 else '#222')
    fig.suptitle('All six held-study confusion matrices: case-balanced, row-normalized');fig.tight_layout();save(fig,'All6_confusion_matrices')
    records=[];fig,axes=plt.subplots(2,3,figsize=(12,7),sharey=True)
    for k,target in enumerate(['regulatory','cytotoxic']):
        for j,source in enumerate(sources):
            ax=axes[k,j]
            for mode,color in zip(modes,['#1679a5','#df7b24']):
                d=predictions[(predictions.kind=='programme')&(predictions.target==target)&(predictions.holdout==source)&(predictions['mode']==mode)]
                rows=d.reference_row.to_numpy(int);truth=d.truth.to_numpy(float);pred=d.predicted.to_numpy(float);error=pred-truth
                patient=groups[rows];counts=pd.Series(patient).value_counts();w=1/pd.Series(patient).map(counts).to_numpy(float)
                values=[]
                for low,high,label in [(0,20,'1-20'),(20,50,'21-50'),(50,100,'51-100'),(100,np.inf,'>100')]:
                    take=(depth[rows]>low)&(depth[rows]<=high);n=int(take.sum());eligible=n>=100 and len(set(patient[take]))>=5
                    rmse=float(np.sqrt(np.average(error[take]**2,weights=w[take]))) if eligible else np.nan
                    records.append(dict(holdout=source,target=target,mode=mode,native_shared_panel_depth=label,cells=n,cases=len(set(patient[take])),diagnostic_eligible=eligible,case_balanced_RMSE=rmse));values.append(rmse)
                ax.plot(range(4),values,'o-',label='native' if mode=='native' else 'depth stress',color=color)
            ax.set(xticks=range(4),xticklabels=['1-20','21-50','51-100','>100'],xlabel='Native shared-panel RNA counts',ylabel='Case-balanced RMSE',title=source+' / '+target)
            if k==0 and j==0:ax.legend(frameon=False)
    fig.suptitle('All studies and programmes: error versus native RNA depth\nFixed descriptive bins; missing insufficient strata retained, no model selection or P');fig.tight_layout();save(fig,'All_programmes_native_depth_error')
    pd.DataFrame(records).to_csv(out/'private_all_native_depth_error_strata.tsv',sep='\t',index=False)
    (out/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Complete depth and confusion diagnostics</title><main style="max-width:1200px;margin:auto;font:16px system-ui"><h1>All held-study confusion and depth diagnostics</h1>'+''.join('<section><img style="width:100%" src="'+n+'.png"><p><a href="'+n+'.pdf">PDF</a> | <a href="'+n+'.svg">SVG</a></p></section>' for n in ['All6_confusion_matrices','All_programmes_native_depth_error'])+'</main>')
    proof=dict(status='ALL_CONFUSION_AND_FIXED_DEPTH_DIAGNOSTICS_RENDERED_NO_REFIT',finished_utc=datetime.now(timezone.utc).isoformat(),producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),no_new_P_or_changed_admission=True)
    (out/'render_technical_public.json').write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);main(p.parse_args())
