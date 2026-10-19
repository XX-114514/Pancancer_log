"""Only source/cancer/lineage aggregates leave the private modelling directory."""
import hashlib,json,os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    control=Path(os.environ['PANCANCER_MV_CONTROL']);out=Path(os.environ['PANCANCER_REPORT_OUTPUT']);assert not out.exists();out.mkdir(mode=0o700)
    audit=Path('/dev/shm/pancancer-real-spatial-mv-background-audit-20261020-v1/private_output/private_terminal.json')
    assert sha(audit)=='bb1ae0a3aef7e458905f71bcbb2cc54b050101ff500b907797a821bd0eae7ff1'
    x=pd.read_csv(control/'private_output/private_complete126_blocked_model_cases.tsv',sep='\t');assert len(x)==126
    contract=json.loads(Path(json.loads((control/'plan.json').read_text())['private_contract']).read_text())
    assert set(x.cancer)<={'LUAD','pleural_mesothelioma'} and set(x.response_lineage)<={'T cell','Macrophage'}
    assert set(x.readout)=={'RAW','NEGATIVE_EXCESS'} and x.dataset_id.nunique()==1
    effects=[]
    for key,b in x.groupby(['dataset_id','cancer','response_lineage','readout']):
        z=b[b.status.eq('ESTIMABLE_BLOCKED_PREDICTION')];row=dict(zip(['public_GSE','cancer','response_lineage','readout'],key))
        row.update(full_family_positions=len(b),estimable_positions=len(z),positive_gain_positions=int(z.gain_R2.gt(0).sum()),positive_full_R2_positions=int(z.multiview_R2.gt(0).sum()),positive_gain_and_full_R2_positions=int((z.gain_R2.gt(0)&z.multiview_R2.gt(0)).sum()))
        for col in ['gain_R2','baseline_R2','multiview_R2']:
            for label,q in [('q25',.25),('median',.5),('q75',.75)]:row[label+'_'+col]=float(z[col].quantile(q)) if len(z) else np.nan
        effects.append(row)
    e=pd.DataFrame(effects);assert len(e)==8;e.to_csv(out/'source_cancer_lineage_readout_model_summaries.tsv',sep='\t',index=False)
    paired=[];coverage=[]
    for i,t in enumerate(contract['targets']):
        f=pd.read_csv(control/('private_output/private_image_%02d_folds.tsv.gz'%i),sep='\t');s=pd.read_csv(control/('private_output/private_image_%02d_field_scope.tsv.gz'%i),sep='\t')
        for status,n in s.status.value_counts().items():coverage.append(dict(public_GSE=t['dataset_id'],cancer=t['cancer'],status=status,technical_fields=int(n)))
        for m in t['proxy_mappings']:
            z=x[x.image_ordinal.eq(i)&x.target_proxy_group_id.eq(m['proxy_group_id'])];assert len(z)==2
            okay=bool(z.status.eq('ESTIMABLE_BLOCKED_PREDICTION').all());ff=f[f.target_proxy_group_id.eq(m['proxy_group_id'])] if len(f) else f
            other=[set(v['signature_indices']) for v in t['proxy_mappings'] if v['major']==m['major'] and v['proxy_group_id']!=m['proxy_group_id']];own=set(m['signature_indices']);union=set().union(*other)
            paired.append(dict(public_GSE=t['dataset_id'],cancer=t['cancer'],response_lineage=m['major'],proxy=m['proxy_group_id'],estimable=okay,
                both_readouts_positive=bool(okay and z.gain_R2.gt(0).all()),all_six_folds_positive=bool(okay and len(ff)==6 and ff.status.eq('ESTIMABLE_BLOCKED_PREDICTION').all() and ff.gain_R2.gt(0).all()),
                overlap=max(len(own&a)/len(own|a) for a in other),covered_by_other_union=len(own&union)/len(own)))
    p=pd.DataFrame(paired);rows=[]
    for key,b in p.groupby(['public_GSE','cancer','response_lineage']):
        rows.append(dict(zip(['public_GSE','cancer','response_lineage'],key))|dict(target_positions=len(b),paired_estimable=int(b.estimable.sum()),both_readouts_positive_gain=int(b.both_readouts_positive.sum()),all_six_fold_gain_positive=int(b.all_six_folds_positive.sum()),median_max_signature_Jaccard=float(b.overlap.median()),median_signature_covered_by_other_union=float(b.covered_by_other_union.median())))
    paired_summary=pd.DataFrame(rows);assert len(paired_summary)==4;paired_summary.to_csv(out/'source_cancer_lineage_paired_and_overlap_summary.tsv',sep='\t',index=False)
    common=[]
    for key,b in p[p.estimable].groupby(['public_GSE','response_lineage','proxy']):
        if b.cancer.nunique()==2:common.append(dict(public_GSE=key[0],response_lineage=key[1],both_cancers_positive=bool(b.both_readouts_positive.all()),all_cancer_folds_positive=bool(b.all_six_folds_positive.all())))
    rows=[]
    for key,b in pd.DataFrame(common).groupby(['public_GSE','response_lineage']):
        rows.append(dict(public_GSE=key[0],response_lineage=key[1],shared_physical_measurements=len(b),both_cancers_both_readouts_gain_positive=int(b.both_cancers_positive.sum()),all_cancers_all_folds_gain_positive=int(b.all_cancer_folds_positive.sum())))
    cross=pd.DataFrame(rows);assert len(cross)==2;cross.to_csv(out/'source_two_cancer_common_measurement_descriptors.tsv',sep='\t',index=False)
    cv=pd.DataFrame(coverage).groupby(['public_GSE','cancer','status'],as_index=False).technical_fields.sum();cv.to_csv(out/'source_cancer_field_scope_counts.tsv',sep='\t',index=False)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,3,figsize=(15.7,5.2),gridspec_kw={'width_ratios':[1,1.2,1.2]})
    groups=[('LUAD','T cell'),('LUAD','Macrophage'),('pleural_mesothelioma','T cell'),('pleural_mesothelioma','Macrophage')]
    labels=['LUAD | T','LUAD | Mac','Mesothelioma | T','Mesothelioma | Mac'];pos=np.arange(4)
    raw=e[e.readout.eq('RAW')].set_index(['cancer','response_lineage']).loc[groups]
    axes[0].barh(pos,raw.full_family_positions,color='#dfe5e8',label='Frozen positions / readout');axes[0].barh(pos,raw.estimable_positions,color='#2a7f9e',label='Estimable / readout')
    for j,row in enumerate(raw.to_dict('records')):axes[0].text(row['full_family_positions']+.4,j,'%d/%d'%(row['estimable_positions'],row['full_family_positions']),va='center')
    axes[0].set_xlim(0,30);axes[0].set_yticks(pos,labels);axes[0].set_xlabel('Correlated measurement positions');axes[0].set_title('Coverage')
    for readout,offset,color,label in [('RAW',-.12,'#347eae','Raw'),('NEGATIVE_EXCESS',.12,'#bd6b31','Negative-background excess')]:
        z=e[e.readout.eq(readout)].set_index(['cancer','response_lineage']).loc[groups]
        for ax,col in [(axes[1],'gain_R2'),(axes[2],'multiview_R2')]:
            mid=z['median_'+col].to_numpy();lo=z['q25_'+col].to_numpy();hi=z['q75_'+col].to_numpy()
            ax.errorbar(mid,pos+offset,xerr=np.vstack([mid-lo,hi-mid]),fmt='o',color=color,capsize=3,label=label)
    for ax in axes[1:]:ax.set_yticks(pos,labels);ax.axvline(0,color='#777777',lw=.8)
    axes[1].set_title('Added opposite-lineage state information');axes[1].set_xlabel('Heldout gain in R²: median and IQR')
    axes[2].set_title('Complete model prediction');axes[2].set_xlabel('Heldout R²: median and IQR')
    for ax in axes:ax.invert_yaxis()
    fig.suptitle('GSE299786 | spatial prediction with background-only controls',fontsize=15,y=.98)
    h0,l0=axes[0].get_legend_handles_labels();h1,l1=axes[1].get_legend_handles_labels()
    fig.legend(h0+h1,l0+l1,loc='upper center',bbox_to_anchor=(.52,.91),ncol=4,frameon=False,fontsize=9)
    fig.text(.025,.035,'Equal FOV test weight; 3 spatial regions and a far-neighbour buffer. IQR is not a confidence interval. No P/FDR or causal claim.\nBaseline retains quality, author subtype and composition; other receiver-state proxies are excluded.\nOne source; reviewed independent patients = 0. Correlated partial proxies are not distinct validated states. Missing positions remain unestimable.',fontsize=9)
    fig.subplots_adjust(left=.12,right=.985,top=.76,bottom=.27,wspace=.45)
    for ext in ['png','svg']:fig.savefig(out/('source_blocked_state_prediction_summary.'+ext),dpi=180)
    plt.close(fig)
    receipt=dict(status='COMPLETE_SOURCE_ONLY_REPORT_PENDING_VISUAL_AND_SCIENTIFIC_REVIEW',completed_UTC=datetime.now(timezone.utc).isoformat(),report_source_sha256=sha(__file__),audit_sha256=sha(audit),files={p.name:dict(bytes=p.stat().st_size,sha256=sha(p)) for p in out.iterdir() if p.is_file()},contains_no_individual_or_programme_ID_gene_weight_matrix=True)
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'status':receipt['status'],'files':len(receipt['files'])}))


if __name__=='__main__':main()
