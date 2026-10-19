"""GSE/cancer/lineage-only capture report; no cells, gene lists or weights."""
import hashlib,json,os,shutil
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    parent=Path('/dev/shm/pancancer-native-state-capture-20261020-v2/private_output')
    audit=Path('/dev/shm/pancancer-native-state-capture-audit-20261020-v2/private_output/private_terminal.json')
    out=Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert out.parent.parent==Path('/dev/shm') and not out.exists() and shutil.disk_usage('/dev/shm').free>512*1024**2
    assert sha(audit)=='95693e8e457d7042c93694563504ad56ca61b2e903c5dee8d20a31071510999c'
    terminal=json.loads((parent/'private_terminal.json').read_text())
    name='private_all_field_proxy_capture_and_quality.tsv.gz';assert sha(parent/name)==terminal['files'][name]['sha256']
    f=pd.read_csv(parent/name,sep='\t',dtype={'sample_id':str,'FOV':str})
    rows=[];points=[];joints=[]
    for (gse,cancer,major),g in f.groupby(['dataset_id','cancer','major'],sort=True):
        q=g[g.native_capture_qualified];v=q.groupby('proxy_group_id').residual_variance_fraction.median()
        rows.append(dict(GSE=gse,cancer=cancer,author_lineage=major,tested_proxy_measurements=g.proxy_group_id.nunique(),
            proxies_with_any_qualified_field=q.proxy_group_id.nunique(),proxies_with_atleast5_qualified_fields=int(q.groupby('proxy_group_id').size().ge(5).sum()),
            technical_fields=len(g[['sample_id','FOV']].drop_duplicates()),fields_with_any_qualified_proxy=len(q[['sample_id','FOV']].drop_duplicates()),
            median_of_proxy_residual_variance_fraction=float(v.median()) if len(v) else np.nan,reviewed_independent_patients=0))
        for value in v:points.append(dict(GSE=gse,cancer=cancer,author_lineage=major,median_residual_variance_fraction=float(value)))
    for (gse,cancer),g in f.groupby(['dataset_id','cancer'],sort=True):
        sets={m:set(map(tuple,g[g.major.eq(m)&g.native_capture_qualified][['sample_id','FOV']].to_numpy())) for m in ['T cell','Macrophage']}
        joints.append(dict(GSE=gse,cancer=cancer,technical_fields=len(g[['sample_id','FOV']].drop_duplicates()),
            fields_with_both_lineages_qualified=len(sets['T cell']&sets['Macrophage']),reviewed_independent_patients=0))
    s=pd.DataFrame(rows);p=pd.DataFrame(points);j=pd.DataFrame(joints);assert len(s)==4 and len(j)==2 and len(p)==47
    status=f.groupby(['dataset_id','cancer','major','status'],sort=True).size().reset_index(name='field_proxy_positions')
    out.mkdir(mode=0o700)
    s.to_csv(out/'source_cancer_lineage_capture_counts.tsv',sep='\t',index=False)
    j.to_csv(out/'source_cancer_joint_capture_counts.tsv',sep='\t',index=False)
    p.to_csv(out/'anonymous_source_proxy_variance_fraction_descriptors.tsv',sep='\t',index=False)
    status.to_csv(out/'source_cancer_lineage_reason_counts.tsv',sep='\t',index=False)
    # Carry the whole input denominator, including unavailable author/geometry sources.
    prior=Path('/dev/shm/pancancer-native-state-probe-match-20261020-v1/private_output')
    a=pd.read_csv(prior/'private_complete_source_target_capture_summary.tsv',sep='\t')
    a.to_csv(out/'all_source_author_target_scope_counts.tsv',sep='\t',index=False)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none'})
    fig,(ax,bx)=plt.subplots(1,2,figsize=(12.2,5.7));colors={'T cell':'#3171a0','Macrophage':'#4f8d69'}
    cancers=['LUAD','pleural_mesothelioma'];labels={'LUAD':'Lung adenocarcinoma','pleural_mesothelioma':'Pleural mesothelioma'}
    xpos=np.arange(2);width=.23
    for k,m in enumerate(['T cell','Macrophage']):
        values=[int(s[s.cancer.eq(c)&s.author_lineage.eq(m)].fields_with_any_qualified_proxy.iloc[0]) for c in cancers]
        bars=ax.bar(xpos+(k-1)*width,values,width,label=m,color=colors[m]);ax.bar_label(bars,padding=3,fontsize=9)
    joint=[int(j[j.cancer.eq(c)].fields_with_both_lineages_qualified.iloc[0]) for c in cancers]
    bars=ax.bar(xpos+width,joint,width,label='Both lineages',color='#c18b54');ax.bar_label(bars,padding=3,fontsize=9)
    totals=[int(j[j.cancer.eq(c)].technical_fields.iloc[0]) for c in cancers]
    ax.set_xticks(xpos,[labels[c]+'\n(all FOVs=%d)'%n for c,n in zip(cancers,totals)]);ax.set(ylabel='Technical fields with qualified capture',ylim=(0,90),title='Matched native programme capture');ax.legend(fontsize=8)
    categories=[(c,m) for c in cancers for m in ['T cell','Macrophage']];rng=np.random.default_rng(2026102003)
    for i,(c,m) in enumerate(categories):
        v=p[p.cancer.eq(c)&p.author_lineage.eq(m)].median_residual_variance_fraction.to_numpy()
        bx.scatter(i+rng.uniform(-.12,.12,len(v)),v,s=24,c=colors[m],alpha=.65,linewidths=0)
        bx.plot([i-.2,i+.2],[np.median(v)]*2,color='black',lw=1.5)
    bx.set_xticks(range(4),['LUAD\nT','LUAD\nMac','Mesothelioma\nT','Mesothelioma\nMac'])
    bx.set(ylabel='Residual variance / original score variance',ylim=(0,1.03),title='After observed panel quality and author subtype')
    fig.suptitle('GSE299786 | 4 native RNA-observed spatial inputs',fontsize=14)
    fig.tight_layout(rect=(0,.23,1,.94))
    fig.text(.035,.035,'Capture: >=30 lineage cells, >=10 cells with 3 signature probes, identifiable full quality/subtype design and score variation.\nFOVs are technical fields, not independent patients. Points are correlated programme measurements, not independent states.\nResidual variation is not necessarily biological. Metadata-only segmentation rows have unknown RNA.\nNo new spatial interaction, clinical effect or full-signature state validation is claimed.',fontsize=9)
    for ext in ['png','svg']:fig.savefig(out/('source_native_state_capture_and_variation.'+ext),dpi=170)
    plt.close(fig)
    report=dict(status='COMPLETE_SOURCE_ONLY_NATIVE_CAPTURE_FIGURE_PENDING_ACTUAL_PIXEL_REVIEW',completed_UTC=datetime.now(timezone.utc).isoformat(),
         producer_sha256=sha(__file__),capture_terminal_sha256=sha(parent/'private_terminal.json'),audit_terminal_sha256=sha(audit),
         source_cancer_lineage_rows=len(s),source_cancer_joint_rows=len(j),anonymous_source_proxy_descriptor_rows=len(p),
         individual_records_probe_definitions_gene_lists_weights_or_matrices_output=False,scientific_outputs_public=False,
         files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()})
    assert sum(v['bytes'] for v in report['files'].values())<2*1024**2
    (out/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['status','source_cancer_lineage_rows','source_cancer_joint_rows','anonymous_source_proxy_descriptor_rows']}))


if __name__=='__main__':main()
