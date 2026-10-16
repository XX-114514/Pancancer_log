"""Private complete evidence table, provenance, diagnostic coverage/forests.

No new hypothesis tests or candidate substitutions. Wait for numerical audit.
"""
import hashlib,json,math,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()

def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);I=R/'NSCLC_automatic_ligand_response_results_20261016_v1'
    launch=json.loads((R/'HNSCC_receiver_followup_job_control_20261016_v2/audit_NSCLC_automatic_ligand_response_results_v1/launch_receipt.json').read_text());deadline=time.monotonic()+100000
    while not (I/'private_independent_audit.json').exists():
        p=Path('/proc')/str(launch['pid'])/'stat';assert p.exists();s=p.read_text().rsplit(')',1)[1].split()
        assert s[19]==str(launch['proc_start_ticks']) and s[0]!='Z' and time.monotonic()<deadline;time.sleep(20)
    audit=json.loads((I/'private_independent_audit.json').read_text());assert audit['status'].startswith('PASS_COMPLETE_UNIQUE_FAMILY_ALL_HC3_')
    path=I/'private_all_unique_ligand_response_results.tsv';assert sha(path)==audit['result_sha256'];t=pd.read_csv(path,sep='\t');assert len(t)==44528
    O=R/'NSCLC_automatic_ligand_response_review_20261016_v1';assert not O.exists();O.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),result_audit_sha256=sha(I/'private_independent_audit.json'),
        scope='All original unique hypotheses, all estimable forest records and all input unavailable reasons, no new P/BH.',
        classes='Source-limited statistical candidate only if original primaryBH<.05, finite pointwiseCI excludes0, >=95% cases LOO estimable and>=95% valid LOO signs agree. BH with incomplete robustness separate; all other evaluated rows insufficient full-family evidence; unavailable retains original reason.',
        joins='Original full-signature clinical family and original two-end coordination are descriptive annotations, NOT jointFDR or validation of the current masked response.',
        provenance='Original CM context means historical programme provenance, not acceptance that native author subtypes map to original172 identities or old CM activities reproduce.',
        figures='All estimated rows, paginated40 per forest; no significance-only omission. Overview full planned/evaluated coverage plus raw/adjusted effects. Pointwise CI and dependent tests labelled.',
        boundary='One NSCLC source. No independent-patient, cross-source, unique-receptor, protein assembly, clinical-prediction or causal claim.')
    (O/'review_contract.json').write_text(json.dumps(plan,indent=2)+'\n')
    V=Path(json.loads((R/'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root'])
    context_path=V/'original_CM_programme_pair_associations_20261015_v2/all_context_coverage_reasons.tsv'
    contexts=pd.read_csv(context_path,sep='\t',keep_default_na=False).set_index('context');assert set(t.context)<=set(contexts.index)
    cp=R/'NSCLC_frozen_signature_clinical_results_20261016_v1';ca=json.loads((cp/'private_independent_audit.json').read_text())
    cpath=cp/'private_all_planned_signature_clinical_results.tsv';assert sha(cpath)==ca['result_sha256']
    clinical=pd.read_csv(cpath,sep='\t').set_index(['context','program','type','sub_cell_type'])
    qp=R/'NSCLC_two_end_coordination_results_20261016_v1';qa=json.loads((qp/'private_independent_audit.json').read_text())
    qpath=qp/'private_all572_two_end_results.tsv';assert sha(qpath)==qa['result_sha256'];coord=pd.read_csv(qpath,sep='\t').set_index('record')
    evaluated=t.sequential_wild_p.notna();bh=evaluated&t.BH_q_complete_unique_family.lt(.05)
    ci=np.isfinite(t.CI_low)&np.isfinite(t.CI_high)&(t.CI_low*t.CI_high).gt(0)
    stable=t.LOO_estimable.ge(.95*t.cases)&t.LOO_same_sign.ge(.95*t.LOO_estimable)
    robust=bh&ci&stable
    t['evidence_class']=np.where(~evaluated,'NOT_ESTIMABLE',np.where(robust,'SOURCE_LIMITED_STATISTICAL_CANDIDATE',np.where(bh,'BH_SUPPORTED_ROBUSTNESS_INCOMPLETE','INSUFFICIENT_COMPLETE_FAMILY_EVIDENCE')))
    annotations=[]
    for r in t.itertuples():
        c=contexts.loc[r.context];receiver_type=r.receiver_type;sender_type='typeB' if receiver_type=='typeA' else 'typeA'
        rc=clinical.loc[(r.context,r.program,receiver_type,r.receiver_author_subtype)]
        sc=clinical.loc[(r.context,r.program,sender_type,r.sender_author_subtype)]
        q=coord.loc[r.parent_pair]
        cm_names=';'.join('2026年10月'+s+'（v20261010）' for s in str(c.CMs).split(';') if s)
        annotations.append(dict(record=r.record,training_gse=c.gse_id,training_stratum=c.stratum,original_CM_provenance=cm_names,
            original_receiver_endpoint=c.endpoint_A if receiver_type=='typeA' else c.endpoint_B,
            original_sender_endpoint=c.endpoint_B if sender_type=='typeB' else c.endpoint_A,
            full_receiver_signature_clinical_beta=rc.beta_MPR_minus_nonMPR,full_receiver_signature_clinical_q=rc.BH_q_full_planned_family,
            full_sender_signature_clinical_beta=sc.beta_MPR_minus_nonMPR,full_sender_signature_clinical_q=sc.BH_q_full_planned_family,
            original_two_end_coordination_rho=q.partial_rank_rho,original_two_end_coordination_q=q.BH_q_full572,
            descriptive_combination_not_joint_FDR=True,native_author_subtype_to_old_node_mapping_validated=False))
    ann=pd.DataFrame(annotations);t=t.merge(ann,on='record',validate='one_to_one');assert len(t)==44528
    t.to_csv(O/'private_all44528_annotated_evidence.tsv.gz',sep='\t',index=False,na_rep='NA',compression='gzip')
    t[t.BH_q_complete_unique_family.lt(.05)].to_csv(O/'private_all_primary_BH_supported_rows.tsv',sep='\t',index=False,na_rep='NA')
    t[t.evidence_class.eq('SOURCE_LIMITED_STATISTICAL_CANDIDATE')].to_csv(O/'private_all_source_limited_statistical_candidates.tsv',sep='\t',index=False,na_rep='NA')
    provenance=contexts.loc[sorted(t.context.unique()),['gse_id','CMs','stratum','endpoint_A','endpoint_B']].reset_index()
    provenance['CMs']=provenance.CMs.map(lambda s:';'.join('2026年10月'+v+'（v20261010）' for v in str(s).split(';') if v))
    provenance.to_csv(O/'private_all_training_context_provenance.tsv',sep='\t',index=False)
    summaries=[]
    for key,f in t.groupby(['context','program','receiver_type'],sort=True):
        summaries.append(dict(context=key[0],program=key[1],receiver_type=key[2],planned=len(f),evaluated=int(f.sequential_wild_p.notna().sum()),
            primary_BH_rows=int(f.BH_q_complete_unique_family.lt(.05).sum()),BY_rows=int(f.BY_q_dependence_sensitivity.lt(.05).sum()),
            robust_rows=int(f.evidence_class.eq('SOURCE_LIMITED_STATISTICAL_CANDIDATE').sum())))
    summary=pd.DataFrame(summaries);assert len(summary)==16 and summary.planned.sum()==len(t)
    summary.to_csv(O/'private_all16_context_direction_summaries.tsv',sep='\t',index=False)
    estimated=t[t.sequential_wild_p.notna()].sort_values(['context','program','receiver_type','record']).copy();assert len(estimated)==audit['all_estimable_points']
    columns=['record','context','program','receiver_type','sender_author_subtype','receiver_author_subtype','ligand','cases','partial_rank_rho','CI_low','CI_high','BH_q_complete_unique_family','BY_q_dependence_sensitivity','raw_spearman','sender_programme_adjusted_rho','pathology_adjusted_rho','evidence_class']
    estimated[columns].to_csv(O/'private_all_estimable_forest_data.tsv',sep='\t',index=False,na_rep='NA')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'svg.fonttype':'none'})
    labels=[f'C{r.context}/{r.program}/{r.receiver_type}' for r in summary.itertuples()]
    fig,axs=plt.subplots(1,2,figsize=(14,6.6));yy=np.arange(len(summary))
    axs[0].barh(yy,summary.planned,color='#d9d9d9',label='All planned hypotheses');axs[0].barh(yy,summary.evaluated,color='#2f627d',label='Evaluated')
    axs[0].set_yticks(yy,labels);axs[0].invert_yaxis();axs[0].set_xlabel('Unique hypotheses');axs[0].set_title('Input coverage by reference programme/direction');axs[0].legend(frameon=False,fontsize=8)
    sig=estimated.BH_q_complete_unique_family.lt(.05)
    axs[1].scatter(estimated.loc[~sig,'raw_spearman'],estimated.loc[~sig,'partial_rank_rho'],s=6,color='#a0a0a0',alpha=.45,rasterized=True)
    axs[1].scatter(estimated.loc[sig,'raw_spearman'],estimated.loc[sig,'partial_rank_rho'],s=13,color='#a83b32',label='Primary BH < .05',rasterized=True)
    axs[1].plot([-1,1],[-1,1],ls='--',lw=.8,color='#777777');axs[1].axhline(0,lw=.6,color='#999999');axs[1].axvline(0,lw=.6,color='#999999')
    axs[1].set(xlim=(-1,1),ylim=(-1,1),xlabel='Unadjusted rank correlation',ylabel='Fully controlled partial rank correlation',title='All evaluated ligand-programme hypotheses')
    if sig.any():axs[1].legend(frameon=False,fontsize=8)
    fig.suptitle(f'GSE243013: {len(t):,} planned / {len(estimated):,} evaluated; one source',fontsize=12)
    fig.tight_layout();fig.savefig(O/'private_coverage_and_adjustment.png',dpi=160);fig.savefig(O/'private_coverage_and_adjustment.svg');plt.close(fig)
    forest=O/'private_all_estimable_forest_pages';forest.mkdir(mode=0o700);pages=[]
    for page,start in enumerate(range(0,len(estimated),40),1):
        f=estimated.iloc[start:start+40];fig,ax=plt.subplots(figsize=(15,11));y=np.arange(len(f))
        for position,r in enumerate(f.itertuples()):
            color='#a83b32' if r.BH_q_complete_unique_family<.05 else '#67727b'
            if pd.notna(r.CI_low) and pd.notna(r.CI_high):ax.hlines(position,r.CI_low,r.CI_high,color=color,lw=1)
            ax.scatter(r.partial_rank_rho,position,color=color,s=17)
        lab=[f'{r.record} C{r.context}/{r.program}/{r.receiver_type}: {r.sender_author_subtype} -> {r.receiver_author_subtype} | {r.ligand} (n={int(r.cases)}, q={r.BH_q_complete_unique_family:.3g})' for r in f.itertuples()]
        ax.set_yticks(y,lab,fontsize=6.3);ax.invert_yaxis();ax.axvline(0,color='#555555',ls='--',lw=.8)
        ax.set(xlim=(-1.03,1.03),xlabel='Partial rank correlation; pointwise 95% case-bootstrap interval',title=f'All evaluated hypotheses, page {page}; red = primary complete-family BH < .05')
        fig.subplots_adjust(left=.56,right=.98,bottom=.08,top=.94);png=forest/f'forest_{page:03d}.png';fig.savefig(png,dpi=140);plt.close(fig)
        pages.append(dict(page=page,file=str(png.relative_to(O)),records=f.record.tolist()))
    assert [r for p in pages for r in p['records']]==estimated.record.tolist()
    (O/'private_all_forest_page_records.json').write_text(json.dumps(pages,indent=2)+'\n')
    fmt=lambda v:'NA' if pd.isna(v) else f'{v:.6g}'
    report=['# 自动配体—接收程序：完整结果与证据边界','',f'UTC: {datetime.now(timezone.utc).isoformat()}','',
        f'完整不重复假设{len(t)}；可估计{len(estimated)}；主BH记录{int(sig.sum())}；满足BH/区间/留一预定规则{int(robust.sum())}；BY敏感性记录{int(estimated.BY_q_dependence_sensitivity.lt(.05).sum())}。',
        '这些是相关假设记录，不能当独立通讯轴数。不可估计项保留输入原因，不能当RNA阴性或通讯不存在。','',
        '## 原CM与参考程序来源','',provenance.to_string(index=False),'',
        '原参考来自2026年10月CM（v20261010）的上下文表达程序，编号同名MCP不能互换。训练分层含CD45pMACS与unsorted，来自同一CRC来源；不是独立癌种/患者复现。当前作者亚群未一一接受为原172节点，也未证明原CM组成活动在肺癌复现。','',
        '## 证据分级','',t.evidence_class.value_counts().to_string(),'','## 所有来源限定统计候选','',
        '完整候选和BH但稳健性不全记录分别见TSV，全部44528条和全部可估计森林均保留。临床注释属于原完整签名，不能当当前剔除基因后结局的临床检验；与原协调的交集没有新的联合FDR。','']
    supported=t[t.evidence_class.eq('SOURCE_LIMITED_STATISTICAL_CANDIDATE')].sort_values(['BH_q_complete_unique_family','record'])
    for r in supported.itertuples():
        report.append(f'- {r.record}: C{r.context}/{r.program}, {r.sender_author_subtype} [{r.ligand}] -> {r.receiver_author_subtype}; n={int(r.cases)}, rho={fmt(r.partial_rank_rho)}, CI=[{fmt(r.CI_low)}, {fmt(r.CI_high)}], BH={fmt(r.BH_q_complete_unique_family)}, BY={fmt(r.BY_q_dependence_sensitivity)}; possible receptors={r.qualified_receptors}; sender-programme adjusted rho={fmt(r.sender_programme_adjusted_rho)}, pathology adjusted rho={fmt(r.pathology_adjusted_rho)}.')
    report+=['','## 决策','',
        '若有预定规则支持候选：只进入独立功能/来源验证规划，逐一保留所有可能受体、原生RNA、额外控制、临床及训练分层限制，不宣布唯一受体或机制。' if len(supported) else '本轮没有满足全部预定规则的配体—程序候选。停止本家族追加阈值、替换配体/评分或缩减多重校正范围；原临床/两端协调仍保留，下一转其独立功能或队列验证。',
        '全家族顺序wild为模型近似，非精确治疗随机实验。病例未全审为独立患者，未知处理/来源混杂及靶RNA相关不等同蛋白信号；病理额外控制可能有collider偏倚。当前仅一个治疗后免疫来源，不作治疗前预测、跨来源验证、正式细胞状态或耐药机制命名。',
        '图为内部诊断：数据/覆盖/区间均来自已审计完整表。绘图可读性需独立图像审阅后再用于正式发表。']
    (O/'PRIVATE_COMPLETE_SCIENTIFIC_REVIEW.md').write_text('\n'.join(report)+'\n')
    for p in [O/'private_coverage_and_adjustment.png',*forest.glob('*.png')]:
        with Image.open(p) as im:assert im.width>=1000 and im.height>=800 and len(im.convert('RGB').getcolors(maxcolors=1000000) or [])>20
    files={str(p.relative_to(O)):sha(p) for p in O.rglob('*') if p.is_file()}
    result=dict(status='COMPLETE_PRIVATE_EVIDENCE_PROVENANCE_CLASSIFICATION_AND_ALL_ESTIMABLE_FORESTS_NO_NEW_TESTS',finished_UTC=datetime.now(timezone.utc).isoformat(),
        code_sha256=sha(__file__),result_audit_sha256=sha(I/'private_independent_audit.json'),clinical_result_sha256=sha(cpath),coordination_result_sha256=sha(qpath),training_context_sha256=sha(context_path),
        complete_hypothesis_rows=len(t),estimated_rows=len(estimated),forest_pages=len(pages),all_estimated_forest_records_covered=True,
        figures_visually_reviewed=False,new_tests_or_P_or_BH=False,receptor_specific_or_causal_axis_validated=False,files=files)
    (O/'private_terminal.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='files'}))

if __name__=='__main__':main()
