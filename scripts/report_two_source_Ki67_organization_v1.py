"""Private two-source descriptive forest, with explicit inference limits."""
import argparse
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def select(rows,context):
    result=[r for r in rows if r['regulatory_context']==context and r['left_marker']=='Ki-67' and r['right_marker']=='Ki-67']
    assert len(result)==1
    return result[0]


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--first',type=Path,required=True);p.add_argument('--first-sha',required=True)
    p.add_argument('--second',type=Path,required=True);p.add_argument('--second-sha',required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    assert hashlib.sha256(a.first.read_bytes()).hexdigest()==a.first_sha
    assert hashlib.sha256(a.second.read_bytes()).hexdigest()==a.second_sha
    first=json.loads(a.first.read_text());combined=json.loads(a.second.read_text());second=combined['analysis']
    assert first['independent_numeric_replay']['all_checks_passed']
    assert second['source_permutation_means_independently_replayed']
    assert second['phenotype_audit']['author_probability_recipes_all_passed']
    assert second['context_audit']['selected_context_accepted']
    a.output.mkdir(mode=0o700,parents=False,exist_ok=False)
    roles=[('author_regulatory','regulatory_Ki67_CD8_Ki67','Regulatory T / CD8'),
        ('author_conventional_CD4','conventional_CD4_Ki67_CD8_Ki67','Conventional CD4 / CD8')]
    fig,axes=plt.subplots(1,2,figsize=(15.5,5.8),sharey=True,sharex=True)
    public_rows=[]
    modes=[('ROI_parent_shuffle',-.1,'#222222','o','Within ROI (main)'),
        ('ROI_parent_grid4_shuffle',.1,'#247c9b','s','Within grid (sensitivity)')]
    for i,(context,comparison,label) in enumerate(roles):
        for mode,offset,color,marker,caption in modes:
            x=select(first['family'][mode],context)
            y=next(r for r in second['family'][mode] if r['comparison']==comparison)
            for ax,row,kind in [(axes[0],x,'case'),(axes[1],y,'material')]:
                if row.get('effect') is None:continue
                value=row['effect'];location=i+offset
                if kind=='case':
                    ci=row.get('case_bootstrap_95_CI')
                    if ci:ax.hlines(location,*ci,color=color,lw=1.6)
                ax.plot(value,location,marker=marker,color=color,ms=6,label=caption if i==0 else None)
                units=row['source_case_units'] if kind=='case' else row['source_materials']
                ax.text(1.03,location,f'n={units}; BH={row["BH"]:.3g}',transform=ax.get_yaxis_transform(),
                    ha='left',va='center',fontsize=8)
            public_rows.append(dict(comparison=comparison,mode=mode,first_source=x,second_source=y))
    bounds=[0.]
    for row in public_rows:
        for name in ['first_source','second_source']:
            record=row[name]
            if record.get('effect') is not None:bounds.append(record['effect'])
            bounds.extend(record.get('case_bootstrap_95_CI') or [])
    lower,upper=min(bounds),max(bounds);margin=.1*max(.1,upper-lower)
    for ax in axes:
        ax.axvline(0,color='#999999',lw=.8);ax.set_ylim(1.5,-.5)
        ax.set_xlim(lower-margin,upper+margin)
        ax.set_xlabel('Conditional Ki-67 rank covariance')
        ax.spines[['top','right']].set_visible(False)
    axes[0].set_yticks(range(2),[r[2] for r in roles]);axes[0].set_title('Source 6838169 (2022): invasive margin\nReviewed case units; 95% case-bootstrap CI')
    axes[1].set_title('Source 5903179 (2022): selected materials\nMaterial units; no patient CI')
    axes[1].legend(frameon=False,loc='lower right',fontsize=8)
    fig.suptitle('2026-10-14 immune-atlas CM02 v1-inspired Ki-67 state exploration',fontsize=13)
    fig.subplots_adjust(left=.19,right=.78,top=.78,bottom=.28,wspace=.65)
    fig.text(.19,.08,'Descriptive comparison across different tissue/therapy contexts; no pooled estimate.\nPositive covariance includes high-high and low-low. No direct communication or clinical inference.\nUnknown repeated materials can affect second-source randomization; patient independence is unreviewed.',fontsize=9)
    fig.savefig(a.output/'two_source_Ki67_organization.png',dpi=180);fig.savefig(a.output/'two_source_Ki67_organization.svg');plt.close(fig)
    differences=[r for rows in second['family'].values() for r in rows if r['comparison']=='regulatory_minus_conventional_shared_materials']
    difference_description=[{k:r.get(k) for k in ['comparison','source_materials','effect','positive_materials','negative_materials','leave_one_material_effect_range']} for r in differences]
    report=dict(status='PRIVATE_TWO_SOURCE_DESCRIPTIVE_REVIEW_NOT_CLINICAL_VALIDATION',source_rows=public_rows,
        difference_description_only=difference_description,
        specificity_inference=dict(formal_equal_nonzero_association_null_not_tested=True,
            shared_CD8_randomization_destroys_both_parent_associations=True,
            original_third_P_is_joint_randomization_not_formal_parent_superiority=True,
            no_priority_promotion_from_third_P_or_BH=True),
        clinical_models_not_run=True,patient_independence_in_second_source_unreviewed=True,
        nuclear_background_and_segmentation_artifacts_not_excluded=True,
        original_CM02_weights_not_projected_or_validated=True,
        source_hashes=dict(first=a.first_sha,second=a.second_sha))
    (a.output/'private_two_source_review.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    lines=['2026-10-14免疫图谱CM02（v1，T_whole_TME_primary_tumor）启发的状态组织探索','']
    for context,comparison,label in roles:
        source1=[select(first['family'][mode],context) for mode,_,_,_,_ in modes]
        source2=[next(r for r in second['family'][mode] if r['comparison']==comparison) for mode,_,_,_,_ in modes]
        for name,records,count_field in [('第一来源',source1,'source_case_units'),('第二来源',source2,'source_materials')]:
            supported=all(r.get('BH') is not None and r['BH']<=.05 and r['effect']>0 for r in records)
            lines.append(f'{label}；{name}：主/网格正向且均BH≤0.05={supported}；单位数={records[0][count_field]}。')
            for caption,record in zip(['主','网格'],records):
                lines.append(f'  {caption}：效应={record.get("effect")}；BH={record.get("BH")}。')
    lines+=['','正协方差同时包含高高与低低，不单独证明Ki-67阳性细胞富集。来源之间的治疗/组织背景需另审。尚无原CM02权重映射或临床验证。',
        '第二来源的父类协方差差值只作描述。原第三置换将两条关系都打散，其零分布不是允许共同非零关系、仅检验两者相等的零假设；不能用该P/BH宣布调节T特异优势。配对病例差异仍待合适单位、区间及检验。',
        '材料独立性尚未审阅；第二来源不提供患者区间或临床模型。同通道/核背景/分割影响尚未排除。',
        '下一步：将高高状态贡献、核/背景敏感性、肿瘤近邻及治疗/取材分层另行固定；先核患者、阶段和终点资格，再决定临床增量。保留原注释和已结束家族。']
    (a.output/'private_two_source_report.txt').write_text('\n'.join(lines)+'\n')
    print(json.dumps(dict(status='PRIVATE_SOURCE_COMPARISON_WRITTEN',no_individual_records=True,
        third_comparison_not_presented_as_formal_superiority=True)))


if __name__=='__main__':main()
