"""Question-specific library/source eligibility from documented preparation.

This is the project's proposed screening policy. It does not claim that source
authors used this table, infer sorting from cell proportions, alter an archive,
or certify patient independence, clinical eligibility or identity annotations.
"""

SORTING={'unknown','unselected','viability_only','CD45_positive','CD45_negative','CD45_fraction_unresolved','lymphocyte_enriched','CD3_positive','T_enriched','CD8_positive','other_lineage_selected','mixed_sorted_fractions'}
MATERIAL={'primary_tumor','tumor_unspecified','metastatic_tumor','normal_adjacent','blood','bone_marrow','lymph_node_unknown','fluid','unknown','mixed','cell_product'}
MODALITY={'scRNA','snRNA','unknown','mixed'}

def qualify_library(review):
    sorting=review.get('experimental_sorting','unknown');material=review.get('material','unknown');modality=review.get('RNA_modality','unknown')
    if sorting not in SORTING or material not in MATERIAL or modality not in MODALITY:raise ValueError('Source/library category requires explicit mapping')
    retention=review.get('retained_data_scope','unknown');product=bool(review.get('cell_product',False)) or material=='cell_product';stimulated=review.get('ex_vivo_stimulation','unknown')=='yes'
    reasons=[]
    if product:reasons.append('CAR/离体制造产品，不代表原位肿瘤组织')
    if stimulated:reasons.append('离体刺激或培养，单独研究实验状态')
    if retention=='selected_nonmalignant_subset':reasons.append('提供的数据已选择非恶性区室，缺少完整TME分母')
    if retention=='selected_other_lineage':reasons.append('提供的数据仅保留指定非T谱系')
    if sorting in ['CD45_positive','lymphocyte_enriched']:reasons.append('免疫/淋巴细胞富集，不能直接比较全TME多谱系比例')
    if sorting in ['CD3_positive','T_enriched','CD8_positive']:reasons.append('T/CD8富集，缺少其他细胞区室的观察机会')
    if sorting=='other_lineage_selected':reasons.append('指定谱系分选，缺少完整TME捕获')
    if sorting=='CD45_negative':reasons.append('CD45阴性分选，免疫区室缺少观察机会')
    if sorting=='CD45_fraction_unresolved':reasons.append('已确认CD45正负分选设计，但本库正负分支未明；不能作为自然全TME比例')
    if sorting=='mixed_sorted_fractions' and not review.get('mixing_weights_reconstructable',False):reasons.append('分选区室混合比例未恢复，不能使用未加权全TME比例')
    if material in ['blood','bone_marrow','fluid']:reasons.append({'blood':'血液/PBMC另建集合','bone_marrow':'骨髓生态另建集合','fluid':'体液/胸腹水另建集合'}[material])
    if material=='normal_adjacent':reasons.append('正常/癌旁作为独立对照集合，不混入肿瘤主发现集')
    if modality=='snRNA':reasons.append('snRNA单独分层，暂不直接混入scRNA主发现集')
    missing=[]
    if sorting=='unknown':missing.append('具体分选/富集方法')
    if modality in ['unknown','mixed']:missing.append('scRNA/snRNA及其文库对应')
    if material in ['unknown','mixed','lymph_node_unknown']:missing.append('明确材料/解剖部位与肿瘤受累')
    if retention=='unknown':missing.append('当前输入文件是否保留完整区室')
    if review.get('mixed_study_library_branch_unresolved',False):missing.append('混合研究中本样本的总组织/富集分支')
    if not review.get('actual_sample_or_library_bridge',False):missing.append('实际样本/物理文库对应')
    whole='排除' if reasons else '待核实' if missing else '条件候选'
    if product:tstatus,tscope='单独分析','CAR-T产品T细胞状态'
    elif stimulated:tstatus,tscope='单独分析','离体刺激/培养T细胞状态'
    elif sorting=='CD45_fraction_unresolved':tstatus,tscope='待核实','先恢复CD45正负文库分支，再限定可比较的T细胞捕获'
    elif material in ['blood','bone_marrow','fluid','normal_adjacent']:
        tstatus,tscope='单独分析',{'blood':'血液T细胞','bone_marrow':'骨髓T细胞','fluid':'体液T细胞','normal_adjacent':'正常/癌旁T细胞'}[material]
    elif sorting=='CD8_positive':tstatus,tscope='条件候选','限定CD8内部，不能代表全部T细胞'
    elif sorting in ['other_lineage_selected','CD45_negative'] or retention=='selected_other_lineage':tstatus,tscope='排除','未提供可比较的全T细胞捕获'
    elif material in ['unknown','mixed','lymph_node_unknown'] or modality in ['unknown','mixed'] or not review.get('actual_sample_or_library_bridge',False):tstatus,tscope='待核实','先核材料、模态与样本/文库桥'
    else:tstatus,tscope='条件候选','同材料、阶段及捕获范围内的T细胞'
    if modality=='snRNA':tscope+='；snRNA分支'
    if sorting=='mixed_sorted_fractions':tscope+='；需保留分选层并核对T亚群捕获'
    if review.get('T_subtype_capture_comparability','unknown')!='accepted' and tstatus=='条件候选':tscope+='；T亚群捕获可比性待核'
    return {'whole_TME_scRNA_main':whole,'whole_TME_primary_reason':reasons[0] if reasons else '；'.join(missing) if missing else '流程候选，仍需身份/QC/来源别名与重复取样资格',
            'whole_TME_all_reasons':reasons,'critical_missing_fields':missing,'T_internal_CM':tstatus,'T_internal_scope':tscope,
            'archive_retained':True,'patient_independence_newly_certified':False,'clinical_eligibility_newly_certified':False,
            'sorting_inferred_from_cell_proportions':False,'policy_is_project_proposal_not_author_protocol_table':True}

def normalize_library_alias(value):
    """Normalize separator spelling while preserving identifier token boundaries.

    BC2 1 and BC21 must remain different. Dataset-specific semantic aliases,
    prefixes, patient numbering and timepoints require separate evidence.
    """
    import re
    return re.sub(r'[^a-z0-9]+','_',str(value).lower()).strip('_')

def resolve_library_aliases(query_aliases,reference_aliases):
    """Return all exact normalized hits; ambiguity is never silently resolved."""
    keys={normalize_library_alias(v) for v in query_aliases}
    hits=sorted(k for k,aliases in reference_aliases.items() if keys & {normalize_library_alias(v) for v in aliases})
    return {'reference_keys':hits,'status':'unique' if len(hits)==1 else 'ambiguous' if hits else 'unmatched'}
