"""Build a versioned, literature-traced review registry; never assign cell labels."""
from __future__ import annotations
import argparse
import copy
import csv
import hashlib
import json
from pathlib import Path

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def words(value):
    return value.split() if value else []

REFERENCES = [
    ('Burclaff2022', 'A Proximal-to-Distal Survey of Healthy Adult Human Small Intestine and Colon Epithelium by Single-Cell Transcriptomics', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC9043569/', 'Human intestine; ISC, TA, Paneth and goblet distinctions; SMOC2 is not an exclusive human ISC marker.'),
    ('HumanColonoid2023', 'Physiological hypoxia improves growth and functional differentiation of human intestinal epithelial organoids', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC9922616/', 'Human colonoids, Fig. 4 absorptive markers; culture context, not an in-vivo pan-cancer classifier.'),
    ('GastricCorpus2023', 'Spatially resolved expression landscape and gene-regulatory network of human gastric corpus epithelium', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC10319429/', 'Human gastric corpus; figure marker lists distinguish surface, chief, parietal and endocrine populations.'),
    ('GastricMucins2004', 'Aberrant epithelial expression of trefoil family factor 2 and mucin 6 in Helicobacter pylori infected gastric antrum, incisura, and body and its association with antralisation', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC1770377/', 'Human biopsy histology; MUC6/TFF2 glandular mucus, including disease-related redistribution; not a unique normal neck identity.'),
    ('Muraro2016', 'A Single-Cell Transcriptome Atlas of the Human Pancreas', 'https://www.hubrecht.eu/app/uploads/2017/11/VanOudenaarden_Key_2016-Muraro_A-single-cell-transcriptome-atlas-of-the-human-pancreas.pdf', 'Human donor pancreas; Fig. 1D and Results pp. 387–388. Isolated islets with exocrine cells; not a TME abundance reference.'),
    ('Kumar2023', 'A spatially resolved single-cell genomic atlas of the adult human breast', 'https://escholarship.org/content/qt81t478n0/qt81t478n0.pdf', 'Normal human breast LumHR, LumSec and basal/myoepithelial; Results/Fig. 3 and spatial validation; KRT5 can occur in LumSec-basal states.'),
    ('Yu2019', 'Single-Cell Transcriptomic Map of the Human and Mouse Bladders', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC6830796/', 'Use human urothelial basal/intermediate/umbrella compartments; UPK-rich differentiated cells do not cover every basal cell.'),
    ('KidneyProtein2014', 'The Kidney Transcriptome and Proteome Defined by Transcriptomics and Antibody-Based Profiling', 'https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0116125', 'Human tissue RNA/protein localization, Table 1 and Fig. 4; not a single-cell clustering algorithm.'),
    ('Lake2023', 'An atlas of healthy and injured cell states and niches in the human kidney', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC10356613/', 'Human nephron segmentation and injury states; Fig. 3d p.589 marker labels in publisher PDF https://www.nature.com/articles/s41586-023-05769-3.pdf; adaptive states may lose mature markers.'),
    ('ThyroidOrganoid2021', 'Adult mouse and human organoids derived from thyroid follicular cells and modeling of Graves’ hyperthyroidism', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC8713972/', 'Human primary tissue/organoid comparisons in Fig. 2C,D and Fig. 4B,D–F; single-cell Fig. 3 is mouse and is not used as a human scRNA validation.'),
    ('FallopianAtlas', 'A cell atlas of the human fallopian tube throughout the menstrual cycle and menopause', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC11698969/', 'Human fallopian secretory epithelium; PAX8/KRT7 and OVGP1, with menstrual/menopause dependence; not all Müllerian organs.'),
    ('GarciaAlonso2021', 'Mapping the temporal and spatial dynamics of the human endometrium in vivo and in vitro', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC8648563/', 'Human endometrium: PAEP secretory phase, SCGB2A2 glandular, FOXJ1/PIFO ciliated; does not support a universal endometrial label from PAEP.'),
    ('Esophagus2022', 'Single-cell RNA-Seq of human esophageal epithelium in homeostasis and allergic inflammation', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC9208762/', 'Human suprabasal differentiation: KRT4/KRT13; an epithelial differentiation direction, not a unique organ barcode.'),
    ('EsophagusMapping2022', 'Mapping the adult human esophagus in vivo and in vitro', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC9720751/', 'Human esophagus: KRT4 mid differentiation, CRNN terminal differentiation; stages need not coexpress in every cell.'),
    ('Cheng2018', 'Transcriptional Programming of Normal and Inflamed Human Epidermis at Single-Cell Resolution', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC6367716/', 'Human epidermis: KRT1/KRT10 spinous and LOR/FLG granular; distinct differentiation stages.'),
]

def curated_specs():
    specs = []
    def add(source, pid, title, context, core, evidence, note, axis='branch', parent='Epithelial'):
        # Each evidence tuple is gene-context-specific; occurrence is not exclusivity.
        blocks = [dict(reference_id=r, genes=words(g), locator=l) for r,g,l in evidence]
        genes = list(dict.fromkeys(g for b in blocks for g in b['genes']))
        specs.append(dict(source_label=source, id=pid, display_name=title, context=context,
                          core_groups=[words(c) for c in core.split(';')], evidence=blocks,
                          genes=genes, note=note, axis=axis, parent=parent))
    add('colonocyte','Epi_colonocyte_like','结肠吸收上皮方向','colon','CA1 CA2;SLC26A3',
        [('HumanColonoid2023','CA1 CA2 SLC26A3','Fig. 4 absorptive marker panel')], '来自人结肠类器官的吸收分化证据；需与原位结肠参考及肠化生对照。')
    add('goblet','Epi_intestinal_goblet_like','肠杯状细胞方向','intestine','MUC2;CLCA1 TFF3',
        [('Burclaff2022','MUC2 CLCA1 TFF3','Goblet cell Results and marker plots')], 'TFF3广泛分泌表达，须结合MUC2/CLCA1；肠化生不等于原发肠来源。')
    add('intestinal_stem','Epi_intestinal_stem_like','肠隐窝干样候选','intestine','LGR5 ASCL2;SLC12A2 RGMB',
        [('Burclaff2022','LGR5 ASCL2 SLC12A2 RGMB OLFM4','ISC/TA/Paneth comparison in Results')], 'OLFM4仅辅助，结肠TA也可表达；SMOC2不进入正向核心，LGR5/ASCL2也须对照Paneth；不是功能性干细胞证明。')
    add('gastric_pit_mucous','Epi_gastric_pit_like','胃表面小凹黏液方向','stomach','MUC5AC TFF1;GKN1 GKN2',
        [('GastricCorpus2023','MUC5AC TFF1 GKN1 GKN2','Main epithelial marker figure')], 'MUC6移至腺体黏液模块；胃型化生和肿瘤可表达该程序。')
    add('gastric_pit_mucous','Epi_gastric_gland_mucous_like','胃腺体黏液方向','stomach','MUC6;TFF2',
        [('GastricMucins2004','MUC6 TFF2','Human biopsy localization and normal/microbial disease comparisons')], '覆盖胃体颈部/胃窦腺体方向，不能强称唯一颈细胞；十二指肠腺和化生亦需排查。')
    add('gastric_chief_parietal','Epi_gastric_chief_like','胃主细胞方向','stomach','PGC PGA3 PGA5;LIPF',
        [('GastricCorpus2023','PGC PGA3 PGA5 LIPF','Main epithelial marker figure')], '与ATP4A/B壁细胞分开。原GIF未纳入；旧符号及人鼠表达差异另核。')
    add('gastric_chief_parietal','Epi_gastric_parietal_like','胃壁细胞方向','stomach','ATP4A;ATP4B',
        [('GastricCorpus2023','ATP4A ATP4B','Main epithelial marker figure')], '胃酸分泌机器方向，不与主细胞消化酶混成一个身份。')
    add('enteroendocrine_neuroendocrine','Epi_neuroendocrine_program','内分泌与神经内分泌分泌程序','','CHGA CHGB;SCG2 SCG3 SCG5',
        [('GastricCorpus2023','CHGA SCG3 SCG5','Endocrine marker figure'),('Muraro2016','CHGB SCG2','Results p.388 endocrine comparisons')], '跨器官程序；须另审内分泌身份和器官。PanNET等肿瘤亦可阳性，不能据此判正常。',axis='state')
    add('pancreatic_acinar','Epi_pancreatic_acinar_like','胰腺腺泡方向','pancreas','PRSS1;PNLIP',
        [('Muraro2016','PRSS1 PNLIP','Fig.1D; Results p.388 coexpression')], '高丰度消化酶易形成环境RNA，需同细胞及来源内证据，孤立PRSS1不足终判。')
    add('pancreatic_ductal','Epi_pancreatic_duct_like','胰腺导管方向','pancreas','KRT19;SPP1',
        [('Muraro2016','KRT19 SPP1','Fig.1D and Results p.388 ductal cluster')], '两者均非胰腺专一；必须审父身份、胆管等对照及胰腺材料/原发部位。')
    add('islet','Epi_pancreatic_islet_like','胰岛内分泌候选','pancreas','INS GCG SST PPY;CHGB SCG2',
        [('Muraro2016','INS GCG SST PPY CHGB SCG2','Fig.1D and Results p.388')], 'INS/GCG/SST/PPY分别指向β/α/δ/PP；此组为候选并集，不要求同一细胞全表达，也不在此强分四类。')
    add('breast_luminal','Epi_breast_LumHR_like','乳腺激素响应管腔方向','breast','ESR1 PGR;ANKRD30A AR',
        [('Kumar2023','ESR1 PGR ANKRD30A AR','Fig.3; smFISH spatial Results pp.7–8')], '正常乳腺LumHR参考；不覆盖全部管腔细胞，不凭受体判肿瘤正常。')
    add('breast_luminal','Epi_breast_LumSec_like','乳腺分泌型管腔方向','breast','LTF SLPI;KRT15 KRT23',
        [('Kumar2023','LTF SLPI KRT15 KRT23 KRT17','Fig.3 spatial/LumSec Results')], '补充原面板漏掉的LumSec；KRT5/KRT15等不应单独把该方向改判基底。')
    add('breast_basal_myoepithelial','Epi_breast_basal_myo_like','乳腺基底肌上皮方向','breast','TP63;ACTA2',
        [('Kumar2023','TP63 ACTA2','Results p.5 basal cells; Fig.3j')], '需独立的上皮/位置证据和周细胞、平滑肌对照；低EPCAM不否定肌上皮。')
    add('prostate_basal_club','Epi_prostate_basal_like','前列腺基底方向','prostate','KRT5 KRT14;TP63',
        [('Henry2018','KRT5 KRT14 TP63','Basal epithelial identity in Results')], '与SCGB1A1 club样方向分开；通用基底程序不唯一指向前列腺。')
    add('prostate_basal_club','Epi_prostate_club_like','前列腺club样候选','prostate','SCGB1A1',
        [('Henry2018','SCGB1A1','OE1/club epithelial population')], '当前仅一个直接核实锚点；必须结合完整DEG和参考映射，不能凭单marker终判，也不能自动当肺细胞。')
    add('urothelial','Epi_urothelial_differentiated_like','分化尿路上皮方向','urothelium','UPK1A UPK1B;UPK2 UPK3A KRT20',
        [('Yu2019','UPK1A UPK1B UPK2 UPK3A KRT20 KRT13','Human epithelial marker table and umbrella/intermediate Results')], 'UPK-rich方向偏伞/中间层，不覆盖所有尿路基底细胞；不是膀胱或上尿路的唯一定位。')
    add('kidney_tubule','Epi_kidney_TAL_like','肾髓袢粗升支方向','kidney','UMOD;SLC12A1',
        [('KidneyProtein2014','UMOD SLC12A1','Table 1 and Fig.4 protein localization'),('Lake2023','UMOD SLC12A1','Fig.3d p.589 marker labels (healthy/altered Visium comparison) / TAL segmentation')], '与近端小管分开；损伤时成熟转运程序可降低。')
    add('kidney_collecting_podocyte','Epi_kidney_principal_like','肾集合系统主细胞方向','kidney','AQP2;AQP3 FXYD4',
        [('KidneyProtein2014','AQP2','Fig.4 collecting duct localization'),('Lake2023','AQP2 AQP3 FXYD4','Fig.3d p.589 marker labels (healthy/altered Visium comparison), PC')], '主细胞候选，不覆盖集合系统所有闰细胞；与足细胞分开。')
    add('kidney_collecting_podocyte','Epi_kidney_podocyte_like','肾足细胞方向','kidney','NPHS2;PODXL',
        [('KidneyProtein2014','NPHS2','Highly kidney-enriched genes / glomerular localization'),('Lake2023','NPHS2 PODXL','Fig.3d p.589 marker labels (healthy/altered Visium comparison), POD')], '特殊上皮，不以EPCAM/KRT8阳性设硬门槛；PODXL也需排查内皮，结合肾小球证据。',parent='*')
    add('thyroid_follicular','Epi_thyroid_follicular_like','甲状腺滤泡方向','thyroid','TG;TPO',
        [('ThyroidOrganoid2021','TG TPO TSHR PAX8','Human Fig.2C,D and Fig.4B,D–F')], 'TG/TPO合成机器为关键，PAX8仅辅助且跨器官共享；本参考的人体证据含bulk/类器官，不冒称人单细胞验证。')
    add('mullerian_secretory','Epi_fallopian_secretory_like','输卵管分泌上皮方向','fallopian_tube','PAX8;OVGP1 KRT7',
        [('FallopianAtlas','PAX8 OVGP1 KRT7','Secretory epithelial populations and OVGP1 histological validation')], '缩窄原Müllerian泛标签；周期/绝经影响OVGP1，不覆盖所有卵巢/内膜细胞；MSLN/WT1不作普适核心。')
    add('endometrial','Epi_endometrial_secretory_like','子宫内膜分泌腺体方向','endometrium','PAEP;SCGB2A2',
        [('GarciaAlonso2021','PAEP SCGB2A2','Secretory-phase and spatial glandular Results')], '特指分泌期腺体方向，PAEP阴性不能排除增殖期或基底层内膜上皮。')
    add('squamous_suprabasal','Epi_squamous_suprabasal_like','鳞状上皮上基底分化方向','squamous_organ','KRT4 KRT13 CRNN',
        [('Esophagus2022','KRT4 KRT13','Human suprabasal epithelial Results'),('EsophagusMapping2022','KRT4 CRNN','Mid/terminal differentiation Results')], '不同分化阶段并列审阅，不能仅由该程序判食管、头颈或宫颈来源。')
    add('keratinocyte','Epi_keratinocyte_differentiated_like','角质形成细胞分化方向','skin_or_keratinizing_epithelium','KRT1 KRT10 LOR FLG',
        [('Cheng2018','KRT1 KRT10 LOR FLG','Spinous/granular differentiation Results')], '棘层与颗粒层标记属于不同阶段；不要求每个细胞同时表达，不覆盖基底角质形成细胞。')
    return specs

REUSED = {
    'alveolar_AT2':'Epi_AT2', 'alveolar_AT1':'Epi_AT1', 'club_secretory':'Epi_club',
    'ciliated':'Epi_ciliated','cholangiocyte':'Epi_cholangiocyte_like',
    'prostate_luminal':'Epi_prostate_luminal_like','kidney_tubule':'Epi_renal_tubular_like',
}
NAMES = {'Epi_AT2':'肺AT2方向','Epi_AT1':'肺AT1方向','Epi_club':'气道club方向',
         'Epi_ciliated':'运动纤毛分化方向','Epi_cholangiocyte_like':'胆管方向',
         'Epi_prostate_luminal_like':'前列腺管腔方向','Epi_renal_tubular_like':'肾近端小管方向'}
SPECIAL = {
    ('intestinal_stem','SMOC2'):('excluded_from_identity_core','Burclaff2022','Human Paneth can show higher SMOC2; not an exclusive ISC marker.'),
    ('club_secretory','MSMB'):('not_imported_into_airway_panel','Henry2018','Prostate luminal context; source note itself acknowledges sharing.'),
    ('gastric_chief_parietal','GIF'):('symbol_and_species_review_pending','','Do not transfer mouse chief-cell Gif logic to human; symbol update not silently applied.'),
    ('mullerian_secretory','MSLN'):('not_imported_as_universal_mullerian_marker','','Shared program, existing Epi_MSLN state remains separate.'),
    ('mullerian_secretory','WT1'):('not_imported_as_universal_mullerian_marker','','Requires organ/cell-specific evidence; not a universal secretory identity.'),
    ('islet','TTR'):('not_imported_as_pan_islet_anchor','','Not required across all endocrine subtypes; original source retained for later gene-context review.'),
}

def build(repo: Path, out: Path):
    basepath=repo/'resources/annotation/v7_panel_repair_v1/registry.json'
    srcpath=repo/'methods/annotation_gene_panel_v1.tsv'
    tracepath=repo/'resources/annotation/v7_panel_repair_v1/gene_reference_trace.tsv'
    with srcpath.open() as handle:
        source=list(csv.DictReader(handle,delimiter='\t'))
    tissue=[r for r in source if r['role']=='tissue' and r['compartment']=='epithelial']
    if len({r['label'] for r in tissue})!=len(tissue): raise ValueError('Duplicate source label')
    base=json.loads(basepath.read_text()); merged=copy.deepcopy(base)
    refs={r['id']:r for r in merged['references']}
    for rid,title,url,scope in REFERENCES:
        if rid in refs: raise ValueError('Unexpected reference ID collision')
        item=dict(id=rid,title=title,url=url,verified_scope=scope)
        refs[rid]=item; merged['references'].append(item)
    panels={p['id']:p for p in merged['panels']}
    # Update only the derived renal copy: previous gap was not a literature claim.
    renal=panels['Epi_renal_tubular_like']
    renal.update(positive_genes=words('SLC34A1 LRP2 CUBN'),reference_ids=['KidneyProtein2014','Lake2023'],
                 evidence_level='curated_literature_informed',
                 interpretation='仅近端小管成熟方向；HAVCR1从身份正向集移出，须作损伤轴独立审阅。人体RNA/蛋白及肾图谱补足旧参考缺口。')
    specs=curated_specs(); evidence={}; maps={r['label']:[] for r in tissue}
    for source_label,pid in REUSED.items(): maps[source_label].append(pid)
    for s in specs:
        if s['id'] in panels: raise ValueError('Panel collision')
        p=dict(id=s['id'],axis=s['axis'],parent_any=[s['parent']],core_groups=s['core_groups'],
               positive_genes=s['genes'],contrast_review_panels=[],expected_lower_genes=[],naming_genes=[],
               required_context=s['context'],reference_ids=list(dict.fromkeys(e['reference_id'] for e in s['evidence'])),
               evidence_level='curated_literature_informed',interpretation=s['note'],auto_final_label=False,
               group_logic='Components are separate review dimensions, never a universal AND/OR expression classifier; aggregate values do not establish same-cell coexpression.',
               negative_evidence_policy='Unmeasured/dropout is not negative. Context mismatch requests review, never rejects a candidate.')
        merged['panels'].append(p); panels[p['id']]=p; maps[s['source_label']].append(p['id'])
        evidence[p['id']]=s['evidence']; NAMES[p['id']]=s['display_name']
    evidence['Epi_renal_tubular_like']=[
        dict(reference_id='KidneyProtein2014',genes=['SLC34A1'],locator='Table 1 and Fig.4 proximal luminal membrane'),
        dict(reference_id='Lake2023',genes=['LRP2','CUBN'],locator='Fig.3d p.589 marker labels (healthy/altered Visium comparison), PT')]
    if set(maps)!=set(REUSED)|{s['source_label'] for s in specs} or any(not v for v in maps.values()):
        raise ValueError('Source coverage mismatch; revise curated mapping explicitly')
    with tracepath.open() as handle:
        oldtrace=list(csv.DictReader(handle,delimiter='\t'))
    oldtrace={(r['panel'],r['gene']):r for r in oldtrace if r['direction']=='positive'}
    target_ids=list(dict.fromkeys(pid for r in tissue for pid in maps[r['label']]))
    traces=[]
    for pid in target_ids:
        for gene in panels[pid]['positive_genes']:
            hits=[dict(reference_id=e['reference_id'],url=refs[e['reference_id']]['url'],locator=e['locator']) for e in evidence.get(pid,[]) if gene in e['genes']]
            old=oldtrace.get((pid,gene)) if pid not in evidence else None
            status='primary_gene_context_verified' if hits else ('inherited_prior_trace' if old else 'gene_context_review_pending')
            traces.append(dict(panel=pid,gene=gene,status=status,primary_hits=hits,prior_trace=old,
                               interpretation='Gene occurrence/context only; not specificity, a validated cutoff or an exact author signature.'))
    lookup={(r['panel'],r['gene']):r for r in traces}
    source_gene_trace=[]
    for r in tissue:
        for gene in r['genes'].split(','):
            canonical=merged['aliases'].get(gene,gene)
            active=[pid for pid in maps[r['label']] if canonical in panels[pid]['positive_genes']]
            core=[pid for pid in active if canonical in sum(panels[pid]['core_groups'],[])]
            entry=dict(source_label=r['label'],original_gene=gene,canonical_gene=canonical,target_panels=active,
                       disposition='retained_core_review' if core else ('retained_auxiliary_review' if active else 'not_imported_pending_gene_context_review'),
                       evidence=[lookup[(pid,canonical)] for pid in active])
            if (r['label'],gene) in SPECIAL:
                action,ref,note=SPECIAL[(r['label'],gene)]
                entry.update(disposition=action,reason=note,review_reference=refs[ref] if ref else None)
            source_gene_trace.append(entry)
    merged.update(resource_version='v7_epithelial_extension_v1',
                  base_resource=dict(path='resources/annotation/v7_panel_repair_v1/registry.json',sha256=sha(basepath)))
    merged['policy'].update(origin_assignment='Context compatibility, not forced organ prediction. Collection site and primary organ remain distinct.',
                            benign_assignment='Never inferred from a normal-reference panel. Use independent CNV/genetic/pathology evidence.',
                            source27_claim='Source file contains 26 tissue rows; no invented 27th row.')
    manifest=dict(resource_version=merged['resource_version'],source_rows=len(source),source_tissue_rows=len(tissue),
                  mapped_target_panels=len(target_ids),new_panels=len(specs),reused_panels=len(REUSED),
                  base_panels=len(base['panels']),merged_panels=len(merged['panels']),
                  merged_unique_positive_genes=len({g for p in merged['panels'] for g in p['positive_genes']}),
                  original_gene_entries=len(source_gene_trace),active_target_gene_entries=len(traces),
                  directly_verified_target_gene_entries=sum(r['status']=='primary_gene_context_verified' for r in traces),
                  original_gene_entries_not_imported=sum(not r['target_panels'] for r in source_gene_trace),
                  input_sha256={'source_tsv':sha(srcpath),'base_registry':sha(basepath),'base_gene_trace':sha(tracepath)},
                  altered_base_panels=['Epi_renal_tubular_like'],final_cell_label_changes=0,
                  author_exact_signature_claim=False,cell_level_annotation_executed=False)
    outputs={'registry_merged.json':merged,'source_tissue_rows.json':tissue,
             'source_to_panels.json':maps,'active_gene_evidence.json':traces,
             'original_gene_disposition.json':source_gene_trace}
    out.mkdir(parents=True,exist_ok=False)
    for name,value in outputs.items(): (out/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')
    rows=['# V7 组织上皮候选模块完整表','',
          f'源文件实际 {len(tissue)} 行组织上皮组合，映射到 {len(target_ids)} 个审阅模块：新增 {len(specs)}、复用 {len(REUSED)}。合并资源共 {len(merged["panels"])} 模块。这里的方向标签不是已经完成的细胞注释。','',
          '“核心组件”用于分别检查表达和DEG证据，斜线表示同组件候选，分号分隔组件；不构成统一布尔分类门槛。纤毛、鳞状、内分泌等跨器官共享。原基因的未纳入项仍完整保存在逐基因处置表。','',
          '| 原组合 | 当前候选模块 | 核心组件 | 全部当前正向基因 | 修订与适用范围 | 文献 |',
          '| --- | --- | --- | --- | --- | --- |']
    for r in tissue:
        for pid in maps[r['label']]:
            p=panels[pid]; links='、'.join(f'[{rid}]({refs[rid]["url"]})' for rid in p['reference_ids'])
            rows.append('| '+' | '.join([r['label'],NAMES[pid]+' (`'+pid+'`)','; '.join('/'.join(c) for c in p['core_groups']),', '.join(p['positive_genes']),p['interpretation'],links])+' |')
    rows+=['','原始组合保留于 [源行快照](source_tissue_rows.json)；[逐基因处置表](ORIGINAL_GENE_REVIEW.md) 覆盖每个原始基因；[当前模块逐基因证据](active_gene_evidence.json) 区分本次直接核实与继承证据。','',
           '## 文献范围','']
    for rid in dict.fromkeys(rid for pid in target_ids for rid in panels[pid]['reference_ids']):
        ref=refs[rid];rows.append(f'- [{rid} — {ref["title"]}]({ref["url"]})：{ref["verified_scope"]}')
    (out/'ALL_PANELS.md').write_text('\n'.join(rows)+'\n')
    review=['# 原始组织上皮基因逐项处置','',
            '保留原始拼写与语境。未纳入表示本次没有足够的该基因—该类型依据或有语义问题，不表示这个基因生物学上无用。继承记录不升级成“本次逐基因已核实”。','',
            '| 原组合 | 原基因 | 当前符号 | 当前模块 | 处置 | 依据或待核实原因 |','| --- | --- | --- | --- | --- | --- |']
    for r in source_gene_trace:
        notes=[]
        for e in r['evidence']:
            if e['primary_hits']:
                notes.extend(f'[{h["reference_id"]}]({h["url"]}) {h["locator"]}' for h in e['primary_hits'])
            elif e['prior_trace']:
                old=e['prior_trace']; notes.append('继承V7 trace：'+old['status']+'；'+old['source_locator'])
        if r.get('reason'):notes.append(r['reason'])
        if r.get('review_reference'):
            ref=r['review_reference'];notes.append(f'[{ref["id"]}]({ref["url"]})')
        review.append('| '+' | '.join([r['source_label'],r['original_gene'],r['canonical_gene'],', '.join(r['target_panels']) or '—',r['disposition'],'；'.join(dict.fromkeys(notes)) or '尚无本次核实的基因—类型直接定位；保留原候选'])+' |')
    (out/'ORIGINAL_GENE_REVIEW.md').write_text('\n'.join(review)+'\n')
    manifest['output_sha256']={p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()}
    (out/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    return manifest

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[2])
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args();print(json.dumps(build(args.repo,args.out),ensure_ascii=False,indent=2))
