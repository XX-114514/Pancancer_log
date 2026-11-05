# V7 组织上皮注释的文献依据与使用方法

组织上皮面板应识别分化方向，再结合材料和原发部位判断来源。正常组织图谱提供参考表达，不能单独证明待注释细胞正常。当前[完整资源](../resources/annotation/v7_epithelial_extension_v1/README.md)将原26组拆分为32项候选，逐个保留原基因与文献依据；所有组合均为项目整理，未冒称作者完整签名或通用分类器。

## 文献支持的关键修订

| 原组合问题 | 新的处理 | 支持范围 |
| --- | --- | --- |
| 肠干细胞四基因一并当独占证据 | LGR5/ASCL2与SLC12A2/RGMB分组件，OLFM4辅助，SMOC2退出身份核心 | [Burclaff 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9043569/)区分人ISC、TA和Paneth；不能从这些marker直接证明功能性干性 |
| 胃主细胞与壁细胞混在一起 | PGC/PGA3/PGA5/LIPF与ATP4A/ATP4B分别审阅 | [人胃体单细胞与空间研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC10319429/)的对应群体marker |
| 乳腺管腔全按ESR1/PGR判 | LumHR、LumSec、基底/肌上皮分别审阅 | [Kumar 2023](https://escholarship.org/content/qt81t478n0/qt81t478n0.pdf)；KRT5等在部分LumSec状态也可表达 |
| 前列腺基底与club混合 | KRT5/KRT14/TP63与SCGB1A1候选分开 | [Henry 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6411034/)正常前列腺；本轮club仅一个核实锚点，必须补完整参考映射 |
| 肾脏两大混合组合 | 近端小管、TAL、集合系统主细胞、足细胞分别审阅 | [Lake 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10356613/)与[Habuka 2014](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0116125)人体肾脏证据 |
| 胰岛四种激素要求共同满足 | INS/GCG/SST/PPY为不同亚型方向的并集，另看分泌机器 | [Muraro 2016](https://www.hubrecht.eu/app/uploads/2017/11/VanOudenaarden_Key_2016-Muraro_A-single-cell-transcriptome-atlas-of-the-human-pancreas.pdf)对应α、β、δ、PP不同群体 |
| 分泌marker当全部子宫内膜 | PAEP/SCGB2A2仅标分泌腺体方向 | [Garcia-Alonso 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8648563/)区分月经周期和位置 |

## 计算和审阅顺序

1. 确定候选父身份，查看全转录组和来源内分布；父类标签仅为上下文。足细胞、肝细胞或肌上皮等不得因泛上皮marker低表达被机械排除。
2. 在谱系/组织相容范围内复核QC、异源marker同细胞共表达与双细胞旁证；必要时重建邻接图和聚类。不能仅凭某个相邻身份panel高分判污染或删除。
3. 基于RNA表达计算群间及近邻群对照DEG；保留来源特异的测量mask和分母。先看富集是否稳定、同细胞是否支持，再结合参考映射。少量环境消化酶、肺表面活性蛋白或激素不能定义整个群。
4. 独立记录父身份、分化/功能程序、器官相容性、良恶性四层信息。采集部位与原发部位分别保留；未知或不符的器官背景进入人工复核，不剔除表达候选，以免把转移或化生强行解释成本地细胞。
5. 群整体不一致时针对混合来源/状态开展子集分析。文献提供审阅逻辑和参照群体，聚类数量需由当前数据的稳定性与可解释性确定，不能按照32个候选预设32群。

现有 `scripts/v7_panel_repair_v1/evidence_engine.py` 兼容新registry，保留全部候选，无argmax、统一分差阈值或自动终判。`each_core_component_has_measured_gene`仅表示组件有测量记录，**即使表达全为零也可能为True**；`min_component_max_detection`是边际检出摘要，不能当同细胞共表达比例。这里未新增文献未验证的数值身份门槛。

## 与TME注释v2共用

远端并发加入179模块的TME v2后，本轮追加[185模块兼容资源](../resources/annotation/v7_epithelial_extension_v1/TME_V2_COMPATIBILITY.md)：保留41个Lodi模块与其他非上皮内容，将19个项目整理上皮组合替换为本次25项，另更新近端小管。沿用TME v2工作流时，把下方registry路径改为`resources/annotation/v7_epithelial_extension_v1/tme_v2_registry_reviewed.json`；下面的137模块示例是V7对照版本。

## 可复现入口

从仓库根目录使用含pandas/numpy的Python环境。构建只读输入，输出目录必须不存在：

```bash
python scripts/v7_epithelial_extension_v1/build.py --out "$NEW_RESOURCE_DIR"
python scripts/v7_epithelial_extension_v1/test_extension.py
python scripts/v7_panel_repair_v1/evidence_engine.py \
  --registry resources/annotation/v7_epithelial_extension_v1/registry_merged.json \
  --expression "$CLUSTER_GENE_EXPRESSION" \
  --de "$RNA_DEG_TABLE" --contexts "$REVIEW_CONTEXT_JSON" \
  --out "$NEW_EVIDENCE_DIR"
```

只有完整声明的DEG对照才加`--de-complete`，不能用于top10截断表。上下文JSON沿用旧引擎字段 `reviewed_parent_ids` 与 `organ_contexts`；对应组织字符串在registry的`required_context`列可查。JSON中的部位相容信息必须经过审阅，不能将表达预测倒填为已知器官。用户现有生产流程尚未自动切换到新registry，本次未在全V7细胞上运行此命令。

## 证据强度与尚待补齐的范围

基因出现在原始研究某类细胞中，不意味着它特异，也不意味着整组是作者定义的分类签名。新表把本次直接核实的语境与继承证据分开；原121条基因记录中42条尚未导入对应模块，有的因语义问题被移出，有的缺本轮直接定位，原值与处置全部保留。

当前面板服务于原26组的纠正，不是所有人体上皮的穷尽目录。例如肠Paneth/BEST4、肾DCT/闰细胞、尿路基底层、增殖期/基底层内膜等仍需独立参考，不能因不匹配本表而一概标为unresolved或异常。允许大类身份明确、亚型暂不命名；32项候选和137项全资源都不代表固定的真实细胞类数。
