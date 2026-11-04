# 泛癌论文原始注释标签完整表提取

取证时间：2026-11-04 22:52:41 UTC / 2026-11-05 06:52:41 UTC+8；本机任务日期2026-11-05 Australia/Melbourne。早期下载的逐文件时间及SHA保留在本地下载收据，未将本次取证时间伪作下载起点。

## 用户要求与起点

用户要求完整读取其他泛癌论文的具体注释结果并输出表格。在前一轮8项研究方法比较基础上，本次读取作者实际发布的补充表、元数据注释列、层级映射和分类器类别，而非以正文示例补成完整列表。

已读仓库AGENTS、README、STATUS、ROADMAP、TODO、相关注释核查日志及RUNTIME_SYNC。已有大量并发修改和未跟踪文件；保留全部内容，不整体stage。未访问数据服务器，不修改V7或当前项目细胞标签，不训练新模型。

## 来源与实际提取范围

|研究|原始来源|实际范围|本次记录数|
|---|---|---|---:|
|Qian 2020|[Supplementary Table S5](https://static-content.springer.com/esm/art%3A10.1038%2Fs41422-020-0355-0/MediaObjects/41422_2020_355_MOESM17_ESM.xlsx)|68个主表型、5个嵌套毛细血管群，保留原名及原表前5基因|73|
|Chu 2023|[Supplementary Tables 3/5/7/8/9/10/11](https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-023-02371-y/MediaObjects/41591_2023_2371_MOESM3_ESM.xlsx)|32一级计数标签；Treg 7、Tfh 5、增殖T 8；另5条非常规T群DEG描述词表|57|
|Kang 2024|[Supplementary Data 2](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-48310-4/MediaObjects/41467_2024_48310_MOESM5_ESM.xlsx)|作者5组NMF状态表：间质23、上皮/神经35、T/NK12、髓系20、B9；排除外部参考signature对照|99|
|Yang 2024|[GSE233236 metadata](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE233nnn/GSE233236/suppl/GSE233236_metadata.tsv.gz)|Annotation列全部20标签及公开新生成批次计数|20|
|Lodi 2025|[作者仓库Master_files](https://github.com/lambrechtslab/Pan-cancer_Lodietal/tree/main/Master_files)|四层注释列唯一值14/17/72/82，共185；另73行作者层级映射|258|
|Luo 2022|[论文补充数据及Source Data](https://doi.org/10.1038/s41467-022-34395-2)、[作者代码](https://github.com/Xiaxy-XuLab/PanCAF)|34个TME数字群、23个上皮数字群、8个fibro解释映射、3个CAF轨迹状态、6个DC二级数字群、TEC/NEC视图|76|
|Ma DC 2025|[原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC12485386/)、[作者代码](https://github.com/matyhoc/DC_altas)|10主标签、7个cDC2嵌套状态、3个LAMP3来源推断视图|20|
|scATOMIC 2023|[论文](https://doi.org/10.1038/s41467-023-37353-8)、作者公开模型对象classes|35个论文参考大类、13个当前公开模型分支的42条类别|77|

最终合计680条分层记录。初次提取672条，补查Luo作者代码增加6个DC和2个内皮视图后更新为680。每条包含研究、层级、作者父类/容器、编号、原始名称、作者解释、可用原表前5基因、可用细胞数、计数范围、条目性质、来源编号、原始位置和版本备注。680不是互斥细胞类型数，父类与子群、重复视图和表达状态不可直接相加。

## 关键事实与解释边界

- Chu 32个一级标签计数合计308,048；其中增殖父类不与8个子群重复计数。原表Tf17/Th17拼写差异保留并标注。
- Yang公开元数据42,099细胞覆盖20个标签；这些计数不代表论文511,847细胞全图谱。没有从标签名称猜造作者5大群映射。
- Lodi公开元数据680,430行，去除Low quality 35,999及Doublets 32,681后611,750。最细列原始82个唯一值含QC，去QC后80；不能将其直接当论文70个共享亚型。73行映射包含70共享亚型与癌/上皮、成纤维、肥大3个粗类。红系等额外类别及Specific_subtype占位值均保留标记。
- Lodi原作者映射中CD8+ CCL5+ Treg放在CD4+ TREG父类，元数据也有该标签。保留原始命名并注明不一致，未擅自纠正；Capillary/capillary大小写差异分列解释。
- Qian的Erythroblast为主表型之一；未强制补造红系发育亚型。C2_SFTPC描述及CHIL3L1/CHI3L1拼写差异保留原文并备注。
- Luo保留所列主注释表、CAF解释层及代码显式给出的DC合并和TEC/NEC映射；数字群名保持数字群名，不把DEG首基因擅改成作者正式名称。代码还有C24重聚类步骤但没有结果标签表，未根据6色调色板猜群数；代码DC父群ident=13与主表编号语境未核定一致，不自动跨对象对应。
- Ma的LAMP3来源视图属于计算推断，非谱系追踪确证；图注总数存在差异，未强行补齐计数。
- scATOMIC列的是参考与分类器词表，不是一个查询图谱的逐细胞结果；当前公开模型不冒称为2023发表时冻结版本。仅用自编R脚本读取公开模型对象的classes字段，未运行分类预测或上游分析脚本。保留Exausted原拼写。

## 执行与修改

本地工作目录为`${WORK_ROOT}/pancancer_annotation_catalogue_20261105`；交付为`${OUTPUT_ROOT}/pancancer_annotation_catalogue_20261105/泛癌论文完整注释结果.xlsx`，两个工作表：论文与读取范围、完整注释表。原始下载、源数据收据、提取脚本和中间JSON保留本地，不将含细胞条码/临床行的原始元数据或完整补充文件提交公共仓库。公开内容仅技术摘要和来源。

执行自编fetch_sources.py分批获取公开资料；build_catalogue.py用pandas读取XLSX/文本并按作者列提取唯一标签；read_model_classes.R仅load公开模型到独立环境并输出classes。使用预置artifact-tool单一build_workbook.mjs创建XLSX，未使用openpyxl写入。verify_workbook.py用zipfile/XML逐字段读回交付文件。

本仓库仅新增本日志及对应daily引用，未修改既有科学结果或当前首页。所有并发修改仍保留。

## 验证与失败记录

提取断言核对8项研究分层计数、Chu细胞总数、Lodi质控差额及映射范围。最终保存后逐项核对680行×13字段与提取JSON完全一致；8项汇总COUNTIFS缓存分别为73/57/99/20/258/76/20/77。原生筛选表A1:M681，冻结首行和前两列，公式无错误；两工作表及来源登记预览已目视检查。交付XLSX SHA256：`3a96300ba38863a94783debb2d912263d583be39ba309fe6f269061e97a8fb5a`。

初始下载域名误写为static-content.springer-cdn.com导致DNS失败，改用正常官方static-content.springer.com后成功；部分PMC/Cell入口返回验证码或403，未绕过访问控制，改用正常公开出版商补充文件、作者仓库与GEO官方FTP。初版提取脚本属性写法无法表示含点列名，改成列名索引后重新运行全部断言。失败未改变原始标签。

尚未验证：未重算作者聚类、DEG或单细胞准确率；原文名词与元数据版本不一致不以本次表格替作者裁决；没有将全部图谱重新统一成本项目标签本体。

仓库轻量validator实际exit0，0 errors / 1项既有warning，无失效本地相对Markdown链接；git diff --check通过。远端main已比本机主工作区更新，使用隔离检出仅准备本任务两份脱敏日志，避免夹带并发修改；同步收据完成后单独保存。建议commit：`Document complete public pan-cancer annotation label extraction`。下一步建议按谱系将原作者身份标签、marker状态和允许粗类终点分别与本项目cluster证据对照；本次未启动该分析。
