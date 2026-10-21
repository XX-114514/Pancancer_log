# 独立NSCLC成像来源的状态、临床与模态资格

## 效果前问题

在原2026-10-14免疫图谱CM02 v1 T_whole_TME_primary_tumor仍为图谱重点的背景下，固定RNA状态临床家族已完成，新来源先核实际测量/病例/治疗/终点，不据旧阴性换阈值。检查Chen等NatImmunol2024 DOI10.1038/s41590-024-01792-2的官方补充表与Zenodo11198494，区分已闭合MERFISH四供者分支、配准RNA/蛋白主面板及小型validation12plex文件。

公开来源：[官方数据](https://zenodo.org/records/11198494)、[论文PubMed](https://pubmed.ncbi.nlm.nih.gov/38503922/)、[作者代码](https://github.com/korsunskylab/nsclc_manuscript_2023)。本轮Nature正文转认证且PMC返回challenge，未绕过，不声称新全文读取。研究摘要的pre-immunotherapy描述不替代逐材料阶段审阅，论文和表总患者数不替代实际细胞文件可测分母。

## 原表与独立审阅

优先复用服务器已保存官方SupplementaryTable1，完整SHA f33f58f219adcecfc5355b1ce8b526dfd36d7eca4555a320ba1deb185d122982。只读B病例表与C/D的Patient#单列，以及P通用染色protocol；不读取C/D已报道细胞频率作新效应。严格整数Excel患者键，不猜Pt前缀/编号，重复/缺失显式计数。完整RECIST类别、非空临床字段、正报告时长及二元事件覆盖，仅来源级汇总；PFS/OS起点和逐标本阶段仍待，不把字段可用当模型可用。

独立stdlib ZIP/XML解析及Decimal数值另核病例/模态键、全部RECIST和事件/时长覆盖，与openpyxl生产结果逐项一致。各模态及临床表的数量范围分开，不将作者汇总表与实际发布细胞文件假定同集合。临床表缺失不补零、未知疗效不按编号推断。

## 实际文件与不同面板

完整官方JSON和每文件最多64KiB前缀只核字段。HTTP206本身不是全文件验证；首次配准前缀收据的ContentRange为null（保存头字段大小写读取方式），不由此推断源未发送，也不宣称范围或完整文件验收。后续小RNA前缀用标准HTTPMessage case-insensitive getter确认起始范围，完整读取另要求200/精确长度/官方MD5。所有前缀人类行只RAM、不保存，不据前缀推患者/细胞数。

官方P协议将panel1RNA、panel2protein和panel4HiPlex分开：主面板Opal通道不能挪用另一面板相同染料的靶标；合并CXCL10/CXCL11探针不拆成两个独立基因读数。主面板缺调节端所需靶标时不将PD1/CD8替成调节身份。validation12plex包含FOXP3/CD4/CD8A/PDCD1/TCF7及组织背景基因，但不继承原完整22基因细胞配方或旧LAG3/TNFRSF9/增殖状态身份。

## 小RNA实际完整资格

ACD_stemimmunity.csv及ACD_tumor.csv按官方固定大小/MD5顺序在RAM全读，仅读取资格，不做新状态拟合/空间P/临床关联。全文CSV的十二基因/坐标有限性、负值和非整数性质、作者type、slide/compartment覆盖与重复位置先审；独立scalar CSV对所有行/靶标/坐标及type/分组数值重放。Exact image/x/y使用原CSV字符串，不从取整或浮点合并推同细胞。图像或slide不是已审核患者，两文件共享位置也不自动证明相同细胞。

这些是作者手工选择的区域，不作随机组织样本或全肿瘤总体。坐标单位、图像/病例/区域层级、type的作者判定及原始计数/变换语义另待核，不能看到正值就给正式细胞状态。主成像巨型文件不因元数据存在而下载。

## 资源与隐私

四项合成检查分别核严格临床键/非阳性时长和已知事件、负/分数RNA保留及缺基因不可填零。完整RNA读前VM free+inactive+speculative须≥768MiB峰值预算加3GiB余量，实际峰值RSS≤768MiB；该VM读数不等同OS保证。默认TLS、原请求级代理，一文件一次完整请求，不重试403或认证，不改全局配置。源raw/临床/表达/坐标/类型行仅RAM，无临时个体文件，来源级科学汇总和解释私有。原权威/身份/已闭家族保持，GSE288199 TableS1用户后补，三项之前拒绝动作不重试。

入口：qualify_Chen_coregistered_clinical_scope_v1.py、audit_Chen_clinical_workbook_XML_v1.py及qualify_Chen_validation_RNA_cells_in_RAM_v1.py。完整输入接受不等于获得直接通讯/临床增量/同终点验证。


实际v1完整资格在RSS断言失败，不接受该输入终态。新v2增量解析与无损紧凑CSV键替代整表/整段Unicode副本，同一输入、科学范围及768MiB/3GiB资源门槛保持。原四项检查及另三项跨块/精确键反例接受；实际v2完整读取另记，不能因代码核通过称真实细胞资格已完成。


实际v2完整输入门槛已接受，详run/私有接受；不将两视图共享坐标键当独立病例。新增[作者公开12plex方法](https://bio-protocol.org/exchange/minidetail?id=20737457&type=30)说明病例按低自发荧光选择，区域按特定趋化/TCF7/CD8信号选择，阳性按手工强度阈值并排AF；由此不把目标marker富集当无偏新发现，也不把非负整数自动解释UMI。与原MERFISH病例重叠及源级测量范围差异私有审阅；当前只接受输入/测量入口，无新空间P/临床效应。
