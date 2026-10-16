# 既有九来源的候选外部输入机会与稳定ID别名核对

## 全范围、层级与限制

复用已发布免疫档案9来源/400文库，对全部原父记录核基因及诊断细胞池的机会，没有表达评分、效应或新下载。预检保留108来源—父记录行、408原组织/时点/治疗/捕获分层行。已有CD8/Treg细身份是候选，不把粗池比较当NSCLC精确作者亚群复现，也不修改270或原172身份。

只按已审患者块计数，未审身份另列文库替代单位；CAR产品、离体刺激、不同组织/时点保持区别。归档临床资格可能早于后续作者表sidecar，不由旧false推断现在没有临床资料。原全部18临床/细胞控制并非能在外部直接重建，未静默删变量拟合。

## 初步字面匹配及失败修复

生产2026-10-16T12:12:03.690397+00:00 /20:12:03.690397 UTC+8启动，12:13:15.650895 UTC完成；源码SHA `b23de7feb89eb8519672b4524760170943ffde584722ba463b7a7f5372151fc0`，终态SHA `a04a053972ac67608a79b00be01527046dde9c49e91d809a974cd9246da514b0`。

审计v1于12:15:43.200344 UTC启动，CSV将True/False自动转成布尔后与原字符串比较失败；仅解析失败，未重启原生产。保留旧日志和源码SHA `f1f5da37aaaea36bf5478e7d4a793f0fd283da3f26a330f2959fe00d45f8534a`。审计v2显式固定原元数据字符串类型，12:17:11.374047 UTC启动、12:17:32.237011 UTC通过；源码SHA `745dca6d45ddd638223080cd9c4bdfa9e6dd3205c7effa00cda6f97ab99f84cf`，审计终态SHA `090d4f93835323ef9af220050178e7daca9552d3f578ea54cca5fceb8ac5a157`。

## 权威别名修正v2与独立审计v3

初步版本把整个别名字符串当单一名称，只是精确名称入口；尚需原生别名解析。已读项目既有`public_gene_identity_contract.py`及`check_HNSCC_frozen_gene_axis_v2.py`：稳定Ensembl ID收集的别名由pipe分隔，排除歧义来源特征，仅唯一物理列可用，冻结基因碰撞全部标未解析，不按权重/效果挑列。

新生产v2于12:23:18.270525 UTC /20:23:18.270525 UTC+8启动，12:23:39.977642 UTC完成；源码SHA `8ead016508d2bb78e4470fd06b03fd0631c06f6fc99a37e4bedb8349ae9b6aec`，终态SHA `d2760e7f0e2240ef4eb6cafff4d20406a50d445e348eab082850eb303aae554c`。新审计v3于12:23:39.768380 UTC启动，12:24:01.055297 UTC验收全部400文库元数据/细胞池、25659基因掩码、108覆盖与408分层；源码SHA `c1874bd86ac8465ab6df4226e254bab73bd75014a23076985480ef11e5ed6776`，终态SHA `41873b55442f8319cda7d9b3d57cddd94731a72ea2eff803c33fb04a45770910`。

36项发送程序的95%测量门槛因准确名称对齐恢复，细胞数量、固定权重、95%门槛和原科学效果均未改变；初步未达门槛不能被解释为RNA不存在。所有具体分子及来源分母仅服务器/授权对话。最终逻辑根`${RUN_ROOT}/existing_immune_cohort_candidate_extension_preflight_20261016_v2`，v1保留。

```bash
python scripts/preflight_existing_immune_cohort_candidate_extension_v2.py
python scripts/audit_existing_immune_candidate_extension_inputs_v3.py
```

来源背景另参考[GSE169246原记录](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169246)及[GSE176021原记录](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176021)，血液/多时点和离体刺激不能混同肿瘤独立病例。另两GEO直开未成功，不据搜索摘要完成新临床映射。

下一[完整原细分单位测量入口](20261016_existing_immune_fine_native_measurement_scope_v1.md)，仍需效果前的精确端点匹配或明确定义的谱系扩展；本预检不是外部复现成功。
