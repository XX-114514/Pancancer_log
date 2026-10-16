# GSE305078官方样本时点和分期sidecar

## 问题与输入

[外部病例机会](20261016_external_fine_pair_case_opportunities_v1.md)发现本来源已有原始表达及32个文库，但原缓存obs只有_index/batch，缺临床字段。GEO官方索引说明全部治疗前采样，原提交单细胞子集32，不能采用原论文更广多组学研究171的分母。来源为[官方GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305078)、[原论文PubMed](https://pubmed.ncbi.nlm.nih.gov/40929245/)。

本轮只从NCBI公开FTP取得小型family SOFT元数据，压缩上限2MiB，没有下载表达/测序或先验。原Zarr、患者身份和临床结局保持。新的sidecar只支持精确采样时点/作者分期和组织，不自动改临床或独立验证资格。

## 原始取得与零匹配诊断

2026-10-16T13:05:06.515283+00:00 /21:05:06.515283 UTC+8启动v1，13:05:10.793635 UTC完成；源码SHA `87468c62464dcbdfd2fefa9ac315f69c8b2af9367571667f7ef448c44bdfd6e4`。官方SOFT SHA `e24c261f5c44157f376c7c1259b779bd77dfe5ef3bbe7eb78dff2ee3c3182d49`，原解析JSON SHA `15da40534708bbf41231ab69f50064960ced6da3fa103c3af62463432fc8e005`。明确overall_design支持全部治疗前，但整字符串与官方title/accession直接相等为0，正确保留映射pending，没有猜时点或疗效。

只读检查原缓存batch格式，所有文库原名均含官方GSM及样本标题的独立字面token；不能把未匹配误判为没有可用来源。未输出实际患者/样本行到公共仓库。

## 双标识修正与独立原文重解析

v2复用已SHA核对的原SOFT，不再次联网下载；要求每原文库含唯一完整GSM token和唯一完整官方title token、二者指向同一记录，全部32双射。Stage/tissue逐条来自原GSM characteristics，不按E/A前缀猜分期；另用不同SOFT分块/正则从原压缩文本独立重读全部title/stage/tissue，和首解析核对。

13:07:01.487748 UTC /21:07:01.487748 UTC+8启动，13:07:01.853615 UTC通过。源码SHA `72675310c6105dcd9c2fb14dbc5c2a755e80a435956110a357dfc261b4af0afd`，终态SHA `f419b0ecb032d6fee9ce2366b237ab3e30b26cdd8d0791152bc44c1e33588058`，sidecar SHA `d443ac6d69af7d7831a499b18b6c23694afe8c205ff87bb2dd0a0f052bd5a3a2`。

原GEO characteristics实际只有tissue/stage，无pCR或反应字段。全治疗前phase只来源于明确all-samples设计声明；不能由此认定治疗组、疗效或真实患者ID已核实。patient_identity_newly_verified/clinical_outcome_adopted均false，未覆盖原authority。

## 全部原细对分期机会及审计

13:08:04 UTC /21:08:04 UTC+8以原全部234家族—细对分别按官方阶段统计，不读表达或效果，完整468分期记录保留。原单位>=20双端及未知患者限制保持；纯元数据集合计算证明SHA `f581811219afeaf5bade181ef04ccc7eb8df6f9a0d0f53e4c1cbdb2ab21516e4`。

独立pivot审计13:08:41.774780 UTC /21:08:41.774780 UTC+8启动，13:08:42.609651 UTC通过全部468记录与32个sidecar标志；源码SHA `f3fd293636481e97ac2732fbc31e8cf9adf5db8238312461ce23897e8fe2271f`，审计SHA `1808df3fdfb91d41a6b606dea97c60239bdebde835c8e553edd2b92b4bfa44b4`。具体阶段/细对数只当前获授权对话及私有目录，公共仓库不刊样本或科学汇总表。

```bash
python scripts/curate_GSE305078_official_sampling_sidecar_v1.py
python scripts/curate_GSE305078_official_sampling_sidecar_v2.py
python scripts/audit_GSE305078_stage_fine_opportunities_v1.py
```

逻辑根`${RUN_ROOT}/GSE305078_official_sampling_sidecar_20261016_v1`及`v2`。下一在[原生测量](20261016_existing_immune_fine_native_signatures_v1.md)独立验收后，效果前固定早/晚阶段、细节点、样本替代单位和可识别控制；未运行新外部模型/P/BH或正式临床验证。
