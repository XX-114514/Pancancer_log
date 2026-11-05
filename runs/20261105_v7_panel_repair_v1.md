# V7全panel修复：资源、算法与回放验证

状态：规则与可执行证据流程已交付；全V7重新注释未完成。任务于2026-11-05进行，初始工作精确秒级时间未记录；服务器合成验证实际时间为06:12:28–06:12:47 UTC / 14:12:28–14:12:47 UTC+8，终端信息于本轮随后取证。

## 输入、输出与方法

旧17谱系162项配置及上一轮2052条panel-gene证据；最新T6的1212条实测marker及top10 RNA DEG；Lodi/Qian/Chu/Ma作者方法、代码/表格及补充原始研究。完整输入SHA见[回放manifest](../inventories/analysis/2026-11-05/v7_panel_repair_v1/replay_v2/manifest.json)。

输出[112项registry](../resources/annotation/v7_panel_repair_v1/README.md)，覆盖409不同基因；旧162逐项迁移；21项来源范围与634条逐基因依据。完整Excel与Markdown显式表同目录。算法取消跨身份/状态argmax和margin，按来源测量、组件、RNA DEG、同细胞支持、命名条件分别报告。[方法](../methods/V7-panel-repair-literature-workflow-20261105.md)明确作者逻辑与项目适配。

## 执行

环境：服务器R环境`scrna_r`，Seurat5.3.1、Matrix1.6.5、jsonlite2.0.0；本机Python3/numpy/pandas执行证据引擎。官方UCell未安装在此次R环境，未使用自制替代。Harmony已安装但本轮合成测试未启用。R精确版本未在该预检中记录，不推测。

```bash
python scripts/v7_panel_repair_v1/test_evidence_engine.py
python scripts/v7_panel_repair_v1/run_smoke.py \
  --rscript Rscript \
  --registry resources/annotation/v7_panel_repair_v1/registry.json \
  --out "${RUN_ROOT}/synthetic_test"
python scripts/v7_panel_repair_v1/replay_aggregates.py \
  --repo "${LOG_REPO}" --out "${RUN_ROOT}/aggregate_replay"
```

合成输入360细胞×2100feature，2来源；seed191，500HVG、20PCA、10图PC、分辨率0.2/0.4、模型基因来源覆盖率1、marker_resolution0.4。先按可审阅清单排除4个测试细胞，再重建356细胞，得到3个图聚类群。这是测试结果，不是预设生物亚型数。

实际退出码fixture/workflow/verify均0。独立验算1227条表达记录、6条Treg同细胞记录和来源IL1RL1不可测分母，最大CP10K绝对误差4.6896e-13。Python15项反例测试通过。旧40群回放4480条、新T6 12群回放1344条；全部`final_label_assigned=false`。工作簿4张表已渲染查看，无公式错误；没有在Excel桌面应用中测试交互。

代码和输入资源SHA、各阶段终端见[服务器测试](../inventories/analysis/2026-11-05/v7_panel_repair_v1/server_smoke_v2.json)与[最终核验](../inventories/analysis/2026-11-05/v7_panel_repair_v1/verification.json)。旧40与T6范围相互重叠，不加总细胞数。旧式log2fc未当Seurat avg_log2FC复用。

## 失败、修正与限制

部分原始网页遭遇404/recaptcha，改用可核实作者PDF、正文或原始研究摘要；不把摘要视为完整分类代码。早期测试字段把“每个组件至少一基因可测”命名过宽，现拆成该字段与全部核心测量完整字段，加入反例。DEG也采用同一版本化alias并拒绝碰撞。

首次独立工作副本快进因受限网络无法补取partial-clone对象失败，随后按已有授权联网快进成功，保留并发TME提交。初轮异步R测试有完成文件但未保存子进程退出码，因此本run以第二轮受监督测试为终端依据。

当前仍有368条模块—基因记录未逐基因核实原始材料；3个器官/神经相关模块保留reference_gap。完整新panel在全V7真实矩阵上尚未全面执行；未做新的CellTypist/scATOMIC、UCell或全谱系来源留出验证。资源与程序测试通过不等于发表级细胞标签完成。

下一步按T/NK/ILC联合边界开始提取完整409基因测量证据，先复核T6 helper/Treg和混合RNA，再进入B/ASC及髓系；每次保留大类、分支、状态、QC和命名证据独立记录。
