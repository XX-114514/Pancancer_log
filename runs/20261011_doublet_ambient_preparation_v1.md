# 双细胞、环境RNA与整合输入准备 v1

状态：方法及环境预检。输入为`${COLLECTION_RUN_ROOT}/objects/qc_<dataset>.zarr`，4队列993,312细胞149文库；新输出写入独立`${PREPROCESS_RUN_ROOT}`。原始计数与旧QC对象保持原样。

输入：CXG_5b9b89d2_2e3e_4ad1_8965_666f64540fea（GSE241783）、GSE286827、GSE292046、GSE294109。GSE286827应用已验收patientfix_v2覆盖；稳定ID不变。

双细胞方法及污染校正参数在环境验证后冻结。每文库需保存scores、flags、移除原因、污染率、校正前后计数及marker检查。缺空液滴或单一谱系限制明确登记。

验收：原始counts与上游完全一致，corrected非负有限且不凭空增加表达；新对象ID唯一且为上游ID子集，所有新过滤可追溯；逐文库成功状态与来源AC状态分开。整合准备保留原始基因覆盖，技术／组织／富集分支，禁止缺失基因被解释为真实零表达。

完整执行、失败和证据见[任务日志](../logs/codex/2026/2026-10-11_preprocess-doublet-ambient.md)。本次尚未进行跨研究整合。

## 21:46 UTC+8 方法冻结与范围更新

初始两个完整文库试运行均经独立全细胞验收。正式运行另立`${PREPROCESS_RUN_ROOT}`=prepare_v2；范围9分支792文库，输入3338473细胞（4QC+5预QC），不能作为最终通过数。参数及限制见[方法](../methods/doublet-ambient-preintegration-v1.md)。运行命令：

```bash
python scripts/run_ambient_preparation.py --run-root "${PREPROCESS_RUN_ROOT}" --collection "${COLLECTION_RUN_ROOT}" --rscript "${SCRNA_RSCRIPT}" --code-git-sha "${CODE_GIT_SHA}" --parallel-datasets 4 --library-workers 4
```

实际启动、终止及独立验收结果将追加；现阶段不声称全量完成。

## 2026-10-12T02:12:48.762758+00:00 恢复记录

prepare_v2因组quota errno122终止，保留部分输出、无全量验收。data4组上限现2048G/已用约1495G；恢复输出为另卷私有`${PREPROCESS_RUN_ROOT}`=prepare_v3，输入原位只读。预算782893303100 bytes+500GiB，小于输出卷约1.598TB实际余量。按原9分支792文库恢复；新增必需参数`--storage-plan "${STORAGE_PLAN}"`。不恢复新下载。

## 2026-10-12T10:38:45.707906+08:00 prepare_v3实际运行

启动02:15:17 UTC/10:15:17 UTC+8；代码cdcea29，4来源×4文库，500GiB预留，存储计划已核实。10:27快照44文库/0失败，275653保留细胞、额外9608；完整独立count/ID/标准化审计将在全队列终态执行。final source AC仍0，未开始全局整合。

```bash
python scripts/run_ambient_preparation.py --run-root "${PREPROCESS_RUN_ROOT}" --collection "${COLLECTION_RUN_ROOT}" --rscript "${SCRNA_RSCRIPT}" --code-git-sha cdcea29b45e989cb8d03235b9a92521611b6dd17 --parallel-datasets 4 --library-workers 4 --storage-plan "${STORAGE_PLAN}"
```

[可读快照](../inventories/runtime/2026-10-12/AMBIENT_PROGRESS_v3_checkpoint2.md)不等于实时状态或最终验收。

## 2026-10-12T10:54:17+08:00 第一批完整来源验收／审计I/O优化

- 实测10:53:50（UTC02:53:50）182/792文库，0失败，944249保留细胞、额外38637个预测双细胞，274673严格技术候选。完整运行仍活跃。
- 3个完整来源独立验收：CAR-T 396593、GSE292046 52525、GSE294109 130243，合计579361细胞/92文库；原始审计02:46:16 UTC和缓存审计02:54:17 UTC均PASS；所有来源/文库科学检查、计数、候选124513及基因覆盖诊断完全一致。报告及SHA匹配证据见runtime/2026-10-12。
- 审计缓存只改读取方式：每次最多一源36GiB稀疏参考，较大来源保持原backed读取，校验数学/门槛不变；原始报告保留。生产代码、双细胞和DecontX参数未改。全量阶段审计脚本升级尚待显式版本登记。
- 数量为技术处理与验收，不等于来源别名排重已完成；final AC仍0，未启动跨来源整合、新下载仍暂停。CAR-T产品及富集snRNA用途限制保留。

## 2026-10-12T11:05:00+08:00 后续整合合同与native初始QC审计冻结

- 3c1cd3fd845b53d2ee21817d7733cad59448bdb1已普通push并独立HTTPS ls-remote核验一致。签名转发首次因Cursor既有remote forwarding端口冲突失败，复用已有SSH master并换唯一的本地签名socket后成功；没有换科学代理端口、复制密钥或重启计算。
- 审计v2在02:59:05 UTC明确登记并应用：保留原v1源文件、最初execution_start不改，最终审计尚未启动；代码SHA51df1b…，与原92文库全量科学检查一致。
- 新整合输入合同分普通组织、富集T/snRNA和CAR-T产品；原始计数/估计层、基因可测覆盖、患者块及作者前处理历史分别核验。CAR-T作者双阳性状态移除3591/30427(11.80%)、高于CD4(4.08%)/CD8(5.06%)；保留全部原始映射，状态关联需敏感性审查。聚合证据无患者行。
- 加强独立审计：对5个native分支的全部原始计数另算floor/mt/作者PASS/作者doublet及全行过滤账本，检查每个非空原始文库均被覆盖。生产不改；缓存native账本和原始文库索引，避免每文库重解压整源账本。新代码SHA37eaef…已语法验证，真实native独立执行仍待下一步，不能提前称验收完成。

## 2026-10-12T11:07:39+08:00 原生初始QC全量复核通过

- 取证03:06:50 UTC/11:06:50 UTC+8：LBCL598289→597962、儿童HGG401253→400899、SPECTRUM927205→880042；合计1926747原始细胞→1878903初始QC保留，全行独立账本检查全部PASS。这不是后续污染/双细胞全量终态。
- 加强审计v3于03:07:39 UTC登记应用，代码37eaef…对应公开20505b1；旧v1/v2脚本和最初运行manifest保留，模型生产参数不变。新增原生全行过滤复算、非空原始文库覆盖门槛，账本每源只读一次。正在另跑已完成来源的全量验收新检查点。
- 全行技术元数据复核4原生CELLxGENE源：SPECTRUM cell，LBCL/HGG/乳腺nucleus，乳腺v2/v3均在各原始文库内唯一；作者是否已污染处理仍unknown，未从整数raw/X推断“未经校正”。
- HNS部分完成40文库中25个污染估计待复核：21个marker损失、8个高UMI损失（可重叠），0不收敛。继续保留raw并阻止采用这些估计，不能把处理成功解释为污染清除通过。
- 时间补正：上一人工节点标题11:05:00不是实际取证时间；20505b1提交实际11:03:42，v2应用10:59:05。实际节点以JSON时间戳/commit时间为准，旧标题作为补录背景保留。
