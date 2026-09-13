# 可执行任务
## Biological Freeze Sprint task board (authoritative; 9 lines)

| Order | Status | Scientific gate |
| --- | --- | --- |
| P0-1 | DONE/FROZEN | Myeloid/DC acceptance rules |
| P0-2 | DONE/FROZEN | Third-GSE validation and `FROZEN_V8_MYELOID_TAXONOMY` |
| P0-3 | DONE/FROZEN | T/NK refinement using the same template |
| P0-4 | DONE (targeted) | Genuine scATOMIC rerun; broader classifier use only if it changes a gate |
| P0-5 | DONE/FROZEN | Fibroblast refinement and FROZEN_V8_FIBROBLAST_TAXONOMY |
| P1-1 | NEXT | Lightweight Endothelial and B/Plasma validation |
| P1-2 | TODO | 10–15 GSE held-out benchmark after rules freeze |
| P1-3/4 | TODO | Freeze identity, verify Copykat `hg20`, then rebuild V8 malignancy once |
| P1-5/FINAL | TODO | One robustness pass → `FROZEN_V8_BIOLOGICAL_REFERENCE` → V8 CoVarNet |


- `V8-003`（done_supporting_only）：`20260913_084700_mydc_full_v1` 已完成 review-only 处理 882,415 scoped
  cells / 1,321 samples / 54 GSE。 终态：1,075 success、246 explicit skipped、0 failed，5,048 cluster rows。完成标准：所有样本 terminal、0 failed、skips 有理由、合并
  cluster/recurrence 表和机器 summary 存在；完成后仍须人工/跨患者审查才能决定任何 V8 写回。
## 2026-09-06 新优先项

- `V8-002`（done_frozen_taxonomy）：Myeloid/DC review-only V3 已完成 2 GSE、4 samples 的
  fail-closed pilot；先冻结 cluster 接受/退回规则，再运行第三 GSE，重点验证 cDC2、pDC、
  SPP1/TREM2 与 unresolved parent fallback。完成 held-out benchmark 和预注册门控前不得写回 V7。

- `REL-001`（done）：冻结并登记 V7 identity/CNV release。
- `REL-002`（done）：登记通讯—空间—临床 evidence-chain release。
- `REL-003`（todo）：在 V7 上重跑依赖身份的 CoVarNet/LIANA，并建立 V5-vs-V7
  下游差异审计；完成前不得把当前通讯 release 称为 V7-derived。
- `REL-004`（done）：登记 V8 full-metadata software candidate 的 rerun2、策略门控和
  外部产物 hash；root audit 进一步记录
  `PASS_DEVELOPMENT_SOFTWARE_CANDIDATE_NOT_FROZEN`、V7 默认和外部记录的 167 项 pytest
  通过；该索引不构成 biological freeze。
- `V8-001`（blocked）：独立 reference 当前为
  `BLOCKED_NO_USABLE_PINNED_CLIENT_AFTER_BINARY_RETRY`。隔离安装和唯一官方 PyPI
  binary-only retry 的失败日志已保存；获得固定 client wheel/cache 或稳定官方连接后，先执行
  metadata-only source/donor/GEO、ontology、license 审查，再进入独立参考与 classifier。
  在此之前不得把软件候选晋升为项目参考或改写 downstream。
- `SUP-001`（todo）：为最终主图建立 `figure-panel-source.tsv`，记录每个 panel
  的 release、源表、筛选、统计量、分母和脚本。
- `SEC-002`（todo）：首次上传新增 releases 前，人工复核 private remote、diff、
  路径脱敏、文件大小和数据治理许可。

| ID | Priority | Status | Task | Dependency | Completion criterion | Related files |
| --- | --- | --- | --- | --- | --- | --- |
| REL-004 | P0 | done | 登记 V8 full-metadata software candidate、root gate 与外部 checksum index | rerun2 summary、statistics、policy、attestation 与 root-audit evidence 可读取 | `PASS_DEVELOPMENT_SOFTWARE_CANDIDATE_NOT_FROZEN`、V7 默认和外部记录 167 pytest 通过均已索引；未形成 biological freeze | `releases/annotation_v8_development_20260912/ROOT_AUDIT_ADDENDUM_20260912.md` |
| V8-001 | P0 | blocked | 为 V8 候选建立独立生物学验证与晋升决策 | 获得可用的固定 Census client wheel/cache 或稳定官方连接；独立真值/参考、实际 classifier 或正交证据、下游影响审计 | 完成 metadata-only source/donor/GEO、ontology、license 审查；新版本化验证 run 明确区分准确率、软件合同和生物学结论；任何晋升均有审查证据 | `runs/20260912_v8_myeloid_dc_review_pilot.md` |
| REC-001 | P0 | done | 核验 GitHub private 远程并完成普通 push | 服务器 GitHub SSH identity | 已通过显式 identity 普通 push，并核对远端 `main` 与本地提交一致 | `logs/codex/2026/2026-09-12_v8-execution-and-cross-gse-pilot.md` |
| REC-002 | P0 | done | 审计并同步统一 CoVarNet 运行的处理、参数和日志索引 | 权威运行目录可读 | 正式 run record、状态、方法、清单、Codex log 和 incident 均通过仓库验证 | `runs/20260725_120105_unified_pancancer_covarnet.md` |
| MON-001 | P0 | in_progress | 监控三个 active extension 并维护可审计快照 | 三个运行目录和日志可读 | 每次状态变化更新 run record/Codex log；不得以进程存在代替产物验证 | `STATUS.md` |
| IMM-001 | P0 | blocked | 完成四队列免疫治疗整合与最终验收 | 统一 11 个 inferCNV reference-insufficient skip 的验收语义 | 保留 19/19、71,398 spots 严格 merge；新 Stage 13 attempt 生成并通过 `FINAL_DELIVERABLE_AUDIT.json` | `runs/20260803_immunotherapy_4cohort_analysis.md` |
| CM-001 | P0 | blocked | 完成恶性识别、CoVarNet、LIANA 与 NicheNet active v2 | `USER002` 组配额余量至少 200 GiB；从零基索引 580 建立新 append-only recovery | Stage 01 全 terminal，active v2 下游重算，`FINAL_AUDIT.json` 通过 | `runs/20260731_malignancy_communication.md` |
| CM-002 | P0 | todo | 核对 CopyKAT 实际 `genome=hg20` 参数是否符合项目意图 | 原执行者或方法依据 | 参数得到书面确认，或以新 attempt 使用纠正值重跑受影响阶段 | `runs/20260731_malignancy_communication.md` |
| DL-001 | P0 | done | 完成 24-GSE 的 141-task recovery retry | 网络、磁盘容量 | final retry summary 存在；127 verified、7 skipped_verified、6 unavailable_upstream、1 failed 均为 terminal | `runs/20260805_external_geo_download.md` |
| DL-003 | P0 | todo | 以新版本化 attempt 恢复 GSE201347 大型 RDS partial | 网络、磁盘容量、retry v3 失败证据 | 文件达到 19,071,156,087 bytes 并通过 size/gzip 完整性检查；旧 partial 和失败记录保留 | `runs/20260805_external_geo_download.md` |
| DL-002 | P1 | todo | 下载后人工复核 core/sidecar 与数据集 guardrail | DL-001 | 每个 GSE 的 GEX、TCR/ATAC/VDJ/spatial、human/xenograft 和重复测量范围明确 | `runs/20260805_external_geo_download.md` |
| GOV-001 | P0 | todo | 对齐源 `stage_status.tsv` 与最终产物状态 | 最终审计文件 | manifest 不再把已完成阶段标记为 `pending`，且保留历史失败/恢复链 | `logs/incidents/2026-07-31_unified-run-record-drift.md` |
| GOV-002 | P0 | blocked | 补齐 `GSE278694` 从最终整合移除的科学理由 | 原执行者决定或原始证据 | ADR/run record 写明理由、影响和是否永久排除 | `logs/incidents/2026-07-31_unified-run-record-drift.md` |
| ANN-001 | P0 | done | 冻结 Phase 08 annotation 输入与注释策略 | P3 integrated model | run record 写明对象、环境、marker/reference 和完成标准 | `runs/20260725_120105_unified_pancancer_covarnet.md` |
| ANN-002 | P0 | done | 执行并验证主类群和亚型注释 | ANN-001 | 标签列、证据表、UMAP、计数表和最终审计齐全 | `methods/project-v3-unified-workflow.md` |
| ANN-003 | P1 | todo | 人工复核低置信度和争议注释 | ANN-002 | 56 low-confidence、10 unresolved、5 ambiguous clusters 有复核结论 | `STATUS.md` |
| EXT-001 | P1 | todo | 决定 `GSE162498` 的 extreme-scale 分支处理 | 资源与分块策略复核 | 形成 ADR 或正式 run，并明确是否进入主整合 | `ROADMAP.md` |
| CNV-001 | P1 | blocked | 核验并完成恶性/CNV 分支及其证据 | P4 annotations；CM-001 组配额阻断解除 | 由 CM-001 的 active v2 run record 和最终审计闭环 | `runs/20260731_malignancy_communication.md` |
| VAL-001 | P1 | done | 对三个免疫治疗队列执行固定 9 模块投影 | 冻结 mapper 和 W | 三个队列 `FINAL_AUDIT.json` 均通过 | `runs/20260725_120105_unified_pancancer_covarnet.md` |
| VAL-002 | P1 | todo | 选择含可靠临床 endpoint 的独立验证队列 | 数据许可与正式 metadata | 预先冻结 endpoint、配对规则和统计模型并完成验证 | `STATUS.md` |
| INV-001 | P1 | todo | 补全最终 31 个 discovery GSE 的逐项记录 | 冻结 final concat manifest | 每行含 dataset ID、角色、状态和权威证据 | `inventories/datasets.tsv` |
| REP-001 | P2 | in_progress | 持续为每次有意义任务写入 Codex log | 无 | 新任务均有日志；正式运行另有 run record | `AGENTS.md`, `templates/` |
| SEC-001 | P1 | todo | push 前进行人工敏感信息审查 | REC-001 | 验证脚本通过且人工确认 diff 不含敏感路径/数据 | `docs/DATA_POLICY.md` |

状态使用：`todo`、`in_progress`、`blocked`、`done`。已完成任务保留原行并链接证据。

## 2026-08-16 07:00 执行门槛复核

- `IMM-001`：允许 skip 的书面验收决定仍缺失，保持 blocked。
- `CM-001`：group-quota 余量约 101.588 GiB，仍未达到表内 200 GiB 依赖，保持 blocked；不得仅因 `df` 可用空间充足而启动。
- `DL-003`：未发现新 recovery；现有 partial 保留 failed，任务仍为 todo。
