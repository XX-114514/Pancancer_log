# ccRCC当前推进入口

当前主线按用户要求沿用原身份，不以新身份复查为前置。

1. `ccrcc_current_identity_forward_v1.py ${RUN_ROOT}`：从已经保存的原身份快照计算完整12组成/108 LR并准备CellChat输入。
2. `audit_ccrcc_current_identity_forward_v1.py ${RUN_ROOT}`：数值核对，不改注释。
3. `run_current_identity_CellChat_v2.py ${RUN_ROOT} ${RSCRIPT}`：未来安全接续入口。当前已运行的v1控制器保留，不并发重启。
4. `render_ccrcc_current_identity_v2.py ${AGGREGATE_DIR}`：只读聚合表绘图，v2修正色条布局，统计量未变。

输入必须是既有已验收7成员、原`spatial_pilot_v1/identity_candidates_PRIVATE.tsv.gz`及已验证的CellChat2.2.0数据库；不从新聚类视图替换身份。代码不包含矩阵、细胞或坐标行，真实路径由调用者提供。

`acquire*`/`prefetch*`/`resume_after_prefetch*`是已保留的下载恢复代码。其他`identity_refinement`/`program_concordance`/`coarse_program_mapping`及相应画图/审核代码为历史诊断旁支，用户已要求后置；不是当前主线的必跑依赖。所有旧输出不覆盖。

[当前方法](../../methods/V4-ccRCC-current-identity-communication-20261109.md) · [结果](../../inventories/analysis/2026-11-09/ccrcc_followup_v1/README.md)。
