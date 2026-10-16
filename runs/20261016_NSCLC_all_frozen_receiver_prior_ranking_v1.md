# NSCLC全部固定接收程序先验：启动节点

## 实际执行

生产2026-10-16T10:33:03.493922+00:00 /18:33:03.493922 UTC+8启动；PID/startticks原收据已核，10:35后取证实际R状态且进度已更新，尚未取到正式终态。生产源码SHA `ebff9a8535b6806b8e0e64edba7bc64c2c42e469ee42583a610c1e2422816169`。环境及私有根沿用[两端run](20261016_NSCLC_two_end_coordination_v1.md)。

方法：[全部固定接收程序自动先验合同](../methods/NSCLC-all-frozen-receiver-prior-ranking-v1.md)。范围24程序、48正负前景、全部1226配体排名及960固定候选槽；无临床或协调结果筛选，无新P/BH。

```bash
python scripts/rank_NSCLC_all_frozen_receiver_priors_v1.py
python scripts/audit_NSCLC_all_frozen_receiver_priors_v1.py
```

独立审计按实际原进程和终态等待，不因观察超时重启。逻辑目录`${RUN_ROOT}/NSCLC_all_frozen_receiver_prior_ranking_20261016_v1`；`plan_before_any_ranking_or_sender_effect.json`、`private_progress.json`及后续`private_terminal.json`/`private_independent_audit.json`分别证明合同、进度、生产和审计。

审计2026-10-16T10:34:21.994954+00:00 /18:34:21.994954 UTC+8实际启动，源码SHA `f3b1daafc8fe9b45eb97111f0c36d0b3d01a027519ab3b315c4695534032bd3c`；生产/审计PID及startticks后续查询均相符，分别R及依赖等待S状态。尚未把等待状态称审计通过。

## 当前限制

本节点为输入/排序运行中，不是已验证候选、受体同细胞共表达或直接互作；候选实际结果/审计另补。原native单位/基因轴使用完整既有验收，无新矩阵下载。科学排序、基因和具体结果仅服务器及获准对话，不进入公共仓库。原身份/9CM和所有失败保留。
