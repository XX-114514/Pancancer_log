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

## 2026-10-16 10:38:46 UTC /18:38:46 UTC+8 全量排序已验收

生产实际终态10:36:53.545381 UTC /18:36:53.545381 UTC+8；独立审计10:38:46.138966 UTC /18:38:46.138966 UTC+8通过全部24背景、58848条手工并列AP/AUC、全部960固定槽及字面参考接口。生产终态SHA `a41336a209b5b37630d15b30d5298978a2c55909cf38a8a79c1cc7568b50f825`，独立审计SHA `a1fc980c77179d1a105ae8c0c27d78089b076c0c6dbf94f09e8e5ac25cdfe7a9`。候选槽SHA `5e1111291c9d63c63075b258d83bc22858b04dfb8ca9e952f5b453391433bc2f`，原接口SHA `06cd774f0df8db311176c11ea0ba76a2465f44735fe82c7e36e7ed825623224d`。

无新临床/通讯效果或P/BH。由此启动[原生同细胞测量及完整输入预检](20261016_NSCLC_prior_joint_RNA_and_response_preflight_v1.md)，不能把先验排序验收称轴验证。
