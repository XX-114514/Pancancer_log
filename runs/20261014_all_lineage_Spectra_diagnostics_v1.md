# 全18大类Spectra稳定性和身份内变化检查

原6block×3seed仍训练；本阶段尚未运行真实诊断。逐完整block接受三个producer completed-seed记录及模型/metadata，然后检查每个factor、scope内三个Hungarian匹配和cycle。

[事前合同](../inventories/analysis/2026-10-14/all_lineage_spectra_diagnostic_plan_v1.json)。已有模型合同的“支持/变异”具体化：同cell至少5签名基因/至少20发现细胞，全部3pair source-residual score rho>=.7；该20为可测支持计数，不是样本级推断资格或状态阳性阈值。权重mincos.9、top50Jaccard.5，所有major结果保留。block_global共享至少2major；只有1major支持为context候选，保留但不叫跨type共享。单major Glial只能context/global或自身lineage，不制造共享。

原node内及来源、RNA/detected深度、同cell原生RNA份额和支持全部列出。签名用匹配3seed L2 weights的median top50；官方factors已包含gene scaling。与旧T/M strict all-permitted-major规则历史并列，不改旧8/173结果。

入口：`python diagnose_all_lineage_spectra_v1.py --block <fixed_block> --plan <frozen_plan>`，环境PANCANCER_PREFLIGHT_ROOT；读取既有模型/原生输入，全程不做cell级p、不按LR/临床筛选、不改172/9CM、不发布正式state。

完整block诊断保存新版本目录；后续独立对所有匹配、major支持和统计复算后冻结候选签名/投影。只有seed文件存在不当作fit完成，须有producer completed记录。
