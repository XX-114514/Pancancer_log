# cNMF无先验状态程序对照 v1

UTC 2026-10-14T01:26:54.724585+00:00/UTC+8 2026-10-14T09:26:54.724585+08:00。[固定合同](../inventories/analysis/2026-10-14/cnmf_state_contrast_plan_v1.json)。与Spectra同12000发现cells，按T/M各6000；官方1.7.0，K10/15/20各100重复，完整600run及consensus，不按任何LR/临床p选择。2000HVGFano选择无外部程序先验，max_NMF_iter1000/defaultdensity.5。

真实官方prepare/get_norm_counts源码核查：externalTPM用于HVG/谱统计，rawcounts基因单位方差用于NMF；native完整RNA分母用于TPM，不能声称NMF本身使用相同Spectra logCP或已消除深度。绝对/normalizedusage并列后续分析。PANCANCER_PREFLIGHT_ROOT/PANCANCER_CNMF_PLAN与隔离venv入口`python train_cnmf_state_contrast_v1.py`。

状态：predeclared_not_started。原输入数据全审接受。zero-HVG-cell或缺run不自动删/跳过，失败保留。完整goal的UCell/全部9CM状态及LR/DIALOGUE/空间仍待；本阶段不晋升正式state。

## 实际运行：2026-10-14 09:37 UTC+8取证

01:32:01 UTC worker123088/child123089已启动，01:37:19两/proc命令身份匹配；164factorization文件/完整600run+consensus未结束。沿同句柄等待真实terminal，不重启。当前完整goal仍active、无正式state。[checkpoint](../inventories/analysis/2026-10-14/state_models_progress_v1/RESULTS.md)。

真实terminal02:03:23 UTC /10:03:23 UTC+8，全600/6共识exit0、进程退出；02:04全对照及02:08独立统计readback完成。原Kstats无全独立重建重演，候选不是state；[全对照](../inventories/analysis/2026-10-14/state_method_comparison_v2/RESULTS.md)。
