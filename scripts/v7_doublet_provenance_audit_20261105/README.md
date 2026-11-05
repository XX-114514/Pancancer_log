# V7 doublet处理历史审计脚本

这四份脚本与本轮实际远程执行的body字节相同；源根由执行器提供，没有嵌入服务器地址、账户或凭据。需要在能读取原Project_v3输入的既有Python3.10/scanpy环境运行，仅输出汇总JSON，不写入输入对象或重跑模型。

```bash
export PROJECT_ROOT=/path/to/Project_v3
python - scripts/v7_doublet_provenance_audit_20261105/audit_V7_doublet_membership_v1.py <<'PY'
import os, sys, runpy
from pathlib import Path
runpy.run_path(sys.argv[1], init_globals={"PROJECT_ROOT": Path(os.environ["PROJECT_ROOT"])})
PY
```

依次可选择：`audit_V7_doublet_membership_v1.py`（全成员/原判定）、`collect_historical_removal_totals.py`（原账本移除统计）、`independent_doublet_provenance_check.py`（完整多重集及原逐细胞账本独立复核）、`final_input_integrity_check.py`（冻结SHA/全分数状态）。更换上述命令的脚本路径即可。

输入是原QC文件及日志中明确命名的版本，不自动搜索“最新”输出。所有54来源的原预测/移除字段存在性由本轮发布校验另外显式检查；原预测False只对有原分数的细胞具有工具阴性含义。当前28,999个未运行记录的False是占位。

[方法](../../methods/V7-doublet-provenance-audit-20261105.md)及[数据结果](../../inventories/analysis/2026-11-05/v7_doublet_provenance_audit_v1/README.md)。公共仓库不含输入细胞条码、患者/样本行或表达矩阵；离线没有原输入时不能宣称已真实执行。
