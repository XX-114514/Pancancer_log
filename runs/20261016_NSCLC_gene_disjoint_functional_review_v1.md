# NSCLC完整基因不重叠功能方向与病例留一描述

## 固定范围与实际终态

接续[全功能测量](20261016_NSCLC_gene_disjoint_receiver_TF_v1.md)。保留全部2568父记录—调控签名槽，1570有输入资格、998缺覆盖。仅描述发送配体与接收靶RNA读数的关系，没有新P/BH、置信区间、阈值筛选或原LR家族改动；同来源筛选后的跟进不当独立确认。

生产2026-10-16T12:02:24.448715+00:00 /20:02:24.448715 UTC+8启动，12:02:50.774029 UTC完成；源码SHA `cb8cff28b9fb61cca04f051aef9033fa8ecdbc5f7fd502f0c1ecad5f326ae06f`，终态SHA `8adfc820d89ce1bd13677ab73f2f4adf022154834718ef57a086fcca21839b1b`。

独立最小二乘实现12:03:08.729765 UTC /20:03:08.729765 UTC+8启动，12:03:50.327917 UTC验收全部点、额外控制和259725病例留一；最大点差4.996e-16。审计源码SHA `c76b61f94ca58acf6c037cdbbcdf42a45563a9eb8916907b28e9dfb6c723fda1`，审计终态SHA `07d644fad739843ae8103d7795e7d0b4630a21de955184992cb42e8f74b20d2f`。

```bash
python scripts/describe_NSCLC_gene_disjoint_TF_response_v1.py
python scripts/audit_NSCLC_gene_disjoint_TF_descriptions_v1.py
```

逻辑根`${RUN_ROOT}/NSCLC_gene_disjoint_TF_response_descriptions_20261016_v1`，所有科学表与病例留一只服务器/授权对话。

## 描述及解释

秩化原发送配体与接收签名，保留原18原尺度标准化控制加截距；另分别加原发送程序或病理响应，及描述接收签名与原masked-MCP的关系。病理控制可能产生collider偏倚，不是新的主要模型。每个原病例删除后重排秩和重核全秩/条件数；失败保持缺失，不删协变量。

调控靶RNA模式不是调控因子蛋白、DNA结合或实际活性。公开短路径的符号也不等于配体一定激活该受体。原完整性/结合/释放/同细胞RNA限制仍保留。完整结果提示需要进一步核签名冗余和来源共同程序，故对全部签名继续[特异性诊断](20261016_NSCLC_regulator_specificity_v1.md)，不只选择方向符合预期的项。

原2026年10月CM（v20261010）、172身份和本来源作者标签均未变更；尚未形成独立验证或整体goal完成。
