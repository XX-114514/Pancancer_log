# V7大类＋marker命名与新冻结

用户要求正式身份标签去掉注释进度词：有可靠亚型才使用亚型，否则使用**大类＋可代表该群身份的marker**。旧V7只作历史追溯，不覆盖；应激、IFN、增殖、同种型、良恶性及复核状态独立存放。不会把DNAJB1/CST3或“缺少亚型”本身造为新细胞类型。

当前为**全量证据审计完成、规则边界收紧及新冻结进行中**的检查点，不能将本页早期建议表当最终身份。最终名称与纳入人数以后续独立冻结审计为准。

## 全量执行与审计

2026-11-06 09:58:45—10:08:02 UTC实际读取全部4,676,787行原始RNA，54来源、490基因、223个审阅分区，109,270条群级marker统计。真实supervisor 10:08:07 terminal0；10:14:00 UTC独立全轴原始计数、完整库UMI、22个正向组件、来源基因轴与全部表达统计审计PASS，最大CP10K误差6.82e-13。QC-v4仍为8,237排除／238暂存／4,668,312保留，未重新删除原对象或覆盖此前QC。

审阅分区保留此前已重新建模群；其他部分按原全局群与旧标签分层用于检查，**不是这次新拟合了223群**，也不能以旧父标签分层制造纯度。单细胞正向组件、独立分类器旁证和逐群marker共同审阅，不能只看群平均或最高分。全轴原始及侧表文件留私有，公开群级汇总和通用脚本。

## 已发现的命名与身份问题

| 旧入口 | 全量新证据 | 处理原则 |
|---|---|---|
| 三个B coarse分区，共111,535个QC保留细胞 | 原两provisional分区的B独立组件比例约1.3%、5.5%；pending分区约4.8%；有T、髓系及DC方向 | 不能整体重命名为B＋marker；正向身份路由和父类纠正分开 |
| B DNAJB1/CST3 | 前轮已证实DNAJB1多数B、CST3父类混合 | 采用B身份marker，状态单列；混合不强命名 |
| 旧pDC与部分Plasma coarse | pDC有LILRA4/CLEC4C/TCF4/IL3RA，也常有MZB1/DERL3/JCHAIN | 共享分泌不能自动doublet，必须用pDC特异组件作对照 |
| 旧Granulocyte | 部分群有真正FCGR3B/CSF3R/CXCR2组合；CSF3R/S100又与单核共享 | 保留真粒细胞证据，收紧以FCGR3B/CEACAM8结构组件为入口，不能由S100单独命名 |
| 旧Macrophage_SERPINE2 | PMEL/MLANA/MITF强，原scATOMIC多为Skin Cancer Cell；另有SOX10/PLP1/MBP共享 | 先核黑色素父身份；SOX10/MBP单独不能把它叫胶质；良恶性仍独立 |
| 旧myCAF某分区 | RGS5/PDGFRB/MYH11与COL1A1/DCN并存 | 需血管壁与成纤维边界复查，不能仅保留功能myCAF名称 |

## 名称与冻结资格分开

“Bcell_MS4A1_CD79A”表达B大类的命名marker，不等于新B亚型，也不要求群内每个细胞恰好检出这两个转录本。命名基因来自可测谱系组件；泛HLA、应激、周期及免疫球蛋白不独自承担父身份。支持大类而不支持细分的群，可以在**大类＋marker粒度**冻结，不以必须细分为门槛。

真正父身份仍混合或信息不足的记录保留完整审计与排查原因；不得通过换名伪造确定身份。新版本正式身份列中不放unresolved/coarse/provisional/pending/review/candidate等进度词；旧词仅可出现在明确标为历史的来源字段。进入冻结集合与当前QC合格是两项独立条件，后者不是发表级身份保证。完整图谱的发表验收也不能由“没有占位词”推导。

## 当前证据文件

- [全轴独立审计](full_axis_diagnostic/independent_profile_audit.json)、[组件汇总](full_axis_diagnostic/all_group_component_summary.tsv)、[既有分类器旁证](full_axis_diagnostic/all_group_reference_top5.tsv)。
- [v1诊断命名建议](v1_diagnostic_name_proposals/all_group_name_proposals.tsv)仅为历史诊断，尚非冻结；22个v1规则及反例保留。
- [v2收紧规则](identity_policy_v2.json)及[23个反例检查](identity_counterexample_results_v2.json)；粒细胞、胶质及ASC组件正在按全量证据重算/独立审计，原490基因原始统计不变。
- [通用代码SHA](code_SHA256.tsv)、[方法](../../../../methods/V7-major-marker-freeze-20261106.md)、[正式运行](../../../../runs/20261106_v7_publication_marker_annotation_v1.md)。

新TME输入与CoVarNet没有因本任务替换；本轮用于旧V7结果的身份重整。未发现的新身份或新亚型不会写成已接受。
