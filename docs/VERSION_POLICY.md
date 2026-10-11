# Annotation and analysis version policy

## 2026-10-11 表达分支与身份版本边界

最新运行状态以[当前快照](../inventories/runtime/2026-10-11/REPORT.md)为准。reviewed V10 core assembly SHA-256为 `41875cd0e904ed46bcae93e5e8df4031e44fc1b71f1f11ff3b5f201c764a53e3`；此索引不是全局最终科学身份晋升声明。

`20261010_uniform_expression_annotation_v1` 用本群表达统计生成单基因描述标签，保留既有major/fine身份；其新CoVarNet/LIANA和CM来源内验证独立登记。表达标签可完成而最终细身份、CNV、生态和临床验证仍未完成。GSE299623最终身份门槛继续保留。历史V7投稿冻结、V8/V9候选及被降级20260923分支均不改写。

## Current authority

V7 is the default cell-identity and malignancy reference from 2026-09-06 onward.
V5/V6 remain historical inputs for already completed downstream analyses.

## Version boundaries

- Identity version: major and fine cell-type definitions.
- Malignancy version: candidate/reference construction, CNV engines and call rule.
- Downstream version: CoVarNet, LIANA, spatial and clinical analyses.
- Manuscript release: curated set of tables, figures and prose.

These axes must be recorded independently. A downstream analysis does not inherit
a newer identity version unless it is actually rerun.

## Increment rule

Create V8 or later when any of these change: cell universe, gene mapping,
candidate/reference definition, marker logic, CNV parameters, malignancy rule or
canonical field semantics. Cosmetic prose or figure-layout changes do not require
a biological annotation version bump, but they do require a dated release update.

## Freeze rule

Frozen artifacts are append-only and verified by hash. Never replace an object
under the same version label. Record corrections in a new run, release or erratum.
