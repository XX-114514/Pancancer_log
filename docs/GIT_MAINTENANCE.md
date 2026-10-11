# Git maintenance and release workflow

## 2026-10-11 授权与公共仓库纠正（优先于历史private条款）

用户本次明确指定 `https://github.com/XX-114514/Pancancer_log` 并授权整理、commit和普通push；已核实目标为public。历史private-only描述保留为当时背景。本次及后续已授权范围内可向该目标正常提交经过审查的代码、汇总和脱敏日志；禁止force push、历史改写、删除和上传患者行/矩阵/凭据。

Cloud开发遵循[工作流](CODEX_CLOUD_DEVELOPMENT.md)及[RUNTIME_SYNC](RUNTIME_SYNC.md)。每个任务必须写Codex日志，执行轻量检查，保留历史。Cloud不能假定能访问数据服务器或本机skill；仓库中的合同与规则是交接依据。缺少真实数据时只交付工具和合成验证，并明确真实执行未验证。

## Repository role

The maintained repository is `${PANCANCER_ROOT}/project-records`, preserving the
history of the existing `Pancancer_log` remote. It is the project control plane
for status, methods, provenance, lightweight releases and manuscript assets. The
adjacent Project_v1-v5 analysis trees are not Git payloads.

## Storage contract

Track:

- Markdown methods, decisions, run records and task logs;
- code and configuration after portability/privacy review;
- aggregate TSV/CSV/JSON tables;
- publication PDF/PNG figures;
- checksums and logical paths for external objects.

Do not track raw human data, cell/patient-level matrices, H5AD/RDS/Zarr/pickle
objects, complete runtime environments, scheduler output, credentials or files
larger than 5 MiB. Git LFS is not installed on the current host, so large objects
must remain external.

## Normal update

```bash
cd "${PANCANCER_ROOT}/project-records"
git status --short
bash scripts/validate_repository.sh
git diff --stat
git diff
git add <explicitly-reviewed-files>
git diff --cached --stat
git diff --cached
git commit -m "type: concise description"
git push origin main
```

Push only after confirming `origin` is the intended private repository and
authentication is available. Never use force push or rewrite frozen history.

## Release rule

A frozen release is append-only. Do not overwrite a V7 file after publication.
Any biological/statistical change creates a new run and version (V8 or later).
A corrected narrative that does not change data may be added as a dated erratum.

Each release should contain:

1. README with scope and version boundaries;
2. machine-readable version metadata;
3. small reports, source tables and figures;
4. external artifact manifest with hashes;
5. source-code provenance;
6. repository validation evidence.

## Before Prism or supplementary export

Use [the Prism and supplementary guide](PRISM_SUPPLEMENTARY_GUIDE.md), verify that
all plotted statistics match the source table, and cite the exact release ID.
Do not mix V7 annotation statistics with V5-derived CoVarNet/LIANA results without
an explicit version statement.
