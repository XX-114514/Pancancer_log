# 已有数据筛选与分区集合执行

## 范围

用户要求优先复用已存数据，判断免疫图谱用途、实际建立合适集合、交付筛选后数据结构，再依据硬盘空间安排后续下载。沿用[分区集合方法](../methods/backed-existing-atlas-collection-v1.md)。科学来源表、具体分母、候选排序和绝对路径只保留获准私有交付；本记录保存技术执行。

## 输入和环境

既有免疫权威对象、公共预处理权威对象及冻结行选择合同；原泛癌全RNA分片与冻结CoVarNet输入参考。合同、发布凭据、选择索引和执行代码分别SHA绑定，参见[技术收据](../inventories/analysis/2026-10-23/existing_immune_collection_v1/technical_execution_receipt.json)。

源运行环境实际为Python3.10.19、NumPy1.26.4、SciPy1.11.4、Zarr2.18.3、AnnData0.11.4。通用代码在dirty工作树按SHA执行，归档时核本机代码与服务器执行代码SHA一致；不是假称事先已经检出本轮最终commit。

## 执行时间线与恢复

- 2026-10-23T07:19 UTC附近：旧探针强制代理导致控制连接关闭；核现成直连与同一服务器/认证配置后仅用于控制，不修改全局SSH/Clash。
- 07:20:55 UTC：只读库存实际exit0，核现有对象及磁盘；随后补组配额、源结构、作者类型与原隔离问题。
- 07:43 UTC附近：先通过真实源环境人工数值测试，首次集合写入后真实读取发现canonical基因ID是DataFrame索引。旧目录保留，没有接受失败结果。
- 修正读取器并将真实索引形态加入反例，在独立attempt_v2执行。控制调用等待超时后未重启原任务，先核同一输出与原进程。
- **07:46:40.324322 UTC /15:46:40.324322 UTC+8**：完整集合终态写出，数值测试实际exit0、全部分区两种视图直接CSR回算和缺失mask检查通过。后续控制核查原源进程不存在；不把超时的控制调用声称为actual exit0。
- **07:51:26.028206 UTC /15:51:26.028206 UTC+8**：独立哈希、来源总量、原发布凭据及旧原生参考入口核查通过，工具实际exit0。
- **08:00:20.760774 UTC /16:00:20.760774 UTC+8**：源软件版本、实际文件规模和最新组配额补核完成，actual exit0。

## 产物与验收

`${RUN_ROOT}/20261023_existing_immune_ecology_collection_v1_attempt_v2`包含config、source coverage、通用读取器/测试、provenance及legacy reference扩展。保存的是小型引用集合；没有新增人类表达下载、完整矩阵复制、联合模型拟合或旧身份修改。

私有工作簿包含现有库存筛选、实际集合分区、旧冻结参考、身份/状态/测量字段合同、容量和下载条件；公式与源汇总独立核对、xlsx结构/错误单元检查及各工作表视觉审阅通过。未决来源维持明确的表达参照、隔离或待完整性审查角色。

## 命令与下一步

```bash
python tests/test_backed_atlas_collection_v1.py -v
python -m py_compile scripts/backed_atlas_collection_v1.py tests/test_backed_atlas_collection_v1.py
bash scripts/validate_repository.sh
```

具体源bootstrap和带来源映射的执行请求留私有，不公开个体行或科学排序。下一先复用现成集合，审研究别名/捕获与亚型状态对应；对部分文件核完整性，预算允许后续传缺段。旧图谱与CM不变，公共push沿用此前拒绝未重试，云端待同步。
