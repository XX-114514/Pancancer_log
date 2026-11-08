# 六份Xenium空间坐标单位修复

问题来自作者预处理中的单位与字段名不一致，以及本项目此前据字段名作出的相等校验。源坐标本身未损坏。

## 根因与证据

[固定版本作者图像处理代码](https://github.com/swarbricklab/P08-breast-tme-spatial-atlas/blob/177e35ffed015b5a6155319c830e03616efc8659/scripts/01_preprocessing/xenium/02_handle_images/01_make_adata_objects_with_images.py)先取原生坐标，除以`microns_per_pixel=0.2125`成为像素；之后把这一已转换数组保存到`xenium_um_spatial`。`obsm['spatial']`则是另外经过逆对齐变换的H&E像素坐标。[初始对象代码](https://github.com/swarbricklab/P08-breast-tme-spatial-atlas/blob/177e35ffed015b5a6155319c830e03616efc8659/scripts/01_preprocessing/xenium/01_make_objects/01a_make_adata_objects.py)说明原生空间坐标来自obs质心。

| 实际字段 | 六份文件中的实际含义 | 本次用途 |
| --- | --- | --- |
| obs[x_centroid,y_centroid] | 原生Xenium微米质心 | 邻域半径、LR边、网格和距离计算的唯一输入 |
| obsm[xenium_um_spatial] | 原生微米÷0.2125，实际原生像素，名称误导 | 乘回0.2125后独立交叉核对，不直接用于微米距离 |
| obsm[spatial] | 对齐后的H&E图像像素 | 图像展示需匹配该图像坐标，不能与原生微米点直接叠放 |

实际六源共3,262,951条坐标全部检查：修正单位解释后的最大绝对误差1.819×10⁻¹² μm；每源固定均匀索引1,000个查询细胞，原生与恢复坐标的最近邻距离最大差2.034×10⁻¹² μm。没有拟合比例/平移来强行匹配，也没有重排条码或改写原生对象。

## 代码改动与拒收条件

新增`native_coordinate_contract_v1.py`，适配器v6在读取RNA矩阵前调用。读取原生obs微米坐标，检查二维、有限及非退化；旧key若本身为真微米允许直接核对，若符合本作者已核实像素定义则×0.2125核对。绝对容差1e-4 μm、相对容差0；其他比例、平移、轴互换、行错位和缺失校验数组一律拒收。每例输出实际单位、换算因子、误差及最终坐标源凭据。

0.2125仅用于这批已核实作者对象，不作为其他Xenium数据的通用比例。没有放宽成任意`allclose`、没有自动回归仿射变换、没有把H&E像素拿来算距离。源文件与作者标签均保留。

此前v5六源均在错误相等断言处终止，主空间结果0/6；不会有一批错误邻域结果被静默覆盖。原错误假设和合成测试只覆盖真微米key，这是原测试遗漏；新回归同时覆盖实际作者像素编码和真微米编码。保留v5及全部失败记录，使用新服务器run重放。

## 验证与执行

本机与科研服务器均通过坐标回归及原生格式端到端合成测试；已验证50μm邻接关系，同时拒收未知缩放/轴/顺序/偏移/NaN/形状/缺失。合成不是病人生物学结论。真实六源坐标审核为全量，聚合凭据见[审核](../inventories/analysis/2026-11-08/v4_coordinate_fix_v1/all6_coordinate_audit_aggregate.json)。

采用[冻结v6](../inventories/analysis/2026-11-08/v4_coordinate_fix_v1/plan_v6.json)，保留原成员/权重、身份门槛、25/50/100/200μm、999置换和全部66假设。现有53Visium/其他空间数据不被这批Xenium的比例更改。

按本地存储约定，新原生输入副本和大型输出放在既有科研服务器的正式项目run。六原生+18映射以普通SSH/rsync复制，逐文件大小/SHA256核验后才写ready凭据；所有本地原件仍保留，无删除。服务器18小时有限流程已部署并启动，等待完整输入；运行与生物学验收分别登记。[状态](../inventories/analysis/2026-11-08/v4_coordinate_fix_v1/runtime_checkpoint_v1.json)是取证快照，传输启动不代表全部接受，服务器启动也不代表六例算完。

后续传输补充：本地上传实测缓慢；默认请求继承的历史临时代理失效，但仅本进程禁用该失效代理后，官方4KB范围请求206/0.98秒通过。新有限直接下载仅使用公开源，逐成员CRC/大小/既有SHA验证；输入与输出使用独立`inputs_direct_v1`/`outputs_direct_v1`，避免与仍活动的原复制任务并发写入。没有更改持久网络配置、没有停止旧上传。此运输改动不改变科学参数。
