# V4血管免疫空间对照有限流程

- 范围：2026-11-06 V4 CM03的乳腺癌组成扩展；六唯一供者，原成员/亚型/权重不改。
- 输入：公开图谱六Xenium原生log1p-UMI与18层级映射；只有完整长度/CRC/SHA接受后读取。公开仓库不含病例行、矩阵、坐标或实际路径。
- 代码基线：4fefaaa36c6d996048bc42e24a1b08350707c41f；当前代码SHA见冻结v5。
- 环境：隔离Python3.9、AnnData0.10.9、NumPy2.0.2、Pandas2.3.3、SciPy1.13.1、h5py3.14.0；绘图隔离Python3.12/matplotlib。
- 计划激活：2026-11-07T13:37 UTC（2026-11-08 Melbourne），18小时有限依赖执行，无自动Git/通知。
- 当前状态/取证：2026-11-07T13:39:33.195033+00:00；18映射通过、原生源待完整接受，0真实病例生物学验收。
- 参数：50UMI/20genes；身份多marker；25/50/100/200μm邻域；100μm区域/血管带；50μmLR边；999条件置换；患者四家族2/10/27/27。
- 产物：私有完整身份/计数旁表、全66假设、身份及3轴位置图；公开仅审阅后的聚合/代码。
- 失败：早期适配器receipt字典语法错误在真实矩阵读取前发现，保留私有失败版本，v3起修正并端到端验证；旧截断下载不接受，分块body长度/Content-Range/CRC/SHA有限恢复。
- 原Visium/通讯队列不动，当前快照不刷新其运行时间。

可复现调用（先在隔离环境验明源清单/计数语义；真实路径从外部注入）：

```bash
python scripts/v4_vascular_controls_20261108/test_vascular_controls_invariants_v2.py
python scripts/v4_vascular_controls_20261108/test_native_vascular_end_to_end_v4.py "${NEW_SYNTHETIC_OUTPUT}"
python scripts/v4_vascular_controls_20261108/test_patient_level_families_v1.py "${NEW_PATIENT_TEST_OUTPUT}"
BREAST_PARENT_WORK="${ACQUISITION_WORK}" VASCULAR_CONTROLS_WORK="${CONTROL_WORK}" VASCULAR_PLOTTING_PYTHON="${PLOT_PYTHON}" VASCULAR_MPL_CACHE="${MPL_CACHE}" python scripts/v4_vascular_controls_20261108/finite_native_vascular_controls_v5.py
```

控制器使用既有逻辑工作清单；需在`${CONTROL_WORK}`部署当前代码、以原名称部署冻结计划，`${ACQUISITION_WORK}`放作者来源/映射清单和完整接收凭据。已存在输出目录会拒绝覆盖；不把本地未下载依赖视作成功。患者桥仅私有。

[方法](../methods/V4-vascular-spatial-controls-20261108.md) · [检查点](../inventories/analysis/2026-11-08/v4_vascular_controls_v1/runtime_checkpoint_v1.json) · [验证](../inventories/analysis/2026-11-08/v4_vascular_controls_v1/software_verification_v1.json)。
