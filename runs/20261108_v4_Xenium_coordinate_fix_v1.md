# V4 Xenium坐标修复重放

- 用户要求：修正空间坐标问题；新大型资料优先保存科研服务器。
- 事件/取证UTC：2026-11-08T00:41:41.638884+00:00；UTC+8：2026-11-08T08:41:41.638884+08:00。
- 根因：作者单位变换后把像素保存至误名um key；原v5未区分单位、六源全部断言失败。
- 输入：同六唯一供者的已CRC/SHA接受原生文件及18映射；不改对象/条码/冻结标签。
- 代码基线：c664a92cb36e82bf294b955ca36e43bde464135f；当前文件SHA见冻结v6。
- 本机验证：Python3.9、AnnData0.10.9、NumPy2.0.2、SciPy1.13.1；服务器验证：Python3.10.19、AnnData0.11.4、NumPy1.26.4、SciPy1.11.4。
- 服务器逻辑run：`${PROJECT_ROOT}/runs/20261108_V4_native_coordinate_fix_v1`；真实路径与复制清单留私有回执。
- 状态：代码/两套合成测试通过，服务器有限阶段已启动；24输入正在逐文件复制校验，未把未完成文件标ready。大型产物仅服务器生成，0本地删除。
- 输出：六源坐标合同、身份审核、12+27+27/患者全66、marker身份及三轴图；计算完成后还需科学/视觉审核。
- 旧失败/运行/参数保持；无自动Git或通知，不改既有Visium/CPDB队列。

复现：

```bash
python scripts/v4_vascular_controls_20261108/test_native_coordinate_contract_v1.py
python scripts/v4_vascular_controls_20261108/test_native_vascular_end_to_end_v5.py "${NEW_SYNTHETIC_OUTPUT}"
python scripts/v4_vascular_controls_20261108/run_coordinate_fixed_cohort_v2.py "${COORDINATE_FIX_RUN_ROOT}"
```

正式执行先部署同SHA scripts、私有config及逐文件大小/SHA通过的inputs。输入清单为六唯一供者，校验凭据写完才创建ready；输出目录已存在则拒绝覆盖。

[方法](../methods/V4-Xenium-coordinate-units-fix-20261108.md) · [当前技术快照](../inventories/analysis/2026-11-08/v4_coordinate_fix_v1/runtime_checkpoint_v1.json)。

运输更新：4KB官方源直连206通过，新服务器直接下载及独立目录runner v2已部署启动；旧上传/旧等待任务保留。默认失效历史代理只在新下载进程中跳过。所有新ready仍需CRC/SHA及坐标合同；未把网络探针当作完整文件接受。
