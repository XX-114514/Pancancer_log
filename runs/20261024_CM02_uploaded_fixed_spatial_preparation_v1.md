# 外部TDLN用户上传输入：固定空间分析准备

## 范围与实际状态

用户已改为本机下载后SFTP上传，当前前缀仍在增长，完整MD5/SHA及真实R任务尚未开始。旧公网输入失败均发生于R前。此次仅准备[上传输入路线](../methods/CM02-user-uploaded-fixed-spatial-input-20261024.md)，原11R组件/拼接R和全部analysis逐项一致，不新增状态、尺度或统计机会。

冻结/验收观察：2026-10-23T13:04:59.501642+00:00 UTC / 2026-10-23T21:04:59.501642+08:00 UTC+8 / 2026-10-24T00:04:59.501642+11:00 Australia/Melbourne (AEDT)。新科学合同SHA：`2fcabbe5544171b5b4546f6bc48642114538ccfb2c9631a254e7668e78279e66`，真实作业未启动。

## 人工执行与修正

现成源Python3.10.19/R4.3.3的人工输入、原R管道和消费检查actual exit0。第一组后台反例暴露主程序退出后后代仍占stdout导致监督器等到超时；修正为就绪读取、轮询工作器、清理专属组后收EOF。另一次人工原位修改未被首遍buffer/元信息区分，已加入同描述符第二遍完整hash。所有失败、代码快照和旧合同保留，没有放宽反例或统计。

修复后五项后台反例从第二个独立SSH连接实际验收：延迟成功发生在首SSH真实结束后；exit7保留；个体canary未落帧；超时工作器已终止；孤立后代无活跃残留。源收据SHA：`b7997c755e0f196675942ac4cbf95b34df81460c374364e6c5877a4f0ae47f0e`。真实pending消费者actual exit3且未建接受输出，部署包小于1MiB。

本机人工R stdout因沙箱Operation not permitted未通过；同链已在真实源R通过，不假称本机R管道接受。源工作器actual wait与监督进程本身OS退出作用域分开；源汇总绑定也不是本机重读真实原生。

## 下一入口（待执行）

完整稳定双hash通过后才允许一次正式任务；只读用户已有文件并保留它，不另建原生缓存。以下命令不是本轮已运行：

```bash
PANCANCER_EXTERNAL_FIXED_TASK="${FROZEN_UPLOADED_TASK}" \
PANCANCER_READONLY_SSH_TARGET="${EXISTING_SOURCE_ALIAS}" \
PANCANCER_REMOTE_PYTHON="${EXISTING_SOURCE_PYTHON}" \
python -B scripts/launch_external_TDLN_detached_RAM_v1.py
```

后续同一控制名采集，uploaded消费者核完整输入/代码/实际退出和原15槽。科学汇总仍私有/获准对话，旧CM/身份/Zarr/已闭家族、临床资料依赖及公共push旧拒绝保持。
