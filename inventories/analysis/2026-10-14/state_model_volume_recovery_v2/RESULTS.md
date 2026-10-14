# 模型主卷容量保护与备用卷恢复v2

原模型v1确实停止：cNMF 05:03:55 UTC在B完成后，Spectra05:09:33在T3seed/Myeloid1seed完成后，均exit1/500+1GiB门槛。此前“其余继续训练”的状态以真实终态纠正，不能把旧RUNNING快照继续当live。

05:34取证主卷组quota可用约479GiB，相对启动约507GiB减少约28GiB，属于共享组总使用变化；未审计到具体归属，不凭此归因单个进程。主库/history不删除、不重写，原输入全部只读。

已有个人data1区可写、无硬quota，复制后free1053.761GiB；恢复输出继续500GiB留空+5GiB预算。301182371bytes/359个完成模型文件逐SHA相同：T3/Myeloid1及B300/3共识，不复制453万fullRNA。新私人目录0700。v2仅output卷/边界恢复/已完成组合跳过，模型科学字段逐项与v1相同。

最初复制包装语法换行失败，服务器无动作；本机嵌套代码先编译/正确转义后完成复制。原模型failed记录、原结果与复制ledger保留。模型副本恢复完成不是finalstate或canonical权威迁移；临时输出需后续完整审计/持久交接。

[恢复参数](../all_lineage_spectra_model_plan_v2.json)、[cNMF参数](../all_lineage_cnmf_model_plan_v2.json)。作业实际启动/进程证明另补。原9CM/172/Zarr保持，全scope/UCell/状态LR/DIALOGUE/空间仍待。

恢复版05:58:48/49 UTC真实启动，06:06 UTC四proc命令一致；Spectra已复用T3/Myeloid1、在Myeloidseed2中，cNMF已复用B300/3、在DC中。actualhandoff/refresh135目录0冲突。详见registration_public_receipt.json。仍无新formalstate/全9CM/LR/DIALOGUE/空间完成。
