# Thrane 辅助来源的表达、坐标和元数据对应

2026-11-20 02:24–02:25 UTC / 10:24–10:25 UTC+8，已下载 spatialGE 作者维护的数据包并完成独立读回。它提供了 4 片的辅助元数据来源，**没有解决原作者患者桥接或 H&E 配准**。本节点用于 2026-11-15 WEHI R4/K15 的 CM12 空间验证资格，不新增生物学检验。

| 切片 | spots | 基因行 | 去除 Ensembl 后缀后的重复 symbol 行 |
| --- | ---: | ---: | ---: |
| ST_mel1_rep2 | 293 | 16148 | 19 |
| ST_mel2_rep1 | 383 | 16831 | 23 |
| ST_mel3_rep1 | 256 | 15653 | 23 |
| ST_mel4_rep2 | 248 | 15991 | 18 |

全部 **19,150,573 个计数值、1,180 个 spot** 与已经验收的作者原始归档对应一致。辅助表仅删去了基因标签中的 Ensembl 后缀，未改变计数；没有自动合并重复 symbol。mapping 中 xpos/ypos 恰好是条码拆出的整数坐标，并非像素或已配准的组织坐标。辅助包的 4 行临床 metadata 与上述四片一一对应，给出 4 个不同 patient 标签，但本轮没有原作者逐患者原件复核，不能据此升级为正式独立患者验证，也不向缺少的另四片外推。

[spatialGE 官方教程](https://fridleylab.github.io/spatialGE/articles/spatial_differential_expression.html)指向该包，明确其含切片元数据。下载自[作者维护仓库](https://github.com/FridleyLab/spatialGE_Data)，相应文件最后变更提交为 `6f53608c0beb59c38ca2a1e28d038ce4344fc322`。Git LFS 声明大小 3,400,344 bytes、SHA256 `c7957972360dd4753d6d49764d98d8a7f63c070ec2de4e4d968f2fc51881f32c`，实际字节及完整 ZIP CRC 均一致。教程和包属于后续方法作者维护来源，不能与 Thrane 原始临床证据等同。

服务器直连下载因 Proxy CONNECT aborted 失败；改用本机流式传输同一公开文件直接存服务器成功，没有改代理配置。初次清单读取误把 Mac 资源旁文件当 CSV 导致 UTF-8 解码错误；正式脚本限定 `melanoma_thrane/` 内容，全部 ZIP 成员仍经过 CRC 检查，保留失败记录。原论文 AACR 正式页及作者稿页本次返回 403，未绕过访问控制。

## 验证与保存

[完整审计](complete.json)记录原始/辅助文件 SHA、所有切片与资格边界；[独立 pandas 回放](independent_replay.json)核对全部数值、行列轴、重复 symbol 和条码坐标，PASS。两个实现均没有使用目标 marker 或临床结局来挑选样本。

原始包和临床行只存服务器 `${RESULTS_ROOT}/CM_validation_20261118_v1/sources/spatialGE_Thrane_v1/`；聚合输出在相应 `spatial_melanoma/spatialGE_bridge_v1/`。公开目录不含表达矩阵、逐 spot 坐标或临床患者行。Python 3.10.19 / scanpy 环境，复现命令：

```sh
python audit_spatialGE_Thrane_bridge.py "${SOURCES_ROOT}" "${OUTPUT_ROOT}"
python replay_spatialGE_Thrane_bridge.py "${SOURCES_ROOT}" "${OUTPUT_ROOT}"
```

现有八片直接 marker 分析已断言每个目标 marker 最多一行，本次重复 symbol 发现未推翻该结果。没有把辅助包替换为新模型输入，也未缩减八片分析范围。下一步仍需原始患者对应资料、条码—图像独立锚点及组织/分割资格；不能用这套整数网格替代配准。
