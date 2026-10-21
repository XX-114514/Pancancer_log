# CM02同输入本地RAM患者与原生邻近测量

## 实际状态

[有界片段资格](20261021_CM02_local_RAM_fallback_preflight_v1.md)通过后，首个完整测量于2026-10-21T07:57:26.296246 UTC /15:57:26.296246 UTC+8启动，但08:00:52.773032 UTC /16:00:52.773032 UTC+8真实失败，未读CSV。失败终态和原进程结束已确认。合并读取/同域名受验证地址选择经实测后，修正版于08:20:58.052225 UTC /16:20:58.052225 UTC+8启动。当前完整CRC/SHA、患者覆盖及图未接受。原远端终态仍不可读、未重启；此为同源计算，非新增证据分母。原CM02时间/范围及权威身份/家族保持。

## 输入与实现

先核先前完整输入资格及其terminal绑定SHA；同一官方ZIP目录必须两次等于原接受摘要，所选CSV完整CRC与SHA256 `45fdf21775bf38fc3cf869ed60d45cb51caac0b1076b2a686572961082589c8e`、字节数和独立CSV/Pandas计数/均值复算必须通过。明确使用缓存固定版本清单，不声明新鲜API前后重放或完整ZIP MD5验算。小型区域/患者桥/通用作者代码保持固定SHA，作者ROI必须对全集双射，未知患者核心区域不外推，零计数与缺失均值分开。

首个入口`...local_v1.py` SHA256 `99435bdb692251315bd443722c767435b483e8c7e0a093307ddc0eded653132a`保留；当前入口`summarize_melanoma_IMC_patient_coverage_local_v3.py` SHA256 `5f16d955a706d4666b029d15a738b654dfaf99269736de59111ea8e7948b44d2`。邻居helper仍为`acbde69dae387ee79d6d8397899681bedf00553dd3da08878abd16b392858471`。十个科学函数与服务器v2的AST相同。实际Python3.13.9/NumPy2.4.1/Pandas3.0.3对当前后端运行原覆盖/节点合同21项、0 skip，含13个纯科学函数AST校验、19/20、ROI非病例及错误输入hash拒绝。精确邻居helper另四项检查通过，没有真实大ROI/跨运行逐数接受声明。

默认安全HTTPS及HTTP/1.1协商，64KiB分块，35秒socket/45秒请求硬上限/三次瞬时重试、原64MiB含失败请求及75分钟range预算。所有个体/表达/坐标仅RAM，Mac RSS按字节，仍执行2GiB相同门槛；不安装软件，不改网络配置或信任库。片段速度只作预算参考，不保证全量。

## 接续与边界

首个节点模块SHA `e8c3193b1176a80c793d2648f2dcdeb8aeb74b2ae831120d06c87a0b163301c1`保留，未因父任务失败而启动。当前节点`...local_v3.py` SHA256 `62bb043a017d66f109e6db883760281d9a41cb314c188bf41e166324c3d80b09`已排入修正工作流，包含于上述21项检查；覆盖真实完成才执行。目录/节点CRC/标签/歧义/不同全集规则保持，不作新鲜API重放声明。

## 失败、修复及严格连接资格

首个完整进程仅成功读取22字节，后续20字节目录探针连续网络超时；本机ZIP源码确认OSError可被包装成BadZipFile，因此不是归档已损坏或生物学阴性。64KiB分块RAM缓存通过五项独立字节/完整合成ZIP CRC/清缓存强制网络重放检查；最大八块，所有物理过读/重试仍占原请求预算。新鲜目录核对前明确清缓存。

合并片段首轮随后遇到证书主机名不匹配，立即停止，未使用错误证书。对当前系统DNS返回的四个官方主机地址仅做完整TLS诊断：两个超时、两个验证成功；此不证明某个运营方根因。新连接仅在当前DNS候选中选择能完全验证原HTTPS主机名者，不用未知IP、不改DNS/SNI/信任/代理、不向不匹配连接发HTTP；外层45秒硬期限不能被吞。五项证书失败/全部失败/硬期限/过期地址/禁止关闭验证检查通过，stdlib和系统未修改。

同一官方route1的新片段于08:18:06.669295 UTC /16:18:06.669295 UTC+8真实通过，目录清缓存后从网络重读一致；片段用时约1.5秒，仅技术资格不代表全量速度保证。旧终态与原PID结束确认后才启动修正版，非因观察超时盲目重启。

覆盖后最多75分钟RAM保留、最多三个另固定自有模块。当前仅索引排队，没有状态/临床效应请求。原生k10、全部类型/十marker/全部类型对及患者门槛保持；不是原RNA状态、物理µm或独立验证接受。下一先核本地同PID/开始身份、完整终态及节点结果，再固定身份和状态组织家族。原完整临床模板信息不足、GSE288199 TableS1待补及受拒动作保持，整体goal未完成。

## 复现

```bash
PANCANCER_IMC_ACCEPTED_MEASUREMENT="${ACCEPTED_FULL_MEMBER_JSON}" \
PANCANCER_RAM_OUTPUT="${LOCAL_RUN_ROOT}/private_output" \
PANCANCER_IMC_TRANSPORT_HELPER="${LOCAL_RUN_ROOT}/scripts/preflight_melanoma_IMC_public_archive_v2.py" \
PANCANCER_IMC_BOUNDED_TRANSPORT="${LOCAL_RUN_ROOT}/scripts/melanoma_IMC_coalesced_ranges_v2.py" \
PANCANCER_RETAIN_RAM_FOR_FROZEN_STEP=1 \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
python "${LOCAL_RUN_ROOT}/scripts/summarize_melanoma_IMC_patient_coverage_local_v3.py"
```

自有源复制后逐SHA核对；真实PID/开始身份、依赖/协议/输入hash与stdout/stderr路径仅私有launch receipt。公共只代码/方法/技术hash，来源效应和科学数量私有。此为单次计算，不创建定时自动化。
