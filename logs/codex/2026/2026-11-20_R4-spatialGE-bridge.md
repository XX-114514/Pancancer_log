# Thrane 辅助来源接收与资格审计

用户要求持续推进 CM12 及其他 R4 模块的空间/临床验证并发布审阅结果。上一轮 CM13 平台审计及发布构成实质进展；本轮开始隔离检出 b560203 干净，原目录并行修改保留。沿当前 R4 WP1/WP4 工作计划，读取状态及方法并查原始来源，不重跑已完成模型。

GSE14520 原始临床页仍要求浏览器验证；没有重复拟合或借二次论文推断终点。空间线检索发现 spatialGE 官方教程与作者维护的 GitHub LFS 数据包。GitHub 插件用于读取目录、LFS 指针和提交元数据；二进制不受该接口支持，使用 HTTPS 流式下载到服务器。

2026-11-20 02:24:55 UTC / 10:24:55 UTC+8 主审计 exit 0；02:25:34 UTC / 10:25:34 UTC+8 独立 pandas 回放 exit 0。四片全部 19,150,573 计数值、1,180 spots 与原始 ZIP 对应，坐标等于条码整数；辅助临床记录四行，但原作者患者桥接仍未验收。来源、SHA 和每片分母见[报告](../../../inventories/analysis/2026-11-20/WEHI_R4_spatialGE_bridge/README.md)。

失败与修正：服务器直接 curl 的 Proxy CONNECT aborted；同一文件经本机流式 SSH 传输成功。初次清单因 __MACOSX 资源文件误读发生 UnicodeDecodeError，正式代码限制内容前缀并保留完整 CRC 校验。AACR 两种原文页面返回 403，没有绕过。

新建两份参数化审计脚本、聚合回执及说明；原始资料在服务器。未读取临床结局用于挑样本/模型，未发布临床行；没有使用 derivative patient 标签宣布原作者验证成功。重复 symbol 问题单列；原直接 marker 脚本已有每个目标最多一行的断言，不需要因此重算。

检查命令包括 git status、定向 rg、GitHub fetch、curl 流式接收、Python 主审计/独立重放及 SHA 检查。下一步获取原作者桥接与图像锚点；四片辅助坐标不含组织变换。建议提交：`Verify Thrane derivative counts and coordinate provenance`。仓库检查、常规 push 和远端 SHA 以实际回执为准。

Validation completed: repository checker exit 0, 0 errors / 1 historical path warning. New manifest, received complete SHA, script syntax, privacy and staged whitespace checks passed.
