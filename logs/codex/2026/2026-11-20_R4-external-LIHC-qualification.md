# R4 CM13 外部肝癌资格审查

用户要求继续实际研究并将审阅结果上传指定公开仓库。开始时隔离检出干净，HEAD 655d086；原工作区并行修改未触碰。读取 README、STATUS、ROADMAP、TODO、AGENTS、RUNTIME_SYNC、R4 工作包计划与组织日志；沿 WP4 查原始临床及表达资格，不重复 TCGA 拟合。

2026-11-20 02:19:00 UTC / 10:19:00 UTC+8，服务器只读输入审计完成，实际 exit 0。02:19:37 UTC / 10:19:37 UTC+8，独立 pandas 读回完成，实际 exit 0。输入为既有 GSE14520 GPL3921 矩阵、官方注释和冻结 R4 签名，输出新增 external_LIHC_qualification_v1；未覆盖旧结果。

执行：定向 find、读取旧处理脚本/矩阵元数据、原始论文及附件检索、官方目录 curl、两份 Python 审计脚本。观察到 445 样本、22,268 探针，225 肿瘤/220 非肿瘤；CM13 主签名 47/60、敏感性 40/60 有精确无歧义探针。120 覆盖行全部独立重放。输入 SHA、时间、输出 SHA 与范围在[结果目录](../../../inventories/analysis/2026-11-20/WEHI_R4_external_LIHC/README.md)。

失败尝试：官方 GEO 附件目录 HTTP 403、PMC 附件浏览器验证；未绕过访问控制。原始论文证明有临床资料，不证明已获得患者桥接。旧 GSE14520 endpoint rules 为空。既有 individual 标签数不等于患者数。未取得临床时间单位/起点/事件编码，未运行 Cox，不称外部验证成功；既往探索史保留。

新增聚合表、审计代码、来源说明和哈希清单，更新 STATUS/TODO 与 daily。无患者行、表达矩阵、原始临床附件或凭据进入公开提交。独立复核资格回执 SHA 为 bc73e09b06faaaba5cd16e22e3e94ca51c2415cc9500a00565c1074dde657f0d。

下一步获取作者临床原件及代码本；核定完整 GSE/平台范围、患者桥、时间原点和预冻结覆盖政策。WP1 仍需独立图像/条码锚点与病理区域资格。建议提交：`Audit CM13 external LIHC measurement and clinical eligibility`。发布前执行全仓库验证、diff/隐私检查及普通 push，发布回执单独保存。

Publication checks: repository validation exit 0; 0 errors and 1 historical internal-path warning. New-file manifest, source syntax, received hashes, privacy review and staged diff check passed. No unrelated files staged.
