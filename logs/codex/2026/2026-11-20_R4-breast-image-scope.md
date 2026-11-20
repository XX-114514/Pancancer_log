# R4 成像资格：恢复后的六片图像检查

用户恢复 active 目标。上一轮因用户暂停只做只读检查，没有新科学进展；本轮从现存状态重新确认，未重启任何旧模型。隔离仓库原 HEAD 8260280，干净；原目录并行工作未触碰。读当前 README/STATUS/ROADMAP/TODO、R4 工作包方法及相关成像代码，按 WP1 检查实际可用图像。

肺癌官方网页可浏览，三个候选文件 HEAD 403，网页 CLI 429；不反复请求。定向查乳腺文件时，自动审批拒绝完整输出私有接收配置，理由为可能含凭据/敏感路径。按更小范围替代，仅列文件名和 H5AD schema，不再读取完整配置；后续安全操作获准，无待批准操作。

六个既有 H5AD 均有 RGB 图像、缩放及细胞面积/核计数，没有多边形。部署主审计与独立范围复核，均最终实际 exit 0。主审计源 SHA 与旧接受回执核对通过，仿射数值和 72 分组指标读回通过。实际时间存于[回执目录](../../../inventories/analysis/2026-11-20/WEHI_R4_breast_image_scope/README.md)，UTC 与 UTC+8 相差八小时。

目视检查完成全部六片总览，点云总体随组织覆盖；未将形态解释为血管/肿瘤诊断。S01 517 个原生质心越界，pc_low_MYH11 分支含 7 ESM1+ EC 和 1 RGS5+ mural；其余五片无越界。无核/多核范围完整保留，不能直接判错误。

失败修正：首次逐位坐标相等断言失败；六片最大差均 2.842e-14 µm，明确容差 1e-10 µm 后通过。未修改已有坐标或数据。新增代码、聚合报告和回执，原始图像/预览及细胞数据不公开。当前没有新 P 值、没有根据结果删除样本。

下一步按固定全部六片/三分支/原参数族，执行图像内单核技术敏感性；同时继续补原始多边形/独立配准。它不能代替病理资格或多患者验证。建议提交：`Audit embedded breast images and segmentation-related analysis scope`。发布执行普通 push，独立核验远端 SHA。

Publication checks completed: repository validator exit 0, 0 errors / 1 historical path warning; source manifest, complete SHA readback, syntax, privacy and staged diff checks passed.
