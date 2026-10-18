# 20261018 GSE305078临床连接资格 v1

## 范围和输入

只读原官方SOFT与一个公开补充文件的格式资格；无模型/基因学习或科学效果测试。[方法](../methods/GSE305078-clinical-link-qualification-20261018.md)。原MCP/CM、Zarr、临床标签及固定统计合同不变。

公开正文来源：[2026年论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC13355034/)。可读方法中的总体治疗分组是后续精确连接的线索，不是当前28样本的已审核分组。原2025全文入口取读未成功；没有从二手评论采纳临床标签。

## 执行和实际结果

- 2026-10-18 12:36:10.139670 UTC /20:36:10.139670 UTC+8：scanpy环境启动PID60468/startticks314207628，源码SHA 2e1c70b8a683cd8a00be99e9a089eca4d0decf26519634ba76183b78183365a6。
- 12:36:17.424945 UTC /20:36:17.424945 UTC+8：终态PUBLIC_SUPPLEMENT_UNAVAILABLE_NO_LABEL_ADOPTED，SHA72981884c349ababdca4c2883fc74c08634a5836d75de5271c88ed6afe24be39。
- 请求预声明公开mmc1.docx链接、上限8MiB，实际HTTP200/text-html/1,817字节；PK ZIP断言正确拒收，未作为DOCX/临床表解析。仅保留小响应，无表达/测序下载、标签采纳或越权访问。

逻辑根`${RUN_ROOT}/GSE305078_public_clinical_link_qualification_20261018_v1`；plan先于网络请求，原SOFT SHA核对后只做元数据资格。

```bash
python scripts/qualify_GSE305078_public_clinical_link_v1.py
```

运行时PANCANCER_PREFLIGHT_ROOT指定真实根；公共日志不复制该根、响应正文或任何样本记录。尚未读取补充文档内容，不能判断图像中的患者映射；精确治疗/pCR连接仍未建立。不需要重启原MCP任务。
