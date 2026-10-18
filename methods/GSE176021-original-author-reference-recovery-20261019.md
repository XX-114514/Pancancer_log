# GSE176021原作者参考标签恢复与临床元数据资格

## 当前问题

既有raw H5AD虽有CD3/CD8作者标签列，但全部为空占位；原作者两个注释文件已缓存，不需重新下载表达。它们是外gzip再包一层原gzip RDS，单层gzfile读取报unknown input format。显式外gzip流再由R gzcon读内层RDS，恢复原560916×7和235851×7表；原值完整保留，不按旧图片或新结果修改作者标签。

## 精确细胞键

作者barcode是原文库命名空间加下划线及10X条码，不是跨文库裸barcode。只允许raw imid + underscore + raw barcode的唯一精确匹配；作者imid与原patient字段必须完全一致。全部作者行登记，未匹配仍缺失，不用模糊匹配、近似条码或首条重复值补标。

当前图谱的source_row_position与既有完整source_rows manifest精确相同，parent_cell_id逐条等于原H5AD该行index；本来源original_native_row_position为空占位、original_cell_id不是此完整原生index，因此不能无审查替代。旁表保存当前全部868494来源行、原始行、原作者行/标签code及全部未匹配作者行索引；原Zarr、当前major/subcluster、原CM与旧统计家族保持。

## 覆盖与独立接受

全部CD3/CD8精确键、患者字段、当前图谱行/ID/标签与全部组织及既有R/E/My捕获池覆盖均另实现全量核对：生产用pandas Index、无引号TSV，审核独立字典映射及R带引号CSV，逐行比较完整旁表，并用独立计数器核全部覆盖。全捕获R/E重叠仅诊断程序测量对象，不重新拟合、检验或命名正式状态。参考标签是作者定义，不能自动升级为本项目正式细胞状态。

当前引用原论文[PMID34290408](https://pubmed.ncbi.nlm.nih.gov/34290408/)和[作者勘误](https://www.nature.com/articles/s41586-021-03893-6)。勘误涉及部分图中名称/位置等，不能凭图片自动交换2024缓存原注释的值；如正式采用这些细分状态，需另核作者版本及原RNA标记。临床MPR、原抗原特异性与当前学习程序彼此分开。

## 论文方法与病理反应资格

从[Europe PMC官方REST](https://europepmc.org/RestfulWebService)完整取得原XML，验证DOI及原Methods：新辅助nivolumab两剂后手术，MPR为残存活肿瘤不超过10%；scRNA使用冻存复苏后高活率CD3分选T细胞。肿瘤/邻正常/LN与脑转移、纵向血液及另设抗原刺激实验分别保留，不能把整个GSE当培养细胞或whole-TME，也不能把手术后关联当基线预测。

原XML明确Supplementary Tables文件，但精确逐病例/材料桥仍需完整原表。服务器CDN TLS失败、官方supplementaryFiles ZIP分块传输失败、已授权Clash经一次回环转发且TLS校验保持的CDN请求仍失败；三份失败回执保留，未采纳任何MPR标签。转发随命令结束，未留下后台代理进程；没有新测序表达下载。论文中的定义/总体组数不当个体标签；没有完整映射，不启动该临床关联。

## 下一阶段

先用完整作者参考身份核原冻结捕获视图及程序的细胞构成/覆盖，逐病例与材料聚合避免细胞伪重复。若需要在作者定义Treg/CD8群重新测量，先另固定完整两程序/原基因轴、全组织及配对病例范围，核变换/训练域及全部资格后再执行；作为新测量问题，不能据原阳性端点选择新对比或重开旧家族。病理反应临床支线待完整原表和精确病例桥。全部科学组成/效应及个体旁表只存服务器及授权汇总对话；public仅源码、方法、维度与SHA。完整互作目标仍active。
