# 52来源全量UCell流式评分容量预检

[冻结计划](../inventories/analysis/2026-10-14/all_source_ucell_stream_preflight_plan_v1.json)。沿用已核对4532933 canonical细胞、52来源、18大类、172原节点和1290样本；核对原始行索引与每个来源真实基因轴。此次仅读取CSR行指针，计算各来源×大类细胞数、非零上界、分块数量与输出容量，不复制RNA、不开始评分。

容量上界为每细胞252个适用程序（45cNMF+最多207Spectra），每cell-programme22bytes：4种float32分数、uint8基因支持、float32 RNA份额和mask。按5000细胞输出块、500细胞rank批次估计；dense排名工作空间采用6倍float64保守估算。输出40GiB预算含8GiB元数据和聚合预留，输出卷保留500GiB。最终字典必须重新检查实际适用程序数量和空间。

主分数计划使用已验证官方pyUCell0.8.0、native全部可测基因、average ties、maxrank1500，1000/3000及cNMF TPM签名敏感性。只有完整字典验收后才进入正式评分；覆盖不足为缺失，真零保留，原身份与CM不变。此处容量通过不等于评分通过或状态发布。

全量评分入口project_all_source_ucell_stream_v1.py已经编写并通过语法检查，尚未执行或数值验收。入口按来源×大类分块读取原始CSR，以source native全部基因为排名分母，保留global row/cell ID/unit index及4种评分、基因支持、RNA份额。先要求完整24scope字典及独立dictionary acceptance，再重新检查空间，最终需所有4532933行恰好一次。此代码准备不能作为UCell已完成的证据。
