# ccRCC既有身份的空间区域分解

目标：定位2026-11-06 V4 CM03大类扩展及既有空间LR候选的贡献，不增加身份审核前置。身份/图/轴/半径/分层沿[上一轮方法](V4-ccRCC-current-identity-communication-20261109.md)，无标签重算或CellChat重拟合。

## 精确期望

LR：在每个身份×floor(x/tile)×floor(y/tile)×floor(log2 RNA深度)层内，求4基因二元profile均值。发送和接收身份互不相同，置换独立；每条边配体/受体双阳性期望为两端相应层的阳性比例乘积。固定近邻边相加即原profile置换的精确均值，不改变原999次P/q。

EC—B/NK边：仅在原B/T/NK免疫细胞的位置×深度层中求该身份比例，沿原EC连接边相加。共享EC：每个EC每层有k个免疫邻居，层大小n，B数b、NK数m。无B概率choose(n−b,k)/choose(n,k)，无NK类似，两者皆无用n−b−m；层间概率相乘后容斥1−P(no B)−P(no NK)+P(neither)。零邻居EC期望0，B/NK互斥不冒充独立Bernoulli。

## 区域与显示

500μm固定网格。LR及EC单边按接收者位置归块；共享EC按EC位置归块。各块保留observed、expected、opportunities和excess，在私有表中保存。公开仅12/108各设置的总体、区域方向汇总及固定期望删一区域范围，无区域坐标行。

删除每块贡献计算100×(全量超额−块超额)/(全量机会−块机会)，不重新构建分层。此为影响诊断，不是jackknife区间或新P。所有有机会块进入删一；≥20机会块用于正向比例和效应图。表达密度用100μm栅格，QC总数<10栅格隐去；各图独立色标，原生微米不做像素换算。

## 验证和范围

独立pandas分组均值核对全72个LR期望；对12/108核区域总和/删一结果；80个共享EC以整数choose核对，并用完整枚举反例证明须考虑不放回互斥。保留原生身份SHA、坐标不变；原P/q只作旧结果，不将本描述结果作为新推断。单供者、相关空间块和参数变体不扩充患者数；没有病理图像，组织类型/血管壁未重建。

[预设设计与实际结果](../inventories/analysis/2026-11-10/ccrcc_localization_v1/README.md)。主脚本`localize_existing_identity_v1.py`，审核`audit_localization_v1.py`，绘图`render_localization_v2.py`（v1保留为早期展示历史）位于`scripts/ccrcc_localization_20261110/`。调用：

```bash
python -B localize_existing_identity_v1.py "${RUN_ROOT}" "${LOCALIZATION_OUTPUT}"
python -B audit_localization_v1.py "${RUN_ROOT}" "${LOCALIZATION_OUTPUT}"
python -B render_localization_v2.py "${LOCALIZATION_OUTPUT}"
```

使用既有scanpy服务器环境（Python3.10.19、NumPy1.26.4、SciPy1.11.4/pandas/Matplotlib），无新安装。细胞/矩阵/坐标/区域原表及完整日志留服务器。

展示v2将效应块的原生500μm网格边界换成与100μm栅格相同的参考原点，并统一物理显示范围；不修改细胞位置或统计。初版仅显示原点分别归零，未作为当前图发布。
