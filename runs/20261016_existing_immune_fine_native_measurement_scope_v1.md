# 既有免疫档案完整T/髓系细分单位与原生测量入口

## 已冻结范围

在[外部输入别名审计](20261016_existing_immune_candidate_extension_preflight_v1.md)通过后，按原节点和物理文库登记全部29个T/髓系节点、6957单位、1745198原细胞，包含不足20细胞的单位，不为增加通过率先合并细分标签。270节点、原172身份与2026年10月CM（v20261010）不变。原生计数层仍为整档2670057细胞×49563物理特征、4586971705个CSR存储条目，尚未本轮读取表达值。

范围生产2026-10-16T12:28:54.423731+00:00 /20:28:54.423731 UTC+8启动、12:29:03.252395 UTC完成；源码SHA `a7720e62e4a393fa2f71a22b7215ebae28083cbf8ebfb1096851bd809aef333b`，终态SHA `24c4af6a2703d4cdb5843d1b8327b2430184b88e75c2cb41510562295e94d7bd`。单位表SHA `a48215e707721d078506e9b632dec50131be81045270f27de9e8abbce622f413`，原细胞—单位映射SHA `28e9c022522ea6c3b23204690e7a689b8f8c3211044d6f30df180ab6a515808f`。

## 资源纠正与拟采用策略

首资源读取匹配quota列表第一文件系统，可能取得不限额卷的零值，不能将其false gate解释为目标数据卷配额耗尽。新v2按真实输出路径的df文件系统精确核对quota；原单位/映射、500GB十进制预留、256MiB输出预算和1GB额外余量不变，没有删除文件、要求用户再次申请配额或缩小科学范围。修正终态以真实后续收据为准。

后续按完整原始计数和来源掩码计算native RNA分母，在内存保留必要的全基因均值工作区，输出固定程序/全部既定调控签名/有限基因QA与完整覆盖；避免新增完整稠密基因矩阵副本。该策略尚未执行表达评分，不构成评分验收。

```bash
python scripts/prepare_existing_immune_fine_native_measurement_scope_v1.py
python scripts/correct_existing_immune_resource_preflight_v2.py
```

逻辑根`${RUN_ROOT}/existing_immune_fine_native_measurement_scope_20261016_v1`及`existing_immune_native_resource_preflight_20261016_v2`。精确NSCLC作者亚群对应仍须效果前审阅；宽谱系扩展、原细身份复现与临床验证分别登记，不互相替代。无新下载、表达效果、标签改写或独立验证声明。

资源v2实际12:31:52.631718 UTC /20:31:52.631718 UTC+8启动，12:31:52.768994 UTC完成并通过原门槛；源码SHA `0ee76ad7b60717926aabcf269e76429aeefa3ba5cac6cfb1ad5b2df68c4d129c`，终态SHA `18b8ea074dbbcb534454682dc40dd8edb8dc86e15bb0e072342a4192038b9426`，修正计划SHA `c3031ee1f3cb16e0f1d988bcd882d7408f4a7e7363f62ca7f210f3ade44e769e`。真实目标文件系统组配额剩余502907535360字节，所需501268435456字节；500GB十进制预留、256MiB输出上限及1GB余量均未降低。资源快照不是后续永久保证，正式写出前再核。旧false来自首文件系统匹配错误，不是实际配额耗尽，原失败记录保留。
