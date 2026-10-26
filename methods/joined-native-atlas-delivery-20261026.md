# 只读全细胞图谱交付接口

本文件描述[全细胞分层注释任务](full-context-atlas-hierarchical-annotation-plan-20261023.md)的交付接口。实际执行与验收以[同任务记录](../runs/20261023_full_context_hierarchical_annotation_v1.md)为准；代码存在不等于真实运行完成。

## 读取方式

服务器上使用该次交付目录中冻结的代码和manifest，继续读取原分区counts。不要把manifest当作自包含Zarr，也不要把输出中的候选身份标志删除后当成新的权威细胞类型。

```python
from joined_native_atlas_delivery_v1 import JoinedNativeAtlas

atlas = JoinedNativeAtlas(delivery_manifest)
source = atlas.sources[0]
metadata = atlas.metadata(source, [0, 1, 0])
panel = atlas.read_counts(source, [0, 1, 0], ["CD3D", "LYZ"])
embedding = atlas.projection(source, "global_overview__shared", [0, 1, 0])
```

每次读取最多4096个来源内位置；顺序和重复位置保留。使用`source`与`native_cell_id`的有序对识别细胞，`native_row`用于绑定底层存储行。

## 字段含义

- `author_label`保留作者证据，`local_major`和`local_subcluster_id`来自来源内分群；`consensus_major/celltype`保留已冻结的逐细胞保护规则输出。
- `confidence`是粗细身份的证据类别组合，不是校准后的概率。`identity_release_status`明确候选层级，`label_version`绑定该次版本。
- `global_subcluster_id`是精确来源内群体的结构ID；展示字段为`celltype_marker_number`。一个marker用于显示，完整本群多基因证据仍在绑定文件。数字不同不自动代表不同共同生物类型。
- SHA单列。fit排除细胞保留在档案中，缺少拟合亚群时不制造编号或借父群marker。
- 原生样本、文库、患者字段以`native__`前缀保留。已审阅的分析块与原生表达单位另列；未知单位不升级为独立患者。
- 六个`RNA_state__`连续读出保留原生表达和测量覆盖。它们不直接确证功能状态，γδ候选也不因存在RNA信号升级为正交确认。
- 原有QC、CNV可评估性、结果、恶性证据版本分别保留，缺字段单列；不通过谱系推断恶性状态。

## 缺失和适用范围

counts必须与`measured`共同读取。未测基因的占位值不代表生物学零，原生数值也不自动证明UMI测量方案。

每个投影包含独立模型适用行轴。`model_domain_available`区分未进入该模型的细胞；缺坐标使用masked数组。原始PCA和候选校正表示并列，raw-only分支不制造校正坐标。拟合暴露、查询适配、几何距离旗标与身份可信度分开。经过资格验收的精确距离覆盖读取时显式应用，原数组和原失败记录保留。

`view_status__*`保留档案、表达参照、完整组织组成、谱系内组成、产品和治疗临床的分别资格。技术读取成功不改变这些资格。尚未加入的科学审阅证据保持unknown；既往已接受的分析仍使用其原冻结范围，不能据新视图待审把它们改写为失败。

## 验收职责

上游完整逐细胞/群体/状态/单位轴的独立审计持续绑定。新交付检查每个来源的计数和mask直接回读、所有精确群体/证据类别代表、每个模型的适用及缺失域、全部精确距离覆盖。计数直接回读为有界探针，不能称全部原生count逐值重算。

定量边界复核使用原已审计的全基因充分统计，分别保留表达量、检测率和固定RNA程序。测量基因集合不同的来源不合并成一个同口径均值；缺测面板保持缺失。本群均值不能证明同细胞共表达，二元谱系共支持不直接等于doublet或污染。

新共同身份、来源设计/重复资格与下游节点粒度仍各有审阅职责。旧带日期与范围的CM继续冻结。
