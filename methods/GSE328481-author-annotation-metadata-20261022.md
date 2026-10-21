# 作者空间身份/区域字段的限定资格

## 问题与范围

2026-10-14免疫图谱CM02 v1（T_whole_TME_primary_tumor）的组织落实支线；此前严格RNA状态家族保持。字段存在不是可用细胞身份，论文转移方法不是当前文件已资格，也不将新身份定义冒充旧RNA模式的精确复现。此阶段没有新空间效应、临床模型或正式重注释。

只选此前效果前已固定的同一最小公开processed试点，保持原URL/长度/ETag/Last-Modified。源数据片段和全部类别/codes仅RAM，客户端不保存H5AD，服务器不新写文件。仅解析anno、anno2、classes、bin50 structural annotation、TLScluster，后续标签校核另允许var基因名词表；不解析表达、坐标、细胞/材料/患者标识、克隆序列或临床值。

## 版本和接受门槛

v1尝试只读旧服务器试点路径，缺失即结束，不自动补下载。v2另冻原URL的metadata-only范围桥至原h5py运行时；保留失败，不把generic AssertionError当已定位HTTP或解析原因。已有schema显示嵌套数值category格式，纯合成重放证实旧reader会TypeError；大heap纯合成未复现读取上限错误，不据猜测扩预算。

v3仅修复跳过数字/嵌套category容器，不能把数值cluster码当生物类型，原字符串词表和资源上限不改。v4另固定生物学标签披露校核：每个词须属于明确细胞/区域词汇或原生唯一基因名，拒绝patient/sample/donor/id以及常见材料/测序编号格式。未知标签不原样输出。标签与源级计数不是患者记录，也不是身份真值；未匹配不当不存在。

范围桥固定64KiB缓存块、最多16MiB总网络字节及1MiB单次h5读取；无自动重试/替代URL。每个请求验证206、Content-Range及原头部。partial metadata读取无法重新验证整个原文件SHA，显式保留；原始全文件指纹仅从此前完整输入继承，不夸大全文件接受。HTTP/认证/长度/签名/预算边界即停。

本地cache/seek/readinto/EOF及cap纯核、原h5py内存H5通过adapter和未知值/序列/标识不输出纯核先于真实metadata字节。分母完整核，缺失与未识别分开，不调用原data矩阵或读取被拒对象。公开作者代码API返回限额时停止，不换路绕过。

## 解释和依赖

真实source terminal/退出码、schema/字段分母和全部未知接受后才说明入口能力。单试点不资格全部来源、逐材料患者独立性、取材/治疗阶段或病理终点；source annotation命名不能自动对应另一个scRNA队列的状态。新空间问题仍须参考/身份/测量/组织位置资格并另冻。科学标签/覆盖和优先级私有，公共只通用代码、纯方法与技术SHA；原CM/身份/Zarr/闭家族及TableS1用户后补安排保持。

参考仅为[公开原研究](https://www.nature.com/articles/s41467-026-75405-x)及[GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE328481)。本轮直接主文页跳认证，只读到公开索引的方法摘录，未据此称阅读全文或绕过认证。
