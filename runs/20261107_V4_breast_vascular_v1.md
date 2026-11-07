# V4乳腺癌血管免疫组成独立分支

逻辑根`${PROJECT_ROOT}/runs/20261107_V4_breast_vascular_v1`。用户将主线收束为血管相关组成—空间组织—募集通讯—治疗关联；状态不是依赖。原V4、53 Visium和25 CPDB队列均保持。

开始取证12:17 UTC／20:17 UTC+8；publication起点8f7b9630，独立新工作树，旧主区15项tracked修改保留。只读旧成像检查12:19:55 UTC确认7/13原生单元数学审核完成：5 Xenium各4变体、首两CosMx112/80变体；身份、MC校准及精确CM03并未一并接受。

新资源12:18–12:22核Zenodo四归档、面板5101名称、六Xenium/18映射成员及作者13代码，固定提交177e35ffed015b5a6155319c830e03616efc8659。检测IGHD/XCL1缺测、Tangram质量语义、原生/图像坐标和CP10K层含义；不把官方记录中约5100与实际CSV5101名称强行相等。

12:23新有限服务器分支启动，12:23:56直连归档URLError后终止，无原生矩阵验收。原失败/代码保持；增加仅loopback、仅两固定公共归档的临时HTTP范围中转，沿既有Tailscale SSH，单个有限远程命令结束即关闭，不改持久proxy/SSH。恢复仍需实际CRC/SHA取证，不仅凭进程活跃宣称下载成功。

I-SPY2复用本机先前12:08:47官方完整来源SHA，未重下载；新增临床分析均在本任务独立目录。12:25:04／20:25:04 UTC+8测量确认248位/69+179、19,134基因；严格三基因panT缺TRAC使24槽均不可估计，CM03另缺NCR3LG1、CM09另缺C1QL4/DMRTA2。尚无可估计效应前增加测量补充v2，用实测CD3D/CD3E作控制，签名/病人/结局/模型门槛保持，缺原程序仍NA。12:26:08生产18/24可估计；12:26:43独立原始键、秩分数、18设计/似然/FDR通过，最大1.19e-10，主21/增量3槽各0 BH。原严格表、原代码及原因表保留。图完整24槽、实际CI与NA均保存，已视觉核对并修正文案/坐标文字重叠。

执行环境：服务器既有scanpy Python3.10及同环境anndata/scipy，单CPU、无GPU。临床生产本机R4.5.2；原始测量/独立核验本机隔离Python3.12/NumPy/Pandas，未新增共享包。新通用代码及绑定SHA见清单，外部作者源码只读取证不整段转载。

脱敏复现命令形状：

```text
${PYTHON} selective_public_zip.py --url ${PUBLIC_ARCHIVE} --output ${PRIVATE_MEMBER_INVENTORY}
${PYTHON} finite_breast_acquisition.py ${BREAST_RUN_ROOT}
BREAST_WORK=${LOCAL_BREAST_WORK} ISPY2_SOURCE=${LOCAL_VERIFIED_SOURCE} V4_SIGNATURES=${FROZEN_V4_SIGNATURES} ${LOCAL_PYTHON} prepare_ISPY2_inputs_v2.py
${LOCAL_RSCRIPT} --vanilla run_ISPY2_logistic.R ${LOCAL_BREAST_WORK}/clinical_v2
BREAST_WORK=${LOCAL_BREAST_WORK} ISPY2_SOURCE=${LOCAL_VERIFIED_SOURCE} V4_SIGNATURES=${FROZEN_V4_SIGNATURES} ${LOCAL_PYTHON} audit_ISPY2_models.py
BREAST_AGGREGATE_ROOT=${REVIEWED_AGGREGATES} ${PLOTTING_PYTHON} render_ISPY2_aggregate.py
```

已有输出目录拒绝覆盖，不能照抄重复运行成功阶段。[当前报告和全部统计](../inventories/analysis/2026-11-07/v4_breast_vascular_v1/README.md)是科学解释权威；空间下载/身份/患者几何/通讯仍须后续接受。暂无新的精确CM03、HEV机制或治疗预测器结论。

补充传输/环境取证：独立SSH单字节范围探针实际206/1字节通过，但4MB块正文读超时，产生的BadZipFile是传输未完而非证明源包损坏。保留各版，新增256KB块/300秒读取/完整块重试，源和科学规则不变；不复用其他线程SSH主连接。服务器未验收新矩阵，不凭活跃进程称成功。为避免传输阻塞所有资格检查，本机先按原字典序取首份对象，3份对应层级映射已完整CRC/SHA通过，原生矩阵仍部分下载；有限依赖阶段只在全部成员CRC/SHA之后读取，不自动Git同步/通知。独立本机原生读取环境Python3.9/AnnData0.10.9/NumPy2.0.2/Pandas2.3.3/SciPy1.13.1/h5py3.14.0，与R临床及服务器环境分开，未改已有环境。

12:44:20 UTC／20:44:20 UTC+8取证：小块服务器传输v4活跃，但新增矩阵0完整/0原生资格结论；本机首份对应三层映射3/3 CRC/SHA通过，矩阵85,983,232 / 442,743,536字节仍部分，完整之后才由有限阶段读取。当前不发布任何新的空间身份、距离或通讯阳性。

较晚传输纠正：服务器临时大文件通道再次实际超时，不能称服务器下载队列健康。全部源/失败保留；本机已启动同六个既定Xenium及18映射的有限接续，复用首份原生下载与已CRC/SHA通过的3份映射，后续15映射不重下前3份；合并映射只硬链接已校验成员、不改变字节。首份资格阶段等完整收据后执行，其余5份逐个完整下载/资格；不做自动Git或通知。科学来源、顺序、身份/距离规则均不改。新增原生生物结论仍为0，当前状态以本机分支实际收据为准。
