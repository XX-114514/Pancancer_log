# 逐文库双细胞与环境RNA处理 v1

## 适用范围

新增公共队列的源内预处理与整合准备，真实计数只在执行服务器。9个可辨识原始文库的计数分支，共792文库；4个已技术QC队列993,312细胞与5个待初始QC队列2,345,161细胞。所有原始计数、旧QC、失败尝试保留。TNBC原始文库映射未确定、CTCL混合来源排重未闭合，继续隔离。

## 双细胞

既有4队列复用已独立验收的逐文库Scrublet分数和1930个移除决策。另用DoubletDetection 4.3.0.post1在floor/mt通过、去双细胞之前的计数独立拟合，保存两方法交集、额外候选与额外移除。正式参数：Louvain、10轮、p_thresh=1e-16、voter_thresh=0.5、标准化scaling、pseudocount=0.1、10000可变基因、30 PCs上限、boost_rate默认0.25、n_jobs=1，seed=20261011+原始文库顺序。依据[作者教程](https://doubletdetection.readthedocs.io/en/stable/tutorial.html)。初次肾癌试运行使用Leiden，软件提示该后端实验性，故正式版本改用教程中的Louvain；试运行结果不能视为正式运行结果。

总双细胞移除上限8%每文库，保留旧决策并在剩余额度中按独立分数移除。超限候选保留flags，相关文库暂停整合放行。<500细胞文库不强制独立双细胞移除、只保留历史标记并注明未充分评估，禁止按固定4%凑数。SPECTRUM复用作者doublet标记；其余新队列缺少已核验方法时登记not_previously_called，不把零标记解释为已确认单细胞。

## 初始QC与作者记录

既有4QC对象不重复应用floor/mt阈值。5个native分支按原始文库评估floor(200 genes/500 UMI)及mt min(25,max(15,median+3MAD))；mt不可测明确登记。转移性乳腺癌全部作者PASS=True，复用该决定，记录通用floor/mt诊断而不无依据重新按通用阈值筛掉作者保留细胞。原始gene覆盖完整保留。

## 环境RNA

4个已QC队列当前均无核验通过的空液滴背景。`raw/X`是原始计数槽名，不能证明包含未过滤液滴。DecontX 1.8.0在每个原始文库中执行filtered-only估计；使用合并的粗谱系诊断分群，不使用T细胞细状态。无足够第二谱系（至少50细胞及1%规模）时，明确not_identifiable_without_empty_droplets，保持计数，禁止伪称去污染完成。

DecontX参数：背景缺失、maxIter=500、delta=(10,10)、estimateDelta=True、convergence=.001、varGenes=5000、确定性seed。谱系诊断由源内粗Leiden cluster的多个marker表达及优势比产生；CAR-T将作者T细胞分类合并成T_NK、单核细胞合并myeloid，unknown单列。诊断分群不等于最终细胞身份。模型为[官方方法说明](https://bioconductor.org/packages/release/bioc/vignettes/decontX/inst/doc/decontX.html)中的无背景分支，不能恢复未在细胞矩阵中出现的环境来源。

正式运行从方法日志提取最后一次iteration和convergence数值，须小于.001才允许采用校正层；达到500轮仍未收敛或缺少该证据时，保留估计并暂停校正放行。独立审计另外解析方法日志，不只相信producer的converged字段。

保持`counts`=校正前原始整数、`counts_decontx`=浮点估计。验证校正不增加计数、非负有限、谱系marker保留。主要谱系marker损失超过50%、文库中位UMI损失>30%或超过20%细胞损失>50%时，校正暂停采用，normalized X退回原始counts并标记复核。单细胞UMI损失>50%另外标记，不进入准备放行掩码；原始细胞和估计层仍保存。此门槛是保守复核规则，不是污染真值。

## 输出与准备验收

完整逐文库Zarr保存原始counts、独立污染估计层、明确来源的normalize_total(1e4)/log1p X、2000源内HVG、原始计数PCA/邻域及源内UMAP诊断。保存doublet表、过滤原因、marker前后统计、参数、污染率、图和逐阶段时间。小文库可以跳过不具可辨识性的模型／embedding，仍保留计数与限制记录。

GSE286827按stable ID应用已验收patientfix_v2；不推断疗效/生存。技术准备放行与source_final_AC分开：来源排重尚未闭合时最终新增AC仍为0，不能自动开始全局整合。CAR-T产品、T富集snRNA、普通肿瘤组织、排序富集文库及正常对照保留各自上下文。

独立验收按所有保留细胞扫描原始counts与上游逐项一致、ID及基因顺序、临床上下文、校正约束、标准化数值与放行门槛。输出技术候选清单与各源基因覆盖交集；正式整合HVG需在适当分支中batch-aware重新选择，源内2000 HVG不直接当全局最终集合。原始基因缺失必须作为不可测信息，不能当真实零表达。

## 资源与运行边界

4源并发、每源4文库进程，总最多16模型进程；BLAS/OMP/NUMBA每进程1线程。实际对象全量保存，不用抽样冒充处理。持久心跳每20秒，失败保留尝试目录、明确非零退出；完成后自动独立验收。全局整合尚未执行。
