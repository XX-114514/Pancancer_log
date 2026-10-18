# 原候选全外部来源家族推断

## 前序与固定参数

全部模型输入于2026-10-18 09:32:57.837581 UTC /17:32:57.837581 UTC+8通过1719ebcf16d4634bb803620c048e028605e1ab076356a255d694bc56b22dd46a，覆盖1040计划、2408 case行与6760文库flag；最大聚合差7.11e-15。完整all-unit主家族和reviewed-only敏感性家族的可估计分母已效果前固定，不能按实际P删除失败或增加来源。

[来源合同](../methods/NSCLC-candidate-external-source-association-v1.md)与私有d4141c94保持。完整双方十QC、两个父比例和全部分类块；至少20单位/df8，原rank cutoff及machine-rank一致性、最小相对奇异值1e-8。源病例响应已原方向定向，不再乘一次方向。19999固定置换/2000bootstrap、1900有效区间、全部留一，种子namespace和全部敏感性不变。

观察P为(count+1)/(19999+1)。每家族m个可估计位置的真实置换尾概率同时99%二项上界为BetaQuantile(1-.01/m,count+1,19999-count)，count=19999取1；对该上界做完整家族BH，同时保留原P的BH/BY。此处明确真实尾概率上界与plus-one观测P，按已冻结外部源码执行，不改历史native计算。

## 实际运行

数学包装源码8b0ff783468d19391b4d804e11dcfb91737e05444c4b4bb1d480566d3a48dbbd先上传，复用已接受rank内核和完整span guard。主源码f57c0c7a938c47cb9e8e20d6ed2306bbc74deca3008f5320986e9551a984b396于2026-10-18 09:37:05.291689 UTC /17:37:05.291689 UTC+8真实启动，两个单BLAS fork workers、64MiB产物、准确quota500GB+1GiB保留。完整结果及审计尚待，不能只看部分来源P或用敏感性救主家族。

独立审核audit_NSCLC_candidate_external_associations_v1.py绑定原启动句柄与正式终态，使用既有accepted OLS oracle加同一full-span guard，完整重放所有真实置换/bootstrap/LOO、全部敏感性及各家族BH/BY/Monte Carlo上界。I/O与整数exact-rank原语共享，但projection/random流/校正独立重算；不是独立生物学队列。审核实际启动/正式接受后续另记。

原不可估计位置、所有负向结果/缺失区间及病例审阅情况保留。效果/数据表/随机产物只server/授权对话，公共Git仅方法/代码/过程和SHA；当前没有直接通讯、临床效益或goal完成声明。

## 全数值正式接受与失败修复（2026-10-18）

生产09:37:50.277735 UTC /17:37:50.277735 UTC+8结束，终态605099dd5741720df181a51a8d4fc976d525ce9d7422f29e069ac39dfea40d67。审核v1于09:38:34启动，在首个效果拟合前因既有oracle无fit入口失败；历史源码e11ea011c848812e5b61a197f67ae2ed240962b044d2ee340a502cf97d97a16d和日志保留。v2仅改为既有calculate入口及说明，未改主要生产、设计、阈值或随机流，也未重新生产主结果。

v2源码af1f9bdfc45ec20a596fca015cdfdd57ee3cdae39506d6849350c6be27d73d85于09:40:46.415359 UTC /17:40:46.415359 UTC+8启动，09:42:03.875279 UTC /17:42:03.875279 UTC+8正式通过。证明dc9e287dea77c2044f03f4b7087c3cc1462436f2d913088315e33736dd0a818b，覆盖全部22可估计模型/1040计划、439978置换/44000bootstrap/736留一；最大effect差5.56e-16、bootstrap差4.27e-13。全部敏感性、BH/BY及Monte Carlo家族上界独立重算；不是独立生物学验证。

## 全结果与实际图像审阅

review源码39054a614edfc7f6c5ac744af4059a30bbf4e90b0c57d729995de300f003225f于09:48:04.193060 UTC /17:48:04.193060 UTC+8启动，09:48:06.700107 UTC /17:48:06.700107 UTC+8完成1020aedd1e548f9b464460276b294e771edd9d92688edca83320c20e5b2bbfbb。全部22结果/失败摘要和13主模型森林图已实际审阅；整图无裁切，材料/时点标签和无有效CI的空心点可读。森林图aa825eeae4f2d925405da4886565136ae720a2c377d1277132e7ada64cc07b04；私有实际审阅凭证15b0bd9983edf88e4689bc0cf0d91d411c98265028d574a45ecfbd4abb4e0b9c，精确登记时间见凭证。实际像素查看发生于09:57:43之前，未记录精确时刻，不虚构事件时间。

family summary的distinct_estimated_GSE实际是dataset标签去重，包含非GSE collection；不可作为独立GSE数引用。main与reviewed-only相同case集合不能算重复验证。原四候选全部保留，包括未形成合格模型、负向及无CI结果；科学效应和新的生物学优先级仅服务器/授权对话，未公开。下一先核实真实发送细胞来源及接收功能，不以RNA关联声明直接通讯或疗效。

## 结果后全节点组成诊断

另见[组成方法](../methods/NSCLC-candidate-external-node-composition-v1.md)。所有22原模型的精确case/library集合均保留，追加双方全部冻结节点比例，包含常数与零列；原全控制及rank/df门槛不变。09:51:43.416775 UTC /17:51:43.416775 UTC+8启动，09:51:46.711279 UTC /17:51:46.711279 UTC+8完成84549504fb5563b1ba744acac7728a59a77b5a3c1b22b5c303d43d0967025fe7；观察09:57至09:59UTC。源码ec127c35201dda65b37b427a8e5dfd81db872c75349fccd9d69ecdf900fd0b0f。全部6760池分母与metadata细节点数逐项一致，原效果复算差<1e-10；6可估计、16完整span剩余df不足均保留。没有新P/CI或主分类变化，没有独立全raw重放。冻结节点是亚型/marker分区，尚非正式状态标签；效应减弱不自动等于技术伪影。
