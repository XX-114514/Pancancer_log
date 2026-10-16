# 肺癌原生MOFA输入与发现

## 真实起点与阶段

上一goal轮为实际进展：完整头颈癌临床收束、新肺癌资格及逐cell metadata重数/固定划分已通过。云main412d16cf4ad7d8761181eaa8a739b82f4545413c于2026-10-16 17:28:46.642257 UTC独立核验。完整目标仍未达成，不能将输入或模型执行替代有价值关系。[方法合同](../methods/NSCLC-native-multicellular-factor-analysis-v1.md)、[资格](20261016_NSCLC_native_multicell_discovery_preflight_v1.md)。

本轮重核真实资格终态、原缓存和软件依赖，选择适配现有原计数的官方样本级MOFA。原一次只读环境检查误把DIALOGUE兄弟目录当当前子目录，FileNotFoundError发生在列目录时，PyPI依赖查询已返回；按现有启动代码定位真实目录，edgeR4.0.16可用，未改任何环境或数据。

## 输入和软件验收

2026-10-16 17:32:09.303016 UTC /2026-10-17 01:32:09.303016 UTC+8启动私有软件安装/合成测试，源618e14371a10f2b581692970905eca9a2d4e7234b34b7531db6db10733a6fbc9；17:35:57.569235完成。固定官方wheel/no-deps，旧环境版本一致，证明9ce1c4267c7b785ee09ebd1f7c3b0361af4116d380c874e1bb1e3a6c3bb468c8。网络等待原句柄/proc状态与ticks多次核实，未因观察超时重启。

17:33:49.481603 UTC输入准备开始，源6afc0c4541c03cbba91f38b7cfbddf4f078a992fec7662a8d23e14a6d00169d8；完整模型/留出合同在新表达模型前写出，SHA2b0b30efc6f22469c96116766b354a614a044337988596c2c351f24f0d0229d2。17:34:56.081928四视图输入完成，原计数/TMM/全部特征/QC/选择和真实R源码各SHA保存。

17:35:45.804425开始独立审计，源1356edbf069d0ac08bdbf3376585f32f15f31c645bc58f914af0877846cca6b6；17:35:57.730803 UTC /01:35:57.730803 UTC+8全部832视图行/31831原生特征/TMM/QR/训练专属选择通过，证明3b6b8086841151de638b8bb6f7f1373fc856debeae2742c392562c1c00607b8d。继承原缓存64单位原count复算范围，不夸称本阶段重扫全部原始细胞。

## 拟合接入失败与标识修正

17:39:10.914995 UTC原发现v1派发，代码54e5d9615f5a5f4e57c7f1b20b055f908dadffb0ae6941df41da4e2fe20ef71b。官方set_data_matrix拒绝不同视图的重复gene feature名称，实际在build/fit之前退出，无真实模型参数；伴随17:39:59.965915审计因生产终止退出。17:40 UTC原/proc及日志已核，不能当等待或阴性科学结果。

新v2仅添加view::gene模型特征命名空间，不改变原基因/矩阵/顺序/病例/参数/种子。原失败plan绑定新v2，旧目录与日志保留。新v2实际启动/终态及完整后验/自身投影审计以私有收据另核，尚未发布留出或临床结论。

## 执行、产物与限制

scanpy Python/OMP及BLAS单线程；PANCANCER_PROJECT_ROOT/PANCANCER_PREFLIGHT_ROOT。软件prepare_NSCLC_MFA_runtime_v1.py，输入prepare_NSCLC_native_MFA_inputs_v1.py及audit_NSCLC_native_MFA_inputs_v1.py，投影NSCLC_native_MFA_math_v1.py，发现run_NSCLC_native_MFA_discovery_v1.py/v2.py及独立audit_NSCLC_native_MFA_discovery_v1.py/v2.py。

版本目录NSCLC_MFA_runtime_20261016_v1、NSCLC_native_MFA_inputs_20261016_v1、NSCLC_native_MFA_discovery_20261016_v1/v2；软件/输入/模型各自终态与SHA、PID/startticks绑定。当前准备完整60留出检验，不能用共享all-view后验作为不同视图分数。真实输出和学习列表server private；公开无患者/科学表/权重。原CM及身份、旧阴性家族/失败保持，goal active。

## 17:46 UTC记录：保存修复与三个完整模型

v2实际17:40:55.536003 UTC /2026-10-17 01:40:55.536003 UTC+8开始，源af87d1838cab2e5fdc08464e8306ebb0475e6fbafb5a9bcf40a1fbcf80d6dc1e。主模型真实Converged，但官方saveParameters对True选项明确deprecated退出，已有部分HDF5权重保留。原job错误输出为空不是仍在运行；通过原/proc缺失及脱敏官方软件日志确认实际终止，未读候选关联。v2审计因生产终止退出。

v3只关闭废弃保存参数选项，W/Z/Tau显式期望、输入/参数/种子保持；17:43:52.055268 UTC启动，源2eee43e0937fa2faa71314df63f037e1ffe10deec90765f76247214147d0af51。三个固定种子17:44:41.408812 UTC完整结束，终态5c35ae3dd7cc6d0b0489be97aac52f94ae9a7903572f7f88235febd2d1c62379。独立v3审计17:44:57.340274开始，源46753e6749056942fbe30443d9c3362503c2ffd078df97d8c02416cdfe079b7d；全部后验/自身投影及v2主参数比对终态待核，不提前称PASS。

完整60槽留出脚本及固定统计helper已准备；helper891024dd8b8d266ea6758259ebce0fee7e6811939f8b3025988cac28a15e8973。只有模型审计PASS和主模型非最大迭代退出才可计算；尚未审阅新关联或临床结局。原全部失败和API检查历史保留，无新表达下载。

## 17:49 UTC取证：完整模型及留出验收

v3审计最初将Tau的view/group层级当直接dataset而TypeError，按实际只读结构修正v4，模型不重跑。v4源b417c9d856a45dddb6d571c801361c9369c8d87c18979073041a2fca0fb00c46，17:47:03.060112 UTC启动、17:47:04.441776 UTC通过全部三个官方模型、24960自身分数、20发现内稳定匹配及保留v2主参数逐元素一致，证明3ebc1087f7d3030c2d6b36caa96f2869a4bcdfb44e5deebbe03651276cb9ad3f。主模型正常收敛，非最大迭代退出；QR投影最大差4.63e-14以内。旧未通过审计导致的门控留出v1退出保留，无P生成。

门控留出v2仅改等待的审计句柄/新输出目录，原全部统计合同不变。17:47:41.333080 UTC /2026-10-17 01:47:41.333080 UTC+8开始，源a80f48970d5143bd1cd35c398cef0bd45b484d8d3aac4368b89040f4011a4ba4；17:48:00.085939完整60槽结束，终态0c3b7f7cf6c8cf5af042295fc7d904f0bdd1c28b93eb3f064cb833c327e4befb。17:48:36.972941独立审计开始，源3e3f130b8f3853959306c0dc78ec010cbd0050155d5ac0648e345b64cf23841c；17:48:44.784824 UTC全部1199940置换/120000bootstrap及留一、完整60家族BH/BY/MCupper通过，证明fcc4a9144139b2ac247b93443dd557fab9381d5f75733b1a7340aba9ed72a081。

完整效应与因子/细胞群支持仅server/获准当前对话审阅，不入public。下一全10因子/四视图的功能和作者亚型混合解释，保留共同刺激/技术背景可能性；留出关联不是直接LR、临床预测或独立患者验证。原CM与完整目标保持、goal active。

## 17:55 UTC完整审阅和下一解释问题

17:53:56.773701 UTC派发全部60关系/10因子汇总及全图，源eed7066b0fcf138af3d42a34c87d91ef1dfc710390cf4fde21067a80ad93d871；17:53:58.100267完成，终态4da97d15671379bd3cd7a04045f05cafe665d1a133b30921e43c9152b490e08d。图像全部格数/标签/BY标记在17:55:10.424990 UTC /2026-10-17 01:55:10.424990 UTC+8实际审阅，回执b16d4e60ef058bee95a1b8547bfedd01c3f0527e0e2abbabe674e654330e1276。无新P或挑选子家族。

下一完整问题PRIVATE_NEXT_ALL_FACTOR_FUNCTION_AND_MIXTURE_SCOPE.json SHA720e394d8a272dd998b88e734b37aeafc929498aab93880f6c968fd7d470bf2c：全部10×四视图×50原既有功能、未选基因RNA及全部作者亚型/技术/人口学诊断，描述性、不新增P/功能FDR。完整60组成敏感性保留原病例/分层。未开始这些新测量或临床关联；不将载荷注释当独立通路验证，不公开学习基因列表。所有失败、原数据和目标保持。
