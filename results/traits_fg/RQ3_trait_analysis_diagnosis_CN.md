# 第三个科学问题分析诊断清单（中文）

## 适用范围

本清单针对当前工作区中与第三个科学问题（RQ3: traits/functional strategies as winners and losers under land-use x warming）相关的三类分析：

1. 单性状三重交互模型
2. 多性状联合两重交互模型
3. 基于 PCA/聚类的功能组（FG）模型

核心判断基于以下文稿与代码/结果：

- 研究计划中 RQ3 的定义：`P(occurrence) ~ LandUse x TempAnom x Traits/FG`
- 单性状脚本：`/Users/dingchenchen/lucc/性状交互作用_出现概率.R`
- 单性状建模脚本：`/Users/dingchenchen/lucc/性状预测.R`
- 多性状联合模型脚本：`/Users/dingchenchen/lucc/筛选性状交互作用_出现概率.R`
- FG 模型脚本：`/Users/dingchenchen/lucc/性状交互作用_FG_Gower1.R`
- 代表性结果目录：
  - `/Users/dingchenchen/lucc/02_ThreeWay_ABS_and_PCT1`
  - `/Users/dingchenchen/lucc/glmmTMB_occurrence_alltraits_AIC_2way_ABS_PCT_medianControls`
  - `/Users/dingchenchen/lucc/glmmTMB_occurrence_FG_GowerClustering_FULL_plusTwoWay`
  - `/Users/dingchenchen/lucc/glmmTMB_occurrence_FG_GowerWinsor_FULL_plusTwoWay`
  - `/Users/dingchenchen/lucc/brms_occurrence_FG3_PCA_rawTraits_FULL_plusTwoWay`

## 一、总体结论

当前第三个科学问题的分析框架已经比较完整，但三类方法的定位应明确区分：

- 单性状三重交互最贴近 RQ3 的原始科学问题，应该作为主分析框架。
- 多性状联合模型更适合作为“控制其他性状后的独立信号筛查”或稳健性分析，不宜替代主分析。
- FG 分析有潜力做机制综合，但只有“分组稳定”的 FG 方案适合作为主文结果；当前 raw Gower 两组版本不建议主用。

## 二、建议保留为主结果

### 1. 单性状三重交互：保留为主分析

原因：

- 与 proposal 中 RQ3 的公式最一致：`Occurrence ~ UI2 x warming x trait`
- 最容易直接回答“哪类性状在何种土地利用下随变暖更敏感/更耐受”
- 解释层次清楚，适合和 hypothesis 一一对应

当前已确认的本地模型例子：

- `/Users/dingchenchen/lucc/brm_OccmodelRSTmean_ub.rds`

该模型的公式为：

- `Occur ~ UI2 * StdTmeanAnomalyRS * RS.rs + (1 || SS) + (1 || SSBS) + (1 || Best_guess_binomial)`

本地核实到的 RS 模型使用了：

- 444,921 条记录
- 2,958 个物种

优点：

- 使用全物种池，不受多变量 complete-case 明显压缩
- brms 后验输出适合直接报告 95% CrI
- 对论文主叙事最友好

### 2. PCA 三组 FG：可作为主文或主文补充中的机制综合结果

推荐对象：

- `/Users/dingchenchen/lucc/brms_occurrence_FG3_PCA_rawTraits_FULL_plusTwoWay`

理由：

- 分组规模相对平衡：FG1 = 1898, FG2 = 514, FG3 = 433 个物种
- 记录层面占比也更合理：
  - FG1 = 216,947
  - FG2 = 66,130
  - FG3 = 149,483
- 性状解释较稳定，能够形成可叙述的生态策略组

对应解释文件：

- `/Users/dingchenchen/lucc/brms_occurrence_FG3_PCA_rawTraits_FULL_plusTwoWay/05_Tables_FG_Interpretation/FG_trait_signatures_and_interpretation.txt`

当前三组概括如下：

- FG1：低扩散、低热耐受、较小范围，偏专性
- FG2：高扩散、较长世代时间、体型偏大
- FG3：分布广、热泛化、高繁殖输出

这类结果比 raw Gower 两组更符合“功能策略组”的生态解释。

## 三、建议作为补充或稳健性分析

### 3. 多性状联合两重交互：适合作为补充，不宜作为主检验

对应脚本：

- `/Users/dingchenchen/lucc/筛选性状交互作用_出现概率.R`

对应结果目录：

- `/Users/dingchenchen/lucc/glmmTMB_occurrence_alltraits_AIC_2way_ABS_PCT_medianControls`

模型公式在脚本中约第 261-266 行：

- `Occur ~ UI2 * StdTmeanAnomalyRS + UI2 * (all traits) + StdTmeanAnomalyRS * (all traits) + random effects`

关键问题：

- 这个模型没有 `UI2 x warming x trait` 三重交互
- 因此它不能直接回答 RQ3 中“trait mediation of the land-use x warming interaction”
- 它更适合回答：
  - 控制其他性状后，哪些 trait 仍然与 UI2 或 warming 相关？
  - 多个 trait 同时进入模型后，单个 trait 的独立贡献是否还存在？

优势：

- 能减少单性状逐个分析造成的“表面显著”
- 能帮助识别哪些 trait 信号在联合模型中仍然稳健

局限：

- 不对应主 hypothesis 的三重交互检验
- 容易被误读为“已经完成 trait mediation 主检验”

建议定位：

- Supplementary robustness analysis
- 或正文一句话说明“we further fitted a multivariate two-way interaction model as a conditional-independence check”

### 4. Winsorized Gower 两组 FG：可作为敏感性分析

对应结果目录：

- `/Users/dingchenchen/lucc/glmmTMB_occurrence_FG_GowerWinsor_FULL_plusTwoWay`

其优点是：

- 对极端值做 winsorization
- 明确加入最小组大小约束
- 得到的分组比 raw Gower 两组更稳定

当前物种数：

- FG1 = 2554
- FG2 = 291

记录占比：

- FG1 = 386,964
- FG2 = 45,596

对应分组选择文件：

- `/Users/dingchenchen/lucc/glmmTMB_occurrence_FG_GowerWinsor_FULL_plusTwoWay/05_Tables_FG_Interpretation/FG_GowerWinsor_method_K_silhouette_min5pct.csv`

可作为说明：

- 当限制极端值并避免超小 cluster 后，FG 结果仍然保留“慢生活史/高体重/高扩散”与“其余物种”之间的区别

但由于仍然只有两组，生态解释比 PCA 三组略粗。

## 四、不建议作为主结果直接使用

### 5. Raw Gower 两组 FG：不建议作为主文结果

对应结果目录：

- `/Users/dingchenchen/lucc/glmmTMB_occurrence_FG_GowerClustering_FULL_plusTwoWay`

对应选择日志：

- `/Users/dingchenchen/lucc/glmmTMB_occurrence_FG_GowerClustering_FULL_plusTwoWay/05_Tables_FG_Interpretation/FG_trait_signatures_and_selection_log.txt`

该版本最严重的问题是分组塌缩：

- FG1 = 2835 个物种
- FG2 = 10 个物种

记录层面也极不平衡：

- FG1 = 427,688
- FG2 = 4,872

FG2 物种包括：

- `Egretta garzetta`
- `Falco peregrinus`
- `Himantopus himantopus`
- `Milvus migrans`
- `Pandion haliaetus`
- `Larus argentatus`

这表明该 FG2 更像“极端性状或离群物种集合”，而不是一个稳定的功能策略组。

原因在代码层面也能定位：

- `/Users/dingchenchen/lucc/性状交互作用_FG_Gower1.R` 第 200-245 行按 silhouette 选 K 和方法
- 但这里没有最小 cluster size 约束

因此：

- silhouette 选出的“最佳两组”不一定生态上合理
- 会把极端值聚成很小的一组
- 后续三重交互显著性很可能部分来自极端小组与主体组的对比，而非稳定策略差异

结论：

- raw Gower 两组结果不宜作为主文证据
- 最多可在方法尝试或敏感性分析中简短提及

## 五、当前代码和结果之间最关键的不一致

### 6. 单性状脚本与实际模型版本不完全一致

问题 1：模型路径仍指向外部 Windows 路径

- `/Users/dingchenchen/lucc/性状交互作用_出现概率.R` 第 789-798 行

这里读取的是：

- `D:/DCC/BIRDLIFE1/...`

这说明当前脚本不是一个完全可在本工作区复现的最终版本。

问题 2：Tmin/Tmax 位置变量是否标准化，脚本之间不一致

在 plotting 脚本中：

- `/Users/dingchenchen/lucc/性状交互作用_出现概率.R` 第 18-25 行写的是 raw 使用

但在建模脚本中：

- `/Users/dingchenchen/lucc/性状预测.R` 第 1227-1239 行实际拟合的是 `Tmin_position.rs` 和 `Tmax_position.rs`

这意味着：

- 图题、文字说明和真实模型尺度可能不一致
- 后续若把结果写进论文，容易被审稿人指出变量定义混乱

问题 3：不同脚本版本混有 BM 结果，但当前 plotting 脚本不再包含 BM

- `/Users/dingchenchen/lucc/02_ThreeWay_ABS_and_PCT1` 中仍有 BM 图
- 但 `/Users/dingchenchen/lucc/性状交互作用_出现概率.R` 当前显式跑的是 RS/HB/TR/HWI/GL/CS/Tmin/Tmax

说明图和脚本来自不同阶段的版本，需统一。

## 六、样本不一致导致的比较风险

不同方法使用的记录/物种集合并不一样：

- 单性状 RS 模型：444,921 记录，2,958 物种
- 多性状联合模型：432,412 记录，2,845 物种
- FG raw-trait 模型：432,560 记录，2,845 物种

缺失主要来自：

- `TR`
- `Tmin_position`
- `Tmax_position`

因此不能直接把不同方法图上的效应强弱做横向比较，否则会混入 species pool 改变带来的差异。

建议：

- 若要正式比较三类方法，最好统一到同一套 2,845 物种 complete-case 数据
- 或者在文中明确说明各分析的样本池不同

## 七、关于 %Δp 图的解释风险

三类分析的 baseline 不一致：

- 单性状三重交互：baseline = PV + temp=0 + trait Q1
- 多性状联合模型：baseline = PV + temp=0 + all traits median
- FG 模型：baseline = PV + temp=0 + FG_base

因此：

- `%Δp` 只能在同一分析内部比较
- 不能跨方法比较绝对大小
- 如果主文需要跨方法综合，建议改用绝对概率、边际 slope 或 contrast

## 八、关于不确定性区间的技术提醒

单性状 brms 模型：

- 后验区间是真正的 posterior CrI

多性状 glmmTMB 与 FG glmmTMB：

- 目前图中的区间是基于固定效应参数协方差矩阵或 link-scale normal approximation 的近似区间

对应代码位置：

- 多性状模型：`/Users/dingchenchen/lucc/筛选性状交互作用_出现概率.R` 第 302-318 行
- FG glmmTMB：`/Users/dingchenchen/lucc/性状交互作用_FG_Gower1.R` 第 412-420 行

因此：

- 可以用于探索和展示趋势
- 但不应与 brms 后验 CrI 当作同等级证据描述

## 九、最终建议排序

### 建议主文保留

1. 单性状三重交互 brms
2. PCA 三组 FG 三重交互

### 建议作为补充或稳健性分析

1. 多性状联合两重交互 glmmTMB
2. Winsorized Gower 两组 FG

### 建议不作为主结果

1. Raw Gower 两组 FG
2. 当前无法完全追溯脚本版本来源的旧图结果

## 十、下一步最值得做的整理

如果准备写论文或做正式汇报，最建议优先完成以下 5 件事：

1. 统一单性状脚本、模型对象和导出图的版本，消除外部路径依赖。
2. 明确 `Tmin_position/Tmax_position` 到底使用 raw 还是标准化版本，并全程统一。
3. 把多性状联合模型重新定义为“补充稳健性分析”，不要与主检验混写。
4. 以 PCA 三组 FG 为主，winsorized Gower 两组为 sensitivity，放弃 raw Gower 两组主解释。
5. 若要跨方法严格比较，统一到同一个 complete-case 物种池再出一版核心图。

## 一句话总结

目前第三个科学问题最可信的主线是：

“用单性状三重交互回答 trait-mediated vulnerability，用稳定的 FG 分组回答机制综合；多性状联合模型只作为控制其他性状后的补充检验；raw Gower 两组结果不宜主用。”
