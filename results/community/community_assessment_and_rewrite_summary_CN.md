# 群落多样性、群落结构与功能群落结构：V2 审查与重写总结

## 一、这一版的定位

这一版不是把你当前目录里的社区层面分析全部推翻重做，而是在保留已有主结果逻辑的基础上，把最需要修正的环节系统收束成一套更适合论文主文、答辩汇报和后续维护的“主版本”。

当前建议的主版本目录为：

- 群落多样性：
  `/Users/dingchenchen/lucc/community_rewrite_v2/diversity_main_nb`
- 群落组成结构：
  `/Users/dingchenchen/lucc/community_rewrite_v2/compositional_main_bray`
- 功能群落结构：
  `/Users/dingchenchen/lucc/community_rewrite_v2/functional_main_bray`

## 二、这次复核后确认的核心问题

### 1. 群落多样性

原主线脚本：
`/Users/dingchenchen/lucc/6_RunLUClimateModel_new.R`

主要问题：

- abundance 用 Gaussian mixed model 基本可以继续保留。
- richness 原先主线依赖 `Poisson GLMM`，而旧诊断已经显示明显 dispersion 问题。
- `GLMERSelect + poly(x,1)` 的写法不利于复核，也不利于答辩中清楚解释。

V2 改进：

- abundance 保留为显式 `lmer`
- richness 改为显式 `glmmTMB::nbinom2`
- 增加统一样本支撑图、系数森林图、预测表和三格式主图

当前 V2 输出见：
`/Users/dingchenchen/lucc/community_rewrite_v2/diversity_main_nb`

### 2. 群落组成结构

原主线脚本：
`/Users/dingchenchen/lucc/群落同质性_dissimilarity1.R`

主要问题：

- 原脚本把 `lmer`、`brms Gaussian(logit D)`、`brms Beta(D)` 以及 9 个指标并列展开，主结果层级过重。
- 从写作和答辩的角度，没有必要让 Bray、Jaccard、Sorensen 三整套指标平行占据主叙事。
- 各土地利用类型的 `pair_clim_mean` 经验支撑范围差异明显，尤其 Urban 更窄，因此预测不能脱离支持区间解释。

V2 改进：

- 主结果只保留 Bray total / turnover / gradient
- 统一保留固定 `PV@0` 对照表达
- 补上 climate support 图、主预测图、差异图和 warming-related coefficient forest

当前 V2 输出见：
`/Users/dingchenchen/lucc/community_rewrite_v2/compositional_main_bray`

### 3. 功能群落结构

原主线脚本：
`/Users/dingchenchen/lucc/功能同质化.R`

主要问题：

- 连续 ATA 空间上的 Bray 分析是合理的主线。
- 但功能版 Jaccard / Sorensen 依赖 `ATA > 0` 的二值化 presence/absence，当前数据中该指标接近退化。
- 这一点已经在数据上确认：
  `site_trait_matrix_ATA01.rds` 中站点层 `rowMeans(ATA > 0)` 的中位数为 `1`。

V2 改进：

- 主叙事仅保留 ATA-space Bray total / turnover / gradient
- 功能 presence/absence 家族不再进入主图
- 补上支撑图、主预测图、PV@0 差异图、warming coefficient forest 和退化诊断表

当前 V2 输出见：
`/Users/dingchenchen/lucc/community_rewrite_v2/functional_main_bray`

## 三、这版重写后的推荐主文结构

### 主文建议保留

1. 群落多样性
   abundance + richness(NB)

2. 群落组成结构
   Bray total / turnover / gradient

3. 功能群落结构
   functional Bray total / turnover / gradient on ATA space

### 建议放补充材料

- compositional Jaccard / Sorensen
- functional Jaccard / Sorensen
- 并列的多套 Bayesian robustness（如后续需要）
- 旧版 GLMERSelect richness 主结果

## 四、这版结果解释上的优点

- 主线更稳定：避免把方法上已不够稳的指标继续当主结论。
- 层级更清晰：主结果与 robustness 不再混放。
- 解释更直接：统一采用更容易讲清楚的固定基线表达。
- 写作更省力：每张主图都已配好中英双语简要描述。

## 五、仍然需要保留的谨慎点

### 1. diversity richness NB

V2 的 NB richness 明显优于 Poisson 主线，但 DHARMa 仍提示 dispersion departure，因此建议写作中避免把它表述成“完全无诊断问题”，更稳妥的说法是：

- NB 主模型显著缓解了原 Poisson richness 的不严谨性
- 当前主结果应以 NB richness 而不是 Poisson richness 为准

### 2. climate support

无论是 compositional 还是 functional，Urban 的气候支持都更窄，因此图中曲线要明确解释为：

- 在经验支持区间内的模型预测
- 不能把边缘端点理解为对全气候空间的强外推

### 3. functional presence/absence

功能版 Jaccard / Sorensen 不是“数学上不能算”，而是当前这套 ATA 数据结构下生态解释已经明显变弱，因此不建议再让它们进入主叙事。

## 六、代码与结果组织方式

这一版统一采用：

- `01_DataSupport`
- `02_Models`
- `03_Plots`
- `04_Tables`
- `05_Figure_Captions`

并保持三格式出图：

- `.png`
- `.pdf`
- `.pptx`

## 七、结论

这次重写之后，社区层面分析的主线已经可以稳定收束为：

1. `diversity_main_nb`
2. `compositional_main_bray`
3. `functional_main_bray`

如果后面继续做答辩版 PPT，我建议就按这三条主线直接组织，不要再把旧版脚本或全部补充指标重新并列拉回主页面。这样最专业、最稳定，也最容易讲清楚。
