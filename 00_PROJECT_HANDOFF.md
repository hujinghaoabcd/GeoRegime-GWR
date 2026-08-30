# 00_PROJECT_HANDOFF

> **跨对话最高优先级入口。任何新的 ChatGPT / Codex 对话在继续本仓库前，先读本文。**

## 1. 项目定位

本仓库 `GeoRegime-GWR` 是 **Geo-Regime Geographically Weighted Regression (GR-GWR)** 方法论文的研究仓库。

它不是成熟软件包，也不追求 pyGWRx 的复杂架构。用途是：

- 独立理解、重构和验证 GR-GWR；
- 保留一个最小基础 GWR 作为比较与内部计算引擎；
- 快速修改模型数学结构与算法实现；
- 运行仿真、消融、基准和实证实验；
- 生成论文表格、图和可复现实验结果；
- 记录所有关键设计变化，确保跨对话不丢状态。

## 2. 最重要原则

### 2.1 当前算法不是最终算法

当前 `src/georegime_gwr/grgwr.py` 是从 pyGWRx 当前思路抽出的 **baseline snapshot**。后续完全允许：

- 删除 ICM；
- 替换 Ward；
- 替换 kNN + MST；
- 直接使用标准空间邻接矩阵 W；
- 改写边界惩罚；
- 改变 regime 数量选择；
- 改变分区特征；
- 改变整个优化框架。

不要因为代码已经存在就把它当作论文最终方法。

### 2.2 研究顺序

每个算法组件都按以下顺序处理：

1. 先说明当前源码实际做什么；
2. 再说明数学形式；
3. 再检索文献判断理论依据和已有先例；
4. 再设计针对性仿真；
5. 最后决定保留、修改还是删除。

### 2.3 所有重要决定必须留痕

任何会改变论文方法定义的决定，都必须同步更新：

- `docs/project/CURRENT_STATUS.md`
- `docs/design/GRGWR_BASELINE_SPEC.md` 或后续正式模型规范
- `docs/decisions/` 中新增 ADR
- 如流程发生重大变化，更新本文

## 3. 标准 GWR 基线已经固定

当前阶段暂时**不研究 MGWR**。

虽然外部参考实现来自 Python 包 `mgwr`，但只使用：

- `mgwr.gwr.GWR`
- `mgwr.sel_bw.Sel_BW`

用途是建立可信的标准 GWR 外部基线，并验证仓库自己的 `BasicGWR`。

### 3.1 Canonical Georgia standard-GWR benchmark

固定规格：

- 159 Georgia counties，1990 Census；
- response: `PctBach`；
- predictors: `PctFB`, `PctBlack`, `PctRural`；
- projected coordinates: `X`, `Y`；
- X 与 y 均按 `ddof=0` 做 z-score；
- adaptive bisquare；
- AICc bandwidth selection；
- standard GWR only。

使用 `mgwr==2.2.1` 的 `mgwr.gwr.GWR` 已在线复现：

- bandwidth = 117；
- RSS = 51.186192；
- ENP = 11.804770；
- AIC = 296.615923；
- AICc = 299.050809；
- BIC = 335.912535；
- R2 = 0.678074；
- adjusted R2 = 0.652080。

特别注意：`mgwr.GWRResults.RSS` 是 location-wise weighted RSS diagnostic；本项目中的模型级总 RSS 固定使用 `sum(resid_response**2)`。

### 3.2 Repository BasicGWR validation — PASSED

仓库内 `src/georegime_gwr/gwr.py::BasicGWR` 已用相同 Georgia 数据、相同标准化、相同 adaptive bisquare bandwidth=117，与 `mgwr.gwr.GWR` 逐点比较。

比较范围：

- 159 × 4 local parameters；
- 159 fitted values；
- 159 residuals；
- 完整 159 × 159 hat matrix；
- RSS；
- trace(S)。

结果：

- max absolute parameter difference = `5.551115123125783e-16`；
- RMSE parameter difference = `9.131398707249937e-17`；
- max absolute fitted difference = `5.551115123125783e-16`；
- max absolute residual difference = `5.551115123125783e-16`；
- max absolute hat-matrix difference = `5.551115123125783e-17`；
- BasicGWR RSS = `51.186191553463985`；
- mgwr.GWR RSS = `51.18619155346399`；
- trace(S) 两者均为 `11.804769716730094`；
- tolerance `1e-8` = **PASS**。

结论：差异只有浮点机器精度量级。当前 `BasicGWR` 的标准 Gaussian adaptive-bisquare GWR 核心计算已经可信，可以作为后续 GR-GWR 研究的基础 GWR 引擎。

关键文件：

- `data/raw/georgia/GData_utm.csv`
- `data/raw/georgia/G_utm.*`
- `experiments/real_data/georgia_gwr_baseline/reproduce.py`
- `results/reproduction/georgia_gwr_baseline/summary.json`
- `experiments/validation/basicgwr_vs_mgwr_georgia.py`
- `results/validation/basicgwr_vs_mgwr_georgia/summary.json`
- `results/validation/basicgwr_vs_mgwr_georgia/pointwise_comparison.csv`
- `.github/workflows/basicgwr-validation.yml`

旧的 Georgia GWR/MGWR 联合复现实验已从当前主分支移除，避免后续研究混入 MGWR。

## 4. 当前 GR-GWR baseline 流程

当前 baseline 仍忠于 pyGWRx 的主要执行顺序：

1. 根据坐标建立 kNN + MST 空间邻接图；
2. 在全域样本上运行普通 GWR，得到每个位置一套局部参数；
3. 去掉截距，对斜率标准化；
4. 将标准化斜率与标准化坐标按 `spatial_constraint_weight` 组合；
5. 使用带空间连通约束的 Ward 层次聚类形成初始 regime；
6. 给定 regime 后，只使用同 regime 样本进行局部 GWR；
7. 可选 ICM 风格逐点调整 regime；
8. 以 `RSS + lambda * graph-cut boundary count` 做整轮 objective guard；
9. 最终按稳定 regime 重新拟合每个位置的局部系数。

## 5. 当前已经识别的高优先级研究问题

### Q1. 空间邻接结构 W
当前代码固定从点坐标构造 `kNN + MST`。但对于 polygon 数据，GIS 中已有 Queen/Rook contiguity；理论上 GR-GWR 更适合直接定义一个标准空间邻接矩阵 `W`，而不是把 kNN+MST 写死为模型定义。

**当前状态：标准 GWR 已验证通过，可以开始正式研究这一问题。**

### Q2. 坐标是否应进入聚类特征
当前同时存在坐标进入 clustering feature 与 adjacency/connectivity 空间约束，可能重复，需要消融验证。

### Q3. ICM 是否必要
当前 ICM 让 regime-restricted GWR 的拟合结果反向修正 regime，但会增加算法复杂度。需要比较无 refinement 与当前 ICM refinement。

### Q4. 第一轮 GWR 系数是否会被真实边界平滑污染
当前用普通 GWR 局部系数发现 regime，而普通 GWR 本身可能跨真实边界借样本。必须通过 sharp-boundary simulation 验证。

### Q5. shared regime 假设
当前所有局部斜率共同产生一套 regime map。不同 covariates 可能具有不同边界，需要 shared-boundary 与 conflicting-boundary DGP 对比。

### Q6. graph-cut penalty 的解释
当前 `B(z) = number of adjacency edges crossing regimes`。它是 graph-cut penalty，不是真实几何边界长度。

## 6. 文献定位（当前阶段结论）

已经确认：

- spatial regimes 不是首创；
- GWR coefficient clustering 不是首创；
- “区内平滑、区间突变”不是首创；
- boundary-aware GWR / multiscale spatially varying coefficient 相关研究已有先例；
- 当前尚未找到与“GWR local relationships -> endogenous contiguous regime -> regime-gated GWR -> LOO-guided iterative regime refinement”完全同构的方法。

论文中的首创性声明必须使用谨慎措辞，例如 `To the best of our knowledge`，并把创新限定在完整耦合框架，而不是单个已有组件。

## 7. 新对话恢复顺序

按以下顺序读取：

1. `00_PROJECT_HANDOFF.md`
2. `ARCHITECTURE_INDEX.md`
3. `docs/project/CURRENT_STATUS.md`
4. `docs/design/GRGWR_BASELINE_SPEC.md`
5. `docs/design/RESEARCH_QUESTIONS.md`
6. `docs/experiments/EXPERIMENT_PLAN.md`
7. `docs/decisions/` 最新 ADR
8. 当前源码与 tests

如果文档互相冲突：

`最新 ADR > CURRENT_STATUS > 当前正式设计规范 > 本交接文档中的旧描述 > README`

## 8. 当前下一步

标准 GWR 基线已经验证完成，不需要继续修 BasicGWR 核心计算。下一阶段：

1. 建立最基础 synthetic DGP；
2. 研究是否把 `W` 从固定 kNN+MST 构造结果提升为显式模型输入；
3. 对 Queen/Rook、kNN、kNN+MST 做比较；
4. 对 coordinates-in-features、Ward、ICM 做逐组件消融；
5. 再决定 GR-GWR v1 paper algorithm。
