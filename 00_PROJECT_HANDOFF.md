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

### 2.1 当前 GR-GWR 算法不是最终算法

当前 `src/georegime_gwr/grgwr.py` 是从 pyGWRx 当前思路抽出的 **baseline snapshot**。后续完全允许：

- 删除 ICM；
- 替换 Ward；
- 替换 kNN + MST；
- 直接使用标准空间邻接矩阵 W；
- 改写边界惩罚；
- 改变 regime 数量选择；
- 改变分区特征；
- 改变 bandwidth policy；
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

## 3. 标准 GWR 基线：端到端验证已通过

当前阶段**不研究 MGWR**。虽然外部参考来自 Python 包 `mgwr`，但只使用标准 GWR 的：

- `mgwr.gwr.GWR`
- `mgwr.sel_bw.Sel_BW`

它们只作为外部验证 oracle。仓库自己的 `BasicGWR` 已经能够独立完成自动带宽搜索和 GWR 拟合。

### 3.1 Canonical Georgia benchmark

固定规格：

- 159 Georgia counties，1990 Census；
- response: `PctBach`；
- predictors: `PctFB`, `PctBlack`, `PctRural`；
- projected coordinates: `X`, `Y`；
- X 与 y 均按 `ddof=0` 做 z-score；
- adaptive bisquare；
- AICc bandwidth selection；
- standard GWR only。

Canonical `mgwr==2.2.1` 结果：

- bandwidth = 117；
- RSS = 51.186192；
- ENP/trace(S) = 11.804770；
- AIC = 296.615923；
- AICc = 299.050809；
- BIC = 335.912535；
- R2 = 0.678074；
- adjusted R2 = 0.652080。

特别注意：`mgwr.GWRResults.RSS` 是 location-wise weighted RSS diagnostic；本项目中的模型级总 RSS 固定使用 `sum(resid_response**2)`。

### 3.2 Repository BasicGWR end-to-end validation — PASSED

现在验证不是“手工给 bandwidth=117 再拟合”，而是：

`BasicGWR(bandwidth="auto")`

自己搜索带宽，再拟合。

结果：

- external `mgwr.sel_bw.Sel_BW` bandwidth = `117`；
- repository BasicGWR bandwidth = `117`；
- bandwidth 完全一致；
- max absolute parameter difference = `5.551115123125783e-16`；
- RMSE parameter difference = `9.131398707249937e-17`；
- max absolute fitted difference = `5.551115123125783e-16`；
- max absolute residual difference = `5.551115123125783e-16`；
- max absolute hat-matrix difference = `5.551115123125783e-17`；
- BasicGWR RSS = `51.186191553463985`；
- mgwr.GWR RSS = `51.18619155346399`；
- trace(S) 两者均为 `11.804769716730094`；
- BasicGWR AICc = `299.0508086830287`；
- mgwr.GWR AICc = `299.0508086830288`；
- tolerance `1e-12` end-to-end validation = **PASS**。

因此标准 GWR 的**自动选带宽 + 最终局部回归**都已经可信。

关键文件：

- `src/georegime_gwr/gwr.py`
- `src/georegime_gwr/bandwidth.py`
- `data/raw/georgia/GData_utm.csv`
- `data/raw/georgia/G_utm.*`
- `experiments/real_data/georgia_gwr_baseline/reproduce.py`
- `experiments/validation/basicgwr_vs_mgwr_georgia.py`
- `results/reproduction/georgia_gwr_baseline/summary.json`
- `results/validation/basicgwr_vs_mgwr_georgia/summary.json`
- `results/validation/basicgwr_vs_mgwr_georgia/pointwise_comparison.csv`
- `results/validation/basicgwr_vs_mgwr_georgia/basicgwr_bandwidth_search_curve.csv`
- `.github/workflows/basicgwr-validation.yml`

## 4. 极重要：两种 bandwidth search policy 不同

后续新对话必须保留这个事实，不得重新误判。

### 4.1 Current PyGWRx policy

当前 pyGWRx 的 adaptive AICc selector 对合法整数 bandwidth 做**穷举扫描**，取完整候选曲线中的最低 AICc。

在 canonical Georgia 数据上：

- exhaustive PyGWRx-style global minimum = `k=116`；
- AICc 约 `298.9856`。

### 4.2 Canonical mgwr 2.2.1 policy

`mgwr==2.2.1` 默认使用**离散 golden-section search**，不是完整整数穷举。

在同一 Georgia 数据上返回：

- `k=117`；
- AICc = `299.0508086830288`。

所以 `116` 与 `117` 的差异来自**搜索算法**，不是 GWR 回归核心实现错误。

### 4.3 Repository decision

仓库同时保留：

- `PyGWRxAdaptiveAICcSelector`：忠实保留 current PyGWRx exhaustive policy；
- `MGWRCompatibleAICcSelector`：本地实现 canonical `mgwr==2.2.1` discrete golden-section policy。

当前：

`BasicGWR(bandwidth="auto")`

使用 `MGWRCompatibleAICcSelector`，因为当前 canonical baseline 的目标是**完整复现论文/官方标准 GWR 实验，包括自动选择 117**。

正式决策：

`docs/decisions/ADR-0003-standard-gwr-bandwidth-search-policy.md`

后续若 GR-GWR 论文主实验选择 exhaustive search，应把它作为新的方法选择明确记录，不能悄悄把 116/117 混为一谈。

## 5. 当前 GR-GWR baseline 流程

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

## 6. 当前高优先级研究问题

### Q1. 空间邻接结构 W
当前代码固定从点坐标构造 `kNN + MST`。对于 polygon 数据，理论上 GR-GWR 更适合直接定义标准空间邻接矩阵 `W`，polygon 可用 Queen/Rook；point data 可用 kNN/Delaunay/distance graph。

### Q2. 坐标是否应进入聚类特征
当前同时存在 coordinates-in-features 与 adjacency/connectivity 约束，可能重复，需要消融。

### Q3. ICM 是否必要
需要比较无 refinement 与 ICM refinement；若收益很小，倾向简化。

### Q4. 第一轮 GWR 系数是否被真实边界平滑污染
必须用 sharp-boundary synthetic DGP 验证。

### Q5. shared regime 假设
不同 covariates 可能具有不同边界，需要 shared-boundary 与 conflicting-boundary DGP。

### Q6. graph-cut penalty 的解释
`B(z) = number of adjacency edges crossing regimes` 是 graph-cut penalty，不是真实几何边界长度。

### Q7. GR-GWR 最终 bandwidth policy
标准 GWR canonical reproduction 已固定为 mgwr-compatible golden search，但 GR-GWR 正式论文算法究竟采用 benchmark-compatible search 还是 exhaustive AICc global minimum，仍未冻结。

## 7. 文献定位（当前阶段结论）

已经确认：

- spatial regimes 不是首创；
- GWR coefficient clustering 不是首创；
- “区内平滑、区间突变”不是首创；
- boundary-aware GWR / multiscale spatially varying coefficient 相关研究已有先例；
- 当前尚未找到与“GWR local relationships -> endogenous contiguous regime -> regime-gated GWR -> LOO-guided iterative regime refinement”完全同构的方法。

首创性声明必须谨慎，例如 `To the best of our knowledge`，并限定在完整耦合框架，而不是单个已有组件。

## 8. 新对话恢复顺序

按以下顺序读取：

1. `00_PROJECT_HANDOFF.md`
2. `ARCHITECTURE_INDEX.md`
3. `docs/project/CURRENT_STATUS.md`
4. `docs/design/GRGWR_BASELINE_SPEC.md`
5. `docs/design/RESEARCH_QUESTIONS.md`
6. `docs/experiments/EXPERIMENT_PLAN.md`
7. `docs/decisions/` 最新 ADR，当前至少读 `ADR-0003-standard-gwr-bandwidth-search-policy.md`
8. 当前源码与 tests

如果文档冲突：

`最新 ADR > CURRENT_STATUS > 当前正式设计规范 > 本交接文档中的旧描述 > README`

## 9. 当前下一步

标准 GWR 已完成端到端验证，不需要继续修基础 GWR。下一阶段：

1. 建立最基础 synthetic DGP；
2. 研究是否把 `W` 从固定 kNN+MST 构造结果提升为显式模型输入；
3. 比较 Queen/Rook、kNN、kNN+MST；
4. 对 coordinates-in-features、Ward、ICM 做逐组件消融；
5. 单独决定 GR-GWR 正式 bandwidth policy；
6. 再冻结 GR-GWR v1 paper algorithm。
