# EXPERIMENT_PLAN

> 当前实验计划服务于“决定算法应该长什么样”，不是为了尽快堆论文结果。

## Stage A. Baseline sanity checks

目标：证明抽离后的研究代码与当前 GR-GWR 思路一致、基础数值行为正常。

- A1. Basic GWR 能恢复简单线性局部关系。
- A2. GR-GWR 在明显两区跳变 synthetic data 上能产生两个连续 regime。
- A3. `refine=False` 与 `refine=True` 都能稳定运行。
- A4. objective history 在接受的迭代上不增加。

## Stage B. Three canonical DGPs

### B1. Globally smooth coefficient field

真值没有机制边界：

`beta(s)` 全局连续平滑。

预期：普通 GWR 应当是强基线；GR-GWR 不应通过人为分区获得虚假的明显优势。

### B2. Piecewise constant regimes

每个 regime 内参数常数，跨区跳变。

预期：传统 spatial regime / clustered regression 应当很强。

### B3. Within-regime smooth + between-regime jump

每个 regime 内：

`beta_r(s)` 连续变化；

跨 regime：

`beta_r(s)` 发生 jump discontinuity。

这是 GR-GWR 的目标 DGP。

主要指标：

- coefficient RMSE / bias；
- prediction RMSE / MAE；
- regime ARI / NMI；
- boundary precision / recall / F1；
- boundary displacement；
- runtime。

## Stage C. Component ablations

### C1. ICM refinement

- GR-GWR no-refinement
- GR-GWR current ICM

判断：ICM 是否值得保留。

### C2. Spatial coordinates in clustering features

- local slopes only
- local slopes + coordinates

并固定同一 adjacency constraint。

判断：坐标是否提供独立信息，还是重复约束。

### C3. Spatial topology W

对 polygon / lattice synthetic data 比较：

- Queen
- Rook
- kNN
- kNN + MST

对 point data 比较：

- kNN
- Delaunay
- distance threshold

判断：正式模型应如何定义空间邻接。

### C4. Boundary penalty

- no penalty
- unweighted graph cut
- weighted graph cut

### C5. Hard vs soft boundary

若 hard gating 在 transition DGP 上明显提高方差，再研究 soft attenuation。

## Stage D. Failure-mode simulations

### D1. Spatial autocorrelated errors, constant coefficients

真值没有 coefficient regime，但误差具有 SAR/CAR 型空间相关。

目的：检查 GR-GWR 是否错误发现机制区。

### D2. Local multicollinearity

检查 GWR coefficient fingerprints 是否被局部共线性驱动。

### D3. Uneven sampling density

检查 graph-cut penalty 与 kNN 构图是否受密度影响。

### D4. Conflicting covariate boundaries

`beta_1(s)` 与 `beta_2(s)` 的真实边界不同。

目的：评估 shared-regime 假设。

### D5. Fuzzy transition zone

没有锐利 jump，而是连续过渡带。

目的：检查 hard regime gating 的代价。

## Stage E. Competitor benchmarks

至少覆盖：

- OLS / global regression；
- GWR；
- MGWR（若实验允许）；
- spatial regime / clustered regression 代表方法；
- piecewise-smooth / clustered varying coefficient 代表方法；
- boundary-aware GWR/MGWR 代表方法（若代码/实现可获得）。

具体对手在正式跑实验前根据可复现代码可用性确定。

## Stage F. Real-data study

真实数据必须用于展示：

- regime map 是否具有空间解释；
- boundary 附近普通 GWR 与 GR-GWR 的系数差异；
- within-regime coefficient surface 是否仍保留局部变化；
- 结果稳定性，而不是只展示更低 RMSE。

真实数据集暂不冻结。

## Reproducibility rules

每个论文实验应记录：

- random seed；
- data-generation parameters；
- bandwidth / K / lambda / gamma；
- adjacency definition；
- code commit SHA；
- raw metrics；
- figure-generation script。

禁止只保留最终图片而丢失原始实验结果。
