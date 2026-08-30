# GR-GWR Baseline Specification

> 本文描述 **当前研究基线实际做什么**，不是最终论文模型定义。

## 1. 输入

给定 `n` 个空间样本：

- 响应变量：`y_i`
- 解释变量：`x_i = (x_{i1}, ..., x_{ip})`
- 空间坐标：`s_i = (u_i, v_i)`

若 `fit_intercept=True`，设计向量为：

`x_i^* = (1, x_{i1}, ..., x_{ip})`

## 2. 空间距离

当前使用欧氏距离：

`d_ij = ||s_i - s_j||_2`

形成距离矩阵 `D = [d_ij]`。

## 3. 当前邻接图

当前 baseline 从坐标自动构造：

`G = G_kNN union G_MST`

- kNN 给每个点局部邻居；
- 对称化后形成无向图；
- MST 用于保证整个图连通。

这是当前实现，不是冻结理论。后续优先研究是否改为显式空间邻接矩阵 `W`。

## 4. 第一轮普通 GWR

对每个位置 `i`：

`w_ij = K(d_ij / b_i)`

然后估计：

`beta_hat_i = (X' W_i X)^(-1) X' W_i y`

因此每个点最终得到一套局部参数：

`beta_hat_i = (beta_0i, beta_1i, ..., beta_pi)`

第一轮 GWR 不知道 regime，因此仍可能跨真实机制边界借样本。

## 5. 分区特征

若有截距，当前代码只取斜率：

`b_i = (beta_1i, ..., beta_pi)`

每一列斜率做 z-score 标准化，得到 `b_tilde_i`。

坐标按列缩放到 `[0, 1]`，得到 `s_tilde_i`。

最终聚类特征：

`f_i = [sqrt(1-gamma) * b_tilde_i, sqrt(gamma) * s_tilde_i]`

其中：

`gamma = spatial_constraint_weight`。

端点：

- `gamma = 0`：只看局部斜率；
- `gamma = 1`：只看坐标；
- 中间值：同时考虑局部关系与位置。

## 6. 初始 regime

当前使用带空间 connectivity 约束的 Ward 层次聚类：

- 特征距离来自 `f_i`；
- 只有空间邻接图允许的区域能够进行层次合并；
- 目标 regime 数由 `n_regimes` 指定；
- 太小的 regime 会被合并。

输出初始标签：

`z_i^(0) in {0, ..., K-1}`。

## 7. Regime-restricted GWR

给定标签 `z` 后，对位置 `i` 仅允许同 regime 样本参与：

`w_ij^GR = K(d_ij / b_i) * I(z_i = z_j)`

于是仍然是“每个点一套局部系数”，而不是“每个 regime 一套固定系数”。

这对应当前 GR-GWR 的核心结构：

- within-regime: local smooth variation
- between-regime: no borrowing across the boundary

## 8. 当前 ICM 风格 refinement

对每个点 `i`，候选 regime 只包括：

- 当前 regime；
- 空间邻接点已经出现的 regime。

对候选 `r`，先把点 `i` 自己排除，用候选 regime 其余样本重新做局部 WLS：

`E_i(r) = (y_i - yhat_i^(-i,r))^2`

再加局部边界不一致惩罚：

`C_i(r) = E_i(r) + lambda * sum_{j in N(i)} I(z_j != r)`

若候选代价严格更小，且源 regime 不会过小/断裂，则允许移动。

## 9. 全局 objective guard

一轮逐点更新后，用新标签重新完整拟合，并计算：

`L(z) = RSS(z) + lambda * B(z)`

其中：

`RSS(z) = sum_i (y_i - yhat_i)^2`

`B(z) = sum_{(i,j) in E} I(z_i != z_j)`

`B(z)` 是 **graph-cut edge count**，不是实际几何边界长度。

若新目标函数变差超过 `tol`，拒绝该轮并停止。

## 10. 最终输出

最终标签重新编号为 `0..K-1`，并重新按最终 regime 做一次完整局部拟合。

主要结果：

- 每个点的 `regime`；
- 每个点的局部回归参数；
- fitted values / residuals；
- objective history；
- regime boundary graph edges。

## 11. 当前基线最需要验证的假设

1. 第一轮 GWR 的局部系数是否足够稳定，可以用作 regime discovery signal？
2. 坐标进入 `f_i` 是否必要，还是与 adjacency constraint 重复？
3. Ward 是否只是方便初始化，还是影响最终结果过大？
4. ICM refinement 是否带来实质增益？
5. hard gating `I(z_i=z_j)` 是否在模糊边界上方差过大？
6. shared regime map 是否足以表达多个 covariates 的空间异质性？
7. `W` 应如何进入正式模型定义？
