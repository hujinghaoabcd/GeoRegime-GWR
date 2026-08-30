# PENDING_GRGWR_METHOD_QUESTIONS

Last updated: 2026-08-30

> 这些问题来自当前 Georgia 真实数据探索。它们是**后续必须解决的问题**，不是已经冻结的算法决策。任何后续对 GR-GWR 方法的修改都应回到本文逐项处理，并在有充分理论、模拟和真实数据证据后再形成 ADR。

## 1. 用普通 GWR 的结果做初始分区是否合理？

当前探索流程先用 ordinary GWR 得到局部系数，再用这些系数识别初始 regime。

这存在一个必须正面回答的理论问题：普通 GWR 本身会跨真实边界借样本，并倾向于把突变关系平滑化；如果真实边界被 pilot GWR 严重平滑，后续基于局部系数的分区可能无法恢复真实边界。

当前暂定解释：ordinary GWR 只作为 **pilot / initialization estimator**，提供局部关系的初始连续场，而不是最终 boundary-aware estimator。

后续必须验证：

- 在已知真实边界的 synthetic DGP 中，pilot GWR 还能保留多少边界信号；
- `GWR -> clustering` 是否足以恢复边界；
- 是否必须加入 `initial segmentation -> regime-restricted GWR -> label refinement` 的迭代修正；
- 完整 GR-GWR 是否确实比普通 GWR 减少 boundary smoothing。

在这些实验完成前，不得声称该问题已经解决。

## 2. 什么才是合适的 local relationship fingerprint？

当前 Georgia 初始探索采用三个标准化局部斜率：

`PctFB`, `PctBlack`, `PctRural`

即：

`f_i = (z beta_PctFB, z beta_PctBlack, z beta_PctRural)`

当前探索中：

- 截距不参与 fingerprint；
- coordinates 不参与 fingerprint；
- 三个斜率按列 z-score (`ddof=0`)；
- 差异使用 Euclidean distance。

这些都只是当前实验规格，不是最终方法定义。

后续必须验证：

- 是否应该排除截距；
- coefficient standardization 是否合理；
- Euclidean distance 是否足够；
- 是否需要考虑局部系数估计不确定性，例如按标准误或协方差加权；
- 是否需要加入其他 local relationship diagnostics。

## 3. coordinates 是否应该进入 clustering features？

当前 Georgia 探索暂时**不把 coordinates 放进 clustering features**。

理由是空间连续性已经由显式邻接矩阵 `W` 负责，如果再把坐标混进系数特征距离，可能重复计算 geography，并混淆：

- `D`：用于 GWR kernel，回答“离多远”；
- `W`：用于 topology / allowed merges，回答“是不是邻居”；
- local coefficient features：用于 relationship similarity，回答“关系是否相似”。

后续必须做 coordinates-in-features vs no-coordinates 的消融后再决定。

## 4. Ward 是否应该保留？

当前初始分区使用 **Queen-constrained Ward agglomerative clustering**。

Ward 的作用是：在只允许空间相邻区域合并的前提下，每一步选择使组内平方和增加最小的合并。

当前 Georgia 实验表明 Queen constraint 可以保证候选 regime 连通，但这不意味着 Ward 是最终最佳分区算法。

后续必须比较或讨论：

- Ward 的统计目标是否与 GR-GWR 的最终目标一致；
- 是否需要直接优化 regression loss + boundary penalty；
- 是否需要替代聚类/区域化算法；
- Ward 是否只适合作为 initializer。

## 5. regime 数量 K 怎么选？

这是当前最容易被审稿人质疑的问题之一。

Georgia 初始探索只是扫描 `K=2..15`。`15` 是探索上限，不是理论上限，也不是最终方法规定。

当前结果中 `K=6` 是一个值得进一步检查的候选，但**尚未冻结**。

重要原则：

- WCSS 随 K 增大必然下降，因此不能单独用 WCSS 选择最终 K；
- elbow 图只能用于探索和缩小候选范围；
- 最终 K 应在完整 GR-GWR 模型层面选择，而不是只在 clustering 层选择。

后续优先考虑的选择证据：

1. spatial / blocked cross-validation prediction error（RMSE、MAE）作为主要泛化证据；
2. AICc / BIC 等模型复杂度指标；
3. regime 连通性与最小样本量约束；
4. partition stability / sensitivity；
5. 若多个 K 的预测性能近似，可用 1-SE rule 选择更小、更简单的 K。

还必须解决 `K_max` 的定义。后续更合理的方向是由最小可接受 regime size 和局部回归可识别性约束决定，而不是人为固定为 15。

## 6. WCSS 的定位

WCSS = within-cluster sum of squares，只衡量当前 fingerprint space 中各 regime 内部的紧凑程度。

它可以：

- 描述不同 K 下的分区紧凑性；
- 观察 elbow；
- 辅助产生候选 K。

它不能：

- 证明某个 K 是最终正确的 regime 数量；
- 证明分区具有真实机制意义；
- 替代完整 GR-GWR 的 out-of-sample validation。

## 7. 当前 Georgia 初始分区证据

当前 experimental specification：

- n = 159 counties；
- ordinary BasicGWR adaptive exhaustive bandwidth = 116；
- fingerprint = 3 standardized local slopes；
- spatial topology = Queen contiguity；
- Queen undirected edges = 431；
- Ward candidate K = 2..15；
- all candidate regimes connected；
- intercept excluded；
- coordinates excluded。

当前 WCSS / size 结果显示：

- K=5: WCSS 54.1594；
- K=6: WCSS 41.1763；regime sizes = 53, 19, 17, 22, 27, 21；
- K=7: WCSS 36.4078，并开始出现 size=7 的较小 regime。

因此 **K=6 只作为 exploratory working candidate**，不将其写成最终选择。

## 8. K=6 分区内重新 GWR：第一轮探索结果

已经完成最简单的 fixed-label experiment：

`ordinary GWR -> fixed K=6 Ward regimes -> each regime fits its own BasicGWR`

本轮严格不做：

- label refinement；
- ICM；
- boundary penalty；
- regime reassignment；
- iterative repartition。

跨 regime 的观测在局部拟合中完全不可借用，即 cross-regime weights = 0。

为了先看结果，本轮每个 regime **独立使用 adaptive exhaustive AICc 自动选择 bandwidth**。这只是实验规格，不是最终 GR-GWR bandwidth policy。

关键发现：6 个 regime 的自动带宽全部等于各自 regime 的样本量：

- R1: n=53, bw=53；
- R2: n=19, bw=19；
- R3: n=17, bw=17；
- R4: n=22, bw=22；
- R5: n=27, bw=27；
- R6: n=21, bw=21。

这意味着在当前真实数据和当前 K=6 划分下，AICc 在每个区内部都倾向于使用尽可能大的局部邻域，而不是更局部的小 bandwidth。这个现象必须后续专门解释。

总体 in-sample 对比：

- ordinary GWR: RSS=51.0829, RMSE=0.5668, MAE=0.3994, R2=0.6787, trace(S)=11.9121, AICc=298.9856；
- K=6 restricted GWR: RSS=43.1001, RMSE=0.5206, MAE=0.3450, R2=0.7289, trace(S)=39.7019, AICc=354.0115。

因此当前结果是：

- 分区限制后 in-sample 拟合误差明显下降；
- 但有效复杂度 trace(S) 大幅上升；
- 按当前整体 Gaussian GWR AICc 计算，AICc 反而显著变差（+55.03）。

这说明“分区后误差下降”本身不能证明模型更好，也进一步说明最终必须使用更严格的 complexity / spatial CV 证据。

### 新增待解决问题：GR-GWR bandwidth 到底怎么定义？

后续必须比较至少：

1. 每个 regime 独立选择 adaptive bandwidth；
2. 所有 regime 共用一个 bandwidth，但禁止跨边界借样本；
3. 使用全局 distance bandwidth 再做 regime masking；
4. 对 regime size / effective degrees of freedom 加更明确的限制；
5. 用 spatial CV 而不是纯 in-sample AICc 选择 bandwidth。

尤其要回答：当自动带宽总是顶到 regime size 上限时，GR-GWR 是否实际上更接近“分区内空间加权回归 / 接近区域回归”，以及这是否符合模型设计目标。

## 9. 接下来的工作顺序

当前已有两层结果：初始 K=6 分区，以及固定 K=6 的最简单 regime-restricted GWR。

下一阶段再逐项研究：

- pilot GWR boundary smoothing；
- K selection；
- per-regime vs shared bandwidth；
- Ward 是否只是 initializer；
- 是否需要 regime-restricted refit + iterative label refinement；
- 如何控制因分区导致的有效复杂度增长；
- spatial / blocked CV 是否支持真实泛化改善。
