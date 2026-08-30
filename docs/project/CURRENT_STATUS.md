# CURRENT_STATUS

Last updated: 2026-08-30

## Current phase

**Phase 1 — GR-GWR component research**

标准 GWR 基线已经完成端到端验证并封存。当前正在逐组件研究 GR-GWR，已经完成显式 `W`、Georgia Queen adjacency、pilot GWR local coefficients、local relationship edge diagnostics 和 Queen-constrained Ward 初始分区探索。

## Completed

- 建立独立仓库 `GeoRegime-GWR`。
- 加入最小 `BasicGWR` 与 `GRGWRBaseline`。
- 明确当前 GR-GWR 只是 baseline snapshot，不是最终论文算法。
- 建立 `00_PROJECT_HANDOFF.md`、CURRENT_STATUS、设计文档、ADR、实验计划和 CI。
- 引入 canonical Georgia 1990 county benchmark（159 counties）。
- 当前基线验证只研究标准 GWR，不运行 MGWR 模型；外部 `mgwr` 包只用 `mgwr.gwr.GWR` 和 `mgwr.sel_bw.Sel_BW` 做验证。
- Canonical Georgia external reference 已复现：bandwidth 117，RSS 51.186192，trace(S) 11.804770，AICc 299.050809，R2 0.678074。
- 给定相同 bandwidth=117 时，仓库 `BasicGWR` 与 `mgwr.gwr.GWR` 的 159×4 局部系数、拟合值、残差和 hat matrix 达到机器精度一致。
- 已从 current PyGWRx 抽离 adaptive exhaustive AICc search。
- 已实现 fixed-distance golden-section AICc search。
- 已实现 `mgwr==2.2.1` adaptive discrete golden-section compatibility mode。
- 三模式 Georgia 自动带宽端到端 CI 已通过。
- 普通测试 CI 已在 Python 3.10 与 3.12 通过。
- 已建立 Georgia Queen contiguity `W`：159 counties、431 undirected edges、无孤立县、整体连通。
- 已绘制 ordinary BasicGWR 的 159×4 local coefficients；研究默认 adaptive exhaustive bandwidth = 116。
- 已建立当前 exploratory local relationship fingerprint：三个标准化局部斜率 `PctFB`, `PctBlack`, `PctRural`；当前不含 intercept、不含 coordinates。
- 已计算 431 条 Queen 邻接边上的 coefficient-fingerprint Euclidean distance。
- 已生成 `K=2..15` 的 Queen-constrained Ward 初始分区；所有候选 regime 均空间连通。

## Standard GWR bandwidth policy — FROZEN BY ADR-0003

三条路径必须明确区分。

### Adaptive research default

```python
BasicGWR(bandwidth="auto")
```

使用 **exhaustive integer AICc search**：

- selector: `PyGWRxAdaptiveAICcSelector`；
- Georgia 搜索范围：8–159；
- selected `k = 116`；
- final AICc = `298.9855995813665`；
- 这是后续论文实验中标准 GWR baseline 的默认 adaptive 策略。

### Fixed-distance research default

```python
BasicGWR(bandwidth="auto", adaptive=False)
```

使用 **continuous golden-section AICc search**：

- selector: `FixedGoldenAICcSelector`；
- fixed bandwidth 是连续距离变量，因此使用连续一维优化；
- 功能已实现，但尚未建立独立外部 fixed-GWR benchmark，不得提前声称已经外部数值验证。

### Historical mgwr reproduction mode

```python
BasicGWR(
    bandwidth="auto",
    adaptive=True,
    search_strategy="mgwr_golden",
)
```

使用 `MGWRCompatibleAICcSelector`：

- Georgia compatibility search range：48–159；
- selected `k = 117`；
- final AICc = `299.0508086830287`；
- external `mgwr.gwr.GWR` AICc = `299.05080868302883`；
- compatibility mode 只用于历史/官方示例复现，不是研究默认值。

正式决策：`docs/decisions/ADR-0003-standard-gwr-bandwidth-search-policy.md`。

## Georgia validation — PASSED

固定规格：

- n = 159；
- response: `PctBach`；
- predictors: `PctFB`, `PctBlack`, `PctRural`；
- projected coordinates: `X`, `Y`；
- X/y z-score (`ddof=0`)；
- adaptive bisquare；
- Gaussian GWR AICc。

实际 CI 结果：

1. external `mgwr.sel_bw.Sel_BW` = **117**；
2. BasicGWR adaptive exhaustive = **116**；
3. BasicGWR explicit `mgwr_golden` = **117**；
4. compatibility 117 与 external GWR 的最大局部参数差 = `9.99e-16`；
5. 最大拟合值差 = `2.66e-15`；
6. 最大 hat-matrix 差 = `5.55e-17`；
7. RSS 差 = `-7.11e-15`；
8. AICc 差 = `-1.14e-13`；
9. `passes_validation = true`；
10. exhaustive 116 的 AICc 明确低于 compatibility 117。

关键证据：

- `experiments/validation/basicgwr_vs_mgwr_georgia.py`
- `results/validation/basicgwr_vs_mgwr_georgia/summary.json`
- `results/validation/basicgwr_vs_mgwr_georgia/strict_exhaustive_bandwidth_curve.csv`
- `results/validation/basicgwr_vs_mgwr_georgia/mgwr_compatible_search_trace.csv`
- `.github/workflows/basicgwr-validation.yml`

## Current GR-GWR baseline flow

当前原始 baseline 代码仍是：

1. 从坐标建立 kNN + MST 空间邻接；
2. 全域普通 GWR；
3. 局部斜率标准化；
4. 斜率 + 坐标组成 clustering features；
5. 空间约束 Ward 得到初始 regime；
6. 只在同 regime 内重新做 GWR；
7. ICM 风格逐点 refinement；
8. `RSS + lambda * graph-cut count` objective guard；
9. 最终重新拟合。

**注意：当前 Georgia component experiments 已经与这份旧 baseline snapshot 有意分离。不要把 baseline 代码里的 kNN+MST、coordinates-in-features、ICM 等当成已确认方案。**

## Current exploratory component choices — NOT FROZEN

Georgia 当前探索规格：

1. `D` 与 `W` 分离：`D` 用于 GWR kernel distance，`W` 用于 topology / allowed merges；
2. polygon `W` 当前使用 Queen contiguity；
3. ordinary GWR 只作为 pilot local-relationship estimator；
4. fingerprint 暂用三个标准化 local slopes；
5. intercept 暂不进入 fingerprint；
6. coordinates 暂不进入 fingerprint；
7. similarity 暂用 standardized-slope Euclidean distance；
8. initial segmentation 暂用 Queen-constrained Ward；
9. `K=2..15` 仅为 exploratory scan；
10. 当前只把 `K=6` 作为 working candidate 做进一步检查，**不是最终 K**。

当前 K=6 初始结果：

- WCSS = `41.17628003831048`；
- regime sizes = `53, 19, 17, 22, 27, 21`；
- all regimes connected = true。

## Pending method questions — MUST RESOLVE LATER

完整记录见：

`docs/project/PENDING_GRGWR_METHOD_QUESTIONS.md`

目前最重要的未解决问题包括：

- 普通 GWR 会平滑真实边界，为什么还能作为分区依据？pilot GWR 是否会丢失过多 boundary signal？
- 是否需要 `initial segmentation -> regime-restricted GWR -> label refinement` 的后续迭代修正？
- local relationship fingerprint 是否应该包含 intercept、估计不确定性或其他信息？
- coordinates 是否应该永远排除，还是仅在显式 `W` 存在时排除？
- Ward 是最终 segmentation 方法还是只应作为 initializer？
- regime 数量 K 如何选择？WCSS 不能单独决定 K。
- `K_max` 应如何由最小 regime size / 可识别性约束定义，而不是人为设 15？
- 最终 K 是否应主要依据 spatial/blocked CV，并结合 AICc/BIC、stability 和 1-SE rule？

这些问题目前都**不得写成已解决或已冻结**。

## Not frozen

以下 GR-GWR 组件全部尚未冻结：

- polygon Queen/Rook 的最终默认；
- point-data `W` 构造；
- coordinates 是否进入 clustering features；
- intercept 是否进入 fingerprint；
- fingerprint metric / uncertainty weighting；
- Ward 是否保留；
- regime 数量选择；
- K search upper bound；
- regime minimum size；
- ICM 是否保留；
- boundary penalty；
- shared vs variable-specific regimes；
- soft boundary / partial borrowing；
- GR-GWR 最终 bandwidth policy；
- 参数选择和统计推断。

## Current highest-priority task

**先单独检查 Georgia 的 exploratory K=6 initial partition 本身是否合理。当前不做 boundary-aware refit，不把 K=6 冻结为最终选择。**

## Next tasks

1. 单独绘制并检查 K=6：空间形状、regime size、各 regime fingerprint 中心与离散程度。
2. 检查是否存在狭长连接、局部碎片、内部 coefficient inconsistency。
3. 在看清初始 K=6 后，再决定是否进入 regime-restricted refit。
4. 后续必须用 synthetic boundary DGP 回答 pilot-GWR boundary smoothing 问题。
5. 后续建立正式 K-selection protocol：spatial/blocked CV + complexity/stability/min-size evidence。
6. 做 coordinates/intercept/fingerprint/Ward 等消融。
7. 在证据充分前，不冻结 GR-GWR v1 paper algorithm。
