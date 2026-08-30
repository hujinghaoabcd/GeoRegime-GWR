# CURRENT_STATUS

Last updated: 2026-08-30

## Current phase

**Phase 1 — GR-GWR component research**

标准 GWR 基线已经完成端到端验证并封存。下一阶段开始逐组件研究 GR-GWR，本阶段最高优先级是空间结构 `W`，然后再研究 clustering features、Ward、ICM 等组件。

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

当前原始 baseline：

1. 从坐标建立 kNN + MST 空间邻接；
2. 全域普通 GWR；
3. 局部斜率标准化；
4. 斜率 + 坐标组成 clustering features；
5. 空间约束 Ward 得到初始 regime；
6. 只在同 regime 内重新做 GWR；
7. ICM 风格逐点 refinement；
8. `RSS + lambda * graph-cut count` objective guard；
9. 最终重新拟合。

## Not frozen

以下 GR-GWR 组件全部尚未冻结：

- kNN + MST 是否保留；
- 是否把标准空间邻接矩阵 `W` 作为显式输入；
- polygon 的 Queen/Rook；
- coordinates 是否继续进入 clustering features；
- Ward 是否保留；
- regime 数量选择；
- ICM 是否保留；
- boundary penalty；
- shared vs variable-specific regimes；
- soft boundary / partial borrowing；
- GR-GWR 最终 bandwidth policy；
- 参数选择和统计推断。

## Current highest-priority task

**基础 GWR 阶段结束。现在进入 GR-GWR 的空间结构设计：是否围绕显式空间邻接矩阵 `W` 重构，而不是把 kNN+MST 写死进模型定义。**

## Next tasks

1. 设计显式 `W` 接口，并明确距离矩阵 `D` 与邻接矩阵 `W` 的不同职责。
2. 建立 globally smooth、piecewise constant、within-regime smooth + between-regime jump 三类基础 DGP。
3. 比较 Queen/Rook、kNN、kNN+MST。
4. 做 coordinates-in-features vs no-coordinates 消融。
5. 做 no-refinement vs ICM-refinement 消融。
6. 再评估 Ward、regime number、boundary penalty 与 GR-GWR 最终 bandwidth policy。
7. 在证据充分前，不冻结 GR-GWR v1 paper algorithm。
