# CURRENT_STATUS

Last updated: 2026-08-30

## Current phase

**Phase 0 — Standard-GWR baseline finalization, then GR-GWR component research**

本仓库是 GR-GWR 方法论文的轻量研究仓库。标准 GWR 的局部回归核心已经与外部 `mgwr.gwr.GWR` 在 Georgia 数据上验证到机器精度；当前正在完成自动带宽搜索策略的显式拆分和最终 CI 验证。完成后进入 GR-GWR 的 `W / clustering / Ward / ICM` 逐组件研究。

## Completed

- 建立独立仓库 `GeoRegime-GWR`。
- 加入最小 `BasicGWR` 与 `GRGWRBaseline`。
- 明确当前 GR-GWR 只是 baseline snapshot，不是最终论文算法。
- 建立 `00_PROJECT_HANDOFF.md`、CURRENT_STATUS、设计文档、ADR、实验计划和 CI。
- 引入 canonical Georgia 1990 county benchmark（159 counties）。
- 当前阶段只研究标准 GWR，不运行 MGWR 模型；外部 `mgwr` 包只用 `mgwr.gwr.GWR` 和 `mgwr.sel_bw.Sel_BW` 做验证。
- Canonical Georgia external reference 已复现：bandwidth 117，RSS 51.186192，trace(S) 11.804770，AICc 299.050809，R2 0.678074。
- 给定相同 bandwidth=117 时，仓库 `BasicGWR` 与 `mgwr.gwr.GWR` 的 159×4 局部系数、拟合值、残差和 hat matrix 达到机器精度一致。
- 已从 current PyGWRx 抽离 adaptive exhaustive AICc search。
- 已实现 fixed-distance golden-section AICc search。
- 已实现 `mgwr==2.2.1` adaptive discrete golden-section compatibility mode。

## Standard GWR bandwidth policy — FROZEN BY ADR-0003

三条路径必须明确区分：

### Adaptive research default

```python
BasicGWR(bandwidth="auto")
```

默认解释为 adaptive bandwidth，并使用 **exhaustive integer AICc search**：

- selector: `PyGWRxAdaptiveAICcSelector`；
- Georgia 全整数候选 AICc 最低点：`k=116`；
- 这是后续论文实验中标准 GWR baseline 的默认策略，因为它消除了离散搜索优化器未访问全局最低候选的影响。

### Fixed-distance research default

```python
BasicGWR(bandwidth="auto", adaptive=False)
```

使用 **continuous golden-section AICc search**：

- selector: `FixedGoldenAICcSelector`；
- fixed bandwidth 是连续距离变量，因此使用连续一维优化；
- 该功能已经实现，但尚未增加独立的外部 fixed-GWR benchmark，不得提前声称其已经外部数值验证。

### Historical mgwr reproduction mode

```python
BasicGWR(
    bandwidth="auto",
    adaptive=True,
    search_strategy="mgwr_golden",
)
```

使用 `MGWRCompatibleAICcSelector`，复现 `mgwr==2.2.1` 的离散 golden-section search。Canonical Georgia 应返回 `k=117`。

**Compatibility mode 不是研究默认值。**

正式决策：`docs/decisions/ADR-0003-standard-gwr-bandwidth-search-policy.md`。

## Georgia validation contract

固定数据规格：

- n = 159；
- response: `PctBach`；
- predictors: `PctFB`, `PctBlack`, `PctRural`；
- projected coordinates: `X`, `Y`；
- X/y z-score (`ddof=0`)；
- adaptive bisquare；
- Gaussian GWR AICc。

最终自动带宽 CI 必须同时满足：

1. external `mgwr.sel_bw.Sel_BW` = 117；
2. BasicGWR adaptive research default = 116；
3. BasicGWR explicit `mgwr_golden` = 117；
4. compatibility-mode 117 的局部参数、拟合、残差、hat matrix、RSS、AICc 与 external `mgwr.gwr.GWR` 在 `1e-12` tolerance 内一致；
5. exhaustive 116 的 AICc 低于 compatibility 117。

关键文件：

- `src/georegime_gwr/gwr.py`
- `src/georegime_gwr/bandwidth.py`
- `experiments/validation/basicgwr_vs_mgwr_georgia.py`
- `results/validation/basicgwr_vs_mgwr_georgia/`
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

**先确认新的三模式 Georgia GWR validation CI 全部通过。**

CI 通过后，基础 GWR 阶段结束，下一问题是：

> GR-GWR 是否应围绕显式空间邻接矩阵 `W` 重构，而不是把 kNN+MST 写死进模型定义？

## Next tasks

1. 确认三模式 Georgia GWR validation CI。
2. 建立 globally smooth、piecewise constant、within-regime smooth + between-regime jump 三类基础 DGP。
3. 设计显式 `W` 接口，比较 Queen/Rook、kNN、kNN+MST。
4. 做 coordinates-in-features vs no-coordinates 消融。
5. 做 no-refinement vs ICM-refinement 消融。
6. 再评估 Ward、regime number、boundary penalty 与 GR-GWR 最终 bandwidth policy。
7. 在证据充分前，不冻结 GR-GWR v1 paper algorithm。
