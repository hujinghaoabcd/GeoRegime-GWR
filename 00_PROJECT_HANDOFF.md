# 00_PROJECT_HANDOFF

> **跨对话最高优先级入口。任何新的 ChatGPT / Codex 对话继续本仓库前，必须先读本文。**

## 1. 项目定位

`GeoRegime-GWR` 是 **Geo-Regime Geographically Weighted Regression (GR-GWR)** 方法论文的轻量研究仓库，不是成熟软件包。

用途：

- 独立理解、重构和验证 GR-GWR；
- 保留可信的最小标准 GWR baseline；
- 快速修改论文算法；
- 做 synthetic / real-data / ablation / benchmark 实验；
- 生成论文结果；
- 通过设计文档、ADR 和交接文档保证跨对话连续开发。

## 2. 最高原则

### 2.1 当前 GR-GWR 不是最终算法

`src/georegime_gwr/grgwr.py` 是 baseline snapshot。后续完全允许删除或替换：

- ICM；
- Ward；
- kNN + MST；
- clustering features；
- boundary penalty；
- regime number selection；
- bandwidth policy；
- 整个优化框架。

不要因为现有代码已经实现就把它当作论文最终定义。

### 2.2 研究顺序

每个组件遵循：

`source behavior -> mathematical statement -> literature check -> targeted simulation -> keep / modify / remove`

### 2.3 方法定义变化必须留痕

至少同步：

- `docs/project/CURRENT_STATUS.md`
- 对应设计规范
- `docs/decisions/` 新增或更新 ADR
- 若影响整体流程，更新本文

## 3. 标准 GWR 基线 — 已验证完成

当前基础验证只研究 **standard GWR**，不运行 MGWR 模型。

外部 `mgwr==2.2.1` 只作为验证 oracle，使用：

- `mgwr.gwr.GWR`
- `mgwr.sel_bw.Sel_BW`

Canonical Georgia benchmark：

- 159 Georgia counties，1990 Census；
- response: `PctBach`；
- predictors: `PctFB`, `PctBlack`, `PctRural`；
- coords: projected `X`, `Y`；
- X 与 y 均 z-score (`ddof=0`)；
- adaptive bisquare；
- Gaussian GWR AICc。

External canonical result：

- bandwidth = 117；
- RSS ≈ 51.186192；
- trace(S) ≈ 11.804770；
- AICc ≈ 299.050809；
- R2 ≈ 0.678074。

给定相同 `k=117` 时，仓库 `BasicGWR` 与 external `mgwr.gwr.GWR` 的局部系数、拟合值、残差和 hat matrix 已验证到机器精度一致。

普通测试 CI 已在 Python 3.10 / 3.12 通过；Georgia 自动带宽端到端 CI 也已通过。

## 4. 极重要：BasicGWR 有三种明确 bandwidth search policy

**不要把 116 与 117 当成实现错误，也不要把 compatibility mode 当成默认研究策略。**

正式依据：`docs/decisions/ADR-0003-standard-gwr-bandwidth-search-policy.md`。

### 4.1 Adaptive research default — exhaustive integer

```python
BasicGWR(bandwidth="auto")
```

默认 adaptive，使用：

`PyGWRxAdaptiveAICcSelector`

Georgia 实际验证：

- search range = 8–159；
- selected `k = 116`；
- AICc = `298.9855995813665`；
- 这是论文研究中标准 GWR baseline 的默认 adaptive 策略。

### 4.2 Fixed-distance research default — golden section

```python
BasicGWR(bandwidth="auto", adaptive=False)
```

使用：

`FixedGoldenAICcSelector`

对连续 fixed-distance bandwidth 做 golden-section AICc optimization。

**该路径已经实现，但还没有独立外部 fixed-GWR benchmark。不要提前声称它已经完成外部数值验证。**

### 4.3 Historical mgwr reproduction — explicit compatibility mode

```python
BasicGWR(
    bandwidth="auto",
    adaptive=True,
    search_strategy="mgwr_golden",
)
```

使用：

`MGWRCompatibleAICcSelector`

复现 `mgwr==2.2.1` 的离散 golden-section search 和 compact-kernel boundary semantics。

Georgia 实际验证：

- search range = 48–159；
- selected `k = 117`；
- AICc = `299.0508086830287`；
- external mgwr AICc = `299.05080868302883`；
- 用于复现论文/官方示例；
- **不是默认研究模式。**

### 4.4 Why 116 vs 117?

Current PyGWRx adaptive exhaustive search 会访问所有整数候选，因此 Georgia 上找到完整离散 AICc 曲线最低点 `116`。

Historical mgwr 2.2.1 默认采用 rounded discrete golden-section，只访问部分候选并返回 `117`。

因此：

- 116 vs 117 = search algorithm difference；
- 不是 standard GWR local WLS implementation difference。

## 5. Georgia automatic-search validation — PASSED

`experiments/validation/basicgwr_vs_mgwr_georgia.py` 已实际验证：

1. external mgwr search = **117**；
2. repository research default exhaustive = **116**；
3. repository explicit `mgwr_golden` = **117**；
4. compatibility-mode 117 与 external `mgwr.gwr.GWR`：
   - max local-parameter difference = `9.99e-16`；
   - max fitted-value difference = `2.66e-15`；
   - max hat-matrix difference = `5.55e-17`；
   - RSS difference = `-7.11e-15`；
   - AICc difference = `-1.14e-13`；
5. exhaustive 116 AICc < compatibility 117 AICc；
6. `passes_validation = true`。

证据：

- `results/validation/basicgwr_vs_mgwr_georgia/summary.json`
- `results/validation/basicgwr_vs_mgwr_georgia/strict_exhaustive_bandwidth_curve.csv`
- `results/validation/basicgwr_vs_mgwr_georgia/mgwr_compatible_search_trace.csv`
- `.github/workflows/basicgwr-validation.yml`

## 6. 当前 GR-GWR baseline 流程

当前 baseline 仍是：

1. kNN + MST adjacency；
2. full-domain ordinary GWR；
3. standardized local slopes；
4. slopes + coordinates clustering features；
5. spatially constrained Ward initial regimes；
6. regime-restricted local GWR；
7. ICM-style pointwise refinement；
8. `RSS + lambda * graph-cut boundary count` objective guard；
9. final refit。

该流程只是研究起点，不是冻结算法。

## 7. 当前未冻结问题

最高优先级包括：

- 是否直接定义显式空间邻接矩阵 `W`；
- polygon 是否优先 Queen/Rook；
- point data 是否用 kNN/Delaunay/distance graph；
- 是否取消强制 MST；
- coordinates 是否进入 clustering features；
- Ward 是否保留；
- regime number；
- ICM 是否保留；
- boundary penalty；
- shared vs variable-specific regimes；
- hard boundary vs partial borrowing；
- GR-GWR 正式 bandwidth policy；
- inference / diagnostics。

## 8. 新对话恢复顺序

必须按顺序读：

1. `00_PROJECT_HANDOFF.md`
2. `ARCHITECTURE_INDEX.md`
3. `docs/project/CURRENT_STATUS.md`
4. `docs/decisions/ADR-0003-standard-gwr-bandwidth-search-policy.md`
5. `docs/design/GRGWR_BASELINE_SPEC.md`
6. `docs/design/RESEARCH_QUESTIONS.md`
7. `docs/experiments/EXPERIMENT_PLAN.md`
8. 其他最新 ADR
9. 当前源码与 tests / CI

冲突优先级：

`最新 ADR > CURRENT_STATUS > 正式设计规范 > handoff 中的旧描述 > README`

## 9. 当前下一步

**基础 GWR 阶段已经结束，不要再从头复核。**

现在进入 GR-GWR 组件研究：

1. 第一优先级：把 `W` 作为显式空间结构重新设计，明确它与距离矩阵 `D` 的不同职责；
2. 建立 globally smooth / piecewise constant / within-regime smooth + between-regime jump DGP；
3. 比较 Queen/Rook、kNN、kNN+MST；
4. 再做 coordinates、Ward、ICM 等消融；
5. 最终依据理论和实验冻结 GR-GWR v1 paper algorithm。
