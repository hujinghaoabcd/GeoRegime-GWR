# CURRENT_STATUS

Last updated: 2026-08-30

## Current phase

**Phase 0 — Research repository extraction and standard-GWR baseline recovery**

目标：把 GR-GWR 从 pyGWRx 的成熟软件包环境中抽离出来，建立一个可快速修改、可做论文实验、可跨对话持续开发的独立研究仓库。**标准 GWR 基线现在已经完成端到端验证，可以进入 GR-GWR 组件研究。**

## Completed

- 建立独立仓库 `GeoRegime-GWR`。
- 加入最小 `BasicGWR`。
- 加入 `GRGWRBaseline`，保留当前 pyGWRx 的主要流程。
- 明确当前实现只是 baseline，不是最终论文算法。
- 建立项目交接、架构索引、设计、ADR 与实验文档体系。
- 加入基础 GWR / GR-GWR invariants tests 与 GitHub Actions。
- Python 3.10 与 Python 3.12 CI 均已通过。
- 引入 canonical PySAL/libpysal Georgia 1990 county benchmark（159 counties），保留 `GData_utm.csv` 与 `G_utm.*` polygon files。
- 当前阶段明确：**不研究 MGWR，只研究标准 GWR**。`mgwr` 只作为外部验证来源，实际只调用 `mgwr.gwr.GWR` 与 `mgwr.sel_bw.Sel_BW`。
- 使用 `mgwr==2.2.1` 成功复现 canonical Georgia standard-GWR 示例：bandwidth 117，RSS 51.186192，ENP 11.804770，AICc 299.050809，R2 0.678074，adjusted R2 0.652080。
- 注意：`mgwr.GWRResults.RSS` 是 location-wise weighted RSS diagnostic；模型级总 RSS 在本项目中固定按 `sum(resid_response**2)` 计算。
- 外部标准-GWR复现代码：`experiments/real_data/georgia_gwr_baseline/reproduce.py`。
- 机器可读结果：`results/reproduction/georgia_gwr_baseline/summary.json`。

## BasicGWR end-to-end validation — PASSED

验证不再使用手工给定的 bandwidth=117。仓库自己的：

`BasicGWR(bandwidth="auto")`

现在会自行完成自动带宽搜索，再完成标准 GWR 拟合。

Canonical Georgia 条件：

- 159 counties；
- `PctBach ~ PctFB + PctBlack + PctRural`；
- X/y z-score (`ddof=0`)；
- adaptive bisquare；
- AICc bandwidth selection。

端到端结果：

- external `mgwr.sel_bw.Sel_BW` selected bandwidth = `117`；
- repository `BasicGWR` selected bandwidth = `117`；
- bandwidth equality = **exact**；
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
- tolerance `1e-12` end-to-end test = **PASS**。

因此标准 GWR 的两层都已验证：

1. 自动带宽搜索；
2. 给定所选带宽后的局部 GWR 数值计算。

关键文件：

- `src/georegime_gwr/gwr.py`
- `src/georegime_gwr/bandwidth.py`
- `experiments/validation/basicgwr_vs_mgwr_georgia.py`
- `results/validation/basicgwr_vs_mgwr_georgia/summary.json`
- `results/validation/basicgwr_vs_mgwr_georgia/pointwise_comparison.csv`
- `results/validation/basicgwr_vs_mgwr_georgia/basicgwr_bandwidth_search_curve.csv`
- `.github/workflows/basicgwr-validation.yml`

## Important bandwidth-search distinction

当前 pyGWRx 与 canonical `mgwr==2.2.1` 的自动搜索策略不是同一个算法，不能混称为相同：

- **current PyGWRx**：对合法整数 adaptive bandwidth 做 exhaustive AICc scan；Georgia 上全曲线最低点为 `k=116`；
- **mgwr 2.2.1 canonical behavior**：离散 golden-section AICc search；Georgia 返回 `k=117`。

因此仓库同时保留：

- `PyGWRxAdaptiveAICcSelector`：忠实保留当前 pyGWRx 穷举策略；
- `MGWRCompatibleAICcSelector`：复现 canonical mgwr 2.2.1 搜索行为。

`BasicGWR(bandwidth="auto")` 当前采用后者，因为当前任务是**完全复现 canonical Georgia standard-GWR 实验，包括自动带宽 117**。

这一决策已记录在：

`docs/decisions/ADR-0003-standard-gwr-bandwidth-search-policy.md`

后续论文正式实验若决定改用 exhaustive global AICc minimum，必须单独记录，不得把 116 与 117 的区别隐藏掉。

## Current GR-GWR baseline flow

当前原始 GR-GWR baseline 流程：

1. 空间邻接结构；
2. 第一轮普通 GWR；
3. 局部斜率标准化；
4. 斜率 + 坐标形成聚类特征；
5. 空间约束 Ward 初始分区；
6. 同 regime 内重新 GWR；
7. ICM 风格逐点 refinement；
8. `RSS + graph-cut penalty` 全局 guard；
9. 最终重新拟合。

## Not frozen

以下全部尚未冻结：

- 是否必须使用 kNN + MST；
- 是否直接以标准空间邻接矩阵 `W` 为模型输入；
- polygon 是否优先 Queen/Rook；
- 是否把坐标继续放入 clustering features；
- 是否继续使用 Ward；
- regime 数量如何选择；
- 是否保留 ICM refinement；
- boundary penalty 的形式；
- shared regime 与 variable-specific regime；
- 是否需要 soft boundary / partial borrowing；
- GR-GWR 最终 bandwidth/search policy；
- 参数选择标准与统计推断。

## Canonical real-data benchmark

Georgia benchmark 固定使用 canonical PySAL/mgwr standard-GWR 规格，不与其他 `GeorgiaEduc` 版本混用：

- response: `PctBach`；
- predictors: `PctFB`, `PctBlack`, `PctRural`；
- projected coordinates: `X`, `Y`；
- z-score X 与 y (`ddof=0`)；
- adaptive bisquare kernel；
- AICc bandwidth selection；
- standard GWR only。

## Current highest-priority question

**With the standard GWR engine now validated end to end, should GR-GWR be reformulated around an explicit spatial adjacency matrix W rather than a hard-coded kNN+MST graph construction?**

当前讨论倾向：理论层面使用标准 `W` 更符合 GIS / spatial statistics 表述；polygon 可用 Queen/Rook；point data 可用 kNN / Delaunay / distance graph；图只是 `W` 的计算表示，不应成为模型定义本身。但尚未冻结，需要实验决定。

## Next tasks

1. 建立三个最基础 DGP：globally smooth、piecewise constant、within-regime smooth + between-regime jump。
2. 设计显式 `W` 接口，并比较 Queen/Rook、kNN、kNN+MST。
3. 建立 coordinates-in-clustering-features vs no-coordinates 消融。
4. 建立 no-refinement vs ICM-refinement 消融。
5. 再评估 Ward、regime number selection、boundary penalty，以及 GR-GWR 最终 bandwidth policy。
6. 在上述结果出来前，不宣布任何 GR-GWR 算法组件为最终方案。
