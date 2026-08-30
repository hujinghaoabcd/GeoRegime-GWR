# CURRENT_STATUS

Last updated: 2026-08-30

## Current phase

**Phase 0 — Research repository extraction and standard-GWR baseline recovery**

目标：把 GR-GWR 从 pyGWRx 的成熟软件包环境中抽离出来，建立一个可快速修改、可做论文实验、可跨对话持续开发的独立研究仓库；当前标准 GWR 基线已经固定并完成数值验证，可以进入 GR-GWR 组件研究。

## Completed

- 建立独立仓库 `GeoRegime-GWR`。
- 加入最小 `BasicGWR`。
- 加入 `GRGWRBaseline`，保留当前 pyGWRx 的主要流程。
- 明确当前实现只是 baseline，不是最终论文算法。
- 建立项目交接、架构索引、设计与实验文档体系。
- 加入基础 GWR / GR-GWR invariants tests 与 GitHub Actions。
- Python 3.10 与 Python 3.12 CI 均已通过。
- 引入 canonical PySAL/libpysal Georgia 1990 county benchmark（159 counties），保留 `GData_utm.csv` 与 `G_utm.*` polygon files。
- 当前阶段明确：**不研究 MGWR，只研究标准 GWR**。`mgwr` 仅作为 Python 包来源，实际只调用 `mgwr.gwr.GWR` 与带宽选择器。
- 使用 `mgwr==2.2.1` 成功复现 canonical Georgia standard-GWR 示例：
  - bandwidth = 117；
  - RSS = 51.186192；
  - ENP = 11.804770；
  - AICc = 299.050809；
  - R2 = 0.678074；
  - adjusted R2 = 0.652080。
- 与 canonical reference 的 bandwidth 完全一致，主要诊断仅存在四舍五入级差异。
- 注意：`mgwr.GWRResults.RSS` 是 location-wise weighted RSS diagnostic；模型级总 RSS 在本项目中固定按 `sum(resid_response**2)` 计算，避免误读包内部属性。
- 外部标准-GWR复现代码：`experiments/real_data/georgia_gwr_baseline/reproduce.py`。
- 机器可读结果：`results/reproduction/georgia_gwr_baseline/summary.json`。
- 159 个县的外部 GWR 局部参数：`results/reproduction/georgia_gwr_baseline/gwr_parameters.csv`。

## BasicGWR numerical validation — PASSED

仓库内 `BasicGWR` 已在 canonical Georgia 数据上，与外部 `mgwr.gwr.GWR` 做逐点数值验证。两者使用完全相同：

- 159 counties；
- `PctBach ~ PctFB + PctBlack + PctRural`；
- X/y z-score (`ddof=0`)；
- adaptive bisquare；
- fixed verified bandwidth = 117。

验证对象包括 159 × 4 local parameters、159 fitted values、159 residuals、完整 159 × 159 hat matrix、RSS 与 trace(S)。

结果：

- max absolute parameter difference = `5.551115123125783e-16`；
- RMSE parameter difference = `9.131398707249937e-17`；
- max absolute fitted difference = `5.551115123125783e-16`；
- max absolute residual difference = `5.551115123125783e-16`；
- max absolute hat-matrix difference = `5.551115123125783e-17`；
- BasicGWR RSS = `51.186191553463985`；
- mgwr.GWR RSS = `51.18619155346399`；
- trace(S) 两者均为 `11.804769716730094`；
- `1e-8` tolerance test = **PASS**。

这些差异处于浮点机器精度量级，因此当前 `BasicGWR` 可视为已数值复现 `mgwr.gwr.GWR` 的标准 Gaussian adaptive-bisquare GWR 核心计算。

关键文件：

- `experiments/validation/basicgwr_vs_mgwr_georgia.py`
- `results/validation/basicgwr_vs_mgwr_georgia/summary.json`
- `results/validation/basicgwr_vs_mgwr_georgia/pointwise_comparison.csv`
- `.github/workflows/basicgwr-validation.yml`

## Current GR-GWR baseline flow

当前对原始 GR-GWR baseline 流程的认识：

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

后续 `BasicGWR`、GR-GWR 及其他对照模型若使用该 benchmark，必须显式记录任何偏离上述规格的地方。

## Current highest-priority question

**With the standard GWR engine now validated, should GR-GWR be reformulated around an explicit spatial adjacency matrix W rather than a hard-coded kNN+MST graph construction?**

当前讨论倾向：

- 理论层面使用标准 `W` 更符合 GIS / spatial statistics 表述；
- polygon 可用 Queen/Rook contiguity；
- point data 可用 kNN / Delaunay / distance graph；
- 图只是 `W` 的计算表示，不应成为模型定义本身。

但该设计尚未正式接受为 v1 paper algorithm，需要通过后续基线与消融实验决定。

## Next tasks

1. 建立三个最基础 DGP：globally smooth、piecewise constant、within-regime smooth + between-regime jump。
2. 设计显式 `W` 接口，并比较 Queen/Rook、kNN、kNN+MST。
3. 建立 coordinates-in-clustering-features vs no-coordinates 消融。
4. 建立 no-refinement vs ICM-refinement 消融。
5. 再评估 Ward、regime number selection 与 boundary penalty。
6. 在上述结果出来前，不宣布任何 GR-GWR 算法组件为最终方案。
