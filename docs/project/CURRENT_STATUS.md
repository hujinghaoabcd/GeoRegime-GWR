# CURRENT_STATUS

Last updated: 2026-08-30

## Current phase

**Phase 0 — Research repository extraction and standard-GWR baseline recovery**

目标：把 GR-GWR 从 pyGWRx 的成熟软件包环境中抽离出来，建立一个可快速修改、可做论文实验、可跨对话持续开发的独立研究仓库；当前先把标准 GWR 基线彻底固定，再进入 GR-GWR 修改。

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

## Current highest-priority task

**Validate repository `BasicGWR` numerically against the trusted external `mgwr.gwr.GWR` Georgia result before changing GR-GWR.**

通过后再继续研究：

- 是否把 GR-GWR 重构为显式空间邻接矩阵 `W`；
- polygon 的 Queen/Rook；
- coordinates-in-features；
- Ward；
- ICM。

## Next tasks

1. 用仓库内 `BasicGWR` 在同一 Georgia 规格下对照外部 `mgwr.gwr.GWR`：
   - bandwidth；
   - 159 × 4 local parameters；
   - fitted values；
   - residuals；
   - hat-matrix based diagnostics。
2. 修正 `BasicGWR` 与标准 GWR 的任何数值/定义差异，直到基线可信。
3. 建立三个最基础 DGP：globally smooth、piecewise constant、within-regime smooth + between-regime jump。
4. 再进入 `W / coordinates-in-features / Ward / ICM` 的逐组件研究。
5. 在上述结果出来前，不宣布任何 GR-GWR 算法组件为最终方案。
