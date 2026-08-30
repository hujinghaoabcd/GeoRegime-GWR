# CURRENT_STATUS

Last updated: 2026-08-30

## Current phase

**Phase 0 — Research repository extraction and baseline recovery**

目标：把 GR-GWR 从 pyGWRx 的成熟软件包环境中抽离出来，建立一个可快速修改、可做论文实验、可跨对话持续开发的独立研究仓库。

## Completed

- 建立独立仓库 `GeoRegime-GWR`。
- 加入最小 `BasicGWR`。
- 加入 `GRGWRBaseline`，保留当前 pyGWRx 的主要流程。
- 明确当前实现只是 baseline，不是最终论文算法。
- 建立项目交接、架构索引、设计与实验文档体系。
- 加入基础 GWR / GR-GWR invariants tests 与 GitHub Actions。
- Python 3.10 与 Python 3.12 CI 均已通过。
- 引入 canonical PySAL/libpysal Georgia 1990 county benchmark（159 counties），保留 `GData_utm.csv` 与 `G_utm.*` polygon files。
- 使用外部 `mgwr==2.2.1` 成功复现 Oshan et al. (2019) Georgia GWR/MGWR 示例：
  - GWR bandwidth = 117；AICc = 299.050809；R2 = 0.678074；
  - MGWR bandwidths = [92, 101, 136, 158]；AICc = 297.120138；R2 = 0.679878；
  - 与历史官方 notebook 的带宽完全一致，主要诊断仅存在四舍五入级差异。
- 外部复现代码：`experiments/real_data/georgia_mgwr_2019/reproduce.py`。
- 机器可读复现结果：`results/reproduction/georgia_mgwr_2019/summary.json`。
- 当前对模型流程已经恢复到以下认识：
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

Georgia benchmark 固定使用 Oshan et al. / PySAL 示例规格，不与其他 `GeorgiaEduc` 版本混用：

- response: `PctBach`；
- predictors: `PctFB`, `PctBlack`, `PctRural`；
- projected coordinates: `X`, `Y`；
- z-score X 与 y (`ddof=0`)；
- adaptive bisquare kernel；
- AICc bandwidth selection。

后续 BasicGWR、GR-GWR 及其他对照模型若使用该 benchmark，必须显式记录任何偏离上述规格的地方。

## Current highest-priority question

**Should GR-GWR be reformulated around an explicit spatial adjacency matrix W rather than a hard-coded kNN+MST graph construction?**

当前讨论倾向：

- 理论层面使用标准 `W` 更符合 GIS / spatial statistics 表述；
- polygon 可用 Queen/Rook contiguity；
- point data 可用 kNN / Delaunay / distance graph；
- 图只是 `W` 的计算表示，不应成为模型定义本身。

但该设计尚未正式接受为 v1 paper algorithm，需要先做基线与消融实验。

## Next tasks

1. 用仓库内 `BasicGWR` 在同一 Georgia 规格下对照外部 `mgwr.GWR`，确认基础 GWR 引擎数值一致性。
2. 建立三个最基础 DGP：
   - globally smooth；
   - piecewise constant；
   - within-regime smooth + between-regime jump。
3. 建立 no-refinement vs ICM-refinement 消融。
4. 建立 coordinates-in-clustering-features vs no-coordinates 消融。
5. 设计显式 `W` 接口，并比较：
   - Queen/Rook；
   - kNN；
   - kNN+MST。
6. 在上述结果出来前，不宣布任何算法组件为最终方案。
