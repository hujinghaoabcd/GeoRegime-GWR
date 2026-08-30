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

## Current highest-priority question

**Should GR-GWR be reformulated around an explicit spatial adjacency matrix W rather than a hard-coded kNN+MST graph construction?**

当前讨论倾向：

- 理论层面使用标准 `W` 更符合 GIS / spatial statistics 表述；
- polygon 可用 Queen/Rook contiguity；
- point data 可用 kNN / Delaunay / distance graph；
- 图只是 `W` 的计算表示，不应成为模型定义本身。

但该设计尚未正式接受为 v1 paper algorithm，需要先做基线与消融实验。

## Next tasks

1. 跑通当前 baseline smoke test。
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
