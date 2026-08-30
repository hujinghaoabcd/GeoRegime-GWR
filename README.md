# GeoRegime-GWR

Research repository for **Geo-Regime Geographically Weighted Regression (GR-GWR)**.

本仓库用于 GR-GWR 方法论文的算法研究、仿真、消融、对比与实证实验。它不是 pyGWRx 的替代品，也不是成熟软件包；当前代码只是一个可读、可修改的研究基线。

## Start here

新的对话或新的研究阶段请按以下顺序恢复项目：

1. `00_PROJECT_HANDOFF.md`
2. `ARCHITECTURE_INDEX.md`
3. `docs/project/CURRENT_STATUS.md`
4. `docs/design/GRGWR_BASELINE_SPEC.md`
5. `docs/design/RESEARCH_QUESTIONS.md`
6. `docs/experiments/EXPERIMENT_PLAN.md`
7. `docs/decisions/` 最新 ADR

## Current baseline

当前基线流程：

1. kNN + MST 空间邻接结构；
2. 全域普通 GWR 得到每个位置的局部参数；
3. 标准化局部斜率，并与标准化坐标组合；
4. 空间连通约束 Ward 得到初始 regime；
5. 仅在同 regime 内进行局部 GWR；
6. 可选 ICM 风格逐点 refinement；
7. `RSS + graph-cut penalty` 做全局 objective guard；
8. 最终重新拟合局部系数。

**以上不是最终论文算法。任何组件都允许通过理论、文献和实验被修改或删除。**

## Minimal code layout

```text
GeoRegime-GWR/
├─ 00_PROJECT_HANDOFF.md
├─ ARCHITECTURE_INDEX.md
├─ pyproject.toml
├─ src/georegime_gwr/
│  ├─ gwr.py
│  └─ grgwr.py
├─ tests/
│  ├─ test_gwr.py
│  └─ test_grgwr.py
├─ experiments/
├─ results/
└─ docs/
   ├─ project/
   ├─ design/
   ├─ experiments/
   ├─ literature/
   └─ decisions/
```

## Install and test

```bash
python -m pip install -e .[dev]
pytest
```

手工 smoke run：

```bash
python scripts/run_smoke.py
```

## Research rule

代码存在不等于理论已经确认。每个关键算法组件都应经过：

**source behavior -> mathematical statement -> literature check -> targeted simulation -> keep/modify/remove decision**。

当模型定义发生变化时，必须同步更新 `CURRENT_STATUS.md`、ADR 和 `SESSION_LOG.md`，确保以后任何对话都能准确接手。
