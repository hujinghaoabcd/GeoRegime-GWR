# ARCHITECTURE_INDEX

本文件是仓库结构索引，不定义最终算法。跨对话恢复时先读 `00_PROJECT_HANDOFF.md`，再读本文件。

## 1. 代码

- `src/georegime_gwr/gwr.py`
  - 最小基础 GWR；
  - 负责距离核、固定/自适应带宽、局部加权最小二乘；
  - 不承担 GR-GWR 的 regime 逻辑。

- `src/georegime_gwr/grgwr.py`
  - 当前 GR-GWR research baseline；
  - 忠于 pyGWRx 现有主要流程，但代码更小、更容易修改；
  - **不是最终论文算法。**

## 2. 项目状态

- `docs/project/CURRENT_STATUS.md`
  - 当前阶段、已完成事项、正在讨论的问题、下一步任务。

## 3. 模型设计

- `docs/design/GRGWR_BASELINE_SPEC.md`
  - 当前 baseline 的逐步数学与代码对应关系。

- `docs/design/RESEARCH_QUESTIONS.md`
  - 当前所有尚未冻结的模型设计问题。

## 4. 实验

- `docs/experiments/EXPERIMENT_PLAN.md`
  - 仿真、消融、对比、实证的计划和优先级。

后续建议目录：

- `experiments/simulations/`
- `experiments/ablations/`
- `experiments/benchmarks/`
- `experiments/case_studies/`
- `results/`

## 5. 决策记录

- `docs/decisions/ADR-0001-research-repository-scope.md`
- 后续每次改变论文方法定义时新增 ADR，不覆盖历史。

## 6. 测试

- `tests/test_gwr.py`
- `tests/test_grgwr.py`

测试目标不是成熟软件包 API 覆盖，而是防止研究过程中的算法修改破坏已知数学性质。

## 7. 文档优先级

发生冲突时：

1. 最新已接受 ADR
2. `docs/project/CURRENT_STATUS.md`
3. 当前正式模型设计规范
4. `00_PROJECT_HANDOFF.md`
5. README

README 只负责项目入口，不作为算法真相来源。
