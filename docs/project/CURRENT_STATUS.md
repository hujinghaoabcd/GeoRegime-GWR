# CURRENT_STATUS

Last updated: 2026-08-30

## Current phase

**Phase 1 — GR-GWR component research**

标准 GWR 基线已经完成端到端验证并封存。当前正在逐组件研究 GR-GWR。已经完成 Georgia Queen W、pilot GWR local coefficients、local relationship edge diagnostics、Queen-constrained Ward 初始分区、K=6 fixed-label regime-restricted GWR、MGWR benchmark 与 K sensitivity。

## Completed

- 建立独立仓库 `GeoRegime-GWR`。
- 加入最小 `BasicGWR` 与 `GRGWRBaseline`；后者仍是旧 baseline snapshot，不代表最终论文算法。
- 建立 `00_PROJECT_HANDOFF.md`、CURRENT_STATUS、设计文档、ADR、实验计划和 CI。
- 引入 canonical Georgia 1990 county benchmark（159 counties）。
- Standard GWR external validation 已通过：给定相同 bandwidth=117 时，仓库 `BasicGWR` 与 external `mgwr.gwr.GWR` 的局部参数、拟合值、残差、hat matrix 达到机器精度一致。
- Adaptive research default 使用 exhaustive integer AICc，Georgia k=116，AICc=298.9855995813665。
- Historical mgwr compatibility mode k=117，AICc=299.0508086830287。
- Fixed-distance golden search 已实现。
- Georgia Queen contiguity W：159 counties、431 undirected edges、无孤立县、整体连通。
- Pilot BasicGWR local coefficient maps 已生成。
- 当前 exploratory relationship fingerprint：`PctFB`, `PctBlack`, `PctRural` 三个标准化 local slopes；不含 intercept，不含 coordinates。
- 431 条 Queen edges 的 coefficient-fingerprint Euclidean distance 已计算。
- Queen-constrained Ward `K=2..15` 初始分区已生成，所有候选 regime 空间连通。
- Focused K=6 初始分区已检查：sizes = 53,19,17,22,27,21；WCSS=41.17628003831048。
- Fixed K=6 regime-restricted GWR 已完成：每个 regime 只借本区样本，cross-regime weights=0；无 label refinement。
- K=6 每区独立 bandwidth 全部顶到 regime size：53,19,17,22,27,21。
- Fixed K=6 模型形式对照已完成：ordinary GWR vs K6 regime OLS vs K6 regime GWR。
- External MGWR benchmark 已完成（`mgwr==2.2.1`）。
- `K=2..15` fixed-partition regime-GWR sensitivity 已完成；当前严格规格可稳定完成 K=2..6。

## Standard GWR bandwidth policy — FROZEN BY ADR-0003

### Adaptive research default

`BasicGWR(bandwidth="auto")`

- exhaustive integer AICc search；
- Georgia selected k=116；
- AICc=298.9855995813665。

### Fixed-distance research default

`BasicGWR(bandwidth="auto", adaptive=False)`

- continuous golden-section AICc search；
- 功能已实现，尚未建立独立 external fixed-GWR benchmark。

### Historical mgwr reproduction mode

`BasicGWR(bandwidth="auto", adaptive=True, search_strategy="mgwr_golden")`

- Georgia selected k=117；
- 与 external `mgwr.gwr.GWR` 数值一致。

正式决策：`docs/decisions/ADR-0003-standard-gwr-bandwidth-search-policy.md`。

## Current exploratory component choices — NOT FROZEN

1. D 与 W 分离：D 负责 GWR kernel distance，W 负责 topology / allowed merges；
2. Georgia polygon W 当前使用 Queen contiguity；
3. ordinary GWR 只作为 pilot local-relationship estimator；
4. fingerprint 暂用三个标准化 local slopes；
5. intercept 暂不进入 fingerprint；
6. coordinates 暂不进入 fingerprint；
7. similarity 暂用 standardized-slope Euclidean distance；
8. initial segmentation 暂用 Queen-constrained Ward；
9. K=2..15 仅为 exploratory scan；
10. 当前继续使用 K=6 作为 working candidate，**不是最终 K**；
11. fixed K=6 后每个 regime 暂独立自动选 bandwidth；
12. AICc、lambda、ICM/refinement、最终 bandwidth policy 均未冻结。

## Current empirical results

### Ordinary GWR / MGWR / K=6 restricted GWR

- ordinary GWR: RMSE=0.56681, MAE=0.39943, R2=0.67872, trace(S)=11.9121, AICc=298.99；
- MGWR: RMSE=0.56579, MAE=0.40098, R2=0.67988, trace(S)=11.3683, AICc=297.12；
- K=6 regime-restricted GWR: RMSE=0.52064, MAE=0.34497, R2=0.72893, trace(S)=39.7019, conditional AICc=354.01。

MGWR bandwidths：Intercept=92, PctFB=101, PctBlack=136, PctRural=158。

这些目前都是 exploratory / in-sample evidence，不能据此宣称 GR-GWR 已优于 MGWR。

### K sensitivity

当前 fixed-partition -> per-regime GWR：

- K=2: RMSE=0.57225, R2=0.67253；
- K=3: RMSE=0.55129, R2=0.69608；
- K=4: RMSE=0.54037, R2=0.70800；
- K=5: RMSE=0.53018, R2=0.71891；
- K=6: RMSE=0.52064, R2=0.72893。

K=7..9 出现 size=7 regime，当前严格 bandwidth search 无法稳定完成；K=10..15 出现 size=4 regime，当前 4 参数局部回归明显不适合作为稳定 regime。

当前只说明 K=6 是可稳定估计候选中样本内拟合最好、且没有极小 regime 的 working candidate；不代表最终 K 已选定。

## Pending method questions — MUST RESOLVE LATER

完整记录：`docs/project/PENDING_GRGWR_METHOD_QUESTIONS.md`

当前重点 OPEN 问题包括：

- pilot GWR boundary smoothing 与真实边界恢复能力；
- K selection 与 Kmax/minimum regime size；
- in-sample K 增大带来的 data-adaptive advantage；
- per-regime vs shared bandwidth；
- 所有 K=5/6 regime bandwidth 顶到最大值的含义；
- MGWR multiscale heterogeneity 与 GR-GWR discontinuous/regime heterogeneity 的正式比较；
- AICc / effective df / partition-selection complexity；
- Ward 是否只是 initializer；
- refinement / ICM / lambda 如何定义；
- spatial / blocked CV、partition stability、synthetic known-boundary validation；
- intercept / coordinates / uncertainty-aware fingerprint 等消融。

这些问题全部保持 OPEN，不在当前阶段强行解决。

## Current algorithm path

为了先看完整算法行为，当前不重构总体思路，继续既有研究链：

`ordinary pilot GWR -> Queen-constrained Ward initial regimes -> fixed-regime restricted GWR -> label refinement -> restricted GWR refit -> iterate / converge`

当前 K=6 仅作为 working candidate 用于跑通这条完整链。

## Current highest-priority task

**在当前 K=6 working partition 上进入第一轮 label refinement。**

要求：

- Queen W 继续作为邻接与连通性约束；
- 不允许 source regime 因单点移动而被切断；
- target label 只从当前 label + Queen 邻居 labels 中选择；
- 不允许产生过小/不可估计 regime；
- refinement 后重新做 regime-restricted GWR；
- 记录 label changes、regime sizes、connectivity、boundary edges、RSS/RMSE/MAE/R2 和 objective change；
- 当前 lambda / objective 仍只是 exploratory baseline，不冻结为论文最终定义。

## Next tasks

1. K=6 第一轮 label refinement；
2. refinement 后 restricted GWR refit；
3. 若仍有有效 label change，则按既有思路继续迭代到稳定；
4. 记录完整 convergence path；
5. 完整链跑通后，再系统解决 K、bandwidth、lambda、AICc、CV、inference 和 synthetic boundary validation。
