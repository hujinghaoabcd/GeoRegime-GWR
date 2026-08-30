# RESEARCH_QUESTIONS

本文件只记录**尚未冻结**的问题。每个问题最终都应通过：理论依据 + 文献 + 仿真 + 消融决定，而不是凭直觉拍板。

## RQ1. 空间拓扑应如何定义？

候选：

- polygon: Queen / Rook contiguity；
- point: kNN / Delaunay / distance threshold；
- network: network adjacency；
- user-supplied spatial weights matrix `W`。

核心问题：正式模型是否应直接以 `W` 为输入，而把自动构图降级为便利功能？

## RQ2. 是否还需要把坐标放进 clustering feature？

当前有两重空间约束：

1. feature 中的坐标；
2. clustering connectivity / regime connectivity。

需要比较：

- coefficients only；
- coefficients + coordinates；
- coefficients only + spatial adjacency constraint。

## RQ3. Ward 是否必要？

当前 Ward 只承担初始化。候选：

- 保留 constrained Ward；
- 直接从区域化算法初始化；
- 多随机初始化；
- 从 coefficient change-point / graph partition 方法初始化；
- 不显式聚类，改为统一目标函数联合估计。

## RQ4. ICM refinement 是否必要？

比较：

- `refine=False`；
- 当前 sequential ICM；
- 其他离散优化；
- 完全联合优化。

如果 ICM 对 boundary recovery / coefficient RMSE / prediction 几乎无增益，应删除以保持方法简洁。

## RQ5. hard boundary 是否过强？

当前：

`z_i != z_j => borrowing weight = 0`

候选：

- hard gating；
- soft cross-boundary attenuation；
- transition-zone model。

## RQ6. 是否允许不同 covariate 拥有不同 regime？

当前 shared regime：

`z_i` 对所有 `beta_k(s)` 相同。

未来可能扩展：

`z_ik`，即变量特异的 regime map。

第一篇论文是否保持 shared regime，取决于 conflicting-boundary simulation。

## RQ7. regime 数 K 如何确定？

当前人工指定 `n_regimes`。

候选：

- conditional AICc；
- spatial CV；
- stability；
- penalized objective；
- hierarchical stopping rule。

## RQ8. 如何避免第一轮 GWR 边界污染？

可能方案：

- 保持当前做法并证明后续步骤可恢复；
- 用更小/多尺度 bandwidth；
- robust coefficient fingerprints；
- bootstrap/stability filtering；
- 完全放弃 GWR-first initialization。

## RQ9. 边界复杂度如何惩罚？

当前是 graph-cut edge count：

`sum I(z_i != z_j)`。

候选：

- unweighted graph cut；
- edge-length weighted cut；
- shared-boundary-length penalty；
- Potts penalty；
- no explicit penalty + connectivity only。

## RQ10. 论文真正的核心主张是什么？

当前最有希望的主张不是“首次发现 spatial regime”，而是：

**让局部关系估计与内生机制区发现形成反馈耦合，从而学习 GWR 应当在哪里停止跨空间边界借样本。**

这一主张仍需通过实验确认。
