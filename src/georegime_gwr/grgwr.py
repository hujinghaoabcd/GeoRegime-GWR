"""Research baseline of Geo-Regime Geographically Weighted Regression.

IMPORTANT
---------
This file is a deliberately readable baseline extracted from the current
pyGWRx GR-GWR idea.  It is NOT the final paper algorithm.  Future experiments
may replace the spatial weights construction, clustering, ICM refinement, or
objective function entirely.  Every major change should be recorded in docs.
"""

from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial.distance import cdist
from sklearn.cluster import AgglomerativeClustering

from .gwr import BasicGWR


@dataclass
class GRGWRResult:
    regimes: np.ndarray
    parameters: np.ndarray
    fitted_values: np.ndarray
    residuals: np.ndarray
    objective_history: list[float]


class GRGWRBaseline:
    """当前 pyGWRx 思路的轻量研究基线。

    当前流程：
    1. kNN + MST 建空间邻接图；
    2. 全域普通 GWR 得到每个位置的局部参数；
    3. 去掉截距、标准化斜率，并与标准化坐标组合为聚类特征；
    4. 空间连通约束 Ward 聚类得到初始 regime；
    5. 仅在同 regime 内重新做局部 GWR；
    6. 可选 ICM 风格逐点 refinement；
    7. 用 RSS + graph-cut penalty 做整轮 objective guard。

    注意：以上每一步都允许在论文研究过程中被替换或删除。
    """

    def __init__(
        self,
        n_regimes=3,
        bandwidth=20,
        kernel="bisquare",
        lambda_boundary=1.0,
        max_iter=10,
        tol=1e-4,
        spatial_constraint_weight=0.5,
        fit_intercept=True,
        n_neighbors=8,
        min_regime_size=None,
        enforce_connectivity=True,
        random_state=42,
        refine=True,
    ):
        self.n_regimes = int(n_regimes)
        self.bandwidth = bandwidth
        self.kernel = kernel
        self.lambda_boundary = float(lambda_boundary)
        self.max_iter = int(max_iter)
        self.tol = float(tol)
        self.spatial_constraint_weight = float(spatial_constraint_weight)
        self.fit_intercept = bool(fit_intercept)
        self.n_neighbors = int(n_neighbors)
        self.min_regime_size = min_regime_size
        self.enforce_connectivity = bool(enforce_connectivity)
        self.random_state = random_state
        self.refine = bool(refine)

        if not 0.0 <= self.spatial_constraint_weight <= 1.0:
            raise ValueError("spatial_constraint_weight must be in [0, 1]")
        if self.n_regimes < 1:
            raise ValueError("n_regimes must be >= 1")

        self._gwr = BasicGWR(bandwidth, kernel, fit_intercept)

    # ------------------------------------------------------------------
    # 1. 当前基线的空间邻接结构：kNN + MST
    # ------------------------------------------------------------------
    def _build_graph(self, coords: np.ndarray) -> None:
        n = coords.shape[0]
        self.distance_matrix_ = cdist(coords, coords)
        k = min(self.n_neighbors, n - 1)
        graph = np.zeros((n, n), dtype=bool)

        for i in range(n):
            neighbours = np.argsort(self.distance_matrix_[i])[1 : k + 1]
            graph[i, neighbours] = True
        graph |= graph.T

        # 当前基线沿用 pyGWRx：MST 仅用于保证图整体连通。
        mst = minimum_spanning_tree(self.distance_matrix_).tocoo()
        for i, j in zip(mst.row, mst.col):
            graph[int(i), int(j)] = True
            graph[int(j), int(i)] = True

        np.fill_diagonal(graph, False)
        self.adjacency_matrix_ = csr_matrix(graph.astype(float))
        self.adjacency_ = tuple(np.flatnonzero(graph[i]) for i in range(n))
        self.edges_ = tuple(
            (int(i), int(j))
            for i in range(n)
            for j in self.adjacency_[i]
            if i < j
        )

    # ------------------------------------------------------------------
    # 2. 第一轮普通 GWR
    # ------------------------------------------------------------------
    def _fit_global_gwr(self) -> np.ndarray:
        n, p = self.X_design_.shape
        parameters = np.zeros((n, p))
        for i in range(n):
            parameters[i], _ = BasicGWR._solve_local(
                self.X_design_,
                self.y_,
                self._gwr._weights(self.distance_matrix_[i]),
            )
        return parameters

    # ------------------------------------------------------------------
    # 3. 分区特征
    # ------------------------------------------------------------------
    @staticmethod
    def _safe_standardize(values: np.ndarray) -> np.ndarray:
        mean = np.mean(values, axis=0)
        scale = np.std(values, axis=0, ddof=0)
        scale = np.where(scale <= 1e-12, 1.0, scale)
        return (values - mean) / scale

    def _clustering_features(self, local_parameters, coords):
        slopes = local_parameters[:, 1:] if self.fit_intercept else local_parameters
        slopes_scaled = self._safe_standardize(slopes)

        coord_range = np.ptp(coords, axis=0)
        coord_range = np.where(coord_range <= 1e-12, 1.0, coord_range)
        coords_scaled = (coords - np.min(coords, axis=0)) / coord_range

        gamma = self.spatial_constraint_weight
        if gamma == 0.0:
            return slopes_scaled
        if gamma == 1.0:
            return coords_scaled
        return np.hstack(
            [np.sqrt(1.0 - gamma) * slopes_scaled,
             np.sqrt(gamma) * coords_scaled]
        )

    @staticmethod
    def _relabel(labels):
        unique = np.unique(labels)
        mapping = {int(v): i for i, v in enumerate(unique)}
        return np.asarray([mapping[int(v)] for v in labels], dtype=int)

    # ------------------------------------------------------------------
    # 4. 初始机制区：空间连通约束 Ward
    # ------------------------------------------------------------------
    def _initial_regimes(self, features):
        n = features.shape[0]
        feasible = max(1, n // self._min_regime_size_)
        k = min(self.n_regimes, feasible)
        if k == 1:
            return np.zeros(n, dtype=int)

        model = AgglomerativeClustering(
            n_clusters=k,
            linkage="ward",
            connectivity=self.adjacency_matrix_,
            compute_full_tree=True,
        )
        labels = model.fit_predict(features)
        return self._merge_small_regimes(self._relabel(labels), features)

    def _merge_small_regimes(self, labels, features):
        labels = self._relabel(labels.copy())
        while True:
            counts = np.bincount(labels)
            small = np.flatnonzero(counts < self._min_regime_size_)
            if small.size == 0 or counts.size == 1:
                return labels

            source = int(small[np.argmin(counts[small])])
            nodes = np.flatnonzero(labels == source)
            adjacent_labels = set()
            for node in nodes:
                adjacent_labels.update(int(labels[j]) for j in self.adjacency_[node])
            adjacent_labels.discard(source)
            candidates = sorted(adjacent_labels) or [
                r for r in range(counts.size) if r != source
            ]
            center = np.mean(features[nodes], axis=0)
            target = min(
                candidates,
                key=lambda r: float(
                    np.linalg.norm(center - np.mean(features[labels == r], axis=0))
                ),
            )
            labels[nodes] = target
            labels = self._relabel(labels)

    # ------------------------------------------------------------------
    # 5. 给定 regime 后，只在同 regime 内做局部 GWR
    # ------------------------------------------------------------------
    def _fit_for_labels(self, labels):
        n, p = self.X_design_.shape
        parameters = np.zeros((n, p))
        fitted = np.zeros(n)

        for regime in range(int(np.max(labels)) + 1):
            indices = np.flatnonzero(labels == regime)
            Xr = self.X_design_[indices]
            yr = self.y_[indices]
            for global_i in indices:
                d = self.distance_matrix_[global_i, indices]
                beta, _ = BasicGWR._solve_local(Xr, yr, self._gwr._weights(d))
                parameters[global_i] = beta
                fitted[global_i] = self.X_design_[global_i] @ beta

        return parameters, fitted

    # ------------------------------------------------------------------
    # 6-7. ICM 风格 refinement 与全局 objective guard
    # ------------------------------------------------------------------
    def _boundary_count(self, labels):
        return int(sum(labels[i] != labels[j] for i, j in self.edges_))

    def _objective(self, labels, fitted):
        rss = float(np.sum((self.y_ - fitted) ** 2))
        return rss + self.lambda_boundary * self._boundary_count(labels)

    def _candidate_error(self, node, regime, labels):
        indices = np.flatnonzero(labels == regime)
        indices = indices[indices != node]
        if indices.size < self.X_design_.shape[1]:
            return float("inf")
        beta, _ = BasicGWR._solve_local(
            self.X_design_[indices],
            self.y_[indices],
            self._gwr._weights(self.distance_matrix_[node, indices]),
        )
        pred = float(self.X_design_[node] @ beta)
        return float((self.y_[node] - pred) ** 2)

    def _can_remove(self, node, labels):
        source = int(labels[node])
        members = np.flatnonzero(labels == source)
        if members.size - 1 < self._min_regime_size_:
            return False
        if not self.enforce_connectivity:
            return True

        remaining = members[members != node]
        if remaining.size <= 1:
            return True
        allowed = set(int(v) for v in remaining)
        visited = {int(remaining[0])}
        stack = [int(remaining[0])]
        while stack:
            current = stack.pop()
            for neighbour in self.adjacency_[current]:
                j = int(neighbour)
                if j in allowed and j not in visited:
                    visited.add(j)
                    stack.append(j)
        return len(visited) == remaining.size

    def _icm_sweep(self, labels, iteration):
        updated = labels.copy()
        seed = (0 if self.random_state is None else int(self.random_state)) + iteration
        order = np.random.default_rng(seed).permutation(labels.size)
        changed = 0

        for node in map(int, order):
            current = int(updated[node])
            if not self._can_remove(node, updated):
                continue
            candidates = set(int(updated[j]) for j in self.adjacency_[node]) | {current}
            costs = {}
            for r in candidates:
                disagreement = int(np.sum(updated[self.adjacency_[node]] != r))
                costs[r] = (
                    self._candidate_error(node, r, updated)
                    + self.lambda_boundary * disagreement
                )
            best = min(costs, key=costs.get)
            if best != current and costs[best] < costs[current] - self.tol:
                updated[node] = best
                changed += 1
        return updated, changed

    # ------------------------------------------------------------------
    # Public fit
    # ------------------------------------------------------------------
    def fit(self, X, y, coords):
        X = np.asarray(X, dtype=float)
        y = np.asarray(y, dtype=float).reshape(-1)
        coords = np.asarray(coords, dtype=float)
        if X.ndim == 1:
            X = X[:, None]
        if X.shape[0] != y.size or coords.shape[0] != y.size:
            raise ValueError("X, y, and coords must have the same rows")

        self.X_ = X
        self.y_ = y
        self.coords_ = coords
        self.X_design_ = (
            np.column_stack([np.ones(X.shape[0]), X]) if self.fit_intercept else X.copy()
        )
        p = self.X_design_.shape[1]
        self._min_regime_size_ = p + 2 if self.min_regime_size is None else int(self.min_regime_size)

        self._build_graph(coords)
        self.initial_gwr_parameters_ = self._fit_global_gwr()
        self.clustering_features_ = self._clustering_features(
            self.initial_gwr_parameters_, coords
        )
        labels = self._initial_regimes(self.clustering_features_)

        parameters, fitted = self._fit_for_labels(labels)
        objective = self._objective(labels, fitted)
        self.objective_history_ = [objective]

        if self.refine:
            for iteration in range(self.max_iter):
                proposed, changed = self._icm_sweep(labels, iteration)
                if changed == 0:
                    break
                p_params, p_fitted = self._fit_for_labels(proposed)
                p_objective = self._objective(proposed, p_fitted)
                if p_objective > objective + self.tol:
                    break
                labels, parameters, fitted, objective = (
                    proposed, p_params, p_fitted, p_objective
                )
                self.objective_history_.append(objective)

        self.regimes_ = self._relabel(labels)
        # 重新按最终标签拟合，确保结果与最终 regime 一致。
        self.parameters_, self.fitted_values_ = self._fit_for_labels(self.regimes_)
        self.residuals_ = self.y_ - self.fitted_values_
        self.result_ = GRGWRResult(
            regimes=self.regimes_.copy(),
            parameters=self.parameters_.copy(),
            fitted_values=self.fitted_values_.copy(),
            residuals=self.residuals_.copy(),
            objective_history=list(self.objective_history_),
        )
        return self
