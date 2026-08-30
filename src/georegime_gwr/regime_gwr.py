"""Unified regime-aware GWR refit used by current GR-GWR experiments.

This module contains the current *experimental refit primitive* for GR-GWR.
It keeps all observations in one model fit, but a focal location may only
borrow observations from its own regime. Adaptive distance scales are defined
within that regime rather than on the full sample.

This is not yet the final paper algorithm. In particular, K selection,
partition learning, complexity control, and the final bandwidth policy remain
open research questions.
"""

from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
from typing import Mapping

import numpy as np
from scipy.spatial.distance import cdist

from .gwr import BasicGWR


@dataclass
class RegimeAwareGWRResult:
    regimes: np.ndarray
    parameters: np.ndarray
    fitted_values: np.ndarray
    residuals: np.ndarray
    hat_matrix: np.ndarray
    local_bandwidths: np.ndarray


class RegimeAwareGWR:
    """一次统一拟合的 regime-aware GWR。

    对位置 i 与观测 j：

        w_ij = K(d_ij / b_i) * I(z_i == z_j)

    其中 adaptive 距离尺度 ``b_i`` 只使用 i 所在 regime 内的距离定义。

    Parameters
    ----------
    bandwidth : "regime_size", int, or mapping
        ``"regime_size"`` 表示每个 regime 的 adaptive k 等于该 regime
        样本量；整数表示所有 regime 使用相同的区内 neighbour count（若
        大于某区样本量则截到该区样本量）；mapping 则允许显式给出
        ``{regime_label: k}``。
    kernel : {"bisquare", "gaussian", "exponential"}
    fit_intercept : bool

    Notes
    -----
    当前 Georgia K=6 基准使用 ``bandwidth="regime_size"``。它已验证与
    之前六个独立 regime GWR（六区带宽均顶到各自样本量）达到机器精度
    等价；因此这里只改变统一模型表述/实现，不降低统计复杂度。
    """

    def __init__(self, bandwidth="regime_size", kernel="bisquare", fit_intercept=True):
        if kernel not in {"bisquare", "gaussian", "exponential"}:
            raise ValueError("kernel must be bisquare, gaussian, or exponential")
        if isinstance(bandwidth, str) and bandwidth != "regime_size":
            raise ValueError("string bandwidth must be 'regime_size'")
        if not isinstance(bandwidth, (str, Integral, Mapping)):
            raise TypeError("bandwidth must be 'regime_size', int, or mapping")
        if isinstance(bandwidth, Integral) and int(bandwidth) < 1:
            raise ValueError("adaptive bandwidth must be >= 1")

        self.bandwidth = bandwidth
        self.kernel = kernel
        self.fit_intercept = bool(fit_intercept)

    def _bandwidth_map(self, regimes: np.ndarray) -> dict[int, int]:
        labels, counts = np.unique(regimes, return_counts=True)
        sizes = {int(r): int(n) for r, n in zip(labels, counts)}

        if self.bandwidth == "regime_size":
            return sizes.copy()
        if isinstance(self.bandwidth, Integral):
            k = int(self.bandwidth)
            return {r: min(k, n) for r, n in sizes.items()}

        result = {}
        for r, n in sizes.items():
            if r not in self.bandwidth:
                raise ValueError(f"missing bandwidth for regime {r}")
            k = int(self.bandwidth[r])
            if k < 1:
                raise ValueError(f"bandwidth for regime {r} must be >= 1")
            result[r] = min(k, n)
        return result

    def _weights(self, distances_i: np.ndarray, same: np.ndarray, k: int) -> np.ndarray:
        d_same = distances_i[same]
        k_eff = min(int(k), d_same.size)
        if k_eff < 1:
            raise ValueError("within-regime adaptive bandwidth must be >= 1")

        bw = float(np.partition(d_same, k_eff - 1)[k_eff - 1])
        if bw <= 1e-12:
            positive = d_same[d_same > 1e-12]
            bw = float(np.min(positive)) if positive.size else 1.0
        # Match the frozen PyGWRx adaptive boundary policy.
        bw = float(np.nextafter(bw, np.inf))

        ratio = d_same / bw
        if self.kernel == "bisquare":
            w_same = np.where(ratio < 1.0, (1.0 - ratio**2) ** 2, 0.0)
        elif self.kernel == "gaussian":
            w_same = np.exp(-0.5 * ratio**2)
        else:
            w_same = np.exp(-ratio)

        weights = np.zeros_like(distances_i, dtype=float)
        weights[same] = w_same
        return weights

    def fit(self, X, y, coords, regimes):
        X = np.asarray(X, dtype=float)
        y = np.asarray(y, dtype=float).reshape(-1)
        coords = np.asarray(coords, dtype=float)
        regimes = np.asarray(regimes).reshape(-1)

        if X.ndim == 1:
            X = X[:, None]
        if coords.ndim != 2 or coords.shape[1] != 2:
            raise ValueError("coords must have shape (n, 2)")
        n = y.size
        if X.shape[0] != n or coords.shape[0] != n or regimes.size != n:
            raise ValueError("X, y, coords, and regimes must have the same rows")

        Xd = np.column_stack([np.ones(n), X]) if self.fit_intercept else X.copy()
        distances = cdist(coords, coords)
        bandwidths = self._bandwidth_map(regimes)

        p = Xd.shape[1]
        parameters = np.empty((n, p), dtype=float)
        fitted = np.empty(n, dtype=float)
        hat = np.zeros((n, n), dtype=float)
        local_bandwidths = np.empty(n, dtype=int)

        for i in range(n):
            label = int(regimes[i])
            same = regimes == regimes[i]
            k_i = int(bandwidths[label])
            weights = self._weights(distances[i], same, k_i)
            beta, C = BasicGWR._solve_local(Xd, y, weights)
            parameters[i] = beta
            fitted[i] = Xd[i] @ beta
            hat[i] = Xd[i] @ C
            local_bandwidths[i] = k_i

        residuals = y - fitted
        self.X_ = X
        self.X_design_ = Xd
        self.y_ = y
        self.coords_ = coords
        self.regimes_ = regimes.copy()
        self.distance_matrix_ = distances
        self.bandwidths_ = bandwidths
        self.local_bandwidths_ = local_bandwidths
        self.parameters_ = parameters
        self.fitted_values_ = fitted
        self.residuals_ = residuals
        self.hat_matrix_ = hat
        self.result_ = RegimeAwareGWRResult(
            regimes=self.regimes_.copy(),
            parameters=parameters.copy(),
            fitted_values=fitted.copy(),
            residuals=residuals.copy(),
            hat_matrix=hat.copy(),
            local_bandwidths=local_bandwidths.copy(),
        )
        return self
