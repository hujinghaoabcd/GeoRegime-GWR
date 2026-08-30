"""Minimal GWR used as the baseline engine for GR-GWR research.

This module is intentionally small. It keeps only the mechanics needed by
paper experiments: Euclidean distance, fixed/adaptive bandwidths, three common
kernels, local weighted least squares, automatic adaptive AICc bandwidth
selection, and fitted local coefficients.
"""

from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral

import numpy as np
from scipy.spatial.distance import cdist

from .bandwidth import MGWRCompatibleAICcSelector


@dataclass
class GWRResult:
    parameters: np.ndarray
    fitted_values: np.ndarray
    residuals: np.ndarray
    hat_matrix: np.ndarray


class BasicGWR:
    """轻量级基础 GWR。

    Parameters
    ----------
    bandwidth : int, float, or "auto"
        Integer = adaptive neighbour-order bandwidth; float = fixed distance;
        ``"auto"`` = adaptive AICc search reproducing the standard-GWR search
        behavior of ``mgwr==2.2.1`` used by the canonical Georgia benchmark.
    kernel : {"bisquare", "gaussian", "exponential"}
    fit_intercept : bool

    Notes
    -----
    这里只服务于 GR-GWR 方法研究，不追求 pyGWRx 的完整公共 API。
    """

    def __init__(self, bandwidth="auto", kernel="bisquare", fit_intercept=True):
        if kernel not in {"bisquare", "gaussian", "exponential"}:
            raise ValueError("kernel must be bisquare, gaussian, or exponential")
        if isinstance(bandwidth, str) and bandwidth.lower() != "auto":
            raise ValueError("string bandwidth must be 'auto'")
        self.bandwidth = bandwidth
        self.kernel = kernel
        self.fit_intercept = fit_intercept

    def _active_bandwidth(self):
        return getattr(self, "bandwidth_", self.bandwidth)

    def _weights(self, distances: np.ndarray) -> np.ndarray:
        bandwidth = self._active_bandwidth()
        if isinstance(bandwidth, Integral):
            k = min(int(bandwidth), distances.size)
            if k < 1:
                raise ValueError("adaptive bandwidth must be >= 1")
            bw = float(np.partition(distances, k - 1)[k - 1])
            if bw <= 1e-12:
                positive = distances[distances > 1e-12]
                bw = float(np.min(positive)) if positive.size else 1.0
            # mgwr 2.2.1 compact-kernel boundary convention.
            bw *= 1.0000001
        else:
            bw = float(bandwidth)
            if bw <= 0:
                raise ValueError("fixed bandwidth must be > 0")

        ratio = distances / bw
        if self.kernel == "bisquare":
            return np.where(ratio < 1.0, (1.0 - ratio**2) ** 2, 0.0)
        if self.kernel == "gaussian":
            return np.exp(-0.5 * ratio**2)
        return np.exp(-ratio)

    @staticmethod
    def _solve_local(X: np.ndarray, y: np.ndarray, weights: np.ndarray):
        """局部加权最小二乘；奇异时按 solve -> ridge -> pinv 回退。"""
        Xw = X * weights[:, None]
        M = Xw.T @ X
        p = X.shape[1]
        try:
            C = np.linalg.solve(M, Xw.T)
            beta = C @ y
            if np.all(np.isfinite(beta)):
                return beta, C
        except np.linalg.LinAlgError:
            pass

        ridge = 1e-6 * (np.trace(M) / max(p, 1) + 1e-12) + 1e-12
        try:
            C = np.linalg.solve(M + ridge * np.eye(p), Xw.T)
            beta = C @ y
            if np.all(np.isfinite(beta)):
                return beta, C
        except np.linalg.LinAlgError:
            pass

        C = np.linalg.pinv(M) @ Xw.T
        return C @ y, C

    def fit(self, X, y, coords):
        X = np.asarray(X, dtype=float)
        y = np.asarray(y, dtype=float).reshape(-1)
        coords = np.asarray(coords, dtype=float)
        if X.ndim == 1:
            X = X[:, None]
        if coords.ndim != 2 or coords.shape[1] != 2:
            raise ValueError("coords must have shape (n, 2)")
        if X.shape[0] != y.size or coords.shape[0] != y.size:
            raise ValueError("X, y, and coords must have the same rows")

        Xd = np.column_stack([np.ones(X.shape[0]), X]) if self.fit_intercept else X.copy()

        if isinstance(self.bandwidth, str) and self.bandwidth.lower() == "auto":
            selector = MGWRCompatibleAICcSelector(kernel=self.kernel)
            search = selector.select(Xd, y, coords)
            self.bandwidth_selector_ = selector
            self.bandwidth_search_ = search
            self.bandwidth_ = int(search.bandwidth)
        else:
            self.bandwidth_ = self.bandwidth
            self.bandwidth_selector_ = None
            self.bandwidth_search_ = None

        distances = cdist(coords, coords)
        n, p = Xd.shape
        parameters = np.zeros((n, p))
        fitted = np.zeros(n)
        hat = np.zeros((n, n))

        for i in range(n):
            beta, C = self._solve_local(Xd, y, self._weights(distances[i]))
            parameters[i] = beta
            fitted[i] = Xd[i] @ beta
            hat[i] = Xd[i] @ C

        self.X_ = X
        self.X_design_ = Xd
        self.y_ = y
        self.coords_ = coords
        self.distance_matrix_ = distances
        self.parameters_ = parameters
        self.fitted_values_ = fitted
        self.residuals_ = y - fitted
        self.hat_matrix_ = hat
        self.result_ = GWRResult(parameters, fitted, self.residuals_, hat)
        return self
