"""Minimal standard GWR used as the baseline engine for GR-GWR research."""

from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral

import numpy as np
from scipy.spatial.distance import cdist

from .bandwidth import (
    FixedGoldenAICcSelector,
    MGWRCompatibleAICcSelector,
    PyGWRxAdaptiveAICcSelector,
)


@dataclass
class GWRResult:
    parameters: np.ndarray
    fitted_values: np.ndarray
    residuals: np.ndarray
    hat_matrix: np.ndarray


class BasicGWR:
    """轻量级标准 GWR。

    Parameters
    ----------
    bandwidth : int, float, or "auto"
        Integer means adaptive neighbour-order bandwidth; float means fixed
        distance. ``"auto"`` activates AICc bandwidth selection.
    kernel : {"bisquare", "gaussian", "exponential"}
    fit_intercept : bool
    adaptive : bool or None, optional
        For ``bandwidth="auto"``, the default ``None`` resolves to ``True``.
        For numeric bandwidths, ``None`` preserves the historical convention:
        integer -> adaptive, float -> fixed.
    search_strategy : str or None, optional
        Automatic-search policy. Defaults to ``"exhaustive"`` for adaptive
        bandwidths and ``"golden_section"`` for fixed bandwidths.

        Adaptive choices:
        - ``"exhaustive"``: current PyGWRx integer exhaustive AICc search;
        - ``"mgwr_golden"``: mgwr 2.2.1 compatibility search.

        Fixed choice:
        - ``"golden_section"``: continuous PyGWRx-style golden-section AICc search.

    Notes
    -----
    The research default is deliberately the strict adaptive exhaustive search.
    ``mgwr_golden`` exists only for exact historical/canonical reproduction.
    """

    def __init__(
        self,
        bandwidth="auto",
        kernel="bisquare",
        fit_intercept=True,
        *,
        adaptive=None,
        search_strategy=None,
    ):
        if kernel not in {"bisquare", "gaussian", "exponential"}:
            raise ValueError("kernel must be bisquare, gaussian, or exponential")
        if isinstance(bandwidth, str) and bandwidth.lower() != "auto":
            raise ValueError("string bandwidth must be 'auto'")
        if adaptive is not None and not isinstance(adaptive, (bool, np.bool_)):
            raise TypeError("adaptive must be bool or None")

        self.bandwidth = bandwidth
        self.kernel = kernel
        self.fit_intercept = bool(fit_intercept)
        self.adaptive = None if adaptive is None else bool(adaptive)
        self.search_strategy = search_strategy

    def _resolve_adaptive(self) -> bool:
        if self.adaptive is not None:
            return self.adaptive
        if isinstance(self.bandwidth, str):
            return True
        return isinstance(self.bandwidth, Integral) and not isinstance(self.bandwidth, (bool, np.bool_))

    def _resolve_search_strategy(self, adaptive: bool) -> str | None:
        if not isinstance(self.bandwidth, str):
            return None
        if self.search_strategy is None:
            return "exhaustive" if adaptive else "golden_section"
        strategy = str(self.search_strategy).strip().lower()
        if adaptive and strategy not in {"exhaustive", "mgwr_golden"}:
            raise ValueError("adaptive search_strategy must be 'exhaustive' or 'mgwr_golden'")
        if not adaptive and strategy != "golden_section":
            raise ValueError("fixed search_strategy must be 'golden_section'")
        return strategy

    def _active_bandwidth(self):
        return getattr(self, "bandwidth_", self.bandwidth)

    def _weights(self, distances: np.ndarray) -> np.ndarray:
        bandwidth = self._active_bandwidth()
        if self.adaptive_:
            if isinstance(bandwidth, (bool, np.bool_)) or not float(bandwidth).is_integer():
                raise ValueError("adaptive bandwidth must be an integer neighbour count")
            k = min(int(bandwidth), distances.size)
            if k < 1:
                raise ValueError("adaptive bandwidth must be >= 1")
            bw = float(np.partition(distances, k - 1)[k - 1])
            if bw <= 1e-12:
                positive = distances[distances > 1e-12]
                bw = float(np.min(positive)) if positive.size else 1.0
            if self.boundary_policy_ == "mgwr":
                bw *= 1.0000001
            else:
                bw = float(np.nextafter(bw, np.inf))
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
        self.adaptive_ = self._resolve_adaptive()
        strategy = self._resolve_search_strategy(self.adaptive_)

        if isinstance(self.bandwidth, str):
            if self.adaptive_ and strategy == "exhaustive":
                selector = PyGWRxAdaptiveAICcSelector(kernel=self.kernel)
                self.boundary_policy_ = "pygwrx"
            elif self.adaptive_ and strategy == "mgwr_golden":
                selector = MGWRCompatibleAICcSelector(kernel=self.kernel)
                self.boundary_policy_ = "mgwr"
            elif not self.adaptive_ and strategy == "golden_section":
                selector = FixedGoldenAICcSelector(kernel=self.kernel)
                self.boundary_policy_ = "fixed"
            else:
                raise RuntimeError("unsupported bandwidth search configuration")

            search = selector.select(Xd, y, coords)
            self.bandwidth_selector_ = selector
            self.bandwidth_search_ = search
            self.bandwidth_ = int(search.bandwidth) if self.adaptive_ else float(search.bandwidth)
            self.search_strategy_ = strategy
        else:
            if self.adaptive_:
                value = float(self.bandwidth)
                if not value.is_integer():
                    raise ValueError("adaptive numeric bandwidth must be an integer")
                self.bandwidth_ = int(value)
                self.boundary_policy_ = "pygwrx"
            else:
                self.bandwidth_ = float(self.bandwidth)
                if self.bandwidth_ <= 0.0:
                    raise ValueError("fixed numeric bandwidth must be > 0")
                self.boundary_policy_ = "fixed"
            self.bandwidth_selector_ = None
            self.bandwidth_search_ = None
            self.search_strategy_ = None

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
