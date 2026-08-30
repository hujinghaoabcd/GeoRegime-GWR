"""Lightweight bandwidth selection for the paper-research GWR baseline.

Two adaptive AICc search policies are intentionally kept separate:

``PyGWRxAdaptiveAICcSelector``
    A faithful extraction of the current PyGWRx policy: scan every valid
    integer neighbour order and choose the global minimum AICc.

``MGWRCompatibleAICcSelector``
    Reproduce the standard-GWR bandwidth search used by ``mgwr==2.2.1``:
    discrete golden-section search with the same initial adaptive interval and
    compact-kernel boundary convention.  This is the policy used by
    ``BasicGWR(bandwidth='auto')`` because the Georgia benchmark must reproduce
    the published/canonical mgwr GWR result exactly, including bandwidth 117.

The distinction is documented rather than hidden because on Georgia the
exhaustive PyGWRx curve has its global minimum at k=116, while mgwr's historical
golden-section procedure returns k=117.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial.distance import cdist


@dataclass(frozen=True)
class BandwidthSearchResult:
    bandwidth: int
    score: float
    search_range: tuple[int, int]
    search_trace: tuple[tuple[int, float], ...]
    strategy: str


def _distance_bandwidth_pygwrx(distances: np.ndarray, k: int) -> float:
    """Current PyGWRx adaptive distance convention."""
    if k < 1 or k > distances.size:
        raise ValueError(f"k must satisfy 1 <= k <= {distances.size}; got {k}")
    bandwidth = float(np.partition(distances, k - 1)[k - 1])
    if bandwidth <= 0.0:
        positive = distances[distances > 0.0]
        if positive.size == 0:
            raise np.linalg.LinAlgError("all distances are zero")
        bandwidth = float(np.min(positive))
    return float(np.nextafter(bandwidth, np.inf))


def _distance_bandwidth_mgwr(distances: np.ndarray, k: int) -> float:
    """mgwr 2.2.1 adaptive compact-kernel convention.

    mgwr's Kernel class multiplies the k-th ordered distance by 1.0000001 so
    that the boundary observation remains inside compact-support kernels.
    """
    if k < 1 or k > distances.size:
        raise ValueError(f"k must satisfy 1 <= k <= {distances.size}; got {k}")
    bandwidth = float(np.partition(distances, k - 1)[k - 1])
    if bandwidth <= 0.0:
        positive = distances[distances > 0.0]
        if positive.size == 0:
            raise np.linalg.LinAlgError("all distances are zero")
        bandwidth = float(np.min(positive))
    return bandwidth * 1.0000001


def _kernel_from_distance_bandwidth(
    distances: np.ndarray,
    distance_bandwidth: float,
    kernel: str,
) -> np.ndarray:
    ratio = distances / distance_bandwidth
    if kernel == "bisquare":
        return np.where(ratio < 1.0, (1.0 - ratio**2) ** 2, 0.0)
    if kernel == "gaussian":
        return np.exp(-0.5 * ratio**2)
    if kernel == "exponential":
        return np.exp(-ratio)
    raise ValueError("kernel must be bisquare, gaussian, or exponential")


def _fit_local_unpenalized(
    X: np.ndarray,
    y: np.ndarray,
    weights: np.ndarray,
    target_row: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    p = X.shape[1]
    positive = weights > 0.0
    if np.count_nonzero(positive) < p:
        raise np.linalg.LinAlgError("too few positive-weight observations")
    Xw_rank = X[positive] * np.sqrt(weights[positive])[:, None]
    if np.linalg.matrix_rank(Xw_rank) < p:
        raise np.linalg.LinAlgError("rank-deficient weighted design")

    Xtw = X.T * weights
    normal = Xtw @ X
    inverse_normal = np.linalg.inv(normal)
    beta = inverse_normal @ Xtw @ y
    hat_row = target_row @ inverse_normal @ Xtw
    return beta, hat_row


def _aicc(y: np.ndarray, fitted: np.ndarray, trace_s: float) -> float:
    """Gaussian GWR AICc used by PyGWRx and matching mgwr for valid fits."""
    n = y.size
    residuals = y - fitted
    rss = max(float(residuals @ residuals), np.finfo(float).tiny)
    denominator = n - 2.0 - float(trace_s)
    if denominator <= 0.0:
        return np.inf
    return float(
        n * np.log(rss / n)
        + n * np.log(2.0 * np.pi)
        + n * (n + float(trace_s)) / denominator
    )


def _score_candidate(
    X: np.ndarray,
    y: np.ndarray,
    distances: np.ndarray,
    k: int,
    kernel: str,
    *,
    boundary_policy: str,
) -> float:
    n = y.size
    fitted = np.zeros(n, dtype=float)
    trace_s = 0.0

    for i in range(n):
        if boundary_policy == "pygwrx":
            bw = _distance_bandwidth_pygwrx(distances[i], k)
        elif boundary_policy == "mgwr":
            bw = _distance_bandwidth_mgwr(distances[i], k)
        else:
            raise ValueError("unknown boundary policy")

        weights = _kernel_from_distance_bandwidth(distances[i], bw, kernel)
        beta, hat_row = _fit_local_unpenalized(X, y, weights, X[i])
        fitted[i] = X[i] @ beta
        trace_s += float(hat_row[i])

    return _aicc(y, fitted, trace_s)


class PyGWRxAdaptiveAICcSelector:
    """Current PyGWRx exhaustive integer adaptive AICc search."""

    def __init__(self, kernel: str = "bisquare") -> None:
        if kernel not in {"bisquare", "gaussian", "exponential"}:
            raise ValueError("kernel must be bisquare, gaussian, or exponential")
        self.kernel = kernel

    def select(self, X, y, coords) -> BandwidthSearchResult:
        X = np.asarray(X, dtype=float)
        y = np.asarray(y, dtype=float).reshape(-1)
        coords = np.asarray(coords, dtype=float)
        n, p = X.shape
        lower = max(p + 1, 2, int(np.ceil(0.05 * n)))
        upper = n
        distances = cdist(coords, coords)

        trace: list[tuple[int, float]] = []
        best_k = None
        best_score = np.inf
        for k in range(lower, upper + 1):
            try:
                score = _score_candidate(
                    X, y, distances, k, self.kernel, boundary_policy="pygwrx"
                )
            except np.linalg.LinAlgError:
                score = np.inf
            trace.append((k, float(score)))
            if np.isfinite(score) and score < best_score:
                best_k, best_score = k, float(score)

        if best_k is None:
            raise RuntimeError("bandwidth selection failed for every candidate")
        result = BandwidthSearchResult(
            bandwidth=int(best_k),
            score=float(best_score),
            search_range=(int(lower), int(upper)),
            search_trace=tuple(trace),
            strategy="pygwrx_exhaustive_integer_aicc",
        )
        self.result_ = result
        return result


class MGWRCompatibleAICcSelector:
    """Reproduce mgwr 2.2.1 standard-GWR discrete golden-section AICc search."""

    def __init__(
        self,
        kernel: str = "bisquare",
        tol: float = 1.0e-6,
        max_iter: int = 200,
    ) -> None:
        if kernel not in {"bisquare", "gaussian", "exponential"}:
            raise ValueError("kernel must be bisquare, gaussian, or exponential")
        self.kernel = kernel
        self.tol = float(tol)
        self.max_iter = int(max_iter)

    def select(self, X, y, coords) -> BandwidthSearchResult:
        X = np.asarray(X, dtype=float)
        y = np.asarray(y, dtype=float).reshape(-1)
        coords = np.asarray(coords, dtype=float)
        if X.ndim != 2:
            raise ValueError("X must be a two-dimensional design matrix")
        if coords.ndim != 2 or coords.shape[1] != 2:
            raise ValueError("coords must have shape (n, 2)")
        if X.shape[0] != y.size or coords.shape[0] != y.size:
            raise ValueError("X, y, and coords must have the same rows")

        n, p = X.shape
        # mgwr Sel_BW._init_section for adaptive GWR: 40 + 2*n_vars ... n.
        a = float(40 + 2 * p)
        c = float(n)
        if a > c:
            raise ValueError("sample size is too small for mgwr-compatible search")

        distances = cdist(coords, coords)
        cache: dict[int, float] = {}
        evaluation_order: list[tuple[int, float]] = []

        def objective(value: float) -> float:
            k = int(np.round(value))
            if k not in cache:
                try:
                    score = _score_candidate(
                        X, y, distances, k, self.kernel, boundary_policy="mgwr"
                    )
                except np.linalg.LinAlgError:
                    score = np.inf
                cache[k] = float(score)
                evaluation_order.append((k, float(score)))
            return cache[k]

        # Exact constants/control flow from mgwr 2.2.1 golden_section.
        delta = 0.38197
        b = a + delta * abs(c - a)
        d = c - delta * abs(c - a)
        opt_val = None
        opt_score = np.inf
        diff = 1.0e9
        iters = 0

        while abs(diff) > self.tol and iters < self.max_iter:
            iters += 1
            b = float(np.round(b))
            d = float(np.round(d))
            score_b = objective(b)
            score_d = objective(d)

            if score_b <= score_d:
                opt_val = b
                opt_score = score_b
                c = d
                d = b
                b = a + delta * abs(c - a)
            else:
                opt_val = d
                opt_score = score_d
                a = b
                b = d
                d = c - delta * abs(c - a)

            diff = score_b - score_d

        if opt_val is None:
            raise RuntimeError("golden-section bandwidth search failed")

        best_k = int(np.round(opt_val))
        result = BandwidthSearchResult(
            bandwidth=best_k,
            score=float(opt_score),
            search_range=(int(40 + 2 * p), int(n)),
            search_trace=tuple(evaluation_order),
            strategy="mgwr_2_2_1_discrete_golden_aicc",
        )
        self.result_ = result
        self.bandwidth_ = best_k
        self.best_score_ = float(opt_score)
        self.search_range_ = result.search_range
        self.search_trace_ = result.search_trace
        return result
