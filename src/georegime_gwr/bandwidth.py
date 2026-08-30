"""Lightweight bandwidth selection for the paper-research GWR baseline.

Three search policies are intentionally explicit:

1. ``PyGWRxAdaptiveAICcSelector``
   Adaptive integer bandwidth; exhaustive AICc scan over every valid ``k``.
   This is the default research policy because it returns the global minimum
   on the discrete candidate domain.

2. ``FixedGoldenAICcSelector``
   Fixed-distance bandwidth; continuous golden-section minimization of AICc.
   This mirrors the current PyGWRx division of labour: discrete exhaustive
   search for adaptive bandwidths, continuous optimization for fixed bandwidths.

3. ``MGWRCompatibleAICcSelector``
   Compatibility mode reproducing the historical ``mgwr==2.2.1`` adaptive
   discrete golden-section search used by the canonical Georgia example.

The policies are kept separate because on canonical Georgia the exhaustive
integer AICc optimum is k=116, whereas mgwr's historical golden-section search
returns k=117.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial.distance import cdist


@dataclass(frozen=True)
class BandwidthSearchResult:
    bandwidth: int | float
    score: float
    search_range: tuple[int | float, int | float]
    search_trace: tuple[tuple[int | float, float], ...]
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
    """mgwr 2.2.1 adaptive compact-kernel convention."""
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
    ratio = distances / float(distance_bandwidth)
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
    """Gaussian GWR AICc using trace(S) as the complexity term."""
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
    bandwidth: int | float,
    kernel: str,
    *,
    policy: str,
) -> float:
    n = y.size
    fitted = np.zeros(n, dtype=float)
    trace_s = 0.0

    for i in range(n):
        if policy == "pygwrx_adaptive":
            bw = _distance_bandwidth_pygwrx(distances[i], int(bandwidth))
        elif policy == "mgwr_adaptive":
            bw = _distance_bandwidth_mgwr(distances[i], int(bandwidth))
        elif policy == "fixed":
            bw = float(bandwidth)
            if bw <= 0.0:
                return np.inf
        else:
            raise ValueError("unknown bandwidth policy")

        weights = _kernel_from_distance_bandwidth(distances[i], bw, kernel)
        beta, hat_row = _fit_local_unpenalized(X, y, weights, X[i])
        fitted[i] = X[i] @ beta
        trace_s += float(hat_row[i])

    return _aicc(y, fitted, trace_s)


def _validate_inputs(X, y, coords):
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float).reshape(-1)
    coords = np.asarray(coords, dtype=float)
    if X.ndim != 2:
        raise ValueError("X must be a two-dimensional design matrix")
    if coords.ndim != 2 or coords.shape[1] != 2:
        raise ValueError("coords must have shape (n, 2)")
    if X.shape[0] != y.size or coords.shape[0] != y.size:
        raise ValueError("X, y, and coords must have the same rows")
    return X, y, coords


class PyGWRxAdaptiveAICcSelector:
    """Exhaustive integer adaptive AICc search extracted from current PyGWRx."""

    def __init__(self, kernel: str = "bisquare") -> None:
        if kernel not in {"bisquare", "gaussian", "exponential"}:
            raise ValueError("kernel must be bisquare, gaussian, or exponential")
        self.kernel = kernel

    def select(self, X, y, coords) -> BandwidthSearchResult:
        X, y, coords = _validate_inputs(X, y, coords)
        n, p = X.shape
        lower = max(p + 1, 2, int(np.ceil(0.05 * n)))
        upper = n
        distances = cdist(coords, coords)

        trace: list[tuple[int, float]] = []
        best_k: int | None = None
        best_score = np.inf
        for k in range(lower, upper + 1):
            try:
                score = _score_candidate(
                    X, y, distances, k, self.kernel, policy="pygwrx_adaptive"
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


class FixedGoldenAICcSelector:
    """Continuous fixed-distance AICc search using golden-section reduction."""

    RESPHI = (np.sqrt(5.0) - 1.0) / 2.0

    def __init__(
        self,
        kernel: str = "bisquare",
        tol: float = 1.0e-4,
        max_iter: int = 100,
    ) -> None:
        if kernel not in {"bisquare", "gaussian", "exponential"}:
            raise ValueError("kernel must be bisquare, gaussian, or exponential")
        self.kernel = kernel
        self.tol = float(tol)
        self.max_iter = int(max_iter)

    def select(self, X, y, coords) -> BandwidthSearchResult:
        X, y, coords = _validate_inputs(X, y, coords)
        distances = cdist(coords, coords)
        positive = distances[distances > 0.0]
        if positive.size == 0:
            raise ValueError("fixed bandwidth search requires distinct coordinates")

        lower = 0.5 * float(np.min(positive))
        upper = 2.0 * float(np.max(distances))
        cache: dict[float, float] = {}

        def objective(value: float) -> float:
            key = float(value)
            if key not in cache:
                try:
                    score = _score_candidate(
                        X, y, distances, key, self.kernel, policy="fixed"
                    )
                except np.linalg.LinAlgError:
                    score = np.inf
                cache[key] = float(score)
            return cache[key]

        a = lower
        b = upper
        # Evaluate endpoints so boundary optima are visible.
        objective(a)
        objective(b)
        c = b - self.RESPHI * (b - a)
        d = a + self.RESPHI * (b - a)
        fc = objective(c)
        fd = objective(d)

        iterations = 0
        while iterations < self.max_iter:
            midpoint = 0.5 * (a + b)
            if (b - a) <= self.tol * (1.0 + abs(midpoint)):
                break
            if fc <= fd:
                b, d, fd = d, c, fc
                c = b - self.RESPHI * (b - a)
                fc = objective(c)
            else:
                a, c, fc = c, d, fd
                d = a + self.RESPHI * (b - a)
                fd = objective(d)
            iterations += 1

        finite = [(bw, score) for bw, score in cache.items() if np.isfinite(score)]
        if not finite:
            raise RuntimeError("fixed golden-section bandwidth search failed")
        best_bw, best_score = min(finite, key=lambda item: (item[1], item[0]))
        result = BandwidthSearchResult(
            bandwidth=float(best_bw),
            score=float(best_score),
            search_range=(float(lower), float(upper)),
            search_trace=tuple(sorted(cache.items(), key=lambda item: item[0])),
            strategy="pygwrx_fixed_golden_aicc",
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
        bw_min: int | None = None,
        bw_max: int | None = None,
    ) -> None:
        if kernel not in {"bisquare", "gaussian", "exponential"}:
            raise ValueError("kernel must be bisquare, gaussian, or exponential")
        self.kernel = kernel
        self.tol = float(tol)
        self.max_iter = int(max_iter)
        self.bw_min = bw_min
        self.bw_max = bw_max

    def select(self, X, y, coords) -> BandwidthSearchResult:
        X, y, coords = _validate_inputs(X, y, coords)
        n, p = X.shape
        # mgwr Sel_BW._init_section for adaptive GWR.
        a = float(40 + 2 * p if self.bw_min is None else self.bw_min)
        c = float(n if self.bw_max is None else self.bw_max)
        if a > c:
            raise ValueError("sample size/search bounds are invalid for mgwr-compatible search")

        distances = cdist(coords, coords)
        cache: dict[int, float] = {}
        evaluation_order: list[tuple[int, float]] = []

        def objective(value: float) -> float:
            k = int(np.round(value))
            if k not in cache:
                try:
                    score = _score_candidate(
                        X, y, distances, k, self.kernel, policy="mgwr_adaptive"
                    )
                except np.linalg.LinAlgError:
                    score = np.inf
                cache[k] = float(score)
                evaluation_order.append((k, float(score)))
            return cache[k]

        delta = 0.38197
        b = a + delta * abs(c - a)
        d = c - delta * abs(c - a)
        opt_val: float | None = None
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
                opt_val, opt_score = b, score_b
                c, d = d, b
                b = a + delta * abs(c - a)
            else:
                opt_val, opt_score = d, score_d
                a, b = b, d
                d = c - delta * abs(c - a)
            diff = score_b - score_d

        if opt_val is None:
            raise RuntimeError("mgwr-compatible golden-section bandwidth search failed")

        best_k = int(np.round(opt_val))
        result = BandwidthSearchResult(
            bandwidth=best_k,
            score=float(opt_score),
            search_range=(int(40 + 2 * p if self.bw_min is None else self.bw_min), int(n if self.bw_max is None else self.bw_max)),
            search_trace=tuple(evaluation_order),
            strategy="mgwr_2_2_1_discrete_golden_aicc",
        )
        self.result_ = result
        return result
