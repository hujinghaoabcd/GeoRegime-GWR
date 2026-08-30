"""Lightweight bandwidth selection extracted from PyGWRx for paper experiments.

Only the standard Gaussian GWR path needed by this research repository is kept
here.  The adaptive AICc selector follows the current PyGWRx policy:

- adaptive bandwidth is an integer neighbour order k;
- automatic range is max(p + 1, 2, ceil(0.05*n)) ... n;
- every valid integer k is evaluated (no approximate optimizer);
- each candidate fits a full unpenalized Gaussian GWR;
- the score is the Gaussian GWR AICc based on RSS and trace(S).

This keeps the research repository small while preserving the exact bandwidth
selection semantics relevant to the Georgia standard-GWR benchmark.
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


def _adaptive_distance_bandwidth(distances: np.ndarray, k: int) -> float:
    """Return the adaptive distance scale using the current PyGWRx rule."""
    if k < 1 or k > distances.size:
        raise ValueError(f"k must satisfy 1 <= k <= {distances.size}; got {k}")

    bandwidth = float(np.partition(distances, k - 1)[k - 1])
    if bandwidth <= 0.0:
        positive = distances[distances > 0.0]
        if positive.size == 0:
            raise np.linalg.LinAlgError("all distances are zero")
        bandwidth = float(np.min(positive))

    # PyGWRx uses the next representable float so a compact kernel includes
    # the k-th boundary neighbour.
    return float(np.nextafter(bandwidth, np.inf))


def _kernel_weights(distances: np.ndarray, k: int, kernel: str) -> np.ndarray:
    bw = _adaptive_distance_bandwidth(distances, k)
    ratio = distances / bw

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
    """Unpenalized local WLS used for bandwidth scoring.

    Invalid/rank-deficient candidates are rejected rather than ridge-regularized,
    matching the current PyGWRx bandwidth selector policy.
    """
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
    """Gaussian GWR AICc formula used by current PyGWRx."""
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


class AdaptiveAICcSelector:
    """Exhaustive adaptive-neighbour AICc selector extracted from PyGWRx."""

    def __init__(self, kernel: str = "bisquare") -> None:
        if kernel not in {"bisquare", "gaussian", "exponential"}:
            raise ValueError("kernel must be bisquare, gaussian, or exponential")
        self.kernel = kernel

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
        lower = max(p + 1, 2, int(np.ceil(0.05 * n)))
        upper = n
        if lower > upper:
            raise ValueError("sample size is too small for adaptive bandwidth selection")

        distances = cdist(coords, coords)
        trace: list[tuple[int, float]] = []
        best_k: int | None = None
        best_score = np.inf

        for k in range(lower, upper + 1):
            fitted = np.zeros(n, dtype=float)
            trace_s = 0.0
            valid = True

            for i in range(n):
                try:
                    weights = _kernel_weights(distances[i], k, self.kernel)
                    beta, hat_row = _fit_local_unpenalized(X, y, weights, X[i])
                except np.linalg.LinAlgError:
                    valid = False
                    break
                fitted[i] = X[i] @ beta
                trace_s += float(hat_row[i])

            score = _aicc(y, fitted, trace_s) if valid else np.inf
            trace.append((k, float(score)))
            if np.isfinite(score) and score < best_score:
                best_k = k
                best_score = float(score)

        if best_k is None:
            raise RuntimeError("bandwidth selection failed for every candidate")

        result = BandwidthSearchResult(
            bandwidth=best_k,
            score=best_score,
            search_range=(lower, upper),
            search_trace=tuple(trace),
        )
        self.result_ = result
        self.bandwidth_ = best_k
        self.best_score_ = best_score
        self.search_range_ = (lower, upper)
        self.search_trace_ = tuple(trace)
        return result
