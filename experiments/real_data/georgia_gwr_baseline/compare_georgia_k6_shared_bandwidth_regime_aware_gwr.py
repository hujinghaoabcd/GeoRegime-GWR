"""Exploratory unified K=6 regime-aware GWR with one shared adaptive bandwidth.

The 159 Georgia counties remain in one model. For focal county i and observation j,

    w_ij^GR = K(d_ij / b_i(k)) * I(z_i == z_j)

where k is one shared adaptive neighbour-count bandwidth for the whole model and
b_i(k) is the ordinary global k-th-neighbour distance at focal location i.
Cross-regime weights are then forced to zero. This is deliberately an experiment,
not a frozen GR-GWR definition.
"""

from __future__ import annotations

from pathlib import Path
import json
import sys

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from georegime_gwr.gwr import BasicGWR

DATA = ROOT / "data" / "raw" / "georgia" / "GData_utm.csv"
LABELS = ROOT / "results" / "spatial" / "georgia_initial_ward_regimes" / "initial_regime_labels_k2_k15.csv"
K6_CURRENT = ROOT / "results" / "real_data" / "georgia_k6_regime_restricted_gwr" / "summary.json"
MGWR_SUMMARY = ROOT / "results" / "real_data" / "georgia_mgwr_benchmark" / "summary.json"
OUT = ROOT / "results" / "real_data" / "georgia_k6_shared_bandwidth_regime_aware_gwr"
PREDICTORS = ["PctFB", "PctBlack", "PctRural"]


def _key(value) -> str:
    try:
        return str(int(float(value)))
    except (TypeError, ValueError):
        return str(value)


def _aicc(y: np.ndarray, fitted: np.ndarray, trace_s: float) -> float:
    n = y.size
    residuals = y - fitted
    rss = max(float(residuals @ residuals), np.finfo(float).tiny)
    denominator = n - 2.0 - float(trace_s)
    if denominator <= 0.0:
        return float("inf")
    return float(
        n * np.log(rss / n)
        + n * np.log(2.0 * np.pi)
        + n * (n + float(trace_s)) / denominator
    )


def _metrics(y: np.ndarray, fitted: np.ndarray, hat: np.ndarray) -> dict[str, float]:
    residuals = y - fitted
    rss = float(residuals @ residuals)
    rmse = float(np.sqrt(np.mean(residuals**2)))
    mae = float(np.mean(np.abs(residuals)))
    tss = float(np.sum((y - np.mean(y)) ** 2))
    trace_s = float(np.trace(hat))
    return {
        "rss": rss,
        "rmse": rmse,
        "mae": mae,
        "r2": float(1.0 - rss / tss),
        "trace_hat": trace_s,
        "aicc_fixed_partition": _aicc(y, fitted, trace_s),
    }


def _kernel_weights(distances: np.ndarray, k: int) -> np.ndarray:
    """Match BasicGWR adaptive bisquare + PyGWRx boundary semantics."""
    k_eff = min(int(k), distances.size)
    bw = float(np.partition(distances, k_eff - 1)[k_eff - 1])
    if bw <= 1e-12:
        positive = distances[distances > 1e-12]
        bw = float(np.min(positive)) if positive.size else 1.0
    bw = float(np.nextafter(bw, np.inf))
    ratio = distances / bw
    return np.where(ratio < 1.0, (1.0 - ratio**2) ** 2, 0.0)


def _fit_shared_masked(
    Xd: np.ndarray,
    y: np.ndarray,
    distances: np.ndarray,
    labels: np.ndarray,
    k: int,
):
    """Fit one regime-aware smoother with one shared adaptive k.

    Candidates are rejected if any local weighted design has fewer than p
    positive same-regime observations or is rank deficient. We do not use ridge
    merely to make an otherwise unidentified bandwidth candidate selectable.
    """
    n, p = Xd.shape
    fitted = np.empty(n, dtype=float)
    params = np.empty((n, p), dtype=float)
    hat = np.zeros((n, n), dtype=float)
    positive_counts = np.empty(n, dtype=int)

    for i in range(n):
        w = _kernel_weights(distances[i], k)
        w = w * (labels == labels[i])
        positive = w > 0.0
        positive_counts[i] = int(np.sum(positive))
        if positive_counts[i] < p:
            return None
        if np.linalg.matrix_rank(Xd[positive]) < p:
            return None
        beta, C = BasicGWR._solve_local(Xd, y, w)
        if not np.all(np.isfinite(beta)):
            return None
        params[i] = beta
        fitted[i] = Xd[i] @ beta
        hat[i] = Xd[i] @ C

    return params, fitted, hat, positive_counts


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(DATA)
    labels_df = pd.read_csv(LABELS)
    df["AreaKey_join"] = df["AreaKey"].map(_key)
    labels_df["AreaKey_join"] = labels_df["AreaKey"].map(_key)
    df = df.merge(labels_df[["AreaKey_join", "K6"]], on="AreaKey_join", how="left", validate="1:1")
    if len(df) != 159 or df["K6"].isna().any():
        raise RuntimeError("Expected complete K=6 labels for all 159 Georgia counties")

    X = df[PREDICTORS].to_numpy(dtype=float)
    y = df["PctBach"].to_numpy(dtype=float)
    coords = df[["X", "Y"]].to_numpy(dtype=float)
    labels = df["K6"].to_numpy(dtype=int)

    Xz = (X - X.mean(axis=0)) / X.std(axis=0, ddof=0)
    yz = (y - y.mean()) / y.std(ddof=0)
    Xd = np.column_stack([np.ones(len(df)), Xz])
    distances = cdist(coords, coords)

    # Ordinary GWR reference under the frozen research-default search.
    ordinary = BasicGWR(bandwidth="auto", kernel="bisquare", fit_intercept=True).fit(Xz, yz, coords)
    ordinary_metrics = _metrics(yz, ordinary.fitted_values_, ordinary.hat_matrix_)

    search_rows = []
    best = None
    for k in range(8, len(df) + 1):
        fit = _fit_shared_masked(Xd, yz, distances, labels, k)
        if fit is None:
            search_rows.append({"k": k, "status": "unidentified"})
            continue
        params, fitted, hat, positive_counts = fit
        m = _metrics(yz, fitted, hat)
        row = {
            "k": k,
            "status": "fit",
            **m,
            "min_positive_same_regime": int(np.min(positive_counts)),
            "median_positive_same_regime": float(np.median(positive_counts)),
            "max_positive_same_regime": int(np.max(positive_counts)),
        }
        search_rows.append(row)
        if np.isfinite(m["aicc_fixed_partition"]) and (
            best is None or m["aicc_fixed_partition"] < best[0]
        ):
            best = (m["aicc_fixed_partition"], k, params, fitted, hat, positive_counts, m)

    if best is None:
        raise RuntimeError("No identifiable shared-bandwidth regime-aware candidate")

    _, best_k, params, fitted, hat, positive_counts, shared_metrics = best
    pd.DataFrame(search_rows).to_csv(OUT / "shared_bandwidth_search.csv", index=False)

    county = pd.DataFrame({
        "AreaKey": df["AreaKey"],
        "regime": labels,
        "fitted_shared_regime_aware": fitted,
        "residual_shared_regime_aware": yz - fitted,
        "positive_same_regime_count": positive_counts,
        "Intercept": params[:, 0],
        "PctFB": params[:, 1],
        "PctBlack": params[:, 2],
        "PctRural": params[:, 3],
    })
    county.to_csv(OUT / "county_results.csv", index=False)

    current = json.loads(K6_CURRENT.read_text(encoding="utf-8"))
    independent_metrics = current["regime_restricted_metrics"]
    mgwr = None
    if MGWR_SUMMARY.exists():
        mgwr = json.loads(MGWR_SUMMARY.read_text(encoding="utf-8")).get("mgwr")

    comparison = pd.DataFrame([
        {"model": "ordinary_gwr", "bandwidth_spec": f"adaptive k={ordinary.bandwidth_}", **ordinary_metrics},
        {"model": "k6_shared_bandwidth_regime_aware_gwr", "bandwidth_spec": f"one shared adaptive k={best_k} + regime mask", **shared_metrics},
        {"model": "k6_independent_regime_gwr", "bandwidth_spec": "six independently selected within-regime bandwidths", **{k: independent_metrics[k] for k in ["rss", "rmse", "mae", "r2", "trace_hat"]}, "aicc_fixed_partition": independent_metrics["aicc"]},
    ])
    if mgwr is not None:
        comparison = pd.concat([
            comparison,
            pd.DataFrame([{"model": "mgwr", "bandwidth_spec": "variable-specific adaptive bandwidths", **{k: mgwr[k] for k in ["rss", "rmse", "mae", "r2", "trace_hat"]}, "aicc_fixed_partition": mgwr["aicc"]}]),
        ], ignore_index=True)
    comparison.to_csv(OUT / "model_comparison.csv", index=False)

    summary = {
        "status": "exploratory_unified_shared_bandwidth_regime_aware_gwr",
        "n_counties": int(len(df)),
        "K": 6,
        "labels_fixed": True,
        "one_model_fit": True,
        "weight_rule": "bisquare(global adaptive k-th-neighbour distance) * I(same K6 regime)",
        "shared_adaptive_bandwidth_k": int(best_k),
        "bandwidth_search_range": [8, int(len(df))],
        "bandwidth_selected_by": "aggregate conditional/fixed-partition Gaussian AICc",
        "same_regime_positive_count_at_best_k": {
            "min": int(np.min(positive_counts)),
            "median": float(np.median(positive_counts)),
            "max": int(np.max(positive_counts)),
        },
        "ordinary_gwr": ordinary_metrics,
        "shared_bandwidth_regime_aware_gwr": shared_metrics,
        "independent_regime_gwr": independent_metrics,
        "mgwr_reference": mgwr,
        "warning": (
            "Exploratory fixed-label comparison only. K6 labels were learned from the same response-derived pilot GWR; "
            "the reported AICc does not charge partition-selection complexity."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(comparison.to_string(index=False))
    print("\n", json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
