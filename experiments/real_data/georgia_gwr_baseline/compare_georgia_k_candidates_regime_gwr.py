"""Compare fixed Ward K candidates under regime-restricted GWR.

Exploratory K-sensitivity test. Each Ward partition K=2..15 is held fixed and
then each regime fits its own BasicGWR. No label refinement or boundary penalty
is applied. K values that contain too-small or non-identifiable regimes are
recorded as infeasible rather than force-fit with ridge regularization.
"""

from __future__ import annotations

from pathlib import Path
import json
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from georegime_gwr.gwr import BasicGWR

DATA = ROOT / "data" / "raw" / "georgia" / "GData_utm.csv"
LABELS = ROOT / "results" / "spatial" / "georgia_initial_ward_regimes" / "initial_regime_labels_k2_k15.csv"
MGWR_SUMMARY = ROOT / "results" / "real_data" / "georgia_mgwr_benchmark" / "summary.json"
OUT = ROOT / "results" / "real_data" / "georgia_k_sensitivity_regime_gwr"
PREDICTORS = ["PctFB", "PctBlack", "PctRural"]
MIN_REGIME_N = 5


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
    return float(n * np.log(rss / n) + n * np.log(2.0 * np.pi) + n * (n + float(trace_s)) / denominator)


def _metrics(y: np.ndarray, fitted: np.ndarray, hat: np.ndarray) -> dict[str, float]:
    residuals = y - fitted
    rss = float(residuals @ residuals)
    rmse = float(np.sqrt(np.mean(residuals**2)))
    mae = float(np.mean(np.abs(residuals)))
    tss = float(np.sum((y - np.mean(y)) ** 2))
    r2 = float(1.0 - rss / tss)
    trace_s = float(np.trace(hat))
    return {
        "rss": rss,
        "rmse": rmse,
        "mae": mae,
        "r2": r2,
        "trace_hat": trace_s,
        "aicc_fixed_partition": _aicc(y, fitted, trace_s),
    }


def _blank_row(K, status, min_n, max_n, sizes, failed_regime=None, failure_reason=None):
    return {
        "K": K,
        "status": status,
        "failed_regime": failed_regime,
        "failure_reason": failure_reason,
        "min_regime_n": min_n,
        "max_regime_n": max_n,
        "regime_sizes": ";".join(str(sizes[r]) for r in sorted(sizes)),
        "rss": np.nan,
        "rmse": np.nan,
        "mae": np.nan,
        "r2": np.nan,
        "trace_hat": np.nan,
        "aicc_fixed_partition": np.nan,
        "all_bandwidths_at_regime_max": np.nan,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(DATA)
    labels = pd.read_csv(LABELS)
    if len(df) != 159 or len(labels) != 159:
        raise RuntimeError("Georgia K-sensitivity experiment expects 159 counties")

    df["AreaKey_join"] = df["AreaKey"].map(_key)
    labels["AreaKey_join"] = labels["AreaKey"].map(_key)
    kcols = [f"K{k}" for k in range(2, 16)]
    df = df.merge(labels[["AreaKey_join", *kcols]], on="AreaKey_join", how="left", validate="1:1")

    X = df[PREDICTORS].to_numpy(dtype=float)
    y = df["PctBach"].to_numpy(dtype=float)
    coords = df[["X", "Y"]].to_numpy(dtype=float)
    Xz = (X - X.mean(axis=0)) / X.std(axis=0, ddof=0)
    yz = (y - y.mean()) / y.std(ddof=0)

    baseline = BasicGWR(bandwidth="auto", kernel="bisquare", fit_intercept=True).fit(Xz, yz, coords)
    baseline_metrics = _metrics(yz, baseline.fitted_values_, baseline.hat_matrix_)

    rows, bw_rows = [], []
    for K in range(2, 16):
        regimes = df[f"K{K}"].to_numpy(dtype=int)
        unique, counts = np.unique(regimes, return_counts=True)
        sizes = {int(r): int(n) for r, n in zip(unique, counts)}
        min_n, max_n = min(sizes.values()), max(sizes.values())

        if min_n < MIN_REGIME_N:
            rows.append(_blank_row(K, "infeasible_small_regime", min_n, max_n, sizes, failure_reason=f"min regime n={min_n} < {MIN_REGIME_N}"))
            continue

        n = len(df)
        fitted = np.full(n, np.nan, dtype=float)
        hat = np.zeros((n, n), dtype=float)
        all_at_max = True
        failed = None

        for regime in sorted(sizes):
            idx = np.flatnonzero(regimes == regime)
            try:
                model = BasicGWR(bandwidth="auto", kernel="bisquare", fit_intercept=True).fit(Xz[idx], yz[idx], coords[idx])
            except (RuntimeError, np.linalg.LinAlgError, ValueError) as exc:
                failed = (int(regime), f"{type(exc).__name__}: {exc}")
                break

            fitted[idx] = model.fitted_values_
            hat[np.ix_(idx, idx)] = model.hat_matrix_
            bw = int(model.bandwidth_)
            all_at_max = all_at_max and (bw == idx.size)
            bw_rows.append({
                "K": K,
                "regime": int(regime),
                "n": int(idx.size),
                "selected_bandwidth": bw,
                "bandwidth_fraction_of_regime": float(bw / idx.size),
                "bandwidth_at_regime_max": bool(bw == idx.size),
            })

        if failed is not None:
            rows.append(_blank_row(K, "infeasible_bandwidth_search", min_n, max_n, sizes, failed_regime=failed[0], failure_reason=failed[1]))
            continue

        m = _metrics(yz, fitted, hat)
        rows.append({
            "K": K,
            "status": "fit",
            "failed_regime": np.nan,
            "failure_reason": None,
            "min_regime_n": min_n,
            "max_regime_n": max_n,
            "regime_sizes": ";".join(str(sizes[r]) for r in sorted(sizes)),
            **m,
            "all_bandwidths_at_regime_max": bool(all_at_max),
        })

    table = pd.DataFrame(rows)
    table.to_csv(OUT / "k2_k15_regime_gwr_comparison.csv", index=False)
    pd.DataFrame(bw_rows).to_csv(OUT / "per_regime_bandwidths.csv", index=False)

    fit_table = table.loc[table["status"] == "fit"].copy()
    summary = {
        "status": "exploratory_fixed_partition_k_sensitivity",
        "candidate_K": list(range(2, 16)),
        "feasible_K": fit_table["K"].astype(int).tolist(),
        "infeasible_K": table.loc[table["status"] != "fit", "K"].astype(int).tolist(),
        "minimum_regime_n_required_here": MIN_REGIME_N,
        "ordinary_gwr_bandwidth": int(baseline.bandwidth_),
        "ordinary_gwr": baseline_metrics,
        "mgwr_reference": json.loads(MGWR_SUMMARY.read_text(encoding="utf-8")).get("mgwr") if MGWR_SUMMARY.exists() else None,
        "best_in_sample_rmse_K": int(fit_table.loc[fit_table["rmse"].idxmin(), "K"]) if len(fit_table) else None,
        "best_in_sample_mae_K": int(fit_table.loc[fit_table["mae"].idxmin(), "K"]) if len(fit_table) else None,
        "best_in_sample_r2_K": int(fit_table.loc[fit_table["r2"].idxmax(), "K"]) if len(fit_table) else None,
        "best_conditional_aicc_K": int(fit_table.loc[fit_table["aicc_fixed_partition"].idxmin(), "K"]) if len(fit_table) else None,
        "warning": "Exploratory in-sample comparison only. Ward labels are learned from the same data; conditional AICc omits partition-selection complexity.",
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(table.to_string(index=False))
    print("\n" + json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
