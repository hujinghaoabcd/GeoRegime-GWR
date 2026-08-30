"""Compare fixed Ward K candidates under regime-restricted GWR.

This is an exploratory K-sensitivity experiment. For each existing Queen-constrained
Ward partition K=2..15, labels are held fixed and each regime fits its own
BasicGWR using only observations from that regime. No label refinement or
boundary penalty is applied.

K values containing a regime with fewer than five observations are marked
infeasible for the current 4-parameter (intercept + 3 slopes) adaptive GWR
search and are not force-fit.
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


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(DATA)
    labels = pd.read_csv(LABELS)
    if len(df) != 159 or len(labels) != 159:
        raise RuntimeError("Georgia K-sensitivity experiment expects 159 counties")

    df["AreaKey_join"] = df["AreaKey"].map(_key)
    labels["AreaKey_join"] = labels["AreaKey"].map(_key)
    kcols = [f"K{k}" for k in range(2, 16)]
    missing = [c for c in kcols if c not in labels.columns]
    if missing:
        raise RuntimeError(f"Missing Ward label columns: {missing}")
    df = df.merge(labels[["AreaKey_join", *kcols]], on="AreaKey_join", how="left", validate="1:1")

    X = df[PREDICTORS].to_numpy(dtype=float)
    y = df["PctBach"].to_numpy(dtype=float)
    coords = df[["X", "Y"]].to_numpy(dtype=float)
    Xz = (X - X.mean(axis=0)) / X.std(axis=0, ddof=0)
    yz = (y - y.mean()) / y.std(ddof=0)

    baseline = BasicGWR(bandwidth="auto", kernel="bisquare", fit_intercept=True).fit(Xz, yz, coords)
    baseline_metrics = _metrics(yz, baseline.fitted_values_, baseline.hat_matrix_)

    rows: list[dict] = []
    bw_rows: list[dict] = []

    for K in range(2, 16):
        regimes = df[f"K{K}"].to_numpy(dtype=int)
        unique, counts = np.unique(regimes, return_counts=True)
        sizes = {int(r): int(n) for r, n in zip(unique, counts)}
        min_n = min(sizes.values())
        max_n = max(sizes.values())

        if min_n < MIN_REGIME_N:
            rows.append({
                "K": K,
                "status": "infeasible_small_regime",
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
            })
            continue

        n = len(df)
        fitted = np.full(n, np.nan, dtype=float)
        hat = np.zeros((n, n), dtype=float)
        all_at_max = True

        for regime in sorted(sizes):
            idx = np.flatnonzero(regimes == regime)
            model = BasicGWR(bandwidth="auto", kernel="bisquare", fit_intercept=True).fit(
                Xz[idx], yz[idx], coords[idx]
            )
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

        if not np.all(np.isfinite(fitted)):
            raise RuntimeError(f"Non-finite fitted values for K={K}")
        m = _metrics(yz, fitted, hat)
        rows.append({
            "K": K,
            "status": "fit",
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
    best_rmse = int(fit_table.loc[fit_table["rmse"].idxmin(), "K"])
    best_mae = int(fit_table.loc[fit_table["mae"].idxmin(), "K"])
    best_r2 = int(fit_table.loc[fit_table["r2"].idxmax(), "K"])
    best_aicc = int(fit_table.loc[fit_table["aicc_fixed_partition"].idxmin(), "K"])

    mgwr_reference = None
    if MGWR_SUMMARY.exists():
        ref = json.loads(MGWR_SUMMARY.read_text(encoding="utf-8"))
        mgwr_reference = ref.get("mgwr")

    summary = {
        "status": "exploratory_fixed_partition_k_sensitivity",
        "candidate_K": list(range(2, 16)),
        "feasible_K": fit_table["K"].astype(int).tolist(),
        "infeasible_K": table.loc[table["status"] != "fit", "K"].astype(int).tolist(),
        "minimum_regime_n_required_here": MIN_REGIME_N,
        "ordinary_gwr_bandwidth": int(baseline.bandwidth_),
        "ordinary_gwr": baseline_metrics,
        "mgwr_reference": mgwr_reference,
        "best_in_sample_rmse_K": best_rmse,
        "best_in_sample_mae_K": best_mae,
        "best_in_sample_r2_K": best_r2,
        "best_conditional_aicc_K": best_aicc,
        "warning": (
            "Exploratory in-sample comparison only. Ward labels are learned from the same data; "
            "AICc treats each partition as fixed and omits partition-selection complexity."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(table.to_string(index=False))
    print("\n", json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
