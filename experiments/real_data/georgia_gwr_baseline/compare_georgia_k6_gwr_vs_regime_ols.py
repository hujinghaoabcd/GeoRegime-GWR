"""Compare ordinary GWR, fixed-K6 regime-restricted GWR, and fixed-K6 regime OLS.

This experiment isolates a specific methodological question: after the current
K=6 segmentation is fixed, is the gain mainly due to blocking cross-boundary
borrowing, or do we still need within-regime local GWR variation?

The K=6 labels are treated as fixed for this exploratory comparison.  The AICc
reported here therefore does NOT account for the data-driven cost of discovering
the partition itself and must not be used as a final model-selection criterion.
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
OUT = ROOT / "results" / "real_data" / "georgia_k6_model_form_comparison"

PREDICTORS = ["PctFB", "PctBlack", "PctRural"]
TERMS = ["Intercept", *PREDICTORS]


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
    r2 = float(1.0 - rss / tss) if tss > 0 else float("nan")
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
    labels = pd.read_csv(LABELS)[["AreaKey", "K6"]].copy()
    df["AreaKey_join"] = df["AreaKey"].map(_key)
    labels["AreaKey_join"] = labels["AreaKey"].map(_key)
    df = df.merge(labels[["AreaKey_join", "K6"]], on="AreaKey_join", how="left", validate="1:1")
    if len(df) != 159 or df["K6"].isna().any():
        raise RuntimeError("Expected complete K6 labels for 159 Georgia counties")

    regimes = df["K6"].to_numpy(dtype=int)
    if sorted(np.unique(regimes).tolist()) != [1, 2, 3, 4, 5, 6]:
        raise RuntimeError(f"Unexpected K6 labels: {np.unique(regimes).tolist()}")

    X = df[PREDICTORS].to_numpy(dtype=float)
    y = df["PctBach"].to_numpy(dtype=float)
    coords = df[["X", "Y"]].to_numpy(dtype=float)
    Xz = (X - X.mean(axis=0)) / X.std(axis=0, ddof=0)
    yz = (y - y.mean()) / y.std(ddof=0)
    Xd = np.column_stack([np.ones(len(df)), Xz])

    # 1) Validated ordinary GWR baseline.
    ordinary = BasicGWR(bandwidth="auto", kernel="bisquare", fit_intercept=True).fit(Xz, yz, coords)
    if int(ordinary.bandwidth_) != 116:
        raise RuntimeError(f"Expected ordinary GWR bandwidth 116; got {ordinary.bandwidth_}")
    ordinary_metrics = _metrics(yz, ordinary.fitted_values_, ordinary.hat_matrix_)

    # 2) Fixed-K6 regime-restricted GWR: each regime independently selects its
    # adaptive bandwidth; observations from other regimes are never used.
    n = len(df)
    p = Xd.shape[1]
    rgwr_fitted = np.full(n, np.nan)
    rgwr_params = np.full((n, p), np.nan)
    rgwr_hat = np.zeros((n, n))
    rgwr_bandwidths: dict[str, int] = {}

    # 3) Fixed-K6 regime OLS: exactly one coefficient vector per regime.
    rols_fitted = np.full(n, np.nan)
    rols_params = np.full((n, p), np.nan)
    rols_hat = np.zeros((n, n))
    regime_rows: list[dict[str, float | int]] = []

    for regime in range(1, 7):
        idx = np.flatnonzero(regimes == regime)

        rgwr = BasicGWR(bandwidth="auto", kernel="bisquare", fit_intercept=True).fit(
            Xz[idx], yz[idx], coords[idx]
        )
        rgwr_bandwidths[str(regime)] = int(rgwr.bandwidth_)
        rgwr_fitted[idx] = rgwr.fitted_values_
        rgwr_params[idx] = rgwr.parameters_
        rgwr_hat[np.ix_(idx, idx)] = rgwr.hat_matrix_

        Xr = Xd[idx]
        yr = yz[idx]
        beta, _, rank, _ = np.linalg.lstsq(Xr, yr, rcond=None)
        if int(rank) != p:
            raise RuntimeError(f"Regime {regime} OLS design rank {rank}, expected {p}")
        pinv = np.linalg.pinv(Xr)
        Hr = Xr @ pinv
        fitted_r = Xr @ beta
        rols_fitted[idx] = fitted_r
        rols_params[idx] = beta
        rols_hat[np.ix_(idx, idx)] = Hr

        rgwr_m = _metrics(yr, rgwr.fitted_values_, rgwr.hat_matrix_)
        rols_m = _metrics(yr, fitted_r, Hr)
        regime_rows.append(
            {
                "regime": regime,
                "n": int(idx.size),
                "rgwr_bandwidth": int(rgwr.bandwidth_),
                "rgwr_trace_hat": rgwr_m["trace_hat"],
                "rgwr_rss": rgwr_m["rss"],
                "rols_trace_hat": rols_m["trace_hat"],
                "rols_rss": rols_m["rss"],
                "rols_rmse": rols_m["rmse"],
                "rols_r2": rols_m["r2"],
            }
        )

    rgwr_metrics = _metrics(yz, rgwr_fitted, rgwr_hat)
    rols_metrics = _metrics(yz, rols_fitted, rols_hat)

    comparison = pd.DataFrame(
        [
            {"model": "ordinary_gwr", "structure": "one global GWR surface", **ordinary_metrics},
            {"model": "k6_regime_restricted_gwr", "structure": "six independent within-regime GWR surfaces", **rgwr_metrics},
            {"model": "k6_regime_ols", "structure": "six regimes, one coefficient vector per regime", **rols_metrics},
        ]
    )
    comparison.to_csv(OUT / "model_comparison.csv", index=False)
    pd.DataFrame(regime_rows).to_csv(OUT / "per_regime_comparison.csv", index=False)

    county = pd.DataFrame({"AreaKey": df["AreaKey"], "regime": regimes})
    for j, term in enumerate(TERMS):
        county[f"ordinary_{term}"] = ordinary.parameters_[:, j]
        county[f"regime_gwr_{term}"] = rgwr_params[:, j]
        county[f"regime_ols_{term}"] = rols_params[:, j]
    county["y_z"] = yz
    county["ordinary_fitted_z"] = ordinary.fitted_values_
    county["regime_gwr_fitted_z"] = rgwr_fitted
    county["regime_ols_fitted_z"] = rols_fitted
    county.to_csv(OUT / "county_model_results.csv", index=False)

    summary = {
        "status": "exploratory_fixed_k6_model_form_comparison",
        "K": 6,
        "labels_fixed": True,
        "partition_selection_complexity_accounted_for": False,
        "standardization": "global X/y z-score, ddof=0",
        "ordinary_gwr_bandwidth": int(ordinary.bandwidth_),
        "regime_gwr_bandwidths": rgwr_bandwidths,
        "ordinary_gwr": ordinary_metrics,
        "k6_regime_restricted_gwr": rgwr_metrics,
        "k6_regime_ols": rols_metrics,
        "interpretation_guard": (
            "AICc values treat the K6 partition as fixed and therefore omit the complexity of discovering labels from y. "
            "Use only as an exploratory conditional comparison."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(comparison.to_string(index=False))
    print("\nPer-regime comparison:\n", pd.DataFrame(regime_rows).to_string(index=False))
    print("\n", json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
