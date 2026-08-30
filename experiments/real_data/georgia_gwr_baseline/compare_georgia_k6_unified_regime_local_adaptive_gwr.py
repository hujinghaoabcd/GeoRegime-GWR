"""Test one unified K=6 regime-aware GWR using regime-local adaptive neighborhoods.

All 159 Georgia counties are fitted in one pass. For focal county i, the
adaptive distance scale is defined only from observations in i's own regime:

    w_ij = K(d_ij / b_i) * I(z_i == z_j)

For this first equivalence experiment, k_i equals the size of i's regime.
Because the current six independent regime GWR fits all selected their regime
size as the optimal adaptive bandwidth, this unified formulation should
numerically reproduce them if the weight geometry is truly equivalent.
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
CURRENT_SUMMARY = ROOT / "results" / "real_data" / "georgia_k6_regime_restricted_gwr" / "summary.json"
CURRENT_COUNTY = ROOT / "results" / "real_data" / "georgia_k6_regime_restricted_gwr" / "county_results.csv"
OUT = ROOT / "results" / "real_data" / "georgia_k6_unified_regime_local_adaptive_gwr"
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
    trace_s = float(np.trace(hat))
    return {
        "rss": rss,
        "rmse": rmse,
        "mae": mae,
        "r2": float(1.0 - rss / tss),
        "trace_hat": trace_s,
        "aicc_fixed_partition": _aicc(y, fitted, trace_s),
    }


def _regime_local_weights(
    distances_i: np.ndarray,
    same_regime: np.ndarray,
    k_within_regime: int,
) -> np.ndarray:
    """Adaptive bisquare scale determined only by same-regime distances."""
    d_same = distances_i[same_regime]
    k = min(int(k_within_regime), d_same.size)
    if k < 1:
        raise ValueError("within-regime adaptive bandwidth must be >= 1")
    bw = float(np.partition(d_same, k - 1)[k - 1])
    if bw <= 1e-12:
        positive = d_same[d_same > 1e-12]
        bw = float(np.min(positive)) if positive.size else 1.0
    bw = float(np.nextafter(bw, np.inf))

    weights = np.zeros_like(distances_i, dtype=float)
    ratio = d_same / bw
    weights[same_regime] = np.where(ratio < 1.0, (1.0 - ratio**2) ** 2, 0.0)
    return weights


def _fit_unified(
    Xd: np.ndarray,
    y: np.ndarray,
    distances: np.ndarray,
    labels: np.ndarray,
):
    n, p = Xd.shape
    params = np.empty((n, p), dtype=float)
    fitted = np.empty(n, dtype=float)
    hat = np.zeros((n, n), dtype=float)
    effective_k = np.empty(n, dtype=int)

    regime_sizes = {int(r): int(np.sum(labels == r)) for r in np.unique(labels)}

    for i in range(n):
        same = labels == labels[i]
        k_i = regime_sizes[int(labels[i])]
        w = _regime_local_weights(distances[i], same, k_i)
        beta, C = BasicGWR._solve_local(Xd, y, w)
        params[i] = beta
        fitted[i] = Xd[i] @ beta
        hat[i] = Xd[i] @ C
        effective_k[i] = k_i

    return params, fitted, hat, effective_k, regime_sizes


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(DATA)
    labels_df = pd.read_csv(LABELS)
    current_county = pd.read_csv(CURRENT_COUNTY)
    current_summary = json.loads(CURRENT_SUMMARY.read_text(encoding="utf-8"))

    df["AreaKey_join"] = df["AreaKey"].map(_key)
    labels_df["AreaKey_join"] = labels_df["AreaKey"].map(_key)
    current_county["AreaKey_join"] = current_county["AreaKey"].map(_key)
    df = df.merge(labels_df[["AreaKey_join", "K6"]], on="AreaKey_join", how="left", validate="1:1")
    df = df.merge(
        current_county[[
            "AreaKey_join",
            "restricted_Intercept",
            "restricted_PctFB",
            "restricted_PctBlack",
            "restricted_PctRural",
            "restricted_fitted_z",
        ]],
        on="AreaKey_join",
        how="left",
        validate="1:1",
    )
    if len(df) != 159 or df["K6"].isna().any():
        raise RuntimeError("Expected complete K=6 Georgia inputs")

    X = df[PREDICTORS].to_numpy(dtype=float)
    y = df["PctBach"].to_numpy(dtype=float)
    coords = df[["X", "Y"]].to_numpy(dtype=float)
    labels = df["K6"].to_numpy(dtype=int)

    Xz = (X - X.mean(axis=0)) / X.std(axis=0, ddof=0)
    yz = (y - y.mean()) / y.std(ddof=0)
    Xd = np.column_stack([np.ones(len(df)), Xz])
    distances = cdist(coords, coords)

    params, fitted, hat, effective_k, regime_sizes = _fit_unified(Xd, yz, distances, labels)
    metrics = _metrics(yz, fitted, hat)

    reference_params = df[[f"restricted_{t}" for t in TERMS]].to_numpy(dtype=float)
    reference_fitted = df["restricted_fitted_z"].to_numpy(dtype=float)
    parameter_difference = params - reference_params
    fitted_difference = fitted - reference_fitted

    current_metrics = current_summary["regime_restricted_metrics"]
    comparison = pd.DataFrame([
        {
            "model": "unified_regime_local_adaptive_gwr",
            "fit_structure": "one 159-county fit; within-regime adaptive distance; k_i=regime size",
            **metrics,
        },
        {
            "model": "six_independent_regime_gwr",
            "fit_structure": "six separate BasicGWR fits; each selected bandwidth=regime size",
            "rss": current_metrics["rss"],
            "rmse": current_metrics["rmse"],
            "mae": current_metrics["mae"],
            "r2": current_metrics["r2"],
            "trace_hat": current_metrics["trace_hat"],
            "aicc_fixed_partition": current_metrics["aicc"],
        },
    ])
    comparison.to_csv(OUT / "model_comparison.csv", index=False)

    county = pd.DataFrame({
        "AreaKey": df["AreaKey"],
        "regime": labels,
        "regime_local_k": effective_k,
        "unified_fitted_z": fitted,
        "reference_independent_fitted_z": reference_fitted,
        "delta_fitted": fitted_difference,
    })
    for j, term in enumerate(TERMS):
        county[f"unified_{term}"] = params[:, j]
        county[f"reference_independent_{term}"] = reference_params[:, j]
        county[f"delta_{term}"] = parameter_difference[:, j]
    county.to_csv(OUT / "county_equivalence.csv", index=False)

    summary = {
        "status": "exploratory_unified_regime_local_adaptive_equivalence_test",
        "n_counties": int(len(df)),
        "K": 6,
        "one_model_fit": True,
        "labels_fixed": True,
        "weight_rule": "bisquare(regime-local adaptive distance) * I(same regime)",
        "bandwidth_rule": "k_i = size of focal county's regime",
        "regime_sizes": regime_sizes,
        "unified_metrics": metrics,
        "independent_regime_gwr_metrics": current_metrics,
        "equivalence": {
            "max_abs_parameter_difference": float(np.max(np.abs(parameter_difference))),
            "rmse_parameter_difference": float(np.sqrt(np.mean(parameter_difference**2))),
            "max_abs_fitted_difference": float(np.max(np.abs(fitted_difference))),
            "rmse_fitted_difference": float(np.sqrt(np.mean(fitted_difference**2))),
            "trace_hat_difference": float(metrics["trace_hat"] - current_metrics["trace_hat"]),
            "rss_difference": float(metrics["rss"] - current_metrics["rss"]),
        },
        "interpretation_guard": (
            "This test changes implementation/formulation, not statistical flexibility. "
            "If equivalent, one unified fit can represent the six current regime GWR fits, "
            "but trace(S) and conditional AICc should remain essentially unchanged."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(comparison.to_string(index=False))
    print("\n", json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
