"""Validate both research and compatibility GWR bandwidth policies on Georgia.

Standard GWR only; MGWR is never imported or run.

Checks
------
1. External mgwr 2.2.1 default adaptive-bisquare AICc search returns 117.
2. BasicGWR research default (exhaustive integer AICc) returns 116.
3. BasicGWR(search_strategy='mgwr_golden') returns 117.
4. Compatibility-mode final GWR matches mgwr.GWR to machine precision.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from mgwr.gwr import GWR
from mgwr.sel_bw import Sel_BW

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from georegime_gwr.gwr import BasicGWR  # noqa: E402

DATA = ROOT / "data" / "raw" / "georgia" / "GData_utm.csv"
OUT_DIR = ROOT / "results" / "validation" / "basicgwr_vs_mgwr_georgia"


def max_abs(a, b) -> float:
    return float(np.max(np.abs(np.asarray(a, dtype=float) - np.asarray(b, dtype=float))))


def rmse_diff(a, b) -> float:
    d = np.asarray(a, dtype=float) - np.asarray(b, dtype=float)
    return float(np.sqrt(np.mean(d * d)))


def gwr_aicc(y: np.ndarray, fitted: np.ndarray, trace_s: float) -> float:
    y = np.asarray(y, dtype=float).reshape(-1)
    fitted = np.asarray(fitted, dtype=float).reshape(-1)
    n = y.size
    rss = float(np.sum((y - fitted) ** 2))
    return float(
        n * np.log(rss / n)
        + n * np.log(2.0 * np.pi)
        + n * (n + trace_s) / (n - 2.0 - trace_s)
    )


def main() -> None:
    df = pd.read_csv(DATA)
    y = df["PctBach"].to_numpy(dtype=float).reshape(-1, 1)
    X = df[["PctFB", "PctBlack", "PctRural"]].to_numpy(dtype=float)
    coords = df[["X", "Y"]].to_numpy(dtype=float)

    Xz = (X - X.mean(axis=0)) / X.std(axis=0, ddof=0)
    yz = (y - y.mean(axis=0)) / y.std(axis=0, ddof=0)

    # External canonical reference: standard GWR only.
    reference_selector = Sel_BW(coords, yz, Xz, fixed=False, kernel="bisquare")
    reference_bandwidth = int(reference_selector.search(criterion="AICc"))
    reference = GWR(
        coords,
        yz,
        Xz,
        reference_bandwidth,
        fixed=False,
        kernel="bisquare",
        constant=True,
        hat_matrix=True,
        n_jobs=1,
    ).fit()

    # Research default: current PyGWRx adaptive exhaustive integer AICc policy.
    strict = BasicGWR(
        bandwidth="auto",
        kernel="bisquare",
        fit_intercept=True,
    ).fit(Xz, yz.reshape(-1), coords)

    # Historical compatibility mode: reproduce mgwr 2.2.1 search behavior.
    compat = BasicGWR(
        bandwidth="auto",
        kernel="bisquare",
        fit_intercept=True,
        adaptive=True,
        search_strategy="mgwr_golden",
    ).fit(Xz, yz.reshape(-1), coords)

    ref_params = np.asarray(reference.params)
    ref_fitted = np.asarray(reference.predy).reshape(-1)
    ref_residuals = np.asarray(reference.resid_response).reshape(-1)
    ref_hat = np.asarray(reference.S)

    strict_aicc = gwr_aicc(yz, strict.fitted_values_, float(np.trace(strict.hat_matrix_)))
    compat_aicc = gwr_aicc(yz, compat.fitted_values_, float(np.trace(compat.hat_matrix_)))

    metrics = {
        "n": int(len(df)),
        "external_mgwr_selected_bandwidth": reference_bandwidth,
        "strict_exhaustive_selected_bandwidth": int(strict.bandwidth_),
        "compat_mgwr_golden_selected_bandwidth": int(compat.bandwidth_),
        "strict_search_strategy": strict.bandwidth_search_.strategy,
        "compat_search_strategy": compat.bandwidth_search_.strategy,
        "strict_search_range": list(strict.bandwidth_search_.search_range),
        "compat_search_range": list(compat.bandwidth_search_.search_range),
        "strict_best_aicc": float(strict.bandwidth_search_.score),
        "compat_search_best_aicc": float(compat.bandwidth_search_.score),
        "strict_final_aicc": strict_aicc,
        "compat_final_aicc": compat_aicc,
        "mgwr_gwr_aicc": float(reference.aicc),
        "strict_better_than_compat_by_aicc": bool(strict_aicc < compat_aicc),
        "compat_bandwidth_matches_mgwr": bool(int(compat.bandwidth_) == reference_bandwidth),
        "parameter_shape": list(compat.parameters_.shape),
        "max_abs_parameter_difference": max_abs(compat.parameters_, ref_params),
        "rmse_parameter_difference": rmse_diff(compat.parameters_, ref_params),
        "max_abs_fitted_difference": max_abs(compat.fitted_values_, ref_fitted),
        "rmse_fitted_difference": rmse_diff(compat.fitted_values_, ref_fitted),
        "max_abs_residual_difference": max_abs(compat.residuals_, ref_residuals),
        "max_abs_hat_matrix_difference": max_abs(compat.hat_matrix_, ref_hat),
        "compat_rss": float(np.sum(compat.residuals_ ** 2)),
        "mgwr_gwr_rss": float(np.sum(ref_residuals ** 2)),
        "rss_difference": float(np.sum(compat.residuals_ ** 2) - np.sum(ref_residuals ** 2)),
        "compat_trace_hat": float(np.trace(compat.hat_matrix_)),
        "mgwr_gwr_trace_hat": float(np.trace(ref_hat)),
        "aicc_difference": float(compat_aicc - float(reference.aicc)),
    }

    names = ["Intercept", "PctFB", "PctBlack", "PctRural"]
    metrics["max_abs_parameter_difference_by_term"] = {
        name: max_abs(compat.parameters_[:, j], ref_params[:, j])
        for j, name in enumerate(names)
    }

    numerical_tolerance = 1e-12
    metrics["numerical_tolerance"] = numerical_tolerance
    metrics["passes_validation"] = bool(
        reference_bandwidth == 117
        and int(strict.bandwidth_) == 116
        and int(compat.bandwidth_) == 117
        and metrics["strict_better_than_compat_by_aicc"]
        and metrics["max_abs_parameter_difference"] < numerical_tolerance
        and metrics["max_abs_fitted_difference"] < numerical_tolerance
        and metrics["max_abs_residual_difference"] < numerical_tolerance
        and metrics["max_abs_hat_matrix_difference"] < numerical_tolerance
        and abs(metrics["rss_difference"]) < numerical_tolerance
        and abs(metrics["aicc_difference"]) < numerical_tolerance
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "summary.json").write_text(
        json.dumps(metrics, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    comparison = pd.DataFrame({
        "AreaKey": df["AreaKey"],
        "compat_fitted": compat.fitted_values_,
        "reference_fitted": ref_fitted,
        "fitted_difference": compat.fitted_values_ - ref_fitted,
    })
    for j, name in enumerate(names):
        comparison[f"compat_{name}"] = compat.parameters_[:, j]
        comparison[f"reference_{name}"] = ref_params[:, j]
        comparison[f"difference_{name}"] = compat.parameters_[:, j] - ref_params[:, j]
    comparison.to_csv(OUT_DIR / "pointwise_comparison.csv", index=False)

    pd.DataFrame(
        strict.bandwidth_search_.search_trace,
        columns=["bandwidth", "aicc"],
    ).to_csv(OUT_DIR / "strict_exhaustive_bandwidth_curve.csv", index=False)

    pd.DataFrame(
        compat.bandwidth_search_.search_trace,
        columns=["bandwidth", "aicc"],
    ).to_csv(OUT_DIR / "mgwr_compatible_search_trace.csv", index=False)

    readme = [
        "# Georgia standard-GWR bandwidth validation",
        "",
        "MGWR is not imported or run.",
        "",
        f"- external mgwr 2.2.1 bandwidth: {reference_bandwidth}",
        f"- strict exhaustive BasicGWR bandwidth: {int(strict.bandwidth_)}",
        f"- mgwr-compatible BasicGWR bandwidth: {int(compat.bandwidth_)}",
        f"- strict final AICc: {strict_aicc:.12f}",
        f"- compatible final AICc: {compat_aicc:.12f}",
        f"- max |compat parameter - mgwr parameter|: {metrics['max_abs_parameter_difference']:.3e}",
        f"- max |compat fitted - mgwr fitted|: {metrics['max_abs_fitted_difference']:.3e}",
        f"- max |compat hat - mgwr hat|: {metrics['max_abs_hat_matrix_difference']:.3e}",
        f"- validation pass: {metrics['passes_validation']}",
    ]
    (OUT_DIR / "README.md").write_text("\n".join(readme) + "\n", encoding="utf-8")

    print(json.dumps(metrics, indent=2))
    if not metrics["passes_validation"]:
        raise SystemExit("BasicGWR bandwidth-policy validation failed")


if __name__ == "__main__":
    main()
