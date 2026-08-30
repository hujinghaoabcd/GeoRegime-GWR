"""Validate and compare standard-GWR bandwidth policies on Georgia.

Standard GWR only; MGWR is never imported or run.

Checks
------
1. External mgwr 2.2.1 default adaptive-bisquare AICc search returns 117.
2. BasicGWR research default (exhaustive integer AICc) returns 116.
3. BasicGWR(search_strategy='mgwr_golden') returns 117.
4. Compatibility-mode final GWR matches mgwr.GWR to machine precision.
5. Persist a descriptive comparison of strict-adaptive, mgwr-compatible adaptive,
   and fixed-distance golden-section fits on the same Georgia data.
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


def mean_abs_diff(a, b) -> float:
    return float(np.mean(np.abs(np.asarray(a, dtype=float) - np.asarray(b, dtype=float))))


def rmse_diff(a, b) -> float:
    d = np.asarray(a, dtype=float) - np.asarray(b, dtype=float)
    return float(np.sqrt(np.mean(d * d)))


def correlation(a, b) -> float:
    a = np.asarray(a, dtype=float).reshape(-1)
    b = np.asarray(b, dtype=float).reshape(-1)
    if np.std(a) <= 1e-15 or np.std(b) <= 1e-15:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


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


def fit_diagnostics(model: BasicGWR, y: np.ndarray) -> dict:
    y = np.asarray(y, dtype=float).reshape(-1)
    residuals = np.asarray(model.residuals_, dtype=float).reshape(-1)
    rss = float(np.sum(residuals**2))
    tss = float(np.sum((y - np.mean(y)) ** 2))
    trace_s = float(np.trace(model.hat_matrix_))
    return {
        "bandwidth": float(model.bandwidth_),
        "aicc": gwr_aicc(y, model.fitted_values_, trace_s),
        "rss": rss,
        "rmse": float(np.sqrt(np.mean(residuals**2))),
        "mae": float(np.mean(np.abs(residuals))),
        "r2": float(1.0 - rss / tss),
        "trace_hat": trace_s,
    }


def pairwise_fit_difference(a: BasicGWR, b: BasicGWR, names: list[str]) -> dict:
    out = {
        "fitted": {
            "mean_abs_difference": mean_abs_diff(a.fitted_values_, b.fitted_values_),
            "rmse_difference": rmse_diff(a.fitted_values_, b.fitted_values_),
            "max_abs_difference": max_abs(a.fitted_values_, b.fitted_values_),
            "correlation": correlation(a.fitted_values_, b.fitted_values_),
        },
        "parameters_overall": {
            "mean_abs_difference": mean_abs_diff(a.parameters_, b.parameters_),
            "rmse_difference": rmse_diff(a.parameters_, b.parameters_),
            "max_abs_difference": max_abs(a.parameters_, b.parameters_),
            "correlation": correlation(a.parameters_, b.parameters_),
        },
        "parameters_by_term": {},
    }
    for j, name in enumerate(names):
        av = a.parameters_[:, j]
        bv = b.parameters_[:, j]
        d = av - bv
        out["parameters_by_term"][name] = {
            "mean_signed_difference": float(np.mean(d)),
            "mean_abs_difference": float(np.mean(np.abs(d))),
            "rmse_difference": float(np.sqrt(np.mean(d**2))),
            "max_abs_difference": float(np.max(np.abs(d))),
            "correlation": correlation(av, bv),
        }
    return out


def main() -> None:
    df = pd.read_csv(DATA)
    y = df["PctBach"].to_numpy(dtype=float).reshape(-1, 1)
    X = df[["PctFB", "PctBlack", "PctRural"]].to_numpy(dtype=float)
    coords = df[["X", "Y"]].to_numpy(dtype=float)

    Xz = (X - X.mean(axis=0)) / X.std(axis=0, ddof=0)
    yz = (y - y.mean(axis=0)) / y.std(axis=0, ddof=0)
    yz1 = yz.reshape(-1)

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

    # Research default: adaptive exhaustive integer AICc policy.
    strict = BasicGWR(
        bandwidth="auto",
        kernel="bisquare",
        fit_intercept=True,
    ).fit(Xz, yz1, coords)

    # Historical compatibility mode: reproduce mgwr 2.2.1 adaptive search.
    compat = BasicGWR(
        bandwidth="auto",
        kernel="bisquare",
        fit_intercept=True,
        adaptive=True,
        search_strategy="mgwr_golden",
    ).fit(Xz, yz1, coords)

    # Fixed-distance research path: continuous golden-section AICc search.
    fixed = BasicGWR(
        bandwidth="auto",
        kernel="bisquare",
        fit_intercept=True,
        adaptive=False,
    ).fit(Xz, yz1, coords)

    ref_params = np.asarray(reference.params)
    ref_fitted = np.asarray(reference.predy).reshape(-1)
    ref_residuals = np.asarray(reference.resid_response).reshape(-1)
    ref_hat = np.asarray(reference.S)

    strict_aicc = gwr_aicc(yz1, strict.fitted_values_, float(np.trace(strict.hat_matrix_)))
    compat_aicc = gwr_aicc(yz1, compat.fitted_values_, float(np.trace(compat.hat_matrix_)))

    names = ["Intercept", "PctFB", "PctBlack", "PctRural"]
    diagnostics = {
        "strict_adaptive_exhaustive": fit_diagnostics(strict, yz1),
        "compat_adaptive_mgwr_golden": fit_diagnostics(compat, yz1),
        "fixed_distance_golden": fit_diagnostics(fixed, yz1),
    }

    strategy_comparison = {
        "strict_vs_compat": pairwise_fit_difference(strict, compat, names),
        "strict_vs_fixed": pairwise_fit_difference(strict, fixed, names),
        "compat_vs_fixed": pairwise_fit_difference(compat, fixed, names),
    }

    metrics = {
        "n": int(len(df)),
        "external_mgwr_selected_bandwidth": reference_bandwidth,
        "strict_exhaustive_selected_bandwidth": int(strict.bandwidth_),
        "compat_mgwr_golden_selected_bandwidth": int(compat.bandwidth_),
        "fixed_golden_selected_bandwidth": float(fixed.bandwidth_),
        "strict_search_strategy": strict.bandwidth_search_.strategy,
        "compat_search_strategy": compat.bandwidth_search_.strategy,
        "fixed_search_strategy": fixed.bandwidth_search_.strategy,
        "strict_search_range": list(strict.bandwidth_search_.search_range),
        "compat_search_range": list(compat.bandwidth_search_.search_range),
        "fixed_search_range": list(fixed.bandwidth_search_.search_range),
        "strict_best_aicc": float(strict.bandwidth_search_.score),
        "compat_search_best_aicc": float(compat.bandwidth_search_.score),
        "fixed_search_best_aicc": float(fixed.bandwidth_search_.score),
        "strict_final_aicc": strict_aicc,
        "compat_final_aicc": compat_aicc,
        "fixed_final_aicc": diagnostics["fixed_distance_golden"]["aicc"],
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
        "fit_diagnostics_by_strategy": diagnostics,
        "pairwise_strategy_differences": strategy_comparison,
    }

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
        "strict_fitted": strict.fitted_values_,
        "compat_fitted": compat.fitted_values_,
        "fixed_fitted": fixed.fitted_values_,
        "strict_minus_compat_fitted": strict.fitted_values_ - compat.fitted_values_,
        "strict_minus_fixed_fitted": strict.fitted_values_ - fixed.fitted_values_,
        "compat_minus_fixed_fitted": compat.fitted_values_ - fixed.fitted_values_,
        "reference_fitted": ref_fitted,
        "compat_minus_reference_fitted": compat.fitted_values_ - ref_fitted,
    })
    for j, name in enumerate(names):
        comparison[f"strict_{name}"] = strict.parameters_[:, j]
        comparison[f"compat_{name}"] = compat.parameters_[:, j]
        comparison[f"fixed_{name}"] = fixed.parameters_[:, j]
        comparison[f"strict_minus_compat_{name}"] = strict.parameters_[:, j] - compat.parameters_[:, j]
        comparison[f"strict_minus_fixed_{name}"] = strict.parameters_[:, j] - fixed.parameters_[:, j]
        comparison[f"reference_{name}"] = ref_params[:, j]
    comparison.to_csv(OUT_DIR / "pointwise_comparison.csv", index=False)

    pd.DataFrame(
        strict.bandwidth_search_.search_trace,
        columns=["bandwidth", "aicc"],
    ).to_csv(OUT_DIR / "strict_exhaustive_bandwidth_curve.csv", index=False)

    pd.DataFrame(
        compat.bandwidth_search_.search_trace,
        columns=["bandwidth", "aicc"],
    ).to_csv(OUT_DIR / "mgwr_compatible_search_trace.csv", index=False)

    pd.DataFrame(
        fixed.bandwidth_search_.search_trace,
        columns=["bandwidth", "aicc"],
    ).to_csv(OUT_DIR / "fixed_golden_search_trace.csv", index=False)

    pd.DataFrame([
        {"strategy": key, **value}
        for key, value in diagnostics.items()
    ]).to_csv(OUT_DIR / "strategy_diagnostics.csv", index=False)

    pairwise_rows = []
    for pair, detail in strategy_comparison.items():
        pairwise_rows.append({"pair": pair, "term": "fitted", **detail["fitted"]})
        pairwise_rows.append({"pair": pair, "term": "parameters_overall", **detail["parameters_overall"]})
        for term, vals in detail["parameters_by_term"].items():
            pairwise_rows.append({"pair": pair, "term": term, **vals})
    pd.DataFrame(pairwise_rows).to_csv(OUT_DIR / "strategy_pairwise_differences.csv", index=False)

    readme = [
        "# Georgia standard-GWR bandwidth validation",
        "",
        "MGWR is not imported or run.",
        "",
        f"- external mgwr 2.2.1 bandwidth: {reference_bandwidth}",
        f"- strict exhaustive BasicGWR bandwidth: {int(strict.bandwidth_)}",
        f"- mgwr-compatible BasicGWR bandwidth: {int(compat.bandwidth_)}",
        f"- fixed-distance golden BasicGWR bandwidth: {float(fixed.bandwidth_):.6f}",
        f"- strict final AICc: {strict_aicc:.12f}",
        f"- compatible final AICc: {compat_aicc:.12f}",
        f"- fixed final AICc: {diagnostics['fixed_distance_golden']['aicc']:.12f}",
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
