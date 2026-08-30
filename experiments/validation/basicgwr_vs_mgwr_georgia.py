"""End-to-end validate BasicGWR against mgwr.gwr.GWR on Georgia data.

This is a standard-GWR-only validation. It checks BOTH stages:

1. automatic adaptive-bisquare AICc bandwidth selection;
2. final GWR fit at the selected bandwidth.

The candidate uses the lightweight PyGWRx-derived bandwidth selector integrated
into ``BasicGWR(bandwidth='auto')``. The external reference uses only
``mgwr.sel_bw.Sel_BW`` and ``mgwr.gwr.GWR``; MGWR is never imported or run.
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

    # Exact same preprocessing as the canonical Georgia GWR example.
    Xz = (X - X.mean(axis=0)) / X.std(axis=0, ddof=0)
    yz = (y - y.mean(axis=0)) / y.std(axis=0, ddof=0)

    # External reference: standard GWR bandwidth selection only.
    reference_selector = Sel_BW(coords, yz, Xz, fixed=False, kernel="bisquare")
    reference_bandwidth = int(reference_selector.search(bw_min=2))

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

    # Candidate: bandwidth is NOT supplied. BasicGWR must recover it itself.
    candidate = BasicGWR(
        bandwidth="auto",
        kernel="bisquare",
        fit_intercept=True,
    ).fit(Xz, yz.reshape(-1), coords)

    candidate_bandwidth = int(candidate.bandwidth_)
    ref_params = np.asarray(reference.params)
    ref_fitted = np.asarray(reference.predy).reshape(-1)
    ref_residuals = np.asarray(reference.resid_response).reshape(-1)
    ref_hat = np.asarray(reference.S)

    candidate_aicc = gwr_aicc(
        yz,
        candidate.fitted_values_,
        float(np.trace(candidate.hat_matrix_)),
    )

    metrics = {
        "n": int(len(df)),
        "external_selected_bandwidth": reference_bandwidth,
        "basicgwr_selected_bandwidth": candidate_bandwidth,
        "bandwidths_identical": bool(candidate_bandwidth == reference_bandwidth),
        "basicgwr_search_range": list(candidate.bandwidth_search_.search_range),
        "basicgwr_search_best_aicc": float(candidate.bandwidth_search_.score),
        "parameter_shape": list(candidate.parameters_.shape),
        "max_abs_parameter_difference": max_abs(candidate.parameters_, ref_params),
        "rmse_parameter_difference": rmse_diff(candidate.parameters_, ref_params),
        "max_abs_fitted_difference": max_abs(candidate.fitted_values_, ref_fitted),
        "rmse_fitted_difference": rmse_diff(candidate.fitted_values_, ref_fitted),
        "max_abs_residual_difference": max_abs(candidate.residuals_, ref_residuals),
        "max_abs_hat_matrix_difference": max_abs(candidate.hat_matrix_, ref_hat),
        "basicgwr_rss": float(np.sum(candidate.residuals_ ** 2)),
        "mgwr_gwr_rss": float(np.sum(ref_residuals ** 2)),
        "rss_difference": float(np.sum(candidate.residuals_ ** 2) - np.sum(ref_residuals ** 2)),
        "basicgwr_trace_hat": float(np.trace(candidate.hat_matrix_)),
        "mgwr_gwr_trace_hat": float(np.trace(ref_hat)),
        "basicgwr_final_aicc": candidate_aicc,
        "mgwr_gwr_aicc": float(reference.aicc),
        "aicc_difference": float(candidate_aicc - float(reference.aicc)),
    }

    names = ["Intercept", "PctFB", "PctBlack", "PctRural"]
    metrics["max_abs_parameter_difference_by_term"] = {
        name: max_abs(candidate.parameters_[:, j], ref_params[:, j])
        for j, name in enumerate(names)
    }

    # Exact integer bandwidth equality is mandatory. Floating-point matrix
    # operations are required to agree to near-machine precision.
    numerical_tolerance = 1e-12
    metrics["numerical_tolerance"] = numerical_tolerance
    metrics["passes_end_to_end_validation"] = bool(
        metrics["bandwidths_identical"]
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
        "basic_fitted": candidate.fitted_values_,
        "reference_fitted": ref_fitted,
        "fitted_difference": candidate.fitted_values_ - ref_fitted,
    })
    for j, name in enumerate(names):
        comparison[f"basic_{name}"] = candidate.parameters_[:, j]
        comparison[f"reference_{name}"] = ref_params[:, j]
        comparison[f"difference_{name}"] = candidate.parameters_[:, j] - ref_params[:, j]
    comparison.to_csv(OUT_DIR / "pointwise_comparison.csv", index=False)

    search_curve = pd.DataFrame(
        candidate.bandwidth_search_.search_trace,
        columns=["bandwidth", "aicc"],
    )
    search_curve.to_csv(OUT_DIR / "basicgwr_bandwidth_search_curve.csv", index=False)

    readme = [
        "# BasicGWR vs mgwr.GWR — Georgia end-to-end validation",
        "",
        "The external `mgwr` package is used only for standard GWR and its bandwidth selector.",
        "MGWR is not imported or run.",
        "",
        f"- external selected bandwidth: {reference_bandwidth}",
        f"- BasicGWR selected bandwidth: {candidate_bandwidth}",
        f"- bandwidths identical: {metrics['bandwidths_identical']}",
        f"- max |parameter difference|: {metrics['max_abs_parameter_difference']:.3e}",
        f"- max |fitted difference|: {metrics['max_abs_fitted_difference']:.3e}",
        f"- max |hat-matrix difference|: {metrics['max_abs_hat_matrix_difference']:.3e}",
        f"- AICc difference: {metrics['aicc_difference']:.3e}",
        f"- end-to-end pass: {metrics['passes_end_to_end_validation']}",
    ]
    (OUT_DIR / "README.md").write_text("\n".join(readme) + "\n", encoding="utf-8")

    print(json.dumps(metrics, indent=2))
    if not metrics["passes_end_to_end_validation"]:
        raise SystemExit("BasicGWR end-to-end validation failed")


if __name__ == "__main__":
    main()
