"""Numerically validate BasicGWR against mgwr.gwr.GWR on canonical Georgia data.

This is a GWR-only validation.  The external ``mgwr`` package is used solely
for its standard ``GWR`` implementation; the MGWR model is not imported or run.

Both implementations receive exactly the same standardized data and the same
already-established adaptive bisquare bandwidth (117).  The script compares
all local parameters, fitted values, residuals, hat matrices, and global RSS.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from mgwr.gwr import GWR

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from georegime_gwr.gwr import BasicGWR  # noqa: E402

DATA = ROOT / "data" / "raw" / "georgia" / "GData_utm.csv"
OUT_DIR = ROOT / "results" / "validation" / "basicgwr_vs_mgwr_georgia"
BANDWIDTH = 117


def max_abs(a, b) -> float:
    return float(np.max(np.abs(np.asarray(a, dtype=float) - np.asarray(b, dtype=float))))


def rmse_diff(a, b) -> float:
    d = np.asarray(a, dtype=float) - np.asarray(b, dtype=float)
    return float(np.sqrt(np.mean(d * d)))


def main() -> None:
    df = pd.read_csv(DATA)

    y = df["PctBach"].to_numpy(dtype=float).reshape(-1, 1)
    X = df[["PctFB", "PctBlack", "PctRural"]].to_numpy(dtype=float)
    coords = df[["X", "Y"]].to_numpy(dtype=float)

    # Exact same preprocessing as canonical Georgia GWR reproduction.
    Xz = (X - X.mean(axis=0)) / X.std(axis=0, ddof=0)
    yz = (y - y.mean(axis=0)) / y.std(axis=0, ddof=0)

    reference = GWR(
        coords,
        yz,
        Xz,
        BANDWIDTH,
        fixed=False,
        kernel="bisquare",
        constant=True,
        hat_matrix=True,
        n_jobs=1,
    ).fit()

    candidate = BasicGWR(
        bandwidth=BANDWIDTH,
        kernel="bisquare",
        fit_intercept=True,
    ).fit(Xz, yz.reshape(-1), coords)

    ref_params = np.asarray(reference.params)
    ref_fitted = np.asarray(reference.predy).reshape(-1)
    ref_residuals = np.asarray(reference.resid_response).reshape(-1)
    ref_hat = np.asarray(reference.S)

    metrics = {
        "n": int(len(df)),
        "bandwidth": BANDWIDTH,
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
    }

    # Per-parameter maxima make it easy to diagnose any future drift.
    names = ["Intercept", "PctFB", "PctBlack", "PctRural"]
    metrics["max_abs_parameter_difference_by_term"] = {
        name: max_abs(candidate.parameters_[:, j], ref_params[:, j])
        for j, name in enumerate(names)
    }

    # Machine-level equality is not required across BLAS implementations, but
    # differences should be many orders of magnitude below scientific relevance.
    tolerance = 1e-8
    metrics["tolerance"] = tolerance
    metrics["passes_tolerance"] = bool(
        metrics["max_abs_parameter_difference"] < tolerance
        and metrics["max_abs_fitted_difference"] < tolerance
        and metrics["max_abs_hat_matrix_difference"] < tolerance
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

    readme = [
        "# BasicGWR vs mgwr.GWR — Georgia numerical validation",
        "",
        "The external `mgwr` package is used only for `mgwr.gwr.GWR`; MGWR is not run.",
        "",
        f"- n: {metrics['n']}",
        f"- adaptive bisquare bandwidth: {BANDWIDTH}",
        f"- max |parameter difference|: {metrics['max_abs_parameter_difference']:.3e}",
        f"- max |fitted difference|: {metrics['max_abs_fitted_difference']:.3e}",
        f"- max |hat-matrix difference|: {metrics['max_abs_hat_matrix_difference']:.3e}",
        f"- BasicGWR RSS: {metrics['basicgwr_rss']:.12f}",
        f"- mgwr.GWR RSS: {metrics['mgwr_gwr_rss']:.12f}",
        f"- pass 1e-8 tolerance: {metrics['passes_tolerance']}",
    ]
    (OUT_DIR / "README.md").write_text("\n".join(readme) + "\n", encoding="utf-8")

    print(json.dumps(metrics, indent=2))
    if not metrics["passes_tolerance"]:
        raise SystemExit("BasicGWR does not match mgwr.GWR within tolerance")


if __name__ == "__main__":
    main()
