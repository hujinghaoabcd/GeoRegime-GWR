"""Compare ordinary GWR, external MGWR, and fixed-K6 regime-restricted GWR.

This is an exploratory external benchmark.  MGWR is not part of GR-GWR; it is
used only as a comparator on the same canonical Georgia data and the same
standardized X/y used by the validated BasicGWR baseline.
"""

from __future__ import annotations

from pathlib import Path
import json
import sys

import numpy as np
import pandas as pd
from mgwr.gwr import MGWR
from mgwr.sel_bw import Sel_BW

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from georegime_gwr.gwr import BasicGWR

DATA = ROOT / "data" / "raw" / "georgia" / "GData_utm.csv"
K6_SUMMARY = ROOT / "results" / "real_data" / "georgia_k6_regime_restricted_gwr" / "summary.json"
OUT = ROOT / "results" / "real_data" / "georgia_mgwr_benchmark"

PREDICTORS = ["PctFB", "PctBlack", "PctRural"]
TERMS = ["Intercept", *PREDICTORS]


def _metrics(y: np.ndarray, fitted: np.ndarray, trace_hat: float, aicc: float) -> dict[str, float]:
    residuals = y - fitted
    rss = float(residuals @ residuals)
    rmse = float(np.sqrt(np.mean(residuals**2)))
    mae = float(np.mean(np.abs(residuals)))
    tss = float(np.sum((y - np.mean(y)) ** 2))
    r2 = float(1.0 - rss / tss)
    return {
        "rss": rss,
        "rmse": rmse,
        "mae": mae,
        "r2": r2,
        "trace_hat": float(trace_hat),
        "aicc": float(aicc),
    }


def _basicgwr_aicc(y: np.ndarray, fitted: np.ndarray, trace_s: float) -> float:
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


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(DATA)
    if len(df) != 159:
        raise RuntimeError(f"Expected 159 Georgia counties; got {len(df)}")

    X = df[PREDICTORS].to_numpy(dtype=float)
    y = df["PctBach"].to_numpy(dtype=float)
    coords = df[["X", "Y"]].to_numpy(dtype=float)

    # Exact same preprocessing as the validated Georgia BasicGWR baseline.
    Xz = (X - X.mean(axis=0)) / X.std(axis=0, ddof=0)
    yz = (y - y.mean()) / y.std(ddof=0)
    yz_col = yz.reshape((-1, 1))

    # Our validated standard-GWR research baseline.
    gwr = BasicGWR(bandwidth="auto", kernel="bisquare", fit_intercept=True).fit(Xz, yz, coords)
    if int(gwr.bandwidth_) != 116:
        raise RuntimeError(f"Expected BasicGWR bandwidth 116; got {gwr.bandwidth_}")
    gwr_trace = float(np.trace(gwr.hat_matrix_))
    gwr_aicc = _basicgwr_aicc(yz, gwr.fitted_values_, gwr_trace)
    gwr_metrics = _metrics(yz, gwr.fitted_values_, gwr_trace, gwr_aicc)

    # External mgwr==2.2.1 MGWR.  Each coefficient, including the intercept,
    # receives its own adaptive bisquare bandwidth via multi-bandwidth search.
    selector = Sel_BW(
        coords,
        yz_col,
        Xz,
        multi=True,
        fixed=False,
        kernel="bisquare",
        constant=True,
    )
    bws = np.asarray(selector.search(multi_bw_min=[2]), dtype=float).reshape(-1)
    if bws.size != len(TERMS):
        raise RuntimeError(f"Expected {len(TERMS)} MGWR bandwidths; got {bws.tolist()}")

    mgwr_model = MGWR(
        coords,
        yz_col,
        Xz,
        selector,
        fixed=False,
        kernel="bisquare",
        constant=True,
    )
    mgwr_results = mgwr_model.fit()

    mgwr_fitted = np.asarray(mgwr_results.predy, dtype=float).reshape(-1)
    mgwr_trace = float(mgwr_results.tr_S)
    mgwr_aicc = float(mgwr_results.aicc)
    mgwr_metrics = _metrics(yz, mgwr_fitted, mgwr_trace, mgwr_aicc)

    params = np.asarray(mgwr_results.params, dtype=float)
    if params.shape != (159, len(TERMS)):
        raise RuntimeError(f"Unexpected MGWR parameter shape: {params.shape}")

    coef = pd.DataFrame({"AreaKey": df["AreaKey"]})
    for j, term in enumerate(TERMS):
        coef[term] = params[:, j]
    coef["fitted_z"] = mgwr_fitted
    coef["residual_z"] = yz - mgwr_fitted
    coef.to_csv(OUT / "mgwr_county_results.csv", index=False)

    bw_table = pd.DataFrame({
        "term": TERMS,
        "bandwidth": bws,
        "bandwidth_fraction_of_n": bws / len(df),
    })
    bw_table.to_csv(OUT / "mgwr_bandwidths.csv", index=False)

    if not K6_SUMMARY.exists():
        raise RuntimeError("K6 regime-restricted GWR summary must exist before MGWR comparison")
    k6 = json.loads(K6_SUMMARY.read_text(encoding="utf-8"))
    k6_metrics = k6["regime_restricted_metrics"]

    comparison = pd.DataFrame([
        {"model": "ordinary_gwr", "bandwidth_spec": "adaptive k=116", **gwr_metrics},
        {
            "model": "mgwr",
            "bandwidth_spec": ";".join(f"{term}={int(round(bw))}" for term, bw in zip(TERMS, bws)),
            **mgwr_metrics,
        },
        {
            "model": "k6_regime_restricted_gwr",
            "bandwidth_spec": "per-regime auto; all selected regime size",
            **{key: float(k6_metrics[key]) for key in ["rss", "rmse", "mae", "r2", "trace_hat", "aicc"]},
        },
    ])
    comparison.to_csv(OUT / "model_comparison.csv", index=False)

    summary = {
        "status": "exploratory_external_mgwr_benchmark",
        "mgwr_version": "2.2.1",
        "n_counties": 159,
        "response": "PctBach",
        "predictors": PREDICTORS,
        "standardization": "global X/y z-score, ddof=0",
        "kernel": "bisquare",
        "adaptive": True,
        "ordinary_gwr_bandwidth": 116,
        "mgwr_bandwidths": {term: int(round(bw)) for term, bw in zip(TERMS, bws)},
        "ordinary_gwr": gwr_metrics,
        "mgwr": mgwr_metrics,
        "k6_regime_restricted_gwr": {key: float(k6_metrics[key]) for key in ["rss", "rmse", "mae", "r2", "trace_hat", "aicc"]},
        "delta_mgwr_minus_gwr": {key: float(mgwr_metrics[key] - gwr_metrics[key]) for key in gwr_metrics},
        "delta_k6_minus_mgwr": {key: float(k6_metrics[key] - mgwr_metrics[key]) for key in ["rss", "rmse", "mae", "r2", "trace_hat", "aicc"]},
        "interpretation_guard": (
            "This is an in-sample exploratory benchmark. MGWR is an external comparator. "
            "The K6 AICc treats the discovered partition as fixed and therefore does not include partition-selection complexity."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(json.dumps(summary, indent=2))
    print("\nModel comparison:\n", comparison.to_string(index=False))
    print("\nMGWR bandwidths:\n", bw_table.to_string(index=False))


if __name__ == "__main__":
    main()
