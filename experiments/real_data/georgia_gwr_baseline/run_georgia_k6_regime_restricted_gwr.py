"""Exploratory K=6 regime-restricted GWR on the Georgia benchmark.

Purpose
-------
This is the simplest post-segmentation experiment requested during method
exploration.  The K=6 Queen-constrained Ward labels are held fixed.  Each
regime then fits its own ordinary BasicGWR using only observations from that
regime.  No label refinement, boundary penalty, ICM, or iterative repartition
is applied.

Important: the per-regime bandwidth policy used here is exploratory, not a
frozen GR-GWR design decision.  Each regime independently selects an adaptive
bisquare bandwidth by the repository's exhaustive integer AICc policy.
"""

from __future__ import annotations

from pathlib import Path
import json
import sys

import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from georegime_gwr.gwr import BasicGWR

DATA = ROOT / "data" / "raw" / "georgia" / "GData_utm.csv"
SHP = ROOT / "data" / "raw" / "georgia" / "G_utm.shp"
LABELS = ROOT / "results" / "spatial" / "georgia_initial_ward_regimes" / "initial_regime_labels_k2_k15.csv"
OUT = ROOT / "results" / "real_data" / "georgia_k6_regime_restricted_gwr"

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
        "aicc": _aicc(y, fitted, trace_s),
    }


def _norm(values: np.ndarray):
    vmin = float(np.nanmin(values))
    vmax = float(np.nanmax(values))
    if vmin < 0.0 < vmax:
        return TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)
    if np.isclose(vmin, vmax):
        eps = 1e-9 if np.isclose(vmin, 0.0) else abs(vmin) * 1e-6
        return Normalize(vmin=vmin - eps, vmax=vmax + eps)
    return Normalize(vmin=vmin, vmax=vmax)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(DATA)
    labels = pd.read_csv(LABELS)[["AreaKey", "K6"]].copy()
    gdf = gpd.read_file(SHP)

    if len(df) != 159 or len(labels) != 159 or len(gdf) != 159:
        raise RuntimeError("Georgia K6 experiment expects exactly 159 counties")

    df["AreaKey_join"] = df["AreaKey"].map(_key)
    labels["AreaKey_join"] = labels["AreaKey"].map(_key)
    df = df.merge(labels[["AreaKey_join", "K6"]], on="AreaKey_join", how="left", validate="1:1")
    if df["K6"].isna().any():
        raise RuntimeError("K6 labels failed to join to Georgia data")
    regimes = df["K6"].to_numpy(dtype=int)
    if sorted(np.unique(regimes).tolist()) != [1, 2, 3, 4, 5, 6]:
        raise RuntimeError(f"Unexpected K6 labels: {np.unique(regimes).tolist()}")

    X = df[PREDICTORS].to_numpy(dtype=float)
    y = df["PctBach"].to_numpy(dtype=float)
    coords = df[["X", "Y"]].to_numpy(dtype=float)

    # Keep the exact same global standardization as the validated standard-GWR baseline.
    Xz = (X - X.mean(axis=0)) / X.std(axis=0, ddof=0)
    yz = (y - y.mean()) / y.std(ddof=0)

    baseline = BasicGWR(bandwidth="auto", kernel="bisquare", fit_intercept=True).fit(Xz, yz, coords)
    if int(baseline.bandwidth_) != 116:
        raise RuntimeError(f"Expected baseline bandwidth 116; got {baseline.bandwidth_}")
    baseline_metrics = _metrics(yz, baseline.fitted_values_, baseline.hat_matrix_)

    n = len(df)
    p = len(TERMS)
    rr_params = np.full((n, p), np.nan, dtype=float)
    rr_fitted = np.full(n, np.nan, dtype=float)
    rr_hat = np.zeros((n, n), dtype=float)
    regime_rows: list[dict[str, float | int]] = []

    for regime in range(1, 7):
        idx = np.flatnonzero(regimes == regime)
        model = BasicGWR(bandwidth="auto", kernel="bisquare", fit_intercept=True).fit(
            Xz[idx], yz[idx], coords[idx]
        )
        rr_params[idx] = model.parameters_
        rr_fitted[idx] = model.fitted_values_
        rr_hat[np.ix_(idx, idx)] = model.hat_matrix_

        m = _metrics(yz[idx], model.fitted_values_, model.hat_matrix_)
        regime_rows.append(
            {
                "regime": regime,
                "n": int(idx.size),
                "selected_bandwidth": int(model.bandwidth_),
                "bandwidth_fraction_of_regime": float(model.bandwidth_ / idx.size),
                **m,
            }
        )

    if not np.all(np.isfinite(rr_params)) or not np.all(np.isfinite(rr_fitted)):
        raise RuntimeError("Non-finite regime-restricted GWR result")

    rr_metrics = _metrics(yz, rr_fitted, rr_hat)
    comparison = pd.DataFrame(
        [
            {"model": "ordinary_gwr", "bandwidth": 116, **baseline_metrics},
            {"model": "k6_regime_restricted_gwr", "bandwidth": "per-regime auto", **rr_metrics},
        ]
    )
    comparison.to_csv(OUT / "model_comparison.csv", index=False)

    regime_summary = pd.DataFrame(regime_rows)
    regime_summary.to_csv(OUT / "regime_gwr_summary.csv", index=False)

    coef = pd.DataFrame({"AreaKey": df["AreaKey"], "regime": regimes})
    for j, term in enumerate(TERMS):
        coef[f"baseline_{term}"] = baseline.parameters_[:, j]
        coef[f"restricted_{term}"] = rr_params[:, j]
        coef[f"delta_{term}"] = rr_params[:, j] - baseline.parameters_[:, j]
    coef["baseline_fitted_z"] = baseline.fitted_values_
    coef["restricted_fitted_z"] = rr_fitted
    coef["baseline_residual_z"] = yz - baseline.fitted_values_
    coef["restricted_residual_z"] = yz - rr_fitted
    coef.to_csv(OUT / "county_results.csv", index=False)

    # Join to polygons for maps.
    shp_key = next((c for c in gdf.columns if str(c).lower() == "areakey"), None)
    if shp_key is None:
        raise RuntimeError("AreaKey not found in Georgia shapefile")
    gdf["AreaKey_join"] = gdf[shp_key].map(_key)
    map_df = gdf[["AreaKey_join", "geometry"]].merge(
        coef.assign(AreaKey_join=coef["AreaKey"].map(_key)),
        on="AreaKey_join",
        how="left",
        validate="1:1",
    )

    # Figure 1: post-segmentation local coefficients.
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    for ax, term in zip(axes.flat, TERMS):
        col = f"restricted_{term}"
        values = map_df[col].to_numpy(dtype=float)
        map_df.plot(
            column=col,
            cmap="coolwarm",
            norm=_norm(values),
            edgecolor="black",
            linewidth=0.25,
            ax=ax,
            legend=True,
            legend_kwds={"shrink": 0.78},
        )
        ax.set_title(term)
        ax.set_axis_off()
    fig.suptitle("Georgia K=6 regime-restricted GWR local coefficients", fontsize=16)
    fig.tight_layout()
    fig.savefig(OUT / "georgia_k6_regime_restricted_coefficients.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "georgia_k6_regime_restricted_coefficients.svg", bbox_inches="tight")
    plt.close(fig)

    # Figure 2: direct before/after diagnostic.
    slope_delta = np.linalg.norm(rr_params[:, 1:] - baseline.parameters_[:, 1:], axis=1)
    map_df["slope_change_norm"] = slope_delta
    map_df["baseline_abs_resid"] = np.abs(yz - baseline.fitted_values_)
    map_df["restricted_abs_resid"] = np.abs(yz - rr_fitted)

    fig, axes = plt.subplots(1, 3, figsize=(18, 6.2))
    map_df.plot(column="regime", categorical=True, cmap="tab10", edgecolor="black", linewidth=0.3, ax=axes[0])
    axes[0].set_title("Fixed K=6 regimes")
    map_df.plot(column="slope_change_norm", cmap="viridis", edgecolor="black", linewidth=0.25, legend=True, ax=axes[1])
    axes[1].set_title("Local slope change: restricted vs ordinary GWR")
    resid_max = float(max(map_df["baseline_abs_resid"].max(), map_df["restricted_abs_resid"].max()))
    # Plot improvement: positive means restricted GWR reduces absolute residual.
    map_df["abs_residual_improvement"] = map_df["baseline_abs_resid"] - map_df["restricted_abs_resid"]
    imp = map_df["abs_residual_improvement"].to_numpy(dtype=float)
    map_df.plot(column="abs_residual_improvement", cmap="coolwarm", norm=_norm(imp), edgecolor="black", linewidth=0.25, legend=True, ax=axes[2])
    axes[2].set_title("Absolute residual improvement (>0 is better)")
    for ax in axes:
        ax.set_axis_off()
    fig.suptitle("Ordinary GWR vs fixed K=6 regime-restricted GWR", fontsize=16)
    fig.tight_layout()
    fig.savefig(OUT / "georgia_k6_gwr_comparison.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "georgia_k6_gwr_comparison.svg", bbox_inches="tight")
    plt.close(fig)

    summary = {
        "status": "exploratory_fixed_k6_regime_restricted_gwr",
        "n_counties": 159,
        "K": 6,
        "labels_fixed": True,
        "label_refinement_applied": False,
        "boundary_penalty_applied": False,
        "regime_restricted_borrowing": True,
        "cross_regime_weights_forced_to_zero": True,
        "standardization": "global X/y z-score, ddof=0, same as validated Georgia baseline",
        "kernel": "bisquare",
        "baseline_bandwidth": 116,
        "per_regime_bandwidth_policy": "independent adaptive exhaustive integer AICc",
        "per_regime_bandwidths": {str(int(r["regime"])): int(r["selected_bandwidth"]) for r in regime_rows},
        "baseline_metrics": baseline_metrics,
        "regime_restricted_metrics": rr_metrics,
        "delta_restricted_minus_baseline": {
            key: float(rr_metrics[key] - baseline_metrics[key])
            for key in ["rss", "rmse", "mae", "r2", "trace_hat", "aicc"]
        },
        "warning": "In-sample exploratory comparison only. K selection, per-regime bandwidth policy, inference, and spatial CV remain unresolved.",
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(json.dumps(summary, indent=2))
    print("\nPer-regime summary:\n", regime_summary.to_string(index=False))


if __name__ == "__main__":
    main()
