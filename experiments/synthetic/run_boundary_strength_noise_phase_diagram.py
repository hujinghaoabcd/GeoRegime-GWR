"""Simulation 3 prototype: boundary-strength x noise phase diagram.

Purpose
-------
Map the empirical operating region of the current GR-GWR prototype.  The
spatial geometry and K=3 partition are held fixed, while the coefficient jump
magnitude and observation noise are varied.  True boundary locations remain
hidden from estimated GR-GWR.

This is deliberately a single-realization diagnostic grid, not final Monte
Carlo evidence.  It is used to locate transition regions worth replicating.
"""

from __future__ import annotations

from pathlib import Path
import json
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

from mgwr.gwr import GWR, MGWR
from mgwr.sel_bw import Sel_BW

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
SYN = ROOT / "experiments" / "synthetic"
for p in (SRC, SYN):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from georegime_gwr import RegimeAwareGWR
import run_canonical_regime_boundary_simulation as core

OUT = ROOT / "results" / "synthetic" / "boundary_strength_noise_phase_diagram"

SEED = 20260903
RHO_X = 0.20
DELTA_VALUES = [0.0, 0.5, 1.0, 1.5, 2.0]
NOISE_VALUES = [0.35, 0.75, 1.25, 2.0]
TERMS = ["Intercept", "X1", "X2", "X3"]

# Simulation-1 coefficient patterns.  Delta scales deviations from their
# across-regime mean: delta=0 => no coefficient discontinuity; delta=1 => the
# canonical Simulation-1 contrast; larger values strengthen the jump.
REFERENCE = np.vstack([
    np.array([0.0, 1.5, -1.0, 0.5]),
    np.array([0.0, -1.0, 1.0, 0.5]),
    np.array([0.0, 0.5, -0.5, -1.0]),
])
REFERENCE_MEAN = REFERENCE.mean(axis=0)
REFERENCE_DEV = REFERENCE - REFERENCE_MEAN


def make_common_random_components(n: int):
    rng = np.random.default_rng(SEED)
    cov = np.full((3, 3), RHO_X, dtype=float)
    np.fill_diagonal(cov, 1.0)
    X_raw = rng.multivariate_normal(np.zeros(3), cov, size=n)
    eps = rng.normal(0.0, 1.0, size=n)
    return X_raw, eps


def raw_beta_for_delta(true: np.ndarray, delta: float) -> np.ndarray:
    table = REFERENCE_MEAN[None, :] + float(delta) * REFERENCE_DEV
    return table[true - 1]


def standardize_scenario(X_raw, y_raw, beta_raw):
    x_mean = X_raw.mean(axis=0)
    x_sd = X_raw.std(axis=0, ddof=0)
    y_mean = float(y_raw.mean())
    y_sd = float(y_raw.std(ddof=0))
    X = (X_raw - x_mean) / x_sd
    y = (y_raw - y_mean) / y_sd

    beta_z = np.empty_like(beta_raw)
    beta_z[:, 1:] = beta_raw[:, 1:] * x_sd[None, :] / y_sd
    beta_z[:, 0] = (
        beta_raw[:, 0]
        + np.sum(beta_raw[:, 1:] * x_mean[None, :], axis=1)
        - y_mean
    ) / y_sd
    return X, y, beta_z, y_sd


def fit_gwr(coords, X, y):
    y_col = y.reshape(-1, 1)
    selector = Sel_BW(coords, y_col, X, fixed=False, kernel="bisquare", constant=True)
    bw = int(round(float(selector.search())))
    res = GWR(coords, y_col, X, bw, fixed=False, kernel="bisquare", constant=True).fit()
    return bw, np.asarray(res.params, dtype=float), np.asarray(res.predy, dtype=float).reshape(-1)


def fit_mgwr(coords, X, y):
    y_col = y.reshape(-1, 1)
    selector = Sel_BW(
        coords, y_col, X,
        multi=True, fixed=False, kernel="bisquare", constant=True,
    )
    bws = np.asarray(
        selector.search(
            multi_bw_min=[2],
            tol_multi=1e-4,
            max_iter_multi=60,
            verbose=False,
        ),
        dtype=float,
    ).reshape(-1)
    res = MGWR(
        coords, y_col, X, selector,
        fixed=False, kernel="bisquare", constant=True,
    ).fit()
    return bws, np.asarray(res.params, dtype=float), np.asarray(res.predy, dtype=float).reshape(-1)


def rmse(y, fitted) -> float:
    return float(np.sqrt(np.mean((y - fitted) ** 2)))


def slope_beta_rmse(params, beta_true) -> float:
    err = params[:, 1:] - beta_true[:, 1:]
    return float(np.sqrt(np.mean(err**2)))


def mean_jump_recovery(params, beta_true, true, edges) -> float:
    mask = core.boundary_edges(true, edges)
    ee = edges[mask]
    ratios = []
    for j in range(1, 4):
        truth = np.abs(beta_true[ee[:, 0], j] - beta_true[ee[:, 1], j])
        est = np.abs(params[ee[:, 0], j] - params[ee[:, 1], j])
        denom = float(np.mean(truth))
        if denom > 1e-12:
            ratios.append(float(np.mean(est) / denom))
    return float(np.mean(ratios)) if ratios else float("nan")


def scenario(delta, sigma, coords, true, edges, adjacency, X_raw, eps):
    beta_raw = raw_beta_for_delta(true, delta)
    signal = beta_raw[:, 0] + np.sum(beta_raw[:, 1:] * X_raw, axis=1)
    y_raw = signal + float(sigma) * eps
    X, y, beta_true, y_sd = standardize_scenario(X_raw, y_raw, beta_raw)

    gwr_bw, gwr_params, gwr_fitted = fit_gwr(coords, X, y)

    initial = core.ward_initial_partition(gwr_params[:, 1:], edges)
    refined, gr_model, history, stop_reason = core.refine_labels(
        initial, X, y, coords, edges, adjacency
    )
    oracle = RegimeAwareGWR(
        bandwidth="regime_size", kernel="bisquare", fit_intercept=True
    ).fit(X, y, coords, true)

    try:
        mgwr_bws, mgwr_params, mgwr_fitted = fit_mgwr(coords, X, y)
        mgwr_error = None
    except Exception as exc:  # keep the grid usable if one multibandwidth search fails
        mgwr_bws = np.full(4, np.nan)
        mgwr_params = np.full_like(gwr_params, np.nan)
        mgwr_fitted = np.full_like(y, np.nan)
        mgwr_error = f"{type(exc).__name__}: {exc}"

    b = core.boundary_scores(true, refined, edges)
    initial_b = core.boundary_scores(true, initial, edges)

    gwr_rmse = rmse(y, gwr_fitted)
    gr_rmse = rmse(y, gr_model.fitted_values_)
    oracle_rmse = rmse(y, oracle.fitted_values_)
    mgwr_rmse = rmse(y, mgwr_fitted) if np.all(np.isfinite(mgwr_fitted)) else float("nan")

    row = {
        "delta": float(delta),
        "noise_sigma": float(sigma),
        "response_sd_raw": float(y_sd),
        "gwr_bandwidth": int(gwr_bw),
        "mgwr_bandwidths": ";".join("nan" if not np.isfinite(v) else str(int(round(v))) for v in mgwr_bws),
        "mgwr_error": mgwr_error,
        "initial_ari": float(adjusted_rand_score(true, initial)),
        "refined_ari": float(adjusted_rand_score(true, refined)),
        "initial_boundary_f1": float(initial_b["f1"]),
        "refined_boundary_f1": float(b["f1"]),
        "refined_boundary_precision": float(b["precision"]),
        "refined_boundary_recall": float(b["recall"]),
        "refined_boundary_edges": int(b["estimated_boundary_edges"]),
        "refinement_iterations": int(max(history["iteration"])) if len(history) else 0,
        "refinement_stop_reason": stop_reason,
        "gwr_rmse": gwr_rmse,
        "mgwr_rmse": mgwr_rmse,
        "grgwr_rmse": gr_rmse,
        "oracle_grgwr_rmse": oracle_rmse,
        "gr_minus_gwr_rmse": float(gr_rmse - gwr_rmse),
        "gr_vs_gwr_rmse_improvement_pct": float((gwr_rmse - gr_rmse) / gwr_rmse * 100.0),
        "gr_minus_mgwr_rmse": float(gr_rmse - mgwr_rmse) if np.isfinite(mgwr_rmse) else float("nan"),
        "gr_vs_mgwr_rmse_improvement_pct": float((mgwr_rmse - gr_rmse) / mgwr_rmse * 100.0) if np.isfinite(mgwr_rmse) else float("nan"),
        "gwr_beta_rmse": slope_beta_rmse(gwr_params, beta_true),
        "mgwr_beta_rmse": slope_beta_rmse(mgwr_params, beta_true) if np.all(np.isfinite(mgwr_params)) else float("nan"),
        "grgwr_beta_rmse": slope_beta_rmse(gr_model.parameters_, beta_true),
        "oracle_beta_rmse": slope_beta_rmse(oracle.parameters_, beta_true),
        "gwr_jump_recovery": mean_jump_recovery(gwr_params, beta_true, true, edges),
        "mgwr_jump_recovery": mean_jump_recovery(mgwr_params, beta_true, true, edges) if np.all(np.isfinite(mgwr_params)) else float("nan"),
        "grgwr_jump_recovery": mean_jump_recovery(gr_model.parameters_, beta_true, true, edges),
        "oracle_jump_recovery": mean_jump_recovery(oracle.parameters_, beta_true, true, edges),
    }
    return row


def heatmap(df, value, title, filename, fmt=".2f"):
    pivot = df.pivot(index="noise_sigma", columns="delta", values=value).sort_index(ascending=False)
    arr = pivot.to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(8.2, 5.5))
    im = ax.imshow(arr, aspect="auto")
    ax.set_xticks(np.arange(len(pivot.columns)), [str(v) for v in pivot.columns])
    ax.set_yticks(np.arange(len(pivot.index)), [str(v) for v in pivot.index])
    ax.set_xlabel("Boundary strength delta")
    ax.set_ylabel("Noise sigma")
    ax.set_title(title)
    for r in range(arr.shape[0]):
        for c in range(arr.shape[1]):
            val = arr[r, c]
            text = "NA" if not np.isfinite(val) else format(val, fmt)
            ax.text(c, r, text, ha="center", va="center")
    fig.colorbar(im, ax=ax, shrink=0.82)
    fig.tight_layout()
    fig.savefig(OUT / f"{filename}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / f"{filename}.svg", bbox_inches="tight")
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    coords, true, edges = core.build_grid()
    adjacency = core.adjacency_from_edges(len(true), edges)
    X_raw, eps = make_common_random_components(len(true))

    rows = []
    for sigma in NOISE_VALUES:
        for delta in DELTA_VALUES:
            print(f"Running delta={delta}, sigma={sigma}", flush=True)
            rows.append(scenario(delta, sigma, coords, true, edges, adjacency, X_raw, eps))

    df = pd.DataFrame(rows).sort_values(["noise_sigma", "delta"]).reset_index(drop=True)
    df.to_csv(OUT / "phase_diagram_results.csv", index=False)

    heatmap(
        df, "gr_vs_gwr_rmse_improvement_pct",
        "GR-GWR RMSE improvement over GWR (%)",
        "rmse_improvement_vs_gwr", ".1f",
    )
    heatmap(
        df, "gr_vs_mgwr_rmse_improvement_pct",
        "GR-GWR RMSE improvement over MGWR (%)",
        "rmse_improvement_vs_mgwr", ".1f",
    )
    heatmap(df, "refined_ari", "Recovered regime ARI", "refined_ari", ".2f")
    heatmap(df, "refined_boundary_f1", "Recovered boundary F1", "refined_boundary_f1", ".2f")
    heatmap(df, "grgwr_beta_rmse", "GR-GWR slope coefficient RMSE", "grgwr_beta_rmse", ".2f")
    heatmap(df, "grgwr_jump_recovery", "GR-GWR mean true-boundary jump recovery", "grgwr_jump_recovery", ".2f")

    positive_vs_gwr = df[df["gr_vs_gwr_rmse_improvement_pct"] > 0]
    positive_vs_mgwr = df[df["gr_vs_mgwr_rmse_improvement_pct"] > 0]
    delta_positive = {}
    for sigma in NOISE_VALUES:
        sub = df[df["noise_sigma"] == sigma].sort_values("delta")
        winners = sub[sub["gr_vs_mgwr_rmse_improvement_pct"] > 0]
        delta_positive[str(sigma)] = None if winners.empty else float(winners.iloc[0]["delta"])

    summary = {
        "status": "single_realization_boundary_strength_noise_phase_diagram",
        "seed": SEED,
        "grid": "25x25",
        "n": 625,
        "K_fixed": 3,
        "true_boundary_locations_hidden": True,
        "lambda_working_baseline": core.LAMBDA_BOUNDARY,
        "delta_values": DELTA_VALUES,
        "noise_values": NOISE_VALUES,
        "common_random_numbers_across_scenarios": True,
        "cells": int(len(df)),
        "cells_grgwr_beats_gwr_rmse": int(len(positive_vs_gwr)),
        "cells_grgwr_beats_mgwr_rmse": int(len(positive_vs_mgwr)),
        "first_delta_beating_mgwr_by_noise": delta_positive,
        "best_refined_ari": float(df["refined_ari"].max()),
        "worst_refined_ari": float(df["refined_ari"].min()),
        "best_boundary_f1": float(df["refined_boundary_f1"].max()),
        "worst_boundary_f1": float(df["refined_boundary_f1"].min()),
        "interpretation_guard": (
            "Single common-random-number realization per cell. Delta=0 has no coefficient discontinuity, "
            "so the latent three-region geometry is not statistically identifiable and its ARI/F1 should not be interpreted as a recovery target. "
            "This grid is diagnostic only; transition cells require Monte Carlo replication."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
