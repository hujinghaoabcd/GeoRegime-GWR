"""Simulation 4 prototype: Georgia topology semi-synthetic GR-GWR experiment.

Uses the real Georgia 159-county UTM coordinates, county polygons and the
previously verified 431-edge Queen graph, while generating synthetic X, y and
known spatially varying regression coefficients.  Thus the geography is real
but the statistical truth is known.

K=4 is fixed for this prototype to isolate recovery of irregular unknown
boundaries from K selection.  True regimes are generated independently of y by
a deterministic multi-source graph growth from four geographically separated
seed counties, guaranteeing Queen-connected irregular regions.
"""
from __future__ import annotations

from collections import deque
from pathlib import Path
import json
import sys

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mgwr.gwr import GWR, MGWR
from mgwr.sel_bw import Sel_BW
from sklearn.metrics import adjusted_rand_score

ROOT = Path(__file__).resolve().parents[2]
SYNTH = Path(__file__).resolve().parent
SRC = ROOT / "src"
for p in (str(SRC), str(SYNTH)):
    if p not in sys.path:
        sys.path.insert(0, p)

import run_canonical_regime_boundary_simulation as base
from georegime_gwr import RegimeAwareGWR

DATA = ROOT / "data" / "raw" / "georgia" / "GData_utm.csv"
SHP = ROOT / "data" / "raw" / "georgia" / "G_utm.shp"
EDGE_FILE = ROOT / "results" / "spatial" / "georgia_queen_w" / "edge_list.csv"
OUT = ROOT / "results" / "synthetic" / "georgia_topology_semisynthetic"

SEED = 20260904
K_TRUE = 4
RHO_X = 0.20
NOISE_SIGMA = 0.75
LAMBDA_BOUNDARY = 0.5
TERMS = ["Intercept", "X1", "X2", "X3"]
SLOPE_TERMS = TERMS[1:]

REGIME_SHIFT = {
    1: np.array([0.0, 1.10, -0.75, 0.40]),
    2: np.array([0.0, -0.90, 0.90, 0.35]),
    3: np.array([0.0, 0.25, -0.10, -1.00]),
    4: np.array([0.0, -0.45, -0.70, 0.95]),
}


def load_geography():
    df = pd.read_csv(DATA)
    edges_df = pd.read_csv(EDGE_FILE)
    gdf = gpd.read_file(SHP)
    if len(df) != 159 or len(gdf) != 159 or len(edges_df) != 431:
        raise RuntimeError("Expected canonical Georgia 159 counties and 431 Queen edges")
    coords = df[["X", "Y"]].to_numpy(dtype=float)
    edges = edges_df[["i", "j"]].to_numpy(dtype=int)
    return df, gdf, coords, edges


def normalized_coords(coords):
    lo = coords.min(axis=0)
    span = np.maximum(coords.max(axis=0) - lo, 1e-12)
    return (coords - lo) / span


def choose_geographic_seeds(coords):
    """Pick counties nearest four normalized-corner targets."""
    z = normalized_coords(coords)
    targets = np.array([[0.12, 0.15], [0.15, 0.85], [0.85, 0.82], [0.87, 0.18]])
    seeds = []
    for target in targets:
        order = np.argsort(np.sum((z - target) ** 2, axis=1))
        seed = next(int(i) for i in order if int(i) not in seeds)
        seeds.append(seed)
    return seeds


def graph_grow_labels(n, edges, seeds):
    """Simultaneous multi-source BFS; each regime is connected by construction."""
    adj = base.adjacency_from_edges(n, edges)
    labels = np.zeros(n, dtype=int)
    q = deque()
    for label, seed in enumerate(seeds, start=1):
        labels[seed] = label
        q.append(seed)
    while q:
        u = q.popleft()
        for v in sorted(adj[u]):
            if labels[v] == 0:
                labels[v] = labels[u]
                q.append(v)
    if np.any(labels == 0):
        raise RuntimeError("Graph growth did not label every county")
    return labels


def true_raw_coefficients(coords, regimes):
    z = normalized_coords(coords)
    x, y = z[:, 0], z[:, 1]
    smooth = np.column_stack([
        0.18 * np.sin(np.pi * x) * np.cos(np.pi * y),
        0.24 * np.sin(np.pi * x) * np.cos(np.pi * y),
        0.22 * (x - 0.5) + 0.10 * np.sin(np.pi * y),
        0.22 * (y - 0.5) + 0.10 * np.cos(np.pi * x),
    ])
    shifts = np.vstack([REGIME_SHIFT[int(r)] for r in regimes])
    return smooth + shifts


def generate_data(coords, true):
    rng = np.random.default_rng(SEED)
    cov = np.full((3, 3), RHO_X, dtype=float)
    np.fill_diagonal(cov, 1.0)
    X_raw = rng.multivariate_normal(np.zeros(3), cov, size=len(true))
    beta_raw = true_raw_coefficients(coords, true)
    signal = beta_raw[:, 0] + np.sum(beta_raw[:, 1:] * X_raw, axis=1)
    y_raw = signal + rng.normal(0.0, NOISE_SIGMA, size=len(true))

    xm = X_raw.mean(axis=0)
    xs = X_raw.std(axis=0, ddof=0)
    ym = float(y_raw.mean())
    ys = float(y_raw.std(ddof=0))
    X = (X_raw - xm) / xs
    y = (y_raw - ym) / ys

    beta_z = np.empty_like(beta_raw)
    beta_z[:, 1:] = beta_raw[:, 1:] * xs[None, :] / ys
    beta_z[:, 0] = (
        beta_raw[:, 0] + np.sum(beta_raw[:, 1:] * xm[None, :], axis=1) - ym
    ) / ys
    return X, y, beta_z, beta_raw


def slope_rmse(params, truth):
    e = params[:, 1:] - truth[:, 1:]
    return float(np.sqrt(np.mean(e * e)))


def fit_rmse(y, fitted):
    return float(np.sqrt(np.mean((y - fitted) ** 2)))


def jump_ratio(params, truth, true, edges):
    ee = edges[base.boundary_edges(true, edges)]
    tj = np.abs(truth[ee[:, 0], 1:] - truth[ee[:, 1], 1:])
    ej = np.abs(params[ee[:, 0], 1:] - params[ee[:, 1], 1:])
    return float(np.mean(ej) / np.mean(tj))


def plot_regimes(gdf, df, true, initial, refined):
    key_col = next(c for c in gdf.columns if str(c).lower() == "areakey")
    gm = gdf.copy()
    gm["_key"] = gm[key_col].astype(float).astype(int).astype(str)
    tmp = pd.DataFrame({
        "_key": df["AreaKey"].astype(float).astype(int).astype(str),
        "true": true,
        "initial": base.align_labels(initial, true),
        "refined": base.align_labels(refined, true),
    })
    gm = gm.merge(tmp, on="_key", how="left", validate="1:1")
    fig, axes = plt.subplots(1, 3, figsize=(16, 6))
    for ax, col, title in zip(
        axes, ["true", "initial", "refined"],
        ["True graph-grown regimes", "Initial Queen-Ward", "Refined GR-GWR"],
    ):
        gm.plot(column=col, categorical=True, cmap="tab10", edgecolor="black", linewidth=0.35, ax=ax)
        ax.set_title(title)
        ax.set_axis_off()
    fig.suptitle("Georgia-topology semi-synthetic regime recovery")
    fig.tight_layout()
    fig.savefig(OUT / "regime_recovery.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "regime_recovery.svg", bbox_inches="tight")
    plt.close(fig)


def plot_x1(gdf, df, beta_true, models):
    key_col = next(c for c in gdf.columns if str(c).lower() == "areakey")
    gm = gdf.copy()
    gm["_key"] = gm[key_col].astype(float).astype(int).astype(str)
    vals = pd.DataFrame({"_key": df["AreaKey"].astype(float).astype(int).astype(str), "True": beta_true[:, 1]})
    for name, params in models.items():
        vals[name] = params[:, 1]
    gm = gm.merge(vals, on="_key", how="left", validate="1:1")
    order = ["True", "GWR", "MGWR", "GR-GWR", "Oracle GR-GWR"]
    vmax = max(float(np.nanmax(np.abs(gm[c].to_numpy(dtype=float)))) for c in order)
    fig, axes = plt.subplots(1, 5, figsize=(22, 5))
    for ax, col in zip(axes, order):
        gm.plot(column=col, cmap="coolwarm", vmin=-vmax, vmax=vmax, edgecolor="black", linewidth=0.25, ax=ax)
        ax.set_title(col)
        ax.set_axis_off()
    fig.suptitle("Georgia topology: X1 true and recovered standardized coefficients")
    fig.tight_layout()
    fig.savefig(OUT / "X1_coefficient_recovery.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "X1_coefficient_recovery.svg", bbox_inches="tight")
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    df, gdf, coords, edges = load_geography()
    seeds = choose_geographic_seeds(coords)
    true = graph_grow_labels(len(df), edges, seeds)
    adj = base.adjacency_from_edges(len(df), edges)
    if not base.all_regimes_connected(true, adj):
        raise RuntimeError("True graph-grown regimes are not connected")
    true_sizes = [int(np.sum(true == r)) for r in range(1, K_TRUE + 1)]
    if min(true_sizes) < 12:
        raise RuntimeError(f"Graph-grown true regime too small: {true_sizes}")

    X, y, beta_true, beta_raw = generate_data(coords, true)
    yc = y[:, None]

    # Standard GWR comparator.
    sel = Sel_BW(coords, yc, X, fixed=False, kernel="bisquare", constant=True)
    gwr_bw = int(round(float(sel.search())))
    gwr = GWR(coords, yc, X, gwr_bw, fixed=False, kernel="bisquare", constant=True).fit()
    gwr_params = np.asarray(gwr.params, dtype=float)
    gwr_fit = np.asarray(gwr.predy, dtype=float).reshape(-1)

    # Estimated GR-GWR. K is fixed, true labels are hidden.
    base.K_TRUE = K_TRUE
    base.LAMBDA_BOUNDARY = LAMBDA_BOUNDARY
    base.MIN_REGIME_N = 6
    base.MAX_REFINEMENT_ITER = 20
    initial = base.ward_initial_partition(gwr_params[:, 1:], edges)
    refined, gr, history, stop = base.refine_labels(initial, X, y, coords, edges, adj)

    # Oracle same refit primitive with true labels.
    oracle = RegimeAwareGWR(bandwidth="regime_size", kernel="bisquare").fit(X, y, coords, true)

    # MGWR comparator.
    msel = Sel_BW(coords, yc, X, multi=True, fixed=False, kernel="bisquare", constant=True)
    bws = np.asarray(msel.search(multi_bw_min=[2], tol_multi=1e-5, max_iter_multi=100, verbose=False), dtype=float)
    mg = MGWR(coords, yc, X, msel, fixed=False, kernel="bisquare", constant=True).fit()
    mg_params = np.asarray(mg.params, dtype=float)
    mg_fit = np.asarray(mg.predy, dtype=float).reshape(-1)

    boundary_rows = []
    for stage, labels in [("initial_ward", initial), ("refined_grgwr", refined)]:
        boundary_rows.append({
            "stage": stage,
            "ari": float(adjusted_rand_score(true, labels)),
            **base.boundary_scores(true, labels, edges),
            "regime_sizes": ";".join(str(int(np.sum(labels == r))) for r in sorted(np.unique(labels))),
        })
    boundary_df = pd.DataFrame(boundary_rows)
    boundary_df.to_csv(OUT / "boundary_recovery.csv", index=False)
    history.to_csv(OUT / "refinement_history.csv", index=False)

    model_params = {
        "GWR": gwr_params,
        "MGWR": mg_params,
        "GR-GWR": gr.parameters_,
        "Oracle GR-GWR": oracle.parameters_,
    }
    rows = []
    for name, params, fitted in [
        ("GWR", gwr_params, gwr_fit),
        ("MGWR", mg_params, mg_fit),
        ("GR-GWR", gr.parameters_, gr.fitted_values_),
        ("Oracle GR-GWR", oracle.parameters_, oracle.fitted_values_),
    ]:
        rows.append({
            "model": name,
            "fit_rmse": fit_rmse(y, fitted),
            "slope_coef_rmse": slope_rmse(params, beta_true),
            "true_boundary_jump_recovery": jump_ratio(params, beta_true, true, edges),
        })
    metrics = pd.DataFrame(rows)
    metrics.to_csv(OUT / "model_recovery_metrics.csv", index=False)

    county = pd.DataFrame({
        "AreaKey": df["AreaKey"],
        "true_regime": true,
        "initial_regime": initial,
        "refined_regime": refined,
    })
    for j, term in enumerate(TERMS):
        county[f"true_{term}"] = beta_true[:, j]
        county[f"gwr_{term}"] = gwr_params[:, j]
        county[f"mgwr_{term}"] = mg_params[:, j]
        county[f"grgwr_{term}"] = gr.parameters_[:, j]
        county[f"oracle_{term}"] = oracle.parameters_[:, j]
    county.to_csv(OUT / "county_results.csv", index=False)

    plot_regimes(gdf, df, true, initial, refined)
    plot_x1(gdf, df, beta_true, model_params)

    summary = {
        "status": "single_run_georgia_topology_semisynthetic_prototype",
        "seed": SEED,
        "n": int(len(df)),
        "queen_edges": int(len(edges)),
        "K_true": K_TRUE,
        "K_fixed_but_true_boundaries_hidden": True,
        "true_regime_generation": "deterministic four-source Queen-graph BFS from geographically separated seeds, independent of y",
        "seed_indices": [int(v) for v in seeds],
        "seed_areakeys": [int(df.iloc[v]["AreaKey"]) for v in seeds],
        "true_regime_sizes": true_sizes,
        "rho_X": RHO_X,
        "noise_sigma_raw": NOISE_SIGMA,
        "lambda_working_baseline": LAMBDA_BOUNDARY,
        "gwr_bandwidth": gwr_bw,
        "mgwr_bandwidths": [int(round(v)) for v in bws],
        "initial_ari": float(boundary_df.loc[boundary_df.stage == "initial_ward", "ari"].iloc[0]),
        "refined_ari": float(boundary_df.loc[boundary_df.stage == "refined_grgwr", "ari"].iloc[0]),
        "initial_boundary_f1": float(boundary_df.loc[boundary_df.stage == "initial_ward", "f1"].iloc[0]),
        "refined_boundary_f1": float(boundary_df.loc[boundary_df.stage == "refined_grgwr", "f1"].iloc[0]),
        "refinement_stop_reason": stop,
        "accepted_refinement_iterations": int(max(0, len(history) - 1)),
        "metrics": rows,
        "warning": "Single realization with K fixed. This tests irregular real geography/topology, not K selection or final Monte Carlo performance."
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print("\nBoundary recovery:\n", boundary_df.to_string(index=False))
    print("\nModel recovery:\n", metrics.to_string(index=False))


if __name__ == "__main__":
    main()
