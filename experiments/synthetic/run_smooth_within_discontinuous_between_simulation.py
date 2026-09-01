"""Simulation 2 prototype: smooth within regimes, discontinuous between regimes.

This deliberately makes the GR-GWR task harder than the canonical piecewise-
constant experiment.  True local coefficients vary smoothly inside each true
regime, while regime-specific offsets create genuine jumps at regime borders.

K=3 is fixed for this prototype, but the true boundary locations and labels are
hidden from estimated GR-GWR.  The current working GR-GWR chain is kept intact:
pilot GWR -> Queen-constrained Ward -> lambda=0.5 refinement -> unified
RegimeAwareGWR.  Oracle GR-GWR receives the true labels and therefore separates
boundary-identification error from within-regime local-estimation error.
"""

from __future__ import annotations

from pathlib import Path
import json
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mgwr.gwr import GWR, MGWR
from mgwr.sel_bw import Sel_BW

ROOT = Path(__file__).resolve().parents[2]
SYNTH = Path(__file__).resolve().parent
SRC = ROOT / "src"
for p in (str(SRC), str(SYNTH)):
    if p not in sys.path:
        sys.path.insert(0, p)

import run_canonical_regime_boundary_simulation as base
from georegime_gwr import RegimeAwareGWR

OUT = ROOT / "results" / "synthetic" / "smooth_within_discontinuous_between"
SEED = 20260902
RHO_X = 0.20
NOISE_SIGMA = 0.75
TERMS = base.TERMS
SLOPE_TERMS = base.SLOPE_TERMS

# Regime offsets are deliberately large enough to create a real discontinuity,
# while the smooth fields below remain strong enough that a piecewise-constant
# regional regression would be misspecified.
REGIME_SHIFT = {
    1: np.array([0.0, 1.10, -0.80, 0.50]),
    2: np.array([0.0, -1.00, 0.90, 0.50]),
    3: np.array([0.0, 0.20, -0.10, -0.90]),
}


def true_raw_coefficients(coords: np.ndarray, regimes: np.ndarray) -> np.ndarray:
    """Construct spatially smooth coefficient fields plus discrete regime shifts."""
    x = coords[:, 0]
    y = coords[:, 1]

    smooth = np.column_stack([
        0.25 * np.sin(np.pi * x) * np.cos(np.pi * y),
        0.35 * np.sin(np.pi * x) * np.cos(np.pi * y),
        0.30 * (x - 0.5) + 0.15 * np.sin(np.pi * y),
        0.30 * (y - 0.5) + 0.15 * np.cos(np.pi * x),
    ])
    shifts = np.vstack([REGIME_SHIFT[int(r)] for r in regimes])
    return smooth + shifts


def generate_data(coords: np.ndarray, true: np.ndarray):
    rng = np.random.default_rng(SEED)
    cov = np.full((3, 3), RHO_X, dtype=float)
    np.fill_diagonal(cov, 1.0)
    X_raw = rng.multivariate_normal(np.zeros(3), cov, size=coords.shape[0])

    beta_raw = true_raw_coefficients(coords, true)
    y_signal = beta_raw[:, 0] + np.sum(beta_raw[:, 1:] * X_raw, axis=1)
    y_raw = y_signal + rng.normal(0.0, NOISE_SIGMA, size=coords.shape[0])

    x_mean = X_raw.mean(axis=0)
    x_sd = X_raw.std(axis=0, ddof=0)
    y_mean = float(y_raw.mean())
    y_sd = float(y_raw.std(ddof=0))
    Xz = (X_raw - x_mean) / x_sd
    yz = (y_raw - y_mean) / y_sd

    # Exact transformation for spatially varying coefficients.
    beta_true_z = np.empty_like(beta_raw)
    beta_true_z[:, 1:] = beta_raw[:, 1:] * x_sd[None, :] / y_sd
    beta_true_z[:, 0] = (
        beta_raw[:, 0] + np.sum(beta_raw[:, 1:] * x_mean[None, :], axis=1) - y_mean
    ) / y_sd

    return Xz, yz, beta_true_z, beta_raw, {
        "x_mean": x_mean.tolist(),
        "x_sd": x_sd.tolist(),
        "y_mean": y_mean,
        "y_sd": y_sd,
        "noise_sigma_raw": NOISE_SIGMA,
    }


def within_regime_gradient_rows(
    model: str,
    params: np.ndarray,
    beta_true: np.ndarray,
    true: np.ndarray,
    edges: np.ndarray,
) -> list[dict[str, float | str]]:
    """Measure whether smooth coefficient changes inside regimes are retained."""
    mask = true[edges[:, 0]] == true[edges[:, 1]]
    ee = edges[mask]
    rows = []
    for j, term in enumerate(SLOPE_TERMS, start=1):
        truth = beta_true[ee[:, 1], j] - beta_true[ee[:, 0], j]
        est = params[ee[:, 1], j] - params[ee[:, 0], j]
        rmse = float(np.sqrt(np.mean((est - truth) ** 2)))
        mae = float(np.mean(np.abs(est - truth)))
        if np.std(truth) > 1e-12 and np.std(est) > 1e-12:
            corr = float(np.corrcoef(truth, est)[0, 1])
        else:
            corr = float("nan")
        rows.append({
            "model": model,
            "term": term,
            "within_edges": int(ee.shape[0]),
            "true_gradient_sd": float(np.std(truth, ddof=0)),
            "estimated_gradient_sd": float(np.std(est, ddof=0)),
            "gradient_rmse": rmse,
            "gradient_mae": mae,
            "gradient_correlation": corr,
        })
    return rows


def plot_true_components(coords, true, beta_raw):
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.1))
    arrays = [beta_raw[:, j] for j in range(4)]
    for ax, term, values in zip(axes, TERMS, arrays):
        vmax = float(np.max(np.abs(values)))
        image = ax.imshow(
            values.reshape(base.GRID_N, base.GRID_N),
            origin="lower", cmap="coolwarm", vmin=-vmax, vmax=vmax,
        )
        ax.set_title(f"True {term}")
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(image, ax=ax, shrink=0.72)
    fig.suptitle("Simulation 2 truth: smooth within regimes + discrete boundary shifts")
    fig.tight_layout()
    fig.savefig(OUT / "true_raw_coefficient_surfaces.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "true_raw_coefficient_surfaces.svg", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    # Reuse the exact same lattice, Queen topology, refinement machinery and
    # diagnostics as Simulation 1.  Redirect base plotting/output to this run.
    base.OUT = OUT
    base.SEED = SEED
    base.RHO_X = RHO_X
    base.NOISE_SIGMA = NOISE_SIGMA
    base.LAMBDA_BOUNDARY = 0.5

    coords, true, edges = base.build_grid()
    adjacency = base.adjacency_from_edges(coords.shape[0], edges)
    X, y, beta_true, beta_raw, scaling = generate_data(coords, true)
    y_col = y.reshape(-1, 1)

    # OLS.
    Xd = np.column_stack([np.ones(X.shape[0]), X])
    beta_ols, *_ = np.linalg.lstsq(Xd, y, rcond=None)
    ols_params = np.tile(beta_ols, (X.shape[0], 1))
    ols_fitted = Xd @ beta_ols

    # Standard GWR.
    gwr_selector = Sel_BW(coords, y_col, X, fixed=False, kernel="bisquare", constant=True)
    gwr_bw = int(round(float(gwr_selector.search())))
    gwr_res = GWR(
        coords, y_col, X, gwr_bw,
        fixed=False, kernel="bisquare", constant=True,
    ).fit()
    gwr_params = np.asarray(gwr_res.params, dtype=float)
    gwr_fitted = np.asarray(gwr_res.predy, dtype=float).reshape(-1)

    # Estimated GR-GWR: K is fixed at 3, but true labels/boundaries are hidden.
    initial = base.ward_initial_partition(gwr_params[:, 1:], edges)
    if not base.all_regimes_connected(initial, adjacency):
        raise RuntimeError("Initial Queen-Ward partition is disconnected")
    refined, gr_model, history, stop_reason = base.refine_labels(
        initial, X, y, coords, edges, adjacency
    )
    if not base.all_regimes_connected(refined, adjacency):
        raise RuntimeError("Refined GR-GWR partition is disconnected")

    # Oracle: identical current refit primitive, but true labels are supplied.
    oracle = RegimeAwareGWR(
        bandwidth="regime_size", kernel="bisquare", fit_intercept=True
    ).fit(X, y, coords, true)

    # MGWR comparator.
    mgwr_selector = Sel_BW(
        coords, y_col, X,
        multi=True, fixed=False, kernel="bisquare", constant=True,
    )
    mgwr_bws = np.asarray(
        mgwr_selector.search(
            multi_bw_min=[2], tol_multi=1e-5, max_iter_multi=100, verbose=False
        ), dtype=float,
    ).reshape(-1)
    mgwr_res = MGWR(
        coords, y_col, X, mgwr_selector,
        fixed=False, kernel="bisquare", constant=True,
    ).fit()
    mgwr_params = np.asarray(mgwr_res.params, dtype=float)
    mgwr_fitted = np.asarray(mgwr_res.predy, dtype=float).reshape(-1)

    model_params = {
        "OLS": ols_params,
        "GWR": gwr_params,
        "MGWR": mgwr_params,
        "GR-GWR": gr_model.parameters_,
        "Oracle GR-GWR": oracle.parameters_,
    }
    model_fitted = {
        "OLS": ols_fitted,
        "GWR": gwr_fitted,
        "MGWR": mgwr_fitted,
        "GR-GWR": gr_model.fitted_values_,
        "Oracle GR-GWR": oracle.fitted_values_,
    }

    metrics_df = pd.DataFrame([
        {"model": name, **base.regression_metrics(y, model_fitted[name])}
        for name in ["OLS", "GWR", "MGWR", "GR-GWR", "Oracle GR-GWR"]
    ])
    metrics_df.to_csv(OUT / "model_fit_metrics.csv", index=False)

    true_bnodes = base.boundary_nodes(true, edges)
    coef_rows, jump_rows, gradient_rows = [], [], []
    for name, params in model_params.items():
        coef_rows.extend(base.coefficient_metrics(name, params, beta_true, true_bnodes))
        jump_rows.extend(base.jump_recovery_rows(name, params, beta_true, true, edges))
        gradient_rows.extend(within_regime_gradient_rows(name, params, beta_true, true, edges))
    coef_df = pd.DataFrame(coef_rows)
    jump_df = pd.DataFrame(jump_rows)
    gradient_df = pd.DataFrame(gradient_rows)
    coef_df.to_csv(OUT / "coefficient_recovery.csv", index=False)
    jump_df.to_csv(OUT / "true_boundary_jump_recovery.csv", index=False)
    gradient_df.to_csv(OUT / "within_regime_gradient_recovery.csv", index=False)

    boundary_rows = []
    for stage, labels in [("initial_ward", initial), ("refined_grgwr", refined)]:
        boundary_rows.append({
            "stage": stage,
            "ari": float(base.adjusted_rand_score(true, labels)),
            **base.boundary_scores(true, labels, edges),
            "regime_sizes": ";".join(
                str(int(np.sum(labels == r))) for r in sorted(np.unique(labels))
            ),
        })
    boundary_df = pd.DataFrame(boundary_rows)
    boundary_df.to_csv(OUT / "boundary_recovery.csv", index=False)
    history.to_csv(OUT / "refinement_history.csv", index=False)

    grid_df = pd.DataFrame({
        "id": np.arange(coords.shape[0]),
        "x": coords[:, 0],
        "y": coords[:, 1],
        "true_regime": true,
        "initial_regime": initial,
        "refined_regime": refined,
        "true_boundary_node": true_bnodes,
    })
    for j, term in enumerate(TERMS):
        grid_df[f"true_raw_{term}"] = beta_raw[:, j]
        grid_df[f"true_z_{term}"] = beta_true[:, j]
        for name, params in model_params.items():
            safe = name.lower().replace(" ", "_").replace("-", "")
            grid_df[f"{safe}_{term}"] = params[:, j]
    grid_df.to_csv(OUT / "grid_results.csv", index=False)

    base.plot_regime_recovery(true, initial, refined, boundary_df)
    base.plot_coefficients(model_params, beta_true)
    plot_true_components(coords, true, beta_raw)

    # True within-regime variation quantifies that this is not a piecewise-
    # constant truth disguised as a regional model.
    within_sd = {}
    for j, term in enumerate(SLOPE_TERMS, start=1):
        within_sd[term] = {
            str(r): float(np.std(beta_true[true == r, j], ddof=0))
            for r in np.unique(true)
        }

    summary = {
        "status": "single_run_smooth_within_discontinuous_between_prototype",
        "seed": SEED,
        "grid": f"{base.GRID_N}x{base.GRID_N}",
        "n": int(coords.shape[0]),
        "K_true": int(base.K_TRUE),
        "K_given_but_boundaries_hidden": True,
        "rho_X": RHO_X,
        "noise_sigma_raw": NOISE_SIGMA,
        "lambda_working_baseline": float(base.LAMBDA_BOUNDARY),
        "truth_design": "smooth coefficient fields within regimes plus regime-specific discontinuous shifts",
        "true_regime_sizes": {str(r): int(np.sum(true == r)) for r in np.unique(true)},
        "true_within_regime_coefficient_sd_standardized": within_sd,
        "gwr_bandwidth": gwr_bw,
        "mgwr_bandwidths": {term: int(round(bw)) for term, bw in zip(TERMS, mgwr_bws)},
        "initial_ari": float(boundary_df.loc[0, "ari"]),
        "refined_ari": float(boundary_df.loc[1, "ari"]),
        "initial_boundary_f1": float(boundary_df.loc[0, "f1"]),
        "refined_boundary_f1": float(boundary_df.loc[1, "f1"]),
        "refinement_stop_reason": stop_reason,
        "refinement_accepted_iterations": int(history["iteration"].max()),
        "scaling": scaling,
        "model_fit": metrics_df.to_dict(orient="records"),
        "interpretation_guard": (
            "Single realization only. K=3 is fixed. The truth contains both continuous within-regime "
            "coefficient variation and discrete boundary jumps. Results are behavioral diagnostics before "
            "Monte Carlo and bandwidth/complexity sensitivity work."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(json.dumps(summary, indent=2))
    print("\nBoundary recovery:\n", boundary_df.to_string(index=False))
    print("\nCoefficient recovery:\n", coef_df.to_string(index=False))
    print("\nBoundary jump recovery:\n", jump_df.to_string(index=False))
    print("\nWithin-regime gradient recovery:\n", gradient_df.to_string(index=False))


if __name__ == "__main__":
    main()
