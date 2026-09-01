"""Canonical synthetic experiment for GR-GWR boundary and coefficient recovery.

Single-run prototype corresponding to Simulation 1 in the paper design.

What is known to the experimenter
---------------------------------
- a 25 x 25 regular lattice (n=625),
- three contiguous true regimes,
- piecewise-constant true regression coefficients.

What is *not* given to estimated GR-GWR
---------------------------------------
- the location of the true boundaries,
- the true regime labels.

For this first prototype K=3 is held fixed on purpose.  K selection is a
separate open problem and is not mixed into the boundary-recovery test.

Compared models
---------------
1. OLS
2. Standard GWR (external mgwr==2.2.1 implementation)
3. MGWR (external mgwr==2.2.1 implementation)
4. Estimated GR-GWR: pilot GWR -> Queen-constrained Ward K=3 -> current
   lambda=0.5 boundary refinement -> unified RegimeAwareGWR refit
5. Oracle GR-GWR: same unified RegimeAwareGWR refit but with true labels

All X and y are globally z-standardized before fitting, matching the Georgia
experiments.  True coefficients are transformed analytically to that same
standardized scale, so coefficient-recovery error has a known ground truth.

This is a diagnostic prototype, not the final Monte Carlo experiment.
"""

from __future__ import annotations

from collections import deque
from pathlib import Path
import json
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.sparse import csr_matrix
from scipy.spatial.distance import cdist
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import adjusted_rand_score

from mgwr.gwr import GWR, MGWR
from mgwr.sel_bw import Sel_BW

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from georegime_gwr import RegimeAwareGWR
from georegime_gwr.gwr import BasicGWR

OUT = ROOT / "results" / "synthetic" / "canonical_regime_boundary"

GRID_N = 25
SEED = 20260901
RHO_X = 0.20
NOISE_SIGMA = 0.75
K_TRUE = 3
LAMBDA_BOUNDARY = 0.5
MAX_REFINEMENT_ITER = 20
MIN_REGIME_N = 6
TERMS = ["Intercept", "X1", "X2", "X3"]
SLOPE_TERMS = TERMS[1:]

# Raw-scale true coefficients: [intercept, X1, X2, X3].
RAW_BETA = {
    1: np.array([0.0, 1.5, -1.0, 0.5]),
    2: np.array([0.0, -1.0, 1.0, 0.5]),
    3: np.array([0.0, 0.5, -0.5, -1.0]),
}


def build_grid() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    axis = np.linspace(0.0, 1.0, GRID_N)
    xx, yy = np.meshgrid(axis, axis)
    coords = np.column_stack([xx.ravel(), yy.ravel()])

    # True regimes: left block + upper-right + lower-right.
    true = np.empty(coords.shape[0], dtype=int)
    x, y = coords[:, 0], coords[:, 1]
    true[x < 0.45] = 1
    true[(x >= 0.45) & (y >= 0.50)] = 2
    true[(x >= 0.45) & (y < 0.50)] = 3

    # Queen graph on the regular lattice.
    edges: list[tuple[int, int]] = []
    for r in range(GRID_N):
        for c in range(GRID_N):
            i = r * GRID_N + c
            for dr, dc in [(0, 1), (1, -1), (1, 0), (1, 1)]:
                rr, cc = r + dr, c + dc
                if 0 <= rr < GRID_N and 0 <= cc < GRID_N:
                    j = rr * GRID_N + cc
                    edges.append((i, j))
    return coords, true, np.asarray(edges, dtype=int)


def adjacency_from_edges(n: int, edges: np.ndarray) -> list[set[int]]:
    adjacency = [set() for _ in range(n)]
    for i, j in edges:
        adjacency[int(i)].add(int(j))
        adjacency[int(j)].add(int(i))
    return adjacency


def connectivity_matrix(n: int, edges: np.ndarray) -> csr_matrix:
    rows = np.concatenate([edges[:, 0], edges[:, 1]])
    cols = np.concatenate([edges[:, 1], edges[:, 0]])
    data = np.ones(rows.size, dtype=float)
    return csr_matrix((data, (rows, cols)), shape=(n, n))


def generate_data(coords: np.ndarray, true: np.ndarray):
    rng = np.random.default_rng(SEED)
    cov = np.full((3, 3), RHO_X, dtype=float)
    np.fill_diagonal(cov, 1.0)
    X_raw = rng.multivariate_normal(np.zeros(3), cov, size=coords.shape[0])

    beta_raw = np.vstack([RAW_BETA[int(r)] for r in true])
    y_signal = beta_raw[:, 0] + np.sum(beta_raw[:, 1:] * X_raw, axis=1)
    y_raw = y_signal + rng.normal(0.0, NOISE_SIGMA, size=coords.shape[0])

    x_mean = X_raw.mean(axis=0)
    x_sd = X_raw.std(axis=0, ddof=0)
    y_mean = float(y_raw.mean())
    y_sd = float(y_raw.std(ddof=0))
    Xz = (X_raw - x_mean) / x_sd
    yz = (y_raw - y_mean) / y_sd

    # Transform the known raw-scale coefficients to the exact standardized
    # coordinate system used for model fitting.
    beta_true_z = np.empty_like(beta_raw)
    beta_true_z[:, 1:] = beta_raw[:, 1:] * x_sd[None, :] / y_sd
    beta_true_z[:, 0] = (
        beta_raw[:, 0] + np.sum(beta_raw[:, 1:] * x_mean[None, :], axis=1) - y_mean
    ) / y_sd

    return Xz, yz, beta_true_z, {
        "x_mean": x_mean.tolist(),
        "x_sd": x_sd.tolist(),
        "y_mean": y_mean,
        "y_sd": y_sd,
        "noise_sigma_raw": NOISE_SIGMA,
    }


def boundary_edges(labels: np.ndarray, edges: np.ndarray) -> np.ndarray:
    return labels[edges[:, 0]] != labels[edges[:, 1]]


def boundary_nodes(labels: np.ndarray, edges: np.ndarray) -> np.ndarray:
    mask = boundary_edges(labels, edges)
    nodes = np.zeros(labels.size, dtype=bool)
    if np.any(mask):
        nodes[np.unique(edges[mask].ravel())] = True
    return nodes


def boundary_scores(true: np.ndarray, pred: np.ndarray, edges: np.ndarray) -> dict[str, float]:
    t = boundary_edges(true, edges)
    p = boundary_edges(pred, edges)
    tp = int(np.sum(t & p))
    fp = int(np.sum(~t & p))
    fn = int(np.sum(t & ~p))
    precision = tp / (tp + fp) if tp + fp else 1.0
    recall = tp / (tp + fn) if tp + fn else 1.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    jaccard = tp / (tp + fp + fn) if tp + fp + fn else 1.0
    return {
        "true_boundary_edges": int(np.sum(t)),
        "estimated_boundary_edges": int(np.sum(p)),
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "precision": float(precision),
        "recall": float(recall),
        "f1": float(f1),
        "jaccard": float(jaccard),
    }


def align_labels(pred: np.ndarray, true: np.ndarray) -> np.ndarray:
    pred_u = np.unique(pred)
    true_u = np.unique(true)
    confusion = np.zeros((pred_u.size, true_u.size), dtype=int)
    for a, pa in enumerate(pred_u):
        for b, tb in enumerate(true_u):
            confusion[a, b] = int(np.sum((pred == pa) & (true == tb)))
    rr, cc = linear_sum_assignment(-confusion)
    mapping = {int(pred_u[r]): int(true_u[c]) for r, c in zip(rr, cc)}
    return np.asarray([mapping[int(v)] for v in pred], dtype=int)


def all_regimes_connected(labels: np.ndarray, adjacency: list[set[int]]) -> bool:
    for r in np.unique(labels):
        nodes = np.flatnonzero(labels == r)
        allowed = set(nodes.tolist())
        if not allowed:
            return False
        seen = {int(nodes[0])}
        q = deque([int(nodes[0])])
        while q:
            u = q.popleft()
            for v in adjacency[u]:
                if v in allowed and v not in seen:
                    seen.add(v)
                    q.append(v)
        if len(seen) != len(allowed):
            return False
    return True


def source_remains_connected(i: int, labels: np.ndarray, adjacency: list[set[int]]) -> bool:
    r = labels[i]
    nodes = np.flatnonzero(labels == r)
    if nodes.size <= 1:
        return False
    remaining = [int(v) for v in nodes if int(v) != i]
    allowed = set(remaining)
    start = remaining[0]
    seen = {start}
    q = deque([start])
    while q:
        u = q.popleft()
        for v in adjacency[u]:
            if v in allowed and v not in seen:
                seen.add(v)
                q.append(v)
    return len(seen) == len(remaining)


def ward_initial_partition(slopes: np.ndarray, edges: np.ndarray) -> np.ndarray:
    z = (slopes - slopes.mean(axis=0)) / slopes.std(axis=0, ddof=0)
    conn = connectivity_matrix(slopes.shape[0], edges)
    cluster = AgglomerativeClustering(
        n_clusters=K_TRUE,
        linkage="ward",
        connectivity=conn,
    )
    return cluster.fit_predict(z).astype(int) + 1


def bisquare_all_training(distances: np.ndarray) -> np.ndarray:
    if distances.size == 0:
        return distances.copy()
    bw = float(np.max(distances))
    if bw <= 1e-12:
        bw = 1.0
    bw = float(np.nextafter(bw, np.inf))
    ratio = distances / bw
    return np.where(ratio < 1.0, (1.0 - ratio**2) ** 2, 0.0)


def loo_cost_for_label(
    i: int,
    target: int,
    labels: np.ndarray,
    Xd: np.ndarray,
    y: np.ndarray,
    distances: np.ndarray,
    adjacency: list[set[int]],
) -> float:
    train = (labels == target)
    train[i] = False
    idx = np.flatnonzero(train)
    if idx.size < Xd.shape[1] + 1:
        return float("inf")
    weights = bisquare_all_training(distances[i, idx])
    beta, _ = BasicGWR._solve_local(Xd[idx], y[idx], weights)
    pred = float(Xd[i] @ beta)
    err2 = float((y[i] - pred) ** 2)
    mismatch = sum(1 for j in adjacency[i] if labels[j] != target)
    return err2 + LAMBDA_BOUNDARY * mismatch


def global_objective(
    X: np.ndarray,
    y: np.ndarray,
    coords: np.ndarray,
    labels: np.ndarray,
    edges: np.ndarray,
):
    model = RegimeAwareGWR(
        bandwidth="regime_size", kernel="bisquare", fit_intercept=True
    ).fit(X, y, coords, labels)
    rss = float(model.residuals_ @ model.residuals_)
    b = int(np.sum(boundary_edges(labels, edges)))
    return rss + LAMBDA_BOUNDARY * b, rss, b, model


def refine_labels(
    initial: np.ndarray,
    X: np.ndarray,
    y: np.ndarray,
    coords: np.ndarray,
    edges: np.ndarray,
    adjacency: list[set[int]],
):
    labels = initial.copy()
    Xd = np.column_stack([np.ones(X.shape[0]), X])
    distances = cdist(coords, coords)
    objective, rss, bcount, model = global_objective(X, y, coords, labels, edges)
    history = [{
        "iteration": 0,
        "changes": 0,
        "rss": rss,
        "boundary_edges": bcount,
        "objective": objective,
        "sizes": ";".join(str(int(np.sum(labels == r))) for r in sorted(np.unique(labels))),
    }]

    stop_reason = "max_iterations"
    for iteration in range(1, MAX_REFINEMENT_ITER + 1):
        before = labels.copy()
        old_objective = objective
        rng = np.random.default_rng(42 + iteration)

        # Only boundary units can have a neighbouring target regime.  Scanning
        # interior units would be mathematically redundant and computationally wasteful.
        candidates = np.flatnonzero(boundary_nodes(labels, edges))
        order = rng.permutation(candidates)
        changes = 0

        for i_raw in order:
            i = int(i_raw)
            current = int(labels[i])
            if int(np.sum(labels == current)) <= MIN_REGIME_N:
                continue
            neighbour_labels = {int(labels[j]) for j in adjacency[i]}
            targets = sorted(neighbour_labels | {current})
            if len(targets) <= 1:
                continue

            current_cost = loo_cost_for_label(i, current, labels, Xd, y, distances, adjacency)
            best = current
            best_cost = current_cost
            for target in targets:
                if target == current:
                    continue
                if not source_remains_connected(i, labels, adjacency):
                    continue
                cost = loo_cost_for_label(i, target, labels, Xd, y, distances, adjacency)
                if cost + 1e-12 < best_cost:
                    best = target
                    best_cost = cost
            if best != current:
                labels[i] = best
                changes += 1

        if changes == 0:
            stop_reason = "no_label_changes"
            break

        new_objective, new_rss, new_bcount, new_model = global_objective(X, y, coords, labels, edges)
        if new_objective > old_objective + 1e-10:
            labels = before
            stop_reason = "global_objective_guard_rejected_sweep"
            break

        objective, rss, bcount, model = new_objective, new_rss, new_bcount, new_model
        history.append({
            "iteration": iteration,
            "changes": changes,
            "rss": rss,
            "boundary_edges": bcount,
            "objective": objective,
            "sizes": ";".join(str(int(np.sum(labels == r))) for r in sorted(np.unique(labels))),
        })

    return labels, model, pd.DataFrame(history), stop_reason


def regression_metrics(y: np.ndarray, fitted: np.ndarray) -> dict[str, float]:
    resid = y - fitted
    rss = float(resid @ resid)
    return {
        "rss": rss,
        "rmse": float(np.sqrt(np.mean(resid**2))),
        "mae": float(np.mean(np.abs(resid))),
        "r2": float(1.0 - rss / np.sum((y - np.mean(y)) ** 2)),
    }


def coefficient_metrics(
    name: str,
    params: np.ndarray,
    beta_true: np.ndarray,
    true_boundary_nodes: np.ndarray,
) -> list[dict[str, float | str | int]]:
    rows = []
    for j, term in enumerate(TERMS):
        err = params[:, j] - beta_true[:, j]
        berr = err[true_boundary_nodes]
        ierr = err[~true_boundary_nodes]
        rows.append({
            "model": name,
            "term": term,
            "rmse_beta": float(np.sqrt(np.mean(err**2))),
            "mae_beta": float(np.mean(np.abs(err))),
            "bias_beta": float(np.mean(err)),
            "boundary_band_rmse_beta": float(np.sqrt(np.mean(berr**2))),
            "interior_rmse_beta": float(np.sqrt(np.mean(ierr**2))),
        })
    return rows


def jump_recovery_rows(
    model: str,
    params: np.ndarray,
    beta_true: np.ndarray,
    true: np.ndarray,
    edges: np.ndarray,
) -> list[dict[str, float | str]]:
    mask = boundary_edges(true, edges)
    ee = edges[mask]
    rows = []
    for j, term in enumerate(SLOPE_TERMS, start=1):
        truth = np.abs(beta_true[ee[:, 0], j] - beta_true[ee[:, 1], j])
        est = np.abs(params[ee[:, 0], j] - params[ee[:, 1], j])
        valid = truth > 1e-12
        ratio = float(np.mean(est[valid]) / np.mean(truth[valid])) if np.any(valid) else float("nan")
        rows.append({
            "model": model,
            "term": term,
            "true_mean_jump": float(np.mean(truth)),
            "estimated_mean_jump": float(np.mean(est)),
            "jump_recovery_ratio": ratio,
        })
    return rows


def plot_regime_recovery(true, initial, refined, boundary_df: pd.DataFrame):
    aligned_initial = align_labels(initial, true)
    aligned_refined = align_labels(refined, true)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.1))
    for ax, labels, title in zip(
        axes,
        [true, aligned_initial, aligned_refined],
        ["True regimes", "Initial Queen-Ward K=3", "Refined GR-GWR regimes"],
    ):
        ax.imshow(labels.reshape(GRID_N, GRID_N), origin="lower", cmap="tab10", vmin=1, vmax=10)
        ax.set_title(title)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle(
        "Canonical boundary recovery | "
        f"initial ARI={boundary_df.loc[0, 'ari']:.3f}, refined ARI={boundary_df.loc[1, 'ari']:.3f}",
        fontsize=14,
    )
    fig.tight_layout()
    fig.savefig(OUT / "regime_recovery.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "regime_recovery.svg", bbox_inches="tight")
    plt.close(fig)


def plot_coefficients(model_params: dict[str, np.ndarray], beta_true: np.ndarray):
    ordered = ["True", "GWR", "MGWR", "GR-GWR", "Oracle GR-GWR"]
    for j, term in enumerate(SLOPE_TERMS, start=1):
        arrays = [beta_true[:, j]] + [model_params[name][:, j] for name in ordered[1:]]
        vmax = max(float(np.max(np.abs(a))) for a in arrays)
        fig, axes = plt.subplots(1, 5, figsize=(20, 4.2))
        image = None
        for ax, name, values in zip(axes, ordered, arrays):
            image = ax.imshow(
                values.reshape(GRID_N, GRID_N),
                origin="lower",
                cmap="coolwarm",
                vmin=-vmax,
                vmax=vmax,
            )
            ax.set_title(name)
            ax.set_xticks([])
            ax.set_yticks([])
        fig.colorbar(image, ax=axes.tolist(), shrink=0.78, pad=0.02)
        fig.suptitle(f"{term}: true and recovered standardized local coefficients", fontsize=14)
        fig.savefig(OUT / f"{term}_coefficient_recovery.png", dpi=300, bbox_inches="tight")
        fig.savefig(OUT / f"{term}_coefficient_recovery.svg", bbox_inches="tight")
        plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    coords, true, edges = build_grid()
    adjacency = adjacency_from_edges(coords.shape[0], edges)
    X, y, beta_true, scaling = generate_data(coords, true)
    y_col = y.reshape(-1, 1)

    # OLS.
    Xd = np.column_stack([np.ones(X.shape[0]), X])
    beta_ols, *_ = np.linalg.lstsq(Xd, y, rcond=None)
    ols_params = np.tile(beta_ols, (X.shape[0], 1))
    ols_fitted = Xd @ beta_ols

    # Standard GWR via the widely used external mgwr implementation.
    gwr_selector = Sel_BW(coords, y_col, X, fixed=False, kernel="bisquare", constant=True)
    gwr_bw = int(round(float(gwr_selector.search())))
    gwr_res = GWR(
        coords, y_col, X, gwr_bw,
        fixed=False, kernel="bisquare", constant=True,
    ).fit()
    gwr_params = np.asarray(gwr_res.params, dtype=float)
    gwr_fitted = np.asarray(gwr_res.predy, dtype=float).reshape(-1)

    # Queen-constrained Ward on pilot GWR slopes; no coordinates or intercept
    # are included in the relationship fingerprint.
    initial = ward_initial_partition(gwr_params[:, 1:], edges)
    if not all_regimes_connected(initial, adjacency):
        raise RuntimeError("Initial Ward partition is unexpectedly disconnected")

    refined, gr_model, history, stop_reason = refine_labels(
        initial, X, y, coords, edges, adjacency
    )
    if not all_regimes_connected(refined, adjacency):
        raise RuntimeError("Refined partition is disconnected")

    # Oracle uses the exact same refit primitive but receives the true labels.
    oracle = RegimeAwareGWR(
        bandwidth="regime_size", kernel="bisquare", fit_intercept=True
    ).fit(X, y, coords, true)

    # External MGWR comparator.
    mgwr_selector = Sel_BW(
        coords, y_col, X,
        multi=True, fixed=False, kernel="bisquare", constant=True,
    )
    mgwr_bws = np.asarray(
        mgwr_selector.search(
            multi_bw_min=[2],
            tol_multi=1e-5,
            max_iter_multi=100,
            verbose=False,
        ),
        dtype=float,
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

    metrics_rows = []
    for name in ["OLS", "GWR", "MGWR", "GR-GWR", "Oracle GR-GWR"]:
        metrics_rows.append({"model": name, **regression_metrics(y, model_fitted[name])})
    metrics_df = pd.DataFrame(metrics_rows)
    metrics_df.to_csv(OUT / "model_fit_metrics.csv", index=False)

    true_bnodes = boundary_nodes(true, edges)
    coef_rows = []
    jump_rows = []
    for name, params in model_params.items():
        coef_rows.extend(coefficient_metrics(name, params, beta_true, true_bnodes))
        jump_rows.extend(jump_recovery_rows(name, params, beta_true, true, edges))
    coef_df = pd.DataFrame(coef_rows)
    jump_df = pd.DataFrame(jump_rows)
    coef_df.to_csv(OUT / "coefficient_recovery.csv", index=False)
    jump_df.to_csv(OUT / "true_boundary_jump_recovery.csv", index=False)

    boundary_rows = []
    for stage, labels in [("initial_ward", initial), ("refined_grgwr", refined)]:
        boundary_rows.append({
            "stage": stage,
            "ari": float(adjusted_rand_score(true, labels)),
            **boundary_scores(true, labels, edges),
            "regime_sizes": ";".join(str(int(np.sum(labels == r))) for r in sorted(np.unique(labels))),
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
        grid_df[f"true_{term}"] = beta_true[:, j]
        for name, params in model_params.items():
            safe = name.lower().replace(" ", "_").replace("-", "")
            grid_df[f"{safe}_{term}"] = params[:, j]
    grid_df.to_csv(OUT / "grid_results.csv", index=False)

    plot_regime_recovery(true, initial, refined, boundary_df)
    plot_coefficients(model_params, beta_true)

    summary = {
        "status": "single_run_canonical_synthetic_prototype",
        "seed": SEED,
        "grid": f"{GRID_N}x{GRID_N}",
        "n": int(coords.shape[0]),
        "K_true": K_TRUE,
        "K_given_but_boundaries_hidden": True,
        "rho_X": RHO_X,
        "noise_sigma_raw": NOISE_SIGMA,
        "lambda_working_baseline": LAMBDA_BOUNDARY,
        "minimum_regime_n_guard": MIN_REGIME_N,
        "gwr_bandwidth": gwr_bw,
        "mgwr_bandwidths": {term: int(round(bw)) for term, bw in zip(TERMS, mgwr_bws)},
        "true_regime_sizes": {str(r): int(np.sum(true == r)) for r in np.unique(true)},
        "initial_ari": float(boundary_df.loc[boundary_df.stage == "initial_ward", "ari"].iloc[0]),
        "refined_ari": float(boundary_df.loc[boundary_df.stage == "refined_grgwr", "ari"].iloc[0]),
        "initial_boundary_f1": float(boundary_df.loc[boundary_df.stage == "initial_ward", "f1"].iloc[0]),
        "refined_boundary_f1": float(boundary_df.loc[boundary_df.stage == "refined_grgwr", "f1"].iloc[0]),
        "refinement_stop_reason": stop_reason,
        "refinement_accepted_iterations": int(history.iteration.max()),
        "scaling": scaling,
        "model_fit": metrics_df.to_dict(orient="records"),
        "interpretation_guard": (
            "Single realization only. K=3 is fixed and known to the experimenter, but true boundary locations are hidden from estimated GR-GWR. "
            "This prototype is for behavioral validation before Monte Carlo, noise-strength, K-selection, and out-of-sample experiments."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(json.dumps(summary, indent=2))
    print("\nBoundary recovery:\n", boundary_df.to_string(index=False))
    print("\nModel fit:\n", metrics_df.to_string(index=False))
    print("\nCoefficient recovery:\n", coef_df.to_string(index=False))
    print("\nTrue-boundary jump recovery:\n", jump_df.to_string(index=False))


if __name__ == "__main__":
    main()
