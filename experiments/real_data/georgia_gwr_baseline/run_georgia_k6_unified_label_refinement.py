"""Run the current K=6 GR-GWR baseline through label refinement to convergence.

Current exploratory chain:

    pilot GWR -> Queen-constrained Ward K=6 -> unified RegimeAwareGWR
    -> Queen-constrained ICM-style label refinement -> unified refit -> repeat

This deliberately follows the existing baseline logic before unresolved method
questions (lambda choice, complexity control, final bandwidth policy, etc.) are
redesigned. The current lambda=1.0 is exploratory only.
"""

from __future__ import annotations

from pathlib import Path
import json
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from georegime_gwr import BasicGWR, RegimeAwareGWR

DATA = ROOT / "data" / "raw" / "georgia" / "GData_utm.csv"
LABELS = ROOT / "results" / "spatial" / "georgia_initial_ward_regimes" / "initial_regime_labels_k2_k15.csv"
EDGES = ROOT / "results" / "spatial" / "georgia_queen_w" / "edge_list.csv"
OUT = ROOT / "results" / "real_data" / "georgia_k6_unified_label_refinement"
PREDICTORS = ["PctFB", "PctBlack", "PctRural"]
TERMS = ["Intercept", *PREDICTORS]
LAMBDA_BOUNDARY = 1.0
MIN_REGIME_SIZE = 6  # p=4, preserving the old p+2 baseline guard.
MAX_ITER = 20
TOL = 1e-10
RANDOM_STATE = 42


def _key(value) -> str:
    try:
        return str(int(float(value)))
    except (TypeError, ValueError):
        return str(value)


def _build_adjacency(n: int, edge_df: pd.DataFrame):
    adjacency = [set() for _ in range(n)]
    edges = []
    for row in edge_df.itertuples(index=False):
        i, j = int(row.i), int(row.j)
        adjacency[i].add(j)
        adjacency[j].add(i)
        edges.append((i, j))
    return tuple(np.asarray(sorted(v), dtype=int) for v in adjacency), tuple(edges)


def _boundary_count(labels: np.ndarray, edges) -> int:
    return int(sum(labels[i] != labels[j] for i, j in edges))


def _all_connected(labels: np.ndarray, adjacency) -> bool:
    for regime in np.unique(labels):
        members = np.flatnonzero(labels == regime)
        if members.size <= 1:
            continue
        allowed = set(map(int, members))
        start = int(members[0])
        visited = {start}
        stack = [start]
        while stack:
            current = stack.pop()
            for neighbour in adjacency[current]:
                j = int(neighbour)
                if j in allowed and j not in visited:
                    visited.add(j)
                    stack.append(j)
        if len(visited) != members.size:
            return False
    return True


def _can_remove(node: int, labels: np.ndarray, adjacency) -> bool:
    source = labels[node]
    members = np.flatnonzero(labels == source)
    if members.size - 1 < MIN_REGIME_SIZE:
        return False

    remaining = members[members != node]
    if remaining.size <= 1:
        return True
    allowed = set(map(int, remaining))
    start = int(remaining[0])
    visited = {start}
    stack = [start]
    while stack:
        current = stack.pop()
        for neighbour in adjacency[current]:
            j = int(neighbour)
            if j in allowed and j not in visited:
                visited.add(j)
                stack.append(j)
    return len(visited) == remaining.size


def _candidate_error(
    node: int,
    regime: int,
    labels: np.ndarray,
    Xd: np.ndarray,
    y: np.ndarray,
    distances: np.ndarray,
) -> float:
    """Leave-one-out prediction cost for assigning node to candidate regime."""
    indices = np.flatnonzero(labels == regime)
    indices = indices[indices != node]
    p = Xd.shape[1]
    if indices.size < MIN_REGIME_SIZE - 1 or indices.size < p:
        return float("inf")
    if np.linalg.matrix_rank(Xd[indices]) < p:
        return float("inf")

    d = distances[node, indices]
    k = indices.size  # current baseline: regime-size adaptive neighborhood.
    bw = float(np.partition(d, k - 1)[k - 1])
    if bw <= 1e-12:
        positive = d[d > 1e-12]
        bw = float(np.min(positive)) if positive.size else 1.0
    bw = float(np.nextafter(bw, np.inf))
    ratio = d / bw
    weights = np.where(ratio < 1.0, (1.0 - ratio**2) ** 2, 0.0)

    beta, _ = BasicGWR._solve_local(Xd[indices], y[indices], weights)
    pred = float(Xd[node] @ beta)
    return float((y[node] - pred) ** 2)


def _metrics(y: np.ndarray, model: RegimeAwareGWR) -> dict[str, float]:
    residuals = model.residuals_
    rss = float(residuals @ residuals)
    return {
        "rss": rss,
        "rmse": float(np.sqrt(np.mean(residuals**2))),
        "mae": float(np.mean(np.abs(residuals))),
        "r2": float(1.0 - rss / np.sum((y - np.mean(y)) ** 2)),
        "trace_hat": float(np.trace(model.hat_matrix_)),
    }


def _regime_sizes(labels: np.ndarray) -> dict[int, int]:
    return {int(r): int(np.sum(labels == r)) for r in np.unique(labels)}


def _fit(X: np.ndarray, y: np.ndarray, coords: np.ndarray, labels: np.ndarray):
    return RegimeAwareGWR(
        bandwidth="regime_size",
        kernel="bisquare",
        fit_intercept=True,
    ).fit(X, y, coords, labels)


def _sweep(
    labels: np.ndarray,
    iteration: int,
    adjacency,
    Xd: np.ndarray,
    y: np.ndarray,
    distances: np.ndarray,
):
    updated = labels.copy()
    order = np.random.default_rng(RANDOM_STATE + iteration).permutation(labels.size)
    proposed_moves = []

    for node in map(int, order):
        current = int(updated[node])
        if not _can_remove(node, updated, adjacency):
            continue

        neighbour_labels = {int(updated[j]) for j in adjacency[node]}
        candidates = sorted(neighbour_labels | {current})
        costs = {}
        errors = {}
        disagreements = {}
        for regime in candidates:
            error = _candidate_error(node, regime, updated, Xd, y, distances)
            disagreement = int(np.sum(updated[adjacency[node]] != regime))
            cost = error + LAMBDA_BOUNDARY * disagreement
            costs[regime] = cost
            errors[regime] = error
            disagreements[regime] = disagreement

        best = min(costs, key=costs.get)
        if best != current and costs[best] < costs[current] - TOL:
            proposed_moves.append({
                "iteration": iteration,
                "node": node,
                "from_regime": current,
                "to_regime": int(best),
                "current_cost": float(costs[current]),
                "best_cost": float(costs[best]),
                "current_loo_error": float(errors[current]),
                "best_loo_error": float(errors[best]),
                "current_boundary_disagreement": int(disagreements[current]),
                "best_boundary_disagreement": int(disagreements[best]),
            })
            updated[node] = best

    return updated, proposed_moves


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(DATA)
    label_df = pd.read_csv(LABELS)
    edge_df = pd.read_csv(EDGES)
    df["AreaKey_join"] = df["AreaKey"].map(_key)
    label_df["AreaKey_join"] = label_df["AreaKey"].map(_key)
    df = df.merge(label_df[["AreaKey_join", "K6"]], on="AreaKey_join", how="left", validate="1:1")
    if len(df) != 159 or df["K6"].isna().any():
        raise RuntimeError("Expected complete Georgia K=6 inputs")

    X = df[PREDICTORS].to_numpy(dtype=float)
    y = df["PctBach"].to_numpy(dtype=float)
    coords = df[["X", "Y"]].to_numpy(dtype=float)
    Xz = (X - X.mean(axis=0)) / X.std(axis=0, ddof=0)
    yz = (y - y.mean()) / y.std(ddof=0)
    Xd = np.column_stack([np.ones(len(df)), Xz])

    labels = df["K6"].to_numpy(dtype=int)
    initial_labels = labels.copy()
    adjacency, edges = _build_adjacency(len(df), edge_df)
    if not _all_connected(labels, adjacency):
        raise RuntimeError("Initial K=6 partition is not Queen-connected")

    model = _fit(Xz, yz, coords, labels)
    distances = model.distance_matrix_
    metrics = _metrics(yz, model)
    boundary = _boundary_count(labels, edges)
    objective = metrics["rss"] + LAMBDA_BOUNDARY * boundary

    history = [{
        "iteration": 0,
        "accepted": True,
        "label_changes_from_previous": 0,
        "rss": metrics["rss"],
        "rmse": metrics["rmse"],
        "mae": metrics["mae"],
        "r2": metrics["r2"],
        "trace_hat": metrics["trace_hat"],
        "boundary_edges": boundary,
        "objective": objective,
        "regime_sizes": ";".join(str(v) for _, v in sorted(_regime_sizes(labels).items())),
    }]
    label_snapshots = {"initial": labels.copy()}
    accepted_move_rows = []
    stop_reason = "max_iter"

    for iteration in range(1, MAX_ITER + 1):
        proposed_labels, proposed_moves = _sweep(
            labels, iteration, adjacency, Xd, yz, distances
        )
        changed = int(np.sum(proposed_labels != labels))
        if changed == 0:
            stop_reason = "no_label_changes"
            break

        if not _all_connected(proposed_labels, adjacency):
            raise RuntimeError(f"Connectivity guard failed after sweep {iteration}")
        sizes = _regime_sizes(proposed_labels)
        if min(sizes.values()) < MIN_REGIME_SIZE:
            raise RuntimeError(f"Minimum regime-size guard failed after sweep {iteration}")

        proposed_model = _fit(Xz, yz, coords, proposed_labels)
        proposed_metrics = _metrics(yz, proposed_model)
        proposed_boundary = _boundary_count(proposed_labels, edges)
        proposed_objective = proposed_metrics["rss"] + LAMBDA_BOUNDARY * proposed_boundary

        accepted = proposed_objective <= objective + TOL
        history.append({
            "iteration": iteration,
            "accepted": accepted,
            "label_changes_from_previous": changed,
            "rss": proposed_metrics["rss"],
            "rmse": proposed_metrics["rmse"],
            "mae": proposed_metrics["mae"],
            "r2": proposed_metrics["r2"],
            "trace_hat": proposed_metrics["trace_hat"],
            "boundary_edges": proposed_boundary,
            "objective": proposed_objective,
            "regime_sizes": ";".join(str(v) for _, v in sorted(sizes.items())),
        })

        if not accepted:
            stop_reason = "global_objective_guard_rejected_sweep"
            break

        previous = labels.copy()
        labels = proposed_labels
        model = proposed_model
        metrics = proposed_metrics
        boundary = proposed_boundary
        objective = proposed_objective
        label_snapshots[f"iter{iteration}"] = labels.copy()

        changed_nodes = set(np.flatnonzero(labels != previous).tolist())
        for row in proposed_moves:
            if row["node"] in changed_nodes:
                row = row.copy()
                row["AreaKey"] = df.iloc[row["node"]]["AreaKey"]
                accepted_move_rows.append(row)

    final_labels = labels.copy()
    final_model = model
    final_metrics = _metrics(yz, final_model)
    final_boundary = _boundary_count(final_labels, edges)
    final_objective = final_metrics["rss"] + LAMBDA_BOUNDARY * final_boundary

    pd.DataFrame(history).to_csv(OUT / "iteration_history.csv", index=False)
    pd.DataFrame(accepted_move_rows).to_csv(OUT / "accepted_label_moves.csv", index=False)

    labels_out = pd.DataFrame({"AreaKey": df["AreaKey"]})
    for name, values in label_snapshots.items():
        labels_out[name] = values
    labels_out["final"] = final_labels
    labels_out.to_csv(OUT / "label_path.csv", index=False)

    county = pd.DataFrame({
        "AreaKey": df["AreaKey"],
        "initial_regime": initial_labels,
        "final_regime": final_labels,
        "regime_changed": initial_labels != final_labels,
        "fitted_z": final_model.fitted_values_,
        "residual_z": final_model.residuals_,
        "regime_local_k": final_model.local_bandwidths_,
    })
    for j, term in enumerate(TERMS):
        county[term] = final_model.parameters_[:, j]
    county.to_csv(OUT / "final_county_results.csv", index=False)

    summary = {
        "status": "exploratory_existing_algorithm_unified_refinement",
        "K": 6,
        "n_counties": int(len(df)),
        "lambda_boundary": LAMBDA_BOUNDARY,
        "minimum_regime_size": MIN_REGIME_SIZE,
        "max_iter": MAX_ITER,
        "random_state": RANDOM_STATE,
        "refit_model": "RegimeAwareGWR(bandwidth='regime_size', kernel='bisquare')",
        "candidate_labels": "current label plus Queen-neighbor labels",
        "candidate_prediction_cost": "leave-one-out squared error under candidate regime",
        "connectivity_guard": True,
        "global_objective": "RSS + lambda * Queen boundary-edge count",
        "initial": {
            "regime_sizes": _regime_sizes(initial_labels),
            "boundary_edges": int(history[0]["boundary_edges"]),
            "rss": float(history[0]["rss"]),
            "rmse": float(history[0]["rmse"]),
            "mae": float(history[0]["mae"]),
            "r2": float(history[0]["r2"]),
            "trace_hat": float(history[0]["trace_hat"]),
            "objective": float(history[0]["objective"]),
        },
        "final": {
            "regime_sizes": _regime_sizes(final_labels),
            "boundary_edges": int(final_boundary),
            **{k: float(v) for k, v in final_metrics.items()},
            "objective": float(final_objective),
            "counties_changed_from_initial": int(np.sum(initial_labels != final_labels)),
            "all_regimes_connected": bool(_all_connected(final_labels, adjacency)),
        },
        "accepted_refinement_iterations": int(len(label_snapshots) - 1),
        "stop_reason": stop_reason,
        "warning": (
            "This deliberately follows the existing exploratory ICM/objective baseline. "
            "lambda=1.0, complexity control, bandwidth policy, and the refinement objective are not frozen."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(pd.DataFrame(history).to_string(index=False))
    print("\n", json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
