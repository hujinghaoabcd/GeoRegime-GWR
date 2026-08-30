"""Sensitivity sweep for the current exploratory K=6 refinement boundary penalty.

This does not change the GR-GWR algorithm. It reruns the existing sequential
Queen-constrained refinement from the same initial K=6 Ward partition for a
small grid of lambda values and records when R3 begins to collapse, how many
boundaries are removed, and how fit/complexity respond.
"""

from __future__ import annotations

from pathlib import Path
import importlib.util
import json

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
BASE_SCRIPT = Path(__file__).with_name("run_georgia_k6_unified_label_refinement.py")
OUT = ROOT / "results" / "real_data" / "georgia_k6_lambda_sensitivity"

# Dense near zero because current standardized LOO error differences are often
# much smaller than one boundary-edge unit under lambda=1.
LAMBDAS = [0.0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]


def _load_base_module():
    spec = importlib.util.spec_from_file_location("k6_refinement_base", BASE_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load base refinement script: {BASE_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _run_one(base, lam, Xz, yz, coords, Xd, initial_labels, adjacency, edges, distances):
    base.LAMBDA_BOUNDARY = float(lam)
    labels = initial_labels.copy()
    model = base._fit(Xz, yz, coords, labels)
    metrics = base._metrics(yz, model)
    boundary = base._boundary_count(labels, edges)
    objective = metrics["rss"] + lam * boundary
    stop_reason = "max_iter"
    accepted_iterations = 0

    for iteration in range(1, base.MAX_ITER + 1):
        proposed_labels, _ = base._sweep(labels, iteration, adjacency, Xd, yz, distances)
        changed = int(np.sum(proposed_labels != labels))
        if changed == 0:
            stop_reason = "no_label_changes"
            break

        if not base._all_connected(proposed_labels, adjacency):
            raise RuntimeError(f"Connectivity failed for lambda={lam} iter={iteration}")
        sizes = base._regime_sizes(proposed_labels)
        if min(sizes.values()) < base.MIN_REGIME_SIZE:
            raise RuntimeError(f"Minimum regime size failed for lambda={lam} iter={iteration}")

        proposed_model = base._fit(Xz, yz, coords, proposed_labels)
        proposed_metrics = base._metrics(yz, proposed_model)
        proposed_boundary = base._boundary_count(proposed_labels, edges)
        proposed_objective = proposed_metrics["rss"] + lam * proposed_boundary

        if proposed_objective > objective + base.TOL:
            stop_reason = "global_objective_guard_rejected_sweep"
            break

        labels = proposed_labels
        model = proposed_model
        metrics = proposed_metrics
        boundary = proposed_boundary
        objective = proposed_objective
        accepted_iterations += 1

    sizes = base._regime_sizes(labels)
    row = {
        "lambda": float(lam),
        "accepted_iterations": accepted_iterations,
        "stop_reason": stop_reason,
        "changed_counties_from_initial": int(np.sum(labels != initial_labels)),
        "r3_to_r4_from_initial": int(np.sum((initial_labels == 3) & (labels == 4))),
        "boundary_edges": int(boundary),
        "rss": float(metrics["rss"]),
        "rmse": float(metrics["rmse"]),
        "mae": float(metrics["mae"]),
        "r2": float(metrics["r2"]),
        "trace_hat": float(metrics["trace_hat"]),
        "objective": float(objective),
    }
    for regime in range(1, 7):
        row[f"R{regime}_n"] = int(sizes.get(regime, 0))
    return row, labels


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    base = _load_base_module()

    df = pd.read_csv(base.DATA)
    label_df = pd.read_csv(base.LABELS)
    edge_df = pd.read_csv(base.EDGES)
    df["AreaKey_join"] = df["AreaKey"].map(base._key)
    label_df["AreaKey_join"] = label_df["AreaKey"].map(base._key)
    df = df.merge(label_df[["AreaKey_join", "K6"]], on="AreaKey_join", how="left", validate="1:1")
    if len(df) != 159 or df["K6"].isna().any():
        raise RuntimeError("Expected complete Georgia K=6 inputs")

    X = df[base.PREDICTORS].to_numpy(dtype=float)
    y = df["PctBach"].to_numpy(dtype=float)
    coords = df[["X", "Y"]].to_numpy(dtype=float)
    Xz = (X - X.mean(axis=0)) / X.std(axis=0, ddof=0)
    yz = (y - y.mean()) / y.std(ddof=0)
    Xd = np.column_stack([np.ones(len(df)), Xz])
    initial_labels = df["K6"].to_numpy(dtype=int)
    adjacency, edges = base._build_adjacency(len(df), edge_df)
    if not base._all_connected(initial_labels, adjacency):
        raise RuntimeError("Initial K=6 partition is not Queen-connected")

    initial_model = base._fit(Xz, yz, coords, initial_labels)
    distances = initial_model.distance_matrix_

    rows = []
    label_table = pd.DataFrame({"AreaKey": df["AreaKey"], "initial": initial_labels})
    for lam in LAMBDAS:
        row, labels = _run_one(
            base, lam, Xz, yz, coords, Xd, initial_labels, adjacency, edges, distances
        )
        rows.append(row)
        label_table[f"lambda_{lam:g}"] = labels

    results = pd.DataFrame(rows)
    results.to_csv(OUT / "lambda_sensitivity.csv", index=False)
    label_table.to_csv(OUT / "lambda_final_labels.csv", index=False)

    initial_metrics = base._metrics(yz, initial_model)
    initial_boundary = base._boundary_count(initial_labels, edges)
    payload = {
        "status": "exploratory_lambda_sensitivity_only",
        "lambda_grid": LAMBDAS,
        "initial_R3_n": int(np.sum(initial_labels == 3)),
        "initial_boundary_edges": int(initial_boundary),
        "initial_rmse": float(initial_metrics["rmse"]),
        "minimum_regime_size": int(base.MIN_REGIME_SIZE),
        "note": (
            "Same initial K=6 partition and same refinement algorithm for every lambda. "
            "This diagnoses penalty sensitivity; it does not select a final lambda."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print(results.to_string(index=False))
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
