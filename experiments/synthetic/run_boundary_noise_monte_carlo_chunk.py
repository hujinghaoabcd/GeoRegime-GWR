"""Paired Monte Carlo worker for representative boundary-strength/noise scenarios.

Each invocation runs a small contiguous chunk of seeds for one (delta, noise)
scenario on the same 25x25 lattice used by Simulation 3.  For every seed the
four methods see exactly the same generated X/y realization:

- standard GWR,
- MGWR,
- estimated GR-GWR (K=3 fixed; boundaries hidden),
- oracle GR-GWR (true labels supplied).

The workflow parallelizes chunks; an independent aggregation script combines
all per-seed rows and computes paired confidence intervals and win rates.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys
import time
import warnings

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score
from mgwr.gwr import GWR, MGWR
from mgwr.sel_bw import Sel_BW

ROOT = Path(__file__).resolve().parents[2]
SYN = ROOT / "experiments" / "synthetic"
SRC = ROOT / "src"
for p in (SYN, SRC):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

import run_boundary_noise_phase_scenario as sim
from georegime_gwr import RegimeAwareGWR


def run_one(delta: float, noise: float, seed: int, coords, true, edges, adj) -> dict:
    X, y, beta_true = sim.generate(delta, noise, seed, coords, true)
    yc = y[:, None]
    identifiable = delta > 1e-12

    t0 = time.perf_counter()
    selector = Sel_BW(coords, yc, X, fixed=False, kernel="bisquare", constant=True)
    gwr_bw = int(round(float(selector.search())))
    gwr = GWR(
        coords, yc, X, gwr_bw,
        fixed=False, kernel="bisquare", constant=True,
    ).fit()
    gwr_runtime = time.perf_counter() - t0
    gwr_params = np.asarray(gwr.params, dtype=float)
    gwr_fit = np.asarray(gwr.predy, dtype=float).reshape(-1)

    t0 = time.perf_counter()
    initial = sim.ward(gwr_params[:, 1:], edges)
    refined, gr, niter, stop = sim.refine(initial, X, y, coords, edges, adj)
    gr_after_pilot_runtime = time.perf_counter() - t0

    t0 = time.perf_counter()
    oracle = RegimeAwareGWR(
        bandwidth="regime_size", kernel="bisquare"
    ).fit(X, y, coords, true)
    oracle_runtime = time.perf_counter() - t0

    mgwr_error = None
    mgwr_runtime = float("nan")
    mgwr_bws = np.full(4, np.nan)
    mgwr_params = np.full_like(gwr_params, np.nan)
    mgwr_fit = np.full_like(y, np.nan)
    t0 = time.perf_counter()
    try:
        mselector = Sel_BW(
            coords, yc, X,
            multi=True, fixed=False, kernel="bisquare", constant=True,
        )
        mgwr_bws = np.asarray(
            mselector.search(
                multi_bw_min=[2], tol_multi=1e-4,
                max_iter_multi=60, verbose=False,
            ),
            dtype=float,
        ).reshape(-1)
        mgwr = MGWR(
            coords, yc, X, mselector,
            fixed=False, kernel="bisquare", constant=True,
        ).fit()
        mgwr_params = np.asarray(mgwr.params, dtype=float)
        mgwr_fit = np.asarray(mgwr.predy, dtype=float).reshape(-1)
        mgwr_runtime = time.perf_counter() - t0
    except Exception as exc:
        mgwr_runtime = time.perf_counter() - t0
        mgwr_error = f"{type(exc).__name__}: {exc}"

    gwr_fit_rmse = sim.fit_rmse(y, gwr_fit)
    gr_fit_rmse = sim.fit_rmse(y, gr.fitted_values_)
    oracle_fit_rmse = sim.fit_rmse(y, oracle.fitted_values_)
    mgwr_ok = np.all(np.isfinite(mgwr_fit)) and np.all(np.isfinite(mgwr_params))
    mgwr_fit_rmse = sim.fit_rmse(y, mgwr_fit) if mgwr_ok else float("nan")

    gwr_coef_rmse = sim.coef_rmse(gwr_params, beta_true)
    gr_coef_rmse = sim.coef_rmse(gr.parameters_, beta_true)
    oracle_coef_rmse = sim.coef_rmse(oracle.parameters_, beta_true)
    mgwr_coef_rmse = sim.coef_rmse(mgwr_params, beta_true) if mgwr_ok else float("nan")

    initial_ari = float(adjusted_rand_score(true, initial)) if identifiable else float("nan")
    refined_ari = float(adjusted_rand_score(true, refined)) if identifiable else float("nan")
    initial_f1 = sim.boundary_f1(true, initial, edges) if identifiable else float("nan")
    refined_f1 = sim.boundary_f1(true, refined, edges) if identifiable else float("nan")

    gwr_jump = sim.jump_ratio(gwr_params, beta_true, true, edges) if identifiable else float("nan")
    gr_jump = sim.jump_ratio(gr.parameters_, beta_true, true, edges) if identifiable else float("nan")
    oracle_jump = sim.jump_ratio(oracle.parameters_, beta_true, true, edges) if identifiable else float("nan")
    mgwr_jump = sim.jump_ratio(mgwr_params, beta_true, true, edges) if identifiable and mgwr_ok else float("nan")

    return {
        "delta": float(delta),
        "noise_sigma": float(noise),
        "seed": int(seed),
        "gwr_bandwidth": int(gwr_bw),
        "mgwr_bw_intercept": float(mgwr_bws[0]),
        "mgwr_bw_x1": float(mgwr_bws[1]),
        "mgwr_bw_x2": float(mgwr_bws[2]),
        "mgwr_bw_x3": float(mgwr_bws[3]),
        "mgwr_error": mgwr_error,
        "initial_ari": initial_ari,
        "refined_ari": refined_ari,
        "initial_boundary_f1": initial_f1,
        "refined_boundary_f1": refined_f1,
        "initial_boundary_edges": int(np.sum(sim.boundary_mask(initial, edges))),
        "refined_boundary_edges": int(np.sum(sim.boundary_mask(refined, edges))),
        "refinement_iterations": int(niter),
        "refinement_stop": stop,
        "gwr_fit_rmse": gwr_fit_rmse,
        "mgwr_fit_rmse": mgwr_fit_rmse,
        "grgwr_fit_rmse": gr_fit_rmse,
        "oracle_fit_rmse": oracle_fit_rmse,
        "gwr_coef_rmse": gwr_coef_rmse,
        "mgwr_coef_rmse": mgwr_coef_rmse,
        "grgwr_coef_rmse": gr_coef_rmse,
        "oracle_coef_rmse": oracle_coef_rmse,
        "gwr_jump_ratio": gwr_jump,
        "mgwr_jump_ratio": mgwr_jump,
        "grgwr_jump_ratio": gr_jump,
        "oracle_jump_ratio": oracle_jump,
        "gr_minus_gwr_fit_rmse": gr_fit_rmse - gwr_fit_rmse,
        "gr_minus_mgwr_fit_rmse": gr_fit_rmse - mgwr_fit_rmse if mgwr_ok else float("nan"),
        "gr_minus_gwr_coef_rmse": gr_coef_rmse - gwr_coef_rmse,
        "gr_minus_mgwr_coef_rmse": gr_coef_rmse - mgwr_coef_rmse if mgwr_ok else float("nan"),
        "gr_fit_better_than_gwr": int(gr_fit_rmse < gwr_fit_rmse),
        "gr_fit_better_than_mgwr": int(gr_fit_rmse < mgwr_fit_rmse) if mgwr_ok else np.nan,
        "gr_coef_better_than_gwr": int(gr_coef_rmse < gwr_coef_rmse),
        "gr_coef_better_than_mgwr": int(gr_coef_rmse < mgwr_coef_rmse) if mgwr_ok else np.nan,
        "runtime_gwr_seconds": float(gwr_runtime),
        "runtime_gr_after_pilot_seconds": float(gr_after_pilot_runtime),
        "runtime_gr_total_seconds": float(gwr_runtime + gr_after_pilot_runtime),
        "runtime_mgwr_seconds": float(mgwr_runtime),
        "runtime_oracle_seconds": float(oracle_runtime),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scenario", required=True)
    ap.add_argument("--delta", type=float, required=True)
    ap.add_argument("--noise", type=float, required=True)
    ap.add_argument("--seed-start", type=int, required=True)
    ap.add_argument("--reps", type=int, default=10)
    ap.add_argument("--chunk", type=int, required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    coords, true, edges = sim.build_grid()
    adj = sim.adjacency_from_edges(len(true), edges)

    rows = []
    for offset in range(args.reps):
        seed = args.seed_start + offset
        print(
            f"scenario={args.scenario} chunk={args.chunk} "
            f"seed={seed} delta={args.delta} noise={args.noise}",
            flush=True,
        )
        try:
            rows.append(run_one(args.delta, args.noise, seed, coords, true, edges, adj))
        except Exception as exc:
            warnings.warn(f"seed {seed} failed: {type(exc).__name__}: {exc}")
            rows.append({
                "delta": args.delta,
                "noise_sigma": args.noise,
                "seed": seed,
                "fatal_error": f"{type(exc).__name__}: {exc}",
            })

    df = pd.DataFrame(rows)
    name = f"mc_{args.scenario}_chunk{args.chunk:02d}.csv"
    df.to_csv(out / name, index=False)
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
