"""Run one cell of the GR-GWR boundary-strength x noise phase diagram.

This is Simulation 3 prototype.  Each workflow-matrix job evaluates one
(delta, noise_sigma) cell on the same 25x25 lattice and true K=3 geometry.
K=3 is fixed here to isolate boundary detectability from K selection.

Compared models: external GWR, external MGWR, estimated GR-GWR, oracle GR-GWR.
Primary criterion is coefficient recovery against known truth; in-sample fit is
reported only as a secondary diagnostic.
"""
from __future__ import annotations

import argparse
from collections import deque
from pathlib import Path
import json
import sys

import numpy as np
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

GRID_N = 25
K_TRUE = 3
RHO_X = 0.20
LAMBDA_BOUNDARY = 0.5
MIN_REGIME_N = 6
MAX_REFINEMENT_ITER = 15

# Common background relationship plus regime contrasts.  delta scales only the
# discontinuous part so the phase diagram isolates boundary strength.
BASE_BETA = np.array([0.0, 0.5, -0.25, 0.25], dtype=float)
CONTRAST = {
    1: np.array([0.0, 1.0, -0.75, 0.25]),
    2: np.array([0.0, -1.0, 0.75, 0.25]),
    3: np.array([0.0, 0.0, 0.0, -0.75]),
}


def build_grid():
    axis = np.linspace(0.0, 1.0, GRID_N)
    xx, yy = np.meshgrid(axis, axis)
    coords = np.column_stack([xx.ravel(), yy.ravel()])
    x, y = coords[:, 0], coords[:, 1]
    true = np.empty(coords.shape[0], dtype=int)
    true[x < 0.45] = 1
    true[(x >= 0.45) & (y >= 0.50)] = 2
    true[(x >= 0.45) & (y < 0.50)] = 3
    edges = []
    for r in range(GRID_N):
        for c in range(GRID_N):
            i = r * GRID_N + c
            for dr, dc in [(0, 1), (1, -1), (1, 0), (1, 1)]:
                rr, cc = r + dr, c + dc
                if 0 <= rr < GRID_N and 0 <= cc < GRID_N:
                    edges.append((i, rr * GRID_N + cc))
    return coords, true, np.asarray(edges, dtype=int)


def adjacency_from_edges(n, edges):
    adj = [set() for _ in range(n)]
    for i, j in edges:
        adj[int(i)].add(int(j)); adj[int(j)].add(int(i))
    return adj


def conn_matrix(n, edges):
    rows = np.concatenate([edges[:, 0], edges[:, 1]])
    cols = np.concatenate([edges[:, 1], edges[:, 0]])
    return csr_matrix((np.ones(rows.size), (rows, cols)), shape=(n, n))


def boundary_mask(labels, edges):
    return labels[edges[:, 0]] != labels[edges[:, 1]]


def boundary_nodes(labels, edges):
    mask = boundary_mask(labels, edges)
    out = np.zeros(labels.size, dtype=bool)
    if np.any(mask):
        out[np.unique(edges[mask].ravel())] = True
    return out


def source_connected_after_remove(i, labels, adj):
    r = labels[i]
    remain = [int(v) for v in np.flatnonzero(labels == r) if int(v) != i]
    if not remain:
        return False
    allowed = set(remain); seen = {remain[0]}; q = deque([remain[0]])
    while q:
        u = q.popleft()
        for v in adj[u]:
            if v in allowed and v not in seen:
                seen.add(v); q.append(v)
    return len(seen) == len(remain)


def generate(delta, noise_sigma, seed, coords, true):
    rng = np.random.default_rng(seed)
    cov = np.full((3, 3), RHO_X); np.fill_diagonal(cov, 1.0)
    X_raw = rng.multivariate_normal(np.zeros(3), cov, size=coords.shape[0])
    beta_raw = np.vstack([BASE_BETA + delta * CONTRAST[int(r)] for r in true])
    signal = beta_raw[:, 0] + np.sum(beta_raw[:, 1:] * X_raw, axis=1)
    y_raw = signal + rng.normal(0.0, noise_sigma, size=coords.shape[0])

    xm = X_raw.mean(axis=0); xs = X_raw.std(axis=0, ddof=0)
    ym = float(y_raw.mean()); ys = float(y_raw.std(ddof=0))
    X = (X_raw - xm) / xs; y = (y_raw - ym) / ys
    beta_z = np.empty_like(beta_raw)
    beta_z[:, 1:] = beta_raw[:, 1:] * xs[None, :] / ys
    beta_z[:, 0] = (beta_raw[:, 0] + np.sum(beta_raw[:, 1:] * xm[None, :], axis=1) - ym) / ys
    return X, y, beta_z


def ward(slopes, edges):
    z = (slopes - slopes.mean(axis=0)) / np.maximum(slopes.std(axis=0, ddof=0), 1e-12)
    return AgglomerativeClustering(n_clusters=K_TRUE, linkage="ward", connectivity=conn_matrix(len(slopes), edges)).fit_predict(z).astype(int) + 1


def bisquare_all(d):
    if d.size == 0: return d.copy()
    bw = float(np.max(d)); bw = 1.0 if bw <= 1e-12 else float(np.nextafter(bw, np.inf))
    q = d / bw
    return np.where(q < 1.0, (1.0 - q*q)**2, 0.0)


def local_cost(i, target, labels, Xd, y, distances, adj):
    mask = labels == target; mask[i] = False
    idx = np.flatnonzero(mask)
    if idx.size < Xd.shape[1] + 1: return float("inf")
    beta, _ = BasicGWR._solve_local(Xd[idx], y[idx], bisquare_all(distances[i, idx]))
    err2 = float((y[i] - Xd[i] @ beta)**2)
    mismatch = sum(1 for j in adj[i] if labels[j] != target)
    return err2 + LAMBDA_BOUNDARY * mismatch


def objective(X, y, coords, labels, edges):
    m = RegimeAwareGWR(bandwidth="regime_size", kernel="bisquare").fit(X, y, coords, labels)
    rss = float(m.residuals_ @ m.residuals_)
    b = int(np.sum(boundary_mask(labels, edges)))
    return rss + LAMBDA_BOUNDARY*b, m


def refine(initial, X, y, coords, edges, adj):
    labels = initial.copy(); Xd = np.column_stack([np.ones(len(y)), X]); dmat = cdist(coords, coords)
    obj, model = objective(X, y, coords, labels, edges)
    accepted = 0; stop = "max_iterations"
    for it in range(1, MAX_REFINEMENT_ITER + 1):
        before = labels.copy(); old_obj = obj; changes = 0
        cand = np.flatnonzero(boundary_nodes(labels, edges))
        order = np.random.default_rng(1000 + it).permutation(cand)
        for ir in order:
            i = int(ir); cur = int(labels[i])
            if np.sum(labels == cur) <= MIN_REGIME_N: continue
            targets = sorted({cur} | {int(labels[j]) for j in adj[i]})
            if len(targets) <= 1: continue
            best = cur; best_cost = local_cost(i, cur, labels, Xd, y, dmat, adj)
            can_remove = source_connected_after_remove(i, labels, adj)
            if not can_remove: continue
            for t in targets:
                if t == cur: continue
                c = local_cost(i, t, labels, Xd, y, dmat, adj)
                if c + 1e-12 < best_cost:
                    best, best_cost = t, c
            if best != cur:
                labels[i] = best; changes += 1
        if changes == 0:
            stop = "no_label_changes"; break
        new_obj, new_model = objective(X, y, coords, labels, edges)
        if new_obj > old_obj + 1e-10:
            labels = before; stop = "global_objective_guard_rejected_sweep"; break
        obj, model = new_obj, new_model; accepted += 1
    return labels, model, accepted, stop


def align(pred, true):
    pu, tu = np.unique(pred), np.unique(true)
    C = np.zeros((len(pu), len(tu)), dtype=int)
    for a,p in enumerate(pu):
        for b,t in enumerate(tu): C[a,b] = np.sum((pred==p)&(true==t))
    rr, cc = linear_sum_assignment(-C); mp = {int(pu[r]):int(tu[c]) for r,c in zip(rr,cc)}
    return np.array([mp[int(v)] for v in pred], dtype=int)


def boundary_f1(true, pred, edges):
    t = boundary_mask(true, edges); p = boundary_mask(pred, edges)
    tp = np.sum(t&p); fp = np.sum(~t&p); fn = np.sum(t&~p)
    precision = tp/(tp+fp) if tp+fp else 1.0; recall = tp/(tp+fn) if tp+fn else 1.0
    return float(2*precision*recall/(precision+recall)) if precision+recall else 0.0


def coef_rmse(params, truth):
    # slopes only: intercept is not part of the boundary-generating contrast question.
    e = params[:,1:] - truth[:,1:]
    return float(np.sqrt(np.mean(e*e)))


def fit_rmse(y, fitted):
    return float(np.sqrt(np.mean((y-fitted)**2)))


def jump_ratio(params, truth, true, edges):
    ee = edges[boundary_mask(true, edges)]
    tj = np.abs(truth[ee[:,0],1:] - truth[ee[:,1],1:])
    ej = np.abs(params[ee[:,0],1:] - params[ee[:,1],1:])
    denom = float(np.mean(tj))
    return float(np.mean(ej)/denom) if denom > 1e-12 else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--delta", type=float, required=True)
    ap.add_argument("--noise", type=float, required=True)
    ap.add_argument("--seed", type=int, default=20260903)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    out = Path(args.output); out.mkdir(parents=True, exist_ok=True)

    coords, true, edges = build_grid(); adj = adjacency_from_edges(len(true), edges)
    X, y, beta_true = generate(args.delta, args.noise, args.seed, coords, true); yc = y[:,None]

    sel = Sel_BW(coords, yc, X, fixed=False, kernel="bisquare", constant=True)
    gwr_bw = int(round(float(sel.search())))
    gwr = GWR(coords, yc, X, gwr_bw, fixed=False, kernel="bisquare", constant=True).fit()
    gwr_params = np.asarray(gwr.params); gwr_fit = np.asarray(gwr.predy).reshape(-1)

    initial = ward(gwr_params[:,1:], edges)
    refined, gr, niter, stop = refine(initial, X, y, coords, edges, adj)
    oracle = RegimeAwareGWR(bandwidth="regime_size", kernel="bisquare").fit(X, y, coords, true)

    msel = Sel_BW(coords, yc, X, multi=True, fixed=False, kernel="bisquare", constant=True)
    bws = np.asarray(msel.search(multi_bw_min=[2], tol_multi=1e-4, max_iter_multi=60, verbose=False), dtype=float)
    mg = MGWR(coords, yc, X, msel, fixed=False, kernel="bisquare", constant=True).fit()
    mg_params = np.asarray(mg.params); mg_fit = np.asarray(mg.predy).reshape(-1)

    identifiable = args.delta > 1e-12
    row = {
        "delta": args.delta, "noise_sigma": args.noise, "seed": args.seed,
        "gwr_bandwidth": gwr_bw, "mgwr_bandwidths": [int(round(v)) for v in bws],
        "initial_ari": float(adjusted_rand_score(true, initial)) if identifiable else None,
        "refined_ari": float(adjusted_rand_score(true, refined)) if identifiable else None,
        "initial_boundary_f1": boundary_f1(true, initial, edges) if identifiable else None,
        "refined_boundary_f1": boundary_f1(true, refined, edges) if identifiable else None,
        "refinement_iterations": niter, "refinement_stop": stop,
        "gwr_coef_rmse": coef_rmse(gwr_params, beta_true),
        "mgwr_coef_rmse": coef_rmse(mg_params, beta_true),
        "grgwr_coef_rmse": coef_rmse(gr.parameters_, beta_true),
        "oracle_coef_rmse": coef_rmse(oracle.parameters_, beta_true),
        "gwr_fit_rmse": fit_rmse(y, gwr_fit),
        "mgwr_fit_rmse": fit_rmse(y, mg_fit),
        "grgwr_fit_rmse": fit_rmse(y, gr.fitted_values_),
        "oracle_fit_rmse": fit_rmse(y, oracle.fitted_values_),
        "gwr_jump_ratio": jump_ratio(gwr_params, beta_true, true, edges) if identifiable else None,
        "mgwr_jump_ratio": jump_ratio(mg_params, beta_true, true, edges) if identifiable else None,
        "grgwr_jump_ratio": jump_ratio(gr.parameters_, beta_true, true, edges) if identifiable else None,
        "oracle_jump_ratio": jump_ratio(oracle.parameters_, beta_true, true, edges) if identifiable else None,
        "note": "delta=0 has no statistically identifiable coefficient boundary; ARI/F1/jump ratios intentionally omitted" if not identifiable else "single-seed phase-diagram prototype",
    }
    tag = f"d{args.delta:g}_n{args.noise:g}".replace(".", "p")
    (out / f"scenario_{tag}.json").write_text(json.dumps(row, indent=2), encoding="utf-8")
    print(json.dumps(row, indent=2))

if __name__ == "__main__":
    main()
